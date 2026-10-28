function [candidate,Q,info]=nav_correct(s,P,obs,limit)
% Поправка общего состояния с учётом заданной зависимости ошибок.
% Возвращает кандидат; пригодность решения определяет отдельный сторож.
assert(isscalar(limit) && isfinite(limit) && limit>0, ...
    'libekf:Threshold','Expected explicit positive squared residual limit.');
paired=iscell(s);
if paired
    assert(numel(s)==2,'libekf:State','Expected two epoch states.');
    states=s; dimension=42;
    assert(strcmp(obs.frame_id,s{1}.frame_id)&&strcmp(obs.frame_id,s{2}.frame_id)&& ...
        obs.t_start==s{1}.t&&obs.t_end==s{2}.t&&s{2}.t>=s{1}.t, ...
        'libekf:Epoch','Observation must refer to both states.');
else
    states={s}; dimension=21;
    assert(strcmp(obs.frame_id,s.frame_id) && obs.t_measure==s.t, ...
        'libekf:Epoch','Observation must refer to state frame and epoch.');
end
assert(any(strcmp(obs.noise_model,{'independent','joint_required'})), ...
    'libekf:Correlation','Unknown dependence model.');
assert(~strcmp(obs.noise_model,'joint_required')||isfield(obs,'state_noise_cov'), ...
    'libekf:Correlation','Correlated prior/measurement requires cross covariance.');
assert(isequal(size(P),[dimension dimension]) && all(isfinite(P(:))) && ...
    norm(P-P.','fro')<=1e-12*max(norm(P,'fro'),realmin), ...
    'libekf:Covariance','Expected finite symmetric state covariance.');
% Проверка в безразмерном масштабе; этот профиль требует положительную определённость.
d=sqrt(diag(P)); assert(all(isfinite(d) & d>0),'libekf:Covariance','Nonpositive variance.');
A=(P./d)./d.';
if paired
    % Совпадающие ошибки сохранённого и текущего состояний дают вырожденную ковариацию.
    eigenvalues=eig((A+A.')/2);
    assert(min(eigenvalues)>=-100*eps*dimension*max(1,max(abs(eigenvalues))), ...
        'libekf:Covariance','Joint covariance has a negative variance direction.');
else
    [~,flag]=chol((A+A.')/2);
    assert(flag==0,'libekf:Covariance','State covariance is not positive definite.');
end
H=obs.H; R=obs.R; n=numel(obs.residual);
assert(isequal(size(H),[n dimension]) && all(isfinite(H(:))) && ...
    isequal(size(R),[n n]) && all(isfinite(R(:))) && ...
    norm(R-R.','fro')<=1e-12*max(norm(R,'fro'),realmin), ...
    'libekf:Covariance','Invalid observation matrices.');
[~,flag]=chol(R); assert(flag==0,'libekf:Covariance','Invalid observation covariance.');
C=zeros(dimension,n);
if isfield(obs,'state_noise_cov')
    C=obs.state_noise_cov;
    assert(isequal(size(C),[dimension n])&&all(isfinite(C(:))), ...
        'libekf:Correlation','Invalid state/measurement cross covariance.');
    assert(~strcmp(obs.noise_model,'independent')||all(C(:)==0), ...
        'libekf:Correlation','Nonzero dependence declared independent.');
    joint=[P C;C.' R]; scales=sqrt(diag(joint)); normalized=(joint./scales)./scales.';
    values=eig((normalized+normalized.')/2);
    assert(min(values)>=-100*eps*(dimension+n)*max(1,max(abs(values))), ...
        'libekf:Correlation','Inconsistent joint state/measurement covariance.');
end
% Проверяем полную зависимость с положением другого аппарата до проекции на дальность.
if isfield(obs,'peer_covariance') || isfield(obs,'peer_cross_cov')
    assert(isfield(obs,'peer_covariance')&&isfield(obs,'peer_cross_cov'), ...
        'libekf:Correlation','Incomplete peer dependence.');
    V=obs.peer_covariance; U=obs.peer_cross_cov;
    assert(isequal(size(V),[3 3])&&isequal(size(U),[dimension 3])&& ...
        all(isfinite(V(:)))&&all(isfinite(U(:)))&& ...
        norm(V-V.','fro')<=1e-12*max(norm(V,'fro'),realmin)&&all(diag(V)>=0), ...
        'libekf:Correlation','Invalid joint peer covariance.');
    joint=[P U;U.' V]; scales=sqrt(diag(joint));
    zero=scales==0;
    assert(all(all(joint(zero,:)==0)), ...
        'libekf:Correlation','Zero variance cannot have nonzero cross covariance.');
    scales(zero)=1; normalized=(joint./scales)./scales.';
    values=eig((normalized+normalized.')/2);
    assert(min(values)>=-100*eps*(dimension+3)*max(1,max(abs(values))), ...
        'libekf:Correlation','Inconsistent joint own/peer covariance.');
end
info=nav_innovation(P,obs); info.accepted=info.nis<=limit;
candidate=s; Q=P;
if ~info.accepted, return; end
% S=D*L*L'*D; решаем систему без обращения матриц.
B=(P*H.'+C)./info.scale.';
K=((B/info.factor.')/info.factor)./info.scale.';
dx=K*obs.residual(:);
for k=1:numel(states)
    states{k}=nav_apply_error(states{k},dx((k-1)*21+(1:21)));
end
if paired, candidate=states; else, candidate=states{1}; end
J=eye(dimension)-K*H;
Q=J*P*J.'+K*R*K.'-J*C*K.'-K*C.'*J.';
% После изменения ориентации переносим ошибки в новые местные оси.
G=eye(dimension);
for k=1:numel(states)
    indices=(k-1)*21+(7:9); v=dx(indices); a=norm(v); V=nav_skew(v);
    if a<1e-4
        b=1/2-a*a/24+a^4/720; c=1/6-a*a/120+a^4/5040;
    else
        b=(1-cos(a))/a^2; c=(a-sin(a))/a^3;
    end
    G(indices,indices)=eye(3)-b*V+c*(V*V);
end
Q=G*Q*G.'; Q=(Q+Q.')/2;
assert(all(isfinite(Q(:))),'libekf:Covariance','Correction covariance overflow.');
end

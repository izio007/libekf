function [candidate,Q,info]=nav_correct(s,P,obs,limit,retained)
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
% Сохранённая опора обновляется также наблюдением, которое её прямо не измеряет.
if nargin>=5
    assert(isstruct(retained)&&isscalar(retained)&& ...
        all(isfield(retained,{'position','covariance','cross_cov','frame_id','t','noise_independent'}))&& ...
        (isequal(retained.noise_independent,true)||isequal(retained.noise_independent,false))&& ...
        ~isfield(obs,'peer_covariance'), ...
        'libekf:Correlation','Retained reference requires an explicit noise dependence model.');
    if retained.noise_independent
        assert(strcmp(obs.noise_model,'independent')&& ...
            (~isfield(obs,'state_noise_cov')||all(obs.state_noise_cov(:)==0))&& ...
            (~isfield(retained,'noise_cross_cov')||all(retained.noise_cross_cov(:)==0)), ...
            'libekf:Correlation','Nonzero dependence declared independent.');
    else
        assert(strcmp(obs.noise_model,'joint_required')&&isfield(obs,'state_noise_cov')&& ...
            isfield(retained,'noise_cross_cov'),'libekf:Correlation','Explicit joint noise covariance required.');
        obs.peer_noise_cross=retained.noise_cross_cov;
        obs.peer_state_noise_cov=obs.state_noise_cov;
    end
    assert(strcmp(retained.frame_id,obs.frame_id)&&isscalar(retained.t)&& ...
        retained.t==states{end}.t,'libekf:Epoch','Retained reference must use current frame and epoch.');
    if isfield(retained,'ids'), obs.peer_ids=retained.ids; end
    obs.peer_position=retained.position; obs.peer_covariance=retained.covariance;
    obs.peer_cross_cov=retained.cross_cov; obs.peer_time=retained.t;
    obs.peer_jacobian=zeros(numel(obs.residual),numel(retained.position)); obs.peer_noise_cov=obs.R;
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
% Сжатие исходных данных в невязку не должно скрывать невозможную зависимость.
if isfield(obs,'input_noise_cov')||isfield(obs,'state_input_cov')
    assert(all(isfield(obs,{'input_noise_cov','state_input_cov','noise_transform'})), ...
        'libekf:Correlation','Incomplete raw input dependence.');
    N=obs.input_noise_cov; U=obs.state_input_cov; T=obs.noise_transform; count=size(N,1);
    assert(isequal(size(N),[count count])&&isreal(N)&&all(isfinite(N(:)))&& ...
        norm(N-N.','fro')<=1e-12*max(norm(N,'fro'),realmin)&&all(diag(N)>=0)&& ...
        isequal(size(U),[dimension count])&&isreal(U)&&all(isfinite(U(:)))&& ...
        isequal(size(T),[n count])&&isreal(T)&&all(isfinite(T(:))), ...
        'libekf:Correlation','Invalid raw input covariance or transform.');
    assert(~strcmp(obs.noise_model,'independent')||all(U(:)==0), ...
        'libekf:Correlation','Raw dependence declared independent.');
    assert(norm(R-T*N*T.','fro')<=1e-12*max(norm(R,'fro'),realmin)&& ...
        norm(C-U*T.','fro')<=1e-12*max(norm(C,'fro'),realmin), ...
        'libekf:Correlation','Projected covariance disagrees with raw inputs.');
    joint=[P U;U.' N]; scales=sqrt(diag(joint)); zero=scales==0;
    assert(~any(any(joint(zero,:)~=0)),'libekf:Correlation','Impossible zero-variance dependence.');
    scales(zero)=1; normalized=(joint./scales)./scales.';
    values=eig(full((normalized+normalized.')/2));
    assert(min(values)>=-100*eps*size(joint,1)*max(1,max(abs(values))), ...
        'libekf:Correlation','Inconsistent joint state/raw-input covariance.');
end
% Проверяем полную зависимость с положением другого аппарата до проекции на дальность.
if isfield(obs,'peer_covariance') || isfield(obs,'peer_cross_cov')
    assert(isfield(obs,'peer_covariance')&&isfield(obs,'peer_cross_cov'), ...
        'libekf:Correlation','Incomplete peer dependence.');
    V=obs.peer_covariance; U=obs.peer_cross_cov; peer_dimension=size(V,1);
    assert(peer_dimension>=3&&mod(peer_dimension,3)==0&&isequal(size(V),[peer_dimension peer_dimension])&&isequal(size(U),[dimension peer_dimension])&& ...
        all(isfinite(V(:)))&&all(isfinite(U(:)))&& ...
        norm(V-V.','fro')<=1e-12*max(norm(V,'fro'),realmin)&&all(diag(V)>=0), ...
        'libekf:Correlation','Invalid joint peer covariance.');
    joint=[P U;U.' V]; scales=sqrt(diag(joint));
    zero=scales==0;
    assert(all(all(joint(zero,:)==0)), ...
        'libekf:Correlation','Zero variance cannot have nonzero cross covariance.');
    scales(zero)=1; normalized=(joint./scales)./scales.';
    values=eig((normalized+normalized.')/2);
    assert(min(values)>=-100*eps*(dimension+peer_dimension)*max(1,max(abs(values))), ...
        'libekf:Correlation','Inconsistent joint own/peer covariance.');
end
if isfield(obs,'peer_jacobian')
    assert(all(isfield(obs,{'peer_covariance','peer_cross_cov','peer_noise_cov','peer_position','peer_time'})), ...
        'libekf:Peer','Incomplete reference update data.');
    Hp=obs.peer_jacobian; N=obs.peer_noise_cov;
    assert(isequal(size(Hp),[n peer_dimension])&&all(isfinite(Hp(:)))&& ...
        isequal(size(N),[n n])&&all(isfinite(N(:)))&& ...
        norm(N-N.','fro')<=1e-12*max(norm(N,'fro'),realmin)&& ...
        isequal(size(obs.peer_position),[peer_dimension 1])&&all(isfinite(obs.peer_position))&& ...
        isscalar(obs.peer_time)&&isfinite(obs.peer_time), ...
        'libekf:Peer','Invalid reference update data.');
    [~,flag]=chol(N); assert(flag==0,'libekf:Covariance','Invalid raw range covariance.');
    W=zeros(dimension,n); D=zeros(peer_dimension,n);
    if isfield(obs,'peer_state_noise_cov'), W=obs.peer_state_noise_cov; end
    if isfield(obs,'peer_noise_cross'), D=obs.peer_noise_cross; end
    assert(isequal(size(W),[dimension n])&&isreal(W)&&all(isfinite(W(:)))&& ...
        isequal(size(D),[peer_dimension n])&&isreal(D)&&all(isfinite(D(:))), ...
        'libekf:Correlation','Invalid joint reference/noise dependence.');
    joint_noise=[P obs.peer_cross_cov W;obs.peer_cross_cov.' obs.peer_covariance D;W.' D.' N];
    scales=sqrt(diag(joint_noise)); zero=scales==0;
    assert(all(all(joint_noise(zero,:)==0)),'libekf:Correlation','Zero variance has nonzero noise dependence.');
    scales(zero)=1; normalized=(joint_noise./scales)./scales.';
    values=eig((normalized+normalized.')/2);
    assert(min(values)>=-100*eps*size(joint_noise,1)*max(1,max(abs(values))), ...
        'libekf:Correlation','Inconsistent joint navigation/reference/noise covariance.');
    expectedR=N+Hp*obs.peer_covariance*Hp.'+Hp*D+D.'*Hp.';
    expectedC=obs.peer_cross_cov*Hp.'+W;
    assert(norm(R-expectedR,'fro')<=1e-12*max(norm(R,'fro'),realmin)&& ...
        norm(C-expectedC,'fro')<=1e-12*max(norm(C,'fro'),realmin), ...
        'libekf:Correlation','Reference and effective observation errors disagree.');
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
info.error_reset=G;
assert(all(isfinite(Q(:))),'libekf:Covariance','Correction covariance overflow.');
% Совместно уточняем опорный объект, чтобы его ошибка не стала независимой снова.
if isfield(obs,'peer_jacobian')
    Hp=obs.peer_jacobian; V=obs.peer_covariance; U=obs.peer_cross_cov;
    assert(isequal(size(Hp),[n peer_dimension])&&all(isfinite(Hp(:)))&& ...
        isequal(size(obs.peer_noise_cov),[n n]),'libekf:Peer','Invalid peer observation geometry.');
    Bp=(U.'*H.'+V*Hp.'+D)./info.scale.';
    Kp=((Bp/info.factor.')/info.factor)./info.scale.';
    prior=[P U;U.' V]; gain=[K;Kp]; A=eye(dimension+peer_dimension)-gain*[H Hp];
    noise_cross=[W;D];
    joint=A*prior*A.'+gain*obs.peer_noise_cov*gain.'- ...
        A*noise_cross*gain.'-gain*noise_cross.'*A.';
    reset=blkdiag(G,eye(peer_dimension)); joint=reset*joint*reset.'; joint=(joint+joint.')/2;
    assert(all(isfinite(joint(:))),'libekf:Covariance','Peer correction overflow.');
    info.peer_update=struct('position',obs.peer_position+Kp*obs.residual(:), ...
        'covariance',joint(dimension+(1:peer_dimension),dimension+(1:peer_dimension)), ...
        'cross_cov',joint(1:dimension,dimension+(1:peer_dimension)), ...
        'frame_id',obs.frame_id,'t',obs.peer_time);
    if isfield(obs,'peer_ids'), info.peer_update.ids=obs.peer_ids; end
end

end

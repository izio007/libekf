function [candidate,Q,info]=nav_correct(s,P,obs,limit)
% Поправка общего состояния по независимому от априорной ошибки наблюдению.
% Возвращает кандидат; пригодность решения определяет отдельный сторож.
assert(isscalar(limit) && isfinite(limit) && limit>0, ...
    'libekf:Threshold','Expected explicit positive squared residual limit.');
assert(strcmp(obs.frame_id,s.frame_id) && obs.t_measure==s.t, ...
    'libekf:Epoch','Observation must refer to state frame and epoch.');
assert(strcmp(obs.noise_model,'independent'), ...
    'libekf:Correlation','Correlated prior/measurement requires a joint model.');
assert(isequal(size(P),[21 21]) && all(isfinite(P(:))) && ...
    norm(P-P.','fro')<=1e-12*max(norm(P,'fro'),realmin), ...
    'libekf:Covariance','Expected finite symmetric state covariance.');
% Проверка в безразмерном масштабе; этот профиль требует положительную определённость.
d=sqrt(diag(P)); assert(all(isfinite(d) & d>0),'libekf:Covariance','Nonpositive variance.');
A=(P./d)./d.'; [~,flag]=chol((A+A.')/2);
assert(flag==0,'libekf:Covariance','State covariance is not positive definite.');
H=obs.H; R=obs.R; n=numel(obs.residual);
assert(isequal(size(H),[n 21]) && all(isfinite(H(:))) && ...
    isequal(size(R),[n n]) && all(isfinite(R(:))) && ...
    norm(R-R.','fro')<=1e-12*max(norm(R,'fro'),realmin), ...
    'libekf:Covariance','Invalid observation matrices.');
[~,flag]=chol(R); assert(flag==0,'libekf:Covariance','Invalid observation covariance.');
info=nav_innovation(P,obs); info.accepted=info.nis<=limit;
candidate=s; Q=P;
if ~info.accepted, return; end
% S=D*L*L'*D; решаем систему без обращения матриц.
B=(P*H.')./info.scale.';
K=((B/info.factor.')/info.factor)./info.scale.';
dx=K*obs.residual(:);
candidate=nav_apply_error(s,dx);
J=eye(21)-K*H;
Q=J*P*J.'+K*R*K.';
% После изменения ориентации переносим ошибки в новые местные оси.
v=dx(7:9); a=norm(v); V=nav_skew(v);
if a<1e-4
    b=1/2-a*a/24+a^4/720; c=1/6-a*a/120+a^4/5040;
else
    b=(1-cos(a))/a^2; c=(a-sin(a))/a^3;
end
G=eye(21); G(7:9,7:9)=eye(3)-b*V+c*(V*V);
Q=G*Q*G.'; Q=(Q+Q.')/2;
assert(all(isfinite(Q(:))),'libekf:Covariance','Correction covariance overflow.');
end

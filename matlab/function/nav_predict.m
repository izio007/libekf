function [out,Q]=nav_predict(s,P,imu,env)
% Инерциальное продолжение в неподвижных местных осях, связанных с Землёй.
% Ускорение — удельная сила, gyro — угловая скорость относительно инерциальных осей.
% noise_density: спектральная матрица белого шума [accel;gyro], а не дисперсия отсчёта.
assert(strcmp(imu.quality,'VALID'),'libekf:Quality','Invalid inertial data.');
assert(strcmp(imu.frame,'BODY') && strcmp(env.frame_id,s.frame_id), ...
    'libekf:Frame','Expected body inertial data and matching Earth-fixed frame.');
assert(strcmp(imu.accel_unit,'m/s^2') && strcmp(imu.gyro_unit,'rad/s'), ...
    'libekf:Units','Invalid inertial units.');
dt=imu.t_end-imu.t_start;
assert(imu.t_start==s.t && isfinite(dt) && dt>0 && ...
    isscalar(env.max_step) && isfinite(env.max_step) && env.max_step>0 && dt<=env.max_step && ...
    isfinite(imu.t_received) && imu.t_received>=imu.t_end, ...
    'libekf:Epoch','Invalid inertial interval, arrival time or step limit.');
assert(ischar(imu.source) && ~isempty(imu.source) && iscell(imu.lineage) && ...
    ~isempty(imu.lineage) && all(cellfun(@(x)ischar(x)&&~isempty(x),imu.lineage)), ...
    'libekf:Lineage','Missing inertial provenance.');
f=vec(imu.accel)-vec(s.ba); w=vec(imu.gyro)-vec(s.bg);
O=nav_skew(vec(env.Omega)); C=nav_rotation(s.q);
g=vec(env.gravity_origin); G=env.gravity_gradient;
assert(isequal(size(G),[3 3]) && all(isfinite(G(:))) && ...
    isscalar(env.gravity_radius) && isfinite(env.gravity_radius) && env.gravity_radius>0 && ...
    norm(s.p)<=env.gravity_radius,'libekf:Gravity','Outside gravity model domain.');
psd(P,21); psd(imu.noise_density,6); psd(env.state_noise_density,21);
% Постоянные показания за короткий интервал; ориентация в середине шага.
Cm=expm(-O*dt/2)*C*expm(nav_skew(w)*dt/2);
out=s; out.q=nav_quaternion(expm(-O*dt)*C*expm(nav_skew(w)*dt));
% При замороженной средней удельной силе линейное движение решается точно.
A=[zeros(3) eye(3) zeros(3,1); G -2*O Cm*f+g; zeros(1,7)];
pv=expm(A*dt)*[s.p(:);s.v(:);1]; out.p=pv(1:3); out.v=pv(4:6);
assert(norm(out.p)<=env.gravity_radius,'libekf:Gravity','Step leaves gravity model domain.');
out.clock=s.clock+s.drift*dt; out.t=imu.t_end;
% Линеаризация местной ошибки около середины шага.
F=zeros(21); F(1:3,4:6)=eye(3); F(4:6,1:3)=G; F(4:6,4:6)=-2*O;
F(4:6,7:9)=-Cm*nav_skew(f); F(4:6,10:12)=-Cm;
F(7:9,7:9)=-nav_skew(w); F(7:9,13:15)=-eye(3); F(20,21)=1;
B=zeros(21,6); B(4:6,1:3)=-Cm; B(7:9,4:6)=-eye(3);
W=B*imu.noise_density*B.'+env.state_noise_density;
% Совместное интегрирование перехода и накопленной погрешности.
E=expm([F W;zeros(21) -F.']*dt); T=E(1:21,1:21);
Q=T*P*T.'+E(1:21,22:42)*T.'; Q=(Q+Q.')/2;
assert(all(isfinite(Q(:))) && all(isfinite([out.p;out.v;out.clock])), ...
    'libekf:Numerics','Inertial propagation overflow.');
end
function x=vec(x)
assert(isnumeric(x)&&numel(x)==3&&all(isfinite(x(:))), ...
    'libekf:Vector','Expected finite 3-vector.'); x=x(:);
end
function psd(M,n)
assert(isequal(size(M),[n n]) && all(isfinite(M(:))) && ...
    norm(M-M.','fro')<=1e-12*max(norm(M,'fro'),realmin), ...
    'libekf:Covariance','Invalid covariance dimensions or symmetry.');
assert(min(eig((M+M.')/2))>=0,'libekf:Covariance','Negative variance direction.');
end

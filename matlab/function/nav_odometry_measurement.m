function obs=nav_odometry_measurement(a,b,m)
% Перемещение или направление движения камеры между двумя состояниями; H имеет 42 столбца.
required={'source','lineage','quality','noise_model','correlation_group','frame_id', ...
    't_start','t_end','t_received','units','z','R','camera_to_body','lever', ...
    'scale_valid','calibration_valid'};
for k=1:numel(required)
    assert(isfield(m,required{k}),'libekf:Packet',['Missing ' required{k}]);
end
mode='metric'; if isfield(m,'mode'), mode=m.mode; end
if any(strcmp(mode,{'pose','direction_rotation'}))
    count=3; first_mode='metric';
    if strcmp(mode,'direction_rotation'), count=2; first_mode='direction'; end
    total=count+3; raw=m.R;
    assert(isnumeric(m.z)&&numel(m.z)==7&&all(isfinite(m.z(:)))&& ...
        iscell(m.units)&&numel(m.units)==total,'libekf:Dimension','Invalid paired camera observation.');
    assert(isequal(size(raw),[total total])&&all(isfinite(raw(:)))&& ...
        norm(raw-raw.','fro')<=1e-12*max(norm(raw,'fro'),realmin), ...
        'libekf:Covariance','Invalid full pose covariance.');
    [~,flag]=chol(raw); assert(flag==0,'libekf:Covariance','Invalid full pose covariance.');
    first=m; first.mode=first_mode; first.z=m.z(1:3);
    first.R=raw(1:count,1:count); first.units=m.units(1:count);
    second=m; second.mode='rotation'; second.z=m.z(4:7);
    second.R=raw(count+1:end,count+1:end); second.units=m.units(count+1:end);
    if isfield(m,'state_noise_cov')
        first=rmfield(first,'state_noise_cov'); second=rmfield(second,'state_noise_cov');
    end
    one=nav_odometry_measurement(a,b,first); two=nav_odometry_measurement(a,b,second);
    obs=one; obs.predicted=[one.predicted;two.predicted];
    obs.residual=[one.residual;two.residual]; obs.H=[one.H;two.H];
    obs.noise_transform=blkdiag(one.noise_transform,two.noise_transform);
    obs.R=obs.noise_transform*raw*obs.noise_transform.'; obs.R=(obs.R+obs.R.')/2;
    obs.units=m.units; obs=attach_dependence(obs,m); return;
end
assert(any(strcmp(mode,{'metric','direction','rotation'})),'libekf:Model','Unknown odometry mode.');
assert(strcmp(m.quality,'VALID') && ...
    (~strcmp(mode,'metric')||isequal(m.scale_valid,true)) && ...
    isequal(m.calibration_valid,true),'libekf:Quality','Valid metric scale and calibration required.');
assert(strcmp(a.frame_id,b.frame_id)&&strcmp(m.frame_id,a.frame_id), ...
    'libekf:Frame','Both states must share a local frame.');
assert(isscalar(m.t_start)&&isscalar(m.t_end)&&isscalar(m.t_received)&& ...
    all(isfinite([m.t_start m.t_end m.t_received]))&& ...
    m.t_start==a.t&&m.t_end==b.t&&b.t>a.t&&m.t_received>=b.t, ...
    'libekf:Epoch','Invalid image pair times.');
assert(ischar(m.source)&&~isempty(m.source)&&iscell(m.lineage)&&~isempty(m.lineage)&& ...
    all(cellfun(@(x)ischar(x)&&~isempty(x),m.lineage))&& ...
    ischar(m.correlation_group)&&~isempty(m.correlation_group)&& ...
    any(strcmp(m.noise_model,{'independent','joint_required'})), ...
    'libekf:Lineage','Missing provenance or dependence declaration.');
units={'m';'m';'m'}; dimension=3; zdimension=3;
if strcmp(mode,'rotation'), units={'rad';'rad';'rad'}; zdimension=4; end
if strcmp(mode,'direction'), units={'1';'1'}; dimension=2; end
assert(iscell(m.units)&&isequal(m.units(:),units), ...
    'libekf:Units','Invalid displacement or tangent-plane units.');
D=m.camera_to_body; l=m.lever(:); z=m.z(:);
assert(isnumeric(D)&&isequal(size(D),[3 3])&&all(isfinite(D(:)))&& ...
    norm(D.'*D-eye(3),'fro')<1e-10&&abs(det(D)-1)<1e-10&& ...
    isnumeric(l)&&numel(l)==3&&all(isfinite(l))&& ...
    isnumeric(z)&&numel(z)==zdimension&&all(isfinite(z)), ...
    'libekf:Geometry','Invalid camera geometry or displacement.');
noise_transform=eye(dimension);
R=m.R;
assert(isequal(size(R),[dimension dimension])&&all(isfinite(R(:)))&& ...
    norm(R-R.','fro')<=1e-12*max(norm(R,'fro'),realmin), ...
    'libekf:Covariance','Invalid displacement covariance.');
[~,flag]=chol(R); assert(flag==0,'libekf:Covariance','Positive covariance required.');
A=nav_rotation(a.q); B=nav_rotation(b.q);
delta=b.p(:)-a.p(:)+B*l-A*l; v=A.'*delta; h=D.'*v;
H=zeros(3,42); H(:,1:3)=-D.'*A.'; H(:,22:24)=D.'*A.';
H(:,7:9)=D.'*nav_skew(v+l); H(:,28:30)=-D.'*A.'*B*nav_skew(l);
if strcmp(mode,'direction')
    for name={'tangent_basis','min_displacement','max_angle'}
        assert(isfield(m,name{1}),'libekf:Packet',['Missing ' name{1}]);
    end
    E=m.tangent_basis; length=norm(h);
    assert(abs(norm(z)-1)<1e-10 && isequal(size(E),[3 2]) && all(isfinite(E(:))) && ...
        norm(E.'*E-eye(2),'fro')<1e-10 && norm(E.'*z)<1e-10 && ...
        isscalar(m.min_displacement)&&isfinite(m.min_displacement)&&m.min_displacement>0 && ...
        isscalar(m.max_angle)&&isfinite(m.max_angle)&&m.max_angle>0&&m.max_angle<pi/2, ...
        'libekf:Geometry','Invalid direction or tangent-plane domain.');
    assert(length>m.min_displacement,'libekf:Geometry','Translation direction undefined.');
    direction=h/length;
    assert(z.'*direction>=cos(m.max_angle),'libekf:Geometry','Outside local angular domain.');
    H=E.'*(eye(3)-direction*direction.')/length*H;
    h=E.'*direction; z=E.'*z;
end
if strcmp(mode,'rotation')
    assert(isfield(m,'max_angle')&&isscalar(m.max_angle)&&isfinite(m.max_angle)&& ...
        m.max_angle>0&&m.max_angle<pi,'libekf:Geometry','Specify local rotation domain.');
    relative=D.'*A.'*B*D; measured=nav_rotation(z);
    errorq=nav_quaternion(relative.'*measured);
    if errorq(1)<0, errorq=-errorq; end
    sine=norm(errorq(2:4)); angle=2*atan2(sine,errorq(1));
    assert(angle<=m.max_angle,'libekf:Geometry','Outside local rotation domain.');
    if sine<1e-12, residual=2*errorq(2:4); else, residual=angle/sine*errorq(2:4); end
    V=nav_skew(residual);
    if angle<1e-4, coefficient=1/12+angle^2/720; ...
    else, coefficient=(1-.5*angle/tan(.5*angle))/angle^2; end
    left=eye(3)-.5*V+coefficient*(V*V);
    right=eye(3)+.5*V+coefficient*(V*V);
    H=zeros(3,42); H(:,7:9)=-left*relative.'*D.'; H(:,28:30)=left*D.';
    noise_transform=right; R=right*R*right.'; R=(R+R.')/2;
    h=nav_quaternion(relative);
else
    residual=z-h;
end
assert(all(isfinite(h))&&all(isfinite(H(:))),'libekf:Geometry','Invalid state geometry.');
obs=struct('noise_transform',noise_transform,'predicted',h,'residual',residual,'H',H,'R',R, ...
    'units',{m.units},'frame_id',m.frame_id,'t_start',a.t,'t_end',b.t, ...
    't_received',m.t_received,'source',m.source,'lineage',{m.lineage}, ...
    'noise_model',m.noise_model,'correlation_group',m.correlation_group);
obs=attach_dependence(obs,m);
end

function obs=attach_dependence(obs,m)
if isfield(m,'state_noise_cov')
    U=m.state_noise_cov;
    assert(isfield(m,'state_noise_frame_id')&&strcmp(m.state_noise_frame_id,obs.frame_id)&& ...
        isequal(size(U),[42 numel(obs.residual)])&&all(isfinite(U(:))), ...
        'libekf:Correlation','Camera cross covariance must cover both epoch errors.');
    obs.state_noise_cov=U*obs.noise_transform.';
end
end

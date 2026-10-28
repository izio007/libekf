function obs = nav_measurement(s, m, env)
% Build one physical observation and analytic local-error Jacobian.
% All working geometry is LOCAL with a named Earth-fixed anchor.
% No state mutation, gating, fusion, or assumption of source independence here.
required={'source','model','z','R','units','frame','frame_id','t_measure','t_received', ...
    'quality','lineage','correlation_group','noise_model'};
for k=1:numel(required)
    assert(isfield(m,required{k}), 'libekf:Packet', ['Missing ' required{k}]);
end
assert(ischar(m.source) && ~isempty(m.source) && ischar(m.model), 'libekf:Packet','Invalid source/model.');
assert(strcmp(m.quality,'VALID'), 'libekf:Quality','Only VALID input enters observation construction.');
assert(strcmp(m.frame,'LOCAL') && strcmp(m.frame_id,s.frame_id) && strcmp(s.frame_id,env.frame_id), ...
    'libekf:Frame','Packet, state and environment must share one local frame ID.');
assert(isscalar(s.t) && isfinite(s.t) && isscalar(m.t_measure) && isfinite(m.t_measure) && ...
    isscalar(m.t_received) && isfinite(m.t_received) && ...
    isfield(env,'time_tolerance') && isscalar(env.time_tolerance) && ...
    isfinite(env.time_tolerance) && env.time_tolerance>=0, 'libekf:Time','Invalid times/tolerance.');
assert(abs(s.t-m.t_measure)<=env.time_tolerance, 'libekf:Epoch','State must be at measurement epoch.');
assert(m.t_received+env.time_tolerance>=m.t_measure,'libekf:Time','Receive time precedes measurement.');
assert(iscell(m.lineage) && ~isempty(m.lineage) && all(cellfun(@(x)ischar(x)&&~isempty(x),m.lineage)), ...
    'libekf:Lineage','Provide original information record IDs.');
assert(ischar(m.correlation_group) && ~isempty(m.correlation_group) && ...
    any(strcmp(m.noise_model,{'independent','joint_required'})), ...
    'libekf:Correlation','Declare correlation group and noise assumption.');
C=nav_rotation(s.q); Rtransform=[]; z=m.z(:);
fields={'p','v','ba','bg','wind'};
for k=1:numel(fields), vector3(s.(fields{k})); end
assert(all(isfinite([s.bh s.clock s.drift])) && ...
    all(cellfun(@(f)isscalar(s.(f)),{'bh','clock','drift'})), 'libekf:State','Invalid augmented state.');
H=zeros(0,21); angular=[];
lever=zeros(3,1);
if isfield(m,'lever'), lever=vector3(m.lever); end
ps=s.p(:)+C*lever;
Jp=zeros(3,21); Jp(:,1:3)=eye(3); Jp(:,7:9)=-C*nav_skew(lever);
switch m.model
    case {'gnss_position','external_position'}
        h=ps; H=Jp; units=repmat({'m'},3,1);
    case {'gnss_velocity','external_velocity','gnss_position_velocity','external_position_velocity','range_rate'}
        if norm(lever)>0
            assert(isfield(m,'omega_ib') && isfield(env,'Omega'), 'libekf:Lever','Nonzero lever needs body rates and Earth rate.');
            omega=vector3(m.omega_ib)-s.bg(:); Omega=vector3(env.Omega);
        else
            omega=zeros(3,1); Omega=zeros(3,1);
        end
        rot=cross(omega,lever);
        va=s.v(:)+C*rot-cross(Omega,C*lever);
        Jv=zeros(3,21); Jv(:,4:6)=eye(3);
        Jv(:,7:9)=-C*nav_skew(rot)+nav_skew(Omega)*C*nav_skew(lever);
        Jv(:,13:15)=C*nav_skew(lever);
        if strcmp(m.model,'range_rate')
            assert(isfield(m,'geometry_at_epoch') && isequal(m.geometry_at_epoch,true), ...
                'libekf:GeometryEpoch','Range-rate geometry must be prepared at observation epoch.');
            d=ps-vector3(m.station); rho=norm(d);
            assert(rho>0,'libekf:Geometry','Zero satellite range.'); u=d/rho;
            dv=va-vector3(m.station_velocity);
            h=u.'*dv+s.drift+m.rate_correction;
            H=((eye(3)-u*u.')*dv/rho).'*Jp+u.'*Jv; H(21)=H(21)+1;
            units={'m/s'};
        elseif any(strcmp(m.model,{'gnss_position_velocity','external_position_velocity'}))
            h=[ps;va]; H=[Jp;Jv]; units=[repmat({'m'},3,1);repmat({'m/s'},3,1)];
        else
            h=va; H=Jv; units=repmat({'m/s'},3,1);
        end
    case {'dme_range','pseudorange'}
        d=ps-vector3(m.station); rho=norm(d);
        assert(rho>0,'libekf:Geometry','Undefined range at station.');
        h=rho; H=d.'/rho*Jp; units={'m'};
        if strcmp(m.model,'pseudorange')
            assert(isfield(m,'geometry_at_epoch') && isequal(m.geometry_at_epoch,true), ...
                'libekf:GeometryEpoch','Satellite transmit/receive geometry must be prepared.');
            h=h+s.clock+m.range_correction; H(20)=H(20)+1;
        end
    case 'vor_azimuth'
        B=basis(m.basis); d=ps-vector3(m.station);
        east=B(:,1).'*d; north=B(:,2).'*d; r2=east^2+north^2;
        assert(r2>0,'libekf:Geometry','Undefined azimuth above station.');
        h=atan2(east,north)+m.azimuth_offset;
        H=(north*B(:,1).'-east*B(:,2).')/r2*Jp;
        units={'rad'}; angular=1;
    case 'baro_height'
        B=basis(env.C_e_l);
        [height,up_e]=nav_height_local(ps,env);
        up=B.'*up_e;
        h=height-m.geoid_offset+s.bh; H=up.'*Jp; H(19)=H(19)+1; units={'m'};
    case 'airdata'
        vb=C.'*(s.v(:)-s.wind(:)); speed=norm(vb); x=vb(1); y=vb(2); z=vb(3);
        r2=x*x+z*z; d=sqrt(r2);
        assert(speed>0 && d>0,'libekf:Geometry','Air angles undefined at zero/vertical lateral velocity.');
        h=[speed;atan2(z,x);atan2(y,d)];
        D=[vb.'/speed; -z/r2 0 x/r2; -y*x/(speed^2*d) d/speed^2 -y*z/(speed^2*d)];
        J=zeros(3,21); J(:,4:6)=C.'; J(:,7:9)=nav_skew(vb); J(:,16:18)=-C.';
        H=D*J; units={'m/s';'rad';'rad'}; angular=[2 3];
    case {'radar_range','radar_terrain','radar_footprint'}
        beam=unit(m.beam_body);
        assert(isfield(m,'coverage_complete') && isequal(m.coverage_complete,true) && ...
            isfield(m,'terrain_quality') && strcmp(m.terrain_quality,'VALID'), ...
            'libekf:TerrainCoverage','Complete, valid terrain coverage of the footprint is required.');
        assert(isfield(m,'nearest_return_valid') && isequal(m.nearest_return_valid,true), ...
            'libekf:RadarFootprint','Nearest-point return law must be validated for this source.');
        hit=nav_terrain_footprint(ps,C*beam,m.half_angle,m.vertices,m.faces,m.edge_margin);
        assert(hit.smooth,'libekf:TerrainNonsmooth','Closest return is not a unique smooth face interior.');
        assert(isscalar(m.max_valid_range) && isfinite(m.max_valid_range) && m.max_valid_range>0 && ...
            m.z<=m.max_valid_range && hit.range<=m.max_valid_range, ...
            'libekf:RadarRange','Radar return exceeds the configured valid range.');
        h=hit.range; H=-hit.normal.'/hit.direction_cosine*Jp;
        if hit.cone_active
            k=hit.normal.'*C*beam; transverse=sqrt(max(0,1-k*k));
            derivative=cos(m.half_angle)-k/transverse*sin(m.half_angle);
            H(7:9)=H(7:9)+hit.plane_distance/hit.direction_cosine^2*derivative*hit.normal.'*C*nav_skew(beam);
        end
        units={'m'};
    case 'ils_angles'
        B=basis(m.basis); d=B.'*(ps-vector3(m.station));
        x=d(1); y=d(2); z=d(3); horizontal=hypot(x,y); distance2=d.'*d;
        assert(x>0 && horizontal>0, 'libekf:Geometry','Outside calibrated forward approach sector.');
        h=[atan2(y,x);atan2(z,horizontal)]-m.angle_reference(:);
        D=[-y/(x*x+y*y) x/(x*x+y*y) 0; ...
            -z*x/(distance2*horizontal) -z*y/(distance2*horizontal) horizontal/distance2];
        H=D*B.'*Jp; units={'rad';'rad'}; angular=[1 2];
    case 'gls_deviation'
        M=m.projection;
        assert(isequal(size(M),[2 3]) && all(isfinite(M(:))) && norm(M*M.'-eye(2),'fro')<1e-10, ...
            'libekf:Geometry','GLS projection needs orthonormal lateral/vertical rows.');
        h=M*(ps-vector3(m.path_origin)); H=M*Jp; units={'m';'m'};
    case 'attitude'
        nav_rotation(m.z); measured=m.z(:); q=s.q(:);
        % Log(q_nominal^{-1} * q_measured), invariant under quaternion sign.
        qe=[q(1)*measured(1)+q(2:4).'*measured(2:4); ...
            q(1)*measured(2:4)-measured(1)*q(2:4)-cross(q(2:4),measured(2:4))];
        if qe(1)<0, qe=-qe; end
        a=2*atan2(norm(qe(2:4)),qe(1));
        assert(isfield(m,'max_angle') && isscalar(m.max_angle) && isfinite(m.max_angle) && ...
            m.max_angle>0 && m.max_angle<pi && a<=m.max_angle, ...
            'libekf:AttitudeDomain','Attitude discrepancy exceeds validated local chart.');
        if a<1e-8, z=2*qe(2:4); coeff=1/12+a*a/720;
        else, z=a*qe(2:4)/norm(qe(2:4)); coeff=(1-0.5*a*cot(a/2))/(a*a); end
        A=nav_skew(z);
        H=zeros(3,21); H(:,7:9)=eye(3)-0.5*A+coeff*A*A;
        % Packet R describes right-local noise at the measured attitude.
        Rtransform=eye(3)+0.5*A+coeff*A*A;
        h=zeros(3,1); units=repmat({'rad'},3,1);
    otherwise
        error('libekf:Model','Unsupported observation model: %s',m.model);
end
if ~strcmp(m.model,'attitude'), z=m.z(:); end
units=units(:);
assert(numel(z)==numel(h) && all(isfinite(z)) && all(isfinite(h)) && all(isfinite(H(:))), ...
    'libekf:Dimension','Invalid observation dimension/value.');
assert(iscell(m.units) && isequal(m.units(:),units), 'libekf:Units','Observation units differ from model.');
R=m.R;
assert(isequal(size(R),[numel(h) numel(h)]) && all(isfinite(R(:))) && ...
    norm(R-R.','fro')<=1e-12*max(norm(R,'fro'),realmin), 'libekf:Covariance','R must be finite and symmetric.');
[~,flag]=chol(R);
assert(flag==0,'libekf:Covariance','This observation profile requires positive definite R.');
if ~isempty(Rtransform), R=Rtransform*R*Rtransform.'; end
residual=z-h;
residual(angular)=atan2(sin(residual(angular)),cos(residual(angular)));
obs=struct('predicted',h,'residual',residual,'H',H,'R',R,'units',{units}, ...
    'source',m.source,'model',m.model,'t_measure',m.t_measure,'t_received',m.t_received, ...
    'frame_id',m.frame_id, ...
    'lineage',{m.lineage},'correlation_group',m.correlation_group,'noise_model',m.noise_model);
end

function v=vector3(v)
assert(isnumeric(v) && numel(v)==3 && all(isfinite(v(:))), 'libekf:Vector','Expected finite 3-vector.'); v=v(:);
end
function u=unit(u)
u=vector3(u); assert(abs(norm(u)-1)<1e-10,'libekf:Geometry','Expected unit direction.');
end
function B=basis(B)
assert(isequal(size(B),[3 3]) && all(isfinite(B(:))) && norm(B.'*B-eye(3),'fro')<1e-10, ...
    'libekf:Geometry','Expected orthonormal geometry basis.');
end

function status=run_tests
% Проверки требований к поведению; один файл, без внешних настроек и отчётов.
base=fileparts(mfilename('fullpath')); addpath(fullfile(base,'function'));
cases={@requirement_source_models,@requirement_coordinates_and_terrain, ...
    @requirement_numerical_invariance,@requirement_broad_radar,@requirement_common_correction, ...
    @requirement_source_provenance,@requirement_inertial_motion,@requirement_watch_intervals, ...
    @requirement_terrain_sequence,@requirement_all_source_fusion};
names={'Показания соответствуют физическим моделям источников', ...
    'Смена координат сохраняет результат; неоднозначность рельефа не скрывается', ...
    'Мера расхождения сохраняется при смене единиц и масштаба', ...
    'Радиовысотомер ищет ближайшую поверхность во всей диаграмме', ...
    'Общая поправка учитывает погрешности и сохраняет состояние при отклонении', ...
    'Повторные и неисправные источники не изменяют общую оценку', ...
    'Инерциальные данные продолжают движение и накапливают погрешность', ...
    'Три интервала контроля различают тревогу, согласие и нехватку истории', ...
    'Привязка использует последовательность и сохраняет неразличимые варианты', ...
    'Все модели уточняют общую оценку независимо от осей прибора'};
checks=repmat(struct('requirement','','passed',false,'error',''),1,numel(cases));
for k=1:numel(cases)
    checks(k).requirement=names{k};
    try
        cases{k}(); checks(k).passed=true;
        fprintf('Выполнено: %s\n',names{k});
    catch ex
        checks(k).error=getReport(ex,'extended','hyperlinks','off');
        fprintf(2,'Не выполнено: %s\n%s\n',names{k},checks(k).error);
    end
end
status=struct('passed',all([checks.passed]),'checks',checks);
end
function requirement_all_source_fusion
[s,e]=fixture;
models={'gnss_position','external_position','gnss_velocity','external_velocity', ...
    'gnss_position_velocity','external_position_velocity','dme_range','pseudorange', ...
    'range_rate','vor_azimuth','baro_height','airdata','radar_range','ils_angles', ...
    'gls_deviation','attitude','radar_footprint'};
sensor0=struct('id',s.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
R=nav_rotation([cos(.4);sin(.4)/sqrt(3)*ones(3,1)]); offset=[10;20;-30];
sensor1=struct('id','rotated_sensor','origin_ecef',e.origin_ecef+e.C_e_l*offset, ...
    'C_e_s',e.C_e_l*R);
adapter=struct('C_e_l',R,'origin_ecef',offset,'frame_id',sensor1.id);
f=struct('state',s,'P',eye(21),'used',{{}},'capacity',100);
for k=1:numel(models)
    m=packet(s,models{k}); m.frame='SENSOR_LOCAL';
    [a,ra]=nav_assimilate(f,m,sensor0,e,1e8);
    assert(ra.accepted,'libekf:AllSources','%s rejected: %s',models{k},ra.reason);
    % Те же физические данные, записанные в иных местных осях прибора.
    original=m; original.frame='ECEF';
    native=nav_packet_from_ecef(original,adapter); native.frame='SENSOR_LOCAL';
    [b,rb]=nav_assimilate(f,native,sensor1,e,1e8);
    assert(rb.accepted,'libekf:AllSources','Rotated %s rejected: %s',models{k},rb.reason);
    close(a.state.p,b.state.p,2e-7); close(a.state.v,b.state.v,2e-7);
    close(nav_rotation(a.state.q),nav_rotation(b.state.q),2e-7);
    close(a.P,b.P,2e-7); close(ra.nis,rb.nis,2e-6);
    assert(isequal(a.used,b.used) && a.state.t==s.t);
end
% Совместная коррекция положения разными независимыми приборами имеет аналитический ответ.
s.q=[1;0;0;0]; f.state=s;
m=packet(s,'gnss_position'); m.frame='SENSOR_LOCAL'; m.lever=zeros(3,1);
m.z=s.p+[0;2;4];
[a,r]=nav_assimilate(f,m,sensor0,e,100); assert(r.accepted);
g=packet(s,'gls_deviation'); g.frame='SENSOR_LOCAL'; g.lever=zeros(3,1);
g.path_origin=zeros(3,1); g.z=s.p(2:3)+[4;2];
[b,r]=nav_assimilate(a,g,sensor0,e,100); assert(r.accepted);
close(b.state.p,s.p+[0;2;2],1e-12);
close(diag(b.P(1:3,1:3)),[.5;1/3;1/3],1e-12);
end
function requirement_terrain_sequence
map=struct('frame_id','terrain','coverage_complete',true,'quality','VALID', ...
    'nearest_return_valid',true,'max_valid_range',2000,'lever',zeros(3,1), ...
    'beam_body',[0;0;-1],'half_angle',pi/3, ...
    'vertices',[-100 -100 0;100 -100 0;0 100 0],'faces',[1 2 3]);
series=struct('t',[0 1 2],'z',[10 10 10],'unit','m','quality','VALID','source','RA', ...
    'lineage',{{'ra0','ra1','ra2'}},'covariance',.01*eye(3));
paths=cell(3,3);
for k=1:3
    paths{1,k}=struct('p',[k;0;10],'q',[1;0;0;0],'t',k-1,'frame_id','terrain');
    paths{2,k}=paths{1,k}; paths{2,k}.p(1)=paths{2,k}.p(1)+10;
    paths{3,k}=paths{1,k};
end
paths{3,3}.p(3)=11;
limits=struct('absolute',9,'difference',2);
r=nav_terrain_sequence(paths,series,map,limits);
assert(strcmp(r.status,'AMBIGUOUS')&&isequal(r.candidates,[1;2]));
close(r.score,[0;0;100],1e-10);
% Первые два измерения ещё не различают третий вариант.
short=series; short.t=series.t(1:2); short.z=series.z(1:2);
short.lineage=series.lineage(1:2); short.covariance=series.covariance(1:2,1:2);
r=nav_terrain_sequence(paths(:,1:2),short,map,limits); assert(numel(r.candidates)==3);
series.z=[20 20 20]; r=nav_terrain_sequence(paths,series,map,limits);
assert(strcmp(r.status,'NO_MATCH')&&isempty(r.candidates));
map.coverage_complete=false;
fails(@()nav_terrain_sequence(paths,series,map,limits),'libekf:TerrainCoverage');
end
function requirement_watch_intervals
c=struct('windows',[1 3 6],'max_gap',1.1,'max_age',1,'capacity',20, ...
    'instant_limit',5,'mean_limit',[.9 .5 .3],'rms_limit',[10 10 10]);
h=[];
for t=0:6
    sample=struct('t',t,'context','independent_check_A','value',0);
    [h,r]=nav_watch(h,sample,t,c);
    if t==0, assert(strcmp(r.status,'NOT_READY')); end
end
assert(strcmp(r.status,'CONSISTENT')&&all(r.ready));
[~,r]=nav_watch(h,[],8,c); assert(strcmp(r.status,'NOT_READY')&&~any(r.ready));
sample.t=7; sample.value=6;
[~,r]=nav_watch(h,sample,7,c); assert(strcmp(r.status,'ALARM')&&r.instant_alarm);
[~,r]=nav_watch([],sample,7,c); assert(strcmp(r.status,'ALARM')&&~any(r.ready));
% Малое устойчивое расхождение выявляется длинной проверкой.
h=[];
for t=0:6
    sample.t=t; sample.value=.4; [h,r]=nav_watch(h,sample,t,c);
end
assert(strcmp(r.status,'ALARM')&&isequal(r.alarm,[false false true]));
sample.t=7; sample.context='changed_sources'; sample.value=0;
[h,r]=nav_watch(h,sample,7,c); assert(strcmp(r.status,'NOT_READY')&&numel(h.t)==1);
fails(@()nav_watch(h,sample,7,c),'libekf:Epoch');
sample.t=10; [h,r]=nav_watch(h,sample,10,c); assert(~all(r.ready));
% Ограниченная память не выдаётся за полноценное покрытие длинной истории.
c.capacity=3; h=[];
for t=0:10
    sample.t=t; [h,r]=nav_watch(h,sample,t,c);
end
assert(numel(h.t)<=3 && ~r.ready(3) && strcmp(r.status,'NOT_READY'));
end
function requirement_inertial_motion
[s,e]=fixture; s.q=[1;0;0;0]; s.ba=zeros(3,1); s.bg=zeros(3,1);
e.Omega=zeros(3,1); e.gravity_origin=zeros(3,1); e.gravity_gradient=zeros(3);
e.gravity_radius=1e5; e.max_step=.1; e.state_noise_density=zeros(21);
u=struct('source','IMU','lineage',{{'imu_interval_100'}},'quality','VALID', ...
    'frame','BODY','accel_unit','m/s^2','gyro_unit','rad/s','t_start',s.t, ...
    't_end',s.t+.1,'t_received',s.t+.11,'accel',[2;0;0],'gyro',zeros(3,1), ...
    'noise_density',zeros(6));
[x,~]=nav_predict(s,eye(21),u,e); dt=u.t_end-u.t_start;
close(x.p,s.p+s.v*dt+[dt^2;0;0],1e-11); close(x.v,s.v+[2*dt;0;0],1e-12);
close(x.clock,s.clock+s.drift*dt,1e-12); assert(x.t==u.t_end);
% Белый шум ускорения даёт известные интегралы погрешности положения и скорости.
u.accel=zeros(3,1); u.noise_density(1:3,1:3)=4*eye(3);
[~,Q]=nav_predict(s,zeros(21),u,e);
close(Q(1:3,1:3),4*dt^3/3*eye(3),1e-12);
close(Q(1:3,4:6),2*dt^2*eye(3),1e-12); close(Q(4:6,4:6),4*dt*eye(3),1e-12);
% Неподвижный на Земле прибор видит её вращение и противодействие тяжести.
s.v=zeros(3,1); e.Omega=[0;0;7.292115e-5]; e.gravity_origin=[0;0;-9.81];
u.gyro=e.Omega; u.accel=-e.gravity_origin;
[x,~]=nav_predict(s,eye(21),u,e);
close(x.p,s.p,1e-11); close(x.v,zeros(3,1),1e-12); close(nav_rotation(x.q),eye(3),1e-12);
u.quality='FAILED'; fails(@()nav_predict(s,eye(21),u,e),'libekf:Quality');
u.quality='VALID'; u.t_start=s.t-.01;
fails(@()nav_predict(s,eye(21),u,e),'libekf:Epoch');
% Полная последовательность через один обработчик, с общим временем и происхождением.
u.t_start=s.t; u.model='inertial';
filter=struct('state',s,'P',eye(21),'used',{{}},'capacity',20);
[a,r]=nav_assimilate(filter,u,[],e,100); assert(r.accepted && a.state.t==u.t_end);
sensor=struct('id',s.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
m=packet(a.state,'gnss_position'); m.frame='SENSOR_LOCAL'; m.lever=zeros(3,1);
m.z=a.state.p+[1;0;0]; m.lineage={'gnss_at_interval_end'};
[b,r]=nav_assimilate(a,m,sensor,e,100);
assert(r.accepted && b.state.p(1)>a.state.p(1) && b.state.p(1)<m.z(1));
assert(b.state.t==u.t_end && numel(b.used)==2);
[c,r]=nav_assimilate(b,u,[],e,100); assert(~r.accepted && isequal(b,c));
m.t_measure=s.t; m.lineage={'late_measurement'};
[c,r]=nav_assimilate(b,m,sensor,e,100); assert(~r.accepted && isequal(b,c));
% Производное решение, использовавшее этот инерциальный отсчёт, не является новым.
m.t_measure=b.state.t; m.lineage=u.lineage; m.model='external_position';
[c,r]=nav_assimilate(b,m,sensor,e,100); assert(~r.accepted && isequal(b,c));
end
function requirement_source_provenance
[s,e]=fixture;
sensor=struct('id',s.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
f=struct('state',s,'P',eye(21),'used',{{}},'capacity',10);
m=packet(s,'gnss_position'); m.frame='SENSOR_LOCAL'; m.lever=zeros(3,1); m.z=s.p+[1;2;3];
[a,r]=nav_assimilate(f,m,sensor,e,100); assert(r.accepted);
% Тот же исходный отсчёт под другим именем системы остаётся тем же свидетельством.
m.source='FMS'; m.model='external_position';
[b,r]=nav_assimilate(a,m,sensor,e,100); assert(~r.accepted && isequal(a,b));
m.lineage={'new_independent_record'}; m.quality='FAILED';
[b,r]=nav_assimilate(a,m,sensor,e,100); assert(~r.accepted && isequal(a,b));
m.quality='VALID'; m.noise_model='joint_required';
[b,r]=nav_assimilate(a,m,sensor,e,100); assert(~r.accepted && isequal(a,b));
m.noise_model='independent'; a.capacity=numel(a.used);
[b,r]=nav_assimilate(a,m,sensor,e,100); assert(~r.accepted && isequal(a,b));
% Здоровый другой источник всё же принимается, когда память доступна.
a.capacity=10;
[b,r]=nav_assimilate(a,m,sensor,e,100); assert(r.accepted && numel(b.used)==2);
end
function requirement_common_correction
[s,e]=fixture; s.q=[1;0;0;0];
m=packet(s,'gnss_position'); m.lever=zeros(3,1); m.z=s.p+[10;0;0]; m.R=9*eye(3);
o=nav_measurement(s,m,e); P=eye(21); P(1:3,1:3)=4*eye(3);
[x,Q,info]=nav_correct(s,P,o,100);
assert(info.accepted); close(x.p,s.p+[40/13;0;0],1e-12);
close(Q(1:3,1:3),(36/13)*eye(3),1e-12); close(Q(4:21,4:21),eye(18),1e-12);
% Второй независимый прибор той же точности: аналитическое объединение трёх оценок.
m.z=s.p+[4;0;0]; o2=nav_measurement(x,m,e);
[y,Q2]=nav_correct(x,Q,o2,100);
close(y.p,s.p+[56/17;0;0],1e-12); close(Q2(1:3,1:3),(36/17)*eye(3),1e-12);
% Связь положения со скоростью переносит поправку и на скорость.
P(1,4)=.5; P(4,1)=.5;
[x,Q]=nav_correct(s,P,o,100); close(x.v,s.v+[5/13;0;0],1e-12);
close(Q(4,4),1-.25/13,1e-12);
[x,Q,info]=nav_correct(s,P,o,.1);
assert(~info.accepted && isequal(x,s) && isequal(Q,P));
o.noise_model='joint_required'; fails(@()nav_correct(s,P,o,100),'libekf:Correlation');
o.noise_model='independent'; o.t_measure=s.t-1;
fails(@()nav_correct(s,P,o,100),'libekf:Epoch');
end
function requirement_source_models
% Independent geometric cases, manifold derivatives and invalid contracts.
[s,e]=fixture;
models={'gnss_position','external_position','gnss_velocity','external_velocity', ...
    'gnss_position_velocity','external_position_velocity','dme_range','pseudorange', ...
    'range_rate','vor_azimuth','baro_height','airdata','radar_range','ils_angles', ...
    'gls_deviation','attitude','radar_footprint'};
for k=1:numel(models)
    m=packet(s,models{k}); o=nav_measurement(s,m,e);
    assert(all(size(o.H)==[numel(o.residual),21]));
    assert(isequal(o.lineage,m.lineage) && strcmp(o.noise_model,m.noise_model));
    % Central differences of residual on the orientation manifold.
    steps=[repmat(.01,1,3),repmat(1e-4,1,3),repmat(1e-6,1,3), ...
        repmat(1e-5,1,3),repmat(1e-6,1,3),repmat(1e-4,1,3),.001,.001,.001];
    numeric=zeros(size(o.H));
    for j=1:21
        dx=zeros(21,1); dx(j)=steps(j);
        a=nav_measurement(nav_apply_error(s,dx),m,e);
        b=nav_measurement(nav_apply_error(s,-dx),m,e);
        numeric(:,j)=-(a.residual-b.residual)/(2*steps(j));
    end
    err=max(abs(numeric-o.H),[],'all');
    assert(err<3e-5, 'libekf:TestJacobian','%s Jacobian error %.6g',models{k},err);
end

% Exact geometric cases do not use the prediction as their expected value.
s.q=[1;0;0;0]; s.p=[1000;400;300];
s.v=[3;4;0]; s.wind=zeros(3,1); s.clock=7; s.drift=2;
m=packet(s,'dme_range'); m.lever=zeros(3,1); m.station=s.p-[3;4;0];
close(nav_measurement(s,m,e).predicted,5,1e-12);
m=packet(s,'pseudorange'); m.lever=zeros(3,1); m.station=s.p+[3;4;0]; m.range_correction=2;
close(nav_measurement(s,m,e).predicted,14,1e-12);
m=packet(s,'range_rate'); m.lever=zeros(3,1); m.station=s.p+[3;4;0];
m.station_velocity=s.v+[0;1;0]; m.rate_correction=.2;
close(nav_measurement(s,m,e).predicted,3,1e-12);
m=packet(s,'vor_azimuth'); m.lever=zeros(3,1); m.station=s.p-[0;4;3]; m.azimuth_offset=0;
close(nav_measurement(s,m,e).predicted,atan2(4,3),1e-12);
m=packet(s,'airdata');
close(nav_measurement(s,m,e).predicted,[5;0;atan2(4,3)],1e-12);
m=packet(s,'ils_angles'); m.lever=zeros(3,1); m.station=s.p-[3;4;12]; m.basis=eye(3);
close(nav_measurement(s,m,e).predicted,[atan2(4,3);atan2(12,5)],1e-12);
m=packet(s,'gls_deviation'); m.lever=zeros(3,1); m.path_origin=[1000;0;0];
close(nav_measurement(s,m,e).predicted,[400;300],1e-12);
s.p=[1000;0;0]; s.bh=8;
m=packet(s,'baro_height'); m.lever=zeros(3,1); m.geoid_offset=3;
close(nav_measurement(s,m,e).predicted,1005,1e-8);
m=packet(s,'radar_range'); m.lever=zeros(3,1); m.beam_body=[-1;0;0];
close(nav_measurement(s,m,e).predicted,1000,1e-10);

% Known 90-degree body-to-ECEF rotation and antenna kinematics.
s.q=[cos(pi/4);0;0;sin(pi/4)]; s.bg=zeros(3,1); e.Omega=zeros(3,1);
m=packet(s,'gnss_position'); m.lever=[1;0;0];
close(nav_measurement(s,m,e).predicted,s.p+[0;1;0],1e-9);
m=packet(s,'gnss_velocity'); m.lever=[1;0;0]; m.omega_ib=[0;0;2];
close(nav_measurement(s,m,e).predicted,s.v+[-2;0;0],1e-12);

% Quaternion sign is a representation choice, not a different attitude.
m=packet(s,'attitude'); a=nav_measurement(s,m,e); m.z=-m.z;
b=nav_measurement(s,m,e); close(a.residual,b.residual,1e-12);

% Azimuth residual must cross +/-pi continuously.
s.q=[1;0;0;0]; m=packet(s,'vor_azimuth'); m.lever=zeros(3,1);
m.station=s.p-[0;-1e-6;-1]; m.z=pi-1e-6; m.azimuth_offset=0;
close(nav_measurement(s,m,e).residual,-2e-6,1e-11);

% Contracts fail explicitly; no silent unit conversion or time reassignment.
m=packet(s,'gnss_position'); bad=m; bad.units={'km';'km';'km'};
fails(@()nav_measurement(s,bad,e),'libekf:Units');
bad=m; bad.t_measure=s.t-1; fails(@()nav_measurement(s,bad,e),'libekf:Epoch');
bad=m; bad.t_received=s.t-1; fails(@()nav_measurement(s,bad,e),'libekf:Time');
bad=m; bad.quality='INVALID'; fails(@()nav_measurement(s,bad,e),'libekf:Quality');
bad=m; bad.R(1,1)=-1; fails(@()nav_measurement(s,bad,e),'libekf:Covariance');
bad=m; bad.R(1,2)=1; fails(@()nav_measurement(s,bad,e),'libekf:Covariance');
bad=m; bad.lineage={}; fails(@()nav_measurement(s,bad,e),'libekf:Lineage');
bad=m; bad.z(1)=NaN; fails(@()nav_measurement(s,bad,e),'libekf:Dimension');
bad=m; bad.lever=[1;0;0]; bad.model='gnss_velocity'; bad=rmfield(bad,'omega_ib');
fails(@()nav_measurement(s,bad,e),'libekf:Lever');
bad=packet(s,'pseudorange'); bad.geometry_at_epoch=false;
fails(@()nav_measurement(s,bad,e),'libekf:GeometryEpoch');

% Pole and equator are analytic ellipsoidal cases.
[h,up]=nav_ellipsoid_height([0;0;e.ellipsoid_a*sqrt(1-e.ellipsoid_e2)+400],e.ellipsoid_a,e.ellipsoid_e2);
close(h,400,1e-8); close(up,[0;0;1],1e-12);
end

function [s,e]=fixture
e=struct('time_tolerance',1e-8,'ellipsoid_a',6378137,'ellipsoid_e2',6.69437999014e-3,'Omega',[0;0;7e-5]);
e.origin_ecef=[e.ellipsoid_a;0;0]; e.C_e_l=eye(3); e.frame_id='test_anchor_1';
s=struct('p',[1000;200;300],'v',[2;100;5], ...
    'q',[.9;.1;-.2;.3]/norm([.9;.1;-.2;.3]),'ba',[.01;-.02;.03], ...
    'bg',[.001;-.002;.003],'wind',[1;2;3],'bh',4,'clock',7,'drift',2,'t',100,'frame_id',e.frame_id);
end
function m=packet(s,model)
m=struct('source',['fixture_' model],'model',model,'z',zeros(3,1),'R',eye(3), ...
    'units',{{'m';'m';'m'}},'frame','LOCAL','frame_id',s.frame_id,'t_measure',s.t,'t_received',s.t+.01, ...
    'quality','VALID','lineage',{{['raw_' model '_100']}},'correlation_group',['group_' model], ...
    'noise_model','independent','lever',[.3;-.2;.1],'omega_ib',[.02;-.03;.04]);
switch model
    case {'gnss_position','external_position'}
        m.z=s.p;
    case {'gnss_velocity','external_velocity'}
        m.z=s.v; m.units=repmat({'m/s'},3,1);
    case {'gnss_position_velocity','external_position_velocity'}
        m.z=[s.p;s.v]; m.R=eye(6); m.units=[repmat({'m'},3,1);repmat({'m/s'},3,1)];
    case {'dme_range','pseudorange','range_rate'}
        m.station=s.p+[300;400;200]; m.station_velocity=[10;20;30];
        m.z=500; m.R=1; m.units={'m'}; m.geometry_at_epoch=true;
        m.range_correction=2; m.rate_correction=.2;
        if strcmp(model,'range_rate'), m.units={'m/s'}; end
    case 'vor_azimuth'
        m.station=s.p-[500;400;100]; m.basis=[0 0 1;1 0 0;0 1 0];
        m.azimuth_offset=.1; m.z=0; m.R=.01; m.units={'rad'};
    case 'baro_height'
        m.geoid_offset=20; m.z=1000; m.R=4; m.units={'m'};
    case 'airdata'
        m.z=[100;0;0]; m.R=diag([1 .001 .001]); m.units={'m/s';'rad';'rad'};
    case {'radar_range','radar_footprint'}
        m.vertices=[0 -10000 -10000;0 20000 -10000;0 -10000 20000]; m.faces=[1 2 3];
        m.edge_margin=.01;
        m.half_angle=pi/3; m.max_valid_range=2000; m.nearest_return_valid=true;
        m.coverage_complete=true; m.terrain_quality='VALID';
        m.beam_body=nav_rotation(s.q).'*[-1;0;0]; m.min_incidence=.2;
        m.z=1000; m.R=1; m.units={'m'};
    case 'ils_angles'
        m.station=s.p-[1000;20;30]; m.basis=eye(3); m.angle_reference=[0;0];
        m.z=[0;0]; m.R=eye(2)*.001; m.units={'rad';'rad'};
    case 'gls_deviation'
        m.path_origin=[0;0;0]; m.projection=[0 1 0;0 0 1];
        m.z=[0;0]; m.R=eye(2); m.units={'m';'m'};
    case 'attitude'
        dx=zeros(21,1); dx(7:9)=[.03;-.02;.01]; measured=nav_apply_error(s,dx);
        m.z=measured.q; m.R=diag([.001 .002 .003]); m.units={'rad';'rad';'rad'}; m.max_angle=.5;
end
end
function close(a,b,tol)
assert(max(abs(a(:)-b(:)))<=tol,'libekf:TestValue','Analytic value mismatch.');
end
function fails(f,id)
try, f(); catch ex, assert(strcmp(ex.identifier,id),'libekf:TestError','Expected %s, got %s',id,ex.identifier); return; end
error('libekf:TestError','Expected failure %s.',id);
end


function requirement_coordinates_and_terrain
env=struct('origin_ecef',[6378137;0;0],'C_e_l',eye(3),'frame_id','common', ...
    'ellipsoid_a',6378137,'ellipsoid_e2',6.69437999014e-3,'time_tolerance',1e-8,'Omega',[0;0;0]);
s=struct('p',[13;24;35],'v',[1;2;3],'q',[1;0;0;0],'ba',zeros(3,1),'bg',zeros(3,1), ...
    'wind',zeros(3,1),'bh',0,'clock',0,'drift',0,'t',100,'frame_id','common');
sensor=struct('id','sensor_A','origin_ecef',env.origin_ecef+[10;20;30], ...
    'C_e_s',[0 -1 0;1 0 0;0 0 1]);
m=struct('source','position_A','model','gnss_position','z',[4;-3;5],'R',diag([1 4 9]), ...
    'units',{{'m';'m';'m'}},'frame','SENSOR_LOCAL','frame_id',sensor.id,'t_measure',100,'t_received',100, ...
    'quality','VALID','lineage',{{'A_100'}},'correlation_group','A','noise_model','independent');
o=nav_observe_source(s,m,sensor,env);
assert(norm(o.residual)<1e-12 && norm(o.predicted-m.z)<1e-12);
assert(norm(o.H(:,1:3)-sensor.C_e_s.','fro')<1e-12 && norm(o.R-m.R,'fro')<1e-12);
assert(isequal(m.z,[4;-3;5]) && strcmp(m.frame,'SENSOR_LOCAL')); % Original untouched.

% Geodetic equator and a nonlocal point use ellipsoid curvature, no flat Earth.
s.p=[1000;0;0]; output=nav_geodetic_output(s,eye(21),env);
assert(abs(output.lat)<1e-12 && abs(output.lon)<1e-12 && abs(output.alt-1000)<1e-8);
assert(strcmp(output.quality.altitude,'UNASSESSED'));
lat=pi/4; lon=pi/2; height=123;
N=env.ellipsoid_a/sqrt(1-env.ellipsoid_e2*sin(lat)^2);
expected=[(N+height)*cos(lat)*cos(lon);(N+height)*cos(lat)*sin(lon); ...
    (N*(1-env.ellipsoid_e2)+height)*sin(lat)];
s.p=expected-env.origin_ecef;
output=nav_geodetic_output(s,eye(21),env);
assert(abs(output.lat-45)<1e-10 && abs(output.lon-90)<1e-10 && abs(output.alt-height)<1e-7);
% Pole is represented explicitly; no divide by vanishing longitude radius.
s.p=[0;0;env.ellipsoid_a*sqrt(1-env.ellipsoid_e2)+100]-env.origin_ecef;
output=nav_geodetic_output(s,eye(21),env);
assert(output.lat==90 && isnan(output.lon) && abs(output.alt-100)<1e-7);
assert(strcmp(output.quality.longitude,'UNDEFINED'));

% Ravine: two slopes, ten metres below the ridge at the centre.
vertices=[0 -10 -10;-10 0 -10;-10 0 10;0 -10 10;0 10 -10;0 10 10];
faces=[1 2 3;1 3 4;2 5 6;2 6 3];
hit=nav_terrain_footprint([2;2;1],[-1;0;0],pi/3,vertices,faces,.01);
assert(norm(hit.point-[-3;7;1])<1e-12 && abs(hit.range-sqrt(50))<1e-12);
s.p=[2;2;1]; m.model='radar_range'; m.frame='LOCAL'; m.frame_id=env.frame_id;
m.z=sqrt(50); m.R=1; m.units={'m'}; m.beam_body=[-1;0;0]; m.nearest_return_valid=true;
m.vertices=vertices; m.faces=faces; m.edge_margin=.01; m.half_angle=pi/3; m.max_valid_range=2000;
m.coverage_complete=true; m.terrain_quality='VALID';
o=nav_measurement(s,m,env); assert(abs(o.predicted-sqrt(50))<1e-12);
% Equal nearest returns from both slopes must retain their ambiguity.
ambiguous=nav_terrain_footprint([2;0;1],[-1;0;0],pi/3,vertices,faces,.01);
assert(~ambiguous.smooth);
fails(@()nav_terrain_footprint([10;100;1],[-1;0;0],pi/3,vertices,faces,.01),'libekf:TerrainCoverage');
m.nearest_return_valid=false; fails(@()nav_measurement(s,m,env),'libekf:RadarFootprint');

% Proper quaternion conversion remains valid at a half turn.
for axis=1:3
    q=zeros(4,1); q(axis+1)=1; C=nav_rotation(q);
    assert(norm(nav_rotation(nav_quaternion(C))-C,'fro')<1e-12);
end
end

function requirement_numerical_invariance
from=struct('origin_ecef',[6378137;0;0],'C_e_l',eye(3),'frame_id','F');
to=struct('origin_ecef',from.origin_ecef+[100;200;300],'C_e_l',[0 -1 0;1 0 0;0 0 1],'frame_id','S');
s=struct('p',[103;204;305],'v',[1;2;3],'q',[1;0;0;0],'ba',zeros(3,1),'bg',zeros(3,1), ...
    'wind',[4;5;6],'bh',0,'clock',0,'drift',0,'t',100,'frame_id','F');
P=diag(1:21); P(1,4)=.1; P(4,1)=.1;
[changed,Q]=nav_change_frame(s,P,from,to);
assert(norm(changed.p-[4;-3;5])<1e-12);
[restored,Pback]=nav_change_frame(changed,Q,to,from);
assert(norm(restored.p-s.p)<1e-12 && norm(restored.v-s.v)<1e-12 && norm(Pback-P,'fro')<1e-12);
assert(norm(nav_rotation(restored.q)-eye(3),'fro')<1e-12 && restored.t==s.t);
% Cross-covariances must transform too, not only three diagonal variances.
assert(abs(Q(2,5)-.1)<1e-12);

obs=struct('H',[1 0;0 1],'R',[2 .5;.5 3],'residual',[1;2]);
prior=[4 1;1 5]; n=nav_innovation(prior,obs);
% Direct analytic inverse of [6,1.5;1.5,8] gives this quadratic value.
expected=(8-6+24)/(48-2.25);
assert(abs(n.nis-expected)<1e-12);
% Change measurement units over 200 orders; normalized result must agree.
D=diag([1e-100,1e100]); scaled=obs;
scaled.H=D*obs.H; scaled.R=D*obs.R*D; scaled.residual=D*obs.residual;
n2=nav_innovation(prior,scaled);
assert(abs(n2.nis-n.nis)<1e-12 && norm(n2.whitened-n.whitened)<1e-12);
% Very small but representable covariance is not a reason for an inverse.
tiny=struct('H',1,'R',1e-300,'residual',sqrt(2e-300));
n3=nav_innovation(1e-300,tiny); assert(abs(n3.nis-1)<1e-12);
end


function requirement_broad_radar
% A slope inside a broad beam can be closer than the nadir intersection.
o=[10;0;0]; b=[-1;0;0]; vertices=[-100 -100 -100;100 100 -100;0 0 200]; faces=[1 2 3];
hit=nav_terrain_footprint(o,b,pi/3,vertices,faces,.01);
assert(abs(hit.range-sqrt(50))<1e-11 && norm(hit.point-[5;5;0])<1e-10 && hit.smooth);
assert(~hit.cone_active);
small=pi/18; constrained=nav_terrain_footprint(o,b,small,vertices,faces,.01);
expected=10/(cos(small)+sin(small));
assert(abs(constrained.range-expected)<1e-11 && constrained.cone_active && constrained.smooth);

env=struct('frame_id','local','time_tolerance',1e-8);
s=struct('p',o,'v',zeros(3,1),'q',[1;0;0;0],'ba',zeros(3,1),'bg',zeros(3,1), ...
    'wind',zeros(3,1),'bh',0,'clock',0,'drift',0,'t',10,'frame_id','local');
m=struct('source','radar','model','radar_footprint','z',expected,'R',1,'units',{{'m'}}, ...
    'frame','LOCAL','frame_id','local','t_measure',10,'t_received',10.01,'quality','VALID', ...
    'lineage',{{'radar_10'}},'correlation_group','radar_map','noise_model','joint_required', ...
    'beam_body',b,'half_angle',small,'vertices',vertices,'faces',faces,'edge_margin',.01, ...
    'max_valid_range',2000,'nearest_return_valid',true,'coverage_complete',true,'terrain_quality','VALID');
obs=nav_measurement(s,m,env); assert(abs(obs.residual)<1e-11);
% Derivative of the continuous cone-boundary optimum, not sampled rays.
for j=[1:3,7:9]
    dx=zeros(21,1); dx(j)=1e-6;
    a=nav_measurement(nav_apply_error(s,dx),m,env); c=nav_measurement(nav_apply_error(s,-dx),m,env);
    numeric=-(a.residual-c.residual)/2e-6;
    assert(abs(numeric-obs.H(j))<1e-7,'libekf:TestJacobian','Footprint Jacobian column %d.',j);
end
bad=m; bad.z=2001; fails(@()nav_measurement(s,bad,env),'libekf:RadarRange');
bad=m; bad.coverage_complete=false; fails(@()nav_measurement(s,bad,env),'libekf:TerrainCoverage');
bad=m; bad.nearest_return_valid=false; fails(@()nav_measurement(s,bad,env),'libekf:RadarFootprint');

% A footprint can see a closer elevated patch away from the beam axis.
v=[vertices; 4 2 -1;4 4 -1;4 3 2]; f=[1 2 3;4 5 6];
closer=nav_terrain_footprint(o,b,pi/3,v,f,.01);
assert(closer.range<hit.range && closer.point(2)>0);
end


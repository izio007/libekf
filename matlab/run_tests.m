function status=run_tests
% Проверки требований к поведению; один файл, без внешних настроек и отчётов.
base=fileparts(mfilename('fullpath')); addpath(fullfile(base,'function'));
cases={@requirement_source_models,@requirement_coordinates_and_terrain, ...
    @requirement_numerical_invariance,@requirement_broad_radar,@requirement_common_correction, ...
    @requirement_source_provenance,@requirement_inertial_motion,@requirement_watch_intervals, ...
    @requirement_terrain_sequence,@requirement_all_source_fusion,@requirement_relative_camera_motion,@requirement_joint_correction,@requirement_visual_ambiguity,@requirement_correlated_measurement,@requirement_source_comparison,@requirement_navigation_outage,@requirement_interval_observations,@requirement_observable_quantities,@requirement_all_source_timing};
names={'Показания соответствуют физическим моделям источников', ...
    'Смена координат сохраняет результат; неоднозначность рельефа не скрывается', ...
    'Мера расхождения сохраняется при смене единиц и масштаба', ...
    'Радиовысотомер ищет ближайшую поверхность во всей диаграмме', ...
    'Общая поправка учитывает погрешности и сохраняет состояние при отклонении', ...
    'Повторные и неисправные источники не изменяют общую оценку', ...
    'Инерциальные данные продолжают движение и накапливают погрешность', ...
    'Три интервала контроля различают тревогу, согласие и нехватку истории', ...
    'Привязка использует последовательность и сохраняет неразличимые варианты', ...
    'Все модели уточняют общую оценку независимо от осей прибора', ...
    'Перемещение камеры учитывает два момента и её размещение на корпусе', ...
    'Совместная поправка сохраняет неопределённость общего положения', ...
    'Альтернативное объяснение кадров сохраняет неоднозначность переноса', ...
    'Зависимость ошибки измерения от оценки учитывается явно', ...
    'Сравнение источников учитывает общие ошибки и передаёт расхождение сторожу', ...
    'Общая оценка продолжается без спутников и отклоняет ложную координату', ...
    'Измерения внутри инерциального интервала учитываются в своё время', ...
    'Исправные источники определяют только наблюдаемые величины', ...
    'Все источники учитывают своё время при пересчёте истории'};
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
models=source_models;
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
% Two-epoch errors remain correlated: the old position is not exact.
[~,Q,T,N]=nav_predict(s,eye(21),u,e);
close(T(1:3,4:6),dt*eye(3),1e-12);
close(T(20,21),dt,1e-12);
close(N(1:3,1:3),4*dt^3/3*eye(3),1e-12);
cross_epoch=T.'; joint=[eye(21) cross_epoch;cross_epoch.' Q];
% Position increment error cancels common initial position error.
D=zeros(3,42); D(:,1:3)=-eye(3); D(:,22:24)=eye(3);
% Initial velocity and accelerometer bias uncertainty also enter the increment.
expected=(dt^2+dt^4/4+4*dt^3/3)*eye(3);
close(D*joint*D.',expected,1e-12);
% Omitting the cross covariance would add two spurious position variances.
close(D*blkdiag(eye(21),Q)*D.'-expected,2*eye(3),1e-12);

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
initial=struct('state',s,'P',eye(21),'used',{{}},'capacity',100);
cases={ ...
    'dme_range',{'station'}; ...
    'pseudorange',{'station','geometry_at_epoch','range_correction'}; ...
    'range_rate',{'station','station_velocity','rate_correction','geometry_at_epoch'}; ...
    'baro_height',{'geoid_offset'}; ...
    'ils_angles',{'basis','station','angle_reference'}; ...
    'gls_deviation',{'projection','path_origin'}; ...
    'radar_footprint',{'beam_body','coverage_complete','terrain_quality','nearest_return_valid', ...
        'vertices','faces','half_angle','edge_margin','max_valid_range'}; ...
    'air_pressures',{'geoid_offset','atmosphere','pressure_valid','calibration_valid','gamma','max_mach','flow_valid'}; ...
    'optical_landmarks',{'calibration_valid','matches_valid','map_valid','camera_to_body','intrinsics','landmarks','min_depth'}; ...
    'magnetic_vector',{'field_valid','calibration_valid','field_gradient','field_radius','field_time','field_origin','field_vector'}; ...
    'rsdn_time_difference',{'geometry_at_epoch','stations','propagation_speed','time_correction','reference_station'}; ...
    'doppler_velocity',{'beams_body','surface_valid','beam_geometry_valid','surface_velocity'}; ...
    'adf_bearing',{'bearing_geometry_valid','station','bearing_offset'}; ...
    'station_range_azimuth',{'basis','station','azimuth_offset'}; ...
    'time_difference',{'geometry_at_epoch','station','reference_station','range_correction','propagation_speed'}};
for c=1:size(cases,1)
    complete=packet(s,cases{c,1}); complete.frame='SENSOR_LOCAL';
    for j=1:numel(cases{c,2})
        incomplete=rmfield(complete,cases{c,2}{j});
        [unchanged,result]=nav_assimilate(initial,incomplete,sensor,e,100);
        assert(~result.accepted && strcmp(result.reason,'libekf:Packet'));
        assert(isequal(unchanged,initial),'libekf:TestAtomic','Incomplete input changed accepted state.');
    end
end

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
models=source_models;
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
    if strcmp(models{k},'magnetic_vector')
        % Tesla-scale derivatives require a physical scale, not a metre-sized tolerance.
        assert(err<1e-12,'libekf:TestJacobian','Magnetic derivative mismatch.');
    end
    assert(err<3e-5, 'libekf:TestJacobian','%s Jacobian error %.6g',models{k},err);
end

% Magnetic vector has an unobservable rotation around a uniform field.
ms=s; ms.q=[1;0;0;0]; m=packet(ms,'magnetic_vector');
m.lever=zeros(3,1); m.field_gradient=zeros(3); m.field_vector=[20;0;40]*1e-6;
o=nav_measurement(ms,m,e); close(o.predicted,[20;0;40]*1e-6,1e-15);
close(o.H(:,1:3),zeros(3),1e-15);
close(o.H(:,7:9)*m.field_vector,zeros(3,1),1e-20);
ms.q=[cos(pi/4);0;0;sin(pi/4)];
close(nav_measurement(ms,m,e).predicted,[0;-20;40]*1e-6,1e-15);
m.field_valid=false; fails(@()nav_measurement(ms,m,e),'libekf:Magnetic');
m.field_valid=true; m.calibration_valid=false;
fails(@()nav_measurement(ms,m,e),'libekf:Magnetic');
m.calibration_valid=true; m.field_time=ms.t-1;
fails(@()nav_measurement(ms,m,e),'libekf:Epoch');
m.field_time=ms.t; m.field_radius=1; ms.p=ms.p+[2;0;0];
fails(@()nav_measurement(ms,m,e),'libekf:Magnetic');

% A uniform axial field gives two angular constraints, no position information.
ms=s; ms.q=[1;0;0;0]; m=packet(ms,'magnetic_vector');
m.lever=zeros(3,1); m.field_gradient=zeros(3); m.field_vector=[0;0;40e-6];
m.z=m.field_vector; m.R=eye(3)*1e-12;
prior=eye(21); prior(7:9,7:9)=eye(3)*.01;
[candidate,Q,info]=nav_correct(ms,prior,nav_measurement(ms,m,e),100);
assert(info.accepted); close(candidate.p,ms.p,1e-12);
close(nav_rotation(candidate.q),eye(3),1e-12);
expected=prior; variance=1/(1/.01+(40e-6)^2/1e-12);
expected(7,7)=variance; expected(8,8)=variance;
close(Q,expected,1e-12);
% An independently specified linear field changes by G*displacement.
m.field_gradient=diag([2 -1 -1])*1e-9;
ms.p=ms.p+[3;4;0];
close(nav_measurement(ms,m,e).predicted,[6e-9;-4e-9;40e-6],1e-15);

% Independent pinhole geometry and simultaneous scale ambiguity.
cs=s; cs.q=[1;0;0;0]; cs.p=zeros(3,1); m=packet(cs,'optical_landmarks');
m.lever=zeros(3,1); m.landmarks=[1 2 10;-2 1 20];
o=nav_measurement(cs,m,e); close(o.predicted,[400;400;240;280],1e-12);
m.landmarks=2*m.landmarks;
close(nav_measurement(cs,m,e).predicted,o.predicted,1e-12);
m.landmarks=[1 2 -10;-2 1 20];
fails(@()nav_measurement(cs,m,e),'libekf:Optical');
m.landmarks=[1 2 10;-2 1 20]; m.matches_valid=false;
fails(@()nav_measurement(cs,m,e),'libekf:Optical');
m.matches_valid=true; m.calibration_valid=false;
fails(@()nav_measurement(cs,m,e),'libekf:Optical');
m.calibration_valid=true; m.map_valid=false;
fails(@()nav_measurement(cs,m,e),'libekf:Optical');

% Prepared distant directions constrain attitude, not position directly.
cs=s; cs.q=[1;0;0;0]; m=packet(cs,'celestial_pixels');
o=nav_measurement(cs,m,e); close(o.predicted,[400;400;160;320],1e-12);
close(o.H(:,1:6),zeros(4,6),1e-12);
cs.p=cs.p+[10000;-2000;500];
close(nav_measurement(cs,m,e).predicted,o.predicted,1e-12);
m.geometry_at_epoch=false; fails(@()nav_measurement(cs,m,e),'libekf:Optical');
m.geometry_at_epoch=true; m.catalog_valid=false;
fails(@()nav_measurement(cs,m,e),'libekf:Optical');

% Peer position is uncertain, not a perfect beacon.
cs=s; m=packet(cs,'cooperative_range'); m.lever=zeros(3,1); m.station=cs.p+[3;4;0];
m.peer_covariance=diag([4 9 1]); m.peer_cross_cov(1,1)=.5; m.noise_model='joint_required';
o=nav_measurement(cs,m,e); close(o.predicted,5,1e-12); close(o.R,8.2,1e-12);
close(o.state_noise_cov(1),.3,1e-12);
m.peer_valid=false; fails(@()nav_measurement(cs,m,e),'libekf:Peer');
m.peer_valid=true; m.peer_time=cs.t-1; fails(@()nav_measurement(cs,m,e),'libekf:Peer');

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

% Разные радиосредства сохраняют физический смысл угла и разности дальностей.
s.q=[1;0;0;0];
m=packet(s,'adf_bearing'); m.lever=zeros(3,1); m.station=s.p+[3;4;12];
close(nav_measurement(s,m,e).predicted,atan2(4,3),1e-12);
s.q=[cos(pi/4);0;0;sin(pi/4)];
close(nav_measurement(s,m,e).predicted,atan2(-3,4),1e-12);
s.q=[1;0;0;0];
m=packet(s,'station_range_azimuth'); m.lever=zeros(3,1);
m.station=s.p-[0;3;4]; m.azimuth_offset=0;
close(nav_measurement(s,m,e).predicted,[5;atan2(3,4)],1e-12);
m=packet(s,'range_difference'); m.lever=zeros(3,1);
m.station=s.p+[3;4;0]; m.reference_station=s.p+[5;12;0];
o=nav_measurement(s,m,e); close(o.predicted,-6,1e-12); assert(o.H(20)==0);
s.clock=s.clock+100; close(nav_measurement(s,m,e).predicted,-6,1e-12);
m.model='time_difference'; m.units={'s'}; m.R=1/m.propagation_speed^2;
close(nav_measurement(s,m,e).predicted,-6/299792458,1e-18);
m.geometry_at_epoch=false; fails(@()nav_measurement(s,m,e),'libekf:GeometryEpoch');

% РСДН: пары с общей опорной станцией обрабатываются совместно.
m=packet(s,'rsdn_time_difference'); m.lever=zeros(3,1);
m.stations=s.p.'+[3 4 0;6 8 0]; m.reference_station=s.p+[5;12;0];
o=nav_measurement(s,m,e);
close(o.predicted,[-6;-2]/299792458,1e-18);
assert(o.R(1,2)>0 && all(o.H(:,20)==0));

% Static pressure remains in pascals and uses the declared atmosphere, not a hidden standard day.
[ps,pe]=fixture; ps.q=[1;0;0;0]; ps.bh=0; ps.p=zeros(3,1);
pm=packet(ps,'static_pressure'); pm.lever=zeros(3,1); pm.geoid_offset=0;
pm.atmosphere.gradient=0;
close(nav_measurement(ps,pm,pe).predicted,100000,1e-7);
ps.p=[9000*log(2);0;0];
po=nav_measurement(ps,pm,pe); close(po.predicted,50000,1e-7);
close(po.H(1),-50000/9000,1e-9); close(po.H(19),po.H(1),1e-12);
ps.p=[1000;0;0]; pm.atmosphere.gradient=-.006;
close(nav_measurement(ps,pm,pe).predicted,100000*.98^(50/9),1e-7);
pm.atmosphere.gradient=1e-12;
close(nav_measurement(ps,pm,pe).predicted,100000*exp(-1/9),1e-6);
pm.pressure_valid=false; fails(@()nav_measurement(ps,pm,pe),'libekf:Atmosphere');
pm.pressure_valid=true; pm.atmosphere.time_end=ps.t-1;
fails(@()nav_measurement(ps,pm,pe),'libekf:Atmosphere');
pm.atmosphere.time_end=ps.t; pm.atmosphere.max_height=999;
fails(@()nav_measurement(ps,pm,pe),'libekf:Atmosphere');

% Impact pressure depends on air-relative speed and the compressible pressure ratio.
[ps,pe]=fixture; ps.p=zeros(3,1); ps.bh=0; ps.q=[1;0;0;0];
ps.wind=[10;20;30]; ps.v=ps.wind+[.5*sqrt(126000);0;0];
pm=packet(ps,'impact_pressure'); pm.lever=zeros(3,1); pm.geoid_offset=0;
pm.atmosphere.gradient=0;
po=nav_measurement(ps,pm,pe); close(po.predicted,100000*(1.05^3.5-1),1e-7);
close(po.H(4:6),-po.H(16:18),1e-12);
ps.v=ps.v+[40;50;60]; ps.wind=ps.wind+[40;50;60];
close(nav_measurement(ps,pm,pe).predicted,po.predicted,1e-7);
ps.v=ps.wind; po=nav_measurement(ps,pm,pe);
close(po.predicted,0,0); close(po.H,zeros(1,21),0);
ps.wind=zeros(3,1); ps.v=[1e-5;0;0];
close(nav_measurement(ps,pm,pe).predicted,.5*(100000/90000)*1e-10,1e-20);
ps.v=ps.wind+[sqrt(126000);0;0];
fails(@()nav_measurement(ps,pm,pe),'libekf:Airflow');
ps.v=ps.wind; pm.flow_valid=false;
fails(@()nav_measurement(ps,pm,pe),'libekf:Airflow');

% Joint pressures and their equivalent static/difference coordinates carry identical information.
[ps,pe]=fixture; ps.p=zeros(3,1); ps.bh=0; ps.q=[1;0;0;0];
ps.wind=zeros(3,1); ps.v=[.5*sqrt(126000);0;0];
pm=packet(ps,'air_pressures'); pm.lever=zeros(3,1); pm.geoid_offset=0;
pm.atmosphere.gradient=0; pm.R=[400 100;100 900];
pm.z=[100010;100000*1.05^3.5+20];
o=nav_measurement(ps,pm,pe);
close(o.predicted,[100000;100000*1.05^3.5],1e-7); close(o.R,pm.R,0);
T=[1 0;-1 1]; transformed=o;
transformed.H=T*o.H; transformed.residual=T*o.residual; transformed.R=T*o.R*T.';
[x,Q,info]=nav_correct(ps,eye(21),o,100); assert(info.accepted);
[y,V,info]=nav_correct(ps,eye(21),transformed,100); assert(info.accepted);
close(x.p,y.p,1e-10); close(x.v,y.v,1e-10); close(Q,V,1e-10);
assert(transformed.R(1,2)==-300); % Shared static pressure error cannot be discarded.
sensor=struct('id',ps.frame_id,'origin_ecef',pe.origin_ecef,'C_e_s',pe.C_e_l);
f=struct('state',ps,'P',eye(21),'used',{{}},'capacity',10);
pm.frame='SENSOR_LOCAL'; pm.lineage={'static_raw','total_raw'};
[f,info]=nav_assimilate(f,pm,sensor,pe,100); assert(info.accepted&&numel(f.used)==2);
pm.model='impact_pressure'; pm.z=diff(pm.z); pm.R=1100; pm.units={'Pa'};
[unchanged,info]=nav_assimilate(f,pm,sensor,pe,100);
assert(~info.accepted&&isequal(unchanged,f));

% Доплеровские проекции относятся к поверхности, а не к воздушному потоку.
s.q=[1;0;0;0]; s.v=[3;4;5];
m=packet(s,'doppler_velocity'); m.lever=zeros(3,1);
close(nav_measurement(s,m,e).predicted,[-3;-1;-7]/sqrt(2),1e-12);
s.wind=[100;200;300];
close(nav_measurement(s,m,e).predicted,[-3;-1;-7]/sqrt(2),1e-12);
s.q=[cos(pi/4);0;0;sin(pi/4)];
close(nav_measurement(s,m,e).predicted,[-1;-7;-9]/sqrt(2),1e-12);
bad=m; bad.surface_valid=false;
fails(@()nav_measurement(s,bad,e),'libekf:Doppler');
s.v=[3;4;0]; s.wind=zeros(3,1);

% Known 90-degree body-to-ECEF rotation and antenna kinematics.
s.q=[cos(pi/4);0;0;sin(pi/4)]; s.bg=zeros(3,1); e.Omega=zeros(3,1);
bad=packet(s,'ils_angles'); bad.angle_reference=0;
fails(@()nav_measurement(s,bad,e),'libekf:Geometry');
bad=packet(s,'baro_height'); bad.geoid_offset=[0;0];
fails(@()nav_measurement(s,bad,e),'libekf:Height');
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
    case 'cooperative_range'
        m.station=s.p+[300;400;200]; m.peer_covariance=eye(3)*.04;
        m.peer_cross_cov=zeros(21,3); m.peer_cross_frame_id=s.frame_id;
        m.peer_time=s.t; m.peer_valid=true; m.z=500; m.R=1; m.units={'m'};
    case {'dme_range','pseudorange','range_rate'}
        m.station=s.p+[300;400;200]; m.station_velocity=[10;20;30];
        m.z=500; m.R=1; m.units={'m'}; m.geometry_at_epoch=true;
        m.range_correction=2; m.rate_correction=.2;
        if strcmp(model,'range_rate'), m.units={'m/s'}; end
    case {'range_difference','time_difference'}
        m.station=s.p+[300;400;200]; m.reference_station=s.p+[-500;100;300];
        m.range_correction=2; m.geometry_at_epoch=true; m.propagation_speed=299792458;
        m.z=0; m.R=1; m.units={'m'};
        if strcmp(model,'time_difference'), m.units={'s'}; m.R=1/m.propagation_speed^2; end
    case 'rsdn_time_difference'
        m.stations=s.p.'+[300 400 200;-400 100 200];
        m.reference_station=s.p+[-500;100;300]; m.propagation_speed=299792458;
        m.time_correction=[2;1]/m.propagation_speed; m.geometry_at_epoch=true;
        m.z=[0;0]; m.R=[1 .5;.5 1]/m.propagation_speed^2; m.units={'s';'s'};
    case 'celestial_pixels'
        m.camera_to_body=eye(3); m.intrinsics=[800;800;320;240];
        directions=[.1 -.2;.2 .1;1 1]; directions=directions./vecnorm(directions);
        m.directions=(nav_rotation(s.q)*directions).'; m.min_forward=.1;
        m.calibration_valid=true; m.matches_valid=true; m.catalog_valid=true; m.geometry_at_epoch=true;
        m.z=[400;400;160;320]; m.R=eye(4); m.units=repmat({'pixel'},4,1);
    case 'optical_landmarks'
        m.camera_to_body=eye(3); m.intrinsics=[800;800;320;240];
        m.landmarks=(s.p+nav_rotation(s.q)*(m.lever+[1 -2;2 1;100 120])).';
        m.min_depth=.1; m.calibration_valid=true; m.matches_valid=true; m.map_valid=true;
        m.z=[328;256;320-40/3;240+20/3]; m.R=eye(4); m.units=repmat({'pixel'},4,1);
    case 'magnetic_vector'
        m.field_vector=[20;5;40]*1e-6; m.field_origin=s.p;
        m.field_gradient=diag([1 -2 1])*1e-10; m.field_radius=10000;
        m.field_time=s.t; m.field_valid=true; m.calibration_valid=true;
        m.z=[20;5;40]*1e-6; m.R=eye(3)*1e-12; m.units=repmat({'T'},3,1);
    case 'adf_bearing'
        m.station=s.p+[300;400;200]; m.bearing_offset=0; m.bearing_geometry_valid=true;
        m.z=0; m.R=.01; m.units={'rad'};
    case 'station_range_azimuth'
        m.station=s.p-[500;400;100]; m.basis=[0 0 1;1 0 0;0 1 0];
        m.azimuth_offset=.1; m.z=[600;0]; m.R=diag([1 .01]); m.units={'m';'rad'};
    case 'vor_azimuth'
        m.station=s.p-[500;400;100]; m.basis=[0 0 1;1 0 0;0 1 0];
        m.azimuth_offset=.1; m.z=0; m.R=.01; m.units={'rad'};
    case {'static_pressure','impact_pressure','air_pressures'}
        m.geoid_offset=20; m.z=100000; m.R=400; m.units={'Pa'};
        m.pressure_valid=true; m.calibration_valid=true;
        m.atmosphere=struct('height',0,'pressure',100000,'temperature',300,'gradient',-.006, ...
            'gravity',10,'gas_constant',300,'min_height',-1000,'max_height',10000, ...
            'time_start',s.t-10,'time_end',s.t+10,'valid',true);
        if any(strcmp(model,{'impact_pressure','air_pressures'}))
            m.gamma=1.4; m.max_mach=.95; m.flow_valid=true; m.z=100;
            if strcmp(model,'air_pressures'), m.z=[100000;100100]; m.R=[400 100;100 900]; m.units={'Pa';'Pa'}; end
        end
    case 'baro_height'
        m.geoid_offset=20; m.z=1000; m.R=4; m.units={'m'};
    case 'doppler_velocity'
        m.beams_body=[1 0 -1;0 1 -1;-1 0 -1]/sqrt(2);
        m.surface_velocity=[1;0;0]; m.surface_valid=true; m.beam_geometry_valid=true;
        m.z=[-3;-1;-7]/sqrt(2); m.R=eye(3); m.units=repmat({'m/s'},3,1);
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

% Output checks apply only to the quantities they actually monitor.
pole=s; s.p=[1000;0;0];
c=struct('windows',[1 3 6],'max_gap',1.1,'max_age',1,'capacity',20, ...
    'instant_limit',5,'mean_limit',[1 1 1],'rms_limit',[2 2 2]);
h=[];
for t=s.t-6:s.t
    [h,good]=nav_watch(h,struct('t',t,'context','height_check','value',0),t,c);
end
control=struct('state_time',s.t,'frame_id',s.frame_id,'now',s.t, ...
    'checks',struct('altitude',{{good}}));
out=nav_geodetic_output(s,eye(21),env,control);
assert(strcmp(out.quality.altitude,'CONSISTENT')&&strcmp(out.quality.latitude,'NOT_READY'));
[~,bad]=nav_watch([],struct('t',s.t,'context','other_height_check','value',6),s.t,c);
control.checks.altitude={good,bad};
out=nav_geodetic_output(s,eye(21),env,control); assert(strcmp(out.quality.altitude,'ALARM'));
[~,unready]=nav_watch([],[],s.t,c); control.checks.altitude={good,unready};
out=nav_geodetic_output(s,eye(21),env,control); assert(strcmp(out.quality.altitude,'NOT_READY'));
control.now=s.t+2;
fails(@()nav_geodetic_output(s,eye(21),env,control),'libekf:WatchOutput');
[~,stale]=nav_watch(h,[],control.now,c); control.checks.altitude={stale};
out=nav_geodetic_output(s,eye(21),env,control); assert(strcmp(out.quality.altitude,'NOT_READY'));
control.state_time=s.t-1;
fails(@()nav_geodetic_output(s,eye(21),env,control),'libekf:WatchOutput');
control.state_time=s.t; control.now=s.t; control.checks.longitude={good};
control.checks.altitude={good};
out=nav_geodetic_output(pole,eye(21),env,control); assert(strcmp(out.quality.longitude,'UNDEFINED'));


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


function requirement_correlated_measurement
[s,~]=fixture; P=eye(21); H=zeros(1,21); H(1)=1;
o=struct('frame_id',s.frame_id,'t_measure',s.t,'noise_model','joint_required', ...
    'H',H,'R',1,'residual',2,'state_noise_cov',zeros(21,1));
o.state_noise_cov(1)=.5;
[x,Q,info]=nav_correct(s,P,o,100); assert(info.accepted);
close(x.p,s.p+[1;0;0],1e-12); close(Q(1,1),.25,1e-12); close(info.S,3,1e-12);
o.state_noise_cov(1)=-.5;
[x,Q,info]=nav_correct(s,P,o,100);
close(x.p,s.p+[1;0;0],1e-12); close(Q(1,1),.75,1e-12); close(info.S,1,1e-12);
o.state_noise_cov(1)=2;
fails(@()nav_correct(s,P,o,100),'libekf:Correlation');
o.state_noise_cov(1)=.5; o.noise_model='independent';
fails(@()nav_correct(s,P,o,100),'libekf:Correlation');

% Repeated ranges share station-position error; compare with the joint linear solution.
[s,e]=fixture; s.q=[1;0;0;0]; P=eye(21); P(1,1)=4;
m=packet(s,'dme_range'); m.lever=zeros(3,1); m.station=s.p-[10;0;0];
m.z=10; m.R=10; m.noise_model='joint_required';
m.state_noise_cov=zeros(21,1); m.state_noise_frame_id=s.frame_id;
o=nav_measurement(s,m,e); [a,Q,info]=nav_correct(s,P,o,100); assert(info.accepted);
close(a.p,s.p,1e-12); close(Q(1,1),20/7,1e-12);
% Station variance 9 and independent range variance 1. Conditioning on the first
% zero residual gives Var(v2|z1)=10-81/14 and Cov(error_x,v2|z1)=-36/14.
m.z=11; m.R=59/14; m.state_noise_cov=zeros(21,1); m.state_noise_cov(1)=-18/7;
m.state_noise_frame_id=s.frame_id;
o=nav_measurement(a,m,e); [b,V,info]=nav_correct(a,Q,o,100); assert(info.accepted);
close(b.p,s.p+[4/27;0;0],1e-12); close(V(1,1),76/27,1e-12);
% Joint observation covariance [10 9;9 10] yields the same posterior.
H=zeros(2,21); H(:,1)=1;
joint=struct('frame_id',s.frame_id,'t_measure',s.t,'noise_model','joint_required', ...
    'H',H,'R',[10 9;9 10],'residual',[0;1],'state_noise_cov',zeros(21,2));
[batch,B,info]=nav_correct(s,P,joint,100); assert(info.accepted);
close(batch.p,b.p,1e-12); close(B,V,1e-12);

% The uncertain-reference model computes the conditional reference errors itself.
m=packet(s,'cooperative_range'); m.lever=zeros(3,1); m.station=s.p-[10;0;0];
m.peer_covariance=diag([9 1 1]); m.peer_cross_cov=zeros(21,3);
m.noise_model='joint_required'; m.z=10; m.R=1;
o=nav_measurement(s,m,e); [a,Q,first_info]=nav_correct(s,P,o,100);
assert(first_info.accepted); ref=first_info.peer_update;
bad_reference=o; bad_reference.peer_noise_cov=2;
fails(@()nav_correct(s,P,bad_reference,100),'libekf:Correlation');
close(ref.covariance(1,1),45/14,1e-12); close(ref.cross_cov(1,1),18/7,1e-12);
m.station=ref.position; m.peer_covariance=ref.covariance; m.peer_cross_cov=ref.cross_cov; m.z=11;
o=nav_measurement(a,m,e); [b,V,second_info]=nav_correct(a,Q,o,100);
assert(second_info.accepted); close(b.p,batch.p,1e-12); close(V,B,1e-12);
o.residual=1e6; [unchanged,unchangedP,rejected]=nav_correct(a,Q,o,1);
assert(~rejected.accepted&&~isfield(rejected,'peer_update')&&isequal(unchanged,a)&&isequal(unchangedP,Q));

% The common handler returns reference errors in common axes, even for rotated sensors.
T=[0 -1 0;1 0 0;0 0 1]; offset=[7;-3;2];
sensor=struct('id','rotated_peer','origin_ecef',e.origin_ecef+e.C_e_l*offset,'C_e_s',e.C_e_l*T);
filter=struct('state',s,'P',P,'used',{{}},'capacity',10);
m=packet(s,'cooperative_range'); m.frame='SENSOR_LOCAL'; m.frame_id=sensor.id;
m.lever=zeros(3,1); m.station=T.'*(s.p-[10;0;0]-offset);
m.peer_covariance=T.'*diag([9 1 1])*T; m.peer_cross_cov=zeros(21,3);
m.peer_cross_frame_id=s.frame_id; m.noise_model='joint_required'; m.R=1; m.z=10;
m.lineage={'reference_range_first'};
[filter,r1]=nav_assimilate(filter,m,sensor,e,100);
assert(r1.accepted&&isfield(r1,'peer_update'));
close(r1.peer_update.position,ref.position,1e-8);
close(r1.peer_update.covariance,ref.covariance,1e-8);
close(r1.peer_update.cross_cov,ref.cross_cov,1e-8);
assert(strcmp(r1.peer_update.frame_id,s.frame_id)&&r1.peer_update.t==s.t);
m.station=T.'*(r1.peer_update.position-offset);
m.peer_covariance=T.'*r1.peer_update.covariance*T;
m.peer_cross_cov=r1.peer_update.cross_cov*T; m.z=11; m.lineage={'reference_range_second'};
[filter,r2]=nav_assimilate(filter,m,sensor,e,100);
assert(r2.accepted&&numel(filter.used)==2);
close(filter.state.p,batch.p,1e-8); close(filter.P,B,1e-8);
% A duplicate range must not emit an update for the reference either.
[same,r2]=nav_assimilate(filter,m,sensor,e,100);
assert(~r2.accepted&&~isfield(r2,'peer_update')&&isequal(same,filter));

% A new independent position observation also conditions the correlated reference.
retained=first_info.peer_update; retained.noise_independent=true;
position=packet(s,'gnss_position'); position.lever=zeros(3,1); position.z=s.p+[1;0;0];
position.R=eye(3); position_obs=nav_measurement(s,position,e);
prior=eye(21); prior(1,1)=20/7;
[corrected,posterior,other]=nav_correct(s,prior,position_obs,100,retained);
assert(other.accepted); close(corrected.p,s.p+[20/27;0;0],1e-12);
close(posterior(1,1),20/27,1e-12);
close(other.peer_update.position,retained.position+[2/3;0;0],1e-12);
close(other.peer_update.covariance(1,1),3/2,1e-12);
close(other.peer_update.cross_cov(1,1),2/3,1e-12);
bad_retained=retained; bad_retained.noise_independent=false;
fails(@()nav_correct(s,prior,position_obs,100,bad_retained),'libekf:Correlation');
bad_retained=retained; bad_retained.t=s.t-1;
fails(@()nav_correct(s,prior,position_obs,100,bad_retained),'libekf:Epoch');

% Prepared retained-reference updates use the same accepted/rejected packet path.
position.frame='SENSOR_LOCAL'; position.retained_reference=retained;
sensor=struct('id',s.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
filter=struct('state',s,'P',prior,'used',{{}},'capacity',10);
[next,accepted]=nav_assimilate(filter,position,sensor,e,100);
assert(accepted.accepted); close(next.state.p,corrected.p,1e-12); close(next.P,posterior,1e-12);
close(accepted.peer_update.position,other.peer_update.position,1e-12);
close(accepted.peer_update.cross_cov,other.peer_update.cross_cov,1e-12);
[same,repeated]=nav_assimilate(next,position,sensor,e,100);
assert(~repeated.accepted&&isequal(same,next)&&~isfield(repeated,'peer_update'));
position.retained_reference.t=s.t-1;
[same,rejected]=nav_assimilate(filter,position,sensor,e,100);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:Epoch')&& ...
    isequal(same,filter)&&~isfield(rejected,'peer_update'));

% Reference error is propagated jointly, not reset at the next epoch.
transition=eye(21); transition(1,4)=2;
retained.cross_cov(4,1)=.1;
motion=struct('A',eye(3),'shift',zeros(3,1),'Q',zeros(3), ...
    't_start',retained.t,'t_end',retained.t+2,'frame_id',retained.frame_id,'noise_independent',true);
propagated=nav_predict_reference(retained,transition,motion);
close(propagated.position,retained.position,0); close(propagated.covariance,retained.covariance,0);
close(propagated.cross_cov(1,1),retained.cross_cov(1,1)+.2,1e-12);
assert(propagated.t==motion.t_end&&~isfield(propagated,'noise_independent'));
motion.A=diag([1 2 1]); motion.shift=[3;0;0]; motion.Q=diag([.5 .2 0]);
propagated=nav_predict_reference(retained,transition,motion);
joint=[prior retained.cross_cov;retained.cross_cov.' retained.covariance];
F=blkdiag(transition,motion.A); expected=F*joint*F.'+blkdiag(zeros(21),motion.Q);
close(propagated.covariance,expected(22:24,22:24),1e-12);
close(propagated.cross_cov,expected(1:21,22:24),1e-12);
close(propagated.position,motion.A*retained.position+[3;0;0],1e-12);
bad_motion=motion; bad_motion.noise_independent=false;
fails(@()nav_predict_reference(retained,transition,bad_motion),'libekf:Correlation');
bad_motion=motion; bad_motion.Q=diag([-1 0 0]);
fails(@()nav_predict_reference(retained,transition,bad_motion),'libekf:Covariance');
bad_motion=motion; bad_motion.t_start=retained.t-1;
fails(@()nav_predict_reference(retained,transition,bad_motion),'libekf:Epoch');

% Impossible peer dependence must be rejected even perpendicular to the measured range.
[s,e]=fixture; m=packet(s,'cooperative_range'); m.lever=zeros(3,1);
m.station=s.p+[10;0;0]; m.z=10; m.peer_covariance=eye(3);
m.peer_cross_cov=zeros(21,3); m.peer_cross_cov(2,2)=2;
m.noise_model='joint_required';
o=nav_measurement(s,m,e);
fails(@()nav_correct(s,eye(21),o,100),'libekf:Correlation');
m.peer_cross_cov(2,2)=.5; o=nav_measurement(s,m,e);
[~,~,info]=nav_correct(s,eye(21),o,100); assert(info.accepted);
m.peer_covariance=zeros(3); m.peer_cross_cov=zeros(21,3);
o=nav_measurement(s,m,e); [~,~,info]=nav_correct(s,eye(21),o,100); assert(info.accepted);
m.peer_cross_cov(2,2)=1e-20; o=nav_measurement(s,m,e);
fails(@()nav_correct(s,eye(21),o,100),'libekf:Correlation');

% The same dependent measurement in rotated sensor axes gives the same correction.
[s,e]=fixture; sensor=struct('id',s.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
f=struct('state',s,'P',eye(21),'used',{{}},'capacity',10);
m=packet(s,'external_position'); m.lever=zeros(3,1); m.frame='SENSOR_LOCAL';
m.z=s.p+[2;0;0]; m.noise_model='joint_required';
m.state_noise_frame_id=s.frame_id; m.state_noise_cov=zeros(21,3); m.state_noise_cov(1,1)=.5;
[normal,info]=nav_assimilate(f,m,sensor,e,100); assert(info.accepted);
close(normal.state.p,s.p+[1;0;0],1e-12); close(normal.P(1,1),.25,1e-12);
turn=nav_rotation([cos(.3);0;0;sin(.3)]);
adapter=struct('C_e_l',turn,'origin_ecef',[10;20;30],'frame_id','alternate');
native=m; native.frame='ECEF'; native=nav_packet_from_ecef(native,adapter); native.frame='SENSOR_LOCAL';
other=struct('id','alternate','origin_ecef',e.origin_ecef+[10;20;30],'C_e_s',turn);
[rotated,info]=nav_assimilate(f,native,other,e,100); assert(info.accepted);
close(rotated.state.p,normal.state.p,1e-11); close(rotated.P,normal.P,1e-11);
native.state_noise_frame_id='wrong_common_axes';
[same,info]=nav_assimilate(f,native,other,e,100); assert(~info.accepted&&isequal(same,f));

% Full visual pose with known dependence on the two-epoch prior.
a=s; a.p=zeros(3,1); a.q=[1;0;0;0]; b=a; b.p=[1;0;0]; b.t=a.t+1;
f=struct('state',b,'reference',a,'P',eye(21)*.01,'jointP',eye(42)*.01, ...
    'used',{{}},'capacity',10);
m=struct('model','visual_odometry','mode','pose','source','camera', ...
    'lineage',{{'image1','image2'}},'quality','VALID','noise_model','joint_required', ...
    'correlation_group','camera','frame_id',a.frame_id,'t_start',a.t,'t_end',b.t, ...
    't_received',b.t,'units',{{'m';'m';'m';'rad';'rad';'rad'}}, ...
    'z',[2;0;0;1;0;0;0],'R',eye(6)*.01,'camera_to_body',eye(3),'lever',zeros(3,1), ...
    'scale_valid',true,'calibration_valid',true,'max_angle',.5, ...
    'state_noise_frame_id',a.frame_id,'state_noise_cov',zeros(42,6));
m.state_noise_cov(22,1)=.002;
[updated,info]=nav_assimilate(f,m,[],[],100); assert(info.accepted);
close(updated.reference.p(1),-.01/.034,1e-12);
close(updated.state.p(1),1+.012/.034,1e-12);
m.state_noise_cov=zeros(21,6);
[same,info]=nav_assimilate(f,m,[],[],100); assert(~info.accepted&&isequal(same,f));


end

function requirement_joint_correction
[a,~]=fixture; a.p=zeros(3,1); a.q=[1;0;0;0]; b=a; b.t=a.t+1; b.p=[1;0;0];
m=struct('source','camera','lineage',{{'first','second'}},'quality','VALID', ...
    'noise_model','independent','correlation_group','camera','frame_id',a.frame_id, ...
    't_start',a.t,'t_end',b.t,'t_received',b.t,'units',{{'m';'m';'m'}}, ...
    'z',[3;0;0],'R',eye(3),'camera_to_body',eye(3),'lever',zeros(3,1), ...
    'scale_valid',true,'calibration_valid',true);
P=[eye(21) .5*eye(21);.5*eye(21) eye(21)];
o=nav_odometry_measurement(a,b,m);
[x,Q,info]=nav_correct({a,b},P,o,100); assert(info.accepted);
close(x{1}.p,[-.5;0;0],1e-12); close(x{2}.p,[1.5;0;0],1e-12);
close(Q([1 22],[1 22]),[.875 .625;.625 .875],1e-12);
close([.5 .5]*Q([1 22],[1 22])*[.5;.5],.75,1e-12);

% Changing Earth-fixed local axes must preserve the complete two-epoch correction.
[~,from]=fixture; to=from; to.frame_id='pair_rotated';
to.C_e_l=nav_rotation([.9;.1;.2;.3]/norm([.9;.1;.2;.3]));
to.origin_ecef=from.origin_ecef+[10;20;30];
[native,Pnative]=nav_change_frame({a,b},P,from,to);
shifted=m; shifted.frame_id=to.frame_id;
onative=nav_odometry_measurement(native{1},native{2},shifted);
close(onative.residual,o.residual,1e-12);
[xn,Qn,info]=nav_correct(native,Pnative,onative,100); assert(info.accepted);
[back,Qback]=nav_change_frame(xn,Qn,to,from);
for j=1:2
    close(back{j}.p,x{j}.p,1e-11);
    close(nav_rotation(back{j}.q),nav_rotation(x{j}.q),1e-12);
end
close(Qback,Q,1e-11);

[x,Q,info]=nav_correct({a,b},P,o,.01);
assert(~info.accepted&&isequal(x,{a,b})&&isequal(Q,P));
% Perfectly shared errors are legal; the relative observation adds no absolute information.
b.p=a.p; P=[eye(21) eye(21);eye(21) eye(21)]; m.z=zeros(3,1);
o=nav_odometry_measurement(a,b,m); [x,Q,info]=nav_correct({a,b},P,o,100);
assert(info.accepted); close(Q,P,1e-12); close(x{1}.p,a.p,1e-12);
P(1,22)=2; P(22,1)=2;
fails(@()nav_correct({a,b},P,o,100),'libekf:Covariance');

% A delayed observation constrains its own epoch, then the correlated current state.
[a,e]=fixture; a.p=zeros(3,1); b=a; b.p=[10;0;0]; b.t=a.t+1;
f=struct('state',b,'P',eye(21),'reference',a, ...
    'jointP',[eye(21) .5*eye(21);.5*eye(21) eye(21)],'used',{{}},'capacity',10);
sensor=struct('id',a.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
g=packet(a,'gnss_position'); g.frame='SENSOR_LOCAL'; g.lever=zeros(3,1);
g.z=[2;0;0]; g.t_received=b.t+.1;
% Age is measured from observation, not from delivery; late data keep their own epoch.
policy=struct('now',g.t_received,'max_age',g.t_received-g.t_measure);
[timed,info]=nav_assimilate(f,g,sensor,e,100,policy); assert(info.accepted);
policy.max_age=policy.max_age-.01;
[same,info]=nav_assimilate(f,g,sensor,e,100,policy);
assert(~info.accepted&&strcmp(info.reason,'libekf:Stale')&&isequal(same,f));
policy.max_age=10; policy.now=g.t_received-.01;
[same,info]=nav_assimilate(f,g,sensor,e,100,policy);
assert(~info.accepted&&strcmp(info.reason,'libekf:Time')&&isequal(same,f));
policy.now=g.t_received+20;
[same,info]=nav_assimilate(f,g,sensor,e,100,policy);
assert(~info.accepted&&strcmp(info.reason,'libekf:Stale')&&isequal(same,f));
[updated,info]=nav_assimilate(f,g,sensor,e,100); assert(info.accepted);
close(timed.P,updated.P,1e-12);
close(updated.reference.p,[1;0;0],1e-12); close(updated.state.p,[10.5;0;0],1e-12);
close(updated.jointP([1 22],[1 22]),[.5 .25;.25 .875],1e-12);
assert(updated.reference.t==a.t&&updated.state.t==b.t);
[same,info]=nav_assimilate(updated,g,sensor,e,100); assert(~info.accepted&&isequal(same,updated));
g.t_measure=a.t+.5; g.lineage={'unretained_epoch'};
[same,info]=nav_assimilate(updated,g,sensor,e,100); assert(~info.accepted&&isequal(same,updated));

% Complete shared path: save reference -> inertial interval -> GNSS -> camera pair.
[a,e]=fixture; a.p=zeros(3,1); a.q=[1;0;0;0]; a.v=[1;0;0];
a.ba=zeros(3,1); a.bg=zeros(3,1);
e.Omega=zeros(3,1); e.gravity_origin=zeros(3,1); e.gravity_gradient=zeros(3);
e.gravity_radius=1e5; e.max_step=.1; e.state_noise_density=eye(21)*1e-8;
f=struct('state',a,'P',eye(21),'reference',a,'jointP',[eye(21) eye(21);eye(21) eye(21)], ...
    'used',{{}},'capacity',100);
u=struct('model','inertial','source','imu','lineage',{{'pair_imu'}},'quality','VALID', ...
    'frame','BODY','accel_unit','m/s^2','gyro_unit','rad/s','t_start',a.t, ...
    't_end',a.t+.1,'t_received',a.t+.1,'accel',zeros(3,1),'gyro',zeros(3,1), ...
    'noise_density',eye(6)*1e-6);
policy=struct('now',u.t_received+1,'max_age',.5);
[same,info]=nav_assimilate(f,u,[],e,100,policy);
assert(~info.accepted&&strcmp(info.reason,'libekf:Stale')&&isequal(same,f));
[f,info]=nav_assimilate(f,u,[],e,100); assert(info.accepted);
sensor=struct('id',a.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
g=packet(f.state,'gnss_position'); g.frame='SENSOR_LOCAL'; g.lever=zeros(3,1);
g.z=f.state.p+[.1;0;0];
[f,info]=nav_assimilate(f,g,sensor,e,100); assert(info.accepted);
assert(f.reference.p(1)>0); close(f.P,f.jointP(22:42,22:42),1e-12);
m.model='visual_odometry'; m.t_start=a.t; m.t_end=f.state.t; m.t_received=f.state.t;
m.z=[.12;0;0]; m.R=eye(3)*.001;
policy=struct('now',m.t_received+1,'max_age',.5);
[same,info]=nav_assimilate(f,m,[],e,100,policy);
assert(~info.accepted&&strcmp(info.reason,'libekf:Stale')&&isequal(same,f));
before=f; [f,info]=nav_assimilate(f,m,[],e,100); assert(info.accepted);
assert(norm(f.state.p-before.state.p)>0); close(f.P,f.jointP(22:42,22:42),1e-12);
[unchanged,info]=nav_assimilate(f,m,[],e,100); assert(~info.accepted&&isequal(unchanged,f));
m.lineage={'new_frame_a','new_frame_b'}; m.z=[1000;0;0];
[unchanged,info]=nav_assimilate(f,m,[],e,.01); assert(~info.accepted&&isequal(unchanged,f));

% Replacing the saved epoch must not create information or erase provenance.
saved=nav_save_reference(f); assert(isequal(saved.state,f.state)&&isequal(saved.P,f.P));
assert(isequal(saved.used,f.used)&&isequal(saved.reference,saved.state));
close(saved.jointP,[f.P f.P;f.P f.P],1e-12);
u.t_start=saved.state.t; u.t_end=u.t_start+.1; u.t_received=u.t_end;
u.lineage={'second_pair_imu'};
[advanced,info]=nav_assimilate(saved,u,[],e,100); assert(info.accepted);
m.t_start=saved.reference.t; m.t_end=advanced.state.t; m.t_received=advanced.state.t;
m.lineage={'new_reference_frame','new_current_frame'};
m.z=advanced.state.p-advanced.reference.p;
[accepted,info]=nav_assimilate(advanced,m,[],e,100); assert(info.accepted);
assert(accepted.reference.t==saved.state.t&&accepted.state.t==u.t_end);
m.lineage={'first','second'};
[unchanged,info]=nav_assimilate(accepted,m,[],e,100);
assert(~info.accepted&&isequal(unchanged,accepted));


end

function requirement_interval_observations
[a,e]=fixture; a.p=zeros(3,1); a.v=[2;0;0]; a.q=[1;0;0;0];
a.ba=zeros(3,1); a.bg=zeros(3,1); e.Omega=zeros(3,1);
e.gravity_origin=zeros(3,1); e.gravity_gradient=zeros(3); e.gravity_radius=1e5;
e.max_step=1; e.state_noise_density=zeros(21);
u=struct('model','inertial','source','imu','lineage',{{'one_imu_interval'}},'quality','VALID', ...
    'frame','BODY','accel_unit','m/s^2','gyro_unit','rad/s','t_start',a.t, ...
    't_end',a.t+1,'t_received',a.t+1,'accel',zeros(3,1),'gyro',zeros(3,1), ...
    'noise_density',zeros(6));
f=struct('state',a,'P',eye(21),'used',{{}},'capacity',20); f=nav_save_reference(f);
sensor=struct('id',a.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
b=a; b.t=a.t+.5; b.p=[1;0;0];
m=packet(b,'gnss_position'); m.frame='SENSOR_LOCAL'; m.lever=zeros(3,1);
m.z=[2;0;0]; m.R=eye(3); m.lineage={'middle_position'};
[x,info]=nav_assimilate_interval(f,u,{m},{sensor},e,100);
assert(info.accepted&&info.observations{1}.accepted&&x.state.t==u.t_end);
% A late observation in an earlier interval reconstructs the same estimate as timely processing.
second=u; second.t_start=u.t_end; second.t_end=u.t_end+1;
second.t_received=second.t_end; second.lineage={'second_imu_interval'};
late=m; late.t_received=second.t_end+.5;
boundary=m; boundary.t_measure=u.t_end; boundary.t_received=second.t_end;
boundary.z=[2;0;0]; boundary.lineage={'boundary_position'};
[without,ri]=nav_replay(f,{u,second},{boundary},{sensor},e,100); assert(ri.accepted);
[replayed,ri]=nav_replay(f,{u,second},{boundary,late},{sensor,sensor},e,100);
assert(ri.accepted&&all(cellfun(@(z)z.accepted,ri.observations)));
[timely,ti]=nav_assimilate_interval(f,u,{m,boundary},{sensor,sensor},e,100); assert(ti.accepted);
[timely,ti]=nav_assimilate(timely,second,[],e,100); assert(ti.accepted);
close(replayed.state.p,timely.state.p,1e-12); close(replayed.jointP,timely.jointP,1e-12);
assert(norm(replayed.state.p-without.state.p)>0&&numel(replayed.used)==4);
[again,ri]=nav_replay(f,{u,second},{boundary,late},{sensor,sensor},e,100);
assert(ri.accepted&&isequal(again,replayed)); % Replaying does not add information twice.
% New late data are admitted by their own age; historical data retain admission.
policy=struct('now',late.t_received,'max_age',late.t_received-late.t_measure);
[again,ri]=nav_replay(f,{u,second},{boundary,late},{sensor,sensor},e,100,{[],policy});
assert(ri.accepted&&ri.observations{2}.accepted&&isequal(again,replayed));
policy.max_age=policy.max_age-.1;
[again,ri]=nav_replay(f,{u,second},{boundary,late},{sensor,sensor},e,100,{[],policy});
assert(ri.accepted&&ri.observations{1}.accepted&&~ri.observations{2}.accepted);
assert(strcmp(ri.observations{2}.reason,'libekf:Stale')&&isequal(again,without));
policy.max_age=10; policy.now=late.t_received-.1;
[again,ri]=nav_replay(f,{u,second},{boundary,late},{sensor,sensor},e,100,{[],policy});
assert(ri.accepted&&strcmp(ri.observations{2}.reason,'libekf:Time')&&isequal(again,without));
policy.now=NaN;
[again,ri]=nav_replay(f,{u,second},{boundary,late},{sensor,sensor},e,100,{[],policy});
assert(ri.accepted&&strcmp(ri.observations{2}.reason,'libekf:Time')&&isequal(again,without));
old=late; old.t_measure=a.t-.5; old.lineage={'before_snapshot'};
[again,ri]=nav_replay(f,{u,second},{boundary,late,old},{sensor,sensor,sensor},e,100);
assert(ri.accepted&&~ri.observations{3}.accepted&&isequal(again,replayed));
broken=second; broken.t_start=broken.t_start+.1;
[same,ri]=nav_replay(f,{u,broken},{late},{sensor},e,100);
assert(~ri.accepted&&isequal(same,f)); % No partially reconstructed result escapes a history gap.
% Camera reference changes are replayed as actions, without inventing measurements.
cam=struct('model','visual_odometry','source','camera','lineage',{{'frame_a','frame_b'}}, ...
    'quality','VALID','noise_model','independent','correlation_group','camera','frame_id',a.frame_id, ...
    't_start',a.t,'t_end',u.t_end,'t_received',u.t_end,'units',{{'m';'m';'m'}}, ...
    'z',[2;0;0],'R',eye(3),'camera_to_body',eye(3),'lever',zeros(3,1), ...
    'scale_valid',true,'calibration_valid',true);
cam2=cam; cam2.t_start=second.t_start; cam2.t_end=second.t_end;
cam2.t_received=second.t_end; cam2.lineage={'frame_c','frame_d'};
save=struct('model','save_reference','t_measure',u.t_end);
[replayed,ri]=nav_replay(f,{u,second},{late,cam,save,cam2},{sensor,[],[],[]},e,100);
camera_policy=struct('now',second.t_end,'max_age',second.t_end-cam.t_end);
[camera_admitted,camera_check]=nav_replay(f,{u,second},{late,cam,save,cam2}, ...
    {sensor,[],[],[]},e,100,{[],camera_policy,[],[]});
assert(camera_check.accepted&&camera_check.observations{2}.accepted&&isequal(camera_admitted,replayed));
camera_policy.max_age=camera_policy.max_age-.01;
[~,camera_check]=nav_replay(f,{u,second},{late,cam,save,cam2}, ...
    {sensor,[],[],[]},e,100,{[],camera_policy,[],[]});
assert(camera_check.accepted&&strcmp(camera_check.observations{2}.reason,'libekf:Stale'));
assert(ri.accepted&&all(cellfun(@(z)z.accepted,ri.observations)));
[timely,ti]=nav_assimilate_interval(f,u,{m,cam},{sensor,[]},e,100); assert(ti.accepted);
timely=nav_save_reference(timely);
[timely,ti]=nav_assimilate_interval(timely,second,{cam2},{[]},e,100); assert(ti.accepted);
close(replayed.state.p,timely.state.p,1e-12); close(replayed.jointP,timely.jointP,1e-12);
assert(replayed.reference.t==u.t_end&&numel(replayed.used)==7);
% A replacement reference does not make an already used frame independent.
cam2.lineage={'frame_b','frame_d'};
[with_repeat,ri]=nav_replay(f,{u,second},{late,cam,save,cam2},{sensor,[],[],[]},e,100);
assert(ri.accepted&&~ri.observations{4}.accepted&&numel(with_repeat.used)==5);
% Reversing equal-time reference/action order cannot silently reinterpret the camera pair.
[~,ri]=nav_replay(f,{u,second},{save,cam},{[],[]},e,100);
assert(ri.accepted&&ri.observations{1}.accepted&&~ri.observations{2}.accepted);
% Independent scalar p-v-bias state transition, conditioning at t=1/2, then to t=1.
F=[1 .5 -.125;0 1 -.5;0 0 1]; Pmid=F*F.'; H=[1 0 0];
K=Pmid*H.'/(H*Pmid*H.'+1); correction=F*K;
expected=F*(Pmid-K*(H*Pmid))*F.';
close(x.state.p(1),2+correction(1),1e-12);
close(x.state.v(1),2+correction(2),1e-12);
close(x.P([1 4 10],[1 4 10]),expected,1e-12);
assert(numel(x.used)==2&&sum(strcmp(x.used,'one_imu_interval'))==1);
close(x.P,x.jointP(22:42,22:42),0);
% Input delivery order cannot move observations to another epoch.
c=m; c.t_measure=a.t+.75; c.t_received=u.t_end; c.lineage={'later_position'}; c.z=[1.5;0;0];
[y,info]=nav_assimilate_interval(f,u,{c,m},{sensor,sensor},e,100); assert(info.accepted);
[z,info]=nav_assimilate_interval(f,u,{m,c},{sensor,sensor},e,100); assert(info.accepted);
close(y.state.p,z.state.p,1e-12); close(y.jointP,z.jointP,1e-12);
% Each source retains its age limit after chronological sorting.
policy=struct('now',u.t_end,'max_age',0,'observation_max_age',[.2 .6]);
[y,info]=nav_assimilate_interval(f,u,{c,m},{sensor,sensor},e,100,policy);
assert(info.accepted&&~info.observations{1}.accepted&&info.observations{2}.accepted);
assert(strcmp(info.observations{1}.reason,'libekf:Stale'));
close(y.state.p,x.state.p,1e-12); close(y.jointP,x.jointP,1e-12);
assert(isequal(y.used,x.used));
policy.observation_max_age=[.6 .2];
[z,info]=nav_assimilate_interval(f,u,{m,c},{sensor,sensor},e,100,policy);
assert(info.accepted&&info.observations{1}.accepted&&~info.observations{2}.accepted);
assert(isequal(y,z));
policy.observation_max_age=-1;
fails(@()nav_assimilate_interval(f,u,{m,c},{sensor,sensor},e,100,policy),'libekf:Time');
% Rejected observations do not stop valid inertial motion or consume source IDs.
u.noise_density=diag([4 4 4 0 0 0]);
bad=m; bad.quality='FAILED';
[y,info]=nav_assimilate_interval(f,u,{bad},{sensor},e,100);
assert(info.accepted&&~info.observations{1}.accepted&&numel(y.used)==1);
[plain,info]=nav_assimilate(f,u,[],e,100); assert(info.accepted);
close(y.state.p,plain.state.p,1e-12); close(y.jointP,plain.jointP,1e-12);
% Rejected data must not change inertial quadrature during a rotating acceleration.
u.accel=[2;1;0]; u.gyro=[0;0;.4];
[plain,info]=nav_assimilate(f,u,[],e,100); assert(info.accepted);
[y,info]=nav_assimilate_interval(f,u,{bad},{sensor},e,100);
assert(info.accepted&&~info.observations{1}.accepted);
close(y.state.p,plain.state.p,1e-12); close(y.state.v,plain.state.v,1e-12);
close(y.jointP,plain.jointP,1e-12);
bad.quality='VALID'; bad.z=[1e6;0;0];
[y,info]=nav_assimilate_interval(f,u,{bad},{sensor},e,100);
assert(info.accepted&&~info.observations{1}.accepted);
close(y.state.p,plain.state.p,1e-12); close(y.jointP,plain.jointP,1e-12);
% Rejected data preceding a good observation must not alter that good update either.
[expected,info]=nav_assimilate_interval(f,u,{c},{sensor},e,100); assert(info.accepted);
[y,info]=nav_assimilate_interval(f,u,{bad,c},{sensor,sensor},e,100);
assert(info.accepted&&~info.observations{1}.accepted&&info.observations{2}.accepted);
close(y.state.p,expected.state.p,1e-12); close(y.jointP,expected.jointP,1e-12);
% Reusing a complete interval is rejected without changing the filter.
[same,info]=nav_assimilate_interval(x,u,{m},{sensor},e,100);
assert(~info.accepted&&isequal(same,x));
end

function requirement_relative_camera_motion
[a,~]=fixture; a.p=zeros(3,1); a.q=[1;0;0;0]; b=a; b.t=a.t+1; b.p=[3;4;0];
m=struct('source','camera','lineage',{{'frame100','frame101'}},'quality','VALID', ...
    'noise_model','independent','correlation_group','camera','frame_id',a.frame_id, ...
    't_start',a.t,'t_end',b.t,'t_received',b.t+.1,'units',{{'m';'m';'m'}}, ...
    'z',[3;4;0],'R',eye(3),'camera_to_body',eye(3),'lever',[1;0;0], ...
    'scale_valid',true,'calibration_valid',true);
o=nav_odometry_measurement(a,b,m); close(o.predicted,[3;4;0],1e-12);
b.p=zeros(3,1); b.q=[cos(pi/4);0;0;sin(pi/4)];
close(nav_odometry_measurement(a,b,m).predicted,[-1;1;0],1e-12);
o=nav_odometry_measurement(a,b,m); numeric=zeros(3,42);
for j=1:42
    d=zeros(21,1); index=mod(j-1,21)+1; d(index)=1e-6;
    if j<=21
        plus=nav_odometry_measurement(nav_apply_error(a,d),b,m);
        minus=nav_odometry_measurement(nav_apply_error(a,-d),b,m);
    else
        plus=nav_odometry_measurement(a,nav_apply_error(b,d),m);
        minus=nav_odometry_measurement(a,nav_apply_error(b,-d),m);
    end
    numeric(:,j)=-(plus.residual-minus.residual)/(2e-6);
end
close(numeric,o.H,1e-8);
a.p=a.p+[100;200;300]; b.p=b.p+[100;200;300];
close(nav_odometry_measurement(a,b,m).predicted,o.predicted,1e-12);
m.scale_valid=false; fails(@()nav_odometry_measurement(a,b,m),'libekf:Quality');
m.scale_valid=true; m.t_start=a.t-1;
fails(@()nav_odometry_measurement(a,b,m),'libekf:Epoch');

% A monocular direction is invariant under positive scaling, not under reversal.
a.p=zeros(3,1); b.p=[3;4;0]; a.q=[1;0;0;0]; b.q=a.q;
m.t_start=a.t; m.mode='direction'; m.scale_valid=false; m.lever=zeros(3,1);
m.z=[.6;.8;0]; m.tangent_basis=[-.8 0;.6 0;0 1]; m.R=eye(2)*.01;
m.units={'1';'1'}; m.min_displacement=.01; m.max_angle=.3;
o=nav_odometry_measurement(a,b,m); close(o.residual,zeros(2,1),1e-12);
close(o.H(:,22:24)*[3;4;0],zeros(2,1),1e-12);
b.p=2*b.p; close(nav_odometry_measurement(a,b,m).residual,o.residual,1e-12);
b.p=[3;4;.2]; o=nav_odometry_measurement(a,b,m);
for j=1:3
    plus=b; minus=b; plus.p(j)=plus.p(j)+1e-6; minus.p(j)=minus.p(j)-1e-6;
    op=nav_odometry_measurement(a,plus,m); om=nav_odometry_measurement(a,minus,m);
    close(-(op.residual-om.residual)/2e-6,o.H(:,21+j),1e-8);
end
b.p=[-3;-4;0]; fails(@()nav_odometry_measurement(a,b,m),'libekf:Geometry');
b.p=zeros(3,1); fails(@()nav_odometry_measurement(a,b,m),'libekf:Geometry');

% Relative rotation remains measurable when camera centres do not translate.
m.mode='rotation'; m.units={'rad';'rad';'rad'}; m.R=eye(3)*.001; m.max_angle=.5;
b.q=[cos(.1);0;0;sin(.1)]; m.z=[cos(.15);0;0;sin(.15)];
o=nav_odometry_measurement(a,b,m); close(o.residual,[0;0;.1],1e-12);
close(o.H(:,1:6),zeros(3,6),1e-12);
for which=1:2
    for j=7:9
        d=zeros(21,1); d(j)=1e-6;
        if which==1
            op=nav_odometry_measurement(nav_apply_error(a,d),b,m);
            om=nav_odometry_measurement(nav_apply_error(a,-d),b,m);
        else
            op=nav_odometry_measurement(a,nav_apply_error(b,d),m);
            om=nav_odometry_measurement(a,nav_apply_error(b,-d),m);
        end
        close(-(op.residual-om.residual)/2e-6,o.H(:,(which-1)*21+j),1e-8);
    end
end
m.z=-m.z; close(nav_odometry_measurement(a,b,m).residual,o.residual,1e-12);
f=struct('state',b,'reference',a,'P',eye(21)*.01,'jointP',eye(42)*.01, ...
    'used',{{}},'capacity',10); m.model='visual_odometry';
[updated,result]=nav_assimilate(f,m,[],[],100); assert(result.accepted);
after=nav_odometry_measurement(updated.reference,updated.state,m);
assert(norm(after.residual)<norm(o.residual));
m.max_angle=.01; fails(@()nav_odometry_measurement(a,b,m),'libekf:Geometry');

% One camera pair carries translation and rotation with their full covariance.
m.max_angle=.5; m.mode='pose'; m.scale_valid=true; b.p=[1;0;0];
m.z=[1.1;0;0;cos(.15);0;0;sin(.15)];
m.units={'m';'m';'m';'rad';'rad';'rad'};
m.R=eye(6)*.01; m.R(1,6)=.002; m.R(6,1)=.002;
o=nav_odometry_measurement(a,b,m); close(o.residual,[.1;0;0;0;0;.1],1e-12);
close(o.R(1,6),.002,1e-12); assert(size(o.H,1)==6);
f.state=b; f.reference=a;
[updated,result]=nav_assimilate(f,m,[],[],100); assert(result.accepted);
assert(numel(updated.used)==2);
[same,result]=nav_assimilate(updated,m,[],[],100); assert(~result.accepted&&isequal(same,updated));
m.mode='direction_rotation'; m.scale_valid=false; m.z(1:3)=[1;0;0];
m.tangent_basis=[0 0;1 0;0 1]; m.units={'1';'1';'rad';'rad';'rad'};
m.R=eye(5)*.01; m.R(1,5)=.002; m.R(5,1)=.002;
o=nav_odometry_measurement(a,b,m); assert(size(o.H,1)==5);
close(o.R(1,5),.002,1e-12);
m.R(1,5)=1; m.R(5,1)=1;
fails(@()nav_odometry_measurement(a,b,m),'libekf:Covariance');



end

function requirement_visual_ambiguity
K=[800 0 320;0 800 240;0 0 1];
[x,y]=meshgrid(linspace(-2,2,6),linspace(-1.5,1.5,5));
world=[x(:).';y(:).';10+2*sin((1:numel(x))*.7)];
training=mod((1:numel(x)).',2)==1; limits=struct('rotation_rms',1e-6,'plane_rms',1e-6);
rotation=nav_rotation([cos(.05);0;sin(.05);0]);
first=K*world; first=(first(1:2,:)./first(3,:)).';
second=K*rotation.'*world; second=(second(1:2,:)./second(3,:)).';
d=nav_visual_ambiguity(first,second,K,training,limits);
assert(d.geometry_ready&&strcmp(d.status,'AMBIGUOUS')&&d.rotation_explains&&d.plane_explains);
% A few bad held-out correspondences must not conceal a dominant rotation explanation.
outliers=second; indices=find(~training); outliers(indices(1),:)=outliers(indices(1),:)+[1000 -1000];
robust=limits; robust.min_support=.9; robust.point_limit=1e-5;
d=nav_visual_ambiguity(first,outliers,K,training,robust);
assert(d.rotation_rms>100&&d.rotation_explains&&strcmp(d.status,'AMBIGUOUS'));
close(d.rotation_support,14/15,1e-12);
outliers(indices(2),:)=outliers(indices(2),:)+[1000 -1000];
d=nav_visual_ambiguity(first,outliers,K,training,robust);
assert(~d.rotation_explains&&~d.plane_explains&&strcmp(d.status,'UNASSESSED'));
close(d.rotation_support,13/15,1e-12);
% Bad fitting matches must not conceal rotation; checking points never select a model.
contaminated=second; fit_indices=find(training);
contaminated(fit_indices(1),:)=contaminated(fit_indices(1),:)+[1000 -1000];
contaminated(indices(1),:)=contaminated(indices(1),:)+[1000 -1000];
sampled=robust; sampled.fit_samples=[2 4 6 8];
d=nav_visual_ambiguity(first,contaminated,K,training,sampled);
assert(d.rotation_explains&&strcmp(d.status,'AMBIGUOUS'));
close(d.rotation_support,14/15,1e-12);
bad_samples=sampled; bad_samples.fit_samples=[2 2 6 8];
fails(@()nav_visual_ambiguity(first,second,K,training,bad_samples),'libekf:Vision');
bad=robust; bad.min_support=1.1;
fails(@()nav_visual_ambiguity(first,second,K,training,bad),'libekf:Vision');

world(3,:)=10; first=K*world; first=(first(1:2,:)./first(3,:)).';
second=K*rotation.'*(world-[1;0;0]); second=(second(1:2,:)./second(3,:)).';
d=nav_visual_ambiguity(first,second,K,training,limits);
assert(d.geometry_ready&&strcmp(d.status,'AMBIGUOUS')&&d.plane_explains);
% A planar scene remains ambiguous despite gross matches in both subsets.
contaminated=second;
contaminated(fit_indices(1),:)=contaminated(fit_indices(1),:)+[1000 -1000];
contaminated(indices(1),:)=contaminated(indices(1),:)+[1000 -1000];
sampled.fit_samples=[2 5 9 13];
d=nav_visual_ambiguity(first,contaminated,K,training,sampled);
assert(d.plane_explains&&strcmp(d.status,'AMBIGUOUS'));
close(d.plane_support,14/15,1e-12);

world(3,:)=10+2*sin((1:numel(x))*.7); first=K*world; first=(first(1:2,:)./first(3,:)).';
second=K*rotation.'*(world-[1;0;0]); second=(second(1:2,:)./second(3,:)).';
d=nav_visual_ambiguity(first,second,K,training,limits);
assert(d.geometry_ready&&~d.rotation_explains&&~d.plane_explains&&strcmp(d.status,'UNASSESSED'));
first(:)=100; second(:)=110; d=nav_visual_ambiguity(first,second,K,training,limits);
assert(~d.geometry_ready&&strcmp(d.status,'UNASSESSED'));
fails(@()nav_visual_ambiguity(first,second,K,true(size(training)),limits),'libekf:Vision');
end

function requirement_source_comparison
a=struct('value',[2;0;0],'covariance',eye(3),'quantity','position','units',{{'m';'m';'m'}}, ...
    'frame_id','local','t',10,'quality','VALID','source','one');
b=a; b.value=zeros(3,1); b.source='two';
cfg=struct('max_age',1,'limit',3,'context','position_one_two_local');
[sample,result]=nav_compare_sources(a,b,.5*eye(3),10,cfg);
close(sample.value,[2;0;0],1e-12); close(result.nis,4,1e-12);
assert(strcmp(result.status,'DISAGREEMENT'));
a.value=a.value+[100;50;20]; b.value=b.value+[100;50;20];
[~,same]=nav_compare_sources(a,b,.5*eye(3),10,cfg); close(same.nis,4,1e-12);
[missing,result]=nav_compare_sources(a,b,.5*eye(3),12,cfg);
assert(isempty(missing)&&strcmp(result.status,'NOT_READY'));
b.quality='FAILED'; [missing,result]=nav_compare_sources(a,b,.5*eye(3),10,cfg);
assert(isempty(missing)&&strcmp(result.status,'NOT_READY'));
b.quality='VALID'; fails(@()nav_compare_sources(a,b,2*eye(3),10,cfg),'libekf:Covariance');
watch=struct('windows',[1 2 3],'max_gap',1,'max_age',1,'capacity',10, ...
    'instant_limit',1.5,'mean_limit',[10 10 10],'rms_limit',[10 10 10]);
[~,alarm]=nav_watch([],sample,10,watch); assert(strcmp(alarm.status,'ALARM'));
end

function requirement_navigation_outage
% Independent straight trajectory and physical measurements through one handler.
[s,e]=fixture; e.C_e_l=[0 0 1;1 0 0;0 1 0]; e.Omega=zeros(3,1);
e.gravity_origin=zeros(3,1); e.gravity_gradient=zeros(3); e.gravity_radius=1e5;
e.max_step=.1; e.state_noise_density=eye(21)*1e-10;
s.p=[5;2;99]; s.v=[1.2;.1;0]; s.q=[1;0;0;0]; s.ba=zeros(3,1); s.bg=zeros(3,1);
s.wind=zeros(3,1); s.bh=0;
P=diag([100*ones(3,1);ones(3,1);1e-4*ones(3,1);1e-8*ones(6,1);ones(3,1);1e-6;1;1]);
f=struct('state',s,'P',P,'used',{{}},'capacity',500);
sensor=struct('id',s.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
stations=[100 0 0;-100 100 0;0 -100 0]; start=s.t;
for k=1:20
    stop=start+k*.1;
    u=struct('model','inertial','source','imu','lineage',{{sprintf('imu_%d',k)}}, ...
        'quality','VALID','frame','BODY','accel_unit','m/s^2','gyro_unit','rad/s', ...
        't_start',f.state.t,'t_end',stop,'t_received',stop, ...
        'accel',zeros(3,1),'gyro',zeros(3,1),'noise_density',eye(6)*1e-6);
    % Decimal time representation may exceed .1 by rounding, not by physical duration.
    e.max_step=.100000000001;
    [f,info]=nav_assimilate(f,u,[],e,100); assert(info.accepted);
    truth=[k*.1;0;100];
    if k<=2 || k==10
        m=packet(f.state,'gnss_position'); m.lever=zeros(3,1); m.frame='SENSOR_LOCAL';
        m.z=truth; m.R=eye(3)*.01; m.lineage={sprintf('satellite_%d',k)};
        if k==10, m.z=m.z+[10000;10000;0]; end
        before=f; [f,info]=nav_assimilate(f,m,sensor,e,100);
        if k==10, assert(~info.accepted&&isequal(f,before)); else, assert(info.accepted); end
    end
    for j=1:3
        m=packet(f.state,'dme_range'); m.lever=zeros(3,1); m.frame='SENSOR_LOCAL';
        m.station=stations(j,:).'; m.z=norm(truth-m.station); m.R=.04;
        m.lineage={sprintf('dme_%d_%d',j,k)};
        [f,info]=nav_assimilate(f,m,sensor,e,100); assert(info.accepted);
    end
    m=packet(f.state,'doppler_velocity'); m.lever=zeros(3,1); m.frame='SENSOR_LOCAL';
    m.surface_velocity=zeros(3,1); m.z=[1;0;-1]/sqrt(2); m.R=eye(3)*.001;
    m.lineage={sprintf('doppler_%d',k)};
    [f,info]=nav_assimilate(f,m,sensor,e,100); assert(info.accepted);
    m=packet(f.state,'baro_height'); m.lever=zeros(3,1); m.frame='SENSOR_LOCAL';
    m.geoid_offset=0; m.z=hypot(e.ellipsoid_a+100,k*.1)-e.ellipsoid_a; m.R=.04;
    m.lineage={sprintf('baro_%d',k)};
    [f,info]=nav_assimilate(f,m,sensor,e,100); assert(info.accepted);
end
assert(norm(f.state.p-[2;0;100])<.1 && norm(f.state.v-[1;0;0])<.1);
assert(f.state.t==start+2 && all(isfinite(f.P(:))));
assert(~ismember('satellite_10',f.used));
end

function requirement_observable_quantities
% Two barometric readings do not separate height from their common bias.
[s,e]=fixture; m=packet(s,'baro_height'); m.lever=zeros(3,1);
o=nav_measurement(s,m,e); H=[o.H;o.H]; R=[4 1;1 4];
targets=zeros(3,21); targets(1,1:3)=o.H(1:3); targets(2,19)=1; targets(3,:)=o.H;
r=nav_observability(H,R,ones(21,1),targets,1e-10);
assert(r.rank==1&&isequal(r.identified,[false;false;true]));
assert(isinf(r.variance(1))&&isinf(r.variance(2))); close(r.variance(3),2.5,1e-10);
% A separately referenced geometric height separates height and barometric bias.
H(2,19)=0; R=diag([4 9]);
r=nav_observability(H,R,ones(21,1),targets,1e-10);
assert(r.rank==2&&all(r.identified)); close(r.variance,[9;13;4],1e-9);
% Air-relative velocity identifies the difference, not ground velocity or wind separately.
H=[eye(3) -eye(3)]; targets=[eye(3) zeros(3);zeros(3) eye(3);H];
r=nav_observability(H,eye(3),ones(6,1),targets,1e-10);
assert(r.rank==3&&isequal(r.identified,[false(6,1);true(3,1)]));
close(r.variance(7:9),ones(3,1),1e-12);
% Changing observation and state units leaves physical targets and uncertainty unchanged.
D=diag([1000 1000 1000 1 1 1]); W=diag([1e3 1e-2 10]);
q=nav_observability(W*H/D,W*W.',diag(D),targets/D,1e-10);
assert(isequal(q.identified,r.identified)); close(q.variance(7:9),r.variance(7:9),1e-12);
% Tiny but nonzero sensitivity must not be reported with a small uncertainty.
r=nav_observability(1e-8,1,1,1,1e-10);
assert(r.identified); close(r.variance,1e16,4);
r=nav_observability(zeros(2,3),eye(2),ones(3,1),eye(3),1e-10);
assert(r.rank==0&&~any(r.identified)&&all(isinf(r.variance)));
fails(@()nav_observability(eye(2),[1 2;2 1],ones(2,1),eye(2),1e-10),'libekf:Covariance');
% Publication cannot hide missing information behind a consistent watchdog.
[s,e]=fixture; m=packet(s,'baro_height'); m.lever=zeros(3,1); o=nav_measurement(s,m,e);
cfg=struct('windows',[1 3 6],'max_gap',1.1,'max_age',1,'capacity',20, ...
    'instant_limit',5,'mean_limit',[1 1 1],'rms_limit',[2 2 2]);
history=[];
for t=s.t-6:s.t
    [history,check]=nav_watch(history,struct('t',t,'context','height','value',0),t,cfg);
end
control=struct('state_time',s.t,'frame_id',s.frame_id,'now',s.t, ...
    'checks',struct('altitude',{{check}}));
control.geometry=struct('H',[o.H;o.H],'R',eye(2), ...
    'scales',ones(21,1),'relative_tolerance',1e-10,'max_variance',[1e-6;1e-6;10]);
out=nav_geodetic_output(s,eye(21)*1e-12,e,control);
assert(strcmp(out.quality.altitude,'UNOBSERVABLE')&&~out.observable(3));
assert(isinf(out.measurement_variance(3))); % A confident prior must not hide absent data.
control.geometry.H(2,19)=0; control.geometry.R=diag([1 100]);
out=nav_geodetic_output(s,eye(21),e,control);
assert(out.observable(3)&&strcmp(out.quality.altitude,'INSUFFICIENT_ACCURACY'));
control.geometry.max_variance(3)=101;
out=nav_geodetic_output(s,eye(21),e,control);
assert(strcmp(out.quality.altitude,'CONSISTENT')&&out.accuracy_sufficient(3));
assert(strcmp(out.quality.latitude,'UNOBSERVABLE'));
[~,alarm]=nav_watch([],struct('t',s.t,'context','height','value',6),s.t,cfg);
control.checks.altitude={alarm}; control.geometry.H(2,:)=o.H;
out=nav_geodetic_output(s,eye(21),e,control);
assert(strcmp(out.quality.altitude,'ALARM')&&~out.observable(3));

end

function models=source_models
models={'gnss_position','external_position','gnss_velocity','external_velocity', ...
    'gnss_position_velocity','external_position_velocity','dme_range','pseudorange', ...
    'range_rate','vor_azimuth','baro_height','airdata','radar_range','ils_angles', ...
    'gls_deviation','attitude','radar_footprint','doppler_velocity', ...
    'adf_bearing','station_range_azimuth','range_difference','time_difference','rsdn_time_difference','magnetic_vector','optical_landmarks','celestial_pixels','cooperative_range','static_pressure','impact_pressure','air_pressures'};
end

function requirement_all_source_timing
[s,e]=fixture; s.ba=zeros(3,1); s.bg=zeros(3,1);
e.Omega=zeros(3,1); e.gravity_origin=zeros(3,1); e.gravity_gradient=zeros(3);
e.gravity_radius=1e5; e.max_step=.125; e.state_noise_density=eye(21)*1e-8;
u=struct('model','inertial','source','imu','lineage',{{'timing_imu'}},'quality','VALID', ...
    'frame','BODY','accel_unit','m/s^2','gyro_unit','rad/s','t_start',s.t, ...
    't_end',s.t+.125,'t_received',s.t+.125,'accel',zeros(3,1),'gyro',zeros(3,1), ...
    'noise_density',eye(6)*1e-6);
f=struct('state',s,'P',eye(21),'used',{{}},'capacity',100);
sensor=struct('id',s.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
first=u; first.t_end=s.t+.0625;
[mid,midP]=nav_predict(s,f.P,first,e);
second=u; second.t_start=mid.t;
[plain,info]=nav_assimilate(f,u,[],e,1e8); assert(info.accepted);
required=fieldnames(u);
for j=1:numel(required)
    incomplete=rmfield(u,required{j});
    [unchanged,failed]=nav_assimilate(f,incomplete,[],e,1e8);
    assert(~failed.accepted&&isequal(unchanged,f));
end
bad=u; bad.t_start=[s.t s.t];
[unchanged,failed]=nav_assimilate(f,bad,[],e,1e8);
assert(~failed.accepted&&isequal(unchanged,f));
models=source_models;
for k=1:numel(models)
    m=packet(mid,models{k}); m.frame='SENSOR_LOCAL'; m.t_received=u.t_end+1;
    at_measurement=f; at_measurement.state=mid; at_measurement.P=midP;
    [reference,info]=nav_assimilate(at_measurement,m,sensor,e,1e8);
    assert(info.accepted,'libekf:Timing','Reference %s rejected: %s',models{k},info.reason);
    [reference.state,reference.P]=nav_predict(reference.state,reference.P,second,e);
    [replayed,info]=nav_replay(f,{u},{m},{sensor},e,1e8);
    assert(info.accepted&&info.observations{1}.accepted,'libekf:Timing', ...
        'Delayed %s rejected',models{k});
    close(replayed.state.p,reference.state.p,1e-10);
    close(replayed.state.v,reference.state.v,1e-10);
    close(nav_rotation(replayed.state.q),nav_rotation(reference.state.q),1e-10);
    close(replayed.P,reference.P,1e-10);
    assert(replayed.state.t==u.t_end&&numel(replayed.used)==2);
    % Each physical model obeys reception and age constraints before entering replay.
    admission=struct('now',m.t_received,'max_age',m.t_received-m.t_measure);
    [admitted,info]=nav_replay(f,{u},{m},{sensor},e,1e8,{admission});
    assert(info.accepted&&info.observations{1}.accepted&&isequal(admitted,replayed));
    admission.max_age=admission.max_age-.01;
    [rejected,info]=nav_replay(f,{u},{m},{sensor},e,1e8,{admission});
    assert(info.accepted&&strcmp(info.observations{1}.reason,'libekf:Stale')&&isequal(rejected,plain));
    admission.max_age=10; admission.now=m.t_received-.01;
    [rejected,info]=nav_replay(f,{u},{m},{sensor},e,1e8,{admission});
    assert(info.accepted&&strcmp(info.observations{1}.reason,'libekf:Time')&&isequal(rejected,plain));

    % Missing contract data reject only the packet; adapters must not throw raw field errors.
    required={'source','model','z','R','units','frame','frame_id','t_measure','t_received', ...
        'quality','lineage','correlation_group','noise_model'};
    for j=1:numel(required)
        incomplete=rmfield(m,required{j});
        [unchanged,failed]=nav_assimilate(at_measurement,incomplete,sensor,e,1e8);
        assert(~failed.accepted&&isequal(unchanged,at_measurement),'libekf:PacketContract', ...
            '%s accepted missing %s',models{k},required{j});
    end
    if any(strcmp(models{k},{'gnss_position','external_position','gnss_velocity','external_velocity', ...
            'gnss_position_velocity','external_position_velocity'}))
        malformed=m; malformed.z=m.z(1:end-1);
        [unchanged,failed]=nav_assimilate(at_measurement,malformed,sensor,e,1e8);
        assert(~failed.accepted&&strcmp(failed.reason,'libekf:Dimension')&&isequal(unchanged,at_measurement));
        malformed=m; malformed.R=m.R(1:end-1,1:end-1);
        [unchanged,failed]=nav_assimilate(at_measurement,malformed,sensor,e,1e8);
        assert(~failed.accepted&&strcmp(failed.reason,'libekf:Dimension')&&isequal(unchanged,at_measurement));
    end
    % A failed source must not affect propagation, regardless of its physical model.
    m.quality='FAILED';
    [replayed,info]=nav_replay(f,{u},{m},{sensor},e,1e8);
    assert(info.accepted&&~info.observations{1}.accepted&&isequal(replayed,plain));
end
end

function status=run_tests
% Проверки требований к поведению; один файл, без внешних настроек и отчётов.
base=fileparts(mfilename('fullpath')); addpath(fullfile(base,'function'));
cases={@requirement_source_models,@requirement_coordinates_and_terrain, ...
    @requirement_numerical_invariance,@requirement_broad_radar,@requirement_common_correction, ...
    @requirement_source_provenance,@requirement_inertial_motion,@requirement_watch_intervals, ...
    @requirement_terrain_sequence,@requirement_all_source_fusion,@requirement_relative_camera_motion,@requirement_joint_correction,@requirement_visual_ambiguity,@requirement_correlated_measurement,@requirement_source_comparison,@requirement_navigation_outage,@requirement_interval_observations,@requirement_observable_quantities,@requirement_all_source_timing,@requirement_timed_rays};
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
    'Все источники учитывают своё время при пересчёте истории', ...
    'Времена строк и погрешности лучей камеры входят в общую поправку'};
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
    % История получает именно принятую модель, а не повторную линеаризацию после поправки.
    model=nav_observe_source(f.state,m,sensor0,e);
    assert(isequal(ra.observation,model)&&isequal(size(ra.error_reset),[21 21])&&ra.state_times==f.state.t);
    [same_state,same_P,correction]=nav_correct(f.state,f.P,ra.observation,1e8);
    close(same_state.p,a.state.p,1e-12); close(same_P,a.P,1e-12);
    close(correction.error_reset,ra.error_reset,1e-12);
    % Общая аддитивная ошибка проходит через физическую модель каждого источника.
    if strcmp(model.noise_model,'independent')&&~isfield(model,'peer_covariance')
        count=numel(model.residual); labels=arrayfun(@(j)sprintf('shared_%d',j),1:count,'UniformOutput',false);
        with_noise=f; with_noise.future_noise=struct('ids',{labels},'units',{model.units}, ...
            'mean',zeros(count,1),'covariance',.1*model.R,'cross_cov',zeros(21,count), ...
            'state_times',s.t,'frame_id',s.frame_id);
        dependent=m; dependent.future_noise_independent=true;
        dependent.shared_noise=struct('ids',{labels},'units',{model.units},'map',eye(count), ...
            'residual_units',{model.units},'residual_frame_id',sensor0.id,'lineage',{{}});
        [actual,accepted]=nav_assimilate(with_noise,dependent,sensor0,e,1e8);
        assert(accepted.accepted,'libekf:AllSources','Shared error in %s rejected: %s',models{k},accepted.reason);
        full=model; full.R=1.1*model.R;
        for field={'input_noise_cov','state_input_cov','noise_transform'}
            if isfield(full,field{1}), full=rmfield(full,field{1}); end
        end
        [expected_state,expected_P]=nav_correct(f.state,f.P,full,1e8);
        close(actual.state.p,expected_state.p,1e-10); close(actual.P,expected_P,1e-10);
        common_history=struct('steps',{{struct('kind',m.model,'time',s.t,'data',accepted)}}, ...
            'initial_times',s.t,'frame_id',s.frame_id,'noise_independent',true, ...
            'initial_noise',with_noise.future_noise);
        observed=nav_observability_sequence(common_history,ones(21,1),eye(21),1e-10);
        expected_geometry=nav_observability(model.H,1.1*model.R,ones(21,1),accepted.error_reset,1e-10);
        assert(isequal(observed.identified,expected_geometry.identified));
        known=observed.identified;
        if any(known), close(observed.variance(known),expected_geometry.variance(known),1e-8); end

    end

    % Определимость из принятой модели совпадает с отдельным расчётом без априорной P.
    history=struct('steps',{{struct('kind',m.model,'time',s.t,'data',ra)}}, ...
        'initial_times',s.t,'frame_id',s.frame_id,'noise_independent',true);
    if isfield(ra.observation,'peer_covariance')||~strcmp(ra.observation.noise_model,'independent')
        fails(@()nav_observability_sequence(history,ones(21,1),eye(21),1e-10),'libekf:Correlation');
    else
        geometry=nav_observability_sequence(history,ones(21,1),eye(21),1e-10);
        expected_geometry=nav_observability(model.H,model.R,ones(21,1),ra.error_reset,1e-10);
        assert(isequal(geometry.identified,expected_geometry.identified));
        known=geometry.identified;
        if any(known), close(geometry.variance(known),expected_geometry.variance(known),1e-8); end
        assert(all(isinf(geometry.variance(~known))));
    end

    bad=m; bad.quality='FAILED';
    [unchanged,failed]=nav_assimilate(f,bad,sensor0,e,1e8);
    assert(~failed.accepted&&isequal(unchanged,f)&&~isfield(failed,'observation'));

    control=struct('state_time',a.state.t,'frame_id',a.state.frame_id,'now',a.state.t+.5, ...
        'checks',struct(),'source_history',a.source_history,'sources_complete',true);
    published=nav_geodetic_output(a.state,a.P,e,control);
    assert(published.age==.5&&published.sources_complete&&numel(published.sources)==1);
    assert(strcmp(published.sources.id,m.source)&&strcmp(published.sources.models{1},m.model)&& ...
        published.sources.age==.5&&~published.sources.excluded);
    control.excluded_sources={m.source}; control=rmfield(control,'sources_complete');
    published=nav_geodetic_output(a.state,a.P,e,control);
    assert(published.sources.excluded&&~published.sources_complete);
    % Прибор не исчезает из происхождения оценки за именем вычислителя.
    derived=m; derived.source='processor';
    derived.input_sources={'physical_sensor','shared_clock','physical_sensor','processor'};
    [tracked,accepted]=nav_assimilate(f,derived,sensor0,e,1e8);
    assert(accepted.accepted&&isequal({tracked.source_history.id}, ...
        {'physical_sensor','processor','shared_clock'}));
    assert(all([tracked.source_history.first_epoch]==s.t)&&all([tracked.source_history.last_epoch]==s.t));
    closed=nav_checkpoint(tracked);
    assert(isequal(closed.source_history,tracked.source_history));
    control.source_history=tracked.source_history; control.excluded_sources={'physical_sensor'};
    published=nav_geodetic_output(tracked.state,tracked.P,e,control);
    assert(isequal({published.sources([published.sources.excluded]).id},{'physical_sensor'}));
    failed_packet=derived; failed_packet.quality='FAILED';
    [unchanged,failed]=nav_assimilate(f,failed_packet,sensor0,e,1e8);
    assert(~failed.accepted&&isequal(unchanged,f));
    future=f; future.source_history=a.source_history; future.source_history.last_epoch=s.t+1;
    [unchanged,failed]=nav_assimilate(future,m,sensor0,e,1e8);
    assert(~failed.accepted&&strcmp(failed.reason,'libekf:SourceHistory')&&isequal(unchanged,future));

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
% LEKF-001/003: все модели одновременно в одном состоянии.
% Нулевые невязки фиксируют точку линеаризации; независимые ошибки
% позволяют получить эталон по сумме информационных матриц.
packets=cell(size(models)); information=eye(21);
for k=1:numel(models)
    m=packet(s,models{k}); m.frame='SENSOR_LOCAL';
    observation=nav_observe_source(s,m,sensor0,e);
    if strcmp(m.model,'attitude'), m.z=s.q;
    elseif strcmp(m.model,'rsdn_arrival_times'), m.z=[0;observation.predicted];
    else, m.z=observation.predicted; end
    observation=nav_observe_source(s,m,sensor0,e);
    assert(norm(observation.residual)<1e-10);
    information=information+observation.H.'*(observation.R\observation.H);
    packets{k}=m;
end
expected=information\eye(21);
for order={1:numel(models),numel(models):-1:1}
    combined=f;
    for k=order{1}
        [combined,accepted]=nav_assimilate(combined,packets{k},sensor0,e,1e8);
        assert(accepted.accepted,'libekf:AllSources','Combined %s: %s',models{k},accepted.reason);
    end
    close(combined.P,expected,2e-9);
    close(combined.state.p,s.p,1e-10); close(combined.state.v,s.v,1e-10);
    close(nav_rotation(combined.state.q),nav_rotation(s.q),1e-10);
    assert(numel(combined.used)==numel(models)&&numel(combined.source_history)==numel(models));
end
% Неисправный прибор не влияет на результат остальных систем.
excluded=7; mixed=f; healthy=f;
for k=1:numel(models)
    m=packets{k};
    if k==excluded, m.quality='FAILED'; end
    [mixed,accepted]=nav_assimilate(mixed,m,sensor0,e,1e8);
    assert(accepted.accepted==(k~=excluded));
    if k~=excluded
        [healthy,accepted]=nav_assimilate(healthy,packets{k},sensor0,e,1e8);
        assert(accepted.accepted);
    end
end
assert(isequal(mixed,healthy));
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
% LEKF-024: поиск по области не приписывает плоскому рельефу горизонтальную привязку.
search=struct('bounds',[-10 10;0 0;-1 1],'step',[10;1;1],'max_candidates',9);
r=nav_terrain_sequence(paths(1,:),series,map,limits,search);
assert(strcmp(r.status,'AMBIGUOUS')&&numel(r.candidates)==3&&~r.search.continuous_verified);
surviving=r.offsets(r.candidates,:);
close(sort(surviving(:,1)),[-10;0;10],1e-12);
close(surviving(:,2:3),zeros(3,2),1e-12);
close(r.predicted,10+r.offsets(:,3)*ones(1,3),1e-12);
% Границы включены даже при шаге, не делящем ширину области.
search.bounds=[-10 10;0 0;0 0]; search.step=[7;1;1];
r=nav_terrain_sequence(paths(1,:),series,map,limits,search);
assert(numel(r.candidates)==4&&min(r.offsets(:,1))==-10&&max(r.offsets(:,1))==10);
assert(r.search.step(1)<=7);
% Нехватка памяти не превращается в неполный поиск с положительным результатом.
search.max_candidates=3;
fails(@()nav_terrain_sequence(paths(1,:),series,map,limits,search),'libekf:TerrainCapacity');
search.max_candidates=9; search.step(1)=0;
fails(@()nav_terrain_sequence(paths(1,:),series,map,limits,search),'libekf:TerrainSearch');
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
    'rsdn_arrival_times',{'geometry_at_epoch','stations','propagation_speed','time_correction','reference_station', ...
        'arrival_time_origin','synchronization_valid','ambiguity_resolved'}; ...
    'doppler_velocity',{'beams_body','surface_valid','beam_geometry_valid','surface_velocity'}; ...
    'doppler_frequency',{'beams_body','surface_valid','beam_geometry_valid','surface_velocity', ...
        'wavelength','frequency_sign','frequency_correction','monostatic','calibration_valid'}; ...
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
% Общие отсчёты допустимы только как условное наблюдение с действительно новой информацией.
% x~N(0,1), z1=x+n1=2, z2=x+n2=4, Var(n1)=Var(n2)=1.
% После z1: E(n1|z1)=1, Var(n1|z1)=1/2, Cov(e,n1|z1)=-1/2.
old=packet(s,'gnss_position'); old.frame='SENSOR_LOCAL'; old.lever=zeros(3,1);
old.z=s.p+[2;0;0]; old.lineage={'conditional_old'};
[conditioned,info]=nav_assimilate(f,old,sensor,e,100); assert(info.accepted);
conditioned.capacity=2;
independent=old; independent.z=s.p+[4;0;0]; independent.lineage={'conditional_new'};
[expected,info]=nav_assimilate(conditioned,independent,sensor,e,100); assert(info.accepted);
combined=old; combined.z=s.p+[3;0;0]; combined.R=.375*eye(3);
combined.lineage={'conditional_old','conditional_new'}; combined.noise_model='joint_required';
combined.state_noise_frame_id=s.frame_id; combined.state_noise_cov=zeros(21,3);
combined.state_noise_cov(1:3,:)=-.25*eye(3);
combined.reuse=struct('lineage',{{'conditional_old'}},'history',{conditioned.used}, ...
    'state_times',conditioned.state.t,'noise_mean',[.5;0;0]);
[actual,info]=nav_assimilate(conditioned,combined,sensor,e,100); assert(info.accepted);
uncentred=nav_observe_source(conditioned.state,combined,sensor,e);
close(info.observation.residual,uncentred.residual-combined.reuse.noise_mean,1e-12);

close(actual.state.p,expected.state.p,1e-12); close(actual.P,expected.P,1e-12);
close(actual.state.p,s.p+[2;0;0],1e-12); close(diag(actual.P(1:3,1:3)),ones(3,1)/3,1e-12);
assert(numel(actual.used)==2);
% Среднее условного шума поворачивается как ошибка, без переноса начала координат.
rotation=[0 -1 0;1 0 0;0 0 1]; offset=[3;4;5];
local=struct('C_e_l',rotation,'origin_ecef',offset,'frame_id','conditional_sensor');
original=combined; original.frame='ECEF'; native=nav_packet_from_ecef(original,local);
native.frame='SENSOR_LOCAL';
rotated_sensor=struct('id',local.frame_id,'C_e_s',e.C_e_l*rotation, ...
    'origin_ecef',e.origin_ecef+e.C_e_l*offset);
[rotated_result,info]=nav_assimilate(conditioned,native,rotated_sensor,e,100);
assert(info.accepted); close(rotated_result.state.p,actual.state.p,1e-12); close(rotated_result.P,actual.P,1e-12);

[unchanged,info]=nav_assimilate(actual,combined,sensor,e,100);
assert(~info.accepted&&isequal(unchanged,actual));
for field={'lineage','history','state_times','noise_mean'}
    bad=combined; bad.reuse=rmfield(bad.reuse,field{1});
    [unchanged,info]=nav_assimilate(conditioned,bad,sensor,e,100);
    assert(~info.accepted&&strcmp(info.reason,'libekf:Reuse')&&isequal(unchanged,conditioned));
end
% Та же условная поправка сохраняет связи с прошлой опорой камеры.
paired=f; paired.reference=s; paired.reference.t=s.t-.1;
paired.jointP=[eye(21) .5*eye(21);.5*eye(21) eye(21)];
[paired,info]=nav_assimilate(paired,old,sensor,e,100); assert(info.accepted);
[paired_expected,info]=nav_assimilate(paired,independent,sensor,e,100); assert(info.accepted);
paired_packet=combined;
paired_packet.state_noise_cov=-.5*paired.jointP(:,22:24);
paired_packet.R=.25*(paired.P(1:3,1:3)+eye(3));
paired_packet.reuse.noise_mean=.5*(old.z-paired.state.p);
paired_packet.reuse.history=paired.used;
paired_packet.reuse.state_times=[paired.reference.t paired.state.t];
[paired_actual,info]=nav_assimilate(paired,paired_packet,sensor,e,100); assert(info.accepted);
close(paired_actual.reference.p,paired_expected.reference.p,1e-12);
close(paired_actual.state.p,paired_expected.state.p,1e-12);
close(paired_actual.jointP,paired_expected.jointP,1e-12);
bad=paired_packet; bad.reuse.state_times(1)=bad.reuse.state_times(1)-.1;
[unchanged,info]=nav_assimilate(paired,bad,sensor,e,100);
assert(~info.accepted&&strcmp(info.reason,'libekf:Reuse')&&isequal(unchanged,paired));
bad=combined; bad.reuse.history={'different_history'};
[unchanged,info]=nav_assimilate(conditioned,bad,sensor,e,100);
assert(~info.accepted&&strcmp(info.reason,'libekf:Reuse')&&isequal(unchanged,conditioned));
bad=combined; bad.noise_model='independent';
[unchanged,info]=nav_assimilate(conditioned,bad,sensor,e,100);
assert(~info.accepted&&strcmp(info.reason,'libekf:Reuse')&&isequal(unchanged,conditioned));

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

% LEKF-027: общая опорная станция связывает разности; ошибка общих часов сокращается.
arrival=packet(s,'rsdn_arrival_times'); arrival.z=[.001;.001003;.001005];
arrival.R=diag([4 9 16])*1e-12;
observed=nav_measurement(s,arrival,e);
close(observed.R,[13 4;4 20]*1e-12,1e-24);
close(observed.residual+observed.predicted,[3;5]*1e-6,1e-18);
common=arrival; common.z=arrival.z+.125; common.R=arrival.R+25e-12*ones(3);
shifted=nav_measurement(s,common,e);
close(shifted.residual,observed.residual,1e-16); close(shifted.R,observed.R,1e-24);
for field={'synchronization_valid','ambiguity_resolved'}
    bad=arrival; bad.(field{1})=false;
    fails(@()nav_measurement(s,bad,e),'libekf:RSDN');
end
% Полный пересчёт исходной взаимной ковариации даёт ту же общую поправку.
joint_arrival=arrival; joint_arrival.noise_model='joint_required';
joint_arrival.state_noise_frame_id=s.frame_id;
joint_arrival.state_noise_cov=zeros(21,3); joint_arrival.state_noise_cov(1,:)=[1 2 3]*1e-8;
prepared=joint_arrival; prepared.model='rsdn_time_difference';
prepared.z=arrival.z(2:end)-arrival.z(1); prepared.units={'s';'s'};
prepared.R=[13 4;4 20]*1e-12; prepared.state_noise_cov=zeros(21,2);
prepared.state_noise_cov(1,:)=[1 2]*1e-8;
[raw_state,raw_P,raw_info]=nav_correct(s,eye(21),nav_measurement(s,joint_arrival,e),1e8);
[prepared_state,prepared_P,prepared_info]=nav_correct(s,eye(21),nav_measurement(s,prepared,e),1e8);
assert(raw_info.accepted&&prepared_info.accepted);
close(raw_state.p,prepared_state.p,1e-12); close(raw_P,prepared_P,1e-12);
% Недопустимая связь с общими часами не скрывается после взятия разностей.
bad=arrival; bad.noise_model='joint_required'; bad.state_noise_frame_id=s.frame_id;
bad.state_noise_cov=zeros(21,3); bad.state_noise_cov(1,:)=1e-3;
bad_observation=nav_measurement(s,bad,e);
close(bad_observation.state_noise_cov,zeros(21,2),1e-15);
fails(@()nav_correct(s,eye(21),bad_observation,1e8),'libekf:Correlation');
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
% Частота отражённого сигнала: знак, удвоенный путь и эквивалентность единиц.
s.q=[1;0;0;0]; s.v=[3;4;5];
frequency=packet(s,'doppler_frequency'); frequency.lever=zeros(3,1);
frequency.beams_body=eye(3); frequency.surface_velocity=zeros(3,1);
frequency.wavelength=.5; frequency.frequency_correction=[1;2;3];
frequency.z=[13;18;23]; frequency.R=[4 1 0;1 9 2;0 2 16];
obs=nav_measurement(s,frequency,e);
close(obs.predicted,[13;18;23],1e-12); close(obs.H(:,4:6),4*eye(3),1e-12);
frequency.frequency_sign=-1;
close(nav_measurement(s,frequency,e).predicted,[-11;-14;-17],1e-12);
frequency.frequency_sign=1;
frequency.z=frequency.z+[2;-3;1];
velocity=frequency; velocity.model='doppler_velocity';
velocity.z=(frequency.z-frequency.frequency_correction)/4;
velocity.R=frequency.R/16; velocity.units=repmat({'m/s'},3,1);
frequency.noise_model='joint_required'; velocity.noise_model='joint_required';
velocity.state_noise_cov=zeros(21,3); velocity.state_noise_cov(4:6,:)=.01*eye(3);
frequency.state_noise_cov=4*velocity.state_noise_cov;
frequency.state_noise_frame_id=s.frame_id; velocity.state_noise_frame_id=s.frame_id;
[sa,Pa,ia]=nav_correct(s,eye(21),nav_measurement(s,frequency,e),100);
[sb,Pb,ib]=nav_correct(s,eye(21),nav_measurement(s,velocity,e),100);
assert(ia.accepted&&ib.accepted); close(sa.v,sb.v,1e-12); close(Pa,Pb,1e-12);
for field={'monostatic','calibration_valid'}
    bad=frequency; bad.(field{1})=false;
    fails(@()nav_measurement(s,bad,e),'libekf:Doppler');
end
bad=frequency; bad.wavelength=0;
fails(@()nav_measurement(s,bad,e),'libekf:Doppler');
bad=frequency; bad.frequency_sign=0;
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
    case {'rsdn_time_difference','rsdn_arrival_times'}
        m.stations=s.p.'+[300 400 200;-400 100 200];
        m.reference_station=s.p+[-500;100;300]; m.propagation_speed=299792458;
        m.time_correction=[2;1]/m.propagation_speed; m.geometry_at_epoch=true;
        m.z=[0;0]; m.R=[1 .5;.5 1]/m.propagation_speed^2; m.units={'s';'s'};
        if strcmp(model,'rsdn_arrival_times')
            m.z=zeros(3,1); m.R=eye(3)*.5/m.propagation_speed^2; m.units={'s';'s';'s'};
            m.arrival_time_origin=s.t; m.synchronization_valid=true; m.ambiguity_resolved=true;
        end
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
    case {'doppler_velocity','doppler_frequency'}
        m.beams_body=[1 0 -1;0 1 -1;-1 0 -1]/sqrt(2);
        m.surface_velocity=[1;0;0]; m.surface_valid=true; m.beam_geometry_valid=true;
        m.z=[-3;-1;-7]/sqrt(2); m.R=eye(3); m.units=repmat({'m/s'},3,1);
        if strcmp(model,'doppler_frequency')
            m.wavelength=.03; m.frequency_sign=1; m.frequency_correction=[1;2;3];
            m.monostatic=true; m.calibration_valid=true;
            m.z=(2/m.wavelength)*m.z+m.frequency_correction;
            m.R=(2/m.wavelength)^2*m.R; m.units=repmat({'Hz'},3,1);
        end
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
% Полные окна проверки только для тревоги не подтверждают выходную величину.
alarm_history=[];
for t=s.t-6:s.t
    [alarm_history,alarm_check]=nav_watch(alarm_history, ...
        struct('t',t,'context','bounded_height','value',0,'alarm_only',true),t,c);
end
assert(all(alarm_check.ready)&&strcmp(alarm_check.status,'NOT_READY'));
control.checks.altitude={alarm_check};
out=nav_geodetic_output(s,eye(21),env,control);
assert(strcmp(out.quality.altitude,'NOT_READY'));
control.checks.altitude={good,alarm_check};
out=nav_geodetic_output(s,eye(21),env,control);
assert(strcmp(out.quality.altitude,'NOT_READY'));
forged=alarm_check; forged.status='CONSISTENT'; control.checks.altitude={forged};
fails(@()nav_geodetic_output(s,eye(21),env,control),'libekf:WatchOutput');
forged=alarm_check; forged.alarm_only=[true false]; control.checks.altitude={forged};
fails(@()nav_geodetic_output(s,eye(21),env,control),'libekf:WatchOutput');
[~,alarm_check]=nav_watch([],struct('t',s.t,'context','bounded_alarm', ...
    'value',6,'alarm_only',true),s.t,c);
control.checks.altitude={good,alarm_check};
out=nav_geodetic_output(s,eye(21),env,control);
assert(strcmp(out.quality.altitude,'ALARM'));
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
[state,environment]=fixture;
% Аналитическая совместная поправка: собственная координата, опора и ошибка прибора.
P=eye(21); H=zeros(1,21); H(1)=1;
U=zeros(21,3); U(1,1)=.2; Uf=zeros(21,1); Uf(1)=.1;
future=struct('ids',{{'calibration'}},'units',{{'m'}},'mean',0,'covariance',3, ...
    'cross_cov',Uf,'measurement_cross_cov',-.3,'state_times',state.t, ...
    'frame_id',state.frame_id,'peer_cross_cov',[.3;0;0],'peer_time',state.t,'peer_ids',{{'peer'}});
o=struct('H',H,'R',3,'residual',3.6,'noise_model','joint_required', ...
    'state_noise_cov',-U(:,1),'frame_id',state.frame_id,'t_measure',state.t, ...
    'units',{{'m'}},'peer_covariance',diag([2 1 1]),'peer_cross_cov',U, ...
    'peer_jacobian',[-1 0 0],'peer_noise_cov',1,'peer_position',zeros(3,1), ...
    'peer_time',state.t,'peer_ids',{{'peer'}},'future_noise',future);
[updated,Q,info]=nav_correct(state,P,o,100);
assert(info.accepted); close(updated.p(1)-state.p(1),.8,1e-12);
close(info.peer_update.position,[-1.8;0;0],1e-12);
close(info.future_noise_update.mean,-.2,1e-12);
% Независимый эталон условной ковариации трёх скалярных величин.
prior=[1 .2 .1;.2 2 .3;.1 .3 3]; response=[.8;-1.8;-.2];
expected=prior-response*response.'/3.6;
actual=[Q(1,1),info.peer_update.cross_cov(1,1),info.future_noise_update.cross_cov(1); ...
    info.peer_update.cross_cov(1,1),info.peer_update.covariance(1,1),info.future_noise_update.peer_cross_cov(1); ...
    info.future_noise_update.cross_cov(1),info.future_noise_update.peer_cross_cov(1),info.future_noise_update.covariance];
close(actual,expected,1e-12);
% Следующая независимая дальность использует обновлённую связь, а не исходную.
next=o; next.residual=0; next.peer_position=info.peer_update.position;
next.peer_covariance=info.peer_update.covariance; next.peer_cross_cov=info.peer_update.cross_cov;
next.R=1+next.peer_covariance(1,1); next.state_noise_cov=-next.peer_cross_cov(:,1);
next.future_noise=info.future_noise_update;
next.future_noise.measurement_cross_cov=-next.future_noise.peer_cross_cov(1,:).';
[~,Q2,second]=nav_correct(updated,Q,next,100);
h=[1 -1 0]; expected2=expected-(expected*h.')*(h*expected)/(1+h*expected*h.');
actual2=[Q2(1,1),second.peer_update.cross_cov(1,1),second.future_noise_update.cross_cov(1); ...
    second.peer_update.cross_cov(1,1),second.peer_update.covariance(1,1),second.future_noise_update.peer_cross_cov(1); ...
    second.future_noise_update.cross_cov(1),second.future_noise_update.peer_cross_cov(1),second.future_noise_update.covariance];
close(actual2,expected2,1e-12);
scaled=o; scaled.future_noise.units={'micrometre'};
scaled.future_noise.covariance=3e12; scaled.future_noise.cross_cov=Uf*1e6;
scaled.future_noise.peer_cross_cov=future.peer_cross_cov*1e6;
scaled.future_noise.measurement_cross_cov=future.measurement_cross_cov*1e6;
[~,scaledP,scaledInfo]=nav_correct(state,P,scaled,100);
close(scaledP,Q,1e-12);
close(scaledInfo.future_noise_update.peer_cross_cov/1e6,info.future_noise_update.peer_cross_cov,1e-12);
bad=o; bad.future_noise=rmfield(future,'peer_cross_cov');
fails(@()nav_correct(state,P,bad,100),'libekf:Correlation');
bad=o; bad.future_noise.peer_ids={'wrong_peer'};
fails(@()nav_correct(state,P,bad,100),'libekf:Correlation');
% Невозможная связь в направлении, не измеряемом дальностью, тоже отклоняется.
bad=o; bad.future_noise.peer_cross_cov(2)=10;
fails(@()nav_correct(state,P,bad,100),'libekf:Correlation');
bad=o; bad.future_noise.peer_time=state.t-1;
fails(@()nav_correct(state,P,bad,100),'libekf:Correlation');
[same,sameP,rejected]=nav_correct(state,P,o,.01);
assert(~rejected.accepted&&isequal(same,state)&&isequal(sameP,P)&&~isfield(rejected,'future_noise_update'));
stored=info.future_noise_update;
unchanged=nav_predict_noise(stored,Q,eye(21),zeros(21),state.t,zeros(21,1));
close(unchanged.peer_cross_cov,stored.peer_cross_cov,0);
fails(@()nav_predict_noise(stored,Q,eye(21),zeros(21),state.t+1,zeros(21,1)), ...
    'libekf:Correlation');
% Условные ошибки повторной калибровки готовятся из сохранённого набора автоматически.
[s,e]=fixture; sensor=struct('id',s.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
initial=struct('state',s,'P',eye(21),'used',{{}},'capacity',20);
initial.future_noise=struct('ids',{{'bias_x','bias_y','bias_z'}},'units',{{'m','m','m'}}, ...
    'mean',zeros(3,1),'covariance',eye(3),'cross_cov',zeros(21,3),'state_times',s.t,'frame_id',s.frame_id);
one=packet(s,'gnss_position'); one.frame='SENSOR_LOCAL'; one.lever=zeros(3,1);
one.z=s.p+[2;0;0]; one.R=eye(3); one.lineage={'calibration','sample_1'};
one.future_noise_independent=true;
one.shared_noise=struct('ids',{initial.future_noise.ids},'units',{initial.future_noise.units}, ...
    'map',eye(3),'residual_units',{{'m','m','m'}},'residual_frame_id',sensor.id,'lineage',{{'calibration'}});
[first_state,check]=nav_assimilate(initial,one,sensor,e,100); assert(check.accepted);
first_record=check;
two=one; two.z=s.p+[4;0;0]; two.lineage={'calibration','sample_2'};
[second_state,check]=nav_assimilate(first_state,two,sensor,e,100); assert(check.accepted);
% R двух измерений = [2 1;1 2]; начальная дисперсия положения = 1.
close(second_state.state.p,s.p+[6/5;0;0],1e-12);
close(second_state.P(1:3,1:3),(3/5)*eye(3),1e-12);
close(second_state.future_noise.mean,[6/5;0;0],1e-12);
close(second_state.future_noise.covariance,(3/5)*eye(3),1e-12);
close(second_state.future_noise.cross_cov,[-.4*eye(3);zeros(18,3)],1e-12);
close(check.observation.conditional_noise_mean,[2/3;0;0],1e-12);
% Проверка определимости восстанавливает исходную общую ошибку, а не считает условные R независимыми.
history=struct('steps',{{struct('kind',one.model,'time',s.t,'data',first_record), ...
    struct('kind',two.model,'time',s.t,'data',check)}}, ...
    'initial_times',s.t,'frame_id',s.frame_id,'noise_independent',true,'initial_noise',initial.future_noise);
target=[eye(3) zeros(3,18)]; geometry=nav_observability_sequence(history,ones(21,1),target,1e-10);
assert(all(geometry.identified)); close(geometry.variance,1.5*ones(3,1),1e-12);
missing=rmfield(history,'initial_noise');
fails(@()nav_observability_sequence(missing,ones(21,1),target,1e-10),'libekf:Correlation');
bad_history=history; bad_history.initial_noise.cross_cov(1,1)=.1;
fails(@()nav_observability_sequence(bad_history,ones(21,1),target,1e-10),'libekf:Correlation');
bad_history=history; bad_history.initial_noise.covariance(3,3)=-.5;
fails(@()nav_observability_sequence(bad_history,ones(21,1),target,1e-10),'libekf:Correlation');
correlated_step=struct('kind','propagation','data',struct('transition',eye(21), ...
    'process_covariance',eye(21),'t_start',s.t,'t_end',s.t+1,'frame_id',s.frame_id, ...
    'future_noise_process_cross_cov',[.1*eye(3);zeros(18,3)]));
bad_history=history; bad_history.steps{end+1}=correlated_step;
fails(@()nav_observability_sequence(bad_history,ones(21,1),target,1e-10),'libekf:Correlation');

assert(numel(second_state.used)==3&&sum(strcmp(second_state.used,'calibration'))==1);
[same,rejected]=nav_assimilate(second_state,two,sensor,e,100);
assert(~rejected.accepted&&isequal(same,second_state));
permuted=first_state; order=[3 1 2];
permuted.future_noise.ids=permuted.future_noise.ids(order);
permuted.future_noise.units=permuted.future_noise.units(order);
permuted.future_noise.mean=permuted.future_noise.mean(order);
permuted.future_noise.covariance=permuted.future_noise.covariance(order,order);
permuted.future_noise.cross_cov=permuted.future_noise.cross_cov(:,order);
[same,check]=nav_assimilate(permuted,two,sensor,e,100); assert(check.accepted);
close(same.state.p,second_state.state.p,1e-12); close(same.P,second_state.P,1e-12);
bad=two; bad.shared_noise.units{1}='s';
[same,rejected]=nav_assimilate(first_state,bad,sensor,e,100);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:Correlation')&&isequal(same,first_state));
bad=two; bad.shared_noise.residual_frame_id='another_sensor_axes';
[same,rejected]=nav_assimilate(first_state,bad,sensor,e,100);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:Correlation')&&isequal(same,first_state));
bad=two; bad.shared_noise.lineage={};
[same,rejected]=nav_assimilate(first_state,bad,sensor,e,100);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:Reuse')&&isequal(same,first_state));
% Смена местного начала и осей сохраняет следующую поправку с общей ошибкой.
to=e; to.frame_id='shared_error_other_axes';
to.C_e_l=e.C_e_l*nav_rotation([.9;.1;-.2;.3]/norm([.9;.1;-.2;.3]));
to.origin_ecef=e.origin_ecef+[300;-200;50];
for paired_case=[false true]
    seed=first_state; if paired_case, seed=nav_save_reference(seed); end
    [expected,ok]=nav_assimilate(seed,two,sensor,e,100); assert(ok.accepted);
    [converted,~,J]=nav_change_frame(seed,[],e,to);
    close(converted.future_noise.cross_cov,J*seed.future_noise.cross_cov,1e-12);
    assert(isequal(converted.future_noise.mean,seed.future_noise.mean)&& ...
        isequal(converted.future_noise.covariance,seed.future_noise.covariance));
    assert(isequal(converted.used,seed.used)&&isequal(converted.source_history,seed.source_history));
    [actual,ok]=nav_assimilate(converted,two,sensor,to,100); assert(ok.accepted);
    expected=nav_change_frame(expected,[],e,to);
    close(actual.state.p,expected.state.p,1e-9); close(actual.P,expected.P,1e-11);
    close(actual.future_noise.mean,expected.future_noise.mean,1e-9);
    close(actual.future_noise.covariance,expected.future_noise.covariance,1e-11);
    close(actual.future_noise.cross_cov,expected.future_noise.cross_cov,1e-11);
    if paired_case
        close(actual.reference.p,expected.reference.p,1e-9); close(actual.jointP,expected.jointP,1e-11);
    end
    restored=nav_change_frame(converted,[],to,e);
    close(restored.state.p,seed.state.p,1e-9); close(restored.P,seed.P,1e-11);
    close(restored.future_noise.cross_cov,seed.future_noise.cross_cov,1e-11);
end
with_reference=rmfield(first_state,'future_noise');
with_reference.reference_bank=struct('ids',{{'peer'}},'valid',true,'position',s.p+[100;20;0], ...
    'covariance',eye(3),'cross_cov',zeros(21,3),'t',s.t,'frame_id',s.frame_id);
converted=nav_change_frame(with_reference,[],e,to);
[~,~,~,expected_reference]=nav_change_frame(with_reference.state,with_reference.P,e,to,with_reference.reference_bank);
assert(isequal(converted.reference_bank,expected_reference));
bad=first_state; bad.future_noise.state_times=s.t-1;
fails(@()nav_change_frame(bad,[],e,to),'libekf:Correlation');

% При поступлении старого пакета условная модель строится по пересчитанной истории.
motion=e; motion.gravity_origin=zeros(3,1); motion.gravity_gradient=zeros(3);
motion.gravity_radius=1e6; motion.max_step=.2; motion.state_noise_density=.01*eye(21);
imu=struct('model','inertial','source','imu','lineage',{{'shared_replay_imu'}},'quality','VALID', ...
    'frame','BODY','accel_unit','m/s^2','gyro_unit','rad/s','t_start',s.t,'t_end',s.t+.1, ...
    't_received',s.t+.1,'accel',zeros(3,1),'gyro',zeros(3,1),'noise_density',.01*eye(6), ...
    'future_noise_independent',true);
one.t_measure=s.t+.025; one.t_received=imu.t_end; one.z=s.p+s.v*.025+[2;0;0];
two.t_measure=s.t+.075; two.t_received=imu.t_end; two.z=s.p+s.v*.075+[4;0;0];
[timely,check]=nav_assimilate_interval(initial,imu,{one,two},{sensor,sensor},motion,100);
assert(check.accepted&&all(cellfun(@(z)z.accepted,check.observations)));
[late,check]=nav_replay(initial,{imu},{two,one},{sensor,sensor},motion,100);
assert(check.accepted&&all(cellfun(@(z)z.accepted,check.observations)));
close(late.state.p,timely.state.p,1e-11); close(late.P,timely.P,1e-11);
close(late.future_noise.mean,timely.future_noise.mean,1e-11);
close(late.future_noise.covariance,timely.future_noise.covariance,1e-11);
close(late.future_noise.cross_cov,timely.future_noise.cross_cov,1e-11);

% Принятое измерение меняет распределение общей ошибки следующего измерения.
[s,e]=fixture; sensor=struct('id',s.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
f=struct('state',s,'P',eye(21),'used',{{}},'capacity',10);
first=packet(s,'gnss_position'); first.frame='SENSOR_LOCAL'; first.lever=zeros(3,1);
first.z=s.p+[2;0;0]; first.R=eye(3); first.lineage={'shared_calibration','first_sample'};
future=struct('ids',{{'next_x','next_y','next_z'}},'units',{{'m','m','m'}}, ...
    'mean',zeros(3,1),'covariance',eye(3),'cross_cov',zeros(21,3), ...
    'measurement_cross_cov',.5*eye(3),'state_times',s.t,'frame_id',s.frame_id);
first.future_noise=future;
[conditioned,info]=nav_assimilate(f,first,sensor,e,100); assert(info.accepted);
future=info.future_noise_update;
close(future.mean,[.5;0;0],1e-12); close(future.covariance,.875*eye(3),1e-12);
close(future.cross_cov,[-.25*eye(3);zeros(18,3)],1e-12);
assert(~isfield(future,'measurement_cross_cov'));
second=rmfield(first,'future_noise'); second.z=s.p+[4;0;0];
second.lineage={'shared_calibration','second_sample'}; second.R=future.covariance;
second.state_noise_cov=future.cross_cov; second.state_noise_frame_id=s.frame_id;
second.noise_model='joint_required';
second.future_noise_independent=false; second.future_noise_cross_cov=future.covariance;
second.reuse=struct('lineage',{{'shared_calibration'}},'history',{conditioned.used}, ...
    'state_times',future.state_times,'noise_mean',future.mean);
[combined,accepted]=nav_assimilate(conditioned,second,sensor,e,100); assert(accepted.accepted);
% Совместное аналитическое решение: R=[1 .5;.5 1], P0=1.
close(combined.state.p,s.p+[12/7;0;0],1e-12);
close(combined.P(1:3,1:3),(3/7)*eye(3),1e-12);
assert(isequal(conditioned.future_noise,future)&&isequal(combined.future_noise,accepted.future_noise_update));
% Связь общей ошибки переносится через движение, промежуточную поправку и смену опоры.
motion=e; motion.gravity_origin=[0;0;-9.8]; motion.gravity_gradient=zeros(3);
motion.gravity_radius=1e6; motion.max_step=.2; motion.state_noise_density=.01*eye(21);
imu=struct('model','inertial','source','imu','lineage',{{'future_noise_imu'}},'quality','VALID', ...
    'frame','BODY','accel_unit','m/s^2','gyro_unit','rad/s','t_start',s.t,'t_end',s.t+.1, ...
    't_received',s.t+.1,'accel',zeros(3,1),'gyro',zeros(3,1),'noise_density',.01*eye(6), ...
    'future_noise_independent',true);
[advanced,ok]=nav_assimilate(conditioned,imu,[],motion,100); assert(ok.accepted);
close(advanced.future_noise.cross_cov,ok.propagation.transition*future.cross_cov,1e-12);
close(advanced.future_noise.covariance,future.covariance,0);
close(advanced.future_noise.mean,future.mean,0); assert(advanced.future_noise.state_times==imu.t_end);
half=imu; half.t_end=s.t+.05;
[middle,Pmiddle,Tmiddle]=nav_predict(conditioned.state,conditioned.P,half,motion);
observation=packet(middle,'external_position'); observation.lever=zeros(3,1); observation.frame='SENSOR_LOCAL';
observation.future_noise_independent=true; observation.t_received=imu.t_end;
observation.z=middle.p;
[through,ok]=nav_assimilate_interval(conditioned,imu,{observation},{sensor},motion,100);
assert(ok.accepted&&ok.observations{1}.accepted);
H=[eye(3) zeros(3,18)]; U=Tmiddle*future.cross_cov; S=H*Pmiddle*H.'+eye(3);
K=Pmiddle*H.'/S; Up=U-K*H*U; Vp=future.covariance-(H*U).'/S*(H*U);
last=imu; last.t_start=middle.t;
[~,~,Tlast]=nav_predict(middle,Pmiddle-K*H*Pmiddle,last,motion);
close(through.future_noise.cross_cov,Tlast*Up,1e-11);
close(through.future_noise.covariance,Vp,1e-11);
[replayed,ok]=nav_replay(conditioned,{imu},{observation},{sensor},motion,100);
assert(ok.accepted); close(replayed.future_noise.cross_cov,through.future_noise.cross_cov,1e-12);
saved=nav_save_reference(through);
close(saved.future_noise.cross_cov,[through.future_noise.cross_cov;through.future_noise.cross_cov],1e-12);
assert(isequal(saved.future_noise.state_times,[imu.t_end;imu.t_end]));
% Для целого шага допускается явно заданная связь нового процессного шума с общим.
correlated=imu; correlated.future_noise_independent=false;
correlated.future_noise_process_cross_cov=[.001*eye(3);zeros(18,3)];
[moved,ok]=nav_assimilate(conditioned,correlated,[],motion,100); assert(ok.accepted);
close(moved.future_noise.cross_cov,advanced.future_noise.cross_cov+correlated.future_noise_process_cross_cov,1e-12);
[same,rejected]=nav_assimilate_interval(conditioned,correlated,{observation},{sensor},motion,100);
assert(~rejected.accepted&&isequal(same,conditioned));
bad=imu; bad.future_noise_process_cross_cov=correlated.future_noise_process_cross_cov;
[same,rejected]=nav_assimilate(conditioned,bad,[],motion,100);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:Correlation')&&isequal(same,conditioned));
bad=rmfield(imu,'future_noise_independent');
[same,rejected]=nav_assimilate(conditioned,bad,[],motion,100);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:Correlation')&&isequal(same,conditioned));
bad=correlated; bad.future_noise_process_cross_cov=100*bad.future_noise_process_cross_cov;
[same,rejected]=nav_assimilate(conditioned,bad,[],motion,100);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:Correlation')&&isequal(same,conditioned));

% Сохранённая опора сохраняет ту же общую ошибку, а не независимую копию.
paired=nav_save_reference(f); with_pair=first;
with_pair.future_noise.cross_cov=zeros(42,3); with_pair.future_noise.state_times=[s.t;s.t];
[~,accepted]=nav_assimilate(paired,with_pair,sensor,e,100); assert(accepted.accepted);
close(accepted.future_noise_update.mean,future.mean,1e-12);
close(accepted.future_noise_update.cross_cov,[future.cross_cov;future.cross_cov],1e-12);
% Единицы будущей ошибки не влияют на полученную физическую поправку.
scaled=first; scaled.future_noise.mean=1e6*scaled.future_noise.mean;
scaled.future_noise.covariance=1e12*scaled.future_noise.covariance;
scaled.future_noise.cross_cov=1e6*scaled.future_noise.cross_cov;
scaled.future_noise.measurement_cross_cov=1e6*scaled.future_noise.measurement_cross_cov;
scaled.future_noise.units={'um','um','um'};
[same,accepted]=nav_assimilate(f,scaled,sensor,e,100); assert(accepted.accepted);
close(same.P,conditioned.P,1e-12); close(accepted.future_noise_update.mean/1e6,future.mean,1e-12);
close(accepted.future_noise_update.covariance/1e12,future.covariance,1e-12);
bad=first; bad.future_noise.measurement_cross_cov=2*eye(3);
[same,rejected]=nav_assimilate(f,bad,sensor,e,100);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:Correlation')&&isequal(same,f));
[same,rejected]=nav_assimilate(f,first,sensor,e,.01);
assert(~rejected.accepted&&isequal(same,f)&&~isfield(rejected,'future_noise_update'));
bad=first; bad.future_noise.state_times=s.t-1;
[same,rejected]=nav_assimilate(f,bad,sensor,e,100);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:Correlation')&&isequal(same,f));

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

% Several references acquire mutual dependence through the shared navigation state.
bank=struct('position',[s.p-[10;0;0];s.p+[20;0;0]], ...
    'covariance',diag([9 1 1 4 1 1]),'cross_cov',zeros(21,6), ...
    'frame_id',s.frame_id,'t',s.t,'noise_independent',true);
bank.cross_cov(1,[1 4])=[1 2]; bankP=eye(21); bankP(1,1)=4;
[~,~,bank_info]=nav_correct(s,bankP,position_obs,100,bank);
assert(bank_info.accepted);
close(bank_info.peer_update.position,bank.position+[.2;0;0;.4;0;0],1e-12);
close(bank_info.peer_update.covariance([1 4],[1 4]),[8.8 -.4;-.4 3.2],1e-12);
close(bank_info.peer_update.cross_cov(1,[1 4]),[.2 .4],1e-12);
motion=struct('A',eye(6),'shift',zeros(6,1),'Q',zeros(6), ...
    't_start',s.t,'t_end',s.t+1,'frame_id',s.frame_id,'noise_independent',true);
next_bank=nav_predict_reference(bank_info.peer_update,eye(21),motion);
close(next_bank.covariance,bank_info.peer_update.covariance,0);
close(next_bank.cross_cov,bank_info.peer_update.cross_cov,0);
% Reordering references must reorder all associated covariance blocks together.
order=[4:6 1:3]; swapped=bank;
swapped.position=bank.position(order); swapped.covariance=bank.covariance(order,order);
swapped.cross_cov=bank.cross_cov(:,order);
[~,~,swapped_info]=nav_correct(s,bankP,position_obs,100,swapped);
close(swapped_info.peer_update.position,bank_info.peer_update.position(order),1e-12);
close(swapped_info.peer_update.covariance,bank_info.peer_update.covariance(order,order),1e-12);

% A range selects its reference by ID and updates the complete correlated set.
bank.ids={'left','right'}; bank.valid=[true true];
m=packet(s,'cooperative_range'); m.lever=zeros(3,1); m.z=11; m.R=1;
m=rmfield(m,{'station','peer_covariance','peer_cross_cov','peer_time','peer_valid','peer_cross_frame_id'});
m.reference_bank=bank; m.peer_id='left'; m.noise_model='joint_required';
o=nav_measurement(s,m,e); [ranged,rangeP,range_info]=nav_correct(s,bankP,o,100);
assert(range_info.accepted&&isequal(range_info.peer_update.ids,bank.ids));
close(ranged.p,s.p+[.25;0;0],1e-12); close(rangeP(1,1),13/4,1e-12);
close(range_info.peer_update.position,bank.position+[-2/3;0;0;1/6;0;0],1e-12);
close(range_info.peer_update.covariance(1,4),4/3,1e-12);
bad=m; bad.peer_id='missing'; fails(@()nav_measurement(s,bad,e),'libekf:Peer');
bad=m; bad.reference_bank.ids={'left','left'}; fails(@()nav_measurement(s,bad,e),'libekf:Peer');
bad=m; bad.reference_bank.valid(1)=false; fails(@()nav_measurement(s,bad,e),'libekf:Peer');
permuted=m; permuted.reference_bank.ids=bank.ids([2 1]);
permuted.reference_bank.position=bank.position(order);
permuted.reference_bank.covariance=bank.covariance(order,order);
permuted.reference_bank.cross_cov=bank.cross_cov(:,order);
o=nav_measurement(s,permuted,e); [same,sameP,permuted_info]=nav_correct(s,bankP,o,100);
close(same.p,ranged.p,1e-12); close(sameP,rangeP,1e-12);
close(permuted_info.peer_update.position,range_info.peer_update.position(order),1e-12);

% A stored set is updated atomically along with the navigation state.
sensor=struct('id',s.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
filter=struct('state',s,'P',bankP,'used',{{}},'capacity',10,'reference_bank',bank);
managed=m; managed=rmfield(managed,'reference_bank'); managed.frame='SENSOR_LOCAL';
[after_range,accepted]=nav_assimilate(filter,managed,sensor,e,100);
assert(accepted.accepted); close(after_range.state.p,ranged.p,1e-12);
close(after_range.reference_bank.position,range_info.peer_update.position,1e-12);
close(after_range.reference_bank.covariance,range_info.peer_update.covariance,1e-12);
assert(isequal(after_range.reference_bank.valid,bank.valid));
position=packet(after_range.state,'gnss_position'); position.lever=zeros(3,1);
position.frame='SENSOR_LOCAL'; position.reference_noise_independent=true;
prepared=after_range.reference_bank; prepared.noise_independent=true;
o=nav_observe_source(after_range.state,position,sensor,e);
[expected,expectedP,expected_info]=nav_correct(after_range.state,after_range.P,o,100,prepared);
[after_position,accepted]=nav_assimilate(after_range,position,sensor,e,100);
assert(accepted.accepted); close(after_position.state.p,expected.p,1e-12); close(after_position.P,expectedP,1e-12);
close(after_position.reference_bank.cross_cov,expected_info.peer_update.cross_cov,1e-12);
position.reference_noise_independent=false;
[same,rejected]=nav_assimilate(after_range,position,sensor,e,100);
assert(~rejected.accepted&&isequal(same,after_range));
cloned=nav_save_reference(after_range);
close(cloned.reference_bank.cross_cov,[after_range.reference_bank.cross_cov;after_range.reference_bank.cross_cov],0);

% Inertial propagation carries the stored set; an absent motion model is rejected.
motion_env=e; motion_env.Omega=zeros(3,1); motion_env.gravity_origin=zeros(3,1);
motion_env.gravity_gradient=zeros(3); motion_env.gravity_radius=1e5;
motion_env.max_step=.1; motion_env.state_noise_density=eye(21)*1e-8;
imu=struct('model','inertial','source','imu','lineage',{{'bank_motion'}},'quality','VALID', ...
    'frame','BODY','accel_unit','m/s^2','gyro_unit','rad/s','t_start',s.t,'t_end',s.t+.1, ...
    't_received',s.t+.1,'accel',zeros(3,1),'gyro',zeros(3,1),'noise_density',eye(6)*1e-8);
[same,rejected]=nav_assimilate(after_range,imu,[],motion_env,100);
assert(~rejected.accepted&&isequal(same,after_range));
imu.reference_motion=struct('A',eye(6),'shift',zeros(6,1),'Q',zeros(6), ...
    't_start',s.t,'t_end',s.t+.1,'frame_id',s.frame_id,'noise_independent',true);
% An impossible stored dependence must not pass an inertial-only step.
broken_filter=after_range; broken_filter.reference_bank.cross_cov(1,1)=1e6;
[unchanged,bad_result]=nav_assimilate(broken_filter,imu,[],motion_env,100);
assert(~bad_result.accepted&&strcmp(bad_result.reason,'libekf:Correlation')&&isequal(unchanged,broken_filter));
broken_filter=after_range; broken_filter.reference_bank.ids={'left','left'};
[unchanged,bad_result]=nav_assimilate(broken_filter,imu,[],motion_env,100);
assert(~bad_result.accepted&&strcmp(bad_result.reason,'libekf:Peer')&&isequal(unchanged,broken_filter));
[~,~,transition]=nav_predict(after_range.state,after_range.P,imu,motion_env);
[moved,accepted]=nav_assimilate(after_range,imu,[],motion_env,100);
assert(accepted.accepted&&moved.reference_bank.t==moved.state.t);
close(moved.reference_bank.position,after_range.reference_bank.position,0);
close(moved.reference_bank.covariance,after_range.reference_bank.covariance,0);
close(moved.reference_bank.cross_cov,transition*after_range.reference_bank.cross_cov,1e-12);
[same,rejected]=nav_assimilate_interval(after_range,imu,{position},{sensor},motion_env,100);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:ReferenceInterval')&&isequal(same,after_range));

% Continuous reference dynamics compose across arbitrary observation epochs.
continuous=struct('F',zeros(6),'drift',[2;0;0;0;0;0],'noise_density',eye(6)*.5, ...
    't_start',s.t,'t_end',s.t+2,'frame_id',s.frame_id,'noise_independent',true);
full=nav_predict_reference(bank,eye(21),continuous);
close(full.position,bank.position+[4;0;0;0;0;0],1e-12);
close(full.covariance,bank.covariance+eye(6),1e-12);
continuous.F(1,1)=.1;
full=nav_predict_reference(bank,eye(21),continuous);
part=continuous; part.t_end=s.t+.7; left=nav_predict_reference(bank,eye(21),part);
part.t_start=part.t_end; part.t_end=continuous.t_end;
split=nav_predict_reference(left,eye(21),part);
close(split.position,full.position,1e-10); close(split.covariance,full.covariance,1e-10);
close(split.cross_cov,full.cross_cov,1e-10);
imu.reference_motion=continuous; imu.reference_motion.t_end=imu.t_end;
first=imu; first.t_end=s.t+.05; first.reference_motion.t_end=first.t_end; first.lineage={'manual_first'};
[mid,ok]=nav_assimilate(after_range,first,[],motion_env,100); assert(ok.accepted);
measurement=packet(mid.state,'gnss_position'); measurement.frame='SENSOR_LOCAL';
measurement.lever=zeros(3,1); measurement.reference_noise_independent=true; measurement.t_received=imu.t_end;
[manual,ok]=nav_assimilate(mid,measurement,sensor,motion_env,100); assert(ok.accepted);
last=imu; last.t_start=first.t_end; last.reference_motion.t_start=last.t_start; last.lineage={'manual_last'};
[manual,ok]=nav_assimilate(manual,last,[],motion_env,100); assert(ok.accepted);
[automatic,ok]=nav_assimilate_interval(after_range,imu,{measurement},{sensor},motion_env,100);
assert(ok.accepted&&ok.observations{1}.accepted);
close(automatic.state.p,manual.state.p,1e-10); close(automatic.P,manual.P,1e-10);
close(automatic.reference_bank.position,manual.reference_bank.position,1e-10);
close(automatic.reference_bank.covariance,manual.reference_bank.covariance,1e-10);
close(automatic.reference_bank.cross_cov,manual.reference_bank.cross_cov,1e-10);

% A late position observation rebuilds all later reference errors from the snapshot.
motion_env.max_step=.100000000001;
late=measurement; late.z=late.z+[.1;0;0]; late.t_received=s.t+.2;
save=struct('model','save_reference','t_measure',imu.t_end);
[timely,ok]=nav_assimilate_interval(after_range,imu,{late,save},{sensor,[]},motion_env,100);
assert(ok.accepted&&all(cellfun(@(x)x.accepted,ok.observations)));
next_imu=imu; next_imu.t_start=imu.t_end; next_imu.t_end=s.t+.2;
next_imu.t_received=next_imu.t_end; next_imu.lineage={'bank_second_interval'};
next_imu.reference_motion.t_start=next_imu.t_start; next_imu.reference_motion.t_end=next_imu.t_end;
first_part=next_imu; first_part.t_end=s.t+.15; first_part.reference_motion.t_end=first_part.t_end;
first_part.lineage={'manual_second_first'};
[at_range,ok]=nav_assimilate(timely,first_part,[],motion_env,100); assert(ok.accepted);
later=managed; later.t_measure=first_part.t_end; later.t_received=next_imu.t_end;
later.lineage={'bank_later_range'};
later.z=norm(at_range.state.p-at_range.reference_bank.position(1:3))+.1;
[expected,ok]=nav_assimilate(at_range,later,sensor,motion_env,100); assert(ok.accepted);
last_part=next_imu; last_part.t_start=first_part.t_end;
last_part.reference_motion.t_start=last_part.t_start; last_part.lineage={'manual_second_last'};
[expected,ok]=nav_assimilate(expected,last_part,[],motion_env,100); assert(ok.accepted);
[rebuilt,ok]=nav_replay(after_range,{imu,next_imu},{later,late,save},{sensor,sensor,[]},motion_env,100);
assert(ok.accepted&&all(cellfun(@(x)x.accepted,ok.observations)));
close(rebuilt.state.p,expected.state.p,1e-10); close(rebuilt.jointP,expected.jointP,1e-10);
close(rebuilt.reference_bank.position,expected.reference_bank.position,1e-10);
close(rebuilt.reference_bank.covariance,expected.reference_bank.covariance,1e-10);
close(rebuilt.reference_bank.cross_cov,expected.reference_bank.cross_cov,1e-10);
[again,ok]=nav_replay(after_range,{imu,next_imu},{later,late,save},{sensor,sensor,[]},motion_env,100);
assert(ok.accepted&&isequal(again,rebuilt));
[without,ok]=nav_replay(after_range,{imu,next_imu},{later,save},{sensor,[]},motion_env,100);
assert(ok.accepted&&norm(without.reference_bank.position-rebuilt.reference_bank.position)>1e-8);

% Changing local axes preserves ranges and complete reference dependencies.
destination=e; destination.frame_id='bank_rotated';
destination.C_e_l=e.C_e_l*[0 -1 0;1 0 0;0 0 1];
destination.origin_ecef=e.origin_ecef+e.C_e_l*[30;-20;10];
[changed,changedP,~,changed_bank]=nav_change_frame(s,bankP,e,destination,bank);
[restored,restoredP,~,restored_bank]=nav_change_frame(changed,changedP,destination,e,changed_bank);
close(restored.p,s.p,1e-8); close(restoredP,bankP,1e-12);
close(restored_bank.position,bank.position,1e-8);
close(restored_bank.covariance,bank.covariance,1e-12);
close(restored_bank.cross_cov,bank.cross_cov,1e-12);
original_packet=m; original_packet.reference_bank=bank;
changed_packet=original_packet; changed_packet.frame_id=destination.frame_id;
changed_packet.reference_bank=changed_bank;
original_observation=nav_measurement(s,original_packet,e);
changed_observation=nav_measurement(changed,changed_packet,destination);
close(changed_observation.predicted,original_observation.predicted,1e-8);
[original_state,originalP,original_info]=nav_correct(s,bankP,original_observation,100);
[changed_state,newP,new_info]=nav_correct(changed,changedP,changed_observation,100);
[expected_state,expectedP,~,expected_bank]=nav_change_frame(original_state,originalP,e,destination,original_info.peer_update);
close(changed_state.p,expected_state.p,1e-8); close(newP,expectedP,1e-8);
close(new_info.peer_update.cross_cov,expected_bank.cross_cov,1e-8);
close(new_info.peer_update.covariance,expected_bank.covariance,1e-8);
[paired_states,pairedP,~,paired_bank]=nav_change_frame( ...
    {cloned.reference,cloned.state},cloned.jointP,e,destination,cloned.reference_bank);
assert(numel(paired_states)==2&&isequal(size(pairedP),[42 42]));
close(paired_bank.cross_cov(1:21,:),paired_bank.cross_cov(22:42,:),1e-12);

% Propagation and frame conversion commute, including the affine origin term.
[~,~,~,rotated_bank,rotated_motion]=nav_change_frame(s,bankP,e,destination,bank,continuous);
predicted_old=nav_predict_reference(bank,eye(21),continuous);
predicted_new=nav_predict_reference(rotated_bank,eye(21),rotated_motion);
[~,~,~,predicted_expected]=nav_change_frame(s,bankP,e,destination,predicted_old);
close(predicted_new.position,predicted_expected.position,1e-8);
close(predicted_new.covariance,predicted_expected.covariance,1e-10);
close(predicted_new.cross_cov,predicted_expected.cross_cov,1e-10);
discrete=struct('A',diag([1.1 1 1 .9 1 1]),'shift',[2;0;0;1;0;0], ...
    'Q',eye(6)*.3,'t_start',s.t,'t_end',s.t+2,'frame_id',s.frame_id,'noise_independent',true);
[~,~,~,rotated_bank,rotated_motion]=nav_change_frame(s,bankP,e,destination,bank,discrete);
predicted_old=nav_predict_reference(bank,eye(21),discrete);
predicted_new=nav_predict_reference(rotated_bank,eye(21),rotated_motion);
[~,~,~,predicted_expected]=nav_change_frame(s,bankP,e,destination,predicted_old);
close(predicted_new.position,predicted_expected.position,1e-8);
close(predicted_new.covariance,predicted_expected.covariance,1e-10);
close(predicted_new.cross_cov,predicted_expected.cross_cov,1e-10);
bad_model=discrete; bad_model.frame_id='wrong';
fails(@()nav_change_frame(s,bankP,e,destination,bank,bad_model),'libekf:Frame');

% Known measurement dependence with both navigation and references is retained.
dependent=packet(s,'gnss_position'); dependent.lever=zeros(3,1); dependent.z=s.p+[1;0;0];
dependent.R=eye(3); dependent.noise_model='joint_required';
dependent.state_noise_cov=zeros(21,3); dependent.state_noise_cov(1,1)=.5;
dependent.state_noise_frame_id=s.frame_id;
retained=bank; retained.noise_independent=false; retained.noise_cross_cov=zeros(6,3);
retained.noise_cross_cov(1,1)=.25;
o=nav_measurement(s,dependent,e); [updated,updatedP,info]=nav_correct(s,bankP,o,100,retained);
assert(info.accepted); close(updated.p,s.p+[.75;0;0],1e-12); close(updatedP(1,1),5/8,1e-12);
close(info.peer_update.position(1),bank.position(1)+5/24,1e-12);
close(info.peer_update.covariance(1,1),839/96,1e-12);
close(info.peer_update.cross_cov(1,1),1/16,1e-12);
impossible=retained; impossible.noise_cross_cov(1,1)=1e6;
fails(@()nav_correct(s,bankP,o,100,impossible),'libekf:Correlation');
filter=struct('state',s,'P',bankP,'used',{{}},'capacity',10,'reference_bank',bank);
dependent.frame='SENSOR_LOCAL'; dependent.reference_noise_independent=false;
dependent.reference_noise_cross_cov=retained.noise_cross_cov;
[through,accepted]=nav_assimilate(filter,dependent,sensor,e,100);
assert(accepted.accepted); close(through.P,updatedP,1e-12);
close(through.reference_bank.cross_cov,info.peer_update.cross_cov,1e-12);
assert(~isfield(through.reference_bank,'noise_cross_cov'));

% Changing measurement axes/scales must transform both noise cross-covariances.
W=[0 -2 0;3 0 0;0 0 .5]; transformed=o;
transformed.residual=W*o.residual; transformed.H=W*o.H; transformed.R=W*o.R*W.';
transformed.state_noise_cov=o.state_noise_cov*W.';
transformed_reference=retained; transformed_reference.noise_cross_cov=retained.noise_cross_cov*W.';
[scaled,scaledP,scaled_info]=nav_correct(s,bankP,transformed,100,transformed_reference);
assert(scaled_info.accepted); close(scaled.p,updated.p,1e-12); close(scaledP,updatedP,1e-12);
close(scaled_info.peer_update.position,info.peer_update.position,1e-12);
close(scaled_info.peer_update.covariance,info.peer_update.covariance,1e-12);
close(scaled_info.peer_update.cross_cov,info.peer_update.cross_cov,1e-12);
% Check the real sensor adapter, including its translated origin.
rotation=[0 -1 0;1 0 0;0 0 1]; offset=[3;-5;2];
rotated_sensor=struct('id','dependent_sensor','origin_ecef',e.origin_ecef+e.C_e_l*offset, ...
    'C_e_s',e.C_e_l*rotation);
rotated_packet=dependent; rotated_packet.frame_id=rotated_sensor.id;
rotated_packet.z=rotation.'*(dependent.z-offset);
rotated_packet.R=rotation.'*dependent.R*rotation;
rotated_packet.state_noise_cov=dependent.state_noise_cov*rotation;
rotated_packet.reference_noise_cross_cov=dependent.reference_noise_cross_cov*rotation;
[rotated_result,accepted]=nav_assimilate(filter,rotated_packet,rotated_sensor,e,100);
assert(accepted.accepted); close(rotated_result.state.p,through.state.p,1e-8);
close(rotated_result.P,through.P,1e-8);
close(rotated_result.reference_bank.position,through.reference_bank.position,1e-8);
close(rotated_result.reference_bank.cross_cov,through.reference_bank.cross_cov,1e-8);

% An independence flag cannot hide a supplied nonzero reference/noise covariance.
contradictory=dependent; contradictory.noise_model='independent';
contradictory.state_noise_cov=zeros(21,3); contradictory.reference_noise_independent=true;
[unchanged,rejected]=nav_assimilate(filter,contradictory,sensor,e,100);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:Correlation')&&isequal(unchanged,filter));
contradictory.reference_noise_cross_cov=zeros(6,3);
[~,accepted]=nav_assimilate(filter,contradictory,sensor,e,100); assert(accepted.accepted);

% Health reports have their own time and provenance; they are not measurements.
health=struct('model','reference_health','source','reference_monitor','lineage',{{'left_failed'}}, ...
    'quality','VALID','peer_id','left','valid',false,'frame_id',s.frame_id, ...
    't_measure',s.t,'t_received',s.t);
[disabled,accepted]=nav_assimilate(after_range,health,[],motion_env,100);
assert(accepted.accepted&&~disabled.reference_bank.valid(1));
assert(isequal(disabled.state,after_range.state)&&isequal(disabled.P,after_range.P));
close(disabled.reference_bank.covariance,after_range.reference_bank.covariance,0);
bad_health=health; bad_health.peer_id='unknown';
[unchanged,rejected]=nav_assimilate(after_range,bad_health,[],motion_env,100);
assert(~rejected.accepted&&isequal(unchanged,after_range));
[unchanged,rejected]=nav_assimilate(disabled,health,[],motion_env,100);
assert(~rejected.accepted&&isequal(unchanged,disabled));
% Недостоверный, просроченный или ещё не полученный отчёт не меняет решение.
bad_health=health; bad_health.quality='INVALID';
[unchanged,rejected]=nav_assimilate(after_range,bad_health,[],motion_env,100);
assert(~rejected.accepted&&isequal(unchanged,after_range));
bad_health=health; bad_health.t_received=s.t+1;
[unchanged,rejected]=nav_assimilate(after_range,bad_health,[],motion_env,100, ...
    struct('now',s.t,'max_age',2));
assert(~rejected.accepted&&isequal(unchanged,after_range));
[unchanged,rejected]=nav_assimilate(after_range,health,[],motion_env,100, ...
    struct('now',s.t+2,'max_age',1));
assert(~rejected.accepted&&isequal(unchanged,after_range));
[boundary,accepted]=nav_assimilate(after_range,health,[],motion_env,100, ...
    struct('now',s.t+1,'max_age',1));
assert(accepted.accepted&&isequal(boundary,disabled));
health.t_measure=s.t+.125; health.t_received=s.t+.2;
[disabled_history,history_check]=nav_replay(after_range,{imu,next_imu},{later,health},{sensor,[]},motion_env,100);
assert(history_check.accepted&&~history_check.observations{1}.accepted&&history_check.observations{2}.accepted);
[health_only,history_check]=nav_replay(after_range,{imu,next_imu},{health},{[]},motion_env,100);
assert(history_check.accepted&&isequal(disabled_history,health_only));
restored_health=health; restored_health.valid=true; restored_health.t_measure=s.t+.14;
restored_health.lineage={'left_restored'};
[restored_history,history_check]=nav_replay(after_range,{imu,next_imu}, ...
    {later,restored_health,health},{sensor,[],[]},motion_env,100);
assert(history_check.accepted&&all(cellfun(@(x)x.accepted,history_check.observations))&& ...
    restored_history.reference_bank.valid(1));

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
[saved,change]=nav_save_reference(f); assert(isequal(saved.state,f.state)&&isequal(saved.P,f.P));
close(change.transition*f.jointP*change.transition.',saved.jointP,1e-12);
assert(isequal(change.input_times,[f.reference.t;f.state.t]));
single=rmfield(f,{'reference','jointP'});
[cloned,change]=nav_save_reference(single);
close(change.transition*single.P*change.transition.',cloned.jointP,1e-12);
assert(isequal(size(change.transition),[42 21])&&change.input_times==single.state.t);
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

% После закрытия истории первый кадр на границе нельзя переиспользовать в новой паре.
closed=nav_checkpoint(saved); dated_imu=u; dated_imu.lineage_times=u.t_end;
[continued,info]=nav_assimilate(closed,dated_imu,[],e,100); assert(info.accepted);
dated_pair=m; dated_pair.lineage={'boundary_frame','after_boundary_frame'};
dated_pair.lineage_times=[m.t_start m.t_end];
[unchanged,info]=nav_assimilate(continued,dated_pair,[],e,100);
assert(~info.accepted&&strcmp(info.reason,'libekf:RetiredHistory')&&isequal(unchanged,continued));
% Пара, оба кадра которой после границы, проходит обычную совместную поправку.
fresh=nav_save_reference(continued);
dated_imu.t_start=fresh.state.t; dated_imu.t_end=fresh.state.t+.1;
dated_imu.t_received=dated_imu.t_end; dated_imu.lineage={'after_boundary_imu'};
dated_imu.lineage_times=dated_imu.t_end;
[fresh,info]=nav_assimilate(fresh,dated_imu,[],e,100); assert(info.accepted);
dated_pair.t_start=fresh.reference.t; dated_pair.t_end=fresh.state.t;
dated_pair.t_received=fresh.state.t; dated_pair.lineage={'fresh_first','fresh_second'};
dated_pair.lineage_times=[dated_pair.t_start dated_pair.t_end];
dated_pair.z=fresh.state.p-fresh.reference.p;
[updated,info]=nav_assimilate(fresh,dated_pair,[],e,100);
assert(info.accepted&&updated.retired_through==closed.retired_through);

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
assert(~ri.accepted&&isequal(same,f)&&isempty(ri.steps)); % Не выдавать части отвергнутой истории.
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
% Индексы относятся к исходному списку, а модели заново получены по хронологии.
events=ri.steps(cellfun(@(z)isfield(z,'packet_index'),ri.steps));
assert(isequal(cellfun(@(z)z.packet_index,events),[1 2 3 4]));
assert(isequal(cellfun(@(z)z.interval_index,events),[1 1 1 2]));
assert(strcmp(events{3}.kind,'save_reference'));
assert(isequal(size(events{3}.data.reference_change.transition),[42 42]));
[again_trace,again_info]=nav_replay(f,{u,second},{late,cam,save,cam2},{sensor,[],[],[]},e,100);
assert(isequal(again_trace,replayed)&&isequaln(again_info.steps,ri.steps));

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
% Общий сдвиг часов не меняет относительный перенос при постоянной скорости.
time_a=a; time_b=b; time_a.v=[3;4;0]; time_b.v=time_a.v;
timed=m; timed.time_covariance=diag([.01 .04]); timed.time_independent=true;
timed.time_linearization=struct('body_rates',zeros(3,2),'epochs',[a.t b.t], ...
    'unit','rad/s','reference','earth','valid',true);
time_obs=nav_odometry_measurement(time_a,time_b,timed);
close(time_obs.time_jacobian,[-time_a.v time_b.v],1e-12);
close(time_obs.R,m.R+.05*(time_a.v*time_a.v.'),1e-12);
timed.time_covariance=.04*ones(2);
close(nav_odometry_measurement(time_a,time_b,timed).R,m.R,1e-12);
% В одном наблюдении времени связывают ошибки направления и поворота.
timed.mode='direction_rotation'; timed.z=[.6;.8;0;1;0;0;0];
timed.tangent_basis=[-.8 0;.6 0;0 1]; timed.min_displacement=.01; timed.max_angle=.3;
timed.units={'1';'1';'rad';'rad';'rad'}; timed.R=eye(5); timed.scale_valid=false;
time_b.v=[1;-2;.5]; timed.time_covariance=[.01 .005;.005 .04];
timed.time_linearization.body_rates=[.1 -.2;.2 .1;.3 .2];
time_obs=nav_odometry_measurement(time_a,time_b,timed);
plain=rmfield(timed,{'time_covariance','time_independent','time_linearization'});
base=nav_odometry_measurement(time_a,time_b,plain);
close(time_obs.R,base.R+time_obs.time_jacobian*timed.time_covariance*time_obs.time_jacobian.',1e-12);
assert(norm(time_obs.R(1:2,3:5),'fro')>0);
for epoch=1:2
    state=time_a; if epoch==2, state=time_b; end
    delta=zeros(21,1); delta(1:3)=state.v*1e-6;
    delta(7:9)=timed.time_linearization.body_rates(:,epoch)*1e-6;
    if epoch==1
        plus=nav_odometry_measurement(nav_apply_error(time_a,delta),time_b,plain);
        minus=nav_odometry_measurement(nav_apply_error(time_a,-delta),time_b,plain);
    else
        plus=nav_odometry_measurement(time_a,nav_apply_error(time_b,delta),plain);
        minus=nav_odometry_measurement(time_a,nav_apply_error(time_b,-delta),plain);
    end
    close(-(plus.residual-minus.residual)/(2e-6),time_obs.time_jacobian(:,epoch),1e-8);
end
[~,Qbase,accepted]=nav_correct({time_a,time_b},eye(42),base,100); assert(accepted.accepted);
[~,Qtime,accepted]=nav_correct({time_a,time_b},eye(42),time_obs,100); assert(accepted.accepted);
assert(trace(Qtime)>trace(Qbase)&&min(eig(Qtime-Qbase))>=-1e-10);
bad_time=timed; bad_time.time_independent=false;
fails(@()nav_odometry_measurement(time_a,time_b,bad_time),'libekf:Correlation');
bad_time=timed; bad_time.time_linearization.unit='deg/s';
fails(@()nav_odometry_measurement(time_a,time_b,bad_time),'libekf:Time');
bad_time=timed; bad_time.time_covariance=[1 2;2 1];
fails(@()nav_odometry_measurement(time_a,time_b,bad_time),'libekf:Covariance');

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
% Соответствия извлекаются из изображений; перенос в пикселях известен независимо.
stream=RandStream('mt19937ar','Seed',17);
texture=conv2(rand(stream,240,320),ones(3)/9,'same');
texture=uint8(255*(texture-min(texture(:)))/(max(texture(:))-min(texture(:))));
shifted=zeros(size(texture),'uint8'); shifted(5:end,7:end)=texture(1:end-4,1:end-6);
frame=struct('image',texture,'source','visible_camera','id','image_a', ...
    't',10,'t_received',10.1,'quality','VALID');
next=frame; next.image=shifted; next.id='image_b'; next.t=10.2; next.t_received=10.3;
calibration=struct('source',frame.source,'valid',true, ...
    'intrinsics',cameraIntrinsics([250 250],[160 120],[240 320]));
options=struct('max_points',1000,'match_ratio',.7,'fit_fraction',.5,'seed',9);
before=rng; matches=nav_visual_matches(frame,next,calibration,options); after=rng;
assert(isequal(before,after)&&strcmp(matches.quality,'UNASSESSED')&&~matches.scale_known);
assert(size(matches.first,1)>=12&&nnz(matches.training)>=8&&nnz(~matches.training)>=4);
delta=matches.second-matches.first;
assert(median(vecnorm(delta-[6 4],2,2))<1);
assert(isequal(matches.lineage,{frame.id,next.id})&&matches.t_start==frame.t&& ...
    matches.t_end==next.t&&matches.t_received==next.t_received&&strcmp(matches.units,'pixel'));
assert(isequal(matches,nav_visual_matches(frame,next,calibration,options)));
close(matches.lineage_times,[frame.t next.t],1e-12);
% Начало экспозиции приводится к её середине; строки сохраняют собственное время.
capture_first=frame; capture_second=next;
capture_first.capture=struct('duration',.02,'readout_time',0,'stamp_reference','first_row_start', ...
    'row_direction',1,'unit','s','valid',true); capture_second.capture=capture_first.capture;
capture_options=options; capture_options.max_points=100;
global_capture=nav_visual_matches(capture_first,capture_second,calibration,capture_options);
close([global_capture.t_start global_capture.t_end],[frame.t next.t]+.01,1e-12);
close(global_capture.point_times,repmat([global_capture.t_start global_capture.t_end],size(global_capture.first,1),1),1e-12);
end_first=capture_first; end_second=capture_second;
end_first.t=end_first.t+.02; end_second.t=end_second.t+.02;
end_first.capture.stamp_reference='last_row_end'; end_second.capture.stamp_reference='last_row_end';
same_capture=nav_visual_matches(end_first,end_second,calibration,capture_options);
close(same_capture.point_times,global_capture.point_times,1e-12);
close(same_capture.lineage_times,global_capture.lineage_times,1e-12);
close(global_capture.lineage_times,[frame.t next.t],1e-12);
capture_first.capture.readout_time=.04; capture_second.capture.readout_time=.04;
rolling_matches=nav_visual_matches(capture_first,capture_second,calibration,capture_options);
close(rolling_matches.point_times(:,1),frame.t+.01+(rolling_matches.first(:,2)-1)/239*.04,1e-10);
reverse_first=capture_first; reverse_second=capture_second;
reverse_first.capture.row_direction=-1; reverse_second.capture.row_direction=-1;
reverse_capture=nav_visual_matches(reverse_first,reverse_second,calibration,capture_options);
close(reverse_capture.point_times(:,1),frame.t+.01+(240-reverse_capture.first(:,2))/239*.04,1e-10);
close(reverse_capture.lineage_times,rolling_matches.lineage_times,1e-12);
early=capture_first; early.t_received=early.t+.01;
fails(@()nav_visual_matches(early,capture_second,calibration,capture_options),'libekf:Time');

% Общая ошибка главной точки камеры сохраняется между точками и кадрами.
noise_options=options; noise_options.max_points=100;
noise_options.localization_covariance=.04*eye(2); noise_options.localization_independent=true;
noise_calibration=calibration; noise_calibration.parameter_covariance=zeros(10);
noise_calibration.parameter_noise_independent=true;
independent_pixels=nav_visual_matches(frame,next,noise_calibration,noise_options);
dimension=4*size(independent_pixels.first,1); assert(dimension>=48);
close(full(independent_pixels.pixel_covariance),.04*eye(dimension),1e-7);
% Ошибка вертикальной координаты строки связана с ошибкой её времени.
joint_options=noise_options; joint_options.capture_covariance=zeros(6);
joint_options.capture_noise_independent=true;
joint=nav_visual_matches(capture_first,capture_second,noise_calibration,joint_options);
rate=.04/239; Cjoint=joint.coordinate_time_covariance;
close(full(Cjoint(2,3)),.04*rate,1e-10); close(full(Cjoint(3,3)),.04*rate^2,1e-12);
count=2*size(joint.first,1); pixel_rows=reshape([3*(1:count)-2;3*(1:count)-1],[],1);
close(full(Cjoint(pixel_rows,pixel_rows)),full(joint.pixel_covariance),1e-10);
joint_options.capture_covariance([1 4],[1 4])=.001;
shared_time=nav_visual_matches(capture_first,capture_second,noise_calibration,joint_options);
time_rows=3:3:size(Cjoint,1);
close(full(shared_time.coordinate_time_covariance(time_rows,time_rows)-Cjoint(time_rows,time_rows)), ...
    .001*ones(count),1e-12);
reverse_joint=nav_visual_matches(reverse_first,reverse_second,noise_calibration,joint_options);
close(full(reverse_joint.coordinate_time_covariance(2,3)),-.04*rate,1e-10);
joint_options.capture_noise_independent=false;
fails(@()nav_visual_matches(capture_first,capture_second,noise_calibration,joint_options),'libekf:Correlation');

noise_calibration.parameter_covariance(3,3)=1;
shared_pixels=nav_visual_matches(frame,next,noise_calibration,noise_options);
common=zeros(dimension,1); common(1:2:end)=1;
close(full(shared_pixels.pixel_covariance-independent_pixels.pixel_covariance),common*common.',1e-6);
close(shared_pixels.calibration_jacobian(:,3),-common,1e-7);
xrows=1:2:dimension;
assert(abs(mean(shared_pixels.pixel_covariance(xrows,xrows),'all')-1-.04/numel(xrows))<1e-6);
bad_calibration=noise_calibration; bad_calibration.parameter_noise_independent=false;
fails(@()nav_visual_matches(frame,next,bad_calibration,noise_options),'libekf:Correlation');
bad_calibration=noise_calibration; bad_calibration.parameter_covariance(1,3)=1;
bad_calibration.parameter_covariance(3,1)=1;
fails(@()nav_visual_matches(frame,next,bad_calibration,noise_options),'libekf:Covariance');

blank=frame; blank.image(:)=0; blank_next=next; blank_next.image(:)=0;
missing=nav_visual_matches(blank,blank_next,calibration,options);
assert(strcmp(missing.quality,'INSUFFICIENT')&&isempty(missing.first));
bad=next; bad.quality='FAILED';
fails(@()nav_visual_matches(frame,bad,calibration,options),'libekf:Quality');
bad=next; bad.t_received=bad.t-1;
fails(@()nav_visual_matches(frame,bad,calibration,options),'libekf:Time');
bad=next; bad.id=frame.id;
fails(@()nav_visual_matches(frame,bad,calibration,options),'libekf:Vision');
bad_calibration=calibration; bad_calibration.source='another_camera';
fails(@()nav_visual_matches(frame,next,bad_calibration,options),'libekf:Quality');
% Направление и поворот проверяются по независимо построенным проекциям 3D-сцены.
points=[4*rand(stream,2,200)-2;4+8*rand(stream,1,200)];
Rtrue=nav_rotation([cos(.04);0;sin(.04);0]); shift=[.5;.15;.03];
K=[250 0 160;0 250 120;0 0 1];
one=K*points; two=K*(Rtrue.'*(points-shift));
pair=matches; pair.first=(one(1:2,:)./one(3,:)).'; pair.second=(two(1:2,:)./two(3,:)).';
pair.training=false(200,1); pair.training(1:2:end)=true;
limits=struct('trials',1000,'confidence',99.9,'fit_distance',.1,'check_distance',.1, ...
    'min_check_support',.9,'min_front_fraction',.9,'seed',123, ...
    'ambiguity',struct('rotation_rms',.1,'plane_rms',.1,'point_limit',.1,'min_support',.9));
before=rng; motion=nav_visual_motion(pair,calibration.intrinsics,limits); after=rng;
unsupported=nav_visual_motion(rolling_matches,calibration.intrinsics,limits);
assert(strcmp(unsupported.quality,'UNSUPPORTED_TIME_MODEL')&&~unsupported.geometry_ready&&isempty(unsupported.rotation));

assert(isequal(before,after)&&motion.geometry_ready&&strcmp(motion.quality,'UNASSESSED'));
dated=pair; dated.lineage_times=[pair.t_start pair.t_end];
dated_motion=nav_visual_motion(dated,calibration.intrinsics,limits);
assert(isequal(dated_motion.lineage_times,dated.lineage_times));
bad=dated; bad.lineage_times(1)=pair.t_start+1;
fails(@()nav_visual_motion(bad,calibration.intrinsics,limits),'libekf:LineageTime');

close(motion.rotation,Rtrue,1e-5); close(motion.direction,shift/norm(shift),1e-5);
assert(~motion.scale_known&&~isfield(motion,'R'));
% Несколько допустимых наборов должны сохранять одну геометрию либо явную неоднозначность.
fitting=find(pair.training); selection_limits=limits;
selection_limits.selection=struct('subsets',{{fitting(1:2:end),fitting(2:2:end)}}, ...
    'rotation_limit',.01,'direction_limit',.05);
stable=nav_visual_motion(pair,calibration.intrinsics,selection_limits);
assert(stable.geometry_ready&&strcmp(stable.selection.status,'CONSISTENT'));
assert(numel(stable.selection.candidates)==3&&strcmp(stable.quality,'UNASSESSED'));
% Две пространственные группы имеют разные движения и обе представлены в проверочной части.
mixed=pair; other_R=nav_rotation([cos(.03);sin(.03);0;0]); other_shift=[-.4;.2;.1];
projection=K*(other_R.'*(points(:,101:end)-other_shift));
mixed.second(101:end,:)=(projection(1:2,:)./projection(3,:)).';
selection_limits.min_check_support=.4;
selection_limits.selection.subsets={fitting(fitting<=100),fitting(fitting>100)};
ambiguous=nav_visual_motion(mixed,calibration.intrinsics,selection_limits);
assert(~ambiguous.geometry_ready&&strcmp(ambiguous.quality,'AMBIGUOUS')&& ...
    strcmp(ambiguous.selection.status,'AMBIGUOUS')&&~isfield(ambiguous,'local_error'));
assert(ambiguous.selection.direction_spread>selection_limits.selection.direction_limit);
assert(nnz(cellfun(@(x)x.geometry_ready,ambiguous.selection.candidates))>=2);
assert(isempty(ambiguous.rotation)&&isempty(ambiguous.direction));
selection_limits.rank_tolerance=1e-10;
with_covariance=nav_visual_motion(mixed,calibration.intrinsics,selection_limits,speye(800)*.04);
assert(strcmp(with_covariance.quality,'AMBIGUOUS')&&~isfield(with_covariance,'local_error'));
assert(nnz(cellfun(@(x)isfield(x,'covariance')&&~isempty(x.covariance),with_covariance.selection.candidates))>=2);

bad_selection=selection_limits; bad_selection.selection.subsets{1}(1)=find(~pair.training,1);
fails(@()nav_visual_motion(pair,calibration.intrinsics,bad_selection),'libekf:Vision');

% Полная погрешность сохраняет связи направления и поворота и масштаб входного шума.
limits.rank_tolerance=1e-10; C=speye(4*size(pair.first,1))*.04;
uncertain=nav_visual_motion(pair,calibration.intrinsics,limits,C);
local=uncertain.local_error;
assert(local.rank==5&&strcmp(local.quality,'LOCAL_LINEAR')&&isequal(size(local.covariance),[5 5]));
assert(norm(local.covariance(1:2,3:5),'fro')>0&&all(eig(local.covariance)>0));
scaled=nav_visual_motion(pair,calibration.intrinsics,limits,4*C);
close(scaled.local_error.covariance,4*local.covariance,1e-10);
common=zeros(size(C,1),1); common(2:4:end)=1;
dependent=nav_visual_motion(pair,calibration.intrinsics,limits,C+.01*(common*common.'));
assert(min(eig(dependent.local_error.covariance-local.covariance))>=-1e-10);
% Якобианы проверяются по независимому пересчёту эпиполярного ограничения.
chosen=uncertain.fitting_indices; a=K\[pair.first(chosen,:).';ones(1,numel(chosen))];
b=K\[pair.second(chosen,:).';ones(1,numel(chosen))];
R=uncertain.rotation; d=uncertain.direction; B=local.tangent_basis; h=1e-6;
for j=1:5
    dp=d; dm=d; Rp=R; Rm=R;
    if j<=2
        dp=d+h*B(:,j); dp=dp/norm(dp); dm=d-h*B(:,j); dm=dm/norm(dm);
    else
        delta=zeros(3,1); delta(j-2)=h;
        Rp=R*expm(nav_skew(delta)); Rm=R*expm(-nav_skew(delta));
    end
    fp=sum(b.*(Rp.'*nav_skew(dp)*a),1).'; fm=sum(b.*(Rm.'*nav_skew(dm)*a),1).';
    close((fp-fm)/(2*h),local.state_jacobian(:,j),1e-8);
end
for j=1:4
    x=pair.first(chosen(1),:); y=pair.second(chosen(1),:); xp=x; xm=x; yp=y; ym=y;
    if j<=2, xp(j)=xp(j)+h; xm(j)=xm(j)-h;
    else, yp(j-2)=yp(j-2)+h; ym(j-2)=ym(j-2)-h; end
    E=R.'*nav_skew(d);
    fp=(K\[yp.';1]).'*E*(K\[xp.';1]); fm=(K\[ym.';1]).'*E*(K\[xm.';1]);
    close((fp-fm)/(2*h),full(local.pixel_jacobian(1,j)),1e-8);
end

% Точки с большей погрешностью имеют меньший вес при уточнении движения.
noisy_pair=pair; noisy_rows=1:4:size(pair.first,1);
noisy_pair.second(noisy_rows,:)=noisy_pair.second(noisy_rows,:)+2*randn(stream,numel(noisy_rows),2);
pixel_variance=.01*ones(4*size(pair.first,1),1);
pixel_variance(4*noisy_rows-1)=4; pixel_variance(4*noisy_rows)=4;
weighted_C=spdiags(pixel_variance,0,numel(pixel_variance),numel(pixel_variance));
weighted_limits=limits; weighted_limits.fit_distance=8; weighted_limits.check_distance=2;
initial_motion=nav_visual_motion(noisy_pair,calibration.intrinsics,weighted_limits);
weighted_limits.refinement=struct('max_iterations',50,'angle_tolerance',1e-9,'pixel_tolerance',1e-7);
refined=nav_visual_motion(noisy_pair,calibration.intrinsics,weighted_limits,weighted_C);
assert(refined.refinement.converged&&refined.geometry_ready&&strcmp(refined.quality,'UNASSESSED'));
assert(norm(refined.direction-shift/norm(shift))<norm(initial_motion.direction-shift/norm(shift)));
rescaled=nav_visual_motion(noisy_pair,calibration.intrinsics,weighted_limits,16*weighted_C);
assert(rescaled.refinement.converged);
close(rescaled.rotation,refined.rotation,1e-8); close(rescaled.direction,refined.direction,1e-8);
close(rescaled.local_error.covariance,16*refined.local_error.covariance,1e-8);
spoiled=noisy_pair; spoiled.second(~pair.training,:)=spoiled.second(~pair.training,:)+[50 80];
heldout=nav_visual_motion(spoiled,calibration.intrinsics,weighted_limits,weighted_C);
assert(heldout.refinement.converged&&~heldout.geometry_ready);
close(heldout.rotation,refined.rotation,0); close(heldout.direction,refined.direction,0);
unfinished_limits=weighted_limits; unfinished_limits.refinement.max_iterations=1;
unfinished=nav_visual_motion(noisy_pair,calibration.intrinsics,unfinished_limits,weighted_C);
assert(~unfinished.refinement.converged&&~unfinished.geometry_ready&&~isfield(unfinished,'local_error'));
bad_pair=pair; bad_pair.t_received=pair.t_start;
fails(@()nav_visual_motion(bad_pair,calibration.intrinsics,limits),'libekf:Time');
bad_pair=pair; bad_pair.lineage={'repeated','repeated'};
fails(@()nav_visual_motion(bad_pair,calibration.intrinsics,limits),'libekf:Lineage');
% Изменение только проверочной части не меняет рассчитанное движение.
corrupted=pair; corrupted.second(~pair.training,:)=corrupted.second(~pair.training,:)+[50 80];
failed=nav_visual_motion(corrupted,calibration.intrinsics,limits);
close(failed.rotation,motion.rotation,0); close(failed.direction,motion.direction,0);
assert(~failed.geometry_ready&&~strcmp(failed.quality,'VALID'));
% При чистом вращении перенос не считается определённым.
two=K*(Rtrue.'*points); rotating=pair; rotating.second=(two(1:2,:)./two(3,:)).';
unresolved=nav_visual_motion(rotating,calibration.intrinsics,limits);
assert(~unresolved.geometry_ready&&unresolved.ambiguity.rotation_explains);

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
% Изоляция сохраняет все объяснения, допустимые по явно заданному числу отказов.
first=struct('value',10,'covariance',1,'quantity','height','units',{{'m'}}, ...
    'frame_id','local','t',10,'quality','VALID','source','solution_a', ...
    'input_sources',{{'a'}},'provenance_complete',true);
second=first; second.value=0; second.source='solution_b'; second.input_sources={'b'};
third=second; third.source='solution_c'; third.input_sources={'c'};
comparison=struct('max_age',1,'limit',3,'context','height_isolation');
[ab_sample,ab]=nav_compare_sources(first,second,0,10,comparison);
[ac_sample,ac]=nav_compare_sources(first,third,0,10,comparison);
[~,bc]=nav_compare_sources(second,third,0,10,comparison);
policy=struct('max_faults',1,'max_hypotheses',10,'max_age',1,'max_span',0);
isolated=nav_isolate_sources({'a','b','c'},{ab,ac,bc},10,policy);
assert(isolated.complete&&strcmp(isolated.status,'ISOLATED')&&isequal(isolated.suspected_sources,{'a'}));
ambiguous=nav_isolate_sources({'a','b','c'},{ab},10,policy);
assert(ambiguous.complete&&strcmp(ambiguous.status,'AMBIGUOUS')&&isempty(ambiguous.suspected_sources));
two_faults=policy; two_faults.max_faults=2;
ambiguous=nav_isolate_sources({'a','b','c'},{ab,ac,bc},10,two_faults);
assert(strcmp(ambiguous.status,'AMBIGUOUS')&&isempty(ambiguous.suspected_sources));
reordered=nav_isolate_sources({'c','a','b'},{bc,ac,ab},10,policy);
assert(isequal(reordered.suspected_sources,isolated.suspected_sources));
% Тревога на F/M/S содержит те же возможные физические причины.
watch_cfg=struct('windows',[1 2 3],'max_gap',1,'max_age',1,'capacity',10, ...
    'instant_limit',5,'mean_limit',[5 5 5],'rms_limit',[5 5 5]);
[~,wab]=nav_watch([],ab_sample,10,watch_cfg); [~,wac]=nav_watch([],ac_sample,10,watch_cfg);
wab.involved_sources=ab.involved_sources; wab.provenance_complete=ab.provenance_complete;
wac.involved_sources=ac.involved_sources; wac.provenance_complete=ac.provenance_complete;
from_watch=nav_isolate_sources({'a','b','c'},{wab,wac},10,policy);
assert(isequal(from_watch.suspected_sources,isolated.suspected_sources));
periodic=wab; periodic.interval=[8 10];
insufficient_period=nav_isolate_sources({'a','b','c'},{periodic,wac},10,policy);
assert(strcmp(insufficient_period.status,'NOT_READY'));
whole_period=policy; whole_period.max_span=2;
covered=nav_isolate_sources({'a','b','c'},{periodic,wac},10,whole_period);
assert(isequal(covered.suspected_sources,{'a'})&&isequal(covered.interval,[8 10]));

% Две неисправности объясняют те же данные иначе; согласующаяся пара не объявляется исправной.
conflict=bc; conflict.status='ALARM';
inconsistent=nav_isolate_sources({'a','b','c'},{ab,ac,conflict},10,policy);
assert(inconsistent.complete&&strcmp(inconsistent.status,'INCONSISTENT')&&isempty(inconsistent.suspected_sources));
small=two_faults; small.max_hypotheses=6;
incomplete=nav_isolate_sources({'a','b','c'},{ab,ac},10,small);
assert(~incomplete.complete&&strcmp(incomplete.status,'NOT_READY')&&isempty(incomplete.suspected_sources));
unknown=ab; unknown.provenance_complete=false;
incomplete=nav_isolate_sources({'a','b','c'},{unknown,ac},10,policy);
assert(~incomplete.complete&&isempty(incomplete.suspected_sources));
stale=nav_isolate_sources({'a','b','c'},{ab,ac},12,policy); assert(strcmp(stale.status,'NOT_READY'));
early=ac; early.t=9.5; early.interval=[9.5 9.5];
mixed=nav_isolate_sources({'a','b','c'},{ab,early},10,policy); assert(strcmp(mixed.status,'NOT_READY'));
% Общая ИНС входит в множество возможных причин каждого вычисленного решения.
first.input_sources={'imu','a'}; second.input_sources={'imu','b'}; third.input_sources={'imu','c'};
third.value=-10;
[~,ab]=nav_compare_sources(first,second,0,10,comparison);
[~,ac]=nav_compare_sources(first,third,0,10,comparison);
[~,bc]=nav_compare_sources(second,third,0,10,comparison);
physical=nav_isolate_sources({'imu','a','b','c'},{ab,ac,bc},10,policy);
assert(physical.complete&&isequal(physical.suspected_sources,{'imu'}));
single=struct('status','ALARM','t',10,'involved_sources',{{'a'}},'provenance_complete',true);
none=policy; none.max_faults=0;
impossible=nav_isolate_sources({'a'},{single},10,none);
assert(strcmp(impossible.status,'INCONSISTENT'));
single_fault=nav_isolate_sources({'a'},{single},10,policy);
assert(isequal(single_fault.suspected_sources,{'a'}));

a=struct('value',[2;0;0],'covariance',eye(3),'quantity','position','units',{{'m';'m';'m'}}, ...
    'frame_id','local','t',10,'quality','VALID','source','one');
b=a; b.value=zeros(3,1); b.source='two';
cfg=struct('max_age',1,'limit',3,'context','position_one_two_local');
% Неизвестная связь не превращается в нулевую ковариацию или подтверждение согласия.
[missing,unavailable]=nav_compare_sources(a,b,[],10,cfg);
assert(isempty(missing)&&strcmp(unavailable.status,'NOT_READY'));
bounded_cfg=cfg; bounded_cfg.unknown_correlation_bound=true;
[bounded,bound]=nav_compare_sources(a,b,[],10,bounded_cfg);
assert(bounded.alarm_only&&strcmp(bound.nis_kind,'LOWER_BOUND')&&strcmp(bound.status,'NOT_READY'));
close(bound.nis,1,1e-12);
for rho=[-.99 -.5 0 .5 .99]
    [~,exact]=nav_compare_sources(a,b,rho*eye(3),10,cfg);
    assert(bound.nis<=exact.nis+1e-12);
end
units=diag([1000 .01 10]); scaled_a=a; scaled_b=b;
scaled_a.value=units*a.value; scaled_b.value=units*b.value;
scaled_a.covariance=units*a.covariance*units.'; scaled_b.covariance=units*b.covariance*units.';
scaled_a.units={'mm';'hm';'dm'}; scaled_b.units=scaled_a.units;
[~,scaled_bound]=nav_compare_sources(scaled_a,scaled_b,[],10,bounded_cfg);
close(scaled_bound.nis,bound.nis,1e-12);
alarm_only_watch=struct('windows',[1 2 3],'max_gap',1,'max_age',1,'capacity',10, ...
    'instant_limit',10,'mean_limit',[.1 .1 .1],'rms_limit',[10 10 10]);
bounded_history=[];
for t=7:10
    entry=bounded; entry.t=t;
    [bounded_history,not_confirmed]=nav_watch(bounded_history,entry,t,alarm_only_watch);
end
assert(all(not_confirmed.ready)&&all(not_confirmed.mean>.1)&& ...
    strcmp(not_confirmed.status,'NOT_READY')&&not_confirmed.alarm_only);
large=a; other=b; large.value=[30;0;0]; large.t=11; other.t=11;
[large_sample,disagreement]=nav_compare_sources(large,other,[],11,bounded_cfg);
assert(strcmp(disagreement.status,'DISAGREEMENT'));
[~,detected]=nav_watch(bounded_history,large_sample,11,alarm_only_watch);
assert(strcmp(detected.status,'ALARM'));
known=a; known.t=11;
[exact_sample,~]=nav_compare_sources(known,other,zeros(3),11,cfg);
assert(~strcmp(exact_sample.context,bounded.context));
[changed_history,not_confirmed]=nav_watch(bounded_history,exact_sample,11,alarm_only_watch);
assert(numel(changed_history.t)==1&&~not_confirmed.alarm_only&&strcmp(not_confirmed.status,'NOT_READY'));

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
% LEKF-012/020: смена реального состава не сохраняет готовность старой истории.
a.value=zeros(3,1); b.value=zeros(3,1);
a.input_sources={'imu','gnss'}; a.provenance_complete=true;
b.input_sources={'gnss','dme'}; b.provenance_complete=true;
[first,shared]=nav_compare_sources(a,b,.5*eye(3),10,cfg);
assert(strcmp(shared.provenance,'SHARED_INPUTS')&&isequal(shared.shared_sources,{'gnss'}));
history=[];
for t=7:10
    check=first; check.t=t; [history,ready]=nav_watch(history,check,t,watch);
end
assert(strcmp(ready.status,'CONSISTENT'));
a.t=11; b.t=11; b.input_sources={'radio'};
[changed,disjoint]=nav_compare_sources(a,b,.5*eye(3),11,cfg);
assert(strcmp(disjoint.provenance,'DISJOINT_DECLARED_INPUTS')&&isempty(disjoint.shared_sources));
[reset,result]=nav_watch(history,changed,11,watch);
assert(strcmp(result.status,'NOT_READY')&&numel(reset.t)==1);
% Порядок аргументов и порядок перечисления приборов не меняют саму проверку.
a.value=[1;2;3]; a.input_sources={'gnss','imu','imu'};
cross=[.2 .1 0;0 .3 0;0 0 .1];
[forward,~]=nav_compare_sources(a,b,cross,11,cfg);
[reverse,~]=nav_compare_sources(b,a,cross.',11,cfg);
assert(strcmp(forward.context,reverse.context)); close(forward.value,reverse.value,1e-12);
unknown=rmfield(b,'provenance_complete');
[~,result]=nav_compare_sources(a,unknown,cross,11,cfg);
assert(strcmp(result.provenance,'UNKNOWN'));
renamed=b; renamed.source='third';
[different,~]=nav_compare_sources(a,renamed,cross,11,cfg);
assert(~strcmp(different.context,forward.context));
bad=b; bad.input_sources={''};
fails(@()nav_compare_sources(a,bad,cross,11,cfg),'libekf:Comparison');

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
% Полный путь: ИНС, три положения, смена опоры, определимость и географический выход.
[s,e]=fixture; s.p=zeros(3,1); s.v=[2;0;0]; s.q=[1;0;0;0];
s.ba=zeros(3,1); s.bg=zeros(3,1); e.Omega=zeros(3,1);
e.gravity_origin=zeros(3,1); e.gravity_gradient=zeros(3); e.gravity_radius=1e5;
e.max_step=1; e.state_noise_density=zeros(21);
imu=struct('model','inertial','source','imu','lineage',{{'history_imu_1'}},'quality','VALID', ...
    'frame','BODY','accel_unit','m/s^2','gyro_unit','rad/s','t_start',s.t, ...
    't_end',s.t+1,'t_received',s.t+1,'accel',zeros(3,1),'gyro',zeros(3,1),'noise_density',zeros(6));
second=imu; second.t_start=s.t+1; second.t_end=s.t+2; second.t_received=s.t+2;
second.lineage={'history_imu_2'};
filter=struct('state',s,'P',eye(21),'used',{{}},'capacity',20);
sensor=struct('id',s.frame_id,'origin_ecef',e.origin_ecef,'C_e_s',e.C_e_l);
packets=cell(1,3);
for j=1:3
    truth=s; truth.t=s.t+j-1; truth.p=s.p+s.v*(j-1);
    packets{j}=packet(truth,'gnss_position'); packets{j}.frame='SENSOR_LOCAL';
    packets{j}.lever=zeros(3,1); packets{j}.z=truth.p; packets{j}.R=eye(3);
    packets{j}.lineage={sprintf('history_position_%d',j)};
end
save=struct('model','save_reference','t_measure',s.t+1);
[state,trace]=nav_replay(filter,{imu,second},{packets{1},packets{2},save,packets{3}}, ...
    {sensor,sensor,[],sensor},e,100);
assert(trace.accepted&&all(cellfun(@(z)z.accepted,trace.observations)));
compiled=struct('steps',{trace.steps},'initial_times',s.t,'frame_id',s.frame_id, ...
    'noise_independent',true,'current_only',true);
target=zeros(6,21); target(:,1:6)=eye(6);
geometry=nav_observability_sequence(compiled,ones(21,1),target,1e-10);
assert(all(geometry.identified)&&geometry.t==s.t+2);
close(geometry.variance,[ones(3,1);6.5*ones(3,1)],1e-10);
control=struct('state_time',state.state.t,'frame_id',s.frame_id,'now',state.state.t,'checks',struct());
control.geometry=struct('history',compiled,'scales',ones(21,1), ...
    'relative_tolerance',1e-10,'max_variance',[1;1;10]);
published=nav_geodetic_output(state.state,state.P,e,control);
assert(all(published.observable)); close(published.measurement_variance(3),1,1e-9);
% Одна общая ошибка положения не уменьшается от повторов, но сокращается в разностях для скорости.
with_common=filter; with_common.future_noise=struct('ids',{{'cx','cy','cz'}},'units',{{'m','m','m'}}, ...
    'mean',zeros(3,1),'covariance',eye(3),'cross_cov',zeros(21,3),'state_times',s.t,'frame_id',s.frame_id);
for j=1:3
    packets{j}.future_noise_independent=true;
    packets{j}.shared_noise=struct('ids',{with_common.future_noise.ids}, ...
        'units',{with_common.future_noise.units},'map',eye(3),'residual_units',{{'m','m','m'}}, ...
        'residual_frame_id',sensor.id,'lineage',{{}});
end
imu.future_noise_independent=true; second.future_noise_independent=true;
[dependent,trace]=nav_replay(with_common,{imu,second},{packets{1},packets{2},save,packets{3}}, ...
    {sensor,sensor,[],sensor},e,100);
assert(trace.accepted&&all(cellfun(@(z)z.accepted,trace.observations)));
common_history=compiled; common_history.steps=trace.steps; common_history.initial_noise=with_common.future_noise;
geometry=nav_observability_sequence(common_history,ones(21,1),target,1e-10);
assert(all(geometry.identified)); close(geometry.variance,[2*ones(3,1);6.5*ones(3,1)],1e-9);
control.geometry.history=common_history;
common_output=nav_geodetic_output(dependent.state,dependent.P,e,control);
assert(all(common_output.observable)); close(common_output.measurement_variance(3),2,1e-9);
% История меняет местные оси после граничного наблюдения и сохранения опоры.
to=e; to.frame_id='history_new_local';
to.C_e_l=e.C_e_l*nav_rotation([.9;.1;-.2;.3]/norm([.9;.1;-.2;.3]));
to.origin_ecef=e.origin_ecef+[300;-200;50];
[changed,changed_trace]=nav_replay(with_common,{imu,second},{packets{1},packets{2},save,packets{3}}, ...
    {sensor,sensor,[],sensor},{e,to},100);
assert(changed_trace.accepted&&all(cellfun(@(z)z.accepted,changed_trace.observations)));
expected=nav_change_frame(dependent,[],e,to);
close(changed.state.p,expected.state.p,1e-9); close(changed.jointP,expected.jointP,1e-10);
close(changed.future_noise.mean,expected.future_noise.mean,1e-9);
close(changed.future_noise.cross_cov,expected.future_noise.cross_cov,1e-10);
assert(sum(cellfun(@(z)strcmp(z.kind,'change_frame'),changed_trace.steps))==1);
frame_history=common_history; frame_history.steps=changed_trace.steps;
geometry=nav_observability_sequence(frame_history,ones(21,1),target,1e-10);
assert(strcmp(geometry.frame_id,to.frame_id));
close(geometry.variance,[2*ones(3,1);6.5*ones(3,1)],1e-9);
changed_control=control; changed_control.frame_id=to.frame_id; changed_control.geometry.history=frame_history;
changed_output=nav_geodetic_output(changed.state,changed.P,to,changed_control);
close(changed_output.measurement_variance,common_output.measurement_variance,1e-9);
contradictory=to; contradictory.frame_id=e.frame_id;
fails(@()nav_replay(with_common,{imu,second},{},{},{e,contradictory},100),'libekf:Frame');
fails(@()nav_replay(with_common,{imu,second},{},{},{e},100),'libekf:Frame');
broken=second; broken.t_start=broken.t_start+.01;
[unchanged,rejected]=nav_replay(with_common,{imu,broken},{},{},{e,to},100);
assert(~rejected.accepted&&isequal(unchanged,with_common)&&isempty(rejected.steps));


unknown=compiled; unknown.noise_independent=false;
fails(@()nav_observability_sequence(unknown,ones(21,1),target,1e-10),'libekf:Correlation');
gap=compiled; first=find(cellfun(@(z)strcmp(z.kind,'propagation'),gap.steps),1);
gap.steps{first}.data.t_start=gap.steps{first}.data.t_start+.01;
fails(@()nav_observability_sequence(gap,ones(21,1),target,1e-10),'libekf:History');
rejected=compiled; first=find(cellfun(@(z)isfield(z.data,'observation'),rejected.steps),1);
rejected.steps{first}.data.accepted=false;
fails(@()nav_observability_sequence(rejected,ones(21,1),target,1e-10),'libekf:History');

% Сохранение копии не создаёт измерительной информации, даже в тот же момент.
clone=struct('H',{{1,zeros(0,2)}},'transitions',{{[1;1]}},'times',[0 0], ...
    'process_covariance',zeros(2),'measurement_covariance',9, ...
    'process_measurement_covariance',zeros(2,1),'frame_id','local');
result=nav_observability_sequence(clone,1,eye(2),1e-10);
assert(all(result.identified)&&result.rank==1); close(result.variance,[9;9],1e-12);
% Удаление наблюдаемой компоненты не делает ненаблюдаемую компоненту известной.
drop=clone; drop.H={[1 0],zeros(0,1)}; drop.transitions={[0 1]};
drop.process_covariance=0; drop.process_measurement_covariance=0;
result=nav_observability_sequence(drop,[1;1],1,1e-10);
assert(~result.identified&&isinf(result.variance));
% Прогноз до выхода без нового наблюдения накапливает ошибку движения.
coast=drop; coast.H={1,zeros(0,1)}; coast.transitions={1}; coast.times=[0 2];
coast.measurement_covariance=1; coast.process_covariance=4;
result=nav_observability_sequence(coast,1,1,1e-10);
assert(result.identified&&result.t==2); close(result.variance,5,1e-12);
coast.process_measurement_covariance=.5;
result=nav_observability_sequence(coast,1,1,1e-10); close(result.variance,4,1e-12);
% Пустой начальный узел допустим; измерение на выходе определяет текущую величину.
coast.H={zeros(0,1),1}; coast.process_measurement_covariance=0;
result=nav_observability_sequence(coast,1,1,1e-10); close(result.variance,1,1e-12);
bad=clone; bad.transitions={eye(2)};
fails(@()nav_observability_sequence(bad,1,eye(2),1e-10),'libekf:Observability');
bad=clone; bad.H={zeros(0,1),zeros(0,2)};
fails(@()nav_observability_sequence(bad,1,eye(2),1e-10),'libekf:Observability');

% Последовательность положений определяет скорость, которую один момент не определяет.
history=struct('H',{{[1 0],[1 0]}},'transitions',{{[1 1;0 1]}},'times',[0 1], ...
    'process_covariance',zeros(2),'measurement_covariance',eye(2), ...
    'process_measurement_covariance',zeros(2),'frame_id','local');
temporal=nav_observability_sequence(history,[1;1],[0 1],1e-10);
assert(temporal.identified&&temporal.rank==2); close(temporal.variance,2,1e-12);
instant=nav_observability([1 0],1,[1;1],[0 1],1e-10); assert(~instant.identified);
% Текущее положение содержит процессную ошибку, общую с поздним измерением.
history.H={1,1}; history.transitions={1}; history.process_covariance=4;
history.measurement_covariance=diag([1 9]); history.process_measurement_covariance=zeros(1,2);
temporal=nav_observability_sequence(history,1,1,1e-10);
close(temporal.variance,45/14,1e-12);
% Связь ошибки движения с первым измерением меняет ошибку прогноза: 1+4-2*.5=4.
history.process_measurement_covariance=[.5 0];
temporal=nav_observability_sequence(history,1,1,1e-10); close(temporal.variance,36/13,1e-12);
centimetres=history; centimetres.H={.01,.01};
centimetres.process_covariance=4e4; centimetres.process_measurement_covariance=[50 0];
same=nav_observability_sequence(centimetres,100,.01,1e-10);
close(same.variance,temporal.variance,1e-12);
bad=history; bad.process_measurement_covariance=[5 0];
fails(@()nav_observability_sequence(bad,1,1,1e-10),'libekf:Covariance');
bad=history; bad.times=[1 0];
fails(@()nav_observability_sequence(bad,1,1,1e-10),'libekf:Observability');
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
% История геометрии подключается к географическому выходу на своём последнем моменте.
temporal_control=control;
temporal_control.geometry=rmfield(control.geometry,{'H','R'});
temporal_control.geometry.history=struct('H',{{control.geometry.H(1,:),control.geometry.H(2,:)}}, ...
    'transitions',{{eye(21)}},'times',[s.t-1 s.t],'frame_id',s.frame_id, ...
    'process_covariance',zeros(21),'measurement_covariance',control.geometry.R, ...
    'process_measurement_covariance',zeros(21,2));
temporal_output=nav_geodetic_output(s,eye(21),e,temporal_control);
assert(isequal(temporal_output.observable,out.observable));
close(temporal_output.measurement_variance(3),out.measurement_variance(3),1e-10);
assert(strcmp(temporal_output.quality.altitude,'CONSISTENT'));
% Географический выход относится к текущему состоянию после удаления старой опоры.
reduced=temporal_control; reduced.geometry.scales=ones(42,1);
reduced.geometry.history.H={[zeros(2,21) control.geometry.H],zeros(0,21)};
reduced.geometry.history.transitions={[zeros(21) eye(21)]};
reduced_output=nav_geodetic_output(s,eye(21),e,reduced);
assert(isequal(reduced_output.observable,out.observable));
close(reduced_output.measurement_variance(3),out.measurement_variance(3),1e-10);

bad_control=temporal_control; bad_control.geometry.history.times(end)=s.t-1;
fails(@()nav_geodetic_output(s,eye(21),e,bad_control),'libekf:WatchOutput');
bad_control=temporal_control; bad_control.geometry.history.frame_id='another_frame';
fails(@()nav_geodetic_output(s,eye(21),e,bad_control),'libekf:WatchOutput');

[~,alarm]=nav_watch([],struct('t',s.t,'context','height','value',6),s.t,cfg);
control.checks.altitude={alarm}; control.geometry.H(2,:)=o.H;
out=nav_geodetic_output(s,eye(21),e,control);
assert(strcmp(out.quality.altitude,'ALARM')&&~out.observable(3));

end

function models=source_models
models={'gnss_position','external_position','gnss_velocity','external_velocity', ...
    'gnss_position_velocity','external_position_velocity','dme_range','pseudorange', ...
    'range_rate','vor_azimuth','baro_height','airdata','radar_range','ils_angles', ...
    'gls_deviation','attitude','radar_footprint','doppler_velocity','doppler_frequency', ...
    'adf_bearing','station_range_azimuth','range_difference','time_difference','rsdn_time_difference','rsdn_arrival_times','magnetic_vector','optical_landmarks','celestial_pixels','cooperative_range','static_pressure','impact_pressure','air_pressures'};
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
[~,~,transition,process_noise]=nav_predict(f.state,f.P,u,e);
close(info.propagation.transition,transition,1e-12);
close(info.propagation.process_covariance,process_noise,1e-12);
close(plain.P,transition*f.P*transition.'+process_noise,1e-12);
assert(info.propagation.t_start==u.t_start&&info.propagation.t_end==u.t_end);

% LEKF-004/020: закрытие истории сохраняет оценку и освобождает только идентификаторы.
seed=f; seed.used={'initial_record'}; seed.capacity=2; seed.excluded_sources={'failed_device'};
checkpoint=nav_checkpoint(seed);
assert(isempty(checkpoint.used)&&checkpoint.retired_through==s.t);
preserved=rmfield(checkpoint,'retired_through'); preserved.used=seed.used;
assert(isequal(preserved,seed));
timed_imu=u; timed_imu.lineage_times=u.t_end;
[closed_plain,accepted]=nav_assimilate(checkpoint,timed_imu,[],e,1e8);
assert(accepted.accepted&&numel(closed_plain.used)==1);
[unchanged,accepted]=nav_assimilate(checkpoint,u,[],e,1e8);
assert(~accepted.accepted&&strcmp(accepted.reason,'libekf:LineageTime')&&isequal(unchanged,checkpoint));
old_imu=timed_imu; old_imu.lineage_times=s.t;
[unchanged,accepted]=nav_assimilate(checkpoint,old_imu,[],e,1e8);
assert(~accepted.accepted&&strcmp(accepted.reason,'libekf:RetiredHistory')&&isequal(unchanged,checkpoint));
% Повторные границы дают продолжать счисление при памяти всего на один новый отсчёт.
bounded=checkpoint; bounded.capacity=1; unbounded=f;
for step=1:12
    next=timed_imu; next.t_start=bounded.state.t; next.t_end=next.t_start+.125;
    next.t_received=next.t_end; next.lineage={sprintf('bounded_imu_%d',step)};
    next.lineage_times=next.t_end;
    [bounded,accepted]=nav_assimilate(bounded,next,[],e,1e8); assert(accepted.accepted);
    [unbounded,accepted]=nav_assimilate(unbounded,next,[],e,1e8); assert(accepted.accepted);
    bounded=nav_checkpoint(bounded);
end
close(bounded.P,unbounded.P,1e-12); close(bounded.state.p,unbounded.state.p,1e-12);
assert(isempty(bounded.used)&&numel(unbounded.used)==12);

required=fieldnames(u);
for j=1:numel(required)
    incomplete=rmfield(u,required{j});
    [unchanged,failed]=nav_assimilate(f,incomplete,[],e,1e8);
    assert(~failed.accepted&&isequal(unchanged,f));
end
bad=u; bad.t_start=[s.t s.t];
[unchanged,failed]=nav_assimilate(f,bad,[],e,1e8);
assert(~failed.accepted&&isequal(unchanged,f));
% Запрет сторожа действует даже при VALID в пакете; другие источники остаются доступны.
report=struct('model','source_health','source','watch','lineage',{{'exclude_source'}}, ...
    'quality','VALID','target_source','test','valid',false,'frame_id',s.frame_id, ...
    't_measure',s.t,'t_received',u.t_end);
models=source_models;
for k=1:numel(models)
    m=packet(mid,models{k}); m.frame='SENSOR_LOCAL';
    % Причина двух тревог становится запретом только при явном ограничении числа отказов.
    alarm=struct('status','ALARM','t',s.t,'interval',[s.t s.t], ...
        'provenance_complete',true,'involved_sources',{{m.source,'aux_a'}});
    second_alarm=alarm; second_alarm.involved_sources={m.source,'aux_b'};
    policy=struct('max_faults',1,'max_hypotheses',8,'max_age',0,'max_span',0);
    isolated=nav_isolate_sources({m.source,'aux_a','aux_b'},{alarm,second_alarm},s.t,policy);
    assert(strcmp(isolated.status,'ISOLATED')&&isolated.complete&&numel(isolated.suspected_sources)==1);
    report.target_source=isolated.suspected_sources{1};
    [admitted,health]=nav_assimilate(f,report,[],e,1e8);
    assert(health.accepted&&ismember(m.source,admitted.excluded_sources));
    close(admitted.P,f.P,0); close(admitted.state.p,f.state.p,0);
    policy.max_faults=2;
    ambiguous=nav_isolate_sources({m.source,'aux_a','aux_b'},{alarm,second_alarm},s.t,policy);
    assert(strcmp(ambiguous.status,'AMBIGUOUS')&&isempty(ambiguous.suspected_sources));
    [blocked,info]=nav_replay(f,{u},{report,m},{[],sensor},e,1e8);
    assert(info.accepted&&info.observations{1}.accepted&&~info.observations{2}.accepted);
    close(blocked.P,plain.P,1e-12); close(blocked.state.p,plain.state.p,1e-12);
    restored=report; restored.valid=true; restored.lineage={'restore_source'};
    restored.t_measure=mid.t;
    [~,info]=nav_replay(f,{u},{report,restored,m},{[],[],sensor},e,1e8);
    assert(info.accepted&&all(cellfun(@(x)x.accepted,info.observations)));
    independent=m; independent.source='other_source';
    [~,info]=nav_replay(f,{u},{report,independent},{[],sensor},e,1e8);
    assert(info.accepted&&all(cellfun(@(x)x.accepted,info.observations)));
    % Имя вычислителя не скрывает отказ прибора, данные которого он использует.
    dependent=independent; dependent.input_sources={m.source,'healthy_auxiliary'};
    [rejected,info]=nav_replay(f,{u},{report,dependent},{[],sensor},e,1e8);
    assert(info.accepted&&~info.observations{2}.accepted&&isequal(rejected,blocked));
    [~,info]=nav_replay(f,{u},{report,restored,dependent},{[],[],sensor},e,1e8);
    assert(info.accepted&&all(cellfun(@(x)x.accepted,info.observations)));
    invalid=independent; invalid.input_sources={''};
    at_measurement=f; at_measurement.state=mid; at_measurement.P=midP;
    [unchanged,info]=nav_assimilate(at_measurement,invalid,sensor,e,1e8);
    assert(~info.accepted&&isequal(unchanged,at_measurement));
end
% Пересчёт от прежнего снимка исключает также производное измерение отказавшего прибора.
suspect=packet(mid,'external_position'); suspect.frame='SENSOR_LOCAL';
suspect.source='position_processor'; suspect.input_sources={'failed_position_sensor'};
suspect.z=suspect.z+[25;0;0];
[contaminated,trace]=nav_replay(f,{u},{suspect},{sensor},e,1e8);
assert(trace.accepted&&trace.observations{1}.accepted&& ...
    any(strcmp({contaminated.source_history.id},'failed_position_sensor')));
assert(norm(contaminated.state.p-plain.state.p)>1);
clean_start=f; clean_start.excluded_sources={'failed_position_sensor'};
[recovered,trace]=nav_replay(clean_start,{u},{suspect},{sensor},e,1e8);
assert(trace.accepted&&~trace.observations{1}.accepted&& ...
    ~any(strcmp({recovered.source_history.id},'failed_position_sensor')));
close(recovered.state.p,plain.state.p,1e-12); close(recovered.P,plain.P,1e-12);
% Запрет в конце истории не отменяет уже принятую поправку.
late_report=report; late_report.target_source='failed_position_sensor';
late_report.t_measure=u.t_end; late_report.t_received=u.t_end;
[still_affected,trace]=nav_replay(f,{u},{suspect,late_report},{sensor,[]},e,1e8);
assert(trace.accepted&&all(cellfun(@(x)x.accepted,trace.observations)));
close(still_affected.state.p,contaminated.state.p,1e-12);
close(still_affected.P,contaminated.P,1e-12);
assert(any(strcmp({still_affected.source_history.id},'failed_position_sensor')));
% Отказ ИНС внутри интервала не допускает продолжения по тем же данным.
report.target_source=u.source; report.t_measure=mid.t;
[unchanged,info]=nav_assimilate_interval(f,u,{report},{[]},e,1e8);
assert(~info.accepted&&strcmp(info.reason,'libekf:SourceHealth')&&isequal(unchanged,f)&&isempty(info.steps));
report.target_source='accelerometer';
dependent_imu=u; dependent_imu.input_sources={'accelerometer','gyroscope'};
[unchanged,info]=nav_assimilate_interval(f,dependent_imu,{report},{[]},e,1e8);
assert(~info.accepted&&strcmp(info.reason,'libekf:SourceHealth')&&isequal(unchanged,f)&&isempty(info.steps));
models=source_models;
for k=1:numel(models)
    m=packet(mid,models{k}); m.frame='SENSOR_LOCAL'; m.t_received=u.t_end+1;
    at_measurement=f; at_measurement.state=mid; at_measurement.P=midP;
    [reference,info]=nav_assimilate(at_measurement,m,sensor,e,1e8);
    assert(info.accepted,'libekf:Timing','Reference %s rejected: %s',models{k},info.reason);
    [reference.state,reference.P]=nav_predict(reference.state,reference.P,second,e);
    % История принятого интервала сохраняет хронологию и воспроизводит погрешность.
    [interval,trace]=nav_assimilate_interval(f,u,{m},{sensor},e,1e8);
    assert(trace.accepted&&numel(trace.steps)==3);
    a=trace.steps{1}.data; b=trace.steps{2}.data; c=trace.steps{3}.data;
    assert(strcmp(trace.steps{1}.kind,'propagation')&&trace.steps{2}.packet_index==1);
    assert(a.t_start==s.t&&a.t_end==mid.t&&c.t_start==mid.t&&c.t_end==u.t_end);
    prior=a.transition*f.P*a.transition.'+a.process_covariance;
    [~,posterior]=nav_correct(mid,prior,b.observation,1e8);
    close(interval.P,c.transition*posterior*c.transition.'+c.process_covariance,1e-10);
    close(interval.P,reference.P,1e-10);
    bad=m; bad.quality='FAILED';
    [~,trace]=nav_assimilate_interval(f,u,{bad},{sensor},e,1e8);
    assert(trace.accepted&&numel(trace.steps)==1);
    assert(trace.steps{1}.data.t_start==s.t&&trace.steps{1}.data.t_end==u.t_end);

    [replayed,info]=nav_replay(f,{u},{m},{sensor},e,1e8);
    assert(info.accepted&&info.observations{1}.accepted,'libekf:Timing', ...
        'Delayed %s rejected',models{k});
    close(replayed.state.p,reference.state.p,1e-10);
    close(replayed.state.v,reference.state.v,1e-10);
    close(nav_rotation(replayed.state.q),nav_rotation(reference.state.q),1e-10);
    close(replayed.P,reference.P,1e-10);
    [~,direct]=nav_assimilate_interval(f,u,{m},{sensor},e,1e8);
    assert(numel(info.steps)==numel(direct.steps));
    for j=1:numel(info.steps)
        assert(info.steps{j}.interval_index==1);
        assert(isequaln(rmfield(info.steps{j},'interval_index'),direct.steps{j}));
    end

    assert(replayed.state.t==u.t_end&&numel(replayed.used)==2);
    % Каждый источник обрабатывается после смены общей системы, в своих исходных осях.
    to=e; to.frame_id='timing_new_local';
    to.C_e_l=e.C_e_l*nav_rotation([.9;.1;-.2;.3]/norm([.9;.1;-.2;.3]));
    to.origin_ecef=e.origin_ecef+[300;-200;50];
    part1=u; part1.t_end=s.t+.03125; part1.lineage={'frame_interval_1'};
    part2=u; part2.t_start=part1.t_end; part2.lineage={'frame_interval_2'};
    shifted=m;
    if isfield(shifted,'peer_cross_frame_id')
        [~,~,J]=nav_change_frame(s,eye(21),e,to);
        shifted.peer_cross_frame_id=to.frame_id; shifted.peer_cross_cov=J*shifted.peer_cross_cov;
    end
    [fixed,accepted]=nav_replay(f,{part1,part2},{m},{sensor},e,1e8); assert(accepted.accepted);
    [changed,accepted]=nav_replay(f,{part1,part2},{shifted},{sensor},{e,to},1e8);
    assert(accepted.accepted&&accepted.observations{1}.accepted,'libekf:Frame', ...
        'Frame history rejected %s: %s',models{k},accepted.observations{1}.reason);
    expected=nav_change_frame(fixed,[],e,to);
    close(changed.state.p,expected.state.p,1e-8); close(changed.P,expected.P,1e-8);

    % Каждая физическая модель сохраняет результат после закрытия старой истории.
    dated=m; dated.lineage_times=m.t_measure;
    [bounded_result,bounded_info]=nav_replay(checkpoint,{timed_imu},{dated},{sensor},e,1e8);
    assert(bounded_info.accepted&&bounded_info.observations{1}.accepted&&numel(bounded_result.used)==2);
    close(bounded_result.P,replayed.P,1e-12); close(bounded_result.state.p,replayed.state.p,1e-12);
    % Новое имя вычисленного пакета не делает старые исходные данные новыми.
    dated.source='derived'; dated.lineage={'old_input'}; dated.lineage_times=s.t;
    [rejected,bounded_info]=nav_replay(checkpoint,{timed_imu},{dated},{sensor},e,1e8);
    assert(bounded_info.accepted&&strcmp(bounded_info.observations{1}.reason,'libekf:RetiredHistory'));
    assert(isequal(rejected,closed_plain));

    inertial_source=replayed.source_history(strcmp({replayed.source_history.id},u.source));
    measured_source=replayed.source_history(strcmp({replayed.source_history.id},m.source));
    assert(inertial_source.first_epoch==u.t_start&&inertial_source.last_epoch==u.t_end);
    assert(measured_source.first_epoch==m.t_measure&&measured_source.last_epoch==m.t_measure);

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

function requirement_timed_rays
% Подготовка движения из ИНС сохраняет зависимость от навигационного состояния.
[first,environment]=fixture; first.v=zeros(3,1); first.q=[1;0;0;0];
second=first; second.t=first.t+.2;
environment.Omega=[0;0;7.292115e-5]; environment.gravity_origin=[0;0;-9.81];
environment.gravity_gradient=zeros(3); environment.gravity_radius=1e6;
raw=struct('accel',[first.ba+[0;0;9.81],second.ba+[0;0;9.81]], ...
    'gyro',[first.bg+environment.Omega,second.bg+environment.Omega], ...
    'epochs',[first.t second.t],'t_received',second.t,'quality','VALID','calibration_valid',true, ...
    'frame','BODY','accel_unit','m/s^2','gyro_unit','rad/s','source','imu', ...
    'lineage',{{'imu_at_first','imu_at_second'}},'lineage_times',[first.t second.t], ...
    'covariance',eye(12)*1e-4,'state_noise_cov',zeros(42,12),'state_noise_frame_id',first.frame_id);
P=eye(42)*1e-3; validity=struct('max_offset',.02,'min_baseline',.01);
prepared=nav_inertial_camera_motion({first,second},P,raw,environment,validity);
close(prepared.motion.acceleration,zeros(3,2),1e-12);
close(prepared.motion.body_rates,[first.bg second.bg],1e-12);
assert(norm(prepared.state_noise_cov,'fro')>0&&strcmp(prepared.noise_model,'joint_required'));
% Численные производные номинально подготовленного входа имеют обратный знак:
% ошибка входа определяется как подготовленное значение минус истинное.
for j=1:42
    d=zeros(21,1); index=mod(j-1,21)+1; d(index)=1e-6;
    plus={first,second}; minus=plus; state_index=1+(j>21);
    plus{state_index}=nav_apply_error(plus{state_index},d);
    minus{state_index}=nav_apply_error(minus{state_index},-d);
    op=nav_inertial_camera_motion(plus,P,raw,environment,validity);
    om=nav_inertial_camera_motion(minus,P,raw,environment,validity);
    value_p=[op.motion.body_rates(:);op.motion.acceleration(:)];
    value_m=[om.motion.body_rates(:);om.motion.acceleration(:)];
    close((value_p-value_m)/(2e-6),-prepared.state_transform(:,j),1e-8);
end
% Полная зависимость исходных ошибок сохраняется, а невозможная отклоняется.
raw.state_noise_cov(1,1)=1e-5;
dependent=nav_inertial_camera_motion({first,second},P,raw,environment,validity);
T=[dependent.state_transform dependent.input_transform]; joint=[P raw.state_noise_cov;raw.state_noise_cov.' raw.covariance];
close(dependent.covariance,T*joint*T.',1e-12);
close(dependent.state_noise_cov,joint(1:42,:)*T.',1e-12);
bad=raw; bad.state_noise_cov(1,1)=1;
fails(@()nav_inertial_camera_motion({first,second},P,bad,environment,validity),'libekf:Covariance');
bad=raw; bad.accel_unit='g';
fails(@()nav_inertial_camera_motion({first,second},P,bad,environment,validity),'libekf:Units');
bad=raw; bad.quality='FAILED';
fails(@()nav_inertial_camera_motion({first,second},P,bad,environment,validity),'libekf:Quality');
[a,e]=fixture; a.p=zeros(3,1); a.q=[1;0;0;0]; a.bg=zeros(3,1); a.v=[2;1;.1];
b=a; b.t=a.t+.2; b.p=[.5;.1;.02]; b.v=[3;-.4;.2]; b.q=[cos(.04);0;sin(.04);0];
rates=[.1 -.2;.2 .1;0 .05]; K=[250 0 160;0 260 120;0 0 1];
D=nav_rotation([cos(.03);0;0;sin(.03)]); lever=[.2;.05;0];
stream=RandStream('mt19937ar','Seed',41); n=12;
landmarks=[2*rand(stream,2,n)-1;4+4*rand(stream,1,n)];
z=zeros(n,6); states={a,b};
for k=1:n
    for j=1:2
        state=states{j}; offset=.018*(k-1)/(n-1)-.009;
        if j==2, offset=-offset; end
        rotation=nav_rotation(state.q)*expm(nav_skew(rates(:,j)*offset));
        position=state.p+state.v*offset+rotation*lever;
        camera_point=D.'*rotation.'*(landmarks(:,k)-position);
        pixel=K*camera_point; z(k,3*j-2:3*j)=[pixel(1:2).'/pixel(3),state.t+offset];
    end
end
motion=struct('body_rates',rates,'epochs',[a.t b.t],'unit','rad/s','reference','earth', ...
    'valid',true,'bias_included',true,'max_offset',.02,'min_baseline',.01);
m=struct('model','visual_odometry','mode','epipolar','source','visible_camera', ...
    'lineage',{{'row_frame_1','row_frame_2'}},'quality','VALID','noise_model','independent', ...
    'correlation_group','camera','frame_id',a.frame_id,'t_start',a.t,'t_end',b.t, ...
    't_received',b.t+.05,'z',z,'intrinsics',K,'camera_to_body',D,'lever',lever, ...
    'scale_valid',false,'calibration_valid',true,'matches_valid',true,'scene_static',true,'row_motion',motion);
m.R=diag([repmat([.04;.04;1e-6;.04;.04;1e-6],n,1);ones(6,1)*1e-4]);
m.units=[repmat({'pixel';'pixel';'s';'pixel';'pixel';'s'},n,1);repmat({'rad/s'},6,1)];
o=nav_odometry_measurement(a,b,m); close(o.residual,zeros(n,1),1e-12);
% LEKF-025: ускорение внутри кадра не должно создавать ложную геометрическую невязку.
accelerated=m; acceleration=[20 -15;-10 12;5 -4];
accelerated.row_motion.acceleration=acceleration;
accelerated.row_motion.acceleration_unit='m/s^2';
accelerated.row_motion.acceleration_frame_id=a.frame_id;
accelerated.row_motion.acceleration_valid=true;
accelerated.R=blkdiag(m.R,.1*eye(6));
accelerated.units=[m.units;repmat({'m/s^2'},6,1)];
for k=1:n
    for j=1:2
        state=states{j}; offset=m.z(k,3*j)-state.t;
        rotation=nav_rotation(state.q)*expm(nav_skew(rates(:,j)*offset));
        position=state.p+state.v*offset+.5*acceleration(:,j)*offset^2+rotation*lever;
        camera_point=D.'*rotation.'*(landmarks(:,k)-position);
        pixel=K*camera_point; accelerated.z(k,3*j-2:3*j-1)=pixel(1:2).'/pixel(3);
    end
end
accelerated_observation=nav_odometry_measurement(a,b,accelerated);
close(accelerated_observation.residual,zeros(n,1),1e-12);
constant=m; constant.z=accelerated.z;
assert(norm(nav_odometry_measurement(a,b,constant).residual)>1e-7);
for j=[3 6 6*n+7:6*n+12]
    plus=accelerated; minus=accelerated;
    if j<=6
        plus.z(1,j)=plus.z(1,j)+1e-6; minus.z(1,j)=minus.z(1,j)-1e-6;
    else
        index=j-6*n-6;
        plus.row_motion.acceleration(index)=plus.row_motion.acceleration(index)+1e-4;
        minus.row_motion.acceleration(index)=minus.row_motion.acceleration(index)-1e-4;
    end
    step=1e-6; if j>6, step=1e-4; end
    op=nav_odometry_measurement(a,b,plus); om=nav_odometry_measurement(a,b,minus);
    close((op.residual-om.residual)/(2*step),accelerated_observation.noise_transform(:,j),1e-8);
end
% Ускорение и его погрешность меняют оси вместе с состоянием.
rotated=e; rotated.frame_id='accelerated_local';
rotated.C_e_l=nav_rotation([.9;.1;.2;.3]/norm([.9;.1;.2;.3]));
rotated.origin_ecef=e.origin_ecef+[10;20;30];
[changed,~]=nav_change_frame({a,b},eye(42),e,rotated);
T=rotated.C_e_l.'*e.C_e_l; input_rotation=blkdiag(eye(6*n+6),T,T);
transformed=accelerated; transformed.frame_id=rotated.frame_id;
transformed.row_motion.acceleration=T*acceleration;
transformed.row_motion.acceleration_frame_id=rotated.frame_id;
transformed.R=input_rotation*accelerated.R*input_rotation.';
transformed_observation=nav_odometry_measurement(changed{1},changed{2},transformed);
close(transformed_observation.residual,accelerated_observation.residual,1e-12);
close(transformed_observation.R,accelerated_observation.R,1e-12);
close(transformed_observation.noise_transform*input_rotation,accelerated_observation.noise_transform,1e-10);
bad=accelerated; bad.row_motion.acceleration_unit='g';
fails(@()nav_odometry_measurement(a,b,bad),'libekf:Acceleration');
bad=accelerated; bad.row_motion.acceleration_valid=false;
fails(@()nav_odometry_measurement(a,b,bad),'libekf:Acceleration');

assert(norm(o.H(:,4:6),'fro')>0&&norm(o.H(:,13:15),'fro')>0);
numeric=zeros(n,42);
for j=1:42
    delta=zeros(21,1); index=mod(j-1,21)+1; delta(index)=1e-6;
    if j<=21
        plus=nav_odometry_measurement(nav_apply_error(a,delta),b,m);
        minus=nav_odometry_measurement(nav_apply_error(a,-delta),b,m);
    else
        plus=nav_odometry_measurement(a,nav_apply_error(b,delta),m);
        minus=nav_odometry_measurement(a,nav_apply_error(b,-delta),m);
    end
    numeric(:,j)=-(plus.residual-minus.residual)/(2e-6);
end
close(numeric,o.H,1e-8);
for j=[1:6 6*n+1 6*n+5]
    plus=m; minus=m;
    if j<=6*n
        row=ceil(j/6); column=mod(j-1,6)+1;
        plus.z(row,column)=plus.z(row,column)+1e-6; minus.z(row,column)=minus.z(row,column)-1e-6;
    else
        plus.row_motion.body_rates(j-6*n)=plus.row_motion.body_rates(j-6*n)+1e-6;
        minus.row_motion.body_rates(j-6*n)=minus.row_motion.body_rates(j-6*n)-1e-6;
    end
    op=nav_odometry_measurement(a,b,plus); om=nav_odometry_measurement(a,b,minus);
    close((op.residual-om.residual)/(2e-6),o.noise_transform(:,j),1e-8);
end
close(o.R,o.noise_transform*m.R*o.noise_transform.',1e-12);
shifted_a=a; shifted_b=b; shifted_a.p=a.p+[100;200;-30]; shifted_b.p=b.p+[100;200;-30];
close(nav_odometry_measurement(shifted_a,shifted_b,m).residual,o.residual,1e-12);
dependent=m; dependent.noise_model='joint_required'; dependent.state_noise_frame_id=a.frame_id;
dependent.state_noise_cov=zeros(42,size(m.R,1)); dependent.state_noise_cov(1,1)=.001;
correlated=nav_odometry_measurement(a,b,dependent);
close(correlated.state_noise_cov,dependent.state_noise_cov*o.noise_transform.',1e-12);
prior=b; prior.p=prior.p+[.1;-.05;0];
filter=struct('state',prior,'P',eye(21),'reference',a,'jointP',eye(42),'used',{{}},'capacity',20);
moving=accelerated; moving.noise_model='joint_required'; moving.state_noise_frame_id=a.frame_id;
moving.state_noise_cov=zeros(42,size(moving.R,1)); moving.state_noise_cov(1,end)=.001;
moving_before=nav_odometry_measurement(a,prior,moving);
[moving_result,moving_info]=nav_assimilate(filter,moving,[],e,1000);
assert(moving_info.accepted);
moving_after=nav_odometry_measurement(moving_result.reference,moving_result.state,moving);
assert(norm(moving_after.residual)<norm(moving_before.residual)&&trace(moving_result.jointP)<trace(filter.jointP));

% Сквозной путь от отдельных новых инерциальных отсчётов до поправки по камере.
camera_env=e; camera_env.Omega=zeros(3,1); camera_env.gravity_origin=zeros(3,1);
camera_env.gravity_gradient=zeros(3); camera_env.gravity_radius=1e6;
camera_raw=raw; camera_raw.accel=[nav_rotation(a.q).'*acceleration(:,1)+a.ba, ...
    nav_rotation(prior.q).'*acceleration(:,2)+prior.ba];
camera_raw.gyro=rates; camera_raw.epochs=[a.t prior.t]; camera_raw.t_received=prior.t;
camera_raw.lineage_times=camera_raw.epochs; camera_raw.state_noise_frame_id=a.frame_id;
camera_raw.state_noise_cov=zeros(42,12);
prepared_motion=nav_inertial_camera_motion({a,prior},filter.jointP,camera_raw,camera_env,validity);
prepared_packet=accelerated; prepared_packet.row_motion=prepared_motion.motion;
prepared_packet.R=blkdiag(accelerated.R(1:6*n,1:6*n),prepared_motion.covariance);
prepared_packet.state_noise_cov=[zeros(42,6*n) prepared_motion.state_noise_cov];
prepared_packet.state_noise_frame_id=a.frame_id; prepared_packet.noise_model='joint_required';
prepared_packet.lineage=[accelerated.lineage(:);prepared_motion.lineage(:)];
prepared_packet.input_sources=prepared_motion.input_sources;
prepared_before=nav_odometry_measurement(a,prior,prepared_packet);
[prepared_result,prepared_info]=nav_assimilate(filter,prepared_packet,[],camera_env,1000);
assert(prepared_info.accepted);
prepared_after=nav_odometry_measurement(prepared_result.reference,prepared_result.state,prepared_packet);
assert(norm(prepared_after.residual)<norm(prepared_before.residual));
% Повторно использованные инерциальные отсчёты не скрываются внутри камеры.
already=filter; already.used=prepared_motion.lineage;
[unchanged,rejected]=nav_assimilate(already,prepared_packet,[],camera_env,1000);
assert(~rejected.accepted&&isequal(unchanged,already));
before=nav_odometry_measurement(a,prior,m);
[updated,info]=nav_assimilate(filter,m,[],e,1000); assert(info.accepted);
after=nav_odometry_measurement(updated.reference,updated.state,m);
assert(norm(after.residual)<norm(before.residual)&&trace(updated.jointP)<trace(filter.jointP));
[unchanged,rejected]=nav_assimilate(updated,m,[],e,1000);
assert(~rejected.accepted&&isequal(unchanged,updated));
bad=m; bad.z(1,3)=a.t+.03;
fails(@()nav_odometry_measurement(a,b,bad),'libekf:Time');
bad=m; bad.scene_static=false;
fails(@()nav_odometry_measurement(a,b,bad),'libekf:Quality');
bad=m; bad.units{3}='ms';
fails(@()nav_odometry_measurement(a,b,bad),'libekf:Units');
bad=m; bad.row_motion.bias_included=false;
fails(@()nav_odometry_measurement(a,b,bad),'libekf:Time');
zero_a=a; zero_b=b; zero_a.v=zeros(3,1); zero_b.v=zeros(3,1); zero_b.p=zero_a.p;
bad=m; bad.lever=zeros(3,1);
fails(@()nav_odometry_measurement(zero_a,zero_b,bad),'libekf:Geometry');
[dependent_update,accepted]=nav_assimilate(filter,dependent,[],e,1000); assert(accepted.accepted);
% Невозможная исходная зависимость отклоняется, даже если исчезает в проекции на невязку.
invisible=null(before.noise_transform); bad=dependent;
bad.state_noise_cov=zeros(size(dependent.state_noise_cov)); bad.state_noise_cov(1,:)=100*invisible(:,1).';
assert(norm(bad.state_noise_cov*before.noise_transform.','fro')<1e-10);
[unchanged,rejected]=nav_assimilate(filter,bad,[],e,1000);
assert(~rejected.accepted&&strcmp(rejected.reason,'libekf:Correlation')&&isequal(unchanged,filter));
inconsistent=nav_odometry_measurement(a,prior,dependent); inconsistent.R=2*inconsistent.R;
fails(@()nav_correct({a,prior},eye(42),inconsistent,1000),'libekf:Correlation');

% Полная поправка сохраняется при смене местного начала и осей, включая общие ошибки.
to=e; to.frame_id='timed_camera_other_axes';
to.C_e_l=e.C_e_l*nav_rotation([.9;.1;-.2;.3]/norm([.9;.1;-.2;.3]));
to.origin_ecef=e.origin_ecef+[300;-200;50];
[states,Pnew,J]=nav_change_frame({a,prior},eye(42),e,to);
rotated_filter=filter; rotated_filter.reference=states{1}; rotated_filter.state=states{2};
rotated_filter.jointP=Pnew; rotated_filter.P=Pnew(22:42,22:42);
for correlated_case=[false true]
    packet=m; expected=updated;
    if correlated_case, packet=dependent; expected=dependent_update; end
    old_observation=nav_odometry_measurement(a,prior,packet);
    packet.frame_id=to.frame_id;
    if correlated_case
        packet.state_noise_frame_id=to.frame_id; packet.state_noise_cov=J*packet.state_noise_cov;
    end
    new_observation=nav_odometry_measurement(states{1},states{2},packet);
    close(new_observation.residual,old_observation.residual,1e-11);
    close(new_observation.H*J,old_observation.H,1e-10);
    close(new_observation.R,old_observation.R,1e-12);
    close(new_observation.noise_transform,old_observation.noise_transform,1e-10);
    [actual,accepted]=nav_assimilate(rotated_filter,packet,[],to,1000); assert(accepted.accepted);
    [converted,covariance]=nav_change_frame({expected.reference,expected.state},expected.jointP,e,to);
    close(actual.reference.p,converted{1}.p,1e-10); close(actual.state.p,converted{2}.p,1e-10);
    close(actual.state.v,converted{2}.v,1e-10);
    close(nav_rotation(actual.state.q),nav_rotation(converted{2}.q),1e-10);
    close(actual.jointP,covariance,1e-10);
end

end

function [out,result]=nav_assimilate_interval(filter,imu,packets,sensors,env,limit,reception)
% Наблюдения внутри одного интервала постоянных показаний ИНС.
% Плотность белого шума относится к непрерывному времени, не к среднему отсчёту.
assert(iscell(packets)&&iscell(sensors)&&numel(packets)==numel(sensors), ...
    'libekf:Packet','Expected matching observation and sensor cells.');
out=filter; result=struct('accepted',false,'reason','','observations',{{}},'steps',{{}});
args={}; if nargin>=7, args={reception}; end
% max_age относится к ИНС; observation_max_age — к исходному порядку наблюдений.
if nargin>=7&&isstruct(reception)&&isfield(reception,'observation_max_age')
    ages=reception.observation_max_age;
    assert(isnumeric(ages)&&isreal(ages)&&numel(ages)==numel(packets)&& ...
        all(isfinite(ages(:)))&&all(ages(:)>=0), ...
        'libekf:Time','One finite nonnegative age limit per observation required.');
end
assert(isfield(imu,'model')&&strcmp(imu.model,'inertial'), ...
    'libekf:Packet','Expected one inertial interval.');
% Проверка исходного интервала, его происхождения и сохранённой опоры.
[whole,check]=nav_assimilate(filter,imu,[],env,limit,args{:});
if ~check.accepted, result.reason=check.reason; return; end
if isempty(packets)
    out=whole; result.accepted=true; result.reason=check.reason;
    result.steps={struct('kind','propagation','data',check.propagation)}; return;
end
if isfield(filter,'future_noise')&&~isequal(imu.future_noise_independent,true)
    result.reason='libekf:Correlation'; return; % Для разбиения нужна временная модель общего шума.
end
if isfield(filter,'future_noise')&&isfield(filter.future_noise,'peer_cross_cov')&& ...
        ~isequal(imu.reference_motion.future_noise_independent,true)
    result.reason='libekf:Correlation'; return; % Связь за целый интервал нельзя копировать на его части.
end
if isfield(filter,'reference_bank')&&~isfield(imu.reference_motion,'F')
    result.reason='libekf:ReferenceInterval'; return;
end
steps={};
candidate=filter; candidate.used=whole.used; % Исходные данные ИНС учитываются ровно один раз.
n=numel(packets); epochs=Inf(1,n); checks=cell(size(packets));
for k=1:n
    checks{k}=struct('accepted',false,'reason','libekf:Epoch','nis',NaN);
    p=packets{k};
    if ~isstruct(p)||~isscalar(p)||~isfield(p,'model'), continue; end
    if strcmp(p.model,'inertial'), continue; end
    field='t_measure'; if strcmp(p.model,'visual_odometry'), field='t_end'; end
    if ~isfield(p,field), continue; end
    t=p.(field);
    if ~isscalar(t)||~isreal(t)||~isfinite(t)||t<imu.t_start||t>imu.t_end, continue; end
    epochs(k)=t;
end
[~,order]=sort(epochs);
try
    for k=order
        if ~isfinite(epochs(k)), continue; end
        [trial,propagation]=advance(candidate,imu,epochs(k),env);
        if strcmp(packets{k}.model,'save_reference')
            % Явное действие истории: снимок не добавляет измерительной информации.
            [trial,reference_change]=nav_save_reference(trial);
            checks{k}=struct('accepted',true,'reason','Опора сохранена','nis',NaN, ...
                'reference_change',reference_change);
        else
            observation_args=args;
            if nargin>=7&&isfield(reception,'observation_max_age')
                policy=reception; policy.max_age=reception.observation_max_age(k);
                observation_args={policy};
            end
            [trial,checks{k}]=nav_assimilate(trial,packets{k},sensors{k},env,limit,observation_args{:});
        end
        % Отвергнутое наблюдение не меняет даже разбиение интегрирования ИНС.
        if checks{k}.accepted
            if ~isempty(propagation)
                steps{end+1}=struct('kind','propagation','data',propagation);
            end
            steps{end+1}=struct('kind',packets{k}.model,'packet_index',k, ...
                'time',epochs(k),'data',checks{k});
            candidate=trial;
        end
    end
    [candidate,propagation]=advance(candidate,imu,imu.t_end,env);
    if ~isempty(propagation), steps{end+1}=struct('kind','propagation','data',propagation); end
    source_history=[]; if isfield(candidate,'source_history'), source_history=candidate.source_history; end
    candidate.source_history=nav_sources(source_history,imu);
    out=candidate; result.accepted=true; result.reason='Интервал обработан по времени измерений';
    result.observations=checks; result.steps=steps;
catch ex
    if startsWith(ex.identifier,'libekf:')
        result.reason=ex.identifier;
    else
        rethrow(ex);
    end
end
end

function [out,record]=advance(in,imu,t,env)
out=in; record=[];
if t==in.state.t, return; end
assert(~isfield(in,'excluded_sources')||~any(strcmp(in.excluded_sources,imu.source)), ...
    'libekf:SourceHealth','Cannot propagate with an excluded inertial source.');
if isfield(in,'excluded_sources')&&isfield(imu,'input_sources')
    assert(~any(ismember(imu.input_sources(:),in.excluded_sources(:))), ...
        'libekf:SourceHealth','Inertial data depend on an excluded source.');
end
part=imu; part.t_start=in.state.t; part.t_end=t;
[out.state,out.P,T,N]=nav_predict(in.state,in.P,part,env);
if isfield(in,'reference')
    cross=in.jointP(1:21,22:42)*T.';
    out.jointP=[in.jointP(1:21,1:21) cross;cross.' out.P];
end
transition=T; process_covariance=N;
if isfield(in,'reference')
    transition=blkdiag(eye(21),T); process_covariance=blkdiag(zeros(21),N);
end
if isfield(in,'future_noise')
    assert(isfield(imu,'future_noise_independent')&&isequal(imu.future_noise_independent,true)&& ...
        (~isfield(imu,'future_noise_process_cross_cov')||all(imu.future_noise_process_cross_cov(:)==0)), ...
        'libekf:Correlation','Splitting correlated process/stored noise requires a temporal model.');
    prior=in.P; epochs=out.state.t;
    if isfield(in,'reference'), prior=in.jointP; epochs=[in.reference.t;out.state.t]; end
    if isfield(in.future_noise,'peer_cross_cov')
        assert(isfield(in,'reference_bank')&&isfield(imu,'reference_motion'),'libekf:Correlation', ...
            'Joint errors require the retained reference motion model.');
        joint_motion=imu.reference_motion; joint_motion.t_start=in.state.t; joint_motion.t_end=t;
        [out.future_noise,out.reference_bank]=nav_predict_noise(in.future_noise,prior,transition,process_covariance, ...
            epochs,zeros(size(transition,1),numel(in.future_noise.mean)),in.reference_bank,joint_motion);
    else
        out.future_noise=nav_predict_noise(in.future_noise,prior,transition,process_covariance, ...
            epochs,zeros(size(transition,1),numel(in.future_noise.mean)));
    end
end
record=struct('transition',transition,'process_covariance',process_covariance, ...
    't_start',part.t_start,'t_end',part.t_end,'frame_id',in.state.frame_id, ...
    'source',imu.source,'lineage',{imu.lineage});
if isfield(in,'future_noise')
    record.future_noise_process_cross_cov=zeros(size(transition,1),numel(in.future_noise.mean));
end
if isfield(in,'reference_bank')&& ...
        ~(isfield(in,'future_noise')&&isfield(in.future_noise,'peer_cross_cov'))
    motion=imu.reference_motion; motion.t_start=in.state.t; motion.t_end=t;
    transition=T; if isfield(in,'reference'), transition=blkdiag(eye(21),T); end
    out.reference_bank=nav_predict_reference(in.reference_bank,transition,motion);
end
end

function [out,result]=nav_assimilate_interval(filter,imu,packets,sensors,env,limit,reception)
% Наблюдения внутри одного интервала постоянных показаний ИНС.
% Плотность белого шума относится к непрерывному времени, не к среднему отсчёту.
assert(iscell(packets)&&iscell(sensors)&&numel(packets)==numel(sensors), ...
    'libekf:Packet','Expected matching observation and sensor cells.');
out=filter; result=struct('accepted',false,'reason','','observations',{{}});
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
if isempty(packets), out=whole; result.accepted=true; result.reason=check.reason; return; end
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
        trial=advance(candidate,imu,epochs(k),env);
        if strcmp(packets{k}.model,'save_reference')
            % Явное действие истории: снимок не добавляет измерительной информации.
            trial=nav_save_reference(trial);
            checks{k}=struct('accepted',true,'reason','Опора сохранена','nis',NaN);
        else
            observation_args=args;
            if nargin>=7&&isfield(reception,'observation_max_age')
                policy=reception; policy.max_age=reception.observation_max_age(k);
                observation_args={policy};
            end
            [trial,checks{k}]=nav_assimilate(trial,packets{k},sensors{k},env,limit,observation_args{:});
        end
        % Отвергнутое наблюдение не меняет даже разбиение интегрирования ИНС.
        if checks{k}.accepted, candidate=trial; end
    end
    candidate=advance(candidate,imu,imu.t_end,env);
    out=candidate; result.accepted=true; result.reason='Интервал обработан по времени измерений';
    result.observations=checks;
catch ex
    if startsWith(ex.identifier,'libekf:')
        result.reason=ex.identifier;
    else
        rethrow(ex);
    end
end
end

function out=advance(in,imu,t,env)
out=in;
if t==in.state.t, return; end
part=imu; part.t_start=in.state.t; part.t_end=t;
[out.state,out.P,T]=nav_predict(in.state,in.P,part,env);
if isfield(in,'reference')
    cross=in.jointP(1:21,22:42)*T.';
    out.jointP=[in.jointP(1:21,1:21) cross;cross.' out.P];
end
end

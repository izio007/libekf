function [out,result]=nav_replay(initial,intervals,packets,sensors,env,limit,admission)
% Повторный расчёт сохранённой истории от неизменённого исходного снимка.
% admission{k}: now/max_age для нового отсчёта; [] для ранее допущенной истории.
% Без admission вызывающая сторона сама проверяет получение и возраст новых данных.
% Это не добавление старых поправок к текущей оценке; исходная история не очищается.
assert(iscell(intervals)&&~isempty(intervals)&&iscell(packets)&&iscell(sensors)&& ...
    numel(packets)==numel(sensors),'libekf:History','Expected inertial history and paired observation cells.');
out=initial; result=struct('accepted',false,'reason','','observations',{{}});
n=numel(packets); epochs=Inf(1,n); checks=cell(size(packets));
if nargin<7, admission=cell(size(packets)); end
assert(iscell(admission)&&numel(admission)==n, ...
    'libekf:Time','One admission policy per observation required.');
for k=1:n
    checks{k}=struct('accepted',false,'reason','libekf:Epoch','nis',NaN);
    p=packets{k};
    if ~isstruct(p)||~isscalar(p)||~isfield(p,'model')||strcmp(p.model,'inertial'), continue; end
    field='t_measure'; if strcmp(p.model,'visual_odometry'), field='t_end'; end
    if ~isfield(p,field), continue; end
    t=p.(field);
    if ~isscalar(t)||~isreal(t)||~isfinite(t), continue; end
    policy=admission{k};
    if ~isempty(policy)
        valid=isstruct(policy)&&isscalar(policy)&&all(isfield(policy,{'now','max_age'}));
        if valid
            valid=isnumeric(policy.now)&&isreal(policy.now)&&isscalar(policy.now)&& ...
                isfinite(policy.now)&&isnumeric(policy.max_age)&&isreal(policy.max_age)&& ...
                isscalar(policy.max_age)&&isfinite(policy.max_age)&&policy.max_age>=0;
        end
        if ~valid, checks{k}.reason='libekf:Time'; continue; end
        if strcmp(p.model,'save_reference')
            % Действие истории не имеет времени получения измерения.
            if t>policy.now, checks{k}.reason='libekf:Time'; continue; end
        else
            if ~isfield(p,'t_received')||~isnumeric(p.t_received)|| ...
                ~isscalar(p.t_received)||~isreal(p.t_received)||~isfinite(p.t_received)|| ...
                t>p.t_received||p.t_received>policy.now
                checks{k}.reason='libekf:Time'; continue;
            end
            if policy.now-t>policy.max_age
                checks{k}.reason='libekf:Stale'; continue;
            end
        end
    end
    epochs(k)=t;
end
candidate=initial; taken=false(1,n);
for k=1:numel(intervals)
    u=intervals{k};
    assert(isstruct(u)&&isscalar(u)&&all(isfield(u,{'t_start','t_end'}))&& ...
        isscalar(u.t_start)&&isscalar(u.t_end)&&all(isfinite([u.t_start u.t_end])), ...
        'libekf:History','Invalid inertial interval in history.');
    selected=find(~taken&epochs>=u.t_start&epochs<=u.t_end);
    [candidate,step]=nav_assimilate_interval(candidate,u,packets(selected),sensors(selected),env,limit);
    if ~step.accepted, result.reason=step.reason; return; end
    checks(selected)=step.observations; taken(selected)=true;
end
out=candidate; result.accepted=true; result.reason='История пересчитана от исходного снимка';
result.observations=checks;
end

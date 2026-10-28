function [out,result]=nav_replay(initial,intervals,packets,sensors,env,limit)
% Повторный расчёт сохранённой истории от неизменённого исходного снимка.
% Вызывающая сторона передаёт только уже полученные и допущенные по возрасту данные.
% Это не добавление старых поправок к текущей оценке; исходная история не очищается.
assert(iscell(intervals)&&~isempty(intervals)&&iscell(packets)&&iscell(sensors)&& ...
    numel(packets)==numel(sensors),'libekf:History','Expected inertial history and paired observation cells.');
out=initial; result=struct('accepted',false,'reason','','observations',{{}});
n=numel(packets); epochs=Inf(1,n); checks=cell(size(packets));
for k=1:n
    checks{k}=struct('accepted',false,'reason','libekf:Epoch','nis',NaN);
    p=packets{k};
    if ~isstruct(p)||~isscalar(p)||~isfield(p,'model')||strcmp(p.model,'inertial'), continue; end
    field='t_measure'; if strcmp(p.model,'visual_odometry'), field='t_end'; end
    if ~isfield(p,field), continue; end
    t=p.(field);
    if isscalar(t)&&isreal(t)&&isfinite(t), epochs(k)=t; end
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

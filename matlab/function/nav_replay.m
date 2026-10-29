function [out,result]=nav_replay(initial,intervals,packets,sensors,env,limit,admission,recovery_exclusions)
% Повторный расчёт сохранённой истории от неизменённого исходного снимка.
% admission{k}: now/max_age для нового отсчёта; [] для ранее допущенной истории.
% Без admission вызывающая сторона сама проверяет получение и возраст новых данных.
% Это не добавление старых поправок к текущей оценке; исходная история не очищается.
assert(iscell(intervals)&&~isempty(intervals)&&iscell(packets)&&iscell(sensors)&& ...
    numel(packets)==numel(sensors),'libekf:History','Expected inertial history and paired observation cells.');
out=initial; result=struct('accepted',false,'reason','','observations',{{}},'steps',{{}});
% Восстановительный пересчёт: запрет действует на всю историю, включая сообщения возврата.
if nargin<8, recovery_exclusions={}; end
assert(iscell(recovery_exclusions)&&all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x), ...
    recovery_exclusions)),'libekf:SourceHealth','Invalid recovery exclusion list.');
if ~isempty(recovery_exclusions)
    if ~isfield(initial,'sources_complete')||~isequal(initial.sources_complete,true)|| ...
            ~isfield(initial,'source_history')|| ...
            (isfield(initial,'invalidated_sources')&&~isempty(initial.invalidated_sources))
        result.reason='libekf:RecoveryInitial'; return;
    end
    history=initial.source_history;
    if ~isstruct(history)||~isfield(history,'id')|| ...
            any(ismember({history.id},recovery_exclusions))
        result.reason='libekf:RecoveryInitial'; return;
    end
    previous={}; if isfield(initial,'excluded_sources'), previous=initial.excluded_sources; end
    assert(iscell(previous)&&all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),previous)), ...
        'libekf:SourceHealth','Invalid snapshot exclusion list.');
    initial.excluded_sources=unique([previous(:);recovery_exclusions(:)],'stable').';
end
n=numel(packets); epochs=Inf(1,n); checks=cell(size(packets));
if nargin<7, admission=cell(size(packets)); end
assert(iscell(admission)&&numel(admission)==n, ...
    'libekf:Time','One admission policy per observation required.');
for k=1:n
    checks{k}=struct('accepted',false,'reason','libekf:Epoch','nis',NaN);
    p=packets{k};
    if ~isstruct(p)||~isscalar(p)||~isfield(p,'model')||strcmp(p.model,'inertial'), continue; end
    if strcmp(p.model,'source_health')&&isfield(p,'target_source')&& ...
            isfield(p,'valid')&&isequal(p.valid,true)&& ...
            any(strcmp(recovery_exclusions,p.target_source))
        checks{k}.reason='libekf:SourceHealth'; continue;
    end
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
% Один контекст на всю историю либо явный контекст каждого интервала.
environments=env; if ~iscell(env), environments=repmat({env},size(intervals)); end
assert(iscell(environments)&&numel(environments)==numel(intervals), ...
    'libekf:Frame','One environment per inertial interval required.');
for j=1:numel(environments)
    frame=environments{j};
    assert(isstruct(frame)&&isscalar(frame)&&all(isfield(frame,{'frame_id','C_e_l','origin_ecef'}))&& ...
        ischar(frame.frame_id)&&isrow(frame.frame_id)&&~isempty(frame.frame_id)&& ...
        isnumeric(frame.origin_ecef)&&isreal(frame.origin_ecef)&&numel(frame.origin_ecef)==3&& ...
        all(isfinite(frame.origin_ecef(:))),'libekf:Frame','Invalid interval frame.');
    nav_quaternion(frame.C_e_l);
    for previous=1:j-1
        known=environments{previous};
        assert(~strcmp(known.frame_id,frame.frame_id)|| ...
            (isequal(known.C_e_l,frame.C_e_l)&&isequal(known.origin_ecef(:),frame.origin_ecef(:))), ...
            'libekf:Frame','One frame identity cannot denote different origins or axes.');
    end
end
active=environments{1};
if ~strcmp(initial.state.frame_id,active.frame_id), result.reason='libekf:Frame'; return; end
candidate=initial; taken=false(1,n); steps={};
for k=1:numel(intervals)
    u=intervals{k};
    assert(isstruct(u)&&isscalar(u)&&all(isfield(u,{'t_start','t_end'}))&& ...
        isscalar(u.t_start)&&isscalar(u.t_end)&&all(isfinite([u.t_start u.t_end])), ...
        'libekf:History','Invalid inertial interval in history.');
    next=environments{k};
    if ~strcmp(active.frame_id,next.frame_id)
        times=candidate.state.t;
        if isfield(candidate,'reference'), times=[candidate.reference.t;candidate.state.t]; end
        try
            [candidate,~,transition]=nav_change_frame(candidate,[],active,next);
        catch ex
            if startsWith(ex.identifier,'libekf:'), result.reason=ex.identifier; return; else, rethrow(ex); end
        end
        change=struct('accepted',true,'transition',transition,'state_times',times, ...
            'from_frame_id',active.frame_id,'to_frame_id',next.frame_id);
        steps{end+1}=struct('kind','change_frame','time',candidate.state.t, ...
            'data',change,'interval_index',k);
    end
    active=next;
    selected=find(~taken&epochs>=u.t_start&epochs<=u.t_end);
    [candidate,step]=nav_assimilate_interval(candidate,u,packets(selected),sensors(selected),active,limit);
    if ~step.accepted, result.reason=step.reason; return; end
    for j=1:numel(step.steps)
        record=step.steps{j}; record.interval_index=k;
        if isfield(record,'packet_index'), record.packet_index=selected(record.packet_index); end
        steps{end+1}=record;
    end
    checks(selected)=step.observations; taken(selected)=true;
end
out=candidate; result.accepted=true; result.reason='История пересчитана от исходного снимка';
result.observations=checks; result.steps=steps;
end

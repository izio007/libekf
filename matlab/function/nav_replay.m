function [out,result]=nav_replay(initial,intervals,packets,sensors,env,limit,admission,recovery_exclusions,recovery_snapshots)
% Повторный расчёт сохранённой истории от неизменённого исходного снимка.
% admission{k}: now/max_age для нового отсчёта; [] для ранее допущенной истории.
% Без admission вызывающая сторона сама проверяет получение и возраст новых данных.
% Это не добавление старых поправок к текущей оценке; исходная история не очищается.
assert(iscell(intervals)&&~isempty(intervals)&&iscell(packets)&&iscell(sensors)&& ...
    numel(packets)==numel(sensors),'libekf:History','Expected inertial history and paired observation cells.');
out=initial; result=struct('accepted',false,'reason','','observations',{{}},'steps',{{}});
n=numel(packets); epochs=Inf(1,n); checks=cell(size(packets));
if nargin<7, admission=cell(size(packets)); end
assert(iscell(admission)&&numel(admission)==n, ...
    'libekf:Time','One admission policy per observation required.');
% Восстановительный пересчёт: запрет действует на всю историю, включая сообщения возврата.
if nargin<8, recovery_exclusions={}; end
assert(iscell(recovery_exclusions)&&all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x), ...
    recovery_exclusions)),'libekf:SourceHealth','Invalid recovery exclusion list.');
first_interval=1;
if nargin>=9
    assert(iscell(recovery_snapshots)&&~isempty(recovery_exclusions), ...
        'libekf:RecoveryInitial','Snapshot selection requires recovery exclusions.');
    required_epoch=Inf;
    for k=1:n
        p=packets{k};
        if isempty(admission{k})||~isstruct(p)||~isscalar(p)||~isfield(p,'model'), continue; end
        field='t_measure'; if strcmp(p.model,'visual_odometry'), field='t_end'; end
        if isfield(p,field)&&isnumeric(p.(field))&&isreal(p.(field))&& ...
                isscalar(p.(field))&&isfinite(p.(field))
            required_epoch=min(required_epoch,p.(field));
        end
    end
    selected=0; newest=-Inf;
    for j=1:numel(recovery_snapshots)
        saved=recovery_snapshots{j};
        if ~recovery_initial(saved,recovery_exclusions)||saved.state.t>required_epoch, continue; end
        boundary=find(cellfun(@(u)isstruct(u)&&isscalar(u)&&isfield(u,'t_start')&& ...
            isequal(u.t_start,saved.state.t),intervals),1);
        if ~isempty(boundary)&&saved.state.t>newest
            selected=j; newest=saved.state.t; first_interval=boundary;
        end
    end
    if selected==0, result.reason='libekf:RecoveryInitial'; return; end
    initial=recovery_snapshots{selected};
    result.snapshot_index=selected; result.snapshot_time=initial.state.t;
end
if ~isempty(recovery_exclusions)
    if ~recovery_initial(initial,recovery_exclusions)
        result.reason='libekf:RecoveryInitial'; return;
    end
    previous={}; if isfield(initial,'excluded_sources'), previous=initial.excluded_sources; end
    assert(iscell(previous)&&all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),previous)), ...
        'libekf:SourceHealth','Invalid snapshot exclusion list.');
    initial.excluded_sources=unique([previous(:);recovery_exclusions(:)],'stable').';
end
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
active=environments{first_interval};
if ~strcmp(initial.state.frame_id,active.frame_id), result.reason='libekf:Frame'; return; end
candidate=initial; taken=false(1,n); steps={};
for k=first_interval:numel(intervals)
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

function valid=recovery_initial(saved,excluded)
valid=false;
if ~isstruct(saved)||~isscalar(saved)||~isfield(saved,'state')|| ...
        ~isstruct(saved.state)||~isscalar(saved.state)||~isfield(saved.state,'t')|| ...
        ~isnumeric(saved.state.t)||~isreal(saved.state.t)||~isscalar(saved.state.t)|| ...
        ~isfinite(saved.state.t)||~isfield(saved,'sources_complete')|| ...
        ~isequal(saved.sources_complete,true)||~isfield(saved,'source_history')|| ...
        (isfield(saved,'invalidated_sources')&&~isempty(saved.invalidated_sources))
    return;
end
try
    history=nav_sources(saved.source_history);
catch ex
    if startsWith(ex.identifier,'libekf:'), return; else, rethrow(ex); end
end
valid=~any(ismember({history.id},excluded))&& ...
    all([history.last_epoch]<=saved.state.t);
end

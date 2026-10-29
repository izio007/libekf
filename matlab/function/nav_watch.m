function [history,result]=nav_watch(history,sample,now,cfg)
% Три интервала контроля одного сопоставимого ряда диагностических проверок.
% sample.value — расхождение, делённое на его погрешность; [] означает отсутствие нового отсчёта.
% context меняется при смене источников, осей или модели нормирования.
w=cfg.windows(:).';
assert(numel(w)==3 && all(isfinite(w)) && all(w>0) && all(diff(w)>0), ...
    'libekf:WatchConfig','Expected increasing F/M/S durations.');
assert(isscalar(cfg.max_gap)&&isfinite(cfg.max_gap)&&cfg.max_gap>0 && ...
    isscalar(cfg.max_age)&&isfinite(cfg.max_age)&&cfg.max_age>=0 && ...
    isscalar(cfg.capacity)&&cfg.capacity>=2&&isfinite(cfg.capacity)&&cfg.capacity==floor(cfg.capacity), ...
    'libekf:WatchConfig','Invalid history limits.');
limits=[cfg.instant_limit cfg.mean_limit(:).' cfg.rms_limit(:).'];
assert(numel(limits)==7 && all(isfinite(limits)) && all(limits>0) && isscalar(now)&&isfinite(now), ...
    'libekf:WatchConfig','Invalid thresholds or current time.');
if isempty(history), history=struct('t',[],'value',[],'context','','now',-Inf); end
if ~isfield(history,'alarm_only'), history.alarm_only=false; end
assert(now>=history.now,'libekf:Epoch','Watch time cannot go backwards.');
history.now=now;
if ~isempty(sample)
    assert(isscalar(sample.t)&&isfinite(sample.t)&&sample.t<=now && ...
        ischar(sample.context)&&~isempty(sample.context) && ...
        isvector(sample.value)&&~isempty(sample.value)&&all(isfinite(sample.value)), ...
        'libekf:WatchSample','Invalid diagnostic check.');
    assert(isempty(history.t)||sample.t>history.t(end),'libekf:Epoch','Repeated or late check.');
    alarm_only=false; if isfield(sample,'alarm_only'), alarm_only=sample.alarm_only; end
    assert(isequal(alarm_only,true)||isequal(alarm_only,false), ...
        'libekf:WatchSample','Invalid alarm-only check mode.');
    if ~strcmp(history.context,sample.context)
        history.t=[]; history.value=[]; history.context=sample.context; history.alarm_only=alarm_only;
    end
    assert(history.alarm_only==alarm_only,'libekf:WatchSample','Changed check purpose requires a new context.');
    assert(isempty(history.value)||size(history.value,1)==numel(sample.value), ...
        'libekf:WatchSample','Changed check dimensions require a new context.');
    history.t(end+1)=sample.t; history.value(:,end+1)=sample.value(:);
end
% Один предшествующий отсчёт нужен для проверки покрытия левой границы.
before=find(history.t<now-w(3),1,'last');
if ~isempty(before), history.t=history.t(before:end); history.value=history.value(:,before:end); end
if numel(history.t)>cfg.capacity
    history.t=history.t(end-cfg.capacity+1:end); history.value=history.value(:,end-cfg.capacity+1:end);
end
result=struct('status','NOT_READY','ready',false(1,3),'alarm',false(1,3), ...
    'mean',NaN(1,3),'rms',NaN(1,3),'instant_alarm',false,'t',now,'age',Inf, ...
    'alarm_only',history.alarm_only);
if isempty(history.t), return; end
result.age=now-history.t(end); fresh=result.age<=cfg.max_age;
result.instant_alarm=fresh && norm(history.value(:,end))>cfg.instant_limit;
for k=1:3
    left=now-w(k); first=find(history.t<=left,1,'last');
    ix=find(history.t>=left & history.t<=now);
    if isempty(ix), continue; end
    V=history.value(:,ix);
    result.mean(k)=norm(mean(V,2)); result.rms(k)=norm(V,'fro')/sqrt(numel(ix));
    if ~isempty(first)
        gaps=diff([history.t(first:end) now]);
        result.ready(k)=fresh && all(gaps<=cfg.max_gap);
    end
    result.alarm(k)=result.ready(k) && ...
        ((~history.alarm_only&&result.mean(k)>cfg.mean_limit(k))||result.rms(k)>cfg.rms_limit(k));
end
if result.instant_alarm||any(result.alarm)
    result.status='ALARM';
elseif all(result.ready)&&~history.alarm_only
    result.status='CONSISTENT';
end
end

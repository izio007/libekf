function result=nav_isolate_sources(sources,checks,now,cfg)
% Условная изоляция по всем допустимым наборам отказов; без голосования.
% Тревога требует хотя бы одной причины из involved_sources соответствующей проверки.
assert(iscell(sources)&&~isempty(sources)&&all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),sources))&& ...
    numel(unique(sources))==numel(sources)&&iscell(checks)&&isscalar(now)&&isreal(now)&&isfinite(now), ...
    'libekf:Isolation','Invalid physical source list, checks or current time.');
n=numel(sources); sources=sources(:);
assert(isstruct(cfg)&&isscalar(cfg)&&all(isfield(cfg,{'max_faults','max_hypotheses','max_age','max_span'}))&& ...
    isscalar(cfg.max_faults)&&isreal(cfg.max_faults)&&isfinite(cfg.max_faults)&& ...
    cfg.max_faults>=0&&cfg.max_faults<=n&&cfg.max_faults==floor(cfg.max_faults)&& ...
    isscalar(cfg.max_hypotheses)&&isreal(cfg.max_hypotheses)&&isfinite(cfg.max_hypotheses)&& ...
    cfg.max_hypotheses>=1&&cfg.max_hypotheses<=flintmax&&cfg.max_hypotheses==floor(cfg.max_hypotheses)&& ...
    isscalar(cfg.max_age)&&isreal(cfg.max_age)&&isfinite(cfg.max_age)&&cfg.max_age>=0&& ...
    isscalar(cfg.max_span)&&isreal(cfg.max_span)&&isfinite(cfg.max_span)&&cfg.max_span>=0, ...
    'libekf:Isolation','Explicit fault bound, search capacity and time limits required.');
result=struct('status','NOT_READY','reason','','suspected_sources',{{}}, ...
    'hypotheses',{{}},'complete',false,'candidate_count',NaN,'t',now,'used_checks',[],'interval',[NaN NaN]);
conflicts={}; periods=zeros(0,2);
for k=1:numel(checks)
    c=checks{k};
    assert(isstruct(c)&&isscalar(c)&&all(isfield(c,{'status','t'}))&& ...
        any(strcmp(c.status,{'ALARM','DISAGREEMENT','CONSISTENT','NOT_READY'})), ...
        'libekf:Isolation','Invalid diagnostic result.');
    if ~any(strcmp(c.status,{'ALARM','DISAGREEMENT'})), continue; end
    assert(isscalar(c.t)&&isreal(c.t)&&isfinite(c.t),'libekf:Isolation','Invalid alarm epoch.');
    if c.t>now||now-c.t>cfg.max_age, continue; end
    period=[c.t c.t]; if isfield(c,'interval'), period=c.interval; end
    assert(isnumeric(period)&&isreal(period)&&numel(period)==2&&all(isfinite(period(:)))&& ...
        period(1)<=period(2)&&period(2)<=c.t,'libekf:Isolation','Invalid diagnostic observation interval.');
    if now-period(2)>cfg.max_age, continue; end

    if ~isfield(c,'provenance_complete')||~isequal(c.provenance_complete,true)
        result.reason='Неизвестен полный состав причин одной из тревог'; return;
    end
    assert(isfield(c,'involved_sources')&&iscell(c.involved_sources)&&~isempty(c.involved_sources)&& ...
        all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),c.involved_sources))&& ...
        all(ismember(c.involved_sources,sources)), ...
        'libekf:Isolation','Alarm refers to an unknown physical source.');
    conflicts{end+1}=find(ismember(sources,c.involved_sources));
    periods(end+1,:)=period(:).'; result.used_checks(end+1)=k;
end
if isempty(conflicts), result.reason='Нет пригодных тревог для изоляции'; return; end
result.interval=[min(periods(:,1)) max(periods(:,2))];
if diff(result.interval)>cfg.max_span
    result.reason='Тревоги относятся к несогласованным временным интервалам'; return;
end
count=1; total=1;
for k=1:cfg.max_faults
    count=round(count*(n-k+1)/k); total=total+count;
    if total>cfg.max_hypotheses
        result.reason='Недостаточно памяти для полного перебора; источник не выделен'; return;
    end
end
feasible={};
for faults=0:cfg.max_faults
    if faults==0, choices=zeros(1,0); else, choices=nchoosek(1:n,faults); end
    for row=1:size(choices,1)
        selected=choices(row,:);
        if all(cellfun(@(edge)any(ismember(edge,selected)),conflicts))
            feasible{end+1}=sources(selected);
        end
    end
end
result.complete=true; result.candidate_count=numel(feasible); result.hypotheses=feasible;
if isempty(feasible)
    result.status='INCONSISTENT'; result.reason='Тревоги противоречат заданной модели числа отказов'; return;
end
forced=feasible{1};
for k=2:numel(feasible), forced=intersect(forced,feasible{k}); end
result.suspected_sources=sources(ismember(sources,forced));
result.status='AMBIGUOUS'; result.reason='Допустимы разные причины; источник не выделен';
if ~isempty(forced)
    result.status='ISOLATED'; result.reason='Источники входят во все объяснения в рамках заданной модели отказов';
end
end

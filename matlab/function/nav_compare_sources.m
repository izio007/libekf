function [sample,result]=nav_compare_sources(a,b,cross_cov,now,cfg)
% Сравнение двух оценок одной физической величины; не определяет виновный источник.
sample=[]; result=struct('status','NOT_READY','reason','','nis',NaN, ...
    'shared_sources',{{}},'involved_sources',{{}},'provenance_complete',false,'t',NaN,'interval',[NaN NaN], ...
    'provenance','UNKNOWN','covariance_model','UNSPECIFIED','nis_kind','UNAVAILABLE');
required={'value','covariance','quantity','units','frame_id','t','quality','source'};
for k=1:numel(required)
    assert(isfield(a,required{k})&&isfield(b,required{k}),'libekf:Comparison','Incomplete source estimate.');
end
assert(isscalar(now)&&isfinite(now)&&isscalar(cfg.max_age)&&isfinite(cfg.max_age)&&cfg.max_age>=0&& ...
    isscalar(cfg.limit)&&isfinite(cfg.limit)&&cfg.limit>0&&ischar(cfg.context)&&~isempty(cfg.context), ...
    'libekf:Comparison','Explicit age, threshold and diagnostic context required.');
assert(isscalar(a.t)&&isscalar(b.t)&&all(isfinite([a.t b.t])), ...
    'libekf:Comparison','Invalid comparison epochs.');
if ~strcmp(a.quality,'VALID')||~strcmp(b.quality,'VALID')
    result.reason='Источник сообщил непригодность'; return;
end
if a.t~=b.t||a.t>now||now-a.t>cfg.max_age
    result.reason='Несогласованное время или устаревшие данные'; return;
end
assert(strcmp(a.quantity,b.quantity)&&strcmp(a.frame_id,b.frame_id)&& ...
    isequal(a.units,b.units)&&ischar(a.source)&&ischar(b.source)&& ...
    ~isempty(a.source)&&~isempty(b.source)&&~strcmp(a.source,b.source), ...
    'libekf:Comparison','Compare distinct sources of the same quantity, units and frame.');
% Единый порядок сохраняет знак ряда при перестановке аргументов.
[~,order]=sort({a.source,b.source});
if order(1)==2, tmp=a; a=b; b=tmp; cross_cov=cross_cov.'; end
[inputs_a,complete_a]=source_inputs(a);
[inputs_b,complete_b]=source_inputs(b);
result.shared_sources=intersect(inputs_a,inputs_b);
result.involved_sources=union(inputs_a,inputs_b);
result.provenance_complete=complete_a&&complete_b; result.t=a.t; result.interval=[a.t a.t];
if ~isempty(result.shared_sources)
    result.provenance='SHARED_INPUTS';
elseif complete_a&&complete_b
    result.provenance='DISJOINT_DECLARED_INPUTS';
end
assert(ischar(a.quantity)&&isrow(a.quantity)&&~isempty(a.quantity)&& ...
    ischar(a.frame_id)&&isrow(a.frame_id)&&~isempty(a.frame_id)&& ...
    iscell(a.units)&&~isempty(a.units)&& ...
    all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),a.units)), ...
    'libekf:Comparison','Invalid quantity, axes or units.');
bounded=isempty(cross_cov);
if isfield(cfg,'unknown_correlation_bound')
    assert(isequal(cfg.unknown_correlation_bound,true)||isequal(cfg.unknown_correlation_bound,false), ...
        'libekf:Comparison','Unknown-correlation bound must be explicitly enabled or disabled.');
end
if bounded
    if ~isfield(cfg,'unknown_correlation_bound')||~cfg.unknown_correlation_bound
        result.reason='Не задана взаимная ковариация ошибок'; return;
    end
    result.covariance_model='UPPER_BOUND'; result.nis_kind='LOWER_BOUND';
else
    result.covariance_model='JOINT'; result.nis_kind='EXACT';
end
% Длины полей исключают совпадение разных имён, содержащих разделители.
tokens=[{cfg.context;result.covariance_model;a.source;b.source;a.quantity;a.frame_id;num2str(numel(a.units))}; ...
    a.units(:);{num2str(complete_a);num2str(numel(inputs_a))};inputs_a(:); ...
    {num2str(complete_b);num2str(numel(inputs_b))};inputs_b(:)];
context='';
for k=1:numel(tokens), context=[context sprintf('%d:',numel(tokens{k})) tokens{k}]; end
x=a.value(:); y=b.value(:); n=numel(x);
assert(n>0&&numel(y)==n&&isreal(x)&&isreal(y)&&all(isfinite([x;y]))&& ...
    isequal(size(a.covariance),[n n])&&isequal(size(b.covariance),[n n])&& ...
    (bounded||isequal(size(cross_cov),[n n])),'libekf:Comparison','Invalid estimate dimensions.');
if bounded
    joint=blkdiag(a.covariance,b.covariance); % Только проверка двух маргинальных ковариаций.
else
    joint=[a.covariance cross_cov;cross_cov.' b.covariance];
end
assert(isreal(joint)&&all(isfinite(joint(:)))&&norm(joint-joint.','fro')<=1e-12*max(norm(joint,'fro'),realmin), ...
    'libekf:Covariance','Invalid comparison covariance.');
d=sqrt(diag(joint)); assert(isreal(d)&&all(d>0&isfinite(d)),'libekf:Covariance','Invalid variances.');
scaled=(joint./d)./d.'; values=eig((scaled+scaled.')/2);
assert(min(values)>=-100*eps*2*n*max(1,max(abs(values))), ...
    'libekf:Covariance','Inconsistent shared source errors.');
if bounded
    S=2*(a.covariance+b.covariance);
else
    S=a.covariance+b.covariance-cross_cov-cross_cov.';
end
S=(S+S.')/2;
if any(diag(S)<=0), result.reason='Вырожденная погрешность сравнения'; return; end
scale=sqrt(diag(S)); normalized=(S./scale)./scale.';
[L,flag]=chol((normalized+normalized.')/2,'lower');
if flag~=0, result.reason='Вырожденная погрешность сравнения'; return; end
value=L\((x-y)./scale); result.nis=value.'*value;
assert(all(isfinite(value))&&isfinite(result.nis),'libekf:Comparison','Comparison overflow.');
sample=struct('t',a.t,'value',value,'context',context,'alarm_only',bounded);
result.status='CONSISTENT'; result.reason='Расхождение в заданных пределах; общие ошибки не исключены';
if bounded
    result.status='NOT_READY'; result.reason='Порог по верхней оценке погрешности не превышен; согласие не подтверждено';
end
if result.nis>cfg.limit
    result.status='DISAGREEMENT'; result.reason='Источники расходятся; причина не определена';
end
end

function [inputs,complete]=source_inputs(source)
inputs={}; complete=false;
if isfield(source,'input_sources')
    inputs=source.input_sources;
    assert(iscell(inputs)&&all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),inputs)), ...
        'libekf:Comparison','Invalid physical input source identities.');
    inputs=unique(inputs(:));
end
if isfield(source,'provenance_complete')
    complete=source.provenance_complete;
    assert(islogical(complete)&&isscalar(complete)&&(~complete||~isempty(inputs)), ...
        'libekf:Comparison','Complete provenance requires physical input sources.');
end
end

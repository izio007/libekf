function [sample,result]=nav_compare_sources(a,b,cross_cov,now,cfg)
% Сравнение двух оценок одной физической величины; не определяет виновный источник.
sample=[]; result=struct('status','NOT_READY','reason','','nis',NaN);
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
x=a.value(:); y=b.value(:); n=numel(x);
assert(n>0&&numel(y)==n&&all(isfinite([x;y]))&& ...
    isequal(size(a.covariance),[n n])&&isequal(size(b.covariance),[n n])&& ...
    isequal(size(cross_cov),[n n]),'libekf:Comparison','Invalid estimate dimensions.');
joint=[a.covariance cross_cov;cross_cov.' b.covariance];
assert(all(isfinite(joint(:)))&&norm(joint-joint.','fro')<=1e-12*max(norm(joint,'fro'),realmin), ...
    'libekf:Covariance','Invalid comparison covariance.');
d=sqrt(diag(joint)); assert(isreal(d)&&all(d>0&isfinite(d)),'libekf:Covariance','Invalid variances.');
scaled=(joint./d)./d.'; values=eig((scaled+scaled.')/2);
assert(min(values)>=-100*eps*2*n*max(1,max(abs(values))), ...
    'libekf:Covariance','Inconsistent shared source errors.');
S=a.covariance+b.covariance-cross_cov-cross_cov.'; S=(S+S.')/2;
if any(diag(S)<=0), result.reason='Вырожденная погрешность сравнения'; return; end
scale=sqrt(diag(S)); normalized=(S./scale)./scale.';
[L,flag]=chol((normalized+normalized.')/2,'lower');
if flag~=0, result.reason='Вырожденная погрешность сравнения'; return; end
value=L\((x-y)./scale); result.nis=value.'*value;
assert(all(isfinite(value))&&isfinite(result.nis),'libekf:Comparison','Comparison overflow.');
sample=struct('t',a.t,'value',value,'context',cfg.context);
result.status='CONSISTENT'; result.reason='Расхождение в заданных пределах; общие ошибки не исключены';
if result.nis>cfg.limit
    result.status='DISAGREEMENT'; result.reason='Источники расходятся; причина не определена';
end
end

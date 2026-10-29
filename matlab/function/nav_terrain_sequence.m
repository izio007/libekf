function result=nav_terrain_sequence(paths,series,map,limits,search)
% Сравнение заданных возможных траекторий со всей последовательностью радиовысот.
% paths — ячейки: строки вариантов, столбцы моментов; в ячейке состояние p,q,t,frame_id.
% series.covariance — полная погрешность последовательности, включая карту и временные связи.
assert(iscell(paths)&&~isempty(paths)&&size(paths,2)==numel(series.t), ...
    'libekf:TerrainSequence','Expected one state per candidate and measurement time.');
t=series.t(:); z=series.z(:); n=numel(t);
assert(n>=2&&numel(z)==n&&all(isfinite([t;z]))&&all(diff(t)>0)&&all(z>=0) && ...
    strcmp(series.unit,'m')&&strcmp(series.quality,'VALID')&& ...
    ischar(series.source)&&~isempty(series.source)&&iscell(series.lineage)&& ...
    numel(series.lineage)==n&&all(cellfun(@(x)ischar(x)&&~isempty(x),series.lineage))&& ...
    numel(unique(series.lineage))==n,'libekf:TerrainSequence','Invalid sequence or repeated source records.');
assert(isequal(map.coverage_complete,true)&&strcmp(map.quality,'VALID')&& ...
    isequal(map.nearest_return_valid,true),'libekf:TerrainCoverage','Unverified terrain coverage or echo model.');
assert(isscalar(map.max_valid_range)&&isfinite(map.max_valid_range)&&map.max_valid_range>0&& ...
    all(z<=map.max_valid_range),'libekf:RadarRange','Sequence exceeds radar range.');
assert(numel(map.lever)==3&&all(isfinite(map.lever))&& ...
    isscalar(limits.absolute)&&isfinite(limits.absolute)&&limits.absolute>0&& ...
    isscalar(limits.difference)&&isfinite(limits.difference)&&limits.difference>=0, ...
    'libekf:TerrainSequence','Invalid lever or explicit comparison limits.');
R=series.covariance;
assert(isequal(size(R),[n n])&&all(isfinite(R(:)))&& ...
    norm(R-R.','fro')<=1e-12*max(norm(R,'fro'),realmin)&&all(diag(R)>0), ...
    'libekf:Covariance','Invalid sequence covariance.');
d=sqrt(diag(R)); A=(R./d)./d.'; [L,flag]=chol((A+A.')/2,'lower');
assert(flag==0,'libekf:Covariance','Sequence covariance must be positive definite.');
count=size(paths,1); offsets=[]; spacing=[];
if nargin>=5
    % Поиск постоянной ошибки положения счисленной траектории в местных осях.
    assert(size(paths,1)==1&&isstruct(search)&&isscalar(search)&& ...
        all(isfield(search,{'bounds','step','max_candidates'})), ...
        'libekf:TerrainSearch','Expected one nominal trajectory and search bounds.');
    bounds=search.bounds; step=search.step(:); capacity=search.max_candidates;
    assert(isnumeric(bounds)&&isreal(bounds)&&isequal(size(bounds),[3 2])&& ...
        all(isfinite(bounds(:)))&&all(bounds(:,2)>=bounds(:,1))&& ...
        isnumeric(step)&&isreal(step)&&numel(step)==3&&all(isfinite(step))&&all(step>0)&& ...
        isnumeric(capacity)&&isreal(capacity)&&isscalar(capacity)&&isfinite(capacity)&& ...
        capacity>=1&&capacity==floor(capacity), ...
        'libekf:TerrainSearch','Invalid metre bounds, spacing or candidate capacity.');
    width=bounds(:,2)-bounds(:,1); cells=ceil(width./step);
    count=prod(cells+1);
    assert(isfinite(count)&&count<=capacity,'libekf:TerrainCapacity', ...
        'Search exceeds candidate capacity; bounds and spacing are not silently reduced.');
    axes=cell(1,3); spacing=zeros(3,1);
    for axis=1:3
        axes{axis}=linspace(bounds(axis,1),bounds(axis,2),cells(axis)+1);
        assert(all(diff(axes{axis})>0),'libekf:TerrainSearch','Spacing is not representable at this origin.');
        if cells(axis)>0, spacing(axis)=width(axis)/cells(axis); end
    end
    [x,y,zoffset]=ndgrid(axes{:}); offsets=[x(:) y(:) zoffset(:)];
end
predicted=NaN(count,n); score=Inf(count,1);
for j=1:count
    for k=1:n
        if nargin>=5
            s=paths{1,k}; s.p=s.p(:)+offsets(j,:).';
        else
            s=paths{j,k};
        end
        assert(s.t==t(k)&&strcmp(s.frame_id,map.frame_id), ...
            'libekf:Epoch','Candidate must match measurement time and map frame.');
        C=nav_rotation(s.q);
        hit=nav_terrain_footprint(s.p(:)+C*map.lever(:),C*map.beam_body(:), ...
            map.half_angle,map.vertices,map.faces,0);
        % Для поиска нужен скаляр дальности; производная на ребре не требуется.
        predicted(j,k)=hit.range;
    end
    if all(isfinite(predicted(j,:)))&&all(predicted(j,:)<=map.max_valid_range)
        u=L\((z-predicted(j,:).')./d); score(j)=u.'*u;
    end
end
best=min(score); survivors=find(isfinite(score)&score<=limits.absolute&score<=best+limits.difference);
status='NO_MATCH';
if numel(survivors)>1, status='AMBIGUOUS'; elseif numel(survivors)==1, status='ONE_CANDIDATE'; end
result=struct('status',status,'candidates',survivors,'score',score,'predicted',predicted, ...
    't',t(end),'lineage',{series.lineage});
if nargin>=5
    result.offsets=offsets;
    result.search=struct('bounds',bounds,'step',spacing,'continuous_verified',false);
end
% ONE_CANDIDATE относится только к предоставленному набору, не ко всей области поиска.
end

function result=nav_visual_ambiguity(first,second,K,training,limits)
% Проверка альтернативных объяснений на отложенных соответствиях двух кадров.
% Входы — координаты после устранения искажений; результат не подтверждает навигацию.
assert(isnumeric(first)&&isnumeric(second)&&isequal(size(first),size(second))&& ...
    size(first,2)==2&&all(isfinite(first(:)))&&all(isfinite(second(:))), ...
    'libekf:Vision','Expected matching finite image coordinates.');
assert(isequal(size(K),[3 3])&&all(isfinite(K(:)))&& ...
    isequal(K(3,:),[0 0 1])&&K(2,1)==0&&K(1,1)>0&&K(2,2)>0, ...
    'libekf:Vision','Invalid pinhole calibration.');
assert(islogical(training)&&numel(training)==size(first,1)&& ...
    nnz(training)>=4&&nnz(~training)>=4,'libekf:Vision','Disjoint fitting and checking subsets required.');
assert(isscalar(limits.rotation_rms)&&isfinite(limits.rotation_rms)&&limits.rotation_rms>0&& ...
    isscalar(limits.plane_rms)&&isfinite(limits.plane_rms)&&limits.plane_rms>0, ...
    'libekf:Vision','Explicit positive checking limits required.');
support_check=isfield(limits,'min_support')||isfield(limits,'point_limit');
if support_check
    assert(all(isfield(limits,{'min_support','point_limit'}))&& ...
        isscalar(limits.min_support)&&isfinite(limits.min_support)&& ...
        limits.min_support>0&&limits.min_support<=1&& ...
        isscalar(limits.point_limit)&&isfinite(limits.point_limit)&&limits.point_limit>0, ...
        'libekf:Vision','Explicit point error and support fraction required.');
end
training=training(:); checking=~training;
first_rays=K\[first.';ones(1,size(first,1))]; second_rays=K\[second.';ones(1,size(second,1))];
p=first_rays(1:2,:).'; q=second_rays(1:2,:).';
first_rays=first_rays./vecnorm(first_rays); second_rays=second_rays./vecnorm(second_rays);
result=struct('status','UNASSESSED','rotation_rms',Inf,'plane_rms',Inf, ...
    'rotation_explains',false,'plane_explains',false,'geometry_ready',false, ...
    'rotation_support',NaN,'plane_support',NaN,'checking_points',nnz(checking));
indices=find(training); samples={indices};
if isfield(limits,'fit_samples')
    a=limits.fit_samples;
    assert(support_check&&isnumeric(a)&&isreal(a)&&size(a,2)==4&& ...
        ~isempty(a)&&all(isfinite(a(:)))&&all(a(:)==floor(a(:)))&& ...
        all(a(:)>=1&a(:)<=numel(indices))&& ...
        all(all(diff(sort(a,2),1,2)>0)), ...
        'libekf:Vision','Four distinct fitting indices per candidate required.');
    for j=1:size(a,1), samples{end+1}=indices(a(j,:)); end
end
ready=false(1,2); names={'rotation','plane'};
for kind=1:2
    best=[]; best_count=-1; best_cost=Inf;
    for j=1:numel(samples)
        M=fit_model(kind,samples{j},first_rays,second_rays,p,q);
        if isempty(M), continue; end
        errors=projection_errors(M,first_rays(:,training),second(training,:),K,kind);
        if support_check
            count=nnz(errors<=limits.point_limit);
            cost=sum(min(errors,limits.point_limit).^2);
        else
            count=0; cost=sum(errors.^2);
        end
        if isempty(best)||count>best_count||(count==best_count&&cost<best_cost)
            best=M; best_count=count; best_cost=cost;
        end
    end
    if isempty(best), continue; end
    % Одно уточнение по согласованным точкам построения, без отложенной части.
    if isfield(limits,'fit_samples')
        errors=projection_errors(best,first_rays(:,training),second(training,:),K,kind);
        agreed=indices(errors<=limits.point_limit);
        if numel(agreed)>=4
            refined=fit_model(kind,agreed,first_rays,second_rays,p,q);
            if ~isempty(refined)
                errors=projection_errors(refined,first_rays(:,training),second(training,:),K,kind);
                count=nnz(errors<=limits.point_limit);
                cost=sum(min(errors,limits.point_limit).^2);
                if count>best_count||(count==best_count&&cost<best_cost), best=refined; end
            end
        end
    end
    ready(kind)=true;
    errors=projection_errors(best,first_rays(:,checking),second(checking,:),K,kind);
    name=names{kind}; value=sqrt(mean(errors.^2));
    result.([name '_rms'])=value;
    explains=value<=limits.([name '_rms']);
    if support_check
        fraction=mean(errors<=limits.point_limit);
        result.([name '_support'])=fraction;
        explains=explains||fraction>=limits.min_support;
    end
    result.([name '_explains'])=explains;
end
result.geometry_ready=all(ready);
if result.rotation_explains||result.plane_explains, result.status='AMBIGUOUS'; end
end

function M=fit_model(kind,indices,a,b,p,q)
M=[];
if kind==1
    spread=svd(a(:,indices));
    if numel(spread)<3||spread(3)<=1e-10*spread(1), return; end
    [U,~,V]=svd(a(:,indices)*b(:,indices).');
    M=V*diag([1 1 det(V*U.')])*U.'; return;
end
A=zeros(2*numel(indices),8); rhs=zeros(2*numel(indices),1);
for j=1:numel(indices)
    i=indices(j); x=p(i,1); y=p(i,2); u=q(i,1); v=q(i,2);
    A(2*j-1,:)=[x y 1 0 0 0 -u*x -u*y]; rhs(2*j-1)=u;
    A(2*j,:)=[0 0 0 x y 1 -v*x -v*y]; rhs(2*j)=v;
end
singular=svd(A);
if numel(singular)<8||singular(8)<=1e-10*singular(1), return; end
h=A\rhs; M=[h(1:3).';h(4:6).';h(7:8).' 1];
end

function errors=projection_errors(M,rays,observed,K,kind)
projected=M*rays;
if kind==1, valid=projected(3,:)>0;
else, valid=abs(projected(3,:))>1e-12; end
errors=Inf(size(observed,1),1);
pixels=K*projected(:,valid); pixels=(pixels(1:2,:)./pixels(3,:)).';
errors(valid)=sqrt(sum((pixels-observed(valid,:)).^2,2));
end

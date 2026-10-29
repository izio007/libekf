function result=nav_visual_motion(matches,camera,limits,pixel_covariance)
% Кандидат относительного движения по исправленным пиксельным соответствиям.
% Масштаб не определяется; условная локальная ковариация необязательна. Поправка фильтра не выдаётся.
assert(isstruct(matches)&&isscalar(matches)&&all(isfield(matches, ...
    {'first','second','training','source','lineage','t_start','t_end','t_received','quality','units'})), ...
    'libekf:Vision','Incomplete image correspondences.');
assert(strcmp(matches.units,'pixel')&&strcmp(matches.quality,'UNASSESSED')&& ...
    isnumeric(matches.first)&&isreal(matches.first)&& ...
    isequal(size(matches.first),size(matches.second))&&size(matches.first,2)==2&& ...
    all(isfinite(matches.first(:)))&&all(isfinite(matches.second(:)))&& ...
    islogical(matches.training)&&numel(matches.training)==size(matches.first,1)&& ...
    nnz(matches.training)>=8&&nnz(~matches.training)>=4, ...
    'libekf:Vision','Disjoint fitting and checking pixel coordinates required.');
assert(ischar(matches.source)&&isrow(matches.source)&&~isempty(matches.source)&& ...
    iscell(matches.lineage)&&numel(matches.lineage)==2&& ...
    all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),matches.lineage))&& ...
    ~strcmp(matches.lineage{1},matches.lineage{2}), ...
    'libekf:Lineage','Two distinct frame identifiers required.');
times=[matches.t_start matches.t_end matches.t_received];
assert(isnumeric(times)&&isreal(times)&&numel(times)==3&&all(isfinite(times))&& ...
    times(1)<times(2)&&times(2)<=times(3), 'libekf:Time','Invalid frame pair times.');
assert(isa(camera,'cameraIntrinsics')||isa(camera,'cameraParameters'), ...
    'libekf:Vision','Pinhole camera calibration required.');
required={'trials','confidence','fit_distance','check_distance','min_check_support', ...
    'min_front_fraction','seed','ambiguity'};
assert(isstruct(limits)&&isscalar(limits)&&all(isfield(limits,required)), ...
    'libekf:Vision','Explicit motion estimation and checking settings required.');
v=[limits.trials limits.confidence limits.fit_distance limits.check_distance ...
    limits.min_check_support limits.min_front_fraction limits.seed];
assert(isnumeric(v)&&isreal(v)&&numel(v)==7&&all(isfinite(v))&& ...
    v(1)>=1&&v(1)==floor(v(1))&&v(2)>0&&v(2)<100&&v(3)>0&&v(4)>0&& ...
    v(5)>0&&v(5)<=1&&v(6)>0&&v(6)<=1&&v(7)>=0&&v(7)<=2^32-1&&v(7)==floor(v(7)), ...
    'libekf:Vision','Invalid motion estimation settings.');
if isfield(limits,'selection')
    covariance_args={}; if nargin>=4, covariance_args={pixel_covariance}; end
    result=check_selection(matches,camera,limits,covariance_args); return;
end
K=[camera.FocalLength(1) camera.Skew camera.PrincipalPoint(1); ...
    0 camera.FocalLength(2) camera.PrincipalPoint(2);0 0 1];
intrinsics=cameraIntrinsics(camera.FocalLength,camera.PrincipalPoint,camera.ImageSize,'Skew',camera.Skew);
result=struct('source',matches.source,'lineage',{matches.lineage}, ...
    't_start',matches.t_start,'t_end',matches.t_end,'t_received',matches.t_received, ...
    'quality','INSUFFICIENT','scale_known',false,'geometry_ready',false, ...
    'rotation',[],'direction',[],'checking_errors',[],'checking_support',0,'front_fraction',0, ...
    'fitting_indices',[]);
% Альтернативы используют то же исходное разделение, без отбора по основной модели.
result.ambiguity=nav_visual_ambiguity(matches.first,matches.second,K,matches.training,limits.ambiguity);
if strcmp(result.ambiguity.status,'AMBIGUOUS'), result.quality='AMBIGUOUS'; end
fit=find(matches.training); check=find(~matches.training);
saved=rng; cleanup=onCleanup(@()rng(saved)); rng(limits.seed,'twister'); %#ok<NASGU>
[E,inliers,status]=estimateEssentialMatrix(matches.first(fit,:),matches.second(fit,:),intrinsics, ...
    'MaxNumTrials',limits.trials,'Confidence',limits.confidence,'MaxDistance',limits.fit_distance);
if status~=0||nnz(inliers)<8, return; end
chosen=fit(inliers); result.fitting_indices=chosen;
[pose,fraction]=estrelpose(E,intrinsics,matches.first(chosen,:),matches.second(chosen,:));
if numel(pose)~=1||~isfinite(fraction)||norm(pose.Translation)==0, return; end
result.rotation=pose.R;
result.direction=pose.Translation(:)/norm(pose.Translation);
result.front_fraction=fraction;
F=K.'\E/K;
x=[matches.first(check,:).';ones(1,numel(check))];
y=[matches.second(check,:).';ones(1,numel(check))];
a=F*x; b=F.'*y; denominator=sum(a(1:2,:).^2,1)+sum(b(1:2,:).^2,1);
errors=Inf(size(denominator)); good=denominator>0;
errors(good)=abs(sum(y(:,good).*a(:,good),1))./sqrt(denominator(good));
result.checking_errors=errors(:);
result.checking_support=mean(errors<=limits.check_distance);
result.quality='REJECTED';
if strcmp(result.ambiguity.status,'AMBIGUOUS')
    result.quality='AMBIGUOUS'; return;
end
result.geometry_ready=result.ambiguity.geometry_ready&&fraction>=limits.min_front_fraction&& ...
    result.checking_support>=limits.min_check_support;
if result.geometry_ready, result.quality='UNASSESSED'; end
if nargin>=4
    assert(isfield(limits,'rank_tolerance')&&isscalar(limits.rank_tolerance)&& ...
        isreal(limits.rank_tolerance)&&isfinite(limits.rank_tolerance)&& ...
        limits.rank_tolerance>0&&limits.rank_tolerance<1, ...
        'libekf:Vision','Explicit local rank tolerance required.');
    count=size(matches.first,1);
    assert(isnumeric(pixel_covariance)&&isreal(pixel_covariance)&& ...
        isequal(size(pixel_covariance),[4*count 4*count])&&all(isfinite(pixel_covariance(:)))&& ...
        norm(pixel_covariance-pixel_covariance.','fro')<=1e-12*max(norm(pixel_covariance,'fro'),realmin)&& ...
        all(diag(pixel_covariance)>0), 'libekf:Covariance','Expected full pixel covariance.');
    d=sqrt(diag(pixel_covariance)); D=spdiags(1./d,0,4*count,4*count);
    [~,flag]=chol(D*pixel_covariance*D);
    assert(flag==0,'libekf:Covariance','Pixel covariance must be positive definite.');
    indices=reshape((4*(chosen(:)-1)+(1:4)).',[],1);
    C=pixel_covariance(indices,indices);
    if isfield(limits,'refinement')
        [R,d,refined]=refine_motion(matches.first(chosen,:),matches.second(chosen,:), ...
            K,result.rotation,result.direction,C,limits.rank_tolerance,limits.refinement);
        result.refinement=refined;
        if ~refined.converged
            result.geometry_ready=false; result.quality='UNASSESSED'; return;
        end
        result.rotation=R; result.direction=d;
        E=R.'*nav_skew(d); F=K.'\E/K;
        a=F*x; b=F.'*y; denominator=sum(a(1:2,:).^2,1)+sum(b(1:2,:).^2,1);
        errors=Inf(size(denominator)); good=denominator>0;
        errors(good)=abs(sum(y(:,good).*a(:,good),1))./sqrt(denominator(good));
        result.checking_errors=errors(:); result.checking_support=mean(errors<=limits.check_distance);
        first_rays=K\[matches.first(chosen,:).';ones(1,numel(chosen))];
        second_rays=K\[matches.second(chosen,:).';ones(1,numel(chosen))];
        front=false(numel(chosen),1);
        for k=1:numel(chosen)
            depths=[first_rays(:,k) -R*second_rays(:,k)]\d;
            front(k)=all(isfinite(depths)&depths>0);
        end
        result.front_fraction=mean(front);
        result.geometry_ready=result.ambiguity.geometry_ready&& ...
            result.front_fraction>=limits.min_front_fraction&&result.checking_support>=limits.min_check_support;
        result.quality='REJECTED'; if result.geometry_ready, result.quality='UNASSESSED'; end
        result.local_error=refined.local_error;
    else
        result.local_error=local_error(matches.first(chosen,:),matches.second(chosen,:), ...
            K,result.rotation,result.direction,C,limits.rank_tolerance);
    end
end
end

function [out,work]=local_error(first,second,K,R,d,C,tolerance)
% Линеаризованная взвешенная модель, условная по выбранным соответствиям.
% Порядок пикселей: x1,y1,x2,y2 для каждой пары; параметры: две компоненты направления, три поворота.
[~,axis]=min(abs(d)); unit=zeros(3,1); unit(axis)=1;
v=cross(d,unit); v=v/norm(v); B=[v cross(d,v)];
E=R.'*nav_skew(d); n=size(first,1);
a=K\[first.';ones(1,n)]; b=K\[second.';ones(1,n)];
pixel_map=K\[eye(2);zeros(1,2)];
J=zeros(n,5); G=sparse(n,4*n); axes=eye(3);
for k=1:n
    for j=1:2, J(k,j)=b(:,k).'*R.'*nav_skew(B(:,j))*a(:,k); end
    for j=1:3, J(k,j+2)=-b(:,k).'*nav_skew(axes(:,j))*E*a(:,k); end
    G(k,4*k-3:4*k)=[b(:,k).'*E*pixel_map a(:,k).'*E.'*pixel_map];
end
V=full(G*C*G.'); V=(V+V.')/2;
scale=sqrt(diag(V)); normalized=(V./scale)./scale.';
[L,flag]=chol(normalized,'lower');
assert(flag==0&&all(isfinite(scale)&scale>0),'libekf:Covariance','Singular epipolar noise.');
A=L\(J./scale); [~,S,W]=svd(A,'econ'); singular=diag(S);
rank=sum(singular>max(tolerance,100*eps*max(size(A)))*max(singular));
work=struct('L',L,'scale',scale,'A',A,'residual',sum(b.*(E*a),1).');
out=struct('rank',rank,'tangent_basis',B,'state_jacobian',J,'pixel_jacobian',G, ...
    'singular_values',singular,'covariance',[],'quality','UNOBSERVABLE');
if rank==5
    factor=W./singular.'; covariance=factor*factor.';
    out.covariance=(covariance+covariance.')/2;
    out.quality='LOCAL_LINEAR';
end
end

function [R,d,info]=refine_motion(first,second,K,R,d,C,rank_tolerance,options)
% Обобщённый метод наименьших квадратов с ошибками координат обоих кадров.
assert(isstruct(options)&&isscalar(options)&&all(isfield(options, ...
    {'max_iterations','angle_tolerance','pixel_tolerance'})), ...
    'libekf:Vision','Explicit refinement convergence settings required.');
v=[options.max_iterations options.angle_tolerance options.pixel_tolerance];
assert(isnumeric(v)&&isreal(v)&&numel(v)==3&&all(isfinite(v))&& ...
    v(1)>=1&&v(1)==floor(v(1))&&all(v(2:3)>0), ...
    'libekf:Vision','Invalid refinement convergence settings.');
y=reshape([first second].',[],1); correction=zeros(size(y));
info=struct('converged',false,'iterations',0,'cost',Inf);
for iteration=1:options.max_iterations
    adjusted=reshape(y+correction,4,[]).';
    [local,work]=local_error(adjusted(:,1:2),adjusted(:,3:4),K,R,d,C,rank_tolerance);
    info.iterations=iteration;
    if local.rank<5, return; end
    rhs=work.L\((work.residual-local.pixel_jacobian*correction)./work.scale);
    step=-(work.A\rhs);
    multiplier=(work.L.'\(rhs+work.A*step))./work.scale;
    candidate=-C*local.pixel_jacobian.'*multiplier;
    if any(~isfinite(step))||any(~isfinite(candidate)), return; end
    d=d+local.tangent_basis*step(1:2); d=d/norm(d);
    R=R*expm(nav_skew(step(3:5)));
    done=norm(step)<=options.angle_tolerance&&max(abs(candidate-correction))<=options.pixel_tolerance;
    correction=candidate;
    if done
        adjusted=reshape(y+correction,4,[]).';
        [local,work]=local_error(adjusted(:,1:2),adjusted(:,3:4),K,R,d,C,rank_tolerance);
        pixel_error=abs(work.residual)./sqrt(full(sum(local.pixel_jacobian.^2,2)));
        if local.rank<5||any(~isfinite(pixel_error))||max(pixel_error)>options.pixel_tolerance, return; end
        info.converged=true; info.local_error=local;
        info.cost=correction.'*(C\correction); return;
    end
end
end

function result=check_selection(matches,camera,limits,covariance_args)
selection=limits.selection;
assert(isstruct(selection)&&isscalar(selection)&&all(isfield(selection, ...
    {'subsets','rotation_limit','direction_limit'}))&&iscell(selection.subsets)&& ...
    ~isempty(selection.subsets), 'libekf:Vision','Explicit fitting subsets and angular limits required.');
angles=[selection.rotation_limit selection.direction_limit];
assert(isnumeric(angles)&&isreal(angles)&&numel(angles)==2&&all(isfinite(angles))&& ...
    all(angles>0&angles<pi),'libekf:Vision','Invalid candidate separation limits.');
base=rmfield(limits,'selection'); fitting=find(matches.training); checking=find(~matches.training);
for k=1:numel(selection.subsets)
    subset=selection.subsets{k};
    assert(isnumeric(subset)&&isreal(subset)&&isvector(subset)&&numel(subset)>=8&& ...
        all(isfinite(subset))&&numel(unique(subset))==numel(subset)&&all(ismember(subset,fitting)), ...
        'libekf:Vision','Candidates must use distinct fitting indices, never held-out points.');
end
result=nav_visual_motion(matches,camera,base,covariance_args{:});
candidates={candidate_summary(result)};
for k=1:numel(selection.subsets)
    subset=selection.subsets{k}(:); indices=[subset;checking];
    part=matches; part.first=matches.first(indices,:); part.second=matches.second(indices,:);
    part.training=[true(numel(subset),1);false(numel(checking),1)];
    args={};
    if ~isempty(covariance_args)
        pixel_indices=reshape((4*(indices-1)+(1:4)).',[],1);
        C=covariance_args{1}; args={C(pixel_indices,pixel_indices)};
    end
    candidate=nav_visual_motion(part,camera,base,args{:});
    candidate.fitting_indices=indices(candidate.fitting_indices);
    candidates{end+1}=candidate_summary(candidate);
end
ready=find(cellfun(@(x)x.geometry_ready,candidates));
rotation_spread=0; direction_spread=0;
for i=1:numel(ready)
    for j=i+1:numel(ready)
        a=candidates{ready(i)}; b=candidates{ready(j)};
        rotation_spread=max(rotation_spread,acos(max(-1,min(1,(trace(a.rotation.'*b.rotation)-1)/2))));
        direction_spread=max(direction_spread,acos(max(-1,min(1,a.direction.'*b.direction))));
    end
end
status='INSUFFICIENT';
if numel(ready)==numel(candidates), status='CONSISTENT'; end
if rotation_spread>selection.rotation_limit||direction_spread>selection.direction_limit|| ...
    any(cellfun(@(x)strcmp(x.quality,'AMBIGUOUS'),candidates))
    status='AMBIGUOUS'; result.quality='AMBIGUOUS';
end
if ~strcmp(status,'CONSISTENT')
    result.geometry_ready=false;
    result.rotation=[]; result.direction=[]; result.fitting_indices=[];
    if isfield(result,'refinement'), result=rmfield(result,'refinement'); end
    if isfield(result,'local_error'), result=rmfield(result,'local_error'); end
end
result.selection=struct('status',status,'rotation_spread',rotation_spread, ...
    'direction_spread',direction_spread,'candidates',{candidates});
end

function summary=candidate_summary(candidate)
summary=struct('rotation',candidate.rotation,'direction',candidate.direction, ...
    'quality',candidate.quality,'geometry_ready',candidate.geometry_ready, ...
    'checking_support',candidate.checking_support,'fitting_indices',candidate.fitting_indices);
if isfield(candidate,'local_error')
    summary.covariance=candidate.local_error.covariance;
    summary.tangent_basis=candidate.local_error.tangent_basis;
end
end

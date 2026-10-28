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
training=training(:); checking=~training;
first_rays=K\[first.';ones(1,size(first,1))]; second_rays=K\[second.';ones(1,size(second,1))];
p=first_rays(1:2,:).'; q=second_rays(1:2,:).';
first_rays=first_rays./vecnorm(first_rays); second_rays=second_rays./vecnorm(second_rays);
result=struct('status','UNASSESSED','rotation_rms',Inf,'plane_rms',Inf, ...
    'rotation_explains',false,'plane_explains',false,'geometry_ready',false);
% Плохое пространственное распределение не считается отсутствием альтернативы.
spread=svd(first_rays(:,training));
if numel(spread)<3||spread(3)<=1e-10*spread(1), return; end
[U,~,V]=svd(first_rays(:,training)*second_rays(:,training).');
rotation=U*diag([1 1 det(U*V.')])*V.';
predicted=rotation.'*first_rays(:,checking);
if all(predicted(3,:)>0)
    pixels=K*predicted; pixels=(pixels(1:2,:)./pixels(3,:)).';
    result.rotation_rms=sqrt(mean(sum((pixels-second(checking,:)).^2,2)));
end
indices=find(training); A=zeros(2*numel(indices),8); rhs=zeros(2*numel(indices),1);
for j=1:numel(indices)
    i=indices(j); x=p(i,1); y=p(i,2); u=q(i,1); v=q(i,2);
    A(2*j-1,:)=[x y 1 0 0 0 -u*x -u*y]; rhs(2*j-1)=u;
    A(2*j,:)=[0 0 0 x y 1 -v*x -v*y]; rhs(2*j)=v;
end
singular=svd(A);
if numel(singular)<8||singular(8)<=1e-10*singular(1), return; end
h=A\rhs; homography=[h(1:3).';h(4:6).';h(7:8).' 1];
projected=homography*[p(checking,:).';ones(1,nnz(checking))];
if all(abs(projected(3,:))>1e-12)
    pixels=K*projected; pixels=(pixels(1:2,:)./pixels(3,:)).';
    result.plane_rms=sqrt(mean(sum((pixels-second(checking,:)).^2,2)));
end
result.geometry_ready=true;
result.rotation_explains=result.rotation_rms<=limits.rotation_rms;
result.plane_explains=result.plane_rms<=limits.plane_rms;
if result.rotation_explains||result.plane_explains, result.status='AMBIGUOUS'; end
end

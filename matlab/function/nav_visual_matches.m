function result=nav_visual_matches(first,second,calibration,options)
% Соответствия двух кадров видимого диапазона до оценки движения.
% Требуется Computer Vision Toolbox. Функция не формирует поправку фильтра.
for frame={first,second}
    f=frame{1};
    assert(isstruct(f)&&isscalar(f)&&all(isfield(f, ...
        {'image','source','id','t','t_received','quality'})), ...
        'libekf:Vision','Incomplete image frame.');
    assert(isa(f.image,'uint8')&&isreal(f.image)&& ...
        (ismatrix(f.image)||(ndims(f.image)==3&&size(f.image,3)==3))&& ...
        size(f.image,1)>=64&&size(f.image,2)>=64, ...
        'libekf:Vision','Expected uint8 grayscale or RGB image.');
    assert(ischar(f.source)&&isrow(f.source)&&~isempty(f.source)&& ...
        ischar(f.id)&&isrow(f.id)&&~isempty(f.id)&&strcmp(f.quality,'VALID'), ...
        'libekf:Quality','Invalid frame identity or quality.');
    assert(isnumeric(f.t)&&isscalar(f.t)&&isreal(f.t)&&isfinite(f.t)&& ...
        isnumeric(f.t_received)&&isscalar(f.t_received)&&isreal(f.t_received)&& ...
        isfinite(f.t_received)&&f.t_received>=f.t, ...
        'libekf:Time','Invalid frame acquisition or arrival time.');
end
assert(strcmp(first.source,second.source)&&~strcmp(first.id,second.id)&& ...
    first.t<second.t&&isequal(size(first.image),size(second.image)), ...
    'libekf:Vision','Distinct ordered frames of one camera required.');
assert(isstruct(calibration)&&isscalar(calibration)&& ...
    all(isfield(calibration,{'intrinsics','source','valid'}))&& ...
    isequal(calibration.valid,true)&&strcmp(calibration.source,first.source), ...
    'libekf:Quality','Valid calibration of this camera required.');
camera=calibration.intrinsics;
assert((isa(camera,'cameraIntrinsics')||isa(camera,'cameraParameters'))&& ...
    isequal(camera.ImageSize,[size(first.image,1) size(first.image,2)]), ...
    'libekf:Vision','Calibration image size differs from input.');
assert(isstruct(options)&&isscalar(options)&& ...
    all(isfield(options,{'max_points','match_ratio','fit_fraction','seed'})), ...
    'libekf:Vision','Explicit matching and partition settings required.');
assert(isscalar(options.max_points)&&isfinite(options.max_points)&& ...
    options.max_points>=12&&options.max_points==floor(options.max_points)&& ...
    isscalar(options.match_ratio)&&isfinite(options.match_ratio)&& ...
    options.match_ratio>0&&options.match_ratio<1&& ...
    isscalar(options.fit_fraction)&&isfinite(options.fit_fraction)&& ...
    options.fit_fraction>0&&options.fit_fraction<1&& ...
    isscalar(options.seed)&&isfinite(options.seed)&&options.seed>=0&& ...
    options.seed<=2^32-1&&options.seed==floor(options.seed), ...
    'libekf:Vision','Invalid matching or partition settings.');
capture_known=isfield(first,'capture')||isfield(second,'capture');
centres=[first.t second.t]; capture=struct('specified',false);
if capture_known
    assert(isfield(first,'capture')&&isfield(second,'capture'),'libekf:Time','Capture metadata required for both frames.');
    [centres(1),span1]=capture_times(first); [centres(2),span2]=capture_times(second);
    assert(centres(1)<centres(2),'libekf:Time','Image exposure centres are not ordered.');
    capture=struct('specified',true,'exposure_duration',[first.capture.duration second.capture.duration], ...
        'readout_time',[first.capture.readout_time second.capture.readout_time], ...
        'row_direction',[first.capture.row_direction second.capture.row_direction], ...
        'intervals',[span1;span2]);
end
result=struct('source',first.source,'lineage',{{first.id,second.id}}, ...
    't_start',centres(1),'t_end',centres(2),'t_received',max(first.t_received,second.t_received), ...
    'units','pixel','first',zeros(0,2),'second',zeros(0,2),'training',false(0,1), ...
    'quality','INSUFFICIENT','scale_known',false,'capture',capture,'point_times',zeros(0,2));
images={im2gray(first.image),im2gray(second.image)};
features=cell(1,2); points=cell(1,2);
for k=1:2
    found=detectORBFeatures(images{k});
    if found.Count==0, return; end
    found=selectStrongest(found,options.max_points);
    [features{k},points{k}]=extractFeatures(images{k},found);
end
pairs=matchFeatures(features{1},features{2},'Unique',true,'MaxRatio',options.match_ratio);
n=size(pairs,1); nfit=floor(n*options.fit_fraction);
if nfit<8||n-nfit<4, return; end
selected_first=double(points{1}(pairs(:,1)).Location);
selected_second=double(points{2}(pairs(:,2)).Location);
if capture_known
    rows=[selected_first(:,2) selected_second(:,2)];
    height=size(first.image,1);
    result.point_times=centres+(rows-(height+1)/2)/(height-1).*capture.readout_time.*capture.row_direction;
end
result.first=double(undistortPoints(selected_first,camera));
result.second=double(undistortPoints(selected_second,camera));
assert(all(isfinite(result.first(:)))&&all(isfinite(result.second(:))), ...
    'libekf:Vision','Nonfinite undistorted coordinates.');
stream=RandStream('mt19937ar','Seed',options.seed);
order=randperm(stream,n); result.training=false(n,1); result.training(order(1:nfit))=true;
result.quality='UNASSESSED';
if isfield(options,'localization_covariance')||isfield(calibration,'parameter_covariance')
    assert(all(isfield(options,{'localization_covariance','localization_independent'}))&& ...
        all(isfield(calibration,{'parameter_covariance','parameter_noise_independent'}))&& ...
        isequal(options.localization_independent,true)&&isequal(calibration.parameter_noise_independent,true), ...
        'libekf:Correlation','Explicit independent localization and calibration noise models required.');
    raw=reshape([selected_first selected_second].',2,[]).';
    [result.pixel_covariance,result.calibration_jacobian,pixel_J]=pixel_uncertainty(raw,camera, ...
        options.localization_covariance,calibration.parameter_covariance);
    result.covariance_model='LOCAL_FIXED_MATCHES';
end
if isfield(options,'capture_covariance')
    assert(capture_known&&isfield(result,'pixel_covariance')&& ...
        isfield(options,'capture_noise_independent')&&isequal(options.capture_noise_independent,true), ...
        'libekf:Correlation','Capture noise requires explicit independence from localization and calibration.');
    result.coordinate_time_covariance=coordinate_time_noise(raw,pixel_J,result.calibration_jacobian, ...
        camera.ImageSize,first.capture,second.capture,options.localization_covariance, ...
        calibration.parameter_covariance,options.capture_covariance);
end
% Отложенные точки не участвовали ни в оценке движения, ни в геометрическом отборе.
% Погрешность движения, различимость и масштаб по числу совпадений не выводятся.
end

function [C,J,pixel_J]=pixel_uncertainty(raw,camera,Cpoint,Ccal)
% Параметры: fx,fy,cx,cy,skew,k1,k2,k3,p1,p2; исправленные пиксели при фиксированном номинальном K.
assert(isnumeric(Cpoint)&&isreal(Cpoint)&&isequal(size(Cpoint),[2 2])&& ...
    all(isfinite(Cpoint(:)))&&norm(Cpoint-Cpoint.','fro')<=1e-12*max(norm(Cpoint,'fro'),realmin), ...
    'libekf:Covariance','Expected localization covariance in pixel squared.');
[~,flag]=chol(Cpoint); assert(flag==0,'libekf:Covariance','Localization covariance must be positive definite.');
assert(isnumeric(Ccal)&&isreal(Ccal)&&isequal(size(Ccal),[10 10])&&all(isfinite(Ccal(:)))&& ...
    norm(Ccal-Ccal.','fro')<=1e-12*max(norm(Ccal,'fro'),realmin)&&all(diag(Ccal)>=0), ...
    'libekf:Covariance','Expected joint calibration covariance.');
scale=sqrt(diag(Ccal)); zero=scale==0;
assert(~any(any(Ccal(zero,:)~=0)),'libekf:Covariance','Zero calibration variance has nonzero dependence.');
normalized=(Ccal(~zero,~zero)./scale(~zero))./scale(~zero).';
assert(isempty(normalized)||min(eig((normalized+normalized.')/2))>=-1e-12, ...
    'libekf:Covariance','Calibration covariance is indefinite.');
radial=zeros(1,3); radial(1:numel(camera.RadialDistortion))=camera.RadialDistortion;
parameters=[camera.FocalLength camera.PrincipalPoint camera.Skew radial camera.TangentialDistortion];
K=[parameters(1) parameters(5) parameters(3);0 parameters(2) parameters(4);0 0 1];
n=size(raw,1); J=zeros(2*n,10); pixel_J=sparse(2*n,2*n);
for j=1:10
    step=eps^(1/3)*max(1,abs(parameters(j))); plus=parameters; minus=parameters;
    plus(j)=plus(j)+step; minus(j)=minus(j)-step;
    derivative=(effective_pixels(raw,plus,K,camera.ImageSize)-effective_pixels(raw,minus,K,camera.ImageSize))/(2*step);
    J(:,j)=reshape(derivative.',[],1);
end
for j=1:2
    step=eps^(1/3)*max(1,max(abs(raw(:,j)))); plus=raw; minus=raw;
    plus(:,j)=plus(:,j)+step; minus(:,j)=minus(:,j)-step;
    derivative=(effective_pixels(plus,parameters,K,camera.ImageSize)-effective_pixels(minus,parameters,K,camera.ImageSize))/(2*step);
    for k=1:n, pixel_J(2*k-1:2*k,2*k-2+j)=derivative(k,:).'; end
end
C=pixel_J*kron(speye(n),Cpoint)*pixel_J.';
if any(Ccal(:)), C=C+J*Ccal*J.'; end
C=(C+C.')/2;
assert(all(isfinite(C(:)))&&all(isfinite(J(:))),'libekf:Covariance','Nonfinite propagated image covariance.');
end

function out=effective_pixels(raw,p,Knom,image_size)
camera=cameraIntrinsics(p(1:2),p(3:4),image_size,'Skew',p(5), ...
    'RadialDistortion',p(6:8),'TangentialDistortion',p(9:10));
corrected=double(undistortPoints(raw,camera));
K=[p(1) p(5) p(3);0 p(2) p(4);0 0 1];
rays=K\[corrected.';ones(1,size(raw,1))]; pixels=Knom*rays;
out=(pixels(1:2,:)./pixels(3,:)).';
end

function [centre,interval]=capture_times(frame)
c=frame.capture;
assert(isstruct(c)&&isscalar(c)&&all(isfield(c, ...
    {'duration','readout_time','stamp_reference','row_direction','unit','valid'}))&& ...
    isequal(c.valid,true)&&strcmp(c.unit,'s')&& ...
    isnumeric(c.duration)&&isscalar(c.duration)&&isreal(c.duration)&&isfinite(c.duration)&&c.duration>=0&& ...
    isnumeric(c.readout_time)&&isscalar(c.readout_time)&&isreal(c.readout_time)&&isfinite(c.readout_time)&&c.readout_time>=0&& ...
    isscalar(c.row_direction)&&any(c.row_direction==[-1 1]), ...
    'libekf:Time','Valid exposure duration, readout time and row direction required.');
half=(c.duration+c.readout_time)/2;
switch c.stamp_reference
    case 'first_row_start', centre=frame.t+half;
    case 'frame_center', centre=frame.t;
    case 'last_row_end', centre=frame.t-half;
    otherwise, error('libekf:Time','Unknown capture timestamp reference.');
end
interval=centre+[-half half];
assert(all(isfinite(interval))&&frame.t_received>=interval(2), ...
    'libekf:Time','Frame received before acquisition completed.');
end

function C=coordinate_time_noise(raw,pixel_J,Jcal,image_size,first,second,Cpoint,Ccal,Ctime)
% Совместные параметры времени: stamp1, exposure1, readout1, stamp2, exposure2, readout2.
assert(isnumeric(Ctime)&&isreal(Ctime)&&isequal(size(Ctime),[6 6])&&all(isfinite(Ctime(:)))&& ...
    norm(Ctime-Ctime.','fro')<=1e-12*max(norm(Ctime,'fro'),realmin)&&all(diag(Ctime)>=0), ...
    'libekf:Covariance','Expected covariance of six capture parameters in seconds squared.');
scale=sqrt(diag(Ctime)); zero=scale==0;
assert(~any(any(Ctime(zero,:)~=0)),'libekf:Covariance','Impossible capture covariance.');
normalized=(Ctime(~zero,~zero)./scale(~zero))./scale(~zero).';
assert(isempty(normalized)||min(eig((normalized+normalized.')/2))>=-1e-12, ...
    'libekf:Covariance','Indefinite capture covariance.');
n=size(raw,1); coordinate_rows=reshape([3*(1:n)-2;3*(1:n)-1],[],1);
A=sparse(3*n,2*n); A(coordinate_rows,:)=pixel_J;
B=zeros(3*n,10); B(coordinate_rows,:)=Jcal;
T=zeros(3*n,6); frames={first,second};
for k=1:n
    frame=mod(k-1,2)+1; c=frames{frame};
    alpha=0;
    if strcmp(c.stamp_reference,'first_row_start'), alpha=.5; end
    if strcmp(c.stamp_reference,'last_row_end'), alpha=-.5; end
    row=c.row_direction*((raw(k,2)-1)/(image_size(1)-1)-.5);
    A(3*k,2*k)=c.row_direction*c.readout_time/(image_size(1)-1);
    T(3*k,3*frame-2:3*frame)=[1 alpha alpha+row];
end
C=A*kron(speye(n),Cpoint)*A.';
if any(Ccal(:)), C=C+B*Ccal*B.'; end
if any(Ctime(:)), C=C+T*Ctime*T.'; end
C=(C+C.')/2;
assert(all(isfinite(C(:))),'libekf:Covariance','Nonfinite coordinate/time covariance.');
end

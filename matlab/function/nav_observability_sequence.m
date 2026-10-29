function out=nav_observability_sequence(history,scales,targets,tolerance)
% Определимость текущих величин по временной последовательности без априорной P.
% x(k+1)=T{k}*x(k)+w(k); y(k)=H{k}*x(k)+n(k).
assert(isstruct(history)&&isscalar(history)&&all(isfield(history,{'H','transitions','times', ...
    'process_covariance','measurement_covariance','process_measurement_covariance','frame_id'})), ...
    'libekf:Observability','Incomplete temporal observation model.');
H=history.H; transitions=history.transitions; times=history.times(:);
assert(iscell(H)&&isvector(H)&&~isempty(H)&&iscell(transitions)&&numel(transitions)==numel(H)-1&& ...
    isnumeric(times)&&isreal(times)&&numel(times)==numel(H)&&all(isfinite(times))&&all(diff(times)>0)&& ...
    ischar(history.frame_id)&&isrow(history.frame_id)&&~isempty(history.frame_id), ...
    'libekf:Observability','Invalid history epochs, transitions or frame.');
n=size(H{1},2); count=numel(H); rows=cellfun(@(h)size(h,1),H); m=sum(rows); nw=n*(count-1);
assert(n>0&&all(rows>0)&&size(targets,2)==n,'libekf:Observability','Invalid state or observation dimensions.');
Q=history.process_covariance; R=history.measurement_covariance; U=history.process_measurement_covariance;
assert(isequal(size(Q),[nw nw])&&isequal(size(R),[m m])&&isequal(size(U),[nw m]), ...
    'libekf:Covariance','Invalid full temporal noise dimensions.');
joint=[Q U;U.' R];
assert(isreal(joint)&&all(isfinite(joint(:)))&& ...
    norm(joint-joint.','fro')<=1e-12*max(norm(joint,'fro'),realmin)&&all(diag(joint)>=0), ...
    'libekf:Covariance','Invalid process/measurement covariance.');
d=sqrt(diag(joint)); zero=d==0;
assert(~any(any(joint(zero,:)~=0)),'libekf:Covariance','Impossible zero-variance dependence.');
d(zero)=1; normalized=(joint./d)./d.'; eigenvalues=eig((normalized+normalized.')/2);
assert(min(eigenvalues)>=-100*eps*(nw+m)*max(1,max(abs(eigenvalues))), ...
    'libekf:Covariance','Indefinite process/measurement covariance.');
Phi=eye(n); B=zeros(n,nw); O=zeros(m,n); F=zeros(m,nw); first=1;
for k=1:count
    hk=H{k};
    assert(isnumeric(hk)&&isreal(hk)&&size(hk,2)==n&&all(isfinite(hk(:))), ...
        'libekf:Observability','Invalid observation Jacobian in history.');
    if k>1
        T=transitions{k-1};
        assert(isnumeric(T)&&isreal(T)&&isequal(size(T),[n n])&&all(isfinite(T(:))), ...
            'libekf:Observability','Invalid error transition.');
        Phi=T*Phi; B=T*B; B(:,(k-2)*n+(1:n))=eye(n);
    end
    ix=first:first+rows(k)-1; first=first+rows(k);
    O(ix,:)=hk*Phi; F(ix,:)=hk*B;
end
effective=F*Q*F.'+R+F*U+U.'*F.'; effective=(effective+effective.')/2;
target_process=targets*B;
noise=struct('covariance',target_process*Q*target_process.', ...
    'cross_cov',target_process*(Q*F.'+U));
out=nav_observability(O,effective,scales,targets*Phi,tolerance,noise);
out.t=times(end); out.frame_id=history.frame_id;
end

function [out,Pout,J,reference_out,motion_out]=nav_change_frame(s,P,from,to,reference,motion)
% Exact Earth-fixed change of origin/axes; no flat-Earth approximation.
% Temporal meaning and right/body-local attitude errors remain unchanged.
reference_out=[]; motion_out=[];
if nargin>=6, motion_out=change_motion(motion,numel(reference.position),from,to); end
if iscell(s)
    assert(numel(s)==2&&isequal(size(P),[42 42])&&all(isfinite(P(:))), ...
        'libekf:Covariance','Expected two states and full joint covariance.');
    [first,~,Jfirst]=nav_change_frame(s{1},P(1:21,1:21),from,to);
    [second,~,Jsecond]=nav_change_frame(s{2},P(22:42,22:42),from,to);
    J=blkdiag(Jfirst,Jsecond); out={first,second}; Pout=J*P*J.';
    if nargin>=5, reference_out=change_reference(reference,J,from,to); end
    return;
end
assert(strcmp(s.frame_id,from.frame_id),'libekf:Frame','Wrong source state frame.');
nav_quaternion(from.C_e_l); nav_quaternion(to.C_e_l);
T=to.C_e_l.'*from.C_e_l;
out=s;
out.p=to.C_e_l.'*(from.origin_ecef(:)-to.origin_ecef(:))+T*s.p(:);
out.v=T*s.v(:); out.wind=T*s.wind(:);
out.q=nav_quaternion(T*nav_rotation(s.q)); out.frame_id=to.frame_id;
assert(isequal(size(P),[21 21]) && all(isfinite(P(:))), 'libekf:Covariance','Expected finite 21x21 covariance.');
J=eye(21); J(1:3,1:3)=T; J(4:6,4:6)=T; J(16:18,16:18)=T;
Pout=J*P*J.';
if nargin>=5, reference_out=change_reference(reference,J,from,to); end
end

function out=change_reference(reference,J,from,to)
assert(isstruct(reference)&&isscalar(reference)&& ...
    all(isfield(reference,{'position','covariance','cross_cov','frame_id'}))&& ...
    strcmp(reference.frame_id,from.frame_id),'libekf:Frame','Invalid reference-set source frame.');
d=numel(reference.position); n=size(J,1);
assert(d>=3&&mod(d,3)==0&&isequal(size(reference.position),[d 1])&& ...
    isreal(reference.position)&&all(isfinite(reference.position))&& ...
    isequal(size(reference.covariance),[d d])&&isreal(reference.covariance)&& ...
    all(isfinite(reference.covariance(:)))&& ...
    isequal(size(reference.cross_cov),[n d])&&isreal(reference.cross_cov)&& ...
    all(isfinite(reference.cross_cov(:))), 'libekf:Covariance','Invalid reference-set dimensions.');
T=to.C_e_l.'*from.C_e_l;
offset=to.C_e_l.'*(from.origin_ecef(:)-to.origin_ecef(:));
K=kron(eye(d/3),T); out=reference;
out.position=K*reference.position+repmat(offset,d/3,1);
out.covariance=K*reference.covariance*K.';
out.cross_cov=J*reference.cross_cov*K.'; out.frame_id=to.frame_id;
assert(all(isfinite(out.position))&&all(isfinite(out.covariance(:)))&&all(isfinite(out.cross_cov(:))), ...
    'libekf:Covariance','Reference frame conversion overflow.');
end

function out=change_motion(model,d,from,to)
assert(isstruct(model)&&isscalar(model)&&isfield(model,'frame_id')&& ...
    strcmp(model.frame_id,from.frame_id)&&d>=3&&mod(d,3)==0, ...
    'libekf:Frame','Invalid reference-motion source frame.');
nav_quaternion(from.C_e_l); nav_quaternion(to.C_e_l);
if isfield(model,'F')
    fields={'F','drift','noise_density'};
    assert(~any(isfield(model,{'A','shift','Q'})),'libekf:Peer','Mixed reference motion models.');
else
    fields={'A','shift','Q'};
end
assert(all(isfield(model,fields)),'libekf:Peer','Incomplete reference motion model.');
M=model.(fields{1}); v=model.(fields{2}); N=model.(fields{3});
assert(isequal(size(M),[d d])&&isreal(M)&&all(isfinite(M(:)))&& ...
    isequal(size(v),[d 1])&&isreal(v)&&all(isfinite(v))&& ...
    isequal(size(N),[d d])&&isreal(N)&&all(isfinite(N(:))), ...
    'libekf:Peer','Invalid reference motion dimensions.');
K=kron(eye(d/3),to.C_e_l.'*from.C_e_l);
b=repmat(to.C_e_l.'*(from.origin_ecef(:)-to.origin_ecef(:)),d/3,1);
out=model; out.frame_id=to.frame_id;
out.(fields{1})=K*M*K.'; out.(fields{3})=K*N*K.';
out.(fields{2})=K*v-out.(fields{1})*b;
if ~isfield(model,'F'), out.shift=out.shift+b; end
assert(all(isfinite(out.(fields{1})(:)))&&all(isfinite(out.(fields{2})(:)))&& ...
    all(isfinite(out.(fields{3})(:))),'libekf:Covariance','Reference motion transform overflow.');
end

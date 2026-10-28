function [out,Pout,J]=nav_change_frame(s,P,from,to)
% Exact Earth-fixed change of origin/axes; no flat-Earth approximation.
% Temporal meaning and right/body-local attitude errors remain unchanged.
if iscell(s)
    assert(numel(s)==2&&isequal(size(P),[42 42])&&all(isfinite(P(:))), ...
        'libekf:Covariance','Expected two states and full joint covariance.');
    [first,~,Jfirst]=nav_change_frame(s{1},P(1:21,1:21),from,to);
    [second,~,Jsecond]=nav_change_frame(s{2},P(22:42,22:42),from,to);
    J=blkdiag(Jfirst,Jsecond); out={first,second}; Pout=J*P*J.'; return;
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
end

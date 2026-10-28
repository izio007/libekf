function q=nav_quaternion(C)
% Rotation matrix to Hamilton scalar-first quaternion; stable near pi.
assert(isequal(size(C),[3 3]) && all(isfinite(C(:))) && ...
    norm(C.'*C-eye(3),'fro')<1e-10 && abs(det(C)-1)<1e-10, ...
    'libekf:Rotation','Expected proper orthonormal rotation.');
d=[1+trace(C);1+2*C(1,1)-trace(C);1+2*C(2,2)-trace(C);1+2*C(3,3)-trace(C)];
[~,i]=max(d); q=zeros(4,1); q(i)=.5*sqrt(max(d(i),0)); divisor=4*q(i);
switch i
    case 1, q(2:4)=[C(3,2)-C(2,3);C(1,3)-C(3,1);C(2,1)-C(1,2)]/divisor;
    case 2, q([1 3 4])=[C(3,2)-C(2,3);C(1,2)+C(2,1);C(1,3)+C(3,1)]/divisor;
    case 3, q([1 2 4])=[C(1,3)-C(3,1);C(1,2)+C(2,1);C(2,3)+C(3,2)]/divisor;
    case 4, q(1:3)=[C(2,1)-C(1,2);C(1,3)+C(3,1);C(2,3)+C(3,2)]/divisor;
end
q=q/norm(q); if q(1)<0, q=-q; end
end

function [h, up, lat, lon] = nav_ellipsoid_height(p, a, e2)
% Geodetic height and ECEF height gradient for an oblate ellipsoid.
% Intended for near-Earth exterior navigation, not the ellipsoid interior.
p=p(:);
assert(numel(p)==3 && all(isfinite(p)) && isscalar(a) && a>0 && ...
    isscalar(e2) && e2>=0 && e2<1, 'libekf:Ellipsoid', 'Invalid ellipsoid or position.');
assert(norm(p)>a/2, 'libekf:EllipsoidDomain', 'Deep ellipsoid interior is outside this navigation profile.');
r=hypot(p(1),p(2));
if r==0
    lat=sign(p(3))*pi/2; lon=NaN;
    up=[0;0;sign(p(3))]; h=abs(p(3))-a*sqrt(1-e2); return;
end
lat=atan2(p(3),r*(1-e2)); converged=false;
for k=1:30
    N=a/sqrt(1-e2*sin(lat)^2);
    next=atan2(p(3)+e2*N*sin(lat),r);
    if abs(next-lat)<1e-14, converged=true; lat=next; break; end
    lat=next;
end
assert(converged,'libekf:EllipsoidConvergence','Geodetic conversion failed.');
up=[cos(lat)*p(1)/r; cos(lat)*p(2)/r; sin(lat)];
lon=atan2(p(2),p(1));
% Projection avoids division by cos(latitude) near a pole.
h=p.'*up-a*sqrt(1-e2*sin(lat)^2);
end

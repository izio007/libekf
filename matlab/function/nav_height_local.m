function [h,up,lat,lon]=nav_height_local(p,env)
% Stable ellipsoidal height relative to a fixed Earth anchor.
% Avoid subtracting two Earth-radius-sized support functions for each sample.
B=env.C_e_l; origin=env.origin_ecef(:); delta=B*p(:);
[~,up,lat,lon]=nav_ellipsoid_height(origin+delta,env.ellipsoid_a,env.ellipsoid_e2);
[h0,n0]=nav_ellipsoid_height(origin,env.ellipsoid_a,env.ellipsoid_e2);
D=diag([1 1 sqrt(1-env.ellipsoid_e2)]);
u=D*up; u0=D*n0; s=norm(u); s0=norm(u0);
denominator=s0*(u.'*u0+s*s0);
if denominator>1e-12
    % Rationalized Gram determinant; use up-n0 before the cross product.
    w=cross(D*(up-n0),u0);
    surfaceGap=-env.ellipsoid_a*(w.'*w)/denominator;
else
    % At nearly antipodal normals the difference is large, not cancelling.
    surfaceGap=env.ellipsoid_a*(u.'*u0/s0-s);
end
h=up.'*delta+h0*(up.'*n0)+surfaceGap;
end

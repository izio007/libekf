function out=nav_geodetic_output(s,P,env)
% Geodetic output only; covariance is NOT a guarantee of integrity.
% lat/lon in degrees, alt is ellipsoidal metres. No implicit MSL conversion.
assert(isequal(size(P),[21 21]) && all(isfinite(P(:))) && ...
    norm(P-P.','fro')<=1e-12*max(norm(P,'fro'),realmin), 'libekf:Covariance','Invalid error covariance.');
assert(strcmp(s.frame_id,env.frame_id),'libekf:Frame','Wrong common navigation anchor.');
nav_quaternion(env.C_e_l); % Proper rotation, not a tangent plane approximation.
[alt,up,lat,lon]=nav_height_local(s.p,env);
N=env.ellipsoid_a/sqrt(1-env.ellipsoid_e2*sin(lat)^2);
M=env.ellipsoid_a*(1-env.ellipsoid_e2)/(1-env.ellipsoid_e2*sin(lat)^2)^1.5;
J=zeros(3,21); J(3,1:3)=up.'*env.C_e_l;
out=struct('lat',lat*180/pi,'lon',lon*180/pi,'alt',alt,'t',s.t, ...
    'units',{{'deg';'deg';'m'}},'altitude_reference','ELLIPSOID', ...
    'quality',struct('latitude','UNASSESSED','longitude','UNASSESSED','altitude','UNASSESSED'));
if isfinite(lon)
    north=[-sin(lat)*cos(lon);-sin(lat)*sin(lon);cos(lat)];
    east=[-sin(lon);cos(lon);0];
    J(1,1:3)=180/pi*north.'*env.C_e_l/(M+alt);
    J(2,1:3)=180/pi*east.'*env.C_e_l/((N+alt)*cos(lat));
else
    % Longitude and its first-order chart are undefined exactly at a pole.
    J(1:2,:)=NaN; out.quality.latitude='CHART_SINGULAR'; out.quality.longitude='UNDEFINED';
end
out.covariance=J*P*J.'; out.jacobian=J;
end

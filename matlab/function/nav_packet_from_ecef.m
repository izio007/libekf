function out=nav_packet_from_ecef(m,env)
% Explicit adapter: global originals -> local geometry and measurement covariance.
% Original packet is not modified. Basis C_e_l maps local vectors to ECEF.
B=env.C_e_l; origin=env.origin_ecef(:);
assert(isequal(size(B),[3 3]) && norm(B.'*B-eye(3),'fro')<1e-10 && ...
    abs(det(B)-1)<1e-10 && numel(origin)==3 && all(isfinite(origin)), ...
    'libekf:Frame','Invalid local Earth anchor.');
assert(strcmp(m.frame,'ECEF'),'libekf:Frame','Adapter expects ECEF originals.');
out=m; out.frame='LOCAL'; out.frame_id=env.frame_id; T=[];
switch m.model
    case {'gnss_position','external_position'}
        out.z=B.'*(m.z(:)-origin); T=B.';
    case {'gnss_velocity','external_velocity'}
        out.z=B.'*m.z(:); T=B.';
    case {'gnss_position_velocity','external_position_velocity'}
        out.z=[B.'*(m.z(1:3)-origin);B.'*m.z(4:6)]; T=blkdiag(B.',B.');
    case 'attitude'
        out.z=nav_quaternion(B.'*nav_rotation(m.z));
        % Right/body-local angular noise covariance is unchanged.
end
if ~isempty(T), out.R=T*m.R*T.'; end
if isfield(m,'state_noise_cov')&&~isempty(T)
    assert(size(m.state_noise_cov,2)==size(T,2),'libekf:Correlation','Wrong cross covariance columns.');
    out.state_noise_cov=m.state_noise_cov*T.';
end
for name={'station','reference_station','path_origin','plane_point','field_origin'}
    if isfield(m,name{1}), out.(name{1})=B.'*(m.(name{1})(:)-origin); end
end
if isfield(m,'station_velocity'), out.station_velocity=B.'*m.station_velocity(:); end
if isfield(m,'surface_velocity'), out.surface_velocity=B.'*m.surface_velocity(:); end
if isfield(m,'basis'), out.basis=B.'*m.basis; end
if isfield(m,'projection'), out.projection=m.projection*B; end
if isfield(m,'plane_normal'), out.plane_normal=B.'*m.plane_normal(:); end
if isfield(m,'vertices'), out.vertices=(B.'*(m.vertices.'-origin)).'; end
if isfield(m,'directions'), out.directions=(B.'*m.directions.').'; end
if isfield(m,'landmarks'), out.landmarks=(B.'*(m.landmarks.'-origin)).'; end
if isfield(m,'stations'), out.stations=(B.'*(m.stations.'-origin)).'; end
if isfield(m,'field_vector'), out.field_vector=B.'*m.field_vector(:); end
if isfield(m,'field_gradient'), out.field_gradient=B.'*m.field_gradient*B; end
if isfield(m,'peer_covariance'), out.peer_covariance=B.'*m.peer_covariance*B; end
if isfield(m,'peer_cross_cov'), out.peer_cross_cov=m.peer_cross_cov*B; end
% Lever, beam and omega_ib are body vectors and must not be rotated here.
end

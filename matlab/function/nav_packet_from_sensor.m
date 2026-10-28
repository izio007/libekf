function out=nav_packet_from_sensor(m, sensor, env)
% Convert sensor-local spatial components/geometry to the common local frame.
% Velocity fields mean Earth-relative physical vectors expressed in the frame,
% not coordinate derivatives in a moving frame. Calibrated body angles/scalars
% retain their native meaning. Sensor axes and Earth anchor are explicit.
assert(isstruct(m)&&isscalar(m)&&all(isfield(m,{'model','z','R','frame','frame_id'})), ...
    'libekf:Packet','Incomplete sensor packet.');
assert(strcmp(m.frame,'SENSOR_LOCAL') && strcmp(m.frame_id,sensor.id), ...
    'libekf:Frame','Wrong sensor-local frame ID.');
B=env.C_e_l; A=sensor.C_e_s;
assert(isequal(size(A),[3 3]) && norm(A.'*A-eye(3),'fro')<1e-10 && abs(det(A)-1)<1e-10, ...
    'libekf:Frame','Sensor axes must be a proper orthonormal frame.');
% This difference is formed once for the anchor, not for every raw measurement.
offset=B.'*(sensor.origin_ecef(:)-env.origin_ecef(:)); T=B.'*A;
% Reuse the geometry transform without constructing large absolute positions.
localenv=struct('C_e_l',T.','origin_ecef',-T.'*offset,'frame_id',env.frame_id);
native=m; native.frame='ECEF'; out=nav_packet_from_ecef(native,localenv);
end

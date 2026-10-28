function obs=nav_observe_source(s,m,sensor,env)
% Native sensor residual and covariance, common-state Jacobian.
% No sensor is declared trustworthy merely because its model was constructed.
local=nav_packet_from_sensor(m,sensor,env);
obs=nav_measurement(s,local,env);
T=env.C_e_l.'*sensor.C_e_s; offset=env.C_e_l.'*(sensor.origin_ecef(:)-env.origin_ecef(:));
switch m.model
    case {'gnss_position','external_position'}
        W=T.'; obs.predicted=W*(obs.predicted-offset);
    case {'gnss_velocity','external_velocity'}
        W=T.'; obs.predicted=W*obs.predicted;
    case {'gnss_position_velocity','external_position_velocity'}
        W=blkdiag(T.',T.'); obs.predicted=W*(obs.predicted-[offset;zeros(3,1)]);
    otherwise
        W=eye(numel(obs.residual));
end
obs.residual=W*obs.residual; obs.H=W*obs.H; obs.R=W*obs.R*W.';
if isfield(obs,'state_noise_cov'), obs.state_noise_cov=obs.state_noise_cov*W.'; end
obs.source_frame_id=sensor.id;
end

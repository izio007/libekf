function out=nav_geodetic_output(s,P,env,control)
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
% Проверки назначаются каждой выходной величине явно; согласие не означает VALID.
if nargin>=4
    assert(isstruct(control)&&isscalar(control)&& ...
        isfield(control,'state_time')&&isequal(control.state_time,s.t)&& ...
        isfield(control,'frame_id')&&strcmp(control.frame_id,s.frame_id)&& ...
        isfield(control,'now')&&isscalar(control.now)&&isfinite(control.now)&& ...
        control.now>=s.t&&isfield(control,'checks')&&isstruct(control.checks), ...
        'libekf:WatchOutput','Control must refer to this state and processing time.');
    out.checked_at=control.now; out.checks=control.checks;
    names={'latitude','longitude','altitude'};
    for k=1:numel(names)
        name=names{k}; status='NOT_READY';
        if isfield(control.checks,name)
            checks=control.checks.(name);
            assert(iscell(checks),'libekf:WatchOutput','Expected a cell of required checks.');
            complete=~isempty(checks); alarm=false;
            for j=1:numel(checks)
                c=checks{j};
                assert(isstruct(c)&&isscalar(c)&&all(isfield(c, ...
                    {'t','ready','alarm','instant_alarm','status'}))&& ...
                    isequal(c.t,control.now)&&islogical(c.ready)&&numel(c.ready)==3&& ...
                    islogical(c.alarm)&&numel(c.alarm)==3&& ...
                    islogical(c.instant_alarm)&&isscalar(c.instant_alarm), ...
                    'libekf:WatchOutput','Checks must be evaluated at publication time.');
                has_alarm=c.instant_alarm||any(c.alarm);
                expected='NOT_READY';
                if has_alarm, expected='ALARM';
                elseif all(c.ready), expected='CONSISTENT'; end
                assert(strcmp(c.status,expected),'libekf:WatchOutput','Contradictory check status.');
                alarm=alarm||has_alarm; complete=complete&&all(c.ready);
            end
            if alarm, status='ALARM'; elseif complete, status='CONSISTENT'; end
        end
        % Неопределённость долготы на полюсе нельзя снять согласованностью сторожа.
        if strcmp(out.quality.(name),'UNASSESSED'), out.quality.(name)=status; end
    end
end
end

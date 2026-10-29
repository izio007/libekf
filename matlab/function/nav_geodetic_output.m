function out=nav_geodetic_output(s,P,env,control)
% Geodetic output only; covariance is NOT a guarantee of integrity.
% lat/lon in degrees, alt is ellipsoidal metres. No implicit MSL conversion.
% Полный снимок автоматически передаёт происхождение и неустранённые отказы.
if isstruct(s)&&isscalar(s)&&isfield(s,'state')
    assert(nargin==4&&isempty(P)&&isfield(s,'P')&&isstruct(control)&&isscalar(control)&& ...
        all(isfield(control,{'now','checks'})),'libekf:WatchOutput', ...
        'Use nav_geodetic_output(filter,[],env,control) with publication time and checks.');
    metadata=struct('state_time',s.state.t,'frame_id',s.state.frame_id);
    for field={'source_history','sources_complete','excluded_sources','invalidated_sources'}
        if isfield(s,field{1}), metadata.(field{1})=s.(field{1}); end
    end
    for field={'state_time','frame_id','source_history','sources_complete','excluded_sources','invalidated_sources'}
        name=field{1};
        if isfield(control,name)
            assert(isfield(metadata,name)&&isequaln(control.(name),metadata.(name)), ...
                'libekf:WatchOutput','Control cannot override snapshot provenance or unresolved faults.');
        elseif isfield(metadata,name)
            control.(name)=metadata.(name);
        end
    end
    out=nav_geodetic_output(s.state,s.P,env,control);
    return;
end
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
out.age=[]; out.sources=nav_sources([]); out.sources_complete=false; out.invalidated_sources={};
% Проверки назначаются каждой выходной величине явно; согласие не означает VALID.
if nargin>=4
    assert(isstruct(control)&&isscalar(control)&& ...
        isfield(control,'state_time')&&isequal(control.state_time,s.t)&& ...
        isfield(control,'frame_id')&&strcmp(control.frame_id,s.frame_id)&& ...
        isfield(control,'now')&&isscalar(control.now)&&isfinite(control.now)&& ...
        control.now>=s.t&&isfield(control,'checks')&&isstruct(control.checks), ...
        'libekf:WatchOutput','Control must refer to this state and processing time.');
    out.checked_at=control.now; out.checks=control.checks;
    out.age=control.now-s.t;
    if isfield(control,'source_history')
        out.sources=nav_sources(control.source_history);
        assert(all([out.sources.last_epoch]<=s.t),'libekf:SourceHistory','Source contribution is newer than the state.');
        excluded={}; if isfield(control,'excluded_sources'), excluded=control.excluded_sources; end
        assert(iscell(excluded)&&all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),excluded)), ...
            'libekf:SourceHistory','Invalid excluded source identities.');
        for j=1:numel(out.sources)
            out.sources(j).age=control.now-out.sources(j).last_epoch;
            out.sources(j).excluded=any(strcmp(excluded,out.sources(j).id));
        end
    end
    if isfield(control,'sources_complete')
        assert(isfield(control,'source_history')&&islogical(control.sources_complete)&&isscalar(control.sources_complete), ...
            'libekf:SourceHistory','Source-history completeness must be explicit.');
        out.sources_complete=control.sources_complete;
    end
    if isfield(control,'invalidated_sources')
        ids=control.invalidated_sources;
        assert(iscell(ids)&&all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),ids)), ...
            'libekf:WatchOutput','Invalid unresolved source influence list.');
        out.invalidated_sources=unique(ids(:));
    end
    if ~isempty(out.sources)&&isfield(out.sources,'excluded')
        out.invalidated_sources=unique([out.invalidated_sources(:); ...
            reshape({out.sources([out.sources.excluded]).id},[],1)]);
    end
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
                alarm_only=false; if isfield(c,'alarm_only'), alarm_only=c.alarm_only; end
                assert(islogical(alarm_only)&&isscalar(alarm_only), ...
                    'libekf:WatchOutput','Invalid diagnostic purpose.');
                has_alarm=c.instant_alarm||any(c.alarm);
                expected='NOT_READY';
                if has_alarm, expected='ALARM';
                elseif all(c.ready)&&~alarm_only, expected='CONSISTENT'; end
                assert(strcmp(c.status,expected),'libekf:WatchOutput','Contradictory check status.');
                alarm=alarm||has_alarm; complete=complete&&all(c.ready)&&~alarm_only;
            end
            if alarm, status='ALARM'; elseif complete, status='CONSISTENT'; end
        end
        % Неопределённость долготы на полюсе нельзя снять согласованностью сторожа.
        if ~isempty(out.invalidated_sources), status='ALARM'; end
        if strcmp(out.quality.(name),'UNASSESSED'), out.quality.(name)=status; end
    end
    if isfield(control,'geometry')
        g=control.geometry;
        assert(isstruct(g)&&isscalar(g)&&all(isfield(g, ...
            {'scales','relative_tolerance','max_variance'}))&& ...
            numel(g.max_variance)==3&&isreal(g.max_variance)&& ...
            all(isfinite(g.max_variance(:))&g.max_variance(:)>0), ...
            'libekf:WatchOutput','Expected observation geometry and explicit output variance limits.');
        usable=all(isfinite(J),2);
        if isfield(g,'history')
            assert(~any(isfield(g,{'H','R'}))&&isstruct(g.history)&&isscalar(g.history)&& ...
                isfield(g.history,'frame_id'), ...
                'libekf:WatchOutput','History geometry must end at the published state and use its axes.');
            geometry=nav_observability_sequence(g.history,g.scales,J(usable,:),g.relative_tolerance);
            assert(geometry.t==s.t&&strcmp(geometry.frame_id,s.frame_id), ...
                'libekf:WatchOutput','History must end at the published epoch and frame.');
        else
            assert(all(isfield(g,{'H','R'})),'libekf:WatchOutput','Missing observation geometry.');
            geometry=nav_observability(g.H,g.R,g.scales,J(usable,:),g.relative_tolerance);
        end
        out.observable=false(3,1); out.measurement_variance=Inf(3,1);
        out.observable(usable)=geometry.identified;
        out.measurement_variance(usable)=geometry.variance;
        out.accuracy_sufficient=out.observable&out.measurement_variance<=g.max_variance(:);
        for k=1:3
            name=names{k};
            % Тревога и особенности координат имеют приоритет над геометрией источников.
            if ~usable(k)||strcmp(out.quality.(name),'ALARM'), continue; end
            if ~out.observable(k)
                out.quality.(name)='UNOBSERVABLE';
            elseif ~out.accuracy_sufficient(k)
                out.quality.(name)='INSUFFICIENT_ACCURACY';
            end
        end
    end
end
end

function obs = nav_measurement(s, m, env)
% Build one physical observation and analytic local-error Jacobian.
% All working geometry is LOCAL with a named Earth-fixed anchor.
% No state mutation, gating, fusion, or assumption of source independence here.
required={'source','model','z','R','units','frame','frame_id','t_measure','t_received', ...
    'quality','lineage','correlation_group','noise_model'};
for k=1:numel(required)
    assert(isfield(m,required{k}), 'libekf:Packet', ['Missing ' required{k}]);
end
assert(ischar(m.source) && ~isempty(m.source) && ischar(m.model), 'libekf:Packet','Invalid source/model.');
assert(strcmp(m.quality,'VALID'), 'libekf:Quality','Only VALID input enters observation construction.');
assert(strcmp(m.frame,'LOCAL') && strcmp(m.frame_id,s.frame_id) && strcmp(s.frame_id,env.frame_id), ...
    'libekf:Frame','Packet, state and environment must share one local frame ID.');
assert(isscalar(s.t) && isfinite(s.t) && isscalar(m.t_measure) && isfinite(m.t_measure) && ...
    isscalar(m.t_received) && isfinite(m.t_received) && ...
    isfield(env,'time_tolerance') && isscalar(env.time_tolerance) && ...
    isfinite(env.time_tolerance) && env.time_tolerance>=0, 'libekf:Time','Invalid times/tolerance.');
assert(abs(s.t-m.t_measure)<=env.time_tolerance, 'libekf:Epoch','State must be at measurement epoch.');
assert(m.t_received+env.time_tolerance>=m.t_measure,'libekf:Time','Receive time precedes measurement.');
assert(iscell(m.lineage) && ~isempty(m.lineage) && all(cellfun(@(x)ischar(x)&&~isempty(x),m.lineage)), ...
    'libekf:Lineage','Provide original information record IDs.');
assert(ischar(m.correlation_group) && ~isempty(m.correlation_group) && ...
    any(strcmp(m.noise_model,{'independent','joint_required'})), ...
    'libekf:Correlation','Declare correlation group and noise assumption.');
C=nav_rotation(s.q); Rtransform=[]; z=m.z(:);
fields={'p','v','ba','bg','wind'};
for k=1:numel(fields), vector3(s.(fields{k})); end
assert(all(isfinite([s.bh s.clock s.drift])) && ...
    all(cellfun(@(f)isscalar(s.(f)),{'bh','clock','drift'})), 'libekf:State','Invalid augmented state.');
H=zeros(0,21); angular=[]; added_variance=0; derived_cross=[];
lever=zeros(3,1);
if isfield(m,'lever'), lever=vector3(m.lever); end
ps=s.p(:)+C*lever;
Jp=zeros(3,21); Jp(:,1:3)=eye(3); Jp(:,7:9)=-C*nav_skew(lever);
switch m.model
    case {'gnss_position','external_position'}
        h=ps; H=Jp; units=repmat({'m'},3,1);
    case {'gnss_velocity','external_velocity','gnss_position_velocity','external_position_velocity','range_rate','doppler_velocity','doppler_frequency'}
        if norm(lever)>0
            assert(isfield(m,'omega_ib') && isfield(env,'Omega'), 'libekf:Lever','Nonzero lever needs body rates and Earth rate.');
            omega=vector3(m.omega_ib)-s.bg(:); Omega=vector3(env.Omega);
        else
            omega=zeros(3,1); Omega=zeros(3,1);
        end
        rot=cross(omega,lever);
        va=s.v(:)+C*rot-cross(Omega,C*lever);
        Jv=zeros(3,21); Jv(:,4:6)=eye(3);
        Jv(:,7:9)=-C*nav_skew(rot)+nav_skew(Omega)*C*nav_skew(lever);
        Jv(:,13:15)=C*nav_skew(lever);
        if any(strcmp(m.model,{'doppler_velocity','doppler_frequency'}))
            require_fields(m,{'beams_body','surface_valid','beam_geometry_valid','surface_velocity'});
            beams=m.beams_body;
            assert(size(beams,2)==3 && size(beams,1)>0 && all(isfinite(beams(:))) && ...
                all(abs(sqrt(sum(beams.^2,2))-1)<1e-10) && ...
                isequal(m.surface_valid,true) && isequal(m.beam_geometry_valid,true), ...
                'libekf:Doppler','Invalid Doppler directions or reflecting surface.');
            vb=C.'*(va-vector3(m.surface_velocity));
            h=beams*vb; H=beams*C.'*Jv;
            H(:,7:9)=H(:,7:9)+beams*nav_skew(vb);
            units=repmat({'m/s'},size(beams,1),1);
            if strcmp(m.model,'doppler_frequency')
                require_fields(m,{'wavelength','frequency_sign','frequency_correction', ...
                    'monostatic','calibration_valid'});
                wavelength=m.wavelength; sign=m.frequency_sign; correction=m.frequency_correction(:);
                assert(isnumeric(wavelength)&&isreal(wavelength)&&isscalar(wavelength)&& ...
                    isfinite(wavelength)&&wavelength>0&&isnumeric(sign)&&isreal(sign)&& ...
                    isscalar(sign)&&any(sign==[-1 1])&&isnumeric(correction)&&isreal(correction)&& ...
                    numel(correction)==size(beams,1)&&all(isfinite(correction))&& ...
                    isequal(m.monostatic,true)&&isequal(m.calibration_valid,true), ...
                    'libekf:Doppler','Invalid wavelength, frequency convention or calibration.');
                scale=2*sign/wavelength;
                h=scale*h+correction; H=scale*H;
                units=repmat({'Hz'},size(beams,1),1);
            end
        elseif strcmp(m.model,'range_rate')
            require_fields(m,{'station','station_velocity','rate_correction','geometry_at_epoch'});
            assert(isfield(m,'geometry_at_epoch') && isequal(m.geometry_at_epoch,true), ...
                'libekf:GeometryEpoch','Range-rate geometry must be prepared at observation epoch.');
            d=ps-vector3(m.station); rho=norm(d);
            assert(rho>0,'libekf:Geometry','Zero satellite range.'); u=d/rho;
            dv=va-vector3(m.station_velocity);
            h=u.'*dv+s.drift+m.rate_correction;
            H=((eye(3)-u*u.')*dv/rho).'*Jp+u.'*Jv; H(21)=H(21)+1;
            units={'m/s'};
        elseif any(strcmp(m.model,{'gnss_position_velocity','external_position_velocity'}))
            h=[ps;va]; H=[Jp;Jv]; units=[repmat({'m'},3,1);repmat({'m/s'},3,1)];
        else
            h=va; H=Jv; units=repmat({'m/s'},3,1);
        end
    case 'cooperative_range'
        peer_indices=1:3;
        if isfield(m,'reference_bank')
            b=m.reference_bank;
            assert(isstruct(b)&&isscalar(b)&&all(isfield(b, ...
                {'ids','position','covariance','cross_cov','frame_id','t','valid'}))&& ...
                isfield(m,'peer_id')&&ischar(m.peer_id), ...
                'libekf:Peer','Incomplete reference set or selected ID.');
            assert(iscell(b.ids)&&~isempty(b.ids)&& ...
                all(cellfun(@(x)ischar(x)&&~isempty(x),b.ids))&& ...
                numel(unique(b.ids))==numel(b.ids)&& ...
                isequal(size(b.position),[3*numel(b.ids) 1])&&isreal(b.position)&&all(isfinite(b.position))&& ...
                islogical(b.valid)&&numel(b.valid)==numel(b.ids)&&strcmp(b.frame_id,s.frame_id), ...
                'libekf:Peer','Invalid reference identities, geometry or frame.');
            selected=find(strcmp(b.ids,m.peer_id));
            assert(numel(selected)==1,'libekf:Peer','Unknown reference ID.');
            assert(~any(isfield(m,{'station','peer_covariance','peer_cross_cov','peer_time','peer_valid','peer_cross_frame_id'})), ...
                'libekf:Peer','Do not specify a reference twice.');
            peer_indices=3*(selected-1)+(1:3);
            m.station=b.position(peer_indices); m.peer_covariance=b.covariance;
            m.peer_cross_cov=b.cross_cov; m.peer_cross_frame_id=b.frame_id;
            m.peer_time=b.t; m.peer_valid=b.valid(selected); peer_position=b.position;
        else
            require_fields(m,{'station'}); peer_position=m.station(:);
        end
        require_fields(m,{'station','peer_covariance','peer_cross_cov','peer_cross_frame_id','peer_time','peer_valid'});
        assert(isequal(m.peer_valid,true)&&isscalar(m.peer_time)&&m.peer_time==s.t, ...
            'libekf:Peer','Peer position must be valid at measurement epoch.');
        assert(~isfield(m,'state_noise_cov'),'libekf:Correlation','Specify peer dependence only once.');
        V=m.peer_covariance; U=m.peer_cross_cov; d=numel(peer_position);
        assert(isequal(size(V),[d d])&&all(isfinite(V(:)))&& ...
            norm(V-V.','fro')<=1e-12*max(norm(V,'fro'),realmin)&& ...
            min(eig((V+V.')/2))>=0,'libekf:Covariance','Invalid peer covariance.');
        assert(any(size(U,1)==[21 42])&&size(U,2)==d&&all(isfinite(U(:)))&& ...
            strcmp(m.peer_cross_frame_id,s.frame_id), ...
            'libekf:Correlation','Invalid common-state/peer covariance.');
        delta=ps-vector3(m.station); distance=norm(delta);
        assert(distance>0,'libekf:Geometry','Undefined peer range.');
        direction=delta/distance; h=distance; H=direction.'*Jp; units={'m'};
        peer_jacobian=zeros(1,d); peer_jacobian(peer_indices)=-direction.';
        added_variance=peer_jacobian*V*peer_jacobian.';
        derived_cross=U*peer_jacobian.';
    case {'dme_range','pseudorange'}
        require_fields(m,{'station'});
        d=ps-vector3(m.station); rho=norm(d);
        assert(rho>0,'libekf:Geometry','Undefined range at station.');
        h=rho; H=d.'/rho*Jp; units={'m'};
        if strcmp(m.model,'pseudorange')
            require_fields(m,{'geometry_at_epoch','range_correction'});
            assert(isfield(m,'geometry_at_epoch') && isequal(m.geometry_at_epoch,true), ...
                'libekf:GeometryEpoch','Satellite transmit/receive geometry must be prepared.');
            h=h+s.clock+m.range_correction; H(20)=H(20)+1;
        end
    case {'vor_azimuth','station_range_azimuth'}
        require_fields(m,{'basis','station','azimuth_offset'});
        B=basis(m.basis); d=ps-vector3(m.station);
        east=B(:,1).'*d; north=B(:,2).'*d; r2=east^2+north^2;
        assert(r2>0,'libekf:Geometry','Undefined azimuth above station.');
        h=atan2(east,north)+m.azimuth_offset;
        H=(north*B(:,1).'-east*B(:,2).')/r2*Jp;
        units={'rad'}; angular=1;
        if strcmp(m.model,'station_range_azimuth')
            rho=norm(d); h=[rho;h]; H=[d.'/rho*Jp;H];
            units={'m';'rad'}; angular=2;
        end
    case {'optical_landmarks','celestial_pixels'}
        require_fields(m,{'calibration_valid','matches_valid','camera_to_body','intrinsics'});
        celestial=strcmp(m.model,'celestial_pixels');
        if celestial
            require_fields(m,{'directions','geometry_at_epoch','catalog_valid','min_forward'});
            points=m.directions; depth=m.min_forward;
            assert(isequal(m.geometry_at_epoch,true)&&isequal(m.catalog_valid,true)&& ...
                size(points,2)==3&&all(abs(sqrt(sum(points.^2,2))-1)<1e-10), ...
                'libekf:Optical','Invalid prepared celestial directions.');
        else
            require_fields(m,{'map_valid','landmarks','min_depth'});
            assert(isequal(m.map_valid,true),'libekf:Optical','Unusable landmark map.');
            points=m.landmarks; depth=m.min_depth;
        end
        assert(isequal(m.calibration_valid,true)&&isequal(m.matches_valid,true), ...
            'libekf:Optical','Unusable camera calibration or correspondence.');
        D=basis(m.camera_to_body);
        assert(abs(det(D)-1)<1e-10,'libekf:Optical','Camera axes must preserve handedness.');
        intr=m.intrinsics(:);
        assert(numel(intr)==4 && all(isfinite(intr)) && all(intr(1:2)>0) && ...
            size(points,2)==3 && size(points,1)>0 && all(isfinite(points(:))) && ...
            isscalar(depth)&&isfinite(depth)&&depth>0, ...
            'libekf:Optical','Invalid camera or landmark geometry.');
        n=size(points,1); h=zeros(2*n,1); H=zeros(2*n,21);
        for j=1:n
            if celestial, vb=C.'*points(j,:).'; else, vb=C.'*(points(j,:).'-ps); end
            vc=D.'*vb;
            assert(vc(3)>depth,'libekf:Optical','Landmark behind or too near camera.');
            rows=2*j-1:2*j;
            h(rows)=[intr(1)*vc(1)/vc(3)+intr(3);intr(2)*vc(2)/vc(3)+intr(4)];
            if celestial, J=zeros(3,21); else, J=-C.'*Jp; end
            J(:,7:9)=J(:,7:9)+nav_skew(vb);
            projection=[intr(1)/vc(3) 0 -intr(1)*vc(1)/vc(3)^2; ...
                0 intr(2)/vc(3) -intr(2)*vc(2)/vc(3)^2];
            H(rows,:)=projection*D.'*J;
        end
        units=repmat({'pixel'},2*n,1);
    case 'magnetic_vector'
        require_fields(m,{'field_valid','calibration_valid','field_gradient','field_radius','field_time','field_origin','field_vector'});
        assert(isfield(m,'field_valid') && isequal(m.field_valid,true) && ...
            isfield(m,'calibration_valid') && isequal(m.calibration_valid,true), ...
            'libekf:Magnetic','Unusable magnetic field or calibration.');
        G=m.field_gradient; radius=m.field_radius;
        assert(isequal(size(G),[3 3]) && all(isfinite(G(:))) && ...
            isscalar(radius)&&isfinite(radius)&&radius>0, ...
            'libekf:Magnetic','Invalid local magnetic model.');
        assert(isscalar(m.field_time)&&isfinite(m.field_time)&& ...
            abs(m.field_time-s.t)<=env.time_tolerance, ...
            'libekf:Epoch','Magnetic model must be prepared at measurement epoch.');
        dp=ps-vector3(m.field_origin);
        assert(norm(dp)<=radius,'libekf:Magnetic','Outside magnetic model domain.');
        field=vector3(m.field_vector)+G*dp;
        assert(norm(field)>0,'libekf:Magnetic','Zero reference magnetic field.');
        h=C.'*field; H=C.'*G*Jp;
        H(:,7:9)=H(:,7:9)+nav_skew(h);
        units=repmat({'T'},3,1);
    case 'adf_bearing'
        require_fields(m,{'bearing_geometry_valid','station','bearing_offset'});
        assert(isequal(m.bearing_geometry_valid,true),'libekf:Bearing','Uncalibrated bearing geometry.');
        vb=C.'*(vector3(m.station)-ps); x=vb(1); y=vb(2); r2=x*x+y*y;
        assert(r2>0,'libekf:Geometry','Relative bearing undefined on body vertical.');
        D=[-y/r2 x/r2 0]; J=-C.'*Jp; J(:,7:9)=J(:,7:9)+nav_skew(vb);
        h=atan2(y,x)+m.bearing_offset; H=D*J; units={'rad'}; angular=1;
    case {'range_difference','time_difference'}
        require_fields(m,{'geometry_at_epoch','station','reference_station','range_correction'});
        assert(isequal(m.geometry_at_epoch,true),'libekf:GeometryEpoch','Unprepared timing geometry.');
        da=ps-vector3(m.station); db=ps-vector3(m.reference_station);
        ra=norm(da); rb=norm(db);
        assert(ra>0&&rb>0,'libekf:Geometry','Undefined range difference at a station.');
        h=ra-rb+m.range_correction; H=(da.'/ra-db.'/rb)*Jp; units={'m'};
        if strcmp(m.model,'time_difference')
            require_fields(m,{'propagation_speed'});
            assert(isscalar(m.propagation_speed)&&isfinite(m.propagation_speed)&&m.propagation_speed>0, ...
                'libekf:Propagation','Invalid propagation speed.');
            h=h/m.propagation_speed; H=H/m.propagation_speed; units={'s'};
        end
    case 'rsdn_time_difference'
        require_fields(m,{'geometry_at_epoch','stations','propagation_speed','time_correction','reference_station'});
        assert(isequal(m.geometry_at_epoch,true),'libekf:GeometryEpoch','Unprepared RSDN timing geometry.');
        stations=m.stations; speed=m.propagation_speed;
        assert(size(stations,2)==3&&size(stations,1)>0&&all(isfinite(stations(:)))&& ...
            isscalar(speed)&&isfinite(speed)&&speed>0,'libekf:Geometry','Invalid RSDN stations.');
        correction=m.time_correction(:);
        assert(numel(correction)==size(stations,1)&&all(isfinite(correction)), ...
            'libekf:Time','Invalid transmitter and propagation corrections.');
        d=ps.'-stations; ranges=sqrt(sum(d.^2,2));
        dr=ps-vector3(m.reference_station); rr=norm(dr);
        assert(all(ranges>0)&&rr>0,'libekf:Geometry','RSDN geometry undefined at a station.');
        h=(ranges-rr)/speed+correction;
        H=(d./ranges-dr.'/rr)*Jp/speed;
        units=repmat({'s'},size(stations,1),1);
    case 'baro_height'
        require_fields(m,{'geoid_offset'});
        assert(isscalar(m.geoid_offset)&&isreal(m.geoid_offset)&&isfinite(m.geoid_offset), ...
            'libekf:Height','Invalid geoid offset.');
        B=basis(env.C_e_l);
        [height,up_e]=nav_height_local(ps,env);
        up=B.'*up_e;
        h=height-m.geoid_offset+s.bh; H=up.'*Jp; H(19)=H(19)+1; units={'m'};
    case {'static_pressure','impact_pressure','air_pressures'}
        require_fields(m,{'geoid_offset','atmosphere','pressure_valid','calibration_valid'});
        a=m.atmosphere;
        require_fields(a,{'height','pressure','temperature','gradient','gravity','gas_constant', ...
            'min_height','max_height','time_start','time_end','valid'});
        parameters=[a.height a.pressure a.temperature a.gradient a.gravity a.gas_constant ...
            a.min_height a.max_height a.time_start a.time_end];
        assert(numel(parameters)==10&&all(isfinite(parameters))&&isreal(parameters)&& ...
            a.pressure>0&&a.temperature>0&&a.gravity>0&&a.gas_constant>0&& ...
            a.min_height<a.max_height&&a.time_start<=s.t&&s.t<=a.time_end&& ...
            isequal(a.valid,true)&&isequal(m.pressure_valid,true)&&isequal(m.calibration_valid,true)&& ...
            isscalar(m.geoid_offset)&&isfinite(m.geoid_offset), ...
            'libekf:Atmosphere','Invalid pressure source or atmosphere domain.');
        [height,up_e]=nav_height_local(ps,env);
        altitude=height-m.geoid_offset+s.bh; delta=altitude-a.height;
        temperature=a.temperature+a.gradient*delta;
        assert(altitude>=a.min_height&&altitude<=a.max_height&&temperature>0&& ...
            numel(m.z)==1+strcmp(m.model,'air_pressures')&&isreal(m.z)&&all(isfinite(m.z(:)))&& ...
            (strcmp(m.model,'impact_pressure')||all(m.z(:)>0)), ...
            'libekf:Atmosphere','Pressure or height outside the atmosphere domain.');
        if a.gradient==0
            exponent=-a.gravity*delta/(a.gas_constant*a.temperature);
        else
            exponent=-a.gravity/(a.gas_constant*a.gradient)*log1p(a.gradient*delta/a.temperature);
        end
        h=a.pressure*exp(exponent);
        assert(isfinite(h)&&h>0,'libekf:Atmosphere','Pressure outside numerical range.');
        H=(env.C_e_l.'*up_e).'*Jp; H(19)=H(19)+1;
        height_jacobian=H;
        H=-a.gravity*h/(a.gas_constant*temperature)*H; units={'Pa'};
        if any(strcmp(m.model,{'impact_pressure','air_pressures'}))
            static_pressure=h; static_jacobian=H;
            require_fields(m,{'gamma','max_mach','flow_valid'});
            assert(isscalar(m.gamma)&&isfinite(m.gamma)&&m.gamma>1&& ...
                isscalar(m.max_mach)&&isfinite(m.max_mach)&&m.max_mach>0&&m.max_mach<1&& ...
                isequal(m.flow_valid,true),'libekf:Airflow','Invalid subsonic flow model.');
            omega=zeros(3,1); Omega=zeros(3,1);
            if norm(lever)>0
                require_fields(m,{'omega_ib'}); require_fields(env,{'Omega'});
                omega=vector3(m.omega_ib)-s.bg(:); Omega=vector3(env.Omega);
            end
            rot=cross(omega,lever);
            air=s.v(:)+C*rot-cross(Omega,C*lever)-s.wind(:);
            Jv=zeros(3,21); Jv(:,4:6)=eye(3); Jv(:,16:18)=-eye(3);
            Jv(:,7:9)=-C*nav_skew(rot)+nav_skew(Omega)*C*nav_skew(lever);
            Jv(:,13:15)=C*nav_skew(lever);
            speed2=air.'*air; mach2=speed2/(m.gamma*a.gas_constant*temperature);
            assert(isfinite(mach2)&&mach2<=m.max_mach^2, ...
                'libekf:Airflow','Speed exceeds the declared subsonic domain.');
            factor=(m.gamma-1)/(2*m.gamma*a.gas_constant);
            increment=factor*speed2/temperature; power=m.gamma/(m.gamma-1);
            ratio=expm1(power*log1p(increment));
            derivative=power*exp((power-1)*log1p(increment));
            H=ratio*H+h*derivative*factor*(2*air.'*Jv/temperature- ...
                speed2*a.gradient/temperature^2*height_jacobian);
            h=h*ratio;
            if strcmp(m.model,'air_pressures')
                h=[static_pressure;static_pressure+h];
                H=[static_jacobian;static_jacobian+H]; units={'Pa';'Pa'};
            end
        end
    case 'airdata'
        vb=C.'*(s.v(:)-s.wind(:)); speed=norm(vb); x=vb(1); y=vb(2); z=vb(3);
        r2=x*x+z*z; d=sqrt(r2);
        assert(speed>0 && d>0,'libekf:Geometry','Air angles undefined at zero/vertical lateral velocity.');
        h=[speed;atan2(z,x);atan2(y,d)];
        D=[vb.'/speed; -z/r2 0 x/r2; -y*x/(speed^2*d) d/speed^2 -y*z/(speed^2*d)];
        J=zeros(3,21); J(:,4:6)=C.'; J(:,7:9)=nav_skew(vb); J(:,16:18)=-C.';
        H=D*J; units={'m/s';'rad';'rad'}; angular=[2 3];
    case {'radar_range','radar_terrain','radar_footprint'}
        require_fields(m,{'beam_body','coverage_complete','terrain_quality','nearest_return_valid', ...
            'vertices','faces','half_angle','edge_margin','max_valid_range'});
        beam=unit(m.beam_body);
        assert(isfield(m,'coverage_complete') && isequal(m.coverage_complete,true) && ...
            isfield(m,'terrain_quality') && strcmp(m.terrain_quality,'VALID'), ...
            'libekf:TerrainCoverage','Complete, valid terrain coverage of the footprint is required.');
        assert(isfield(m,'nearest_return_valid') && isequal(m.nearest_return_valid,true), ...
            'libekf:RadarFootprint','Nearest-point return law must be validated for this source.');
        hit=nav_terrain_footprint(ps,C*beam,m.half_angle,m.vertices,m.faces,m.edge_margin);
        assert(hit.smooth,'libekf:TerrainNonsmooth','Closest return is not a unique smooth face interior.');
        assert(isscalar(m.max_valid_range) && isfinite(m.max_valid_range) && m.max_valid_range>0 && ...
            m.z<=m.max_valid_range && hit.range<=m.max_valid_range, ...
            'libekf:RadarRange','Radar return exceeds the configured valid range.');
        h=hit.range; H=-hit.normal.'/hit.direction_cosine*Jp;
        if hit.cone_active
            k=hit.normal.'*C*beam; transverse=sqrt(max(0,1-k*k));
            derivative=cos(m.half_angle)-k/transverse*sin(m.half_angle);
            H(7:9)=H(7:9)+hit.plane_distance/hit.direction_cosine^2*derivative*hit.normal.'*C*nav_skew(beam);
        end
        units={'m'};
    case 'ils_angles'
        require_fields(m,{'basis','station','angle_reference'});
        assert(numel(m.angle_reference)==2&&isreal(m.angle_reference)&&all(isfinite(m.angle_reference(:))), ...
            'libekf:Geometry','ILS needs two calibrated reference angles.');
        B=basis(m.basis); d=B.'*(ps-vector3(m.station));
        x=d(1); y=d(2); z=d(3); horizontal=hypot(x,y); distance2=d.'*d;
        assert(x>0 && horizontal>0, 'libekf:Geometry','Outside calibrated forward approach sector.');
        h=[atan2(y,x);atan2(z,horizontal)]-m.angle_reference(:);
        D=[-y/(x*x+y*y) x/(x*x+y*y) 0; ...
            -z*x/(distance2*horizontal) -z*y/(distance2*horizontal) horizontal/distance2];
        H=D*B.'*Jp; units={'rad';'rad'}; angular=[1 2];
    case 'gls_deviation'
        require_fields(m,{'projection','path_origin'});
        M=m.projection;
        assert(isequal(size(M),[2 3]) && all(isfinite(M(:))) && norm(M*M.'-eye(2),'fro')<1e-10, ...
            'libekf:Geometry','GLS projection needs orthonormal lateral/vertical rows.');
        h=M*(ps-vector3(m.path_origin)); H=M*Jp; units={'m';'m'};
    case 'attitude'
        nav_rotation(m.z); measured=m.z(:); q=s.q(:);
        % Log(q_nominal^{-1} * q_measured), invariant under quaternion sign.
        qe=[q(1)*measured(1)+q(2:4).'*measured(2:4); ...
            q(1)*measured(2:4)-measured(1)*q(2:4)-cross(q(2:4),measured(2:4))];
        if qe(1)<0, qe=-qe; end
        a=2*atan2(norm(qe(2:4)),qe(1));
        assert(isfield(m,'max_angle') && isscalar(m.max_angle) && isfinite(m.max_angle) && ...
            m.max_angle>0 && m.max_angle<pi && a<=m.max_angle, ...
            'libekf:AttitudeDomain','Attitude discrepancy exceeds validated local chart.');
        if a<1e-8, z=2*qe(2:4); coeff=1/12+a*a/720;
        else, z=a*qe(2:4)/norm(qe(2:4)); coeff=(1-0.5*a*cot(a/2))/(a*a); end
        A=nav_skew(z);
        H=zeros(3,21); H(:,7:9)=eye(3)-0.5*A+coeff*A*A;
        % Packet R describes right-local noise at the measured attitude.
        Rtransform=eye(3)+0.5*A+coeff*A*A;
        h=zeros(3,1); units=repmat({'rad'},3,1);
    otherwise
        error('libekf:Model','Unsupported observation model: %s',m.model);
end
if ~strcmp(m.model,'attitude'), z=m.z(:); end
units=units(:);
assert(numel(z)==numel(h) && all(isfinite(z)) && all(isfinite(h)) && all(isfinite(H(:))), ...
    'libekf:Dimension','Invalid observation dimension/value.');
assert(iscell(m.units) && isequal(m.units(:),units), 'libekf:Units','Observation units differ from model.');
R=m.R;
assert(isequal(size(R),[numel(h) numel(h)]) && all(isfinite(R(:))) && ...
    norm(R-R.','fro')<=1e-12*max(norm(R,'fro'),realmin), 'libekf:Covariance','R must be finite and symmetric.');
[~,flag]=chol(R);
assert(flag==0,'libekf:Covariance','This observation profile requires positive definite R.');
R=R+added_variance;
if ~isempty(Rtransform), R=Rtransform*R*Rtransform.'; end
residual=z-h;
residual(angular)=atan2(sin(residual(angular)),cos(residual(angular)));
obs=struct('predicted',h,'residual',residual,'H',H,'R',R,'units',{units}, ...
    'source',m.source,'model',m.model,'t_measure',m.t_measure,'t_received',m.t_received, ...
    'frame_id',m.frame_id, ...
    'lineage',{m.lineage},'correlation_group',m.correlation_group,'noise_model',m.noise_model);
if ~isempty(derived_cross)
    obs.state_noise_cov=derived_cross;
    obs.peer_covariance=m.peer_covariance;
    obs.peer_cross_cov=m.peer_cross_cov;
    obs.peer_jacobian=peer_jacobian;
    obs.peer_noise_cov=m.R;
    obs.peer_position=peer_position;
    if isfield(m,'reference_bank'), obs.peer_ids=m.reference_bank.ids; end
    obs.peer_time=m.peer_time;
end
if isfield(m,'state_noise_cov')
    U=m.state_noise_cov;
    assert(isfield(m,'state_noise_frame_id')&&strcmp(m.state_noise_frame_id,s.frame_id)&& ...
        any(size(U,1)==[21 42])&&size(U,2)==numel(h)&&all(isfinite(U(:))), ...
        'libekf:Correlation','Cross covariance rows must use common-state error axes.');
    if ~isempty(Rtransform), U=U*Rtransform.'; end
    obs.state_noise_cov=U;
end

end

function v=vector3(v)
assert(isnumeric(v) && numel(v)==3 && all(isfinite(v(:))), 'libekf:Vector','Expected finite 3-vector.'); v=v(:);
end
function u=unit(u)
u=vector3(u); assert(abs(norm(u)-1)<1e-10,'libekf:Geometry','Expected unit direction.');
end
function B=basis(B)
assert(isequal(size(B),[3 3]) && all(isfinite(B(:))) && norm(B.'*B-eye(3),'fro')<1e-10, ...
    'libekf:Geometry','Expected orthonormal geometry basis.');
end

function require_fields(packet,names)
for k=1:numel(names)
    assert(isfield(packet,names{k}),'libekf:Packet',['Missing model field ' names{k}]);
end
end

function [out,result]=nav_assimilate(filter,packet,sensor,env,limit,reception)
% Инерциальный интервал либо наблюдение из местных осей прибора в общую оценку.
% reception (для обработки в реальном времени): now, max_age в секундах.
% Без reception выполняется расчёт по времени измерений без проверки текущего возраста.
% filter: state, P, used (исходные ID, включая использованные при инициализации), capacity.
assert(iscell(filter.used) && all(cellfun(@(x)ischar(x)&&~isempty(x),filter.used)) && ...
    isscalar(filter.capacity) && isfinite(filter.capacity) && ...
    filter.capacity>=numel(filter.used) && filter.capacity==floor(filter.capacity), ...
    'libekf:History','Invalid source history or capacity.');
out=filter;
result=struct('accepted',false,'reason','','nis',NaN);
try
    % Проверяем физическую модель и договор измерения до изменения состояния.
    assert(isstruct(packet)&&isscalar(packet)&&isfield(packet,'model')&& ...
        ischar(packet.model)&&~isempty(packet.model),'libekf:Packet','Missing physical model.');
    if isfield(filter,'retired_through')||isfield(packet,'lineage_times')
        epoch_field='t_measure';
        if any(strcmp(packet.model,{'inertial','visual_odometry'})), epoch_field='t_end'; end
        assert(all(isfield(packet,{'lineage','lineage_times',epoch_field}))&& ...
            iscell(packet.lineage)&&~isempty(packet.lineage)&& ...
            isnumeric(packet.lineage_times)&&isreal(packet.lineage_times)&& ...
            isvector(packet.lineage_times)&&numel(packet.lineage_times)==numel(packet.lineage)&& ...
            all(isfinite(packet.lineage_times(:)))&&isnumeric(packet.(epoch_field))&& ...
            isreal(packet.(epoch_field))&&isscalar(packet.(epoch_field))&&isfinite(packet.(epoch_field))&& ...
            all(packet.lineage_times(:)<=packet.(epoch_field)), ...
            'libekf:LineageTime','Provide valid original record epochs, not packet arrival times.');
        if isfield(filter,'retired_through')
            boundary=filter.retired_through;
            assert(isnumeric(boundary)&&isreal(boundary)&&isscalar(boundary)&&isfinite(boundary)&& ...
                boundary<=filter.state.t,'libekf:History','Invalid retired history boundary.');
            assert(all(packet.lineage_times(:)>boundary),'libekf:RetiredHistory', ...
                'Observation reuses original records from closed history.');
        end
    end
    if nargin>=6
        assert(isstruct(reception)&&isscalar(reception)&& ...
            isfield(reception,'now')&&isfield(reception,'max_age')&& ...
            isscalar(reception.now)&&isreal(reception.now)&&isfinite(reception.now)&& ...
            isscalar(reception.max_age)&&isreal(reception.max_age)&& ...
            isfinite(reception.max_age)&&reception.max_age>=0, ...
            'libekf:Time','Invalid processing time or permitted age.');
        epoch_field='t_measure';
        if any(strcmp(packet.model,{'inertial','visual_odometry'})), epoch_field='t_end'; end
        assert(isfield(packet,epoch_field)&&isfield(packet,'t_received'), ...
            'libekf:Time','Missing observation or arrival time.');
        epoch=packet.(epoch_field); received=packet.t_received;
        assert(isscalar(epoch)&&isreal(epoch)&&isfinite(epoch)&& ...
            isscalar(received)&&isreal(received)&&isfinite(received)&& ...
            epoch<=received&&received<=reception.now, ...
            'libekf:Time','Observation is not available at processing time.');
        assert(reception.now-epoch<=reception.max_age, ...
            'libekf:Stale','Observation age exceeds the source limit.');
    end
    paired=isfield(filter,'reference');
    if paired
        assert(isfield(filter,'jointP')&&isequal(size(filter.jointP),[42 42])&& ...
            all(isfinite(filter.jointP(:)))&& ...
            isequal(filter.P,filter.jointP(22:42,22:42))&& ...
            strcmp(filter.reference.frame_id,filter.state.frame_id)&& ...
            filter.reference.t<=filter.state.t,'libekf:State','Invalid saved reference state.');
        J=filter.jointP; d=sqrt(diag(J));
        assert(all(isreal(d)&isfinite(d)&d>0)&& ...
            norm(J-J.','fro')<=1e-12*max(norm(J,'fro'),realmin), ...
            'libekf:Covariance','Invalid joint covariance.');
        scaled=(J./d)./d.'; ev=eig((scaled+scaled.')/2);
        assert(min(ev)>=-100*eps*42*max(1,max(abs(ev))), ...
            'libekf:Covariance','Indefinite joint covariance.');
        reference=filter.reference;
    end
    odometry=strcmp(packet.model,'visual_odometry');
    inertial=strcmp(packet.model,'inertial');
    source_health=strcmp(packet.model,'source_health');
    health=strcmp(packet.model,'reference_health')||source_health;
    excluded={};
    if isfield(filter,'excluded_sources'), excluded=filter.excluded_sources; end
    assert(iscell(excluded)&&all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),excluded)), ...
        'libekf:SourceHealth','Invalid excluded source list.');
    assert(~isfield(packet,'source')||~any(strcmp(excluded,packet.source)), ...
        'libekf:SourceHealth','Source excluded by a diagnostic report.');
    if isfield(packet,'input_sources')
        dependencies=packet.input_sources;
        assert(iscell(dependencies)&&all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),dependencies)), ...
            'libekf:SourceHealth','Invalid input source identities.');
        assert(~any(ismember(dependencies(:),excluded(:))), ...
            'libekf:SourceHealth','Observation depends on an excluded source.');
    end
    assert(~inertial||~isfield(packet,'retained_reference'),'libekf:Peer', ...
        'Retained reference correction is not an inertial propagation model.');
    stored_noise=isfield(filter,'future_noise');
    if stored_noise
        future=filter.future_noise; prior=filter.P; epochs=filter.state.t;
        if paired, prior=J; epochs=[filter.reference.t;filter.state.t]; end
        assert(~isfield(packet,'future_noise')&&isfield(future,'state_times')&& ...
            isfield(future,'frame_id')&&isequal(future.state_times(:),epochs(:))&& ...
            strcmp(future.frame_id,filter.state.frame_id), ...
            'libekf:Correlation','Stored errors cannot be overridden or used at another epoch/frame.');
        future=nav_predict_noise(future,prior,eye(size(prior)),zeros(size(prior)),epochs, ...
            zeros(size(prior,1),numel(future.mean)));
    end
    stored_bank=isfield(filter,'reference_bank');
    if stored_bank
        bank=filter.reference_bank;
        assert(isstruct(bank)&&isscalar(bank)&&all(isfield(bank,{'t','frame_id'}))&& ...
            isscalar(bank.t)&&bank.t==filter.state.t&&strcmp(bank.frame_id,filter.state.frame_id), ...
            'libekf:Peer','Stored reference set must match current state.');
        prior=filter.P; if paired, prior=J; end
        check_bank(bank,prior);
        assert(~isfield(packet,'reference_bank')&&~isfield(packet,'retained_reference'), ...
            'libekf:Peer','Stored reference set cannot be overridden by a packet.');
        if ~inertial&&~health
            if strcmp(packet.model,'cooperative_range')
                packet.reference_bank=bank;
            else
                assert(isfield(packet,'reference_noise_independent')&& ...
                    (isequal(packet.reference_noise_independent,true)||isequal(packet.reference_noise_independent,false)), ...
                    'libekf:Correlation','Reference noise dependence must be explicit.');
                bank.noise_independent=packet.reference_noise_independent;
                if ~bank.noise_independent
                    assert(isfield(packet,'reference_noise_cross_cov'),'libekf:Correlation','Missing reference/noise covariance.');
                end
                if isfield(packet,'reference_noise_cross_cov')
                    bank.noise_cross_cov=packet.reference_noise_cross_cov;
                end
                packet.retained_reference=bank;
            end
        end
    end
    if health
        target='peer_id'; if source_health, target='target_source'; end
        assert((stored_bank||source_health)&&all(isfield(packet,{'source','lineage','quality',target,'valid', ...
            't_measure','t_received','frame_id'})), 'libekf:PeerHealth','Incomplete health report.');
        assert(ischar(packet.source)&&~isempty(packet.source)&& ...
            iscell(packet.lineage)&&~isempty(packet.lineage)&& ...
            all(cellfun(@(x)ischar(x)&&~isempty(x),packet.lineage))&& ...
            strcmp(packet.quality,'VALID')&&ischar(packet.(target))&&isrow(packet.(target))&&~isempty(packet.(target))&& ...
            islogical(packet.valid)&&isscalar(packet.valid)&&strcmp(packet.frame_id,filter.state.frame_id)&& ...
            isnumeric(packet.t_measure)&&isscalar(packet.t_measure)&&isreal(packet.t_measure)&& ...
            packet.t_measure==filter.state.t&&isnumeric(packet.t_received)&& ...
            isscalar(packet.t_received)&&isreal(packet.t_received)&&isfinite(packet.t_received)&& ...
            packet.t_received>=packet.t_measure, 'libekf:PeerHealth','Invalid reference health report.');
        if source_health
            selected=strcmp(excluded,packet.target_source);
            if packet.valid
                excluded(selected)=[];
            elseif ~any(selected)
                excluded{end+1}=packet.target_source;
            end
        else
            selected=find(strcmp(bank.ids,packet.peer_id));
            assert(numel(selected)==1,'libekf:PeerHealth','Unknown reference health ID.');
            bank.valid(selected)=packet.valid;
        end
        candidate=filter.state; Q=filter.P; ids=unique(packet.lineage(:));
    elseif inertial
        [candidate,Q,T,N]=nav_predict(filter.state,filter.P,packet,env);
        if paired
            cross=J(1:21,22:42)*T.';
            J=[J(1:21,1:21) cross;cross.' Q];
        end
        if stored_bank
            assert(isfield(packet,'reference_motion'),'libekf:Peer','Reference motion model required.');
            transition=T; if paired, transition=blkdiag(eye(21),T); end
            bank=nav_predict_reference(bank,transition,packet.reference_motion);
            assert(bank.t==candidate.t,'libekf:Epoch','Reference and navigation propagation end times differ.');
        end
        if stored_noise
            assert(isfield(packet,'future_noise_independent')&& ...
                (isequal(packet.future_noise_independent,true)||isequal(packet.future_noise_independent,false)), ...
                'libekf:Correlation','Declare new process/stored-noise dependence.');
            transition=T; process=N; epochs=candidate.t; prior=filter.P;
            if paired
                transition=blkdiag(eye(21),T); process=blkdiag(zeros(21),N);
                prior=filter.jointP; epochs=[filter.reference.t;candidate.t];
            end
            W=zeros(size(transition,1),numel(future.mean));
            if isfield(packet,'future_noise_process_cross_cov'), W=packet.future_noise_process_cross_cov; end
            assert((packet.future_noise_independent&&all(W(:)==0))|| ...
                (~packet.future_noise_independent&&isfield(packet,'future_noise_process_cross_cov')), ...
                'libekf:Correlation','Missing or contradictory process/stored-noise covariance.');
            future=nav_predict_noise(future,prior,transition,process,epochs,W);
        end
        ids=unique(packet.lineage(:));
    elseif odometry
        assert(paired,'libekf:State','Odometry requires a saved reference state.');
        obs=nav_odometry_measurement(reference,filter.state,packet);
        ids=unique(obs.lineage(:));
    else
        assert(isfield(packet,'t_measure')&&isscalar(packet.t_measure)&&isfinite(packet.t_measure), ...
            'libekf:Time','Missing measurement epoch.');
        at_reference=paired&&packet.t_measure==reference.t&&reference.t<filter.state.t;
        if at_reference
            obs=nav_observe_source(reference,packet,sensor,env);
        else
            obs=nav_observe_source(filter.state,packet,sensor,env);
        end
        ids=unique(obs.lineage(:));
        if paired
            if at_reference, obs.H=[obs.H zeros(size(obs.H,1),21)];
            else, obs.H=[zeros(size(obs.H,1),21) obs.H]; end
            obs.t_start=reference.t; obs.t_end=filter.state.t;
        end
    end
    if isfield(packet,'shared_noise')
        assert(stored_noise&&~inertial&&~health&&~isfield(packet,'reuse'), ...
            'libekf:Correlation','Shared errors require the stored model and automatic conditioning.');
        obs=shared_observation(obs,future,packet);
    elseif stored_noise&&~inertial&&~health
        assert(isfield(packet,'future_noise_independent')&& ...
            (isequal(packet.future_noise_independent,true)||isequal(packet.future_noise_independent,false)), ...
            'libekf:Correlation','Declare measurement/stored-noise dependence.');
        D=zeros(numel(future.mean),numel(obs.residual));
        if isfield(packet,'future_noise_cross_cov'), D=packet.future_noise_cross_cov; end
        assert((packet.future_noise_independent&&all(D(:)==0))|| ...
            (~packet.future_noise_independent&&isfield(packet,'future_noise_cross_cov')), ...
            'libekf:Correlation','Missing or contradictory measurement/stored-noise covariance.');
        future.measurement_cross_cov=D; obs.future_noise=future;
    end
    if isfield(packet,'future_noise')
        assert(~inertial&&~health,'libekf:Correlation','Future-noise conditioning requires an observation.');
        obs.future_noise=packet.future_noise;
    end
    overlap=intersect(ids,filter.used(:));
    if ~isempty(overlap)&&isfield(packet,'shared_noise')
        assert(isequal(sort(packet.shared_noise.lineage(:)),overlap(:)), ...
            'libekf:Reuse','Declared common errors do not cover all reused records.');
        epochs=filter.state.t; if paired, epochs=[reference.t;filter.state.t]; end
        packet.reuse=struct('lineage',{overlap},'history',{filter.used}, ...
            'state_times',epochs,'noise_mean',zeros(numel(obs.residual),1));
    end

    if ~isempty(overlap)
        fresh=ids(~ismember(ids,filter.used));
        if ~isfield(packet,'reuse')||isempty(fresh)
            result.reason='Повторное использование исходного отсчёта'; return;
        end
        % Условное наблюдение содержит новую информацию и явно связано с текущей историей.
        reuse=packet.reuse;
        assert(~inertial&&~health&&strcmp(obs.noise_model,'joint_required')&& ...
            isfield(obs,'state_noise_cov')&&isstruct(reuse)&&isscalar(reuse)&& ...
            all(isfield(reuse,{'lineage','history','state_times','noise_mean'})), ...
            'libekf:Reuse','Conditional measurement mean and covariance required.');
        assert(iscell(reuse.lineage)&&iscell(reuse.history)&& ...
            all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),reuse.lineage))&& ...
            all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),reuse.history))&& ...
            isequal(sort(reuse.lineage(:)),overlap(:))&& ...
            isequal(reuse.history(:),filter.used(:)), ...
            'libekf:Reuse','Shared records or conditioning history differ from the current filter.');
        epochs=filter.state.t; if paired, epochs=[reference.t;filter.state.t]; end
        mean_noise=reuse.noise_mean;
        assert(isnumeric(reuse.state_times)&&isreal(reuse.state_times)&& ...
            isequal(reuse.state_times(:),epochs(:))&&isnumeric(mean_noise)&&isreal(mean_noise)&& ...
            isvector(mean_noise)&&numel(mean_noise)==numel(obs.residual)&&all(isfinite(mean_noise(:))), ...
            'libekf:Reuse','Invalid conditioning epochs or residual noise mean.');
        obs.residual=obs.residual-mean_noise(:);
        ids=fresh; % Общие исходные записи не добавляются второй раз.
    elseif isfield(packet,'reuse')
        error('libekf:Reuse','Conditional reuse declared without a shared record.');
    end
    if numel(filter.used)+numel(ids)>filter.capacity
        result.reason='Исчерпана память происхождения измерений'; return;
    end
    if ~inertial&&~health
        retained_args={};
        if isfield(packet,'retained_reference'), retained_args={packet.retained_reference}; end
        if paired
            [states,J,info]=nav_correct({reference,filter.state},J,obs,limit,retained_args{:});
            reference=states{1}; candidate=states{2}; Q=J(22:42,22:42);
        else
            [candidate,Q,info]=nav_correct(filter.state,filter.P,obs,limit,retained_args{:});
        end
        result.nis=info.nis;
        if ~info.accepted
            result.reason='Превышено допустимое расхождение'; return;
        end
    end
    if ~health
        source_history=[];
        if isfield(filter,'source_history')
            source_history=nav_sources(filter.source_history);
            assert(all([source_history.last_epoch]<=filter.state.t), ...
                'libekf:SourceHistory','Source history is newer than the input state.');
        end
        source_history=nav_sources(source_history,packet);
    end
    if ~inertial&&~health&&isfield(info,'peer_update'), result.peer_update=info.peer_update; end
    if ~inertial&&~health&&isfield(info,'future_noise_update')
        result.future_noise_update=info.future_noise_update; out.future_noise=info.future_noise_update;
    elseif stored_noise
        out.future_noise=future;
    end
    if stored_bank
        if ~inertial&&~health
            update=info.peer_update;
            for name={'position','covariance','cross_cov','t','frame_id'}
                bank.(name{1})=update.(name{1});
            end
            for name={'noise_independent','noise_cross_cov'}
                if isfield(bank,name{1}), bank=rmfield(bank,name{1}); end
            end
        end
        out.reference_bank=bank;
    end
    if paired, out.reference=reference; out.jointP=J; end
    out.state=candidate; out.P=Q; out.used=[filter.used(:);ids];
    if ~health, out.source_history=source_history; end
    if inertial
        transition=T; process_covariance=N;
        if paired, transition=blkdiag(eye(21),T); process_covariance=blkdiag(zeros(21),N); end
        result.propagation=struct('transition',transition,'process_covariance',process_covariance, ...
            't_start',packet.t_start,'t_end',packet.t_end,'frame_id',filter.state.frame_id, ...
            'source',packet.source,'lineage',{packet.lineage});
        if stored_noise, result.propagation.future_noise_process_cross_cov=W; end
    elseif ~health
        result.observation=obs;
        result.error_reset=info.error_reset;
        result.state_times=filter.state.t;
        if paired, result.state_times=[filter.reference.t;filter.state.t]; end
    end
    result.accepted=true; result.reason='Оценка обновлена; достоверность проверяет сторож';
    if health, result.reason='Признак исправности опоры обновлён'; end
    if source_health
        out.excluded_sources=excluded;
        result.reason='Допуск источника обновлён по сообщению диагностики';
    end
catch ex
    if startsWith(ex.identifier,'libekf:')
        result.reason=ex.identifier;
    else
        rethrow(ex);
    end
end
end

function check_bank(bank,P)
assert(all(isfield(bank,{'ids','valid','position','covariance','cross_cov'}))&& ...
    iscell(bank.ids)&&~isempty(bank.ids)&& ...
    all(cellfun(@(x)ischar(x)&&~isempty(x),bank.ids))&& ...
    numel(unique(bank.ids))==numel(bank.ids)&& ...
    islogical(bank.valid)&&numel(bank.valid)==numel(bank.ids), ...
    'libekf:Peer','Invalid stored reference identities or health flags.');
d=3*numel(bank.ids); n=size(P,1);
assert(any(n==[21 42])&&isequal(size(P),[n n])&&isreal(P)&&all(isfinite(P(:)))&& ...
    isequal(size(bank.position),[d 1])&&isreal(bank.position)&&all(isfinite(bank.position))&& ...
    isequal(size(bank.covariance),[d d])&&isreal(bank.covariance)&&all(isfinite(bank.covariance(:)))&& ...
    isequal(size(bank.cross_cov),[n d])&&isreal(bank.cross_cov)&&all(isfinite(bank.cross_cov(:))), ...
    'libekf:Covariance','Invalid stored reference matrices.');
joint=[P bank.cross_cov;bank.cross_cov.' bank.covariance];
assert(norm(joint-joint.','fro')<=1e-12*max(norm(joint,'fro'),realmin)&&all(diag(joint)>=0), ...
    'libekf:Covariance','Invalid stored joint covariance.');
scales=sqrt(diag(joint)); zero=scales==0;
assert(all(all(joint(zero,:)==0)),'libekf:Correlation','Zero variance has nonzero reference dependence.');
scales(zero)=1; normalized=(joint./scales)./scales.';
values=eig((normalized+normalized.')/2);
assert(min(values)>=-100*eps*(n+d)*max(1,max(abs(values))), ...
    'libekf:Correlation','Stored navigation/reference covariance is inconsistent.');
end


function obs=shared_observation(obs,future,packet)
% R пакета относится только к новой независимой части шума в его модели.
shared=packet.shared_noise; m=numel(obs.residual); n=size(obs.H,2); b=numel(future.mean);
residual_frame=obs.frame_id; if isfield(obs,'source_frame_id'), residual_frame=obs.source_frame_id; end
assert(isstruct(shared)&&isscalar(shared)&&all(isfield(shared,{'ids','units','map','residual_units','residual_frame_id','lineage'}))&& ...
    isfield(packet,'future_noise_independent')&&isequal(packet.future_noise_independent,true)&& ...
    ~isfield(packet,'future_noise_cross_cov')&&~isfield(obs,'peer_covariance')&& ...
    strcmp(obs.noise_model,'independent')&& ...
    (~isfield(obs,'state_noise_cov')||all(obs.state_noise_cov(:)==0))&& ...
    (~isfield(obs,'state_input_cov')||all(obs.state_input_cov(:)==0)), ...
    'libekf:Correlation','Shared model requires a fresh noise part independent of state and stored errors.');
assert(iscell(shared.ids)&&~isempty(shared.ids)&& ...
    all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),shared.ids))&& ...
    numel(unique(shared.ids))==numel(shared.ids)&&iscell(shared.units)&& ...
    numel(shared.units)==numel(shared.ids)&&iscell(shared.residual_units)&& ...
    isequal(shared.residual_units(:),obs.units(:))&&strcmp(shared.residual_frame_id,residual_frame)&& ...
    iscell(shared.lineage)&& ...
    all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),shared.lineage))&& ...
    numel(unique(shared.lineage))==numel(shared.lineage)&&all(ismember(shared.lineage,packet.lineage)), ...
    'libekf:Correlation','Invalid common-error identities, units or provenance.');
[found,indices]=ismember(shared.ids(:),future.ids(:)); stored_units=future.units(:);
assert(all(found)&&isequal(shared.units(:),stored_units(indices(:)))&& ...
    isnumeric(shared.map)&&isreal(shared.map)&&isequal(size(shared.map),[m numel(indices)])&& ...
    all(isfinite(shared.map(:))),'libekf:Correlation','Unknown common error or invalid residual map.');
M=zeros(m,b); M(:,indices)=shared.map;
% Сохраняем полную входную модель, чтобы проекция не скрыла невозможную связь.
raw=obs.R; raw_state=zeros(n,m); transform=eye(m);
if any(isfield(obs,{'input_noise_cov','state_input_cov','noise_transform'}))
    assert(all(isfield(obs,{'input_noise_cov','state_input_cov','noise_transform'})), ...
        'libekf:Correlation','Incomplete fresh input noise model.');
    raw=obs.input_noise_cov; raw_state=obs.state_input_cov; transform=obs.noise_transform;
end
obs.input_noise_cov=blkdiag(raw,future.covariance);
obs.state_input_cov=[raw_state future.cross_cov]; obs.noise_transform=[transform M];
obs.shared_noise_model=struct('ids',{future.ids},'units',{future.units}, ...
    'map',M,'fresh_covariance',obs.R);
obs.conditional_noise_mean=M*future.mean(:);
obs.residual=obs.residual-obs.conditional_noise_mean;
obs.R=obs.R+M*future.covariance*M.'; obs.R=(obs.R+obs.R.')/2;
obs.state_noise_cov=future.cross_cov*M.'; obs.noise_model='joint_required';
future.measurement_cross_cov=future.covariance*M.'; obs.future_noise=future;
end

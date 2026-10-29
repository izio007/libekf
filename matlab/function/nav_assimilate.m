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
    assert(~inertial||~isfield(packet,'retained_reference'),'libekf:Peer', ...
        'Retained reference correction is not an inertial propagation model.');
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
        if ~inertial
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
    if inertial
        [candidate,Q,T]=nav_predict(filter.state,filter.P,packet,env);
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
    if any(ismember(ids,filter.used))
        result.reason='Повторное использование исходного отсчёта'; return;
    end
    if numel(filter.used)+numel(ids)>filter.capacity
        result.reason='Исчерпана память происхождения измерений'; return;
    end
    if ~inertial
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
    if ~inertial&&isfield(info,'peer_update'), result.peer_update=info.peer_update; end
    if stored_bank
        if ~inertial
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
    result.accepted=true; result.reason='Оценка обновлена; достоверность проверяет сторож';
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

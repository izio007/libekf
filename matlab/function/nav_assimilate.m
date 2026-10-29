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
    if inertial
        [candidate,Q,T]=nav_predict(filter.state,filter.P,packet,env);
        if paired
            cross=J(1:21,22:42)*T.';
            J=[J(1:21,1:21) cross;cross.' Q];
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

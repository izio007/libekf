function [out,result]=nav_assimilate(filter,packet,sensor,env,limit)
% Инерциальный интервал либо наблюдение из местных осей прибора в общую оценку.
% filter: state, P, used (исходные ID, включая использованные при инициализации), capacity.
assert(iscell(filter.used) && all(cellfun(@(x)ischar(x)&&~isempty(x),filter.used)) && ...
    isscalar(filter.capacity) && isfinite(filter.capacity) && ...
    filter.capacity>=numel(filter.used) && filter.capacity==floor(filter.capacity), ...
    'libekf:History','Invalid source history or capacity.');
out=filter;
result=struct('accepted',false,'reason','','nis',NaN);
try
    % Проверяем физическую модель и договор измерения до изменения состояния.
    assert(isfield(packet,'model'),'libekf:Packet','Missing physical model.');
    inertial=strcmp(packet.model,'inertial');
    if inertial
        [candidate,Q]=nav_predict(filter.state,filter.P,packet,env);
        ids=unique(packet.lineage(:));
    else
        obs=nav_observe_source(filter.state,packet,sensor,env);
        ids=unique(obs.lineage(:));
    end
    if any(ismember(ids,filter.used))
        result.reason='Повторное использование исходного отсчёта'; return;
    end
    if numel(filter.used)+numel(ids)>filter.capacity
        result.reason='Исчерпана память происхождения измерений'; return;
    end
    if ~inertial
        [candidate,Q,info]=nav_correct(filter.state,filter.P,obs,limit);
        result.nis=info.nis;
        if ~info.accepted
            result.reason='Превышено допустимое расхождение'; return;
        end
    end
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

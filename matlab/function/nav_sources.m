function history=nav_sources(history,packet)
% Сведения об источниках, фактически внесших принятые данные; не оценка исправности.
fields={'id','models','first_epoch','last_epoch'};
if isempty(history), history=struct('id',{},'models',{},'first_epoch',{},'last_epoch',{}); end
assert(isstruct(history)&&all(isfield(history,fields))&&numel(fieldnames(history))==4, ...
    'libekf:SourceHistory','Invalid source contribution history.');
for k=1:numel(history)
    h=history(k);
    assert(ischar(h.id)&&isrow(h.id)&&~isempty(h.id)&&iscell(h.models)&&~isempty(h.models)&& ...
        all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),h.models))&& ...
        isnumeric(h.first_epoch)&&isreal(h.first_epoch)&&isscalar(h.first_epoch)&&isfinite(h.first_epoch)&& ...
        isnumeric(h.last_epoch)&&isreal(h.last_epoch)&&isscalar(h.last_epoch)&&isfinite(h.last_epoch)&& ...
        h.first_epoch<=h.last_epoch, 'libekf:SourceHistory','Invalid source identity or time span.');
end
assert(numel(unique({history.id}))==numel(history),'libekf:SourceHistory','Duplicate source identity.');
if nargin>=2
    assert(isstruct(packet)&&isscalar(packet)&&all(isfield(packet,{'source','model'}))&& ...
        ischar(packet.source)&&isrow(packet.source)&&~isempty(packet.source)&& ...
        ischar(packet.model)&&isrow(packet.model)&&~isempty(packet.model), ...
        'libekf:SourceHistory','Invalid accepted source packet.');
    if any(strcmp(packet.model,{'inertial','visual_odometry'}))
        first=packet.t_start; last=packet.t_end;
    else
        first=packet.t_measure; last=first;
    end
    assert(isscalar(first)&&isscalar(last)&&isreal(first)&&isreal(last)&& ...
        isfinite(first)&&isfinite(last)&&first<=last,'libekf:SourceHistory','Invalid contribution epoch.');
    contributors={packet.source};
    if isfield(packet,'input_sources')
        assert(iscell(packet.input_sources)&& ...
            all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),packet.input_sources)), ...
            'libekf:SourceHistory','Invalid physical input source identities.');
        contributors=unique([contributors(:);packet.input_sources(:)]);
    end
    for j=1:numel(contributors)
        id=contributors{j}; index=find(strcmp({history.id},id));
        if isempty(index)
            history(end+1)=struct('id',id,'models',{{packet.model}},'first_epoch',first,'last_epoch',last);
        else
            history(index).models=unique([history(index).models(:);{packet.model}]);
            history(index).first_epoch=min(history(index).first_epoch,first);
            history(index).last_epoch=max(history(index).last_epoch,last);
        end
    end
end
[~,order]=sort({history.id}); history=history(order);
end

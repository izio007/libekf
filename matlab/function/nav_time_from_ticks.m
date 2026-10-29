function seconds=nav_time_from_ticks(ticks,origin,tick_seconds)
% Целочисленные метки одного часового отсчёта -> секунды от общего начала.
% Сначала вычитаем в uint64; абсолютная метка не проходит через double.
assert(isa(ticks,'uint64')&&isvector(ticks)&&~isempty(ticks)&& ...
    isa(origin,'uint64')&&isscalar(origin)&&isa(tick_seconds,'double')&& ...
    isreal(tick_seconds)&&isscalar(tick_seconds)&&isfinite(tick_seconds)&&tick_seconds>0, ...
    'libekf:Time','Expected original uint64 ticks, uint64 origin and positive seconds per tick.');
negative=ticks<origin; difference=zeros(size(ticks),'uint64');
difference(negative)=origin-ticks(negative);
difference(~negative)=ticks(~negative)-origin;
assert(all(difference(:)<=uint64(flintmax)),'libekf:TimePrecision', ...
    'Choose a closer common time origin before converting to double.');
seconds=double(difference)*tick_seconds;
assert(all(isfinite(seconds(:)))&&all(eps(seconds(:))<=tick_seconds), ...
    'libekf:TimePrecision','Relative seconds cannot resolve the original clock tick.');
seconds(negative)=-seconds(negative);
end

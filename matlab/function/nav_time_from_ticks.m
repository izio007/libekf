function seconds=nav_time_from_ticks(ticks,origin,tick_seconds)
% Целочисленные метки одного часового отсчёта -> секунды от общего начала.
% Сначала вычитаем в uint64; абсолютная метка не проходит через double.
% Допускаются также десятичные строки из файла: char, string или cellstr.
ticks=decimal_ticks(ticks); origin=decimal_ticks(origin);
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

function ticks=decimal_ticks(value)
if isa(value,'uint64'), ticks=value; return; end
if ischar(value)
    assert(isrow(value)&&~isempty(value),'libekf:Time','One decimal timestamp per string required.');
    value={value};
elseif isstring(value)
    assert(~any(ismissing(value(:))),'libekf:Time','Missing timestamp.');
    value=cellstr(value);
end
assert(iscell(value)&&isvector(value)&&~isempty(value),'libekf:Time', ...
    'Expected uint64 ticks or decimal timestamp strings.');
ticks=zeros(size(value),'uint64'); ten=uint64(10); maximum=intmax('uint64');
for k=1:numel(value)
    text=value{k};
    assert(ischar(text)&&isrow(text)&&~isempty(text)&&all(text>='0'&text<='9'), ...
        'libekf:Time','Timestamp must contain decimal digits only.');
    number=uint64(0);
    for j=1:numel(text)
        digit=uint64(text(j)-'0');
        assert(number<=idivide(maximum-digit,ten,'floor'),'libekf:TimePrecision', ...
            'Decimal timestamp exceeds uint64.');
        number=number*ten+digit;
    end
    ticks(k)=number;
end
end

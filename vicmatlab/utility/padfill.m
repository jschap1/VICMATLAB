% PadFill
%
% Fills missing data in time series (t2, y2). 
% t1 must include t2 (currently)

function yprime = padfill(t1, t2, y2, fillval)

if t1(1) < t2(1) && t1(end) > t2(end)
    yprime = fillval*ones(length(t1), 1);
    [a,b] = ismember(t2, t1);
    [c,d] = ismember(t1, t2);
    yprime(c) = y2;
else
    error('t1 does not include t2')
end

% figure
% plot(t1, yprime)

return
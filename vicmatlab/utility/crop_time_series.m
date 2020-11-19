% Crop time series
%
% Crops two times series to match

function [y1,y2,t] = crop_time_series(y1, y2, t1, t2)
 
if t1(1) < t2(1) & t1(end) < t2(end)
    i1 = find(t1 == t2(1));
    i2 = find(t2 == t1(end));
    t = t1(i1:end);
    y1 = y1(i1:end);
    y2 = y2(1:i2);
    
elseif t1(1) > t2(1) & t1(end) > t2(end)
    i1 = find(t2 == t1(1));
    i2 = find(t1 == t2(end));
    t = t1(1:i2);
    y1 = y1(1:i2);
    y2 = y2(i1:end);
    
elseif t1(end) < t2(1)
    error('Time series do not overlap')
    
elseif t1(1) < t2(1) & t1(end) > t2(end)
    % t2 is contained within t1
    t = t2;
    i1 = find(t1 == t2(1));
    i2 = find(t1 == t2(end));
    y1 = y1(i1:i2);
    
elseif t1(1) > t2(1) & t1(end) < t2(end)
    % t1 is contained within t2
    t = t1;
    i1 = find(t2 == t1(1));
    i2 = find(t2 == t1(end));
    y2 = y2(i1:i2);
    
end


firsttime = max(t1(1), t2(1));
lasttime = min(t1(end), t2(end));




return
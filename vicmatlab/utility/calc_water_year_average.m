% Calculate water year average
%
% Calculate the average value of a time series variable on each day of the
% water year. Drops the last day of a leap year.

function [var_avg, dowy] = calc_water_year_average(t, var)

[t, var] = check_water_year(t, var);

dowy = 1:365;

% Choosing an arbitrary year with 365 days
dowy = datetime(2000, 10 ,1):datetime(2001,9,30);

var_avg = NaN(365,1);
for i=1:365
    ind = month(t) == month(dowy(i)) & day(t) == day(dowy(i));
    var_avg(i) = mean(var(ind));
end

return
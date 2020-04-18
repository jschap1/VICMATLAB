% Annual average
%
% Calculates the annual average values (on a water-year basis) of the input
% variable, which should be input as a time series

function [t_out, var_out] = annual_average(timevector, var, category)

%% Input checking

needs_clipping = 0;

% Check to make sure the inputs are on a water-year basis
if month(timevector(1)) ~= 10 || day(timevector(1)) ~= 1
    disp('Must be on a water year basis')
    disp('First day is not Oct. 1. Clipping.')
    needs_clipping = 1;
end
if month(timevector(end)) ~= 9 || day(timevector(end)) ~= 30
    disp('Must be on a water year basis')
    disp('Last day is not Sept. 30. Clipping.')
    needs_clipping = 1;
end

% If the inputs are not on a water-year basis, do something to ensure they
% are. Try clipping off the ends.

if needs_clipping

    t1 = find(month(timevector) == 10 & day(timevector)==1);
    t2 = find(month(timevector) == 9 & day(timevector)==30);
    t_first = t1(1);
    t_last = t2(end);

    timevector = timevector(t_first:t_last);
    var = var(t_first:t_last);

end

%% Calculate the average

t_out = unique(year(timevector));
t_out = t_out(2:end);
nyears = length(t_out);

t1 = find(month(timevector) == 10 & day(timevector)==1);
t2 = find(month(timevector) == 9 & day(timevector)==30);

var_out = zeros(nyears, 1);
for t=1:nyears
    
    switch category
        case 'flux'
            var_out(t) = sum(var(t1(t):t2(t)));
        case 'store'
            var_out(t) = mean(var(t1(t):t2(t)));
    end      
    
end

return
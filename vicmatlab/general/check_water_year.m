function [tv_out, var_out] = check_water_year(timevector, var_in)

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

if needs_clipping

    t1 = find(month(timevector) == 10 & day(timevector)==1);
    t2 = find(month(timevector) == 9 & day(timevector)==30);
    t_first = t1(1);
    t_last = t2(end);

    tv_out = timevector(t_first:t_last);
    var_out = var_in(t_first:t_last);

else
    % no clipping needed
    tv_out = timevector;
    var_out = var_in;
end

return
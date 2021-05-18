% Utility function for aggregating daily to monthly data
%
% t1 must be a column vector (nt by 1)

function [t_m, q_m] = hourly_to_daily(t1, q1, method)

YMD = unique([year(t1), month(t1), day(t1)], 'rows');
nt_m = size(YMD, 1);
q_m = zeros(nt_m, 1);
t_m = datetime(YMD(:,1), YMD(:,2), YMD(:,3));

for t=1:nt_m
    ind = year(t1) == YMD(t,1) & month(t1) == YMD(t,2) & day(t1) == YMD(t,3);
    try
        if strcmp(method, 'mean')
            q_m(t) = mean(q1(ind));
        elseif strcmp(method, 'sum')
            q_m(t) = sum(q1(ind));
        end
    catch
        continue;
    end
end

% figure
% plot(t, q_m)

return
% Utility function for aggregating daily to monthly data
%
% t1 must be a column vector (nt by 1)

function [t_m, q_m] = daily_to_monthly(t1, q1, method)

YM = unique([year(t1), month(t1)], 'rows');
nt_m = size(YM, 1);
q_m = zeros(nt_m, 1);
t_m = datetime(YM(:,1), YM(:,2), 15);

for t=1:nt_m
    ind = year(t1) == YM(t,1) & month(t1) == YM(t,2);
    if strcmp(method, 'mean')
        q_m(t) = mean(q1(ind));
    elseif strcmp(method, 'sum')
        q_m(t) = sum(q1(ind));
    end
end

return
% Utility function for aggregating monthly to annual data
%
% t1 must be a column vector (nt by 1)
% Function works on a water-year basis

function [t2, q2] = monthly_to_annual(t1, q1, method)

water_year = 1;

if water_year
    
    disp('Calculating on a water year basis')
    WY = CY2WY(t1);
    Y = unique(WY);
    nt_y = size(Y, 1);
    q2 = zeros(nt_y, 1);
    t2 = datetime(Y(:,1), 6, 15);
    for t=1:nt_y
        ind = WY == Y(t,1);
        if strcmp(method, 'mean')
            q2(t) = mean(q1(ind));
        elseif strcmp(method, 'sum')
            q2(t) = sum(q1(ind));
        end
    end

else
    disp('Calculating on a calendar year basis')
    Y = unique(year(t1));
    nt_y = size(Y, 1);
    q2 = zeros(nt_y, 1);
    t2 = datetime(Y(:,1), 6, 15);
    for t=1:nt_y
        ind = year(t1) == Y(t,1);
        if strcmp(method, 'mean')
            q2(t) = mean(q1(ind));
        elseif strcmp(method, 'sum')
            q2(t) = sum(q1(ind));
        end
    end    
end    

return
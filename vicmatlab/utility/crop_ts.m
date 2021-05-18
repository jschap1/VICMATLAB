% Crop time series
%
% INPUTS
% y = variable
% t = datetime vector
% d1, d2 = initial and final datetimes

function [y2, t2] = crop_ts(y, t, d1, d2)

i1 = find(t == d1);
i2 = find(t == d2);
y2 = y(i1:i2);
t2 = t(i1:i2);

return
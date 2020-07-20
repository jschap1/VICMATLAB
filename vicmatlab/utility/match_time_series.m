% Match time series
%
% This function takes two or more time series with different start/end
% dates as inputs, and outputs clipped time series with matching start/end dates. 
%
% INPUTS
% ts_struct = structure containing time series (t1, y1), (t2, y2), etc.
% names must be as follows:
% ts_struct.t.t1, .t2, .t3
% ts_struct.y.y1, .y2, .y3
%
% OUTPUTS
% t_out = datetime vector
% y_out = matrix containing the time series variables
%
% WORK IN PROGRESS..

function [t_out, y_out] = match_time_series(start_date, end_date, t, y)

% nseries = length(varargin);

[~, i1] = ismember(start_date, t);
[~, i2] = ismember(end_date, t);

t_out = t(i1:i2);
y_out = y(i1:i2);

return
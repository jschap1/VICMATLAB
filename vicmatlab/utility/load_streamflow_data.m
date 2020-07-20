% Load streamflow data
%
% Loads streamflow data from a text file with YMDQ fields

function [t, q] = load_streamflow_data(filename, q_col, format)

% filename = '/home/jschap/Documents/ESSD/data/naturalized_flow_usgs_09380000.txt'

dat = readmatrix(filename);

if strcmp(format, 'daily')
    t = datetime(dat(:,1), dat(:,2), dat(:,3));
elseif strcmp(format, 'monthly')
    t = datetime(dat(:,1), dat(:,2), 15);
end

q = dat(:,q_col);

return
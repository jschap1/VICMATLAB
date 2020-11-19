% Function for reading streamflow data from USGS
%
%
function [t_usgs, q_usgs] = read_usgs(fname, start_date, end_date, varargin)

% Set default values for optional arguments

numvarargs = length(varargin);
if numvarargs > 3
    error('subset_forcings requires at most one optional input');
end
optargs = {0};
optargs(1:numvarargs) = varargin;
[plotflag] = optargs{:};

%% Process data

dat = readmatrix(fname);
usgs = struct();
usgs.t = datetime(dat(:,1), dat(:,2), 15);
usgs.Q_cms = dat(:,4);

% USGS
i1 = find(usgs.t == start_date);
i2 = find(usgs.t == end_date);
t_usgs = usgs.t(i1:i2);
q_usgs = usgs.Q_cms(i1:i2);

if plotflag
    figure, plot(t_usgs, q_usgs)
end

return
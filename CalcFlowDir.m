function [fd, datapoints] = CalcFlowDir(X, Y)

% INPUTS
% X, a vector of lon coordinates
% Y, a vector of lat coordinates
%
% OUTPUTS
% Flow directions for each grid cell (assumes one coordinate per grid cell)
% Datapoints, a set of lat/lon coordinates that are in the grid cell

numcells = length(X);
u = NaN(numcells,1);
v = NaN(numcells,1);

% assumes [X,Y] are listed in downstream direction
for i=1:numcells-1
    u(i) = X(i+1) - X(i);
    v(i) = Y(i+1) - Y(i);
end

theta = atan2(v,u)*180/pi; % units of degrees

% Assign flow directions using VIC modeling numbering convention
fd = NaN(numcells,1);
fd(theta >= 112.5 & theta <157.5) = 8;
fd(theta >= 157.5 | theta <= -157.5) = 7;
fd(theta <= -112.5 & theta > -157.5) = 6;
fd(theta <= -67.5 & theta > -112.5) = 5;
fd(theta <= -22.5 & theta > -67.5) = 4;
fd(theta <= 22.5 & theta > -22.5) = 3;
fd(theta < 67.5 & theta >= 22.5) = 2;
fd(theta >= 67.5 & theta <112.5) = 1;

% Remove NaNs
nan_inds = isnan(fd);
fd(nan_inds) = [];
datapoints = [X(~nan_inds); Y(~nan_inds)]';

return

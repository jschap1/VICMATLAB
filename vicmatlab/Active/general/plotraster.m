% Plot Raster
%
% Revised 4/4/2020 JRS
%
% Optional arguments (and defaults):
% xtext = 'Lon'
% ytext = 'Lat'
% inan = 1
% fontsize = 18
% target_axis = gca
%
% Usage
% plotraster(x, y, r, title, xlab, ylab, inan, fontsize, target_axis)
% plotraster(x,y,r,title, xlab, ylab, inan)
% plotraster(x,y,r,title)

function im = plotraster(x, y, r, titletext, varargin)

numvarargs = length(varargin);
if numvarargs > 5
    error('The max number of optional arguments is 5')
end

optargs = {'Lon', 'Lat', 1, 18, gca};
optargs(1:numvarargs) = varargin;
[xtext, ytext, inan, fontsize, target_axis] = optargs{:};

if length(x) > 2
    xx = zeros(2,1);
    xx(1) = min(x);
    xx(2) = max(x);
    x = xx;
end

if length(y) > 2
    yy = zeros(2,1);
    yy(1) = min(y);
    yy(2) = max(y);
    y = yy;
end

% inan = 1;
if inan
    im = imagescnan(x,y,r, 'Parent', target_axis);
else
    im = imagesc(x, y, r, 'Parent', target_axis);
end

title(target_axis, titletext)
xlabel(target_axis, xtext)
ylabel(target_axis, ytext)

% cmap = colormap;
% cmap(1,:) = [1,1,1];
% colormap(cmap);
colorbar(target_axis);
% 

set(target_axis, 'ydir', 'normal')
set(target_axis, 'fontsize', fontsize)

return
% Upscale raster
%
%
% INPUTS
% r = input raster
% lons = list of input longitudes
% lats = list of input latitudes
% newres, oldres = new and old resolutions
% method = linear, nearest, etc.

function [newr, newlons, newlats] = upscale_raster(r, lons, lats, newres, oldres, method)

remove_nan = 1; % flag for removing nan values from output

[m1, n1] = size(r); % original dimensions

f = newres/oldres; % upscaling (or downscaling) factor

% Upscale (or downscale) (method from Walter Roberson in MATLAB Answers)
minlat = min(lats);
maxlat = max(lats);
minlon = min(lons);
maxlon = max(lons);

oldlats = linspace(minlat, maxlat, m1);
oldlons = linspace(minlon, maxlon, n1)';

% find which method works best and use that 

newlats1 = (minlat - oldres/2):newres:(maxlat + oldres/2); % include outside edge
newlons1 = (minlon - oldres/2):newres:(maxlon + oldres/2);
newlons1 = newlons1';

newlats2 = (minlat - oldres/2) + newres:newres:(maxlat + oldres/2) - newres; % exclude outside edge
newlons2 = (minlon - oldres/2) + newres:newres:(maxlon + oldres/2) - newres; 
newlons2 = newlons2';

a1 = abs(newlons1(1) - minlon) + abs(newlons1(end) - maxlon);
a2 = abs(newlons2(1) - minlon) + abs(newlons2(end) - maxlon);
ii = find(min([a1,a2]));

b1 = abs(newlats1(1) - minlat);
b2 = abs(newlats1(end) - maxlat);
jj = find(min([b1,b2]));

if ii==1
    newlons = newlons1;
elseif ii==2
    newlons = newlons2;
end

if jj==1
    newlats = newlats1;
elseif jj==2
    newlats = newlats2;
end

newr = interp2(oldlats, oldlons, r', newlats, newlons, method);
newr = newr';

if remove_nan
    newr(isnan(newr))=0;
end

return
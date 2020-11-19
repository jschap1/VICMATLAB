function inflows = GetInflowPoints(river, chainage, labels, facc, flowthres)

% Calculates inflow points for creating BCI file.
% Part of the LISFLOOD input file preparation workflow
%
% I would use Kostas' Python function to do this, but it is thrown off by
% my rasters. Perhaps it would work if I created all the inputs in GRASS.
%
% I am holding off on making this code for now. I will TRY TO GET KOSTAS
% CODE TO WORK, INSTEAD.
%
% Example usage:
%
% cd /Users/jschapMac/Desktop/Tuolumne/TuoSub/GIS
% addpath TempGIS
% 
% [facc, R] = geotiffread('facc.tif');
% [river, r] = arcgridread('river.asc');
% labels = geotiffread('labels.tif');
% chainage = geotiffread('chainage.tif');
% 
% facc = double(facc);
% facc(facc<0) = NaN;
% 
% flowthres = 100;
% 
% inflows = GetInflowPoints(river, chainage, labels, facc, flowthres)

% upstream inflows are where chainage = 0
inflows = find(chainage == 0);







return
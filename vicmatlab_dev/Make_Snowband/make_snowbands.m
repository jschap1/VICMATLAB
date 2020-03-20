% Make snow band/elevation band file
%
% For now, uses the soil parameter file to make a snowband file with just
% one snowband.
%
% Todo:
% Add functionality to create an actual snowband file that uses high(er)
% resolution elevation data.
%
% INPUTS
% soils = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/soils/global_soils_1_16.txt');
% n = 5; number of elevation bands
% savename = 'global_snowbands_1_16.txt';

function snowbands = make_snowbands(soils, n, savename)

lat = soils(:,3);
lon = soils(:,4);

nlat = length(unique(lat));
nlon = length(unique(lon));
ncells = length(lat);

rasterSize = [nlat, nlon];
latlim = [min(lat), max(lat)];
lonlim = [min(lon), max(lon)];

R = georefcells(latlim,lonlim,rasterSize);

snowbands = zeros(ncells, 3*n+1);
snowbands(:,1) = soils(:,2); % grid cell ID

snowbands(:,2) = 1;
snowbands(:,n+2) = soils(:,18); % elevation
snowbands(:,2*n+2) = 1;

% save the soil parameter file
dlmwrite(savename, snowbands)

% save a geotiff version (for example)
% area_fraction = xyz2grid(lon, lat, snowbands(:, 2));
% geotiffwrite('nveg_map.tif', flipud(area_fraction), R)

return
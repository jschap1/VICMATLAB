% Plot elevation bands
%
% INPUTS
% elevband, elevation band file
% soils, soil parameter file (just need the first 5 columns)
% layer = which elevation band to plot

% elevband = load('/Volumes/HD3/SWOTDA/Data/IRB/VIC/IRB_elevation_bands/mysnowbands.txt');
% soils = load('/Volumes/HD3/SWOTDA/Data/IRB/VIC/IRB_elevation_bands/soils.SB');

function plot_elevband(elevband, soils, layer)

lat = soils(:,3);
lon = soils(:,4);

nbands = (size(elevband,2)-1)/3;
band = xyz2grid(lon, lat, elevband(:,nbands+1+layer));

imagescnan(band)

return


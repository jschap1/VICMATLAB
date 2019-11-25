% Plot elevation bands
%
% INPUTS
% elevband, elevation band file
% soils, soil parameter file (just need the first 5 columns)
% layer = which elevation band to plot

% elevband = load('/Volumes/HD3/SWOTDA/Data/IRB/VIC/IRB_elevation_bands/mysnowbands.txt');
% soils = load('/Volumes/HD3/SWOTDA/Data/IRB/VIC/IRB_elevation_bands/soils.SB');

% soils = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/soils/soils_3L_MERIT_latest.txt');
% elevband = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/snowbands_MERIT_latest.txt');

function plot_elevband(elevband, soils, layer, savename)

lat = soils(:,3);
lon = soils(:,4);

ncells = size(soils,1);
gridID = elevband(:,1);
eblat = zeros(ncells,1);
eblon = zeros(ncells,1);
for k=1:ncells
    eblat(k) = lat(ismember(gridID, soils(k,2)));
end

nbands = (size(elevband,2)-1)/3;
band = xyz2grid(lon, lat, elevband(:,nbands+1+layer));

% imagesc(band)

R = makerefmat(min(lon), min(lat), 0.0625, 0.0625);
geotiffwrite(savename, band, R);

return


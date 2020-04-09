% Plot dsmax
%
% This code should be generalized to crop any soil parameter to a basin
% mask and plot it

% Plot dsmax

% basin_mask = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/GIS/dem.tif';
basin_mask = '/Volumes/HD4/SWOTDA/Data/UpperMiss/umrb_mask.tif';
outname = '/Volumes/HD4/SWOTDA/Data/UpperMiss/out_wy1993/figures/Livneh_Comparison/l15_dsmax.mat';

% L2015
[ws, Rvg, lon1, lat1] = geotiffread2('/Volumes/HD3/VICParametersCONUS/soil_data/Ws.tif');
% [ds, Rvg, lon1, lat1] = geotiffread2('/Volumes/HD3/VICParametersCONUS/soil_data/Ds.tif');
% [dsmax, Rvg, lon1, lat1] = geotiffread2('/Volumes/HD3/VICParametersCONUS/soil_data/Dsmax.tif');

% VICGlobal
% [dsmax, Rvg, lon1, lat1] = geotiffread2('/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Figures/tifs/dsmax.tif');

% [mask1, Rmask, lon_mask, lat_mask] = geotiffread2(basin_mask);

% figure
% plotraster(lon_mask, lat_mask, mask1, 'Elevation (m)', '', '')

figure
plotraster(lon1, lat1, ws, 'Ws', '', '')

% Mask out the part of dsmax that we want to keep

[mask1, R1] = geotiffread(basin_mask);
mask1 = flipud(mask1);
R1mat = georefobj2mat(R1, 'LL');
[masklon1, masklat1] = pixcenters(R1mat, size(mask1));
[masklon, masklat, ~] = grid2xyz(masklon1', masklat1', mask1);
ncells = length(masklon);

figure
plotraster(masklon1, masklat1, mask1, 'Elevation (m)', '', '')

lat_ind = zeros(ncells,1);
lon_ind = zeros(ncells,1);
for k=1:ncells
        [~, lat_ind(k)] = min(abs(masklat(k) - lat1));
        [~, lon_ind(k)] = min(abs(masklon(k) - lon1));
end    

B = NaN(ncells, 1); % the subsetted data
for k=1:ncells      
    B(k) = ws(lat_ind(k),lon_ind(k)); % swapping the indices makes it work...   
end

ws_sub = xyz2grid(masklon, masklat, B);
figure, plotraster(masklon, masklat, ws_sub, 'ws', '', '')

save(outname, 'dsmax_sub', 'masklon', 'masklat');
% save(fullfile(outdir, 'vg_dsmax.mat'), 'dsmax_sub', 'masklon', 'masklat');

livneh = load(fullfile(outdir, 'l15_dsmax.mat'));
vicglobal = load(fullfile(outdir, 'vg_dsmax.mat'));

figure
subplot(2,1,1)
plotraster(livneh.masklon, livneh.masklat, livneh.dsmax_sub, 'Dsmax (L15) (mm/day)', 'Lon', 'Lat')
subplot(2,1,2)
plotraster(vicglobal.masklon, vicglobal.masklat, vicglobal.dsmax_sub, 'Dsmax (VG) (mm/day)', 'Lon', 'Lat')

%% Checking out some other soil parameters

ds_livneh = geotiffread2('/Volumes/HD3/VICParametersCONUS/soil_data/Ds.tif');
ws_livneh = geotiffread2('/Volumes/HD3/VICParametersCONUS/soil_data/Ws.tif');

ds_vg = geotiffread2('/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Figures/tifs/ds.tif');
ws_vg = geotiffread2('/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Figures/tifs/ws.tif');

%% Changing soil parameter values to get a better match

soils_tuo = load('/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_7/Classic_VG_modified/upptuo_soils.txt');

soils_tuo(:,6) = 0.01; % ds
soils_tuo(:,7) = 8.88; % dsmax
soils_tuo(:,8) = 0.75; % ws

write_soils(5, soils_tuo, '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_7/Classic_VG_modified/upptuo_soils.txt', '3l')




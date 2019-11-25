% Make Soil File
%
% Script that uses HWSD data to make a soil parameter file for VIC 4 or 
% VIC 5 (classic)
%
% Dependencies
% addpath(genpath('/Users/jschap/Documents/MATLAB/from_file_exchange'))
% myregrid.m, soil_classification.m, pedotransfer_table.txt
%
% Updated 3/12/2019 JRS - added fs_active to base set of soil parameters
% since it's required to run VIC, even if not using frozen soils
%
% Updated 4/1/2019 - changed to MERIT DEM, three soil layers
%
% Updated 4/4/2019 - complete global coverage (formerly just covered HWSD
% region) (Well, actually from -60 S to 90 N)
%
% Updated 4/10/2019 - filled in missing HWSD soil data using PDE-based
% inpainting method
%
% Updated 4/15/2019 - put lat/lons on a grid consistent with the downscaled
% MERRA-2 data
%
% Updated 9/13/2019 - took out extraneous code, made it easier to run in a
% consistent manner

%% INPUTS

addpath(genpath('/Users/jschap/Documents/MATLAB/from_file_exchange'))
addpath(genpath('/Volumes/HD3/SWOTDA/Codes'));

res = 1/16;
grid_decimal = 5; % number of decimal points for forcing file names

elevname = '/Volumes/HD2/MERIT/DEM/Merged_1_16/merged_merit_dem_1_16.tif';
maskname = '/Volumes/HD2/MERIT/DEM/Merged_1_16/merit_mask_1_16.tif';
hwsd_dir = '/Volumes/HD3/HWSD/HWSD_1247/data';
pedo_name = '/Users/jschap/Documents/Codes/VICMATLAB/pedotransfer_table.txt';
quartz_table = '/Users/jschap/Documents/Codes/VICMATLAB/peters-lidard_quartz_values';
prec_names = dir(fullfile('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec', '*tif'));
temp_names = dir(fullfile('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg', '*tif'));
prec_dir = '/Volumes/HD3/WorldClim/wc2.0_2.5m_prec';
tavg_dir = '/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg';

outdir = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_4/';
savename = 'soils_3L_MERIT.txt';

%% Part 1 - load elevation and land fraction data

% load elevation data from MERIT
[elev, Rdem] = geotiffread(elevname);

% target lat lon for regridding (locations at cell centers)
Rdem_mat = georefobj2mat(Rdem, 'LL');
[target_lon, target_lat] = pixcenters(Rdem_mat, size(elev));

elev = flipud(elev);
elev(elev==-9999) = NaN;

% load land cover mask
[landmask, ~] = geotiffread(maskname);
landmask = flipud(landmask);
landmask = logical(landmask);
landcells = find(landmask);
ncells = length(landcells);

% plot MERIT land mask
figure, subplot(2,1,1)
imagesc(target_lon, target_lat, landmask)
xlabel('Longitude')
ylabel('Latitude')
title('MERIT coverage mask')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

% plot MERIT elevation map
subplot(2,1,2)
imagesc(target_lon, target_lat, elev), colorbar
xlabel('Longitude')
ylabel('Latitude')
title('MERIT elevations (m)')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

% compute slope
[~, slope, ~, ~] = gradientm(elev, Rdem);

% find missing values of slope
missing_slope = isnan(slope) & landmask==1;
figure
plotraster(target_lon, target_lat, missing_slope, 'Missing slope locations', 'Lon', 'Lat')

% fill missing values of slope
slope1 = slope;
slope1 = fillmissing(slope, 'nearest');
slope1(~landmask) = NaN;
slope = slope1;

% Convert slope from degrees to radians
slope = pi*slope/180;

% plot slope
figure
imagesc([target_lon(1), target_lon(end)], [target_lat(1), target_lat(end)], slope)
xlabel('Longitude')
ylabel('Latitude')
title('Slope')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

% calculate lonlat for soil parameter file
xyz_upscaled = raster2xyz(target_lon', target_lat', ones(size(elev)));
lonlat = xyz_upscaled(landcells,1:2);

% lonlat2 = [];
% [lonlat2(:,1), lonlat2(:,2), lonlat2(:,3)] = grid2xyz(target_lon', target_lat', ones(size(elev)));
% lonlat2 = lonlat2(landcells,1:2);
% isequal(lonlat2, lonlat);
% these methods are equivalent

% Save slope raster
geotiffwrite(fullfile(outdir, 'slope.tif'), flipud(slope), Rdem)
disp(['Saved slope raster to ' fullfile(outdir, 'slope.tif')]);

%% Part 2 - load HWSD data

% metadata = ncinfo('../HWSD/HWSD_1247/data/T_SAND.nc4');

% lat, lon for HWSD 0.05 degree grid
hwsd_lat = ncread(fullfile(hwsd_dir, 'T_SAND.nc4'), 'lat');
hwsd_lon = ncread(fullfile(hwsd_dir, 'T_SAND.nc4'), 'lon');

% soil texture
t_sand = ncread(fullfile(hwsd_dir, 'T_SAND.nc4'), 'T_SAND'); 
s_sand = ncread(fullfile(hwsd_dir, 'S_SAND.nc4'), 'S_SAND'); 
t_clay = ncread(fullfile(hwsd_dir, 'T_CLAY.nc4'), 'T_CLAY');
s_clay = ncread(fullfile(hwsd_dir, 'S_CLAY.nc4'), 'S_CLAY');

% organic fraction
t_org = ncread(fullfile(hwsd_dir, 'T_OC.nc4'), 'T_OC');
s_org = ncread(fullfile(hwsd_dir, 'S_OC.nc4'), 'S_OC');

% bulk density
t_bd = ncread(fullfile(hwsd_dir, 'T_BULK_DEN.nc4'), 'T_BULK_DEN');
s_bd = ncread(fullfile(hwsd_dir, 'S_BULK_DEN.nc4'), 'S_BULK_DEN');

% regrid to desired resolution
[hwsd_lons, hwsd_lats] = ndgrid(hwsd_lon, hwsd_lat);
[rglons, rglats] = ndgrid(target_lon, target_lat);

method1 = 'linear';
t_sand_rg = myregrid(hwsd_lons, hwsd_lats, rglons, rglats, t_sand, method1)';
t_clay_rg = myregrid(hwsd_lons, hwsd_lats, rglons, rglats, t_clay, method1)';
s_sand_rg = myregrid(hwsd_lons, hwsd_lats, rglons, rglats, s_sand, method1)';
s_clay_rg = myregrid(hwsd_lons, hwsd_lats, rglons, rglats, s_clay, method1)';
t_bd_rg = myregrid(hwsd_lons, hwsd_lats, rglons, rglats, t_bd, method1)';
s_bd_rg = myregrid(hwsd_lons, hwsd_lats, rglons, rglats, s_bd, method1)';
t_org_rg = myregrid(hwsd_lons, hwsd_lats, rglons, rglats, t_org, method1)';
s_org_rg = myregrid(hwsd_lons, hwsd_lats, rglons, rglats, s_org, method1)';

%% Fill in missing data

% takes 4.3 minutes per function call
censor1 = 1;
t_sand_filled = inpaint_hwsd(t_sand_rg, landmask, censor1); 
t_clay_filled = inpaint_hwsd(t_clay_rg, landmask, censor1);
s_sand_filled = inpaint_hwsd(s_sand_rg, landmask, censor1);
s_clay_filled = inpaint_hwsd(s_clay_rg, landmask, censor1);
t_bd_filled = inpaint_hwsd(t_bd_rg, landmask, censor1);
s_bd_filled = inpaint_hwsd(s_bd_rg, landmask, censor1);
t_org_filled = inpaint_hwsd(t_org_rg, landmask, censor1);
s_org_filled = inpaint_hwsd(s_org_rg, landmask, censor1);

% Checking the effect of censoring the data (negligible)
% addpath('/Users/jschap/Documents/MATLAB')
% t_sand_filled = censor(t_sand_filled, 0, 100);
% t_clay_filled = censor(t_clay_filled, 0, 100);
% s_sand_filled = censor(s_sand_filled, 0, 100);
% s_clay_filled = censor(s_clay_filled, 0, 100);
% 
% soiltexture.censored.topsoil.sand = t_sand_filled;
% soiltexture.censored.subsoil.sand = s_sand_filled;
% soiltexture.censored.topsoil.clay = t_clay_filled;
% soiltexture.censored.subsoil.clay = s_clay_filled;
% 
% figure
% subplot(2,2,1)
% plotraster(target_lon, target_lat, soiltexture.censored.topsoil.sand, 'Sand (topsoil)','','')
% subplot(2,2,2)
% plotraster(target_lon, target_lat, soiltexture.censored.topsoil.clay, 'Clay (topsoil)','','')
% subplot(2,2,3)
% plotraster(target_lon, target_lat, soiltexture.censored.subsoil.sand, 'Sand (subsoil)','','')
% subplot(2,2,4)
% plotraster(target_lon, target_lat, soiltexture.censored.subsoil.clay, 'Clay (subsoil)','','')

% Save the inpainted HWSD data
geotiffwrite(fullfile(outdir, 'percent_sand_topsoil_inpainted.tif'), flipud(t_sand_filled), Rdem)
geotiffwrite(fullfile(outdir, 'percent_clay_topsoil_inpainted.tif'), flipud(t_clay_filled), Rdem)
geotiffwrite(fullfile(outdir, 'percent_sand_subsoil_inpainted.tif'), flipud(s_sand_filled), Rdem)
geotiffwrite(fullfile(outdir, 'percent_clay_subsoil_inpainted.tif'), flipud(s_clay_filled), Rdem)
geotiffwrite(fullfile(outdir, 'bulk_density_topsoil_inpainted.tif'), flipud(t_bd_filled), Rdem)
geotiffwrite(fullfile(outdir, 'bulk_density_subsoil_inpainted.tif'), flipud(s_bd_filled), Rdem)
disp(['Saved inpainted, interpolated HWSD data to ' outdir]);

% Plot bulk density before and after inpainting

missing_bulk_dens = isnan(t_bd_rg) & landmask==1;
t_bd_rd_copy = t_bd_rg;
t_bd_rd_copy(missing_bulk_dens) = 0;
t_bd_diff = t_bd_filled - t_bd_rd_copy;

figure
subplot(3,1,1)
plotraster(target_lon, target_lat, 1000*t_bd_rg, 'Topsoil bulk density (kg/m^3) (original)', 'Lon', 'Lat')
subplot(3,1,2)
plotraster(target_lon, target_lat, 1000*t_bd_filled, 'Topsoil bulk density (kg/m^3) (filled)', 'Lon', 'Lat')
subplot(3,1,3)
plotraster(target_lon, target_lat, 1000*t_bd_diff, 'Difference (filled - original)', 'Lon', 'Lat')

%% Part 3

% Classify soils using USGS soil triangle (this step takes a few minutes)
t_percent_sand = t_sand_filled(landcells);
t_percent_clay = t_clay_filled(landcells);
[~, ~, t_sc_int, ~] = soil_classification(t_percent_sand/100, t_percent_clay/100, 'usda', 0);

s_percent_sand = s_sand_filled(landcells);
s_percent_clay = s_clay_filled(landcells);
[~, ~, s_sc_int, ~] = soil_classification(s_percent_sand/100, s_percent_clay/100, 'usda', 0);
% sc_int is the USGS soil class, which is used in the lookup table
 
figure 
subplot(2,2,1), hist(t_percent_sand), title('Percent Sand (topsoil)')
subplot(2,2,2), hist(s_percent_sand), title('Percent Sand (subsoil)')
subplot(2,2,3), hist(t_sc_int), title('USDA class (topsoil)')
subplot(2,2,4), hist(s_percent_sand), title('Percent Sand (subsoil)')

topsoil_soil_class = NaN(size(t_sand_rg));
topsoil_soil_class(landcells) = t_sc_int;

subsoil_soil_class = NaN(size(t_sand_rg));
subsoil_soil_class(landcells) = s_sc_int;

figure, 
subplot(2,1,1)
plotraster(target_lon, target_lat, topsoil_soil_class, 'Topsoil','','')
subplot(2,1,2)
plotraster(target_lon, target_lat, subsoil_soil_class, 'Subsoil','','')

% Rsoilclass = makerefmat(min(lonlat(:,1)), min(lonlat(:,2)), res, res);
geotiffwrite(fullfile(outdir, 'topsoil_soil_class.tif'), flipud(topsoil_soil_class), Rdem)
geotiffwrite(fullfile(outdir, 'subsoil_soil_class.tif'), flipud(subsoil_soil_class), Rdem)
disp(['Saved USDA soil class data to ' outdir]);

% % ------------------------------------------------------
% % Save soil textures as geotiffs (how to do this?)
% 
% M1 = single(landmask);
% M2 = single(landmask);
% 
% % this assignment does not work for some unknown reason...
% M1(landcells) = t_sc_int;
% M2(landcells) = s_sc_int;
% 
% % landmask(landcells)
% % 
% % M1 = xyz2grid(lonlat(:,1), lonlat(:,2), t_sc_int);
% % M2 = xyz2grid(lonlat(:,1), lonlat(:,2), s_sc_int);
% 
% % append rows to the soil texture matrix to account for the missing
% % southern latitudes
% % nmissingrows = size(t_sand_rg, 1) - size(M1, 1);
% % M1 = [M1; NaN(nmissingrows, size(M1, 2))];
% % M2 = [M2; NaN(nmissingrows, size(M2, 2))];
% 
% R = makerefmat(target_lon(1), target_lat(1), res, res);
% 
% geotiffwrite('topsoil_texture.tif', M1, R)
% geotiffwrite('subsoil_texture.tif', M2, R)
% 
% figure
% imagesc(target_lon, target_lat, M2)
% xlabel('Longitude')
% ylabel('Latitude')
% title('USDA Soil Textures (subsoil)')
% set(gca, 'ydir', 'normal')
% set(gca, 'fontsize', 18)
% 
% % Not sure how to get legend to appear as desired:
% % legend('sand','loamy-sand','sandy-loam','silty-loam',...
% %     'silt','loam','sandy-clay-loam','silty-clay-loam',...
% %     'clay','sandy-clay','silty-clay', 'clay');

% ------------------------------------------------------

% Import pedotransfer table
T = dlmread(pedo_name, '\t', 1, 0);

% Estimate hydraulic parameters from soil texture

ksat1 = zeros(ncells,1);
expt1 = zeros(ncells,1);
% bulk_dens1 = zeros(ncells,1);
wcr_fract1 = zeros(ncells,1);
wpwp_fract1 = zeros(ncells,1);
porosity1 = zeros(ncells,1);

ksat2 = zeros(ncells,1);
expt2 = zeros(ncells,1);
% bulk_dens2 = zeros(ncells,1);
wcr_fract2 = zeros(ncells,1);
wpwp_fract2 = zeros(ncells,1);
porosity2 = zeros(ncells,1);

for sc=1:12
    ksat1(t_sc_int == sc) = 240*T(sc, 6); % mm/day
    expt1(t_sc_int == sc) = 3 + 2*T(sc, 7);
%     bulk_dens1(t_sc_int == sc) = 1000*T(sc, 2); % kg/m^3
    wcr_fract1(t_sc_int == sc) = T(sc, 3)/T(sc, 5); % fraction of maximum moisture (double check this calculation, might be incorrect)
    wpwp_fract1(t_sc_int == sc) = T(sc, 4)/T(sc, 5); % fraction of maximum moisture
    porosity1(t_sc_int == sc) = T(sc,5);
    
    ksat2(s_sc_int == sc) = 240*T(sc, 6);
    expt2(s_sc_int == sc) = 3 + 2*T(sc, 7);
%     bulk_dens2(s_sc_int == sc) = 1000*T(sc, 2);
    wcr_fract2(s_sc_int == sc) = T(sc, 3)/T(sc, 5); 
    wpwp_fract2(s_sc_int == sc) = T(sc, 4)/T(sc, 5);
    porosity2(s_sc_int == sc) = T(sc,5);
    disp(sc)
end

% Plot porosity
porosity_map1 = NaN(size(t_sand_rg));
porosity_map1(landcells) = porosity1;

porosity_map2 = NaN(size(t_sand_rg));
porosity_map2(landcells) = porosity2;

figure
plotraster(target_lon, target_lat, porosity_map1, 'Porosity (topsoil)', '','')

% Save porosity map
geotiffwrite(fullfile(outdir, 'porosity_topsoil.tif'), flipud(porosity_map1), Rdem)
geotiffwrite(fullfile(outdir, 'porosity_subsoil.tif'), flipud(porosity_map2), Rdem)

% Note: T is for topsoil, S is for subsoil;
% exponent n is from the Brooks-Corey relationship

%% Part 4

% Average temperature. One option is the ISP GCM output.
% Another is http://www.worldclim.org/formats1

% Downloaded 2.5 minute resolution monthly average temperature data from
% WorldClim, which is a project of Fick and Hijmanns
% http://worldclim.org/version2

% Load monthly average temperature
[temperature1, Rworldclim] = geotiffread(fullfile(tavg_dir, temp_names(1).name));
avgT.m1 = flipud(temperature1);
nan_ind = find(avgT.m1 <= -9999);

for m=1:12
    avgT.(['m' num2str(m)]) = flipud(geotiffread(fullfile(tavg_dir, temp_names(m).name)));
end

sum1 = 0;
for m=1:12
    avgT.(['m' num2str(m)])(nan_ind) = NaN;
    sum1 = sum1 + avgT.(['m' num2str(m)]);
end
avgT.y = sum1/12;

figure, imagesc(avgT.y)
set(gca, 'ydir', 'normal')

% Regrid
Rworldclim_mat = georefobj2mat(Rworldclim, 'LL');
[Tlon, Tlat] = pixcenters(Rworldclim_mat, size(temperature1));
[Tlons, Tlats] = ndgrid(Tlon, Tlat);
avgT_rg = myregrid(Tlons, Tlats, rglons, rglats, avgT.y', 'linear');
avgT_rg = avgT_rg';

ann_T = avgT_rg;
ann_T(~landmask) = NaN;
figure, plotraster(target_lon, target_lat, ann_T, 'Temperatures','','')
set(gca, 'ydir', 'normal')

% find missing values of temperature
missing_T = isnan(ann_T) & landmask==1;
figure, plotraster(target_lon, target_lat, missing_T, 'Missing temperature data','','')
% OK, so fillmissing puts spurious values of temperature around the
% coastlines. Probably better not to have it, or to write a better version
% of the fillmissing function for my specific application.

% fill missing values of temperature
ann_T = fillmissing(ann_T, 'nearest');
ann_T(~landmask) = NaN;

figure, plotraster(target_lon, target_lat, ann_T, 'Interpolated temperatures','','')

% Save the average temperature map
% R = makerefmat(target_lon(1), target_lat(1), res, res);
geotiffwrite(fullfile(outdir, 'worldclim_tavg.tif'), flipud(ann_T), Rdem)
disp(['Saved annual average temperature map to ' fullfile(outdir, 'worldclim_tavg.tif')])

% Also compute July_Tavg

JT_rg = myregrid(Tlons, Tlats, rglons, rglats, avgT.m7', 'linear');
JT_rg = JT_rg';

% fill missing values of temperature
JT_rg = fillmissing(JT_rg, 'nearest');
JT_rg(~landmask) = NaN;

geotiffwrite(fullfile(outdir, 'worldclim_July_tavg.tif'), flipud(JT_rg), Rdem)
disp(['Saved annual average temperature map to ' fullfile(outdir, 'worldclim_July_tavg.tif')])

%% Part 5 

% quartz content - use a pedotransfer function
% this one is from "Estimation of Periodic Water Balance Components...",
% cited on the UW VIC website (2016)

quartz1 = zeros(ncells, 1);
quartz2 = zeros(ncells, 1);
quartz_ptf = readtable(quartz_table);
quartz_vals = quartz_ptf.Quartz;

for sc=1:12
    quartz1(s_sc_int == sc) = quartz_vals(sc);
    quartz2(t_sc_int == sc) = quartz_vals(sc);    
end

%% Part 6

% off_gmt parameter

% tz = shaperead('/Volumes/HD3/TimeZones/ne_10m_time_zones/ne_10m_time_zones.shp');
% convert to raster using gdal/ogr

% gdal_rasterize -te -179.9750 -89.9750 179.9750 89.9750 -ts 5760 2880 -a zone /Volumes/HD3/TimeZones/ne_10m_time_zones/ne_10m_time_zones.shp time_zones.tif
% gdal_rasterize -te -179.9750 -59.9750 179.9750 89.9750 -ts 5760 2320 -a zone /Volumes/HD3/TimeZones/ne_10m_time_zones/ne_10m_time_zones.shp time_zones_above_60.tif

tz = geotiffread('/Volumes/HD3/TimeZones/time_zones_above_60.tif');
tz = flipud(tz);
figure, imagesc(target_lon, target_lat, tz)
set(gca, 'ydir', 'normal')

tz_mask = tz;
tz_mask(~landmask) = NaN;
figure, imagesc(target_lon, target_lat, tz_mask)
set(gca, 'ydir', 'normal')

%% Part 7 

% annual average precipitation

% Again, using the WorldClim data

avgP.m1 = flipud(geotiffread(fullfile(prec_dir, prec_names(1).name)));
nan_ind = find(avgP.m1 <= -9999);

for m=1:12
    avgP.(['m' num2str(m)]) = flipud(geotiffread(fullfile(prec_dir, prec_names(m).name)));
end

sum1 = 0;
for m=1:12
    avgP.(['m' num2str(m)])(nan_ind) = NaN;
    sum1 = sum1 + avgP.(['m' num2str(m)]);
end
avgP.y = sum1;

figure, plotraster(target_lon, target_lat, avgP.y, 'Average PPT', '', '')

% Regrid
avgP_rg = myregrid(Tlons, Tlats, rglons, rglats, double(avgP.y'), 'linear');
avgP_rg = avgP_rg';

ann_P = avgP_rg;
ann_P(~landmask) = NaN;
figure, plotraster(target_lon, target_lat, ann_P, 'Average PPT (regridded)', '', '')

% fill missing values of precipitation
ann_P = fillmissing(ann_P, 'nearest');
ann_P(~landmask) = NaN;
figure, plotraster(target_lon, target_lat, ann_P, 'Average PPT (interpolated)', '', '')

% Save the annual precipitation map
R = makerefmat(target_lon(1), target_lat(1), res, res);
geotiffwrite(fullfile(outdir, 'worldclim_prec.tif'), flipud(ann_P), R)
disp(['Saved annual average precipitation map to ' fullfile(outdir, 'worldclim_prec.tif')])

%% Part 8 

% organic soil properties

%% Part 9 

% residual soil moisture

%% Part 10

% Frozen soil parameters

%% Make the soil parameter file

varnames = get_soil_var_names('3L-no-org-frost-msds');
nvars = length(varnames);

soils = zeros(ncells, nvars); % this is about 1.5 GB, so pretty big

soils(:,1) = 1; % run_cell
soils(:,2) = 1:ncells; % grid_cell
soils(:,3) = lonlat(:,2); % lat
soils(:,4) = lonlat(:,1); % lon
soils(:,5) = 0.2; % b_infilt
soils(:,6) = 0.001; % Ds

% calculate mean ksat for the soil column
ksat_bar = ksat1.*0.1 + ksat1*0.2 + ksat1*0.7;
dsmax = ksat_bar.*slope(landmask); % Dsmax

% quality control on dsmax - may want to do this, unclear if necessary
% set high values equal to a maximum threshold value
% dsmax(dsmax>1e4) = 1e4;
soils(:,7) = dsmax;

soils(:,8) = 0.9; % Ws
soils(:,9) = 2; % c

soils(:,10) = expt1; % expt1
soils(:,11) = expt1; % expt2
soils(:,12) = expt2; % expt3

soils(:,13) = ksat1; % Ksat1
soils(:,14) = ksat1; % Ksat2
soils(:,15) = ksat2; % Ksat3

soils(:,16) = -999; % phi_s1
soils(:,17) = -999; % phi_s2
soils(:,18) = -999; % phi_s3

soils(:,19) = 1000.*wcr_fract1.*0.3.*porosity1; % init_moist1
soils(:,20) = 1000.*wcr_fract1.*0.7.*porosity1; % init_moist2
soils(:,21) = 1000.*wcr_fract2.*0.7.*porosity2; % init_moist3

soils(:,22) = elev(landmask); % elev

% Thickness of each soil layer (m)
soils(:,23) = 0.1;
soils(:,24) = 0.2;
soils(:,25) = 0.7;

soils(:,26) = ann_T(landmask); % avg_T !!!
soils(:,27) = 4; % dp

soils(:,28) = 0.32*expt1 + 4.3; % bubble1
soils(:,29) = 0.32*expt1 + 4.3; % bubble2
soils(:,30) = 0.32*expt2 + 4.3; % bubble3

soils(:,31) = quartz1; % quartz1
soils(:,32) = quartz1; % quartz2
soils(:,33) = quartz2; % quartz3

soils(:,34) = 1000*t_bd_filled(landmask); % bulk_dens1
soils(:,35) = 1000*t_bd_filled(landmask); % bulk_dens2
soils(:,36) = 1000*s_bd_filled(landmask); % bulk_dens3

soils(:,37) = 2685; % soil_dens1
soils(:,38) = 2685; % soil_dens2
soils(:,39) = 2685; % soil_dens3

soils(:,40) = tz(landmask); % off_gmt

soils(:,41) = wcr_fract1; % wcr_fract1
soils(:,42) = wcr_fract1; % wcr_fract2
soils(:,43) = wcr_fract2; % wcr_fract3

soils(:,44) = wpwp_fract1; % wpwp_fract1
soils(:,45) = wpwp_fract1; % wpwp_fract2
soils(:,46) = wpwp_fract2; % wpwp_fract3

soils(:,47) = 0.001; % rough
soils(:,48) = 0.0005; % snow_rough

soils(:,49) = ann_P(landmask); % annual_prec

soils(:,50) = 0; % resid_moist1 (leaving it as 0...)
soils(:,51) = 0; % resid_moist2
soils(:,52) = 0; % resid_moist3

soils(:,53) = 1; % fs_active

soils(:,54) = JT_rg(landmask); % july_Tavg

write_soils(grid_decimal, soils, fullfile(outdir, savename), '3l')

%% Make geotiffs from the soil parameter file

setup = '3L-no-org-frost-msds';
goutdir = fullfile(outdir, 'Figures/tifs');
mkdir(goutdir)
convert_soil_parameters(soils, setup, goutdir, maskname)

%% Plot soil parameters

tifnames = dir(fullfile(goutdir, '*.tif'));
poutdir = fullfile(outdir, 'Figures/pngs');
mkdir(poutdir)

for k=1:length(tifnames)
    
    tmpname = strsplit(tifnames(k).name, '.tif');
    varname1 = tmpname{1};
    
    var1 = flipud(geotiffread(fullfile(goutdir, tifnames(k).name)));
    f = figure('visible', 'off');
    plotraster(target_lon, target_lat, var1, varname1, '','')
    saveas(f, fullfile(poutdir, [varname1 '.png']));
    
end
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

%% INPUTS

addpath(genpath('/Users/jschap/Documents/MATLAB/from_file_exchange'))
addpath(genpath('/Volumes/HD3/SWOTDA/Codes'));

res = 1/16;
grid_decimal = 5; % number of decimal points for forcing file names

% target lat lon for regridding (locations at cell centers) (old)
% target_lon2 = -180+res/2:0.0625:180-res/2;
% target_lat2 = -60+res/2:0.0625:85-res/2;

outdir = '/Volumes/HD3/VICParametersGlobal/Global_1_16/soils/global_3L-2/';
savename = 'soils_3L_MERIT.txt';

%% Part 1 - load elevation and land fraction data

% load elevation data from MERIT
[elev, Rdem] = geotiffread('/Volumes/HD3/VICParametersGlobal/Global_1_16/merged_merit_dem_1_16.tif');

% target lat lon for regridding (locations at cell centers)
Rdem_mat = georefobj2mat(Rdem);
[target_lon, target_lat] = pixcenters(Rdem_mat, size(elev));

elev = flipud(elev);
elev(elev==-9999) = NaN;

% load land cover mask
[landmask, Rmask] = geotiffread('/Volumes/HD3/VICParametersGlobal/Global_1_16/landmask/merit_mask_1_16.tif');
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

subplot(2,1,2)
imagesc(target_lon, target_lat, elev), colorbar
xlabel('Longitude')
ylabel('Latitude')
title('MERIT elevations (m)')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

% compute slope
[~, slope, ~, ~] = gradientm(elev, Rmask);

% find missing values of slope
missing_slope = isnan(slope) & landmask==1;
figure
plotraster(target_lon, target_lat, missing_slope, 'Missing slope locations', 'Lon', 'Lat')

% fill missing values of slope
slope1 = slope;
slope1 = fillmissing(slope, 'nearest');
slope1(~landmask) = NaN;
slope = slope1;

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

% Perform surgery on soil parameter file....
% soils1 = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_2/soils_3L_MERIT.txt');
% soils1(:,3) = lonlat(:,2); % lat
% soils1(:,4) = lonlat(:,1); % lon

% Save slope raster
geotiffwrite(fullfile(outdir, 'slope.tif'), flipud(slope), Rmask)

%% Part 2 - load HWSD data

% metadata = ncinfo('../HWSD/HWSD_1247/data/T_SAND.nc4');

% lat, lon for HWSD 0.05 degree grid
lat = ncread('../HWSD/HWSD_1247/data/T_SAND.nc4', 'lat');
lon = ncread('../HWSD/HWSD_1247/data/T_SAND.nc4', 'lon');

% soil texture
t_sand = ncread('../HWSD/HWSD_1247/data/T_SAND.nc4', 'T_SAND'); 
s_sand = ncread('../HWSD/HWSD_1247/data/S_SAND.nc4', 'S_SAND'); 
t_clay = ncread('../HWSD/HWSD_1247/data/T_CLAY.nc4', 'T_CLAY');
s_clay = ncread('../HWSD/HWSD_1247/data/S_CLAY.nc4', 'S_CLAY');

% organic fraction
t_org = ncread('../HWSD/HWSD_1247/data/T_OC.nc4', 'T_OC');
s_org = ncread('../HWSD/HWSD_1247/data/S_OC.nc4', 'S_OC');

% bulk density
t_bd = ncread('../HWSD/HWSD_1247/data/T_BULK_DEN.nc4', 'T_BULK_DEN');
s_bd = ncread('../HWSD/HWSD_1247/data/S_BULK_DEN.nc4', 'S_BULK_DEN');

% regrid to desired resolution
[lons, lats] = ndgrid(lon, lat);
[rglons, rglats] = ndgrid(target_lon, target_lat);

method1 = 'linear';
t_sand_rg = myregrid(lons, lats, rglons, rglats, t_sand, method1)';
t_clay_rg = myregrid(lons, lats, rglons, rglats, t_clay, method1)';
s_sand_rg = myregrid(lons, lats, rglons, rglats, s_sand, method1)';
s_clay_rg = myregrid(lons, lats, rglons, rglats, s_clay, method1)';
t_bd_rg = myregrid(lons, lats, rglons, rglats, t_bd, method1)';
s_bd_rg = myregrid(lons, lats, rglons, rglats, s_bd, method1)';
t_org_rg = myregrid(lons, lats, rglons, rglats, t_org, method1)';
s_org_rg = myregrid(lons, lats, rglons, rglats, s_org, method1)';

% t_sand_rg_ngb = myregrid(lons, lats, rglons, rglats, t_sand, 'nearest')';
% t_sand_rg_bil = myregrid(lons, lats, rglons, rglats, t_sand, 'linear')';
% t_sand_rg_dif = t_sand_rg_bil - t_sand_rg_ngb;
% figure
% subplot(3,1,1)
% plotraster(lon, lat, t_sand_rg_ngb, 'Nearest Neighbor','','')
% subplot(3,1,2)
% plotraster(lon, lat, t_sand_rg_bil, 'Bilinear','','')
% subplot(3,1,3)
% plotraster(lon, lat, t_sand_rg_dif, 'Difference','','')

%% Fill in missing data

% mask of missing data (bulk density)
missing_bulk_dens = isnan(t_bd_rg_orig) & landmask==1;
figure, imagesc(missing_bulk_dens)

% mask of missing data (soil texture)
missing_sand_percent = isnan(t_sand_rg_orig) & landmask==1;
figure, imagesc(missing_sand_percent)

% fill missing values
% t_sand_rg_filled = t_sand_rg;
% filled_values = fillmissing(t_sand_rg(landmask), 'nearest');
% t_sand_rg_filled(landmask) = filled_values;

% inpaint missing values
% nodatavalue = 9999;
% t_sand_orig = t_sand_rg;
% t_sand_rg(~landmask) = nodatavalue;
% t_sand_rg_filled = inpaint_nans(t_sand_rg);
% t_sand_rg_filled(~landmask) = NaN;

figure, subplot(2,1,1); imagesc(t_sand_orig), colorbar, title('original')
subplot(2,1,2); imagesc(t_sand_rg_filled), colorbar, title('inpainted')

% inpaintn does not allow any NaNs to stay missing. This means that I need
% to fill all the missing values, then mask out the ones I wish to keep as
% NaN

t_sand_rg = inpaint_hwsd(t_sand_rg, landmask); % takes ~4 minutes
t_clay_rg = inpaint_hwsd(t_clay_rg, landmask);
s_sand_rg = inpaint_hwsd(s_sand_rg, landmask);
s_clay_rg = inpaint_hwsd(s_clay_rg, landmask);
t_bd_rg = inpaint_hwsd(t_bd_rg, landmask);
s_bd_rg = inpaint_hwsd(s_bd_rg, landmask);
t_org_rg = inpaint_hwsd(t_org_rg, landmask);
s_org_rg = inpaint_hwsd(s_org_rg, landmask);

% plot before, after, and difference for the filled images
figure

subplot(2,3,1)
imagesc(target_lon, target_lat, t_sand_rg_orig), colorbar
title('Percent sand (topsoil) - original'), xlabel('Lon'), ylabel('Lat')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

subplot(2,3,2)
imagesc(target_lon, target_lat, t_sand_rg), colorbar
title('Percent sand (topsoil) - filled'), xlabel('Lon'), ylabel('Lat')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

subplot(2,3,3)
imagesc(target_lon, target_lat, missing_sand_percent), colorbar
title('Locations of missing values'), xlabel('Lon'), ylabel('Lat')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

subplot(2,3,4)
imagesc(target_lon, target_lat, 1000*t_bd_rg_orig), colorbar
title('Bulk density (topsoil, kg/m^3) - original'), xlabel('Lon'), ylabel('Lat')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

subplot(2,3,5)
imagesc(target_lon, target_lat, 1000*t_bd_rg), colorbar
title('Bulk density (topsoil, kg/m^3) - filled'), xlabel('Lon'), ylabel('Lat')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

subplot(2,3,6)
imagesc(target_lon, target_lat, missing_bulk_dens), colorbar
title('Locations of missing values'), xlabel('Lon'), ylabel('Lat')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)


% testing the "fillmissing" method
% landmask1 = [1,1,1,1,1,0;
%     1,1,1,1,0,1;
%     1,1,1,1,0,1;
%     0,1,1,1,0,0];
% sand1 = [2,2,3,7,9,4;
%     3, NaN, 3,3,3,0;
%     7, NaN, 3,4,4,9;
%     9, NaN, NaN, NaN, NaN, NaN];
% sand_filled = fillmissing(sand1, 'nearest');
% sand_filled(~landmask1) = NaN;
% not optimal because it chooses the nearest neighbor by column, so it does
% not take into account the 2D nature of the data
% but it does basically fill in the missing values, just could be more
% accurate. Probably OK for a dataset that needs to be calibrated anyway.
% what I really want is some sort of bilinear interpolation
%
% Try the inpainting method from John D'Arrico

% sand_filled = inpaint_nans(sand1); % takes ~4 minutes
% sand_filled(~landmask1) = NaN;
% sand_filled(sand_filled>100) = 100;
% sand_filled(sand_filled<0) = 0;

% t_sand_rg = inpaintn(t_sand_rg); % takes about 2 minutes
% t_sand_rg(~landmask) = NaN;
% t_sand_rg(t_sand_rg<0) = 0; % keep the range realistic
% t_sand_rg(t_sand_rg>100) = 100;

% t_clay_rg = inpaintn(t_clay_rg);
% t_clay_rg(~landmask) = NaN;
% t_clay_rg(t_clay_rg<0) = 0;
% t_clay_rg(t_clay_rg>100) = 100;
% 
% s_sand_rg = inpaintn(s_sand_rg);
% s_sand_rg(~landmask) = NaN;
% s_sand_rg(s_sand_rg<0) = 0;
% s_sand_rg(s_sand_rg>100) = 100;
% 
% s_clay_rg = inpaintn(s_clay_rg);
% s_clay_rg(~landmask) = NaN;
% s_clay_rg(s_clay_rg<0) = 0;
% s_clay_rg(s_clay_rg>100) = 100;
% 
% t_bd_rg = inpaintn(t_bd_rg);
% t_bd_rg(~landmask) = NaN;
% t_bd_rg(t_bd_rg<0) = 0;
% 
% s_bd_rg = inpaintn(s_bd_rg);
% s_bd_rg(~landmask) = NaN;
% s_bd_rg(s_bd_rg<0) = 0;
% 
% t_org_rg = inpaintn(t_org_rg);
% t_org_rg(~landmask) = NaN;
% t_org_rg(t_org_rg<0) = 0;
% 
% s_org_rg = inpaintn(s_org_rg);
% s_org_rg(~landmask) = NaN;
% s_org_rg(s_org_rg<0) = 0;

% write out a mask showing where the missing data were filled in
% that is, where there are MERIT data but no HWSD data
% It might vary somewhat depending on the specific HWSD variable

% output masks showing where missing data were filled in
% R = makerefmat(target_lon(1), target_lat(1), res, res);
% geotiffwrite('/Volumes/HD3/VICParametersGlobal/Global_1_16/missing_bd.tif', missing_bulk_dens, R)
% geotiffwrite('/Volumes/HD3/VICParametersGlobal/Global_1_16/missing_texture.tif', missing_sand_percent, R)

% figure, imagesc(target_lon, target_lat, t_sand_rg_nom)
% sum(isnan(t_sand_rg(landmask)))/sum(~isnan(t_sand_rg(landmask))); % about 5% of land cells are missing HWSD data
% t_sand_rg(~landmask) = -9999;
% t_sand_rg_nom = fillmissing(t_sand_rg, 'nearest');
% t_sand_rg_nom(~landmask) = NaN;

% figure, imagesc(target_lon, target_lat, t_sand_rg)

% Make a land mask from HWSD. Use the field "ISSOIL"

% is_soil = ncread('../HWSD/HWSD_1247/data/ISSOIL.nc4', 'ISSOIL');
% 
% % 1 - 5.7e6 pixels -> indicates a soil mapping unit is soil
% % 0 - 4.8e5 pixels -> indicates a soil mapping unit is not soil
% % NaN - 2.0e7 pixels -> ocean
% 
% is_soil_rg = myregrid(lons, lats, rglons, rglats, is_soil, 'linear');
% is_soil_rg = is_soil_rg';

% figure, imagesc(target_lon, target_lat, is_soil_rg)
% set(gca, 'ydir', 'normal')

% Since we've upscaled the soil mask, some pixels are partially soil and
% partially not soil. Need to account for this. Simplest thing to do:
% assume that all partial soil units are soil.

% is_soil_rg(is_soil_rg>0) = 1;

% test = is_soil_rg;
% test(is_soil_rg>0) = 1;
% figure, imagesc(test)

% One more plot for the ESSD paper:

t_bd_rg_filled = inpaint_hwsd(t_bd_rg, landmask);

missing_bulk_dens = isnan(t_bd_rg) & landmask==1;
t_bd_rd_copy = t_bd_rg;
t_bd_rd_copy(missing_bulk_dens) = 0;

t_bd_diff = t_bd_rg_filled - t_bd_rd_copy;

figure
subplot(3,1,1)
plotraster(target_lon, target_lat, t_bd_rg, 'Topsoil bulk density (kg/m^3) (original)', 'Lon', 'Lat')
subplot(3,1,2)
plotraster(target_lon, target_lat, t_bd_rg_filled, 'Topsoil bulk density (kg/m^3) (filled)', 'Lon', 'Lat')
subplot(3,1,3)
plotraster(target_lon, target_lat, t_bd_diff, 'Difference (filled - original)', 'Lon', 'Lat')

%% Part 3

% Classify soils using USGS soil triangle (this step takes a few minutes)
t_percent_sand = t_sand_rg(landcells);
t_percent_clay = t_clay_rg(landcells);
[~, ~, t_sc_int, ~] = soil_classification(t_percent_sand/100, t_percent_clay/100, 'usda', 0);

s_percent_sand = s_sand_rg(landcells);
s_percent_clay = s_clay_rg(landcells);
[~, ~, s_sc_int, ~] = soil_classification(s_percent_sand/100, s_percent_clay/100, 'usda', 0);
% sc_int is the USGS soil class, which is used in the lookup table 

% ------------------------------------------------------
% Save soil textures as geotiffs (how to do this?)

M1 = single(landmask);
M2 = single(landmask);

% this assignment does not work for some unknown reason...
M1(landcells) = t_sc_int;
M2(landcells) = s_sc_int;

% landmask(landcells)
% 
% M1 = xyz2grid(lonlat(:,1), lonlat(:,2), t_sc_int);
% M2 = xyz2grid(lonlat(:,1), lonlat(:,2), s_sc_int);

% append rows to the soil texture matrix to account for the missing
% southern latitudes
% nmissingrows = size(t_sand_rg, 1) - size(M1, 1);
% M1 = [M1; NaN(nmissingrows, size(M1, 2))];
% M2 = [M2; NaN(nmissingrows, size(M2, 2))];

R = makerefmat(target_lon(1), target_lat(1), res, res);

geotiffwrite('topsoil_texture.tif', M1, R)
geotiffwrite('subsoil_texture.tif', M2, R)

figure
imagesc(target_lon, target_lat, M2)
xlabel('Longitude')
ylabel('Latitude')
title('USDA Soil Textures (subsoil)')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

% Not sure how to get legend to appear as desired:
% legend('sand','loamy-sand','sandy-loam','silty-loam',...
%     'silt','loam','sandy-clay-loam','silty-clay-loam',...
%     'clay','sandy-clay','silty-clay', 'clay');

% ------------------------------------------------------

% Import pedotransfer table
T = dlmread('/Users/jschap/Documents/Codes/VICMATLAB/pedotransfer_table.txt', '\t', 1, 0);

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
end

% Note: T is for topsoil, S is for subsoil;
% exponent n is from the Brooks-Corey relationship

%% Part 4

% Average temperature. One option is the ISP GCM output.
% Another is http://www.worldclim.org/formats1

% Downloaded 2.5 minute resolution monthly average temperature data from
% WorldClim, which is a project of Fick and Hijmanns
% http://worldclim.org/version2

% Load monthly average temperature

avgT.m1 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_01.tif'));
nan_ind = find(avgT.m1 <= -9999);

avgT.m1 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_01.tif'));
avgT.m2 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_02.tif'));
avgT.m3 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_03.tif'));
avgT.m4 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_04.tif'));
avgT.m5 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_05.tif'));
avgT.m6 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_06.tif'));
avgT.m7 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_07.tif'));
avgT.m8 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_08.tif'));
avgT.m9 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_09.tif'));
avgT.m10 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_10.tif'));
avgT.m11 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_11.tif'));
avgT.m12 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_tavg/wc2.0_2.5m_tavg_12.tif'));

sum1 = 0;
for k=1:12
    avgT.(['m' num2str(k)])(nan_ind) = NaN;
    sum1 = sum1 + avgT.(['m' num2str(k)]);
end
avgT.y = sum1/12;

figure, imagesc(avgT.y)
set(gca, 'ydir', 'normal')

% Regrid
Tres = 0.0416667;
Tlat = -90:Tres:90-0.5*Tres; % 2.5 minute resolution
Tlon = -180:Tres:180-0.5*Tres;
[Tlons, Tlats] = ndgrid(Tlon, Tlat);
avgT_rg = myregrid(Tlons, Tlats, rglons, rglats, avgT.y', 'linear');
avgT_rg = avgT_rg';

ann_T = avgT_rg;
ann_T(~landmask) = NaN;
figure, imagesc(target_lon, target_lat, ann_T)
set(gca, 'ydir', 'normal')

% find missing values of temperature
missing_T = isnan(ann_T) & landmask==1;
figure, imagesc(missing_T)

% fill missing values of temperature
ann_T = fillmissing(ann_T, 'nearest');
ann_T(~landmask) = NaN;

% size(mean([avgT(1).m; avgT(2).m]))
% 
% avgT.y = mean([avgT.m]);
% 
% figure, imagesc(target_lon, target_lat, avgT(7).m)
% set(gca, 'ydir', 'normal')

% Calculate annual average temperature
% avgT_ann = mean(avgT.m1);

% Get indices that are in IRB

% Domain boundaries:
% min      max
% x 66.15833 82.45000
% y 24.02500 37.08333

% Save the average temperature map
R = makerefmat(target_lon(1), target_lat(1), res, res);
geotiffwrite('worldclim_tavg.tif', ann_T, R)

% Also compute July_Tavg

JT_rg = myregrid(Tlons, Tlats, rglons, rglats, avgT.m7', 'linear');
JT_rg = JT_rg';

% fill missing values of temperature
JT_rg = fillmissing(JT_rg, 'nearest');
JT_rg(~landmask) = NaN;

%% Part 5 

% quartz content - use a pedotransfer function
% this one is from "Estimation of Periodic Water Balance Components...",
% cited on the UW VIC website (2016)

quartz1 = zeros(ncells, 1);
quartz2 = zeros(ncells, 1);
quartz_ptf = readtable('peters-lidard_quartz_values');
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

avgP.m1 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_01.tif'));
nan_ind = find(avgP.m1 <= -9999);

avgP.m1 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_01.tif'));
avgP.m2 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_02.tif'));
avgP.m3 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_03.tif'));
avgP.m4 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_04.tif'));
avgP.m5 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_05.tif'));
avgP.m6 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_06.tif'));
avgP.m7 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_07.tif'));
avgP.m8 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_08.tif'));
avgP.m9 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_09.tif'));
avgP.m10 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_10.tif'));
avgP.m11 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_11.tif'));
avgP.m12 = flipud(geotiffread('/Volumes/HD3/WorldClim/wc2.0_2.5m_prec/wc2.0_2.5m_prec_12.tif'));

sum1 = 0;
for k=1:12
    avgP.(['m' num2str(k)])(nan_ind) = NaN;
    sum1 = sum1 + avgP.(['m' num2str(k)]);
end
avgP.y = sum1;

figure, imagesc(avgP.y)
set(gca, 'ydir', 'normal')

% Regrid
Tres = 0.0416667;
Tlat = -90:Tres:90-0.5*Tres; % 2.5 minute resolution
Tlon = -180:Tres:180-0.5*Tres;
[Tlons, Tlats] = ndgrid(Tlon, Tlat);
avgP_rg = myregrid(Tlons, Tlats, rglons, rglats, double(avgP.y'), 'linear');
avgP_rg = avgP_rg';

ann_P = avgP_rg;
ann_P(~landmask) = NaN;
figure, imagesc(target_lon, target_lat, ann_P)
set(gca, 'ydir', 'normal')

% fill missing values of precipitation
ann_P = fillmissing(ann_P, 'nearest');
ann_P(~landmask) = NaN;

% Save the annual precipitation map
R = makerefmat(target_lon(1), target_lat(1), res, res);
geotiffwrite('worldclim_prec.tif', ann_P, R)

%% Part 8 

% organic soil properties

% organic = fraction of soil layer that is organic, from HWSD
% bdo = bulk density of the organic portion of soil, ???
% sdo = soil density of the organic portion of soil, normally 1300 kg/m3

% See Farouki (1981) for a description of the impact of organic content on
% soil thermal properties.

% Several of these parameters (also bulk density of total soil) are
% available from HWSD)
%
% Carbon content (kg C/m2)
% Area-weighted carbon content (kg C/m2)
% Organic carbon (% weight) -> can directly use this for "organic" soil
% parameter in the VIC soil parameter file

% These are not worth worrying about right now.
% I don't even think the VIC 5 image driver has the organic soils module
% built-in. I haven't seen much information about it in the literature.

%% Part 9 

% residual soil moisture

% The amount of soil moisture that cannot be removed from the soil by
% drainage or evapotranspiration (volume of residual soil moisture content / total volume of soil)
%
% Can simply be set to zero; though there might be an advantage in terms of
% model accuracy when using realistic values.
%
% Skipping this for now

%% Part 10

% Frozen soil parameters

% If using the frozen soil model, parameters must be given.
% Here are some.

%   'fs_active','frost_slope','msds'

% frost_slope = slope of the uniform distribution of soil temperature
% (optional, needed for SPATIAL_FROST) (deg. C)

% msds = maximum slope of the snow depth distribution (optional, needed for
% SPATIAL_SNOW) (m)

%% Flip the rasters to the correct direction
% 
% figure, imagesc(target_lon, target_lat, slope)
% set(gca, 'ydir', 'normal')
% 
% elev = flipud(elev);
% 
% geotiffwrite('elev.tif', elev, R)


%% Make the soil parameter file

% Two-layer soil parameter file
% varnames = {'run_cell','grid_cell','lat','lon','b_infilt','ds','dsmax', ... % all 
%     'ws','c','expt1','expt2','ksat1','ksat2','phi_s1','phi_s2', ...
%     'init_moist1','init_moist2','elev','depth1','depth2','avg_T', ...
%     'dp','bubble1','bubble2','quartz1','quartz2','bulk_dens1','bulk_dens2', ...
%     'soil_dens1','soil_dens2','organic1','organic2','bdo1','bdo2','sdo1','sdo2', ...
%     'off_gmt','wcr_fract1','wcr_fract2','wpwp_fract1','wpwp_fract2','rough','snow_rough', ...
%     'annual_prec','resid_moist1','resid_moist2','fs_active','frost_slope','msds','july_Tavg'};

% varnames = {'run_cell','grid_cell','lat','lon','b_infilt','ds','dsmax', ... % no organic, no frozen soil, no July_Tavg
%     'ws','c','expt1','expt2','ksat1','ksat2','phi_s1','phi_s2', ...
%     'init_moist1','init_moist2','elev','depth1','depth2','avg_T', ...
%     'dp','bubble1','bubble2','quartz1','quartz2','bulk_dens1','bulk_dens2', ...
%     'soil_dens1','soil_dens2', ...
%     'off_gmt','wcr_fract1','wcr_fract2','wpwp_fract1','wpwp_fract2','rough','snow_rough', ...
%     'annual_prec','resid_moist1','resid_moist2'};

% Three-layer soil parameter file (assuming top two layers are identical
% given limited data)
% varnames = {'run_cell','grid_cell','lat','lon','b_infilt','ds','dsmax', ... % all 
%     'ws','c','expt1','expt2','expt3','ksat1','ksat2','ksat3','phi_s1','phi_s2','phi_s3', ...
%     'init_moist1','init_moist2','init_moist3','elev','depth1','depth2','depth3','avg_T', ...
%     'dp','bubble1','bubble2','bubble3','quartz1','quartz2','quartz3','bulk_dens1','bulk_dens2','bulk_dens3', ...
%     'soil_dens1','soil_dens2','soil_dens3','organic1','organic2','organic3','bdo1','bdo2','bdo3','sdo1','sdo2','sdo3', ...
%     'off_gmt','wcr_fract1','wcr_fract2','wcr_fract3','wpwp_fract1','wpwp_fract2','wpwp_fract3','rough','snow_rough', ...
%     'annual_prec','resid_moist1','resid_moist2','resid_moist3','fs_active','frost_slope','msds','july_Tavg'};

varnames = {'run_cell','grid_cell','lat','lon','b_infilt','ds','dsmax', ... % three layers, no organic, no frost slope or msds
    'ws','c','expt1','expt2','expt3','ksat1','ksat2','ksat3','phi_s1','phi_s2','phi_s3', ...
    'init_moist1','init_moist2','init_moist3','elev','depth1','depth2','depth3','avg_T', ...
    'dp','bubble1','bubble2','bubble3','quartz1','quartz2','quartz3','bulk_dens1','bulk_dens2','bulk_dens3', ...
    'soil_dens1','soil_dens2','soil_dens3', ...
    'off_gmt','wcr_fract1','wcr_fract2','wcr_fract3','wpwp_fract1','wpwp_fract2','wpwp_fract3','rough','snow_rough', ...
    'annual_prec','resid_moist1','resid_moist2','resid_moist3','fs_active','july_Tavg'};

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
% ksat_simple_mean =  mean([ksat1, ksat2, ksat2],2);
% figure, subplot(2,1,1), plot(ksat_simple_mean);
% subplot(2,1,2), plot(ksat_bar);
% isequal(ksat_simple_mean, ksat_bar);
% mean(ksat_bar)
% mean(ksat_simple_mean)

dsmax = ksat_bar.*slope(landmask); % Dsmax

% There are some very high values of dsmax
% The max value is 9200*6; where there is exceptionally high conductivity and
% also a very steep slope

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

soils(:,16) = -99; % phi_s1
soils(:,17) = -99; % phi_s2
soils(:,18) = -99; % phi_s3

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

soils(:,34) = 1000*t_bd_rg(landmask); % bulk_dens1
soils(:,35) = 1000*t_bd_rg(landmask); % bulk_dens2
soils(:,36) = 1000*s_bd_rg(landmask); % bulk_dens3

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

% Use if frozen soils = TRUE
% soils(:,42) = NaN; % frost_slope
% soils(:,43) = NaN; % msds

% soils(:,44) = JT_rg(landmask); % july_Tavg
% noJT = find(isnan(JT_rg(landmask)));
% soils(noJT,44) = nanmean(JT_rg(landmask));

% (Use if organic=TRUE)
% soils(:,31) = t_org_rg(landmask); % organic1
% soils(:,32) = s_org_rg(landmask); % organic2
% soils(:,33) = NaN; % bdo1
% soils(:,34) = NaN; % bdo2
% soils(:,35) = 1300; % sdo1
% soils(:,36) = 1300; % sdo2

% soils(:,37) = NaN; % off_gmt
% soils(:,38) = wcr_fract1; % wcr_fract1
% soils(:,39) = wcr_fract2; % wcr_fract2
% soils(:,40) = wpwp_fract1; % wpwp_fract1
% soils(:,41) = wpwp_fract2; % wpwp_fract2
% soils(:,42) = NaN; % rough
% soils(:,43) = NaN; % snow_rough
% soils(:,44) = NaN; % annual_prec
% soils(:,45) = NaN; % resid_moist1
% soils(:,46) = NaN; % resid_moist2
% soils(:,47) = NaN; % fs_active
% soils(:,48) = NaN; % frost_slope
% soils(:,49) = NaN; % msds

%% Write out the soil parameter file

% should probably switch to using write_soils()

% soils(isnan(soils)) = -99; % nodata value

fstring = ['%.' num2str(grid_decimal) 'f'];

% This is used for the Livneh (2013) soil parameter file
% fspec = ['%d %d ' fstring ' ' fstring ' %.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %d\n'];

% This is used for VIC-2L w no optional variables (41 columns in the soil parameter file).
% fspec = ['%d %d ' fstring ' ' fstring ' ' '%.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %d %d %.3f %.3f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d\n'];

% This is used VIC-3L w no optional variables (54 columns in the soil parameter file).
fspec = ['%d %d ' fstring ' ' fstring ' ' '%.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %d %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %.2f\n'];

fID = fopen(fullfile(outdir, savename),'w');
% fID = fopen(savename,'w');
fprintf(fID, fspec, soils1');
fclose(fID);
display(['Soils data saved as ' fullfile(outdir, savename)])

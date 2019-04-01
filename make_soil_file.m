% Make Soil File
%
% Script that uses HWSD data to make a soil parameter file for VIC 4 or 5
%
% Dependencies
% addpath(genpath('/Users/jschap/Documents/MATLAB/from_file_exchange'))
% myregrid.m, soil_classification.m, pedotransfer_table.txt
%
% To do
%
% 1. Incorporate better elevation data (SRTM, it's under "static data" on Elqui/hoffman2)
%    Also, the HWSD has elevation data and land mask files available on its
%    website, under "terrain"
% 2. Fill gaps in coverage (by interpolation or finding different data)
% 3. Write up procedure of making the soil parameter file
%
% Updated 3/12/2019 JRS - added fs_active to base set of soil parameters
% since it's required to run VIC, even if not using frozen soils

%% INPUTS

addpath(genpath('/Users/jschap/Documents/MATLAB/from_file_exchange'))

res = 1/16;
grid_decimal = 5; % number of decimal points for forcing file names
savename = 'mysoilparam.txt';

%% Part 1 - load HWSD data

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
target_lon = min(lon):0.0625:max(lon);
target_lat = min(lat):0.0625:max(lat);
[lons, lats] = ndgrid(lon, lat);
[rglons, rglats] = ndgrid(target_lon, target_lat);

t_sand_rg = myregrid(lons, lats, rglons, rglats, t_sand, 'linear')';
t_clay_rg = myregrid(lons, lats, rglons, rglats, t_clay, 'linear')';
s_sand_rg = myregrid(lons, lats, rglons, rglats, s_sand, 'linear')';
s_clay_rg = myregrid(lons, lats, rglons, rglats, s_clay, 'linear')';
t_bd_rg = myregrid(lons, lats, rglons, rglats, t_bd, 'linear')';
s_bd_rg = myregrid(lons, lats, rglons, rglats, s_bd, 'linear')';
t_org_rg = myregrid(lons, lats, rglons, rglats, t_org, 'linear')';
s_org_rg = myregrid(lons, lats, rglons, rglats, s_org, 'linear')';

% figure, imagesc(target_lon, target_lat, t_sand_rg)

% Make a land mask from HWSD. Use the field "ISSOIL"

is_soil = ncread('../HWSD/HWSD_1247/data/ISSOIL.nc4', 'ISSOIL');

% 1 - 5.7e6 pixels -> indicates a soil mapping unit is soil
% 0 - 4.8e5 pixels -> indicates a soil mapping unit is not soil
% NaN - 2.0e7 pixels -> ocean

is_soil_rg = myregrid(lons, lats, rglons, rglats, is_soil, 'linear');
is_soil_rg = is_soil_rg';

% figure, imagesc(target_lon, target_lat, is_soil_rg)
% set(gca, 'ydir', 'normal')

% Since we've upscaled the soil mask, some pixels are partially soil and
% partially not soil. Need to account for this. Simplest thing to do:
% assume that all partial soil units are soil.

is_soil_rg(is_soil_rg>0) = 1;

% test = is_soil_rg;
% test(is_soil_rg>0) = 1;
% figure, imagesc(test)

%% Part 2 - load elevation and land fraction data

% (Preprocessing - make geotiffs from TBASE data using R)

% Load elevation data
[elev, R] = geotiffread('/Volumes/HD3/SWOTDA/Data/JISAO/tbase_elev_0.25.tif');

elat = -89.875:0.25:89.875;
elon = -179.875:0.25:179.875;

% elat = ncread('./Data/JISAO/elev.0.25-deg.nc', 'lat'); 
% elon = ncread('./Data/JISAO/elev.0.25-deg.nc', 'lon');
% elon = 180 - elon;

elev=fliplr(elev');
elev(elev==-9999) = NaN; % nodata values

% Regrid elevation data
[elons, elats] = ndgrid(elon, elat);
elev_rg = myregrid(elons, elats, rglons, rglats, elev, 'linear');
elev_rg = elev_rg';

% Repeat for land cover fraction data
% [fract, R] = geotiffread('/Volumes/HD3/SWOTDA/Data/JISAO/tbase_fract_0.25.tif');
% fract=fliplr(fract');
% fract(fract==-9999) = NaN; % nodata values
% fract_rg = myregrid(elons, elats, rglons, rglats, fract);
% fract_rg = fract_rg';

% figure
% subplot(1,3,1)
% imagesc(target_lon, target_lat, t_sand_rg)
% set(gca, 'ydir', 'normal')
% title('Percent sand (topsoil)')
% subplot(1,3,2)
% imagesc(target_lon, target_lat, elev_rg)
% set(gca, 'ydir', 'normal')
% title('Elevation')
% 
% subplot(1,3,3)
% imagesc(target_lon, target_lat, is_soil_rg)
% set(gca, 'ydir', 'normal')
% title('Soil mapping units')

% imagesc(target_lon, target_lat, fract_rg)
% set(gca, 'ydir', 'normal')
% title('Land cover fraction')

% Make a land mask

% hwsd_landmask = geotiffread('/Volumes/HD3/HWSD/hwsd_landmask_1_16.tif');
% hwsd_landmask(hwsd_landmask == 255) = NaN;
% hwsd_landmask = flipud(hwsd_landmask);
% figure, imagesc(hwsd_landmask)
% this should be very similar to the ISSOIL mask bc both are from the same data source, ultimately
% SRTM has limited spatial coverage. Should really use a more global
% dataset like MERIT DEM or even GTOPO30.
% 
% For now, let's try SRTM where it is available, and GTOPO for the rest of
% the world
%
% Note: we have elevation data for the HMA region saved in 1 degree by 1
% degree tiles on Elqui: /Volumes/elqui_hd5/PROJECTS/SWE_REANALYSIS/HMA/INPUT_DATA

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Constructing a 1/16 degree global elevation map of the world from SRTM





% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Load a landmask from MODIS

% gsw_landmask = geotiffread('/Volumes/HD3/VICParametersGlobal/High_Resolution/landmask_from_GSW.tif');
% gsw_landmask(gsw_landmask<0) = NaN;
% figure, imagesc(gsw_landmask)


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

landmask = is_soil_rg; 
landmask(isnan(landmask)) = 0;
landmask = logical(landmask);
landcells = find(landmask); % where HWSD has data..., these are considered to be land..., should fill in missing data
ncells = length(landcells);

% landmask = fract_rg>0;
% ncells = sum(landmask(:));

% sanity check - what percent of Earth is covered by water?
% 100-100*ncells/(5760*2880) % About 64 percent by my calc, but it should be about 71 percent, in reality... Could be bc of fraction vs. whole cell

% Exclude Antarctica (below 60 S latitude)
% landmask(:, target_lat<=-60) = 0;
% ncells = sum(landmask(:));
% 
% land_ind = find(landmask);

figure
imagesc(target_lon, target_lat, landmask)
xlabel('Longitude')
ylabel('Latitude')
title('Land mask from HWSD')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)

elev1 = elev_rg;
elev1(~landmask) = NaN;

R = makerefmat(-179.975, -89.975, res, res);

% arcgridwrite('./Data/JISAO/elev1.asc', target_lon, target_lat, elev1')
% geotiffwrite('./Data/JISAO/elev1.tif', elev1, R);

% computing slope in Matlab
[~, slope, ~, ~] = gradientm(elev1, R); % there are a lot of missing values

% Address missing values of slope where there is a land pixel:
noslope = find(isnan(slope(landmask)));
% slope(isnan(slope(landmask))) = 0;

% sum(isnan(slope(landmask'))) % about 10% of them are missing
% missing_ind = find(isnan(slope(landmask')));
% slope(missing_ind) = 1e-4;
% replace missing values with small values

% computing slope in gdal
% gdaldem slope -s 111120 -p -compute_edges  N34E069.hgt N34E069_slope.tif
% noslope = find(isnan(slope(landmask)));

% figure, imagesc(target_lon, target_lat, elev1)
% 
% [slope, R] = geotiffread('./Data/JISAO/slope.tif');
% slope(slope==-9999) = NaN;

% calculate lonlat for soil parameter file
xyz_upscaled = raster2xyz(target_lon, target_lat, ones(size(t_clay_rg)));
lonlat = xyz_upscaled(landcells,1:2);

%% Part 3

% Problem: there are soil data pixels with nodata because the 1/4 degree
% land mask from the DEM classifies some non-land pixels from the HWSD data
% as land. 
%
% Potential solution: use the HWSD data to derive the land mask - however,
% sometimes there are NaN data even though there is land, just because
% there are no soil data available there. Need the land mask used to
% derive the HWSD dataset.
% 
% Done, using ISSOILS.nc4 file

% Classify soils using USGS soil triangle
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

%% Part 5 

% quartz content - use a pedotransfer function

quartz1 = zeros(ncells, 1);
quartz2 = zeros(ncells, 1);
quartz_ptf = [0.95, 0.85, 0.69, 0.19, 0.05, 0.41, 0.61, 0.09, 0.25, 0.50, 0.08, 0.25];
for sc=1:12
    quartz1(s_sc_int == sc) = quartz_ptf(sc);    
    quartz2(t_sc_int == sc) = quartz_ptf(sc);    
end

%% Part 6

% off_gmt parameter

% tz = shaperead('/Volumes/HD3/TimeZones/ne_10m_time_zones/ne_10m_time_zones.shp');
% convert to raster using gdal/ogr

% gdal_rasterize -te -179.9750 -89.9750 179.9750 89.9750 -ts 5760 2880 -a zone /Volumes/HD3/TimeZones/ne_10m_time_zones/ne_10m_time_zones.shp time_zones.tif

tz = geotiffread('/Volumes/HD3/TimeZones/time_zones.tif');
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

%% Make the soil parameter file

% Two-layer soil parameter file
% varnames = {'run_cell','grid_cell','lat','lon','b_infilt','ds','dsmax', ... % all 
%     'ws','c','expt1','expt2','ksat1','ksat2','phi_s1','phi_s2', ...
%     'init_moist1','init_moist2','elev','depth1','depth2','avg_T', ...
%     'dp','bubble1','bubble2','quartz1','quartz2','bulk_dens1','bulk_dens2', ...
%     'soil_dens1','soil_dens2','organic1','organic2','bdo1','bdo2','sdo1','sdo2', ...
%     'off_gmt','wcr_fract1','wcr_fract2','wpwp_fract1','wpwp_fract2','rough','snow_rough', ...
%     'annual_prec','resid_moist1','resid_moist2','fs_active','frost_slope','msds','july_Tavg'};

varnames = {'run_cell','grid_cell','lat','lon','b_infilt','ds','dsmax', ... % no organic, no frozen soil, no July_Tavg
    'ws','c','expt1','expt2','ksat1','ksat2','phi_s1','phi_s2', ...
    'init_moist1','init_moist2','elev','depth1','depth2','avg_T', ...
    'dp','bubble1','bubble2','quartz1','quartz2','bulk_dens1','bulk_dens2', ...
    'soil_dens1','soil_dens2', ...
    'off_gmt','wcr_fract1','wcr_fract2','wpwp_fract1','wpwp_fract2','rough','snow_rough', ...
    'annual_prec','resid_moist1','resid_moist2'};

nvars = length(varnames);

soils = zeros(ncells, nvars); % this is about 1.5 GB, so pretty big

soils(:,1) = 1; % run_cell
soils(:,2) = 1:ncells; % grid_cell
soils(:,3) = lonlat(:,2); % lat
soils(:,4) = lonlat(:,1); % lon
soils(:,5) = 0.2; % b_infilt
soils(:,6) = 0.001; % Ds

dsmax = mean([ksat1, ksat2],2).*slope(landmask); % Dsmax
% noslope = find(isnan(slope(landmask)));
% soils(noslope,7) = 0; 
soils(:,7) = fillmissing(dsmax, 'nearest');
% filling missing values with the nearest non-missing value
%
% There are some very high values of dsmax
% The max value is 9200*6; where there is exceptionally high conductivity and
% also a very steep slope

soils(:,8) = 0.9; % Ws
soils(:,9) = 2; % c
soils(:,10) = expt1; % expt1
soils(:,11) = expt2; % expt2
soils(:,12) = ksat1; % Ksat1
soils(:,13) = ksat2; % Ksat2
soils(:,14) = -99; % phi_s1
soils(:,15) = -99; % phi_s2
soils(:,16) = 1000.*wcr_fract1.*0.3.*porosity1; % init_moist1
soils(:,17) = 1000.*wcr_fract2.*0.7.*porosity2; % init_moist2

% soils(:,18) = elev1(landmask); % elev
% noelev = find(isnan(elev1(landmask)));
% soils(noelev,18) = 0;
soils(:,18) = fillmissing(elev1(landmask), 'nearest'); % elev

soils(:,19) = 30/100; % depth1 (m)
soils(:,20) = 100/100; % depth2 (m)

soils(:,21) = fillmissing(avgT_rg(landmask), 'nearest'); % avg_T
% soils(:,21) = avgT_rg(landmask); % avg_T
% noT = find(isnan(avgT_rg(landmask)));
% soils(noT,21) = 15; % if no data, assume they're 15 degrees C

soils(:,22) = 4; % dp
soils(:,23) = 0.32*expt1 + 4.3; % bubble1
soils(:,24) = 0.32*expt2 + 4.3; % bubble2
soils(:,25) = quartz1; % quartz1
soils(:,26) = quartz2; % quartz2

soils(:,27) = fillmissing(1000*t_bd_rg(landmask), 'nearest'); % bulk_dens1
soils(:,28) = fillmissing(1000*s_bd_rg(landmask), 'nearest'); % bulk_dens2

% soils(:,27) = bulk_dens1; % bulk_dens1
% soils(:,28) = bulk_dens2; % bulk_dens2

soils(:,29) = 2685; % soil_dens1
soils(:,30) = 2685; % soil_dens2
soils(:,31) = tz(landmask); % off_gmt
soils(:,32) = wcr_fract1; % wcr_fract1
soils(:,33) = wcr_fract2; % wcr_fract2
soils(:,34) = wpwp_fract1; % wpwp_fract1
soils(:,35) = wpwp_fract2; % wpwp_fract2
soils(:,36) = 0.001; % rough
soils(:,37) = 0.0005; % snow_rough

% soils(:,38) = avgP_rg(landmask); % annual_prec
soils(:,38) = fillmissing(avgP_rg(landmask), 'nearest');
% noprecip = find(isnan(avgT_rg(landmask)));
% soils(noprecip,21) = -99; 

soils(:,39) = 0; % resid_moist1 (leaving it as 0...)
soils(:,40) = 0; % resid_moist2
soils(:,41) = 0; % fs_active

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

% soils(:,50) = JT_rg(landmask); % july_Tavg
% noJT = find(isnan(JT_rg(landmask)));
% soils(noJT,50) = nanmean(JT_rg(landmask));

%% Write out the soil parameter file

% soils(isnan(soils)) = -99; % nodata value

fstring = ['%.' num2str(grid_decimal) 'f'];

% This is used for the Livneh (2013) soil parameter file
% fspec = ['%d %d ' fstring ' ' fstring ' %.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %.3f %.3f %d %d %d %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d %d\n'];

% This is used for VIC w no optional variables (41 columns in the soil parameter file).
fspec = ['%d %d ' fstring ' ' fstring ' ' '%.4f %.4f %.4f %.4f %d %.3f %.3f %.3f %.3f %d %d %.3f %.3f %.2f %.2f %.2f %d %d %.3f %.3f %.3f %.3f %.2f %.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f %.2f %d %d %d %d\n'];

fID = fopen(savename,'w');
fprintf(fID, fspec, soils');
fclose(fID);
display(['Soils data saved as ' savename])

%% Scrap

% Get land pixels
% land_ind = find(~isnan(t_sand_rg));
% land_mask = ~isnan(t_sand_rg);
% ncells = length(land_ind);
% %

% Do some operations in GDAL to get the domain in a reasonable size and
% format, such as a geotiff file for the study area extent
% gdal_translate -a_srs epsg:4326 ./HWSD_RASTER/hwsd.bil ./HWSD_RASTER/hwsd_geog.tif

% Domain: [66.158333, 24.025, 82.45, 37.083333] % [xll, yll, xur, yur]
% gdal_translate

% load soil parameter data for the region of interest
% soils = geotiffread('../HWSD/HWSD_RASTER/hwsd_irb.tif');

% Each code in the raster links to the attribute database
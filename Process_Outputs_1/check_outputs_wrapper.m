% VIC outputs workflow
%
% Perform various tasks to process the VIC model outputs and analyze the
% results. Can handle VIC Classic or Image drivers.
%
% Last revised 10/7/2019 JRS

clearvars -except soils_vg
addpath(genpath('/Users/jschap/Documents/Codes/VICMATLAB'))

%% Classic Driver --------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

% 1-degree test
% wb_out_dir = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/34N_75E/Classic_Driver_Tests/Raw_WB/daily/wb';
% eb_out_dir = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/34N_75E/Classic_Driver_Tests/Raw_WB/daily/eb';
% results_dir = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/34N_75E/Classic_Driver_Tests/Processed_WB/daily';

% full IRB
wb_out_dir = '/Volumes/HD4/SWOTDA/Data/IRB/Classic/Raw_WB/wb';
eb_out_dir = '/Volumes/HD4/SWOTDA/Data/IRB/Classic/Raw_WB/eb';
results_dir = '/Volumes/HD4/SWOTDA/IRB/Classic/Processed_WB';
figdir = '/Volumes/HD4/SWOTDA/IRB/Classic/Processed_WB/Figures';

% tuolumne test
% wb_out_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_VICGlobal/L2013/Raw_EB_FS/wb';
% eb_out_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_VICGlobal/L2013/Raw_EB_FS/eb';
% results_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_VICGlobal/L2013/Processed_EB_FS';

% UMRB
% wb_out_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_L2015/L2013/Raw_EB_FS_1980-2011/wb';
% eb_out_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_L2015/L2013/Raw_EB_FS_1980-2011/eb';
% results_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_L2015/L2013/Processed_1980-2011';

timestep_out = 'daily';

info = get_vic_run_metadata(wb_out_dir, eb_out_dir, timestep_out);

mkdir(results_dir)
save(fullfile(results_dir, 'vic_run_metadata.mat'), 'info');

%% Plot a time series for a given location

lat = 47; lon = -96;
figure; plot_one_grid_cell(lat, lon, info, 'eb', 'OUT_AIR_TEMP')

grid on
hold on

lat = 38; lon = -87;
plot_one_grid_cell(lat, lon, info, 'eb', 'OUT_AIR_TEMP')

legend('cold cell','hot cell','location','northwest')

%% Calculate time-average maps and area-average time series

parpool(1);
[wb_avg_map, wb_avg_ts, wb_sum_ts, ~, ~] = readVIC_ds(info.wb_out_dir, length(info.wbvars), info.ncells, info.nt);
[eb_avg_map, eb_avg_ts, eb_sum_ts, ~, ~] = readVIC_ds(info.eb_out_dir, length(info.ebvars), info.ncells, info.nt);

% A = load('/Volumes/HD4/SWOTDA/Data/UMRB/Classic_Livneh_met_L15/Raw/Processed/readVIC_outputs.mat');
% wb_avg_map = A.avg_map;
% wb_avg_ts = A.avg_ts;
% wb_sum_ts = A.sum_ts;

%% Assemble the data into an easy to deal with format

OUTPUTS = make_outputs_struct(info, wb_avg_ts, wb_avg_map, eb_avg_ts, eb_avg_map);
mkdir(results_dir)
save(fullfile(results_dir, 'vic_outputs_summarized_1980_daily.mat'), 'OUTPUTS');
% OUTPUTS_VG = load('/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_VICGlobal/Processed_WB_SB/vic_outputs_summarized_1999_daily.mat');
% OUTPUTS_VG = OUTPUTS_VG.OUTPUTS;

%% Plots

OUTPUTS_VG = load('/Volumes/HD4/SWOTDA/Data/UMRB/Classic_VG_met_L15/Processed/vic_outputs_summarized_1999_daily.mat');
OUTPUTS_VG = OUTPUTS_VG.OUTPUTS;

OUTPUTS = load('/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_L2015/L2013/Processed_1980-2011/vic_outputs_summarized_1999_daily.mat');
OUTPUTS = OUTPUTS.OUTPUTS;

% figdir = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/Classic/Figures';
% figdir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Classic_L2015/Figures';
mkdir(figdir)

plot_spatial_avg_ts(OUTPUTS, figdir);
plot_time_avg_maps(OUTPUTS, figdir);

figure, 
plotraster(OUTPUTS_VG.lon, OUTPUTS_VG.lat, OUTPUTS_VG.EB.maps.OUT_AIR_TEMP, 'Temperature','','')

figure, 
plotraster(OUTPUTS_VG.lon, OUTPUTS_VG.lat, OUTPUTS_VG.WB.maps.OUT_SWE, 'SWE','','')

figure
plotraster(OUTPUTS_VG.lon, OUTPUTS_VG.lat, OUTPUTS_VG.WB.maps.OUT_BASEFLOW, 'Baseflow','','')

% plot_difference_ts(OUTPUTS, OUTPUTS_VG, figdir)
% plot_difference_maps(OUTPUTS, OUTPUTS_VG, figdir, 0)

% Write out GeoTiffs
R = makerefmat(min(OUTPUTS.lon), min(OUTPUTS.lat), 1/16, 1/16);
geotiffwrite(fullfile(figdir, 'average_precipitation.tif'), flipud(OUTPUTS.WB.maps.OUT_PREC), R)
geotiffwrite(fullfile(figdir, 'average_evaporation.tif'), flipud(OUTPUTS.WB.maps.OUT_EVAP), R)
geotiffwrite(fullfile(figdir, 'average_runoff.tif'), flipud(OUTPUTS.WB.maps.OUT_RUNOFF), R)
geotiffwrite(fullfile(figdir, 'average_baseflow.tif'), flipud(OUTPUTS.WB.maps.OUT_BASEFLOW), R)
geotiffwrite(fullfile(figdir, 'average_temperature.tif'), flipud(OUTPUTS.EB.maps.OUT_AIR_TEMP), R)

%% Image driver --------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% The rotation of these images is correct based on comparison with 
% WorldClim precipitation and temperature patterns -- 10/14/2019 JRS

wkdir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Image_VICGlobal/L2013/EB_FS_SB';
fluxfile = fullfile(wkdir, 'fluxes.1999-01-01.nc');
info_image = get_vic_run_metadata_image(fluxfile);

% Compute time average maps and spatially-average time series
[avg_ts, avg_map] = calc_average_image(fluxfile, info_image);

% Put it all in a structure just like for classic mode
OUTPUTS_IM = make_outputs_struct_image(info_image, avg_ts, avg_map);
timevector = OUTPUTS_VG.time;

save(fullfile(wkdir, 'vic_outputs_summarized_1999_daily.mat'), 'OUTPUTS_IM', 'timevector')

% Plots
figdir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Image_VICGlobal/L2013/Figures';
mkdir(figdir)
fontsize = 18;
height1 = 500;
width1 = 700;

% Add snow data to OUTPUTS structure
snowfile = fullfile(wkdir, 'snow.1999-01-01.nc');
info_image_snow = get_vic_run_metadata_image(snowfile);
[~, avg_map_snow] = calc_average_image(snowfile, info_image_snow);
OUTPUTS_IM.avg_map_snow = avg_map_snow;
save(fullfile(wkdir, 'vic_outputs_summarized_1999_daily.mat'), 'OUTPUTS_IM')

figure, 
plotraster(OUTPUTS_IM.lon, OUTPUTS_IM.lat, OUTPUTS_IM.avg_map.OUT_AIR_TEMP,'Temperature','','')

figure, 
plotraster(OUTPUTS_IM.lon, OUTPUTS_IM.lat, OUTPUTS_IM.avg_map.OUT_PREC,'Precipitation','','')

figure, 
plotraster(OUTPUTS_IM.lon, OUTPUTS_IM.lat, fliplr(flipud(OUTPUTS_IM.avg_map.OUT_PREC)),'Precipitation','','')

figure, 
plotraster(OUTPUTS_IM.lon, OUTPUTS_IM.lat, OUTPUTS_IM.avg_map.OUT_R_NET,'Rnet','','')

% Time average maps
% Water balance variables
% Fluxes
f5 = figure;
set(f5, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
plotraster(OUTPUTS_IM.lon, OUTPUTS_IM.lat, OUTPUTS_IM.avg_map.OUT_PREC, 'Precipitation (mm)', 'Lon', 'Lat')
subplot(2,2,2)
plotraster(OUTPUTS_IM.lon, OUTPUTS_IM.lat, OUTPUTS_IM.avg_map.OUT_EVAP, 'Evaporation (mm)', 'Lon', 'Lat')
subplot(2,2,3)
plotraster(OUTPUTS_IM.lon, OUTPUTS_IM.lat, OUTPUTS_IM.avg_map.OUT_RUNOFF, 'Runoff (mm)', 'Lon', 'Lat')
subplot(2,2,4)
plotraster(OUTPUTS_IM.lon, OUTPUTS_IM.lat, OUTPUTS_IM.avg_map.OUT_BASEFLOW, 'Baseflow (mm)', 'Lon', 'Lat')
% saveas(f5, fullfile(figdir, 'water_balance_fluxes_maps.png'))

% Area-average time series
% Water balance variables
% Fluxes
f6 = figure;
upper1 = 60;
set(f6, 'Position',  [100, 100, 100+width1, 100+height1])
subplot(2,2,1)
jsplot(timevector, OUTPUTS_IM.avg_ts.OUT_PREC, 'Precipitation (mm)', 'Time', 'Precipitation', fontsize);
grid on
ylim([0,upper1])
subplot(2,2,2)
jsplot(timevector, OUTPUTS_IM.avg_ts.OUT_EVAP, 'Evaporation (mm)', 'Time', 'Evaporation', fontsize);
grid on
ylim([0,upper1])
subplot(2,2,3)
jsplot(timevector, OUTPUTS_IM.avg_ts.OUT_RUNOFF, 'Runoff (mm)', 'Time', 'Runoff', fontsize);
grid on
ylim([0,upper1])
subplot(2,2,4)
jsplot(timevector, OUTPUTS_IM.avg_ts.OUT_BASEFLOW, 'Baseflow (mm)', 'Time', 'Baseflow', fontsize);
grid on
ylim([0,upper1])
% saveas(f6, fullfile(figdir, 'water_balance_fluxes_ts.png'))

% Write out GeoTiffs
figdir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Image_VICGlobal/L2013/Figures';
R = makerefmat(min(OUTPUTS_IM.lon), min(OUTPUTS_IM.lat), 1/16,1/16);
geotiffwrite(fullfile(figdir, 'average_temperature.tif'), OUTPUTS_IM.avg_map.OUT_AIR_TEMP, R)
geotiffwrite(fullfile(figdir, 'average_precipitation.tif'), OUTPUTS_IM.avg_map.OUT_PREC, R)
geotiffwrite(fullfile(figdir, 'average_evaporation.tif'), OUTPUTS_IM.avg_map.OUT_EVAP, R)
geotiffwrite(fullfile(figdir, 'average_runoff.tif'), OUTPUTS_IM.avg_map.OUT_RUNOFF, R)
geotiffwrite(fullfile(figdir, 'average_baseflow.tif'), OUTPUTS_IM.avg_map.OUT_BASEFLOW, R)
geotiffwrite(fullfile(figdir, 'average_swe.tif'), OUTPUTS_IM.avg_map_snow.OUT_SWE, R)
% no need to flipud these; has to do w their origin as NetCDF data

%% Plots for report

OUTPUTS_VG = load('/Volumes/HD4/SWOTDA/Data/Tuolumne/Classic_VICGlobal/Processed_EB_SB/vic_outputs_summarized_1999_daily.mat');
OUTPUTS_VG = OUTPUTS_VG.OUTPUTS;

OUTPUTS_L15 = load('/Volumes/HD4/SWOTDA/Data/Tuolumne/Classic_L2015/Processed_EB_SB/vic_outputs_summarized_1999_daily.mat');
OUTPUTS_L15 = OUTPUTS_L15.OUTPUTS;

OUTPUTS_IM = load('/Volumes/HD4/SWOTDA/Data/Tuolumne/Image_VICGlobal/Results_EB_SB_L/vic_outputs_summarized_1999_daily.mat');
OUTPUTS_IM = OUTPUTS_IM.OUTPUTS_IM;

%%
figure, subplot(4,1,1)
fontsize = 18;
timevector = OUTPUTS_VG.time;
jsplot(timevector, OUTPUTS_VG.WB.ts.OUT_RUNOFF, 'Runoff', 'Time', 'Runoff (mm)', fontsize);
grid on
ylim([0,20])
hold on
plot(timevector, OUTPUTS_L15.WB.ts.OUT_RUNOFF, 'k');
% plot(timevector, OUTPUTS_IM.avg_ts.OUT_RUNOFF, 'r');
% legend('VICGlobal (Classic)', 'L2015', 'VICGlobal (Image)', 'Location', 'northeast')
legend('VICGlobal (Classic)', 'L2015', 'Location', 'northeast')

subplot(4,1,2)
jsplot(timevector, OUTPUTS_VG.WB.ts.OUT_BASEFLOW, 'Baseflow', 'Time', 'Baseflow (mm)', fontsize);
grid on
ylim([0,20])
hold on
plot(timevector, OUTPUTS_L15.WB.ts.OUT_BASEFLOW, 'k');
% plot(timevector, OUTPUTS_IM.avg_ts.OUT_BASEFLOW, 'r');

subplot(4,1,3)
jsplot(timevector, OUTPUTS_VG.WB.ts.OUT_PREC, 'Precipitation', 'Time', 'Precipitation (mm)', fontsize);
grid on
ylim([0,30])
hold on
plot(timevector, OUTPUTS_L15.WB.ts.OUT_PREC, 'k');
% plot(timevector, OUTPUTS_IM.avg_ts.OUT_PREC, 'r');

subplot(4,1,4)
jsplot(timevector, OUTPUTS_VG.EB.ts.OUT_SENSIBLE, 'Sensible heat flux', 'Time', 'H (W/m^2)', fontsize);
grid on
ylim([-50,150])
hold on
plot(timevector, OUTPUTS_L15.EB.ts.OUT_SENSIBLE, 'k');
% plot(timevector, OUTPUTS_IM.avg_ts.OUT_SENSIBLE, 'r');

%% SWE and net radiation difference maps (L2015 vs. VG)

lon = OUTPUTS_VG.lon;
lat = OUTPUTS_VG.lat;
percent_diff = 0;

f1 = figure(1);
subplot(3,2,1)
plotraster(lon, lat, OUTPUTS_VG.WB.maps.OUT_SWE, 'VICGlobal Classic (SWE, mm)', 'Lon','Lat');
% caxis([0,60])
subplot(3,2,3)
plotraster(lon, lat, OUTPUTS_L15.WB.maps.OUT_SWE, 'L2015 (SWE, mm)', 'Lon','Lat')
subplot(3,2,5)
plotdiff_map(lon, lat, OUTPUTS_VG.WB.maps.OUT_SWE, OUTPUTS_L15.WB.maps.OUT_SWE, ...
    'L2015 - VICGlobal Classic (SWE, mm)','Lon','Lat', percent_diff)

subplot(3,2,2)
plotraster(lon, lat, OUTPUTS_VG.EB.maps.OUT_R_NET, 'VICGlobal Classic (Rnet, W/m^2)', 'Lon','Lat');
subplot(3,2,4)
plotraster(lon, lat, OUTPUTS_L15.EB.maps.OUT_R_NET, 'L2015 (Rnet, W/m^2)', 'Lon','Lat')
subplot(3,2,6)
plotdiff_map(lon, lat, OUTPUTS_VG.EB.maps.OUT_R_NET, OUTPUTS_L15.EB.maps.OUT_R_NET, ...
    'L2015 - VICGlobal Classic (Rnet, W/m^2)','Lon','Lat', percent_diff)

% SWE and net radiation difference maps (VG classic vs. VG image)
f2 = figure(2);
subplot(3,2,1)
plotraster(lon, lat, OUTPUTS_VG.WB.maps.OUT_SWE, 'VICGlobal Classic (SWE, mm)', 'Lon','Lat');
% caxis([0,60])
subplot(3,2,3)
plotraster(lon, lat, flipud(OUTPUTS_IM.avg_map_snow.OUT_SWE), 'VICGlobal Image (SWE, mm)', 'Lon','Lat')
subplot(3,2,5)
plotdiff_map(lon, lat, OUTPUTS_VG.WB.maps.OUT_SWE, flipud(OUTPUTS_IM.avg_map_snow.OUT_SWE), ...
    'Image - Classic (SWE, mm)','Lon','Lat', percent_diff)

subplot(3,2,2)
plotraster(lon, lat, OUTPUTS_VG.EB.maps.OUT_R_NET, 'VICGlobal Classic (Rnet, W/m^2)', 'Lon','Lat');
subplot(3,2,4)
plotraster(lon, lat, flipud(OUTPUTS_IM.avg_map.OUT_R_NET), 'VICGlobal Image (Rnet, W/m^2)', 'Lon','Lat')
subplot(3,2,6)
plotdiff_map(lon, lat, OUTPUTS_VG.EB.maps.OUT_R_NET, flipud(OUTPUTS_IM.avg_map.OUT_R_NET), ...
    'Image - Classic (Rnet, W/m^2)','Lon','Lat', percent_diff)

%%%%%%%%%%%%%% Not used %%%%%%%%%%%%%%
% SWE and net radiation difference maps (VG classic vs. VG image)
orig_lon = unique(OUTPUTS_IM.lon);
orig_lat = unique(OUTPUTS_IM.lat);
rglon = unique(OUTPUTS_VG.lon);
rglat = unique(OUTPUTS_VG.lat);
[orig_lons, orig_lats] = ndgrid(orig_lon, orig_lat);
[rglons, rglats] = ndgrid(rglon, rglat);

f2 = figure(2);
subplot(3,2,1)
swe_image = fliplr(OUTPUTS_IM.avg_map_snow.OUT_SWE);
swe_image_resampled = myregrid(orig_lons, orig_lats, rglons, rglats, swe_image', 'linear')';
plotraster(lon, lat, OUTPUTS_VG.WB.maps.OUT_SWE, 'VICGlobal Classic (SWE, mm)', 'Lon','Lat');
caxis([0,60])
subplot(3,2,3)
plotraster(lon, lat, swe_image_resampled, 'VICGlobal Image (SWE, mm)', 'Lon','Lat')
caxis([0,60])
subplot(3,2,5)
plotdiff_map(lon, lat, OUTPUTS_VG.WB.maps.OUT_SWE, swe_image_resampled, ...
    'Image - VICGlobal Classic (SWE, mm)','Lon','Lat', percent_diff)
caxis([0,60])

subplot(3,2,2)
rnet_image = fliplr(OUTPUTS_IM.avg_map.OUT_R_NET);
rnet_image_resampled = myregrid(orig_lons, orig_lats, rglons, rglats, rnet_image', 'linear')';
plotraster(lon, lat, OUTPUTS_VG.EB.maps.OUT_R_NET, 'VICGlobal Classic (Rnet, W/m^2)', 'Lon','Lat');
caxis([0,120])
subplot(3,2,4)
plotraster(lon, lat, rnet_image_resampled, 'VICGlobal Image (Rnet, W/m^2)', 'Lon','Lat')
caxis([0,120])
subplot(3,2,6)
plotdiff_map(lon, lat, OUTPUTS_VG.EB.maps.OUT_R_NET, rnet_image_resampled, ...
    'Image - Classic (Rnet, W/m^2)','Lon','Lat', percent_diff)
caxis([0,120])
%%%%%%%%%%%%%% Not used %%%%%%%%%%%%%%


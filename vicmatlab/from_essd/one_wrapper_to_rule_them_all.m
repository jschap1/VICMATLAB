% Wrapper for UCRB-Livneh comparison for ESSD paper
%
% Written 1/28/2020 JRS

%% Compare summed runoff and baseflow between Livneh and VICGLobal simulations

livneh = load('/Volumes/HD4/SWOTDA/Data/Colorado/L15/livneh_q.mat');
indir1 = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/processed';
indir2 = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_2000-2011_VGMOD/processed';
vicglobal = load(fullfile(indir1, 'summed_q_standard.mat'));
vicglobal_mod = load(fullfile(indir2, 'summed_q_standard.mat'));

% report only a few years of simulation data so it is easier to see how L15
% and VICGlobal results compare
i1 = find(livneh.timevector == datetime(2000,10,1));
i2 = find(livneh.timevector == datetime(2011,9,30));

j1 = find(vicglobal.timevector == datetime(2000,10,1));
j2 = find(vicglobal.timevector == datetime(2011,9,30));

mean_area = 37.5; % km^2
vicglobal.summed_discharge = mean_area*(vicglobal.runoff_sum1 + vicglobal.baseflow_sum1)*1000/(24*3600); % m^3/s
vicglobal_mod.summed_discharge = mean_area*(vicglobal_mod.runoff_sum1 + vicglobal_mod.baseflow_sum1)*1000/(24*3600); % m^3/s
livneh.summed_discharge_cms = livneh.summed_discharge*37.5*1000/3600/24;

% Calculate RMSE and NSE

y = livneh.summed_discharge_cms(i1:i2);
y_hat = vicglobal.summed_discharge(j1:j2)';

rmse = myRMSE(y,y_hat);
nse = myNSE(y,y_hat);

y_hat2 = vicglobal_mod.summed_discharge';

rmse2 = myRMSE(y,y_hat2);
nse2 = myNSE(y,y_hat2);

% Also include naturalized flow data
gage_table = readtable('/Volumes/HD4/SWOTDA/Data/Colorado/gage/naturalized_flow_usgs_09380000.txt');
gage = struct();
gage.time = datetime(str2double(gage_table{:,1}), gage_table{:,2}, 15);
gage.q = gage_table{:,4};
k1 = find(gage.time == datetime(2000,10,15));
k2 = find(gage.time == datetime(2011,9,15));

%%


% % Convert color code to 1-by-3 RGB array (0~1 each)
% str1 = '#a6cee3'; % light blue
% str2 = '#fdae61'; % orange
% str3 = '#1f78b4'; % dark blue % 1f78b4; 2c7bb6
% color1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;
% color2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;
% color3 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255;

lw = 2;
fs = 18;
figure
plot(vicglobal.timevector(j1:j2), vicglobal.summed_discharge(j1:j2), 'linewidth', lw, 'color', [0.3010, 0.7450, 0.9330]) % light blue
title('Upper Colorado discharge')
xlabel('Time') 
ylabel('Q (m^3/s)')
grid on
hold on
plot(vicglobal_mod.timevector, vicglobal_mod.summed_discharge, 'linewidth', lw, 'color', 'magenta') % red
plot(livneh.timevector(i1:i2), livneh.summed_discharge_cms(i1:i2), 'linewidth', lw, 'color', [0.8500, 0.3250, 0.0980]) % orange
plot(gage.time(k1:k2), gage.q(k1:k2), 'linewidth', lw, 'color', 'black') % black
legend('VICGlobal', 'VICGlobal (mod)','Livneh','Gage', 'Location','Best')
set(gca, 'fontsize', fs)

%% Redoing the plot with monthly data

q1 = vicglobal.summed_discharge(j1:j2);
q2 = livneh.summed_discharge_cms(i1:i2);
q3 = vicglobal_mod.summed_discharge;

t1 = vicglobal.timevector(j1:j2);
t2 = livneh.timevector(i1:i2);
t3 = vicglobal_mod.timevector;

% Aggregate to monthly
[t1_m, q1_m] = daily_to_monthly(t1, q1, 'mean');
[t2_m, q2_m] = daily_to_monthly(t2', q2, 'mean');
[t3_m, q3_m] = daily_to_monthly(t3, q3, 'mean');

% Calculate GoF
y = q2_m; % livneh
y_hat = q3_m; % vicglobal

addpath('/Volumes/HD3/SWOTDA/Calibration/CalVIC')
myRMSE(y, y_hat)
myNSE(y, y_hat)

figure
plot(t1_m, q1_m, 'linewidth', lw, 'color', [0.3010, 0.7450, 0.9330]) % light blue
title('Upper Colorado discharge')
xlabel('Time') 
ylabel('Q (m^3/s)')
grid on
hold on
% plot(t3_m, q3_m, 'linewidth', lw, 'color', 'magenta') % red
plot(t2_m, q2_m, 'linewidth', lw, 'color', [0.8500, 0.3250, 0.0980]) % orange
plot(gage.time(k1:k2), gage.q(k1:k2), 'linewidth', lw, 'color', 'black') % black
legend('VICGlobal','L2013','Gage', 'Location','Best')
% legend('VICGlobal', 'VICGlobal (mod)','Livneh','Gage', 'Location','Best')
set(gca, 'fontsize', fs)


%% Investigate soil parameter differences


%% Change VICGlobal soil parameters to match Livneh parameters

soils = load('/Volumes/HD4/SWOTDA/Data/Colorado/colo_soils_vg.txt');
soils_out = soils;
outname = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_soils_vgmod.txt';

soils_out(:,5) = 0.088; % b
soils_out(:,6) = 0.019; % ds
soils_out(:,7) = 8.82; % dsmax (mm/day)
soils_out(:,8) = 0.81; % ws
soils_out(:,25) = 1.4; % layer 3 soil thickness (m)

% soils(1,23:25)

write_soils(5, soils_out, outname, '3l')

%% Run a simulation for UCRB with modified soil parameters

%% Compare simulated SWE

% Extract Livneh SWE time series data over the Upper Colorado Basin

addpath('/Users/jschap/Documents/Research/VICGlobal/Codes/Livneh_Comparison')

start_date = datetime(2000, 10, 1); % must start on the first day of the month
end_date = datetime(2011, 9, 30);
outname = '/Volumes/HD4/SWOTDA/Data/Colorado/L15/swe_sub.mat';
varname = 'SWE';
swe_dir_livneh = '/Volumes/HD3/Livneh_2015/Fluxes';
basin_mask = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif';
[swe_sub_livneh, ~, masklon, masklat] = subset_livneh(start_date, end_date, basin_mask, swe_dir_livneh, outname, varname);

% Livneh subsetting: do it in R. Easier that way. Use --> subset_livneh.R
tifdir = '/Users/jschap/Documents/Research/VICGlobal/Data/uppcololivnehswe/';
tifnames = dir([tifdir '*.tif']);
for tt=1:length(tifnames)
    if tt==1
        [swe_crop, Rswe, swe_lon, swe_lat] = geotiffread2(fullfile(tifdir, tifnames(tt).name));
        [nlat, nlon] = size(swe_crop);
        monthly_swe_crop_livneh = zeros(nlat, nlon, length(tifnames));
    else
        swe_crop = geotiffread2(fullfile(tifdir, tifnames(tt).name));
    end
    monthly_swe_crop_livneh(:,:,tt) = swe_crop;
end

save('/Users/jschap/Documents/Research/VICGlobal/Data/monthly_swe_crop_livneh.mat', 'monthly_swe_crop_livneh', 'nlat', 'nlon', 'tifdir', 'Rswe', 'swe_lon', 'swe_lat')

% % Yeah, each one is (was) the same. That is the problem.
% figure, subplot(2,1,1)
% plotraster(swe_lon, swe_lat, swe_crop, '1', '', '')
% subplot(2,1,2)
% plotraster(swe_lon, swe_lat, swe_crop2, '2', '', '')

% figure, plotraster(swe_lon, swe_lat, swe_crop,'','','')

% Extract VICGlobal SWE data over the Upper Colorado Basin
wb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/wb';
eb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/eb';
swe_col = 27;
[timevector, swe_sub_vg, swe_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, swe_col);

%% Plot basin-average Livneh SWE time series

fs = 20;
figure
plot(start_date:end_date, mean(swe_sub_vg,1), 'linewidth', lw)
title('Upper Colorado SWE')
xlabel('Time')
ylabel('SWE (mm)')
grid on
hold on
plot(start_date:end_date, mean(swe_sub_livneh,1), 'linewidth', lw)
legend('VICGlobal', 'Livneh', 'Location','Best')
set(gca, 'fontsize', fs)

y = mean(swe_sub_livneh,1)';
y_hat = mean(swe_sub_vg,1)';
myRMSE(y,y_hat)

%% Plot time-average SWE over the Upper Colorado basin

% Or potentially average SWE for a particular day of the year?

[mask1, R1, lon1, lat1] = geotiffread2(basin_mask);
nanmask = mask1;
nanmask(nanmask==0) = NaN;

swe_map_vg = xyz2grid(info.lon, info.lat, mean(swe_vg, 2));
swe_map_livneh_sub = flipud(xyz2grid(masklon, masklat, mean(swe_sub_livneh, 2)));

figure, subplot(1,2,1)
plotraster(lon1, lat1, swe_map_vg.*nanmask, 'WY 2000-2011 average SWE (mm)', 'Lon', 'Lat')
subplot(1,2,2)
plotraster(lon1, lat1, swe_map_livneh_sub, 'WY 2000-2011 average SWE (mm)', 'Lon', 'Lat')


%% Plot SWE figure for paper

elevdata = load('/Users/jschap/Documents/Research/VICGlobal/Data/upp_colo_elev.mat');

output_vg = load('/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/processed/vic_outputs_summarized_daily.mat');
output_vg = output_vg.OUTPUTS;

% Basically, I need an image mode counterpart to check_outputs_wrapper
% output_livneh = load('');
% [swe_sub_livneh, masklon, masklat]
% swe_mean_livneh = mean(swe_sub_livneh,2);
% swe_map_livneh = xyz2grid(masklon, masklat, swe_mean_livneh);

% Masking to basin boundary
% [mask1, R1] = geotiffread(basin_mask);
% mask1 = flipud(mask1);
% clipmask = double(mask1);
% clipmask(clipmask==0) = NaN;
% % R1mat = georefobj2mat(R1, 'LL');
% % [masklon1, masklat1] = pixcenters(R1mat, size(mask1));
% % [masklon, masklat, ~] = grid2xyz(masklon1', masklat1', mask1);
% % ncells = length(masklon);
% clipmask = double(mask1);
% clipmask(clipmask==0) = NaN;
% clipped_soil_parameter = cropped_soil_parameter.*clipmask;

%%

nanmask = double(nanmask);
nanmask(nanmask==0) = NaN;

figure
subplot(2,3,1)
plotraster(fliplr(xcoord), fliplr(ycoord), nanmask.*elevdata.cropped_soil_parameter, 'Elevation (m)', 'Lon', 'Lat')

subplot(2,3,2)
swe_map_vg = rot90(nanmask, 2).*output_vg.WB.maps.OUT_SWE;
plotraster(xcoord, ycoord, swe_map_vg, 'VICGlobal SWE (mm)', 'Lon', 'Lat')

% subplot(2,3,3)
% plotraster(fliplr(xcoord), ycoord, flipud(nanmask), 'Basin mask', 'Lon', 'Lat')

subplot(2,3,3)
annual_average_swe_crop_livneh = mean(monthly_swe_crop_livneh, 3);
swe_map_livneh = nanmask.*annual_average_swe_crop_livneh;
plotraster(swe_lon, swe_lat, swe_map_livneh, 'L2013 SWE (mm)', 'Lon', 'Lat')

% figure, plotraster(xcoord, ycoord, annual_average_swe_crop_livneh, '', '','')

subplot(2,3,[4 6])
plot(start_date:end_date, mean(swe_sub_vg,1), 'linewidth', lw)
title('Basin average SWE')
xlabel('Time')
ylabel('SWE (mm)')
grid on
hold on
plot(start_date:end_date, mean(swe_sub_livneh,1), 'linewidth', lw)
legend('VICGlobal', 'Livneh', 'Location','Best')
set(gca, 'fontsize', fs)

% subplot([],1)

%% Load evaporation data for VG and L15

% Extract VICGlobal evaporation data over the Upper Colorado Basin

wb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/wb';
eb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/eb';
evap_col = 4;
[timevector, evap_sub_vg, evap_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, evap_col);

save('/Users/jschap/Documents/Research/VICGlobal/Data/monthly_et_vicglobal.mat', 'evap_sub_vg', 'evap_vg')

vicglobal.evap_sum1 = nansum(evap_sub_vg,1);
sum(vicglobal.evap_sum1(j1:j2))/7833/11 % mm per grid cell and water year

% Livneh subsetting: do it in R. Easier that way. Use --> subset_livneh.R
tifdir = '/Users/jschap/Documents/Research/VICGlobal/Data/uppcololivnehevap/';
tifnames = dir([tifdir '*.tif']);
for tt=1:length(tifnames)
    if tt==1
        [evap_crop, ~, ~, ~] = geotiffread2(fullfile(tifdir, tifnames(tt).name));
        [nlat, nlon] = size(evap_crop);
        monthly_evap_crop_livneh = zeros(nlat, nlon, length(tifnames));
    else
        evap_crop = geotiffread2(fullfile(tifdir, tifnames(tt).name));
    end
    monthly_evap_crop_livneh(:,:,tt) = evap_crop;
end

vicglobal_evap = nanmean(evap_sub_vg,1); % mm
[t_evap, vicglobal_evap_monthly] = daily_to_monthly(timevector, vicglobal_evap, 'mean');

i1 = find(t_evap == datetime(2000,10,15));
i2 = find(t_evap == datetime(2011,9,15));

figure, plotraster(swe_lon, swe_lat, nanmean(monthly_evap_crop_livneh,3), 'ET Livneh','','')

livneh_totalET = squeeze(nanmean(nanmean(monthly_evap_crop_livneh,1),2));

sum(livneh_totalET) % mm/day average for the month, summed over all grid cells, and water years
sum(livneh_totalET)*30/11 % mm/water year, assuming 30 days in a month

%%

figure

subplot(2,2,1)
evap_map_vg = rot90(nanmask, 2).*output_vg.WB.maps.OUT_EVAP;
plotraster(xcoord, ycoord, evap_map_vg, 'VICGlobal ET (mm)', 'Lon', 'Lat')
caxis([0.4,1.2])

% subplot(2,3,3)
% plotraster(fliplr(xcoord), ycoord, flipud(nanmask), 'Basin mask', 'Lon', 'Lat')

subplot(2,2,2)
annual_average_evap_crop_livneh = mean(monthly_evap_crop_livneh, 3);
evap_map_livneh = nanmask.*annual_average_evap_crop_livneh;
plotraster(swe_lon, swe_lat, evap_map_livneh, 'L2013 ET (mm)', 'Lon', 'Lat')
caxis([0.4,2])

subplot(2,2,[3 4])
plot(t_evap(i1:i2), vicglobal_evap_monthly(i1:i2), 'linewidth', lw)
title('Basin average ET')
xlabel('Time')
ylabel('ET (mm)')
grid on
hold on
plot(t_evap(i1:i2), livneh_totalET, 'linewidth', lw)
legend('VICGlobal ET', 'L2013 ET', 'Location','Best')
set(gca, 'fontsize', fs)

%% Load soil moisture data for VG and L15

% Extract VICGlobal evaporation data over the Upper Colorado Basin
sm_1_col = 21;
[timevector, sm1_sub_vg, sm1_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, sm_1_col);
sm_2_col = 22;
[timevector, sm3_sub_vg, sm3_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, sm_2_col);
sm_3_col = 23;
[timevector, sm3_sub_vg, sm3_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, sm_3_col);

save('/Users/jschap/Documents/Research/VICGlobal/Data/owtrta_backup_02032020.mat')
% load('/Users/jschap/Documents/Research/VICGlobal/Data/owtrta_backup_02032020.mat')

% Livneh subsetting: do it in R. Easier that way. Use --> subset_livneh.R
tifdir = '/Users/jschap/Documents/Research/VICGlobal/Data/uppcololivneh_soil_moisture_1/';
tifnames = dir([tifdir '*.tif']);
for tt=1:length(tifnames)
    if tt==1
        monthly_sm1_crop_livneh = zeros(nlat, nlon, length(tifnames));
    end
    sm1_crop = geotiffread2(fullfile(tifdir, tifnames(tt).name));
    monthly_sm1_crop_livneh(:,:,tt) = sm1_crop;
end

vicglobal_sm1 = nanmean(sm1_sub_vg,1); % mm
[t_sm, vicglobal_sm1_monthly] = daily_to_monthly(timevector, vicglobal_sm1, 'mean');

i1 = find(t_sm == datetime(2000,10,15));
i2 = find(t_sm == datetime(2011,9,15));

figure, plotraster(swe_lon, swe_lat, nanmean(monthly_sm1_crop_livneh,3), 'SM1 Livneh','','')

livneh_sm1 = squeeze(nanmean(nanmean(monthly_sm1_crop_livneh,1),2));

%%

figure

subplot(2,2,1)
evap_map_vg = rot90(nanmask, 2).*output_vg.WB.maps.OUT_SOIL_MOIST_0;
plotraster(xcoord, ycoord, evap_map_vg, 'VG Soil Moisture 1 (mm)', 'Lon', 'Lat')
% caxis([0,115])

subplot(2,2,2)
annual_average_sm1_crop_livneh = mean(monthly_sm1_crop_livneh, 3);
sm1_map_livneh = nanmask.*annual_average_sm1_crop_livneh;
plotraster(swe_lon, swe_lat, sm1_map_livneh, 'L15 Soil Moisture 1 (mm)', 'Lon', 'Lat')
% caxis([0,115])

subplot(2,2,[3 4])
plot(t_evap(i1:i2), vicglobal_sm1_monthly(i1:i2), 'linewidth', lw)
title('UCRB SM1 comparison')
xlabel('Time')
ylabel('Units (mm)')
grid on
hold on
plot(t_evap(i1:i2), livneh_sm1, 'linewidth', lw)
legend('VICGlobal SM1', 'Livneh SM1', 'Location','Best')
set(gca, 'fontsize', fs)

%% Layer 2 soil moisture

tifdir = '/Users/jschap/Documents/Research/VICGlobal/Data/uppcololivneh_soil_moisture_2/';
tifnames = dir([tifdir '*.tif']);
for tt=1:length(tifnames)
    if tt==1
        monthly_sm3_crop_livneh = zeros(nlat, nlon, length(tifnames));
    end
    sm3_crop = geotiffread2(fullfile(tifdir, tifnames(tt).name));
    monthly_sm3_crop_livneh(:,:,tt) = sm3_crop;
end

vicglobal_sm3 = nanmean(sm3_sub_vg,1); % mm
[t_sm, vicglobal_sm3_monthly] = daily_to_monthly(timevector, vicglobal_sm3, 'mean');

i1 = find(t_sm == datetime(2000,10,15));
i2 = find(t_sm == datetime(2011,9,15));

figure, plotraster(swe_lon, swe_lat, nanmean(monthly_sm3_crop_livneh,3), 'sm3 Livneh','','')

livneh_sm3 = squeeze(nanmean(nanmean(monthly_sm3_crop_livneh,1),2));

figure

subplot(2,2,1)
sm3_map_vg = rot90(nanmask, 2).*output_vg.WB.maps.OUT_SOIL_MOIST_1;
plotraster(xcoord, ycoord, sm3_map_vg, 'VG Soil Moisture 2 (mm)', 'Lon', 'Lat')
% caxis([0,115])

subplot(2,2,2)
annual_average_sm3_crop_livneh = mean(monthly_sm3_crop_livneh, 3);
sm3_map_livneh = nanmask.*annual_average_sm3_crop_livneh;
plotraster(swe_lon, swe_lat, sm3_map_livneh, 'L15 Soil Moisture 2 (mm)', 'Lon', 'Lat')
% caxis([0,115])

subplot(2,2,[3 4])
plot(t_evap(i1:i2), vicglobal_sm3_monthly(i1:i2), 'linewidth', lw)
title('UCRB sm3 comparison')
xlabel('Time')
ylabel('Units (mm)')
grid on
hold on
plot(t_evap(i1:i2), livneh_sm3, 'linewidth', lw)
legend('VICGlobal sm3', 'Livneh sm3', 'Location','Best')
set(gca, 'fontsize', fs)

%% Layer 3 soil moisture

tifdir = '/Users/jschap/Documents/Research/VICGlobal/Data/uppcololivneh_soil_moisture_3/';
tifnames = dir([tifdir '*.tif']);
for tt=1:length(tifnames)
    if tt==1
        monthly_sm3_crop_livneh = zeros(nlat, nlon, length(tifnames));
    end
    sm3_crop = geotiffread2(fullfile(tifdir, tifnames(tt).name));
    monthly_sm3_crop_livneh(:,:,tt) = sm3_crop;
end

vicglobal_sm3 = nanmean(sm3_sub_vg,1); % mm
[t_sm, vicglobal_sm3_monthly] = daily_to_monthly(timevector, vicglobal_sm3, 'mean');

i1 = find(t_sm == datetime(2000,10,15));
i2 = find(t_sm == datetime(2011,9,15));

figure, plotraster(swe_lon, swe_lat, nanmean(monthly_sm3_crop_livneh,3), 'SM3 Livneh','','')

livneh_sm3 = squeeze(nanmean(nanmean(monthly_sm3_crop_livneh,1),2));

figure

subplot(2,2,1)
sm3_map_vg = rot90(nanmask, 2).*output_vg.WB.maps.OUT_SOIL_MOIST_2;
plotraster(xcoord, ycoord, sm3_map_vg, 'VG Soil Moisture 3 (mm)', 'Lon', 'Lat')
% caxis([0,115])

subplot(2,2,2)
annual_average_sm3_crop_livneh = mean(monthly_sm3_crop_livneh, 3);
sm3_map_livneh = nanmask.*annual_average_sm3_crop_livneh;
plotraster(swe_lon, swe_lat, sm3_map_livneh, 'L15 Soil Moisture 3 (mm)', 'Lon', 'Lat')
% caxis([0,115])

subplot(2,2,[3 4])
plot(t_evap(i1:i2), vicglobal_sm3_monthly(i1:i2), 'linewidth', lw)
title('UCRB SM3 comparison')
xlabel('Time')
ylabel('Units (mm)')
grid on
hold on
plot(t_evap(i1:i2), livneh_sm3, 'linewidth', lw)
legend('VICGlobal SM3', 'Livneh SM3', 'Location','Best')
set(gca, 'fontsize', fs)

%% Water balance calculations

livneh.runoff_sum1 = sum(livneh.runoff_sub, 2);
livneh.baseflow_sum1 = sum(livneh.baseflow_sub, 2);

vicglobal.runoff_sum1 % mm
vicglobal.baseflow_sum1 % mm

i1 = find(livneh.timevector == datetime(2000,10,1));
i2 = find(livneh.timevector == datetime(2011,9,30));

j1 = find(vicglobal.timevector == datetime(2000,10,1));
j2 = find(vicglobal.timevector == datetime(2011,9,30));

figure, subplot(2,1,1)
plot(vicglobal.timevector(j1:j2), vicglobal.runoff_sum1(j1:j2), 'linewidth', lw)
hold on
plot(livneh.timevector(i1:i2), livneh.runoff_sum1(i1:i2), 'linewidth', lw)
legend('VICGlobal','Livneh')
ylabel('Runoff (mm)')
xlabel('Time')
set(gca, 'fontsize', fs)

subplot(2,1,2)
plot(vicglobal.timevector(j1:j2), vicglobal.baseflow_sum1(j1:j2), 'linewidth', lw)
hold on
plot(livneh.timevector(i1:i2), livneh.baseflow_sum1(i1:i2), 'linewidth', lw)
legend('VICGlobal','Livneh')
ylabel('Baseflow (mm)')
xlabel('Time')
set(gca, 'fontsize', fs)

% Over the WY2001-2011 time period, the total runoff and baseflow are:
sum(vicglobal.baseflow_sum1(j1:j2)) % mm per all grid cells and times
sum(vicglobal.baseflow_sum1(j1:j2))/7833 % mm per grid cell over all times
(sum(vicglobal.baseflow_sum1(j1:j2))/7833)/10 % mm per grid cell and water year

(sum(vicglobal.runoff_sum1(j1:j2))/7833)/10 % mm per grid cell and water year

(sum(livneh.runoff_sum1(i1:i2))/7833)/10 % mm per grid cell and water year
(sum(livneh.baseflow_sum1(i1:i2))/7833)/10 % mm per grid cell and water year

sum(vicglobal_evap)/(7833*10)

%% Precipitation

output_vg.WB.maps.OUT_PREC; % time average precipiation over 1980-2011 simulaiton
nanmean(output_vg.WB.maps.OUT_PREC(:))*365 % mm/year on average

%% Load precipitation data for VG and L15

% Extract VICGlobal PPT data over the Upper Colorado Basin
addpath('/Users/jschap/Documents/Research/VICGlobal/Codes')
addpath(genpath('/Users/jschap/Documents/Codes/VICMATLAB'))
wb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/wb';
eb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/eb';
ppt_col = 13;
[timevector, ppt_sub_vg, ppt_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, ppt_col);

vicglobal.ppt1_sum1 = nansum(ppt_sub_vg,1);
sum(vicglobal.ppt1_sum1(j1:j2))/7833/11 % mm per grid cell and water year

% Livneh PPT

info_nc = ncinfo('/Volumes/HD3/Livneh_2013/MetNC/prec.mon.mean.nc');
ppt_livneh = ncread('/Volumes/HD3/Livneh_2013/MetNC/prec.mon.mean.nc', 'prec');
lon_livneh = ncread('/Volumes/HD3/Livneh_2013/MetNC/prec.mon.mean.nc', 'lon');
lat_livneh = ncread('/Volumes/HD3/Livneh_2013/MetNC/prec.mon.mean.nc', 'lat');
size(ppt_livneh)
[mask1, R1, mask_lon, mask_lat] = geotiffread2(basin_mask);
% Presumably it is the same as VICGlobal PPT, unless there was some other
% modification that I was unaware of. Another possibility is that there was
% a difference introduced when I downscaled the PPT data using
% MT-CLIM/VIC-4 prior to running the model in VIC-5.

%%

figure

subplot(2,2,1)
evap_map_vg = rot90(nanmask, 2).*output_vg.WB.maps.OUT_EVAP;
plotraster(xcoord, ycoord, evap_map_vg, 'VICGlobal Evap (mm)', 'Lon', 'Lat')
% caxis([0.4,2])

% subplot(2,3,3)
% plotraster(fliplr(xcoord), ycoord, flipud(nanmask), 'Basin mask', 'Lon', 'Lat')

subplot(2,2,2)
annual_average_evap_crop_livneh = mean(monthly_evap_crop_livneh, 3);
evap_map_livneh = nanmask.*annual_average_evap_crop_livneh;
plotraster(swe_lon, swe_lat, evap_map_livneh, 'Livneh Total ET (mm)', 'Lon', 'Lat')
% caxis([0.4,2])

subplot(2,2,[3 4])
plot(t_evap(i1:i2), vicglobal_evap_monthly(i1:i2), 'linewidth', lw)
title('UCRB ET comparison')
xlabel('Time')
ylabel('Units (mm)')
grid on
hold on
plot(t_evap(i1:i2), livneh_totalET, 'linewidth', lw)
legend('VICGlobal Evaporation', 'Livneh Total ET', 'Location','Best')
set(gca, 'fontsize', fs)

%% Load canopy evaporation data

% Extract VICGlobal canopy evaporation data over the Upper Colorado Basin
addpath('/Users/jschap/Documents/Research/VICGlobal/Codes')
wb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/wb';
eb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/eb';
evapcan_col = 7;
[timevector, evapcan_sub_vg, evapcan_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, evapcan_col);

save('/Users/jschap/Documents/Research/VICGlobal/Data/monthly_evapcan_vicglobal.mat', 'evapcan_sub_vg', 'evapcan_vg')

vicglobal.evapcan_sum1 = nansum(evapcan_sub_vg,1);
sum(vicglobal.evapcan_sum1(j1:j2))/7833/11 % mm per grid cell and water year

% % Livneh subsetting: do it in R. Easier that way. Use --> subset_livneh.R
% tifdir = '/Users/jschap/Documents/Research/VICGlobal/Data/uppcololivnehevap/';
% tifnames = dir([tifdir '*.tif']);
% for tt=1:length(tifnames)
%     if tt==1
%         [evap_crop, ~, ~, ~] = geotiffread2(fullfile(tifdir, tifnames(tt).name));
%         [nlat, nlon] = size(evap_crop);
%         monthly_evap_crop_livneh = zeros(nlat, nlon, length(tifnames));
%     else
%         evap_crop = geotiffread2(fullfile(tifdir, tifnames(tt).name));
%     end
%     monthly_evap_crop_livneh(:,:,tt) = evap_crop;
% end

vicglobal_evapcan = nanmean(evapcan_sub_vg,1); % mm
[t_evap, vicglobal_evapcan_monthly] = daily_to_monthly(timevector, vicglobal_evapcan, 'mean');

i1 = find(t_evap == datetime(2000,10,15));
i2 = find(t_evap == datetime(2011,9,15));

% figure, plotraster(swe_lon, swe_lat, nanmean(monthly_evap_crop_livneh,3), 'ET Livneh','','')

% livneh_totalET = squeeze(nanmean(nanmean(monthly_evap_crop_livneh,1),2));

% sum(livneh_totalET) % mm/day average for the month, summed over all grid cells, and water years
% sum(livneh_totalET)*30/11 % mm/water year, assuming 30 days in a month

%% Plot canopy evaporation data

figure

subplot(2,2,1)
evapcan_map_vg = rot90(nanmask, 2).*output_vg.WB.maps.OUT_EVAP_CANOP;
plotraster(xcoord, ycoord, evapcan_map_vg, 'VG Canopy Evaporation (mm)', 'Lon', 'Lat')
% caxis([0,115])

subplot(2,2,2)
% annual_average_sm3_crop_livneh = mean(monthly_sm3_crop_livneh, 3);
% sm3_map_livneh = nanmask.*annual_average_sm3_crop_livneh;
% plotraster(swe_lon, swe_lat, sm3_map_livneh, 'L15 Soil Moisture 2 (mm)', 'Lon', 'Lat')
% caxis([0,115])

subplot(2,2,[3 4])
plot(t_evap(i1:i2), vicglobal_evapcan_monthly(i1:i2), 'linewidth', lw)
title('UCRB Evap_{can} comparison')
xlabel('Time')
ylabel('Units (mm)')
grid on
hold on
% plot(t_evap(i1:i2), livneh_sm3, 'linewidth', lw)
% legend('VICGlobal sm3', 'Livneh sm3', 'Location','Best')
set(gca, 'fontsize', fs)

%% Load PET data

% Extract VICGlobal PET data over the Upper Colorado Basin
addpath('/Users/jschap/Documents/Research/VICGlobal/Codes')
wb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/wb';
eb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/eb';
pet_col = 12;
[timevector, pet_sub_vg, pet_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, pet_col);

save('/Users/jschap/Documents/Research/VICGlobal/Data/monthly_pet_vicglobal.mat', 'pet_sub_vg', 'pet_vg')

vicglobal.pet_sum1 = nansum(pet_sub_vg,1);
sum(vicglobal.pet_sum1(j1:j2))/7833/11 % mm per grid cell and water year

% Livneh subsetting: do it in R. Easier that way. Use --> subset_livneh.R
tifdir = '/Users/jschap/Documents/Research/VICGlobal/Data/uppcololivneh_pet_short/';
tifnames = dir([tifdir '*.tif']);
for tt=1:length(tifnames)
    if tt==1
        [pet_crop, ~, ~, ~] = geotiffread2(fullfile(tifdir, tifnames(tt).name));
        [nlat, nlon] = size(pet_crop);
        monthly_pet_crop_livneh = zeros(nlat, nlon, length(tifnames));
    else
        pet_crop = geotiffread2(fullfile(tifdir, tifnames(tt).name));
    end
    monthly_pet_crop_livneh(:,:,tt) = pet_crop;
end

vicglobal_pet = nanmean(pet_sub_vg,1); % mm
[t_evap, vicglobal_pet_monthly] = daily_to_monthly(timevector, vicglobal_pet, 'mean');

i1 = find(t_evap == datetime(2000,10,15));
i2 = find(t_evap == datetime(2011,9,15));

% figure, plotraster(swe_lon, swe_lat, nanmean(monthly_evap_crop_livneh,3), 'ET Livneh','','')

livneh_pet = squeeze(nanmean(nanmean(monthly_pet_crop_livneh,1),2));

sum(livneh_pet) % mm/day average for the month, summed over all grid cells, and water years
sum(livneh_pet)*30/11 % mm/water year, assuming 30 days in a month

%% Plot PET data

figure

subplot(2,2,1)
pet_map_vg = rot90(nanmask, 2).*output_vg.WB.maps.OUT_PET;
plotraster(xcoord, ycoord, pet_map_vg, 'VG PET (mm)', 'Lon', 'Lat')
caxis([1.5,4])

subplot(2,2,2)
annual_average_pet_crop_livneh = mean(monthly_pet_crop_livneh, 3);
pet_map_livneh = nanmask.*annual_average_pet_crop_livneh;
plotraster(swe_lon, swe_lat, pet_map_livneh, 'L2013 PET (mm)', 'Lon', 'Lat')
caxis([1.5,4])

subplot(2,2,[3 4])
plot(t_evap(i1:i2), vicglobal_pet_monthly(i1:i2), 'linewidth', lw)
title('UCRB PET comparison')
xlabel('Time')
ylabel('Units (mm)')
grid on
hold on
plot(t_evap(i1:i2), livneh_pet, 'linewidth', lw)
legend('VICGlobal PET', 'L2013 PET', 'Location','Best')
set(gca, 'fontsize', fs)



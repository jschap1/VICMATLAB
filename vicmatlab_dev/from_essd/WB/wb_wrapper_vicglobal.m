% Water balance wrapper
%
% This wrapper is for analyzing the water balance from VIC modeling
% Written 2/26/2020 JRS
%
% Currently in use for the Upper Colorado River Basin, with the VICGlobal
% results

%% Setup

addpath('/Users/jschap/Documents/Research/VICGlobal/Codes/Livneh_Comparison')
addpath('/Users/jschap/Documents/Research/VICGlobal/Codes/WB')
addpath('/Users/jschap/Documents/Codes/VICMATLAB/Utility')
addpath('/Users/jschap/Documents/Research/VICGlobal/Codes')
addpath(genpath('/Users/jschap/Documents/Research/VICMATLAB/'))
addpath('/Users/jschap/Documents/Codes/VICMATLAB/Process_Outputs_1')
addpath('/Users/jschap/Documents/Codes/VICMATLAB/v5')

start_date = datetime(2000, 10, 1);
end_date = datetime(2011, 9, 30);
timevector_daily = (start_date:end_date)';
[timevector_monthly, ~] = daily_to_monthly(timevector_daily, ones(size(timevector_daily)), 'mean');

basin_mask = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif';
[mask1, Rmask, masklon, masklat] = geotiffread2(basin_mask);
nanmask = double(mask1);
nanmask(nanmask==0) = NaN;
[nlon, nlat] = size(mask1);

wb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/wb';
eb_dir_vicglobal = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/eb';
saveloc = '/Users/jschap/Documents/Research/VICGlobal/Data/vicglobal_results';

%% SWE

% indir = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/processed';
% dat = load(fullfile(indir, 'vic_outputs_summarized_daily.mat'));

% Extract VICGlobal SWE data over the Upper Colorado Basin (takes about 40
% minutes to run on Atlantic). Best practice is to run in a separate Matlab
% window to get some crude parallelization.
swe_col = 27;
[timevector, swe_sub_vg, swe_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, swe_col);
save(fullfile(saveloc, 'swe_vg.mat'), 'timevector', 'swe_sub_vg', 'swe_vg', 'info');

%% Soil moisture

t1 = 0.1; % m, depth of soil layer 1
t2 = 0.2; % m, depth of soil layer 2 
t3 = 0.7; % m, depth of soil layer 3

sm_1_col = 21;
[timevector, sm1_sub_vg, sm1_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, sm_1_col);
save(fullfile(saveloc, 'sm1_vg.mat'), 'timevector', 'sm1_sub_vg', 'sm1_vg', 'info');
sm_2_col = 22;
[timevector, sm2_sub_vg, sm2_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, sm_2_col);
save(fullfile(saveloc, 'sm2_vg.mat'), 'timevector', 'sm2_sub_vg', 'sm2_vg', 'info');
sm_3_col = 23;
[timevector, sm3_sub_vg, sm3_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, sm_3_col);
save(fullfile(saveloc, 'sm3_vg.mat'), 'timevector', 'sm3_sub_vg', 'sm3_vg', 'info');

%% ET

evap_col = 4;
[timevector, evap_sub_vg, evap_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, evap_col);
save(fullfile(saveloc, 'evap_vg.mat'), 'timevector', 'evap_sub_vg', 'evap_vg', 'info');

%% Precipitation

ppt_col = 13;
[timevector, ppt_sub_vg, ppt_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, ppt_col);
save(fullfile(saveloc, 'ppt_vg.mat'), 'timevector', 'ppt_sub_vg', 'ppt_vg', 'info');

%% Runoff

runoff_col = 5;
[timevector, runoff_sub_vg, runoff_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, runoff_col);
save(fullfile(saveloc, 'runoff_vg.mat'), 'timevector', 'runoff_sub_vg', 'runoff_vg', 'info');

%% Baseflow

baseflow_col = 6;
[timevector, baseflow_sub_vg, baseflow_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, baseflow_col);
save(fullfile(saveloc, 'baseflow_vg.mat'), 'timevector', 'baseflow_sub_vg', 'baseflow_vg', 'info');

%% Canopy moisture

wdew_col = 17;
[timevector, wdew_sub_vg, wdew_vg, info] = load_vic_output(wb_dir_vicglobal, eb_dir_vicglobal, basin_mask, wdew_col);
save(fullfile(saveloc, 'wdew_vg.mat'), 'timevector', 'wdew_sub_vg', 'wdew_vg', 'info');

%% Make water balance object

out_wdew = load(fullfile(saveloc, 'wdew_vg.mat'));
timevector = out_wdew.timevector;

j1 = find(timevector == datetime(2000,10,1));
j2 = find(timevector == datetime(2011,9,30));

% figure, plot(out1.wdew_sub_vg(:,j1:j2))

WB = struct();
WB.daily.time = timevector(j1:j2);
WB.daily.WDEW = mean(out_wdew.wdew_sub_vg(:,j1:j2))';
[WB.monthly.time, WB.monthly.WDEW] = daily_to_monthly(WB.daily.time, WB.daily.WDEW, 'mean');

figure, plot(WB.daily.time, WB.daily.WDEW);
hold on
plot(WB.monthly.time, WB.monthly.WDEW);

out_et = load(fullfile(saveloc, 'evap_vg.mat'));
WB.daily.ET = nanmean(out_et.evap_sub_vg(:,j1:j2))';
[~, WB.monthly.ET] = daily_to_monthly(WB.daily.time, WB.daily.ET, 'sum');

out_runoff = load(fullfile(saveloc, 'runoff_vg.mat'));
WB.daily.RUNOFF = nanmean(out_runoff.runoff_sub_vg(:,j1:j2))';
[~, WB.monthly.RUNOFF] = daily_to_monthly(WB.daily.time, WB.daily.RUNOFF, 'sum');

out_baseflow = load(fullfile(saveloc, 'baseflow_vg.mat'));
WB.daily.BASEFLOW = nanmean(out_baseflow.baseflow_sub_vg(:,j1:j2))';
[~, WB.monthly.BASEFLOW] = daily_to_monthly(WB.daily.time, WB.daily.BASEFLOW, 'sum');

out_precip = load(fullfile(saveloc, 'ppt_vg.mat'));
WB.daily.PRECIP = nanmean(out_precip.ppt_sub_vg(:,j1:j2))';
[~, WB.monthly.PRECIP] = daily_to_monthly(WB.daily.time, WB.daily.PRECIP, 'sum');

out_swe = load(fullfile(saveloc, 'swe_vg.mat'));
WB.daily.SWE = nanmean(out_swe.swe_sub_vg(:,j1:j2))';
[~, WB.monthly.SWE] = daily_to_monthly(WB.daily.time, WB.daily.SWE, 'mean');

out_sm1 = load(fullfile(saveloc, 'sm1_vg.mat'));
WB.daily.MOIST1 = nanmean(out_sm1.sm1_sub_vg(:,j1:j2))';
[~, WB.monthly.MOIST1] = daily_to_monthly(WB.daily.time, WB.daily.MOIST1, 'mean');

out_sm2 = load(fullfile(saveloc, 'sm2_vg.mat'));
WB.daily.MOIST2 = nanmean(out_sm2.sm2_sub_vg(:,j1:j2))';
[~, WB.monthly.MOIST2] = daily_to_monthly(WB.daily.time, WB.daily.MOIST2, 'mean');

out_sm3 = load(fullfile(saveloc, 'sm3_vg.mat'));
WB.daily.MOIST3 = nanmean(out_sm3.sm3_sub_vg(:,j1:j2))';
[~, WB.monthly.MOIST3] = daily_to_monthly(WB.daily.time, WB.daily.MOIST3, 'mean');

% calculate total storage
WB.monthly.S = WB.monthly.SWE + WB.monthly.WDEW + WB.monthly.MOIST1 + WB.monthly.MOIST2 + WB.monthly.MOIST3;
nmonth = length(WB.monthly.S);

% calculate precipitation as residual
WB.monthly.P_resid = zeros(nmonth, 1);
for t=2:nmonth
    WB.monthly.P_resid(t) = WB.monthly.ET(t) + WB.monthly.BASEFLOW(t) + WB.monthly.RUNOFF(t) + WB.monthly.S(t) - WB.monthly.S(t-1);
end

figure, plot(WB.monthly.time, WB.monthly.P_resid), hold on,
plot(WB.monthly.time, WB.monthly.PRECIP);

%% Calculate values for water balance

mean(WB.monthly.S)

%% Calculate water-year average change in storage

for i=1:11
   storage_change(i) = WB.monthly.S(12*i) - WB.monthly.S(12*(i-1)+1);
end
 
%% Plot water balance time series (subplots)
% 
% figure
% lw = 1.5;
% fs = 18;
% 
% subplot(4,1,1)
% plot(WB.monthly.time, [0; diff(WB.monthly.S)], 'linewidth', lw)
% grid on
% title('\Delta S (mm)')
% set(gca, 'fontsize', fs)
% xlabel('Time (months)')
% ylabel('\Delta S (mm)')
% 
% subplot(4,1,2)
% plot(WB.monthly.time, WB.monthly.ET, 'linewidth', lw)
% grid on
% title('ET (mm)')
% set(gca, 'fontsize', fs)
% xlabel('Time (months)')
% ylabel('ET (mm)')
% 
% subplot(4,1,3)
% plot(WB.monthly.time, WB.monthly.BASEFLOW + WB.monthly.RUNOFF, 'linewidth', lw)
% grid on
% title('Q (mm)')
% set(gca, 'fontsize', fs)
% xlabel('Time (months)')
% ylabel('Q (mm)')
% 
% subplot(4,1,4)
% plot(WB.monthly.time, WB.monthly.P_resid, 'linewidth', lw)
% % plot(WB.monthly.time, WB.monthly.PRECIP, 'linewidth', lw)
% grid on
% title('P (mm)')
% set(gca, 'fontsize', fs)
% xlabel('Time (months)')
% ylabel('P (mm)')
% 
% %% Plot water balance time series (one figure)
% 
% figure
% plot(WB.time, WB.P, 'linewidth', lw)
% hold on
% plot(WB.time, WB.BASEFLOW + WB.RUNOFF, 'linewidth', lw)
% plot(WB.time, WB.ET, 'linewidth', lw)
% plot(WB.time, [0; diff(WB.S)], 'linewidth', lw)
% legend('P', 'Q','ET','\Delta S','Location','Best')
% xlabel('Time (months)')
% ylabel('Value (mm)')
% title('Upper Colorado Water Balance')
% grid on
% set(gca, 'fontsize', fs)
% 
% 
% 
% %% Calculate water-year averages
% 
% WB.WY = struct();
% [WB.WY.time, WB.WY.P] = monthly_to_annual(WB.time, WB.P, 'sum');
% [~, WB.WY.S] = monthly_to_annual(WB.time, WB.S, 'sum');
% [~, WB.WY.Q] = monthly_to_annual(WB.time, WB.RUNOFF + WB.BASEFLOW, 'sum');
% [~, WB.WY.ET] = monthly_to_annual(WB.time, WB.ET, 'sum');
% 
% figure
% plot(WB.WY.time, WB.WY.P, 'linewidth', lw)
% hold on
% plot(WB.WY.time, WB.WY.Q, 'linewidth', lw)
% plot(WB.WY.time, WB.WY.ET, 'linewidth', lw)
% plot(WB.WY.time, [0; diff(WB.WY.S)], 'linewidth', lw)
% legend('P', 'Q','ET','\Delta S','Location','Best')
% xlabel('Time (months)')
% ylabel('Value (mm)')
% title('Upper Colorado Water Balance')
% grid on
% set(gca, 'fontsize', fs)

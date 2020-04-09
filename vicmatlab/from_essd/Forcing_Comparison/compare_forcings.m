% Comparing the different iterations of the Livneh forcing data to ensure
% that the correct version is being used for the VICGlobal validation
%
% Written 2/27/2020 JRS
%
% A check so we are sure that we are comparing VIC simulations that are
% consistent with one another. 
 
%% Header/setup

addpath('/Users/jschap/Documents/Codes/VICMATLAB/Utility')
addpath('/Users/jschap/Documents/Research/VICGlobal/Codes/Forcing_Comparison')

% choose a grid cell in the Upper Colorado River Basin
target_lat = 35.59375;
target_lon = -106.84375;

% choose a range of dates to look at
start_date = datetime(2000, 10, 1);
end_date = datetime(2011, 9, 30);
time.day = (start_date:end_date)';
[time.month, ~] = daily_to_monthly(time.day, ones(size(time.day)), 'mean');
[time.year, ~] = monthly_to_annual(time.month, ones(size(time.month)), 'mean', 0);

nday = length(time.day);
nmonth = length(time.month);
nyear = length(time.year);
res = 1/16;

%% A) Livneh et al. (2013) Metnc folder

% Adapting code from VICMATLAB/v4/vicinputworkflow.m

forcingpath = '/Volumes/HD3/Livneh_2013/MetNC';
forc_A = process_L13_forcing_data_A(forcingpath, start_date, end_date, target_lon, target_lat);

%%

%%

%%

%% F) ASCII forcing data converted from (A) Livneh et al. (2013)

forcingpath = 1;
forc_B = 1;

%%

%% Plot each variable

figure
lw = 1.5;
fs = 18;
legendtext = ['L13 Metnc (A)'];

subplot(2,2,1)
plot(forc_A.time_daily, forc_A.prec, 'linewidth', lw)
title('Precipitation'), xlabel('Time')
ylabel('PPT (mm/day)'), grid on, hold on
legend(legendtext)
set(gca, 'fontsize', fs)

subplot(2,2,2)
plot(forc_A.time_daily, forc_A.tmax, 'linewidth', lw)
title('Tmax'), xlabel('Time')
ylabel('Tmax (deg. C)'), grid on, hold on
legend(legendtext)
set(gca, 'fontsize', fs)

subplot(2,2,3)
plot(forc_A.time_daily, forc_A.tmin, 'linewidth', lw)
title('Tmin'), xlabel('Time')
ylabel('Tmin (deg. C)'), grid on, hold on
legend(legendtext)
set(gca, 'fontsize', fs)

subplot(2,2,4)
plot(forc_A.time_daily, forc_A.wind, 'linewidth', lw)
title('Wind'), xlabel('Time')
ylabel('Wind (m/s)'), grid on, hold on
legend(legendtext)
set(gca, 'fontsize', fs)


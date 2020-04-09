% Water balance wrapper
%
% This wrapper is for analyzing the water balance from VIC modeling
% Written 2/24/2020 JRS
%
% Currently in use for the Upper Colorado River Basin, with the Livneh et
% al. (2013) results. 

%% Setup

addpath('/Users/jschap/Documents/Research/VICGlobal/Codes/Livneh_Comparison')
addpath('/Users/jschap/Documents/Research/VICGlobal/Codes/WB')

start_date = datetime(2000, 10, 1);
end_date = datetime(2011, 9, 30);
timevector_daily = (start_date:end_date)';
[timevector_monthly, ~] = daily_to_monthly(timevector_daily, ones(size(timevector_daily)), 'mean');

basin_mask = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif';
[mask1, Rmask, masklon, masklat] = geotiffread2(basin_mask);
nanmask = double(mask1);
nanmask(nanmask==0) = NaN;
[nlon, nlat] = size(mask1);

%% SWE

tifdir_swe = '/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/SWE/';
SWE = nanmask.*read_var_from_geotiffs(tifdir_swe, [nlon, nlat]);

figure, plotraster(masklon, masklat, SWE(:,:,103), 'April 2009 average SWE (mm)', '', '')
% timevector_monthly(103)

%% Soil moisture

t1 = 0.1; % m, depth of soil layer 1
t2 = 0.2; % m, depth of soil layer 2 
t3 = 0.7; % m, depth of soil layer 3

% Subset the Livneh et al. (2013) results to the basin. 
% Use subset_livneh.R.
% Results are saved in tifdir. 

tifdir1 = '/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/SM1';
tifdir2 = '/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/SM2';
tifdir3 = '/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/SM3';

SM.layer1 = nanmask.*read_var_from_geotiffs(tifdir1, [nlon, nlat]);
SM.layer2 = nanmask.*read_var_from_geotiffs(tifdir2, [nlon, nlat]);
SM.layer3 = nanmask.*read_var_from_geotiffs(tifdir3, [nlon, nlat]);

%% ET

% The L13 fluxes are reported in mm/s. Thus, a conversion is needed to
% mm/day.

tifdir_et = '/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/ET';
ET = nanmask.*read_var_from_geotiffs(tifdir_et, [nlon, nlat]);

ndays_in_month = eomday(year(timevector_monthly),month(timevector_monthly));
nmonth = size(ET,3);
for t=1:nmonth
    ET(:,:,t) = ndays_in_month(t).*ET(:,:,t);
end
ET = ET.*3600*24; % unit conversion

figure, plotraster(masklon, masklat, ET(:,:,103), 'April 2009 average ET (mm)', '', '')

%% Runoff

tifdir_runoff = '/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/runoff';
RUNOFF = nanmask.*read_var_from_geotiffs(tifdir_runoff, [nlon, nlat]);
for t=1:nmonth
    RUNOFF(:,:,t) = ndays_in_month(t).*RUNOFF(:,:,t);
end
RUNOFF = RUNOFF.*3600*24; % unit conversion 

figure, plotraster(masklon, masklat, RUNOFF(:,:,103), 'April 2009 average runoff (mm)', '', '')

%% Baseflow

tifdir_baseflow = '/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/baseflow';
BASEFLOW = nanmask.*read_var_from_geotiffs(tifdir_baseflow, [nlon, nlat]);
for t=1:nmonth
    BASEFLOW(:,:,t) = ndays_in_month(t).*BASEFLOW(:,:,t);
end
BASEFLOW = BASEFLOW.*3600*24; % unit conversion

figure, plotraster(masklon, masklat, BASEFLOW(:,:,103), 'April 2009 average baseflow (mm)', '', '')

%% Canopy moisture

tifdir_w = '/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/wdew';
WDEW = nanmask.*read_var_from_geotiffs(tifdir_w, [nlon, nlat]);

figure, plotraster(masklon, masklat, WDEW(:,:,103), 'April 2009 average canopy moisture (mm)', '', '')

%% Precipitation

tifdir_prec = '/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/PREC';
PREC = nanmask.*read_var_from_geotiffs(tifdir_prec, [nlon, nlat]);
for t=1:nmonth
    PREC(:,:,t) = ndays_in_month(t).*PREC(:,:,t);
end
PREC = PREC.*3600*24; % unit conversion

figure, plotraster(masklon, masklat, PREC(:,:,103), 'April 2009 average precipitation (mm)', '', '')

%% Plot storage

% STOR = struct();
% STOR.map = SWE + SM.layer1 + SM.layer2 + SM.layer3;
% STOR.mapavg = nanmean(STOR.map, 3);
% STOR.ts = nanmean_space(STOR.map);
% 
% figure
% plot(timevector_monthly, STOR.ts)
% xlabel('Time')
% ylabel('Storage')
% title('Basin average storage (mm)')
% set(gca, 'fontsize', 20)
% 
% figure
% plotraster(masklon, masklat, STOR.mapavg, 'Time average storage (mm)', 'Lon', 'Lat')
% set(gca, 'fontsize', 20)

%% Make water balance object

WB = struct();
WB.time = timevector_monthly;

WB.SWE = nanmean_space(SWE);
WB.WDEW = nanmean_space(WDEW);
WB.MOIST1 = nanmean_space(SM.layer1);
WB.MOIST2 = nanmean_space(SM.layer2);
WB.MOIST3 = nanmean_space(SM.layer3);

WB.ET = nanmean_space(ET);
WB.BASEFLOW = nanmean_space(BASEFLOW);
WB.RUNOFF = nanmean_space(RUNOFF);
WB.P = nanmean_space(PREC);

% calculate total storage
% WB.S = WB.SWE + WB.MOIST1 + WB.MOIST2 + WB.MOIST3;
WB.S = WB.SWE + WB.WDEW + WB.MOIST1 + WB.MOIST2 + WB.MOIST3;

% calculate precipitation as residual
WB.P_resid = zeros(size(WB.S));
for t=2:nmonth
    WB.P_resid(t) = WB.ET(t) + WB.BASEFLOW(t) + WB.RUNOFF(t) + WB.S(t) - WB.S(t-1);
end

figure, plot(WB.time, WB.P), grid on
hold on, plot(WB.time, WB.P_resid)
legend('Precip', 'Precip (residual)', 'Location', 'Best');

%% Plot water balance time series (subplots)

figure
lw = 1.5;
fs = 18;

subplot(4,1,1)
plot(WB.time, [0; diff(WB.S)], 'linewidth', lw)
grid on
title('\Delta S (mm)')
set(gca, 'fontsize', fs)
xlabel('Time (months)')
ylabel('\Delta S (mm)')

subplot(4,1,2)
plot(WB.time, WB.ET, 'linewidth', lw)
grid on
title('ET (mm)')
set(gca, 'fontsize', fs)
xlabel('Time (months)')
ylabel('ET (mm)')

subplot(4,1,3)
plot(WB.time, WB.BASEFLOW + WB.RUNOFF, 'linewidth', lw)
grid on
title('Q (mm)')
set(gca, 'fontsize', fs)
xlabel('Time (months)')
ylabel('Q (mm)')

subplot(4,1,4)
plot(WB.time, WB.P, 'linewidth', lw)
grid on
title('P (mm)')
set(gca, 'fontsize', fs)
xlabel('Time (months)')
ylabel('P (mm)')

%% Plot water balance time series (one figure)

figure
plot(WB.time, WB.P, 'linewidth', lw)
hold on
plot(WB.time, WB.BASEFLOW + WB.RUNOFF, 'linewidth', lw)
plot(WB.time, WB.ET, 'linewidth', lw)
plot(WB.time, [0; diff(WB.S)], 'linewidth', lw)
legend('P', 'Q','ET','\Delta S','Location','Best')
xlabel('Time (months)')
ylabel('Value (mm)')
title('Upper Colorado Water Balance')
grid on
set(gca, 'fontsize', fs)

%% Calculate water-year average change in storage

for i=1:11
   storage_change(i) = WB.S(12*i) - WB.S(12*(i-1)+1);
end

%% Calculate water-year averages (nah, don't bother with it)

WB.WY = struct();
[WB.WY.time, WB.WY.P] = monthly_to_annual(WB.time, WB.P, 'sum');
[~, WB.WY.S] = monthly_to_annual(WB.time, WB.S, 'sum');
[~, WB.WY.Q] = monthly_to_annual(WB.time, WB.RUNOFF + WB.BASEFLOW, 'sum');
[~, WB.WY.ET] = monthly_to_annual(WB.time, WB.ET, 'sum');

figure
plot(WB.WY.time, WB.WY.P, 'linewidth', lw)
hold on
plot(WB.WY.time, WB.WY.Q, 'linewidth', lw)
plot(WB.WY.time, WB.WY.ET, 'linewidth', lw)
plot(WB.WY.time, [0; diff(WB.WY.S)], 'linewidth', lw)
legend('P', 'Q','ET','\Delta S','Location','Best')
xlabel('Time (months)')
ylabel('Value (mm)')
title('Upper Colorado Water Balance')
grid on
set(gca, 'fontsize', fs)

%% Runoff (I don't think this is working...)

livneh = load('/Volumes/HD4/SWOTDA/Data/Colorado/L15/livneh_q.mat');
i1 = find(livneh.timevector == start_date);
i2 = find(livneh.timevector == end_date);
Qd = livneh.runoff_sub(i1:i2,:);

npix = size(Qd,1); % this code block is inefficient
Qd_monthly = zeros(npix, nmonth);
for k=1:npix
    [~, Qd_monthly(k,:)] = daily_to_monthly(livneh.timevector(i1:i2)', Qd(k,:), 'sum');
    if mod(k,10)==0
        disp(k)
    end
end

RUNOFF.map_daily = array2map_vic(flipud(Qd), nanmask);


RUNOFF.map_monthly = array2map_vic(Qd_monthly, nanmask);

figure, plotraster(masklon, masklat, RUNOFF.map_monthly(:,:,103), 'April 2009 average runoff (mm)','','')

figure, plotraster(masklon, masklat, RUNOFF.map_daily(:,:,103), 'Average runoff (mm)','','')
caxis([0, 0.1])

%% Baseflow (meh, probably need to redo it)
tic
Qb = livneh.baseflow_sub(i1:i2,:);

Qb_monthly = zeros(npix, nmonth);
for k=1:npix
    [~, Qb_monthly(k,:)] = daily_to_monthly(livneh.timevector(i1:i2)', Qb(k,:), 'sum');
    if mod(k,10)==0
        disp(k)
    end
end

BASEFLOW.map_daily = array2map_vic(Qb, nanmask);
BASEFLOW.map_monthly = array2map_vic(Qb_monthly, nanmask);

figure, plotraster(masklon, masklat, BASEFLOW.map_monthly(:,:,103), 'April 2009 average baseflow (mm)','','')


toc
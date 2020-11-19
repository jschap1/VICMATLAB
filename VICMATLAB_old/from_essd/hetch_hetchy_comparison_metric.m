% Routing comparison for Tuolumne River Basin
%
% Comparing simulated flows to naturalized flow data at La Grange Dam, near
% Hetch Hetchy.
%
% 11/13/2019 JRS

%% Load gauge data

data_dir = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/gauge_data';
gaugeID = 'TLG';
ngages = length(gaugeID);
gage = struct();

fname = 'tuolumne_la_grange_discharge_TLG_naturalized.txt';
gage.(['gauge_' gaugeID]) = dlmread(fullfile(data_dir, fname), '\t', 1, 0); 

% remove bad (negative) values (why are there negative values?)
gage.q = gage.(['gauge_' gaugeID])(:,4);
gage.(['gauge_' gaugeID])(gage.q<0,4) = NaN;
gage.q = gage.(['gauge_' gaugeID])(:,4);
gage.q = gage.q*0.028316847; % convert from cfs to m^3/s

gage.time = (datetime(1982,10,1):days(1):datetime(2011,12,31))';

%% Load routing model outputs (VICGlobal)

indir = '/Volumes/HD3/SWOTDA/tlg_rout/rout_out_vg';
fname = 'tlg_station_cfs.txt';
rout = struct();
aa = load(fullfile(indir, fname));
rout.vg.time = datetime(aa(:,1), aa(:,2), aa(:,3));
rout.vg.q = aa(:,4); 
rout.vg.q = rout.vg.q*0.028316847; % convert from cfs to m^3/s

rout.vg.q = rout.vg.q(rout.vg.time>=gage.time(1));
rout.vg.time = rout.vg.time(rout.vg.time>=gage.time(1));


%% Load routing model outputs (Livneh)

indir = '/Volumes/HD3/SWOTDA/tlg_rout/rout_out_livneh';
fname = 'tlg_station_cfs.txt';
aa = load(fullfile(indir, fname));
rout.livneh.time = datetime(aa(:,1), aa(:,2), aa(:,3));
rout.livneh.q = aa(:,4);
rout.livneh.q = rout.livneh.q*0.028316847; % convert from cfs to m^3/s

rout.livneh.q = rout.livneh.q(rout.livneh.time>=gage.time(1));
rout.livneh.time = rout.livneh.time(rout.livneh.time>=gage.time(1));

%% Plot all three on one axis

linewidth = 1;

figure
title('La Grange Dam')
xlabel('Time')
ylabel('Discharge (m^3/s)')
plot(rout.livneh.time, rout.livneh.q, 'LineWidth', linewidth, 'Color', 'Blue')
hold on
plot(rout.vg.time, rout.vg.q, 'LineWidth', linewidth, 'Color', 'Red')
plot(gage.time, gage.q, 'LineWidth', linewidth, 'Color', 'Black')
legend('L2013','VICGlobal','Gage','Location','NW')
set(gca, 'fontsize', 18)

%% Log plot

figure
title('La Grange Dam')
xlabel('Time')
ylabel('Discharge (m^3/s)')
plot(rout.livneh.time, log(rout.livneh.q), 'LineWidth', linewidth, 'Color', 'Blue')
hold on
plot(rout.vg.time, log(rout.vg.q), 'LineWidth', linewidth, 'Color', 'Red')
plot(gage.time, log(gage.q), 'LineWidth', linewidth, 'Color', 'Black')
legend('L2013','VICGlobal','Gage','Location','NW')
set(gca, 'fontsize', 18)

%% Plot a 30-day moving average of the discharge data

linewidth = 1;

figure
title('La Grange Dam')
xlabel('Time')
ylabel('Discharge (m^3/s)')
plot(rout.livneh.time, rout.livneh.q, 'LineWidth', linewidth, 'Color', 'Blue')
hold on
plot(rout.vg.time, rout.vg.q, 'LineWidth', linewidth, 'Color', 'Red')
plot(gage.time, gage.q, 'LineWidth', linewidth, 'Color', 'Black')
legend('L2013','VICGlobal','Gage','Location','NW')
set(gca, 'fontsize', 18)

save('/Users/jschap/Documents/Research/VICGlobal/Data/upper_tuolumne_discharge.mat', 'rout','gage')

%% Calculate goodness of fit statistics

addpath('/Volumes/HD3/SWOTDA/Calibration')

NSE = struct();
NSE.livneh = myNSE(gage.q, rout.livneh.q);
NSE.vg = myNSE(gage.q, rout.vg.q);

RMSE = struct();
RMSE.livneh = myRMSE(gage.q, rout.livneh.q);
RMSE.vg = myRMSE(gage.q, rout.vg.q);

%% Annual average flows

doy = struct();
doy.vic = day(rout.livneh.time,'dayofyear');
doy.gage = day(gage.time, 'dayofyear');

% Calculate annual average flow

rout.livneh.meanq = zeros(366,1);
rout.vg.meanq = zeros(366,1);
gage.meanq = zeros(366,1);

for i=1:366
    rout.livneh.meanq(i) = mean(rout.livneh.q(doy.vic==i));
    rout.vg.meanq(i) = mean(rout.vg.q(doy.vic==i));
    gage.meanq(i) = nanmean(gage.q(doy.gage==i));
end

%% Plot average flows three on one axis

% Arrange by DOWY, not DOY
aa = struct();
% day(datetime(1980,10,1), 'dayofyear')
aa.gage = cy2wy2(gage.meanq, 275);
aa.vg = cy2wy2(rout.vg.meanq, 275);
aa.livneh = cy2wy2(rout.livneh.meanq, 275);

linewidth = 2;
figure
plot(1:366, aa.livneh, 'LineWidth', linewidth, 'Color', 'Blue')
hold on
plot(1:366, aa.vg, 'LineWidth', linewidth, 'Color', 'Red')
plot(1:366, aa.gage, 'LineWidth', linewidth, 'Color', 'Black')
legend('L2013','VICGlobal','Gage','Location','NE')
title('Annual average flow')
xlabel('Day of water year')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', 18)
xlim([1,366])

%% Figure for Dongyue

linewidth = 2;
figure
plot(1:366, aa.vg, 'LineWidth', linewidth, 'Color', 'Red')
hold on
plot(1:366, aa.gage, 'LineWidth', linewidth, 'Color', 'Black')
legend('VICGlobal','Gage','Location','NE')
title('Annual average flow')
xlabel('Day of water year')
ylabel('Discharge (m^3/s)')
set(gca, 'fontsize', 18)
xlim([1,366])


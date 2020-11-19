% Plotting water balance for one grid cell
%
% Results are from the VICGlobal UCRB simulation
%
% 2/26/2020 JRS

wbfile = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_1980-2011_VG/wb/wb_35.59375_-105.84375.txt';
dat = dlmread(wbfile, '\t', 3, 0);

% dat = dat(1:15,:);

WB = struct();
WB.time = datetime(dat(:,1), dat(:,2), dat(:,3));

WB.runoff = dat(:,5);
WB.baseflow = dat(:,6);
WB.evap = dat(:,4);
WB.precip = dat(:,13);

WB.swe = dat(:,27);
WB.sm1 = dat(:,21);
WB.sm2 = dat(:,22);
WB.sm3 = dat(:,23);
WB.wdew = dat(:,17);

% Calculate storage

WB.storage = WB.swe + WB.sm1 + WB.sm2 + WB.sm3 + WB.wdew;

% Calculate precipitation backwards as a check

ndays = length(WB.time);
P = zeros(ndays, 1);
for i=2:(ndays)
    P(i) = WB.runoff(i) + WB.baseflow(i) + WB.evap(i) + WB.storage(i) - WB.storage(i-1);
end

% This checks out. This is a valid way to back-calculate precipitation, using
% the water balance equation.

%% Plot water balance

figure
lw = 1.5;
fs = 18;

subplot(4,1,1)
% plot(WB.time, WB.storage, 'linewidth', lw)
plot(WB.time, [0; diff(WB.storage)], 'linewidth', lw)
grid on
title('\Delta S (mm)')
set(gca, 'fontsize', fs)
xlabel('Time (months)')
ylabel('\Delta S (mm)')

subplot(4,1,2)
plot(WB.time, WB.evap, 'linewidth', lw)
grid on
title('ET (mm)')
set(gca, 'fontsize', fs)
xlabel('Time (months)')
ylabel('ET (mm)')

subplot(4,1,3)
plot(WB.time, WB.baseflow + WB.runoff, 'linewidth', lw)
grid on
title('Q (mm)')
set(gca, 'fontsize', fs)
xlabel('Time (months)')
ylabel('Q (mm)')

subplot(4,1,4)
plot(WB.time, WB.precip, 'linewidth', lw)
grid on
title('P (mm)')
set(gca, 'fontsize', fs)
xlabel('Time (months)')
ylabel('P (mm)')
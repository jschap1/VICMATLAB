% Compare summed runoff and baseflow

% indir1 = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_7/Classic_L15/out/wb';
% indir2 = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_7/Classic_VG/out/wb';

indir1 = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_2000-2011_VGMOD/wb';

basin_mask = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif';
% basin_mask = '/Users/jschap/Documents/Research/Glaciers/Skagit/skagit_mask.tif';
[mask1, R1, lon1, lat1] = geotiffread2(basin_mask);

% indir1 = '/Users/jschap/Documents/Research/Glaciers/Skagit/out/wb';

fluxnames1 = dir(fullfile(indir1, '*.txt'));
fluxnames2 = dir(fullfile(indir2, '*.txt'));

dat = dlmread(fullfile(indir1, fluxnames1(1).name), '\t', 3, 0);
timevector = datetime(dat(:,1), dat(:,2), dat(:,3));
nt = length(timevector);

fluxnames1(mask1(:)==0) = [];
ncells = length(fluxnames1);

runoff_col = 5;
baseflow_col = 6;

runoff1 = NaN(ncells, nt);
baseflow1 = NaN(ncells, nt);

for k=1:ncells
    dat = dlmread(fullfile(indir1, fluxnames1(k).name), '\t', 3, 0);
    runoff1(k,:) = dat(:,runoff_col);
    baseflow1(k,:) = dat(:,baseflow_col);
    if mod(k,1000) == 0
        disp(k)
    end
end
runoff_sum1 = nansum(runoff1, 1);
baseflow_sum1 = nansum(baseflow1, 1);

% outdir = '/Users/jschap/Documents/Research/Glaciers/Skagit/out/processed';
outdir = '/Volumes/HD4/SWOTDA/Data/Colorado/EB_2000-2011_VGMOD/processed';
% outdir = '/Volumes/HD4/SWOTDA/Data/UpperMiss/out_wy1993/processed';
save(fullfile(outdir, 'summed_q_standard.mat'), 'nt','ncells','runoff1','baseflow1','runoff_sum1','baseflow_sum1', 'timevector')

%%
runoff2 = NaN(ncells, nt);
baseflow2 = NaN(ncells, nt);

for k=1:ncells
    dat = dlmread(fullfile(indir2, fluxnames2(k).name), '\t', 3, 0);
    runoff2(k,:) = dat(:,runoff_col);
    baseflow2(k,:) = dat(:,baseflow_col);
end
runoff_sum2 = sum(runoff2, 1);
baseflow_sum2 = sum(baseflow2, 1);

%% Convert summed discharge units

mean_area = 37.5; % km^2
summed_discharge = mean_area*(runoff_sum1 + baseflow_sum1)*1000/(24*3600); % m^3/s

figure, 
jsplot(timevector, summed_discharge, 'Summed discharge (UCRB, VICGlobal)', 'Time', 'Q (m^3/s)', 18)
grid on
hold on

% livneh = load('/Volumes/HD4/SWOTDA/Data/Colorado/L15/livneh_q.mat');

plot(livneh.timevector, livneh.summed_discharge*37.5*1000/3600/24)

%%
fs = 18;
upp = 6e4;

figure

subplot(3,1,1)
plot(timevector, runoff_sum1 + baseflow_sum1)
hold on
% plot(timevector, runoff_sum2 + baseflow_sum2)
title('Runoff + baseflow')
xlabel('Time')
ylabel('Runoff + baseflow (mm)')
% ylim([0,2e5])
ylim([0,upp])
grid on
set(gca, 'fontsize', fs)

subplot(3,1,2)
plot(timevector, baseflow_sum1)
hold on
% plot(timevector, baseflow_sum2)
title('Baseflow')
xlabel('Time')
ylabel('Baseflow (mm)')
% ylim([0,2e5])
ylim([0,upp])
grid on
set(gca, 'fontsize', fs)

subplot(3,1,3)
plot(timevector, runoff_sum1)
hold on
% plot(timevector, runoff_sum2)
title('Runoff')
xlabel('Time')
ylabel('Runoff (mm)')
% ylim([0,2e5])
ylim([0,upp])
grid on
set(gca, 'fontsize', fs)

% Runoff is comparable between the VICGlobal and L15 simulations.
%
% Baseflow is much higher for the VICGlobal parameters than for the Livneh
% et al. (2013) parameters. 
%
% These differences are probably most closely related to soil parameters,
% but vegetation cover could play a role, as well.
%
% A longer simulation period and comparison with gauge data is needed to
% really evaluate how well each simulation performs.
%
% It would be cool to compare soil parameters among the two Tuolumne
% setups. I will do this in my research notebook (LaTeX file)

%%

figure
subplot(2,1,1)
plot(timevector, baseflow_sum1, 'linewidth', 2)
hold on
plot(timevector, baseflow_sum2, 'linewidth', 2)
title('Baseflow')
xlabel('Time')
ylabel('Baseflow (mm)')
legend('L15','VG')
set(gca, 'fontsize', fs)

subplot(2,1,2)
plot(timevector, runoff_sum1, 'linewidth', 2)
hold on
plot(timevector, runoff_sum2, 'linewidth', 2)
title('Runoff')
xlabel('Time')
ylabel('Runoff (mm)')
legend('L15','VG')
set(gca, 'fontsize', fs)

%% Comparing VICGlobal and L2015 summed runoff and baseflow

% Water year 1993

infile1 = '/Volumes/HD4/SWOTDA/Data/UpperMiss/out_wy1993/processed/livneh_q.mat'; % L2015
infile2 = '/Volumes/HD4/SWOTDA/Data/UpperMiss/out_wy1993/processed/summed_q.mat'; % VICGlobal

Q1 = load(infile1);
Q2 = load(infile2);

fs = 18;

% Summed runoff plus baseflow
figure
plot(Q1.timevector, Q1.summed_discharge, 'linewidth', 2)
hold on
plot(Q1.timevector, Q2.runoff_sum1 + Q2.baseflow_sum1, 'linewidth', 2)
title('Baseflow + Runoff')
xlabel('Time')
ylabel('Baseflow + Runoff (mm)')
legend('L15','VG')
set(gca, 'fontsize', fs)

% Runoff
figure
plot(Q1.timevector, nansum(Q1.runoff_sub,2), 'linewidth', 2)
hold on
plot(Q1.timevector, Q2.runoff_sum1, 'linewidth', 2)
title('Runoff')
xlabel('Time')
ylabel('Runoff (mm)')
legend('L15','VG')
set(gca, 'fontsize', fs)

% Baseflow
figure
plot(Q1.timevector, nansum(Q1.baseflow_sub,2), 'linewidth', 2)
hold on
plot(Q1.timevector, Q2.baseflow_sum1, 'linewidth', 2)
title('Baseflow')
xlabel('Time')
ylabel('Baseflow (mm)')
legend('L15','VG')
set(gca, 'fontsize', fs)


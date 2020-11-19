% Loads, processes, and makes plots for VIC routing model outputs (Lohmann
% routing model)
%
% Updated 2/12/2019 JRS

%% Inputs

% Location of routing output files
rout_out_dir = './FDT/Rout/smaller/Output_mansnapped/CLSM_1959-1982';

% Station locations (latlon)
% stnlocs = './FDT/Rout/smaller/gages_snapped.txt';
stnlocs = './Data/Gauges/irb_stations_rgm_snap2.txt';

% Station location input file from VIC
stnvic = './FDT/Rout/smaller/stnloc_mansnapped.txt';

% prefix = 'IN001'; % name of gauge station/routing model output file prefix

% Path to VICMATLAB codes
addpath(genpath('/Users/jschap/Documents/Codes/VICMATLAB/'))

units = 'cfs'; % mm or cfs
timestep = 'day'; % day, month, or year
   
invisible = 0; % flag to turn on/off plotting
saveflag = 1;
saveloc = './FDT/Rout/smaller/Output_mansnapped/Figures';

%% Load routing results

T = readtable(stnlocs);
[nlocs,~] = size(T);
prefix = cell(nlocs, 1);

fID = fopen(stnvic, 'r');
for k=1:2*nlocs
    tmp = fgetl(fID);
    if mod(k, 2)~=0
        ind = strfind(tmp, 'IRB');
        prefix{(k+1)/2} = tmp(ind:ind+4);
%         prefix{(k+1)/2} = tmp(3:7);
    end
end
fclose(fID);

for k=1:nlocs
    rr = LoadVIC4Results(rout_out_dir, prefix{k}, units, timestep);
    
    switch timestep
        case 'day'
            ROUT.(strtrim(prefix{k})).time = datetime([rr(:,1), rr(:,2), rr(:,3)]);
        case 'month'
            ROUT.(strtrim(prefix{k})).time = datetime([rr(:,1), rr(:,2), 0]);
        case 'year'
            ROUT.(strtrim(prefix{k})).time = datetime([rr(:,1), 0, 0]);
    end
    ROUT.(strtrim(prefix{k})).discharge = rr(:,end);    
    ROUT.(strtrim(prefix{k})).location = T{k,1:2};
    
end

save('./FDT/Rout/smaller/Output_mansnapped/CLSM_1959-1982/ROUT.mat', 'ROUT');

%% Plot each routed streamflow time series

gaugenames = fieldnames(ROUT);

for k=1:nlocs
    figure
    plot(ROUT.IRB1.time, ROUT.(gaugenames{k}).discharge)
    title(['Gauge ' gaugenames{k}])
    xlabel('Time')
    ylabel('Discharge (cfs)')
end

% Plot inflow at Jinnah Reservoir

% Construct a monthly discharge time series from Tarbela Dam for comparison
t_in = [14874	16953	20106	37663	79487	136076	295251	268003	80936	41419	33966	23209]; % cfs
t1 = repmat(t_in(1), 31, 1);
t2 = repmat(t_in(2), 28, 1);
t3 = repmat(t_in(3), 31, 1);
t4 = repmat(t_in(4), 30, 1);
t5 = repmat(t_in(5), 31, 1);
t6 = repmat(t_in(6), 30, 1);
t7 = repmat(t_in(7), 31, 1);
t8 = repmat(t_in(8), 31, 1);
t9 = repmat(t_in(9), 30, 1);
t10 = repmat(t_in(10), 31, 1);
t11 = repmat(t_in(11), 30, 1);
t12 = repmat(t_in(12), 31, 1);
tarbela_inflow = vertcat(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12);
tarbela_inflow = repmat(tarbela_inflow, 3, 1);
tarbela_inflow(end+1) = tarbela_inflow(end);

% Convert to m^3/s
jinnah_inflow = ROUT.(gaugenames{k}).discharge*(12/39.37)^3;
tarbela_inflow = tarbela_inflow*(12/39.37)^3;

figure, k = 12; hold on;
plot(ROUT.IRB1.time, jinnah_inflow, 'b', 'linewidth', 2)
plot(ROUT.IRB1.time, tarbela_inflow, 'r', 'linewidth', 2)
% plot(ROUT.IRB1.time, repmat(mean(tarbela_inflow), 1096, 1), 'r', 'linewidth', 2)
title('Streamflow Routing')
xlabel('Time')
ylabel('Discharge (m^3/s)')
legend('Predicted inflow at Jinnah Barrage','Measured Inflow at Tarbela Dam (2015)', 'location', 'southoutside')
set(gca, 'fontsize', 18)


%% Scrap


%% Plots

figure
plot(ROUT.time, ROUT.ts)
titletext = ['Discharge time series (' timestep ') at ' prefix];
title(titletext); 
xlabel('time'); ylabel(['Q (' units ')'])
set(gca, 'FontSize', 14)

if saveflag
    saveas(gcf, fullfile(saveloc, [prefix '_' timestep '_ts.png']));
    % savefig(gcf, fullfile(saveloc, [prefix '_' timestep '_ts.fig']));
end

%% Load gage data

gage_data = dlmread('./Data/Gauges/irb_gage_data.txt', '\t', 2, 0);
gtime = datetime(gage_data(:,1), gage_data(:,2), gage_data(:,3));
pandoh_inflow = gage_data(:,8); % cfs
bhakra_inflow = gage_data(:,5); % cfs

%% Pandoh Dam

start_time = datetime(2013, 6, 15);
end_time = datetime(2014, 11, 30);

figure
subplot(2,1,1); hold on
[estim, ~] = plot_ts(start_time, end_time, 0.0283168*ROUT.IRB11.discharge, ROUT.IRB11.time);
[truth, ~] = plot_ts(start_time, end_time, 0.0283168*pandoh_inflow, gtime);
set(gca, 'fontsize', 18)
xlabel('Time')
ylabel('Discharge (m^3/s)')
legend('Estimated (GLDAS)','Measured (Gauge)')
title('Pandoh Dam Inflow')

%% Bhakra Dam

subplot(2,1,2); hold on
[estim, ~] = plot_ts(start_time, end_time, 0.0283168*ROUT.IRB5.discharge, ROUT.IRB5.time);
[truth, ~] = plot_ts(start_time, end_time, 0.0283168*bhakra_inflow, gtime);
set(gca, 'fontsize', 18)
xlabel('Time')
ylabel('Discharge (m^3/s)')
legend('Estimated (GLDAS)','Measured (Gauge)')
title('Bhakra Dam Inflow')

%% Validate against a USGS gauge

% Load USGS gauge data
gaugeID = 11290000; % 11290000 (non-ref), 11274630 (ref), 11303000 (non-ref)
fname = ['usgs' num2str(gaugeID) '.txt'];
fID = fopen(fullfile('/Users/jschapMac/Documents/Data/StreamGauges',fname));
fstring = '%s %d %s %f %s';
usgs_data = textscan(fID, fstring);
fclose(fID);

obs = usgs_data{4};

% Assuming the USGS data and the routing model output use the same start
% time and timestep, the ROUT.time vector can be reused. Otherwise, we
% would have to convert the usgs_data{3} entry to a datetime array. (That
% entry contains the timestamps for the USGS streamflow measurements.)

figure, hold on
plot(ROUT.time, ROUT.ts)
plot(ROUT.time, obs)
titletext = ['USGS gauge ' num2str(gaugeID) ' vs. routing model output'];
title(titletext); 
xlabel('time'); ylabel(['Q (' units ')'])
legend('ROUT','USGS');
set(gca, 'FontSize', 14)
hold off

if saveflag
    saveas(gcf, fullfile(saveloc, 'compare_RO.png'));
    savefig(gcf, fullfile(saveloc, 'compare_RO.fig'));
end

%%
% Compare results from different routing model runs

fd = load('/Users/jschapMac/Desktop/Tuolumne/Tuolumne8/Rout_Results/fd/ROUT_OUT.mat');
fd_Wu = load('/Users/jschapMac/Desktop/Tuolumne/Tuolumne8/Rout_Results/fd_Wu/ROUT_OUT.mat');
fd_Wu_Mu = load('/Users/jschapMac/Desktop/Tuolumne/Tuolumne8/Rout_Results/fd_Wu_Mu/ROUT_OUT.mat');

% Plot results on one figure
figure, hold on
plot(fd.ROUT.time, fd.ROUT.ts)
plot(fd.ROUT.time, fd_Wu.ROUT.ts)
plot(fd.ROUT.time, fd_Wu_Mu.ROUT.ts)
legend('1','2','3')

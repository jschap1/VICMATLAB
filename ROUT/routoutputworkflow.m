% Loads, processes, and makes plots for VIC routing model outputs (Lohmann
% routing model)
%
% Updated 2/12/2019 JRS

%% Inputs

% Location of routing output files
cd /Volumes/HD3/SWOTDA

% Station locations (latlon)
stnlocs = './Data/UMRB/stations.txt';

% Station location input file from VIC
stnvic = './Data/UMRB/rout_dl/umrb.stations';

% prefix = 'IN001'; % name of gauge station/routing model output file prefix

% Path to VICMATLAB codes
addpath('/Users/jschap/Documents/Codes/VICMATLAB')

units = 'cfs'; % mm or cfs
timestep = 'day'; % day, month, or year
   
invisible = 1; % flag to turn on/off plotting
saveflag = 1;
saveloc = './Outputs/ROUT_UMRB/plots';

%% Load routing results

T = readtable(stnlocs);
[nlocs,~] = size(T);
prefix = cell(nlocs, 1);

fID = fopen(stnvic, 'r');
for k=1:2*nlocs
    tmp = fgetl(fID);
    if mod(k, 2)~=0
        prefix{(k+1)/2} = tmp(3:7);
    end
end
fclose(fID);

cd /Volumes/HD3/SWOTDA/Outputs/ROUT_UMRB/Raw
for k=1:nlocs
    rr = LoadVIC4Results(prefix{k}, units, timestep);
    
    switch timestep
        case 'day'
            ROUT.(prefix{k}).time = datetime([rr(:,1), rr(:,2), rr(:,3)]);        
        case 'month'
            ROUT.(prefix{k}).time = datetime([rr(:,1), rr(:,2), 0]);
        case 'year'
            ROUT.(prefix{k}).time = datetime([rr(:,1), 0, 0]);
    end
    ROUT.(prefix{k}).discharge = rr(:,end);    
    ROUT.(prefix{k}).location = T{k,1:2};
    
end
cd /Volumes/HD3/SWOTDA/

save('ROUT.mat', 'ROUT');




%% Scrap

% 
% %% Plots
% 
% figure
% plot(ROUT.time, ROUT.ts)
% titletext = ['Discharge time series (' timestep ') at ' prefix];
% title(titletext); 
% xlabel('time'); ylabel(['Q (' units ')'])
% set(gca, 'FontSize', 14)
% 
% if saveflag
%     saveas(gcf, fullfile(saveloc, [prefix '_' timestep '_ts.png']));
%     % savefig(gcf, fullfile(saveloc, [prefix '_' timestep '_ts.fig']));
% end
% 
% %% Validate against a USGS gauge
% 
% % Load USGS gauge data
% gaugeID = 11290000; % 11290000 (non-ref), 11274630 (ref), 11303000 (non-ref)
% fname = ['usgs' num2str(gaugeID) '.txt'];
% fID = fopen(fullfile('/Users/jschapMac/Documents/Data/StreamGauges',fname));
% fstring = '%s %d %s %f %s';
% usgs_data = textscan(fID, fstring);
% fclose(fID);
% 
% obs = usgs_data{4};
% 
% % Assuming the USGS data and the routing model output use the same start
% % time and timestep, the ROUT.time vector can be reused. Otherwise, we
% % would have to convert the usgs_data{3} entry to a datetime array. (That
% % entry contains the timestamps for the USGS streamflow measurements.)
% 
% figure, hold on
% plot(ROUT.time, ROUT.ts)
% plot(ROUT.time, obs)
% titletext = ['USGS gauge ' num2str(gaugeID) ' vs. routing model output'];
% title(titletext); 
% xlabel('time'); ylabel(['Q (' units ')'])
% legend('ROUT','USGS');
% set(gca, 'FontSize', 14)
% hold off
% 
% if saveflag
%     saveas(gcf, fullfile(saveloc, 'compare_RO.png'));
%     savefig(gcf, fullfile(saveloc, 'compare_RO.fig'));
% end
% 
% %%
% % Compare results from different routing model runs
% 
% fd = load('/Users/jschapMac/Desktop/Tuolumne/Tuolumne8/Rout_Results/fd/ROUT_OUT.mat');
% fd_Wu = load('/Users/jschapMac/Desktop/Tuolumne/Tuolumne8/Rout_Results/fd_Wu/ROUT_OUT.mat');
% fd_Wu_Mu = load('/Users/jschapMac/Desktop/Tuolumne/Tuolumne8/Rout_Results/fd_Wu_Mu/ROUT_OUT.mat');
% 
% % Plot results on one figure
% figure, hold on
% plot(fd.ROUT.time, fd.ROUT.ts)
% plot(fd.ROUT.time, fd_Wu.ROUT.ts)
% plot(fd.ROUT.time, fd_Wu_Mu.ROUT.ts)
% legend('1','2','3')

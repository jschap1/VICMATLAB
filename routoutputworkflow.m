% % Generic workflow for loading and processing VIC results, and generating
% figures.
%
% Loads VIC results (flux, snow) and optionally routing model results and 
% arranges them into a nicely formatted structure. 
%
% Allows easy creation of time series plots and maps for fluxes
%
% Should be run from the same directory where the VIC results are located.
%
% Dependencies:
% LoadVICResults
% ProcessVICFluxResults
% ProcessVICSnowResults
% GetCoords

%% Inputs

% Location of routing output files
cd /Users/jschapMac/Desktop/Tuolumne2/RoutOutputs

% Path to VICMATLAB codes
addpath('/Users/jschapMac/Desktop/VIC/VICMATLAB')

prefix = 'PRPT '; % name of gauge station/routing model output file prefix
units = 'cfs'; % mm or cfs
timestep = 'day'; % day, month, or year
   
invisible = 1; % flag to turn on/off plotting
saveflag = 1;
saveloc = '/Users/jschapMac/Desktop/Tuolumne2/VICOutputs/wb2006-2011/Plots';

%% Load routing results

% Pull results from routing model output files into Matlab
routresults = LoadVICResults(prefix, units, timestep);

% Parse the results
switch timestep
    case 'day'
        ROUT.time = datetime([routresults(:,1), routresults(:,2), routresults(:,3)]);        
    case 'month'
        ROUT.time = datetime([routresults(:,1), routresults(:,2), 0]);;
    case 'year'
        ROUT.time = datetime([routresults(:,1), 0, 0]);
end
ROUT.ts = routresults(:,end);      
save('ROUT.mat', 'ROUT');

%% Plots

figure
plot(ROUT.time, ROUT.ts)
titletext = ['Discharge time series (' timestep ') at ' prefix];
title(titletext); 
xlabel('time'); ylabel(['Q (' units ')'])
set(gca, 'FontSize', 14)

if saveflag
    saveas(gcf, fullfile(saveloc, [prefix '_' timestep '_ts.png']));
    savefig(gcf, fullfile(saveloc, [prefix '_' timestep '_ts.fig']));
end

%% Validate against a USGS gauge

% Load USGS gauge data
gaugeID = 11290000; % 11290000 (non-ref), 11274630 (ref), 11303000 (non-ref)
fname = ['usgs' num2str(gaugeID) '.txt'];
fID = fopen(fullfile('/Users/jschapMac/Documents/HydrologyData/StreamGauges',fname));
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


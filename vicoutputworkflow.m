% Generic workflow for processing VIC results. 
%
% Loads VIC results (flux, snow) and optionally routing model results and 
% arranges them into a nicely formatted structure. 
%
% Allows easy creation of time series plots and maps for fluxes
%
% Must be run from the same directory where the VIC results are located.
%
% Dependencies:
% LoadVICResults
% ProcessVICFluxResults, GetDateTime
% GetCoords
% ProcessVICFluxResultsMaps

%% Inputs

% Path to VICMATLAB codes
addpath('/Users/jschapMac/Desktop/VIC/VICMATLAB')

rout = 0; % Specify whether or not to process routing results
if rout == 1
    prefix = 'STEHE'; % Provide info about routing files
    units = 'mm';
    timestep = 'daily';
    [gridcells, fluxresults, snowresults, routresults] = LoadVICResults(rout, prefix, units, timestep);
else
    % Pull results from VIC output files into Matlab
    [gridcells, fluxresults, snowresults] = LoadVICResults();
end

% Provide info about the VIC model run
precision = 5;
nlayers = 3;
run_type = 'WATER_BALANCE';
rec_interval = 'daily';

%%
% Process VIC flux results, arranging them into a structure with entries 
% for each grid cell
FLUXES = ProcessVICFluxResults(gridcells, fluxresults, nlayers, run_type, rec_interval);

[lat, lon] = GetCoords(gridcells, precision);
FLUXES.lat = lat;
FLUXES.lon = lon;

%% Plot flux time series

saveflag = 0;
figure

subplot(2,1,1)
plot(FLUXES.time, FLUXES.ts.fluxes_48_1875_120_6875.prec);
subplot(2,1,2)
plot(FLUXES.time, FLUXES.ts.fluxes_48_1875_120_6875.runoff);

% Save the time series of FLUXES.time
timevectorpath = '/Users/jschapMac/Desktop/Stehekin_4_2/results/vic/default';
timevector = FLUXES.time;
save(fullfile(timevectorpath, 'timevector'), 'timevector')

titletext = 'Daily wdew for (48.1875, 120.6875)';
xlabel('Time')
ylabel('wdew (mm)')
title(titletext)

if saveflag
    saveas(gcf, fullfile(saveloc, 'WDEW.png'));
    savefig(gcf, fullfile(saveloc, 'WDEW.fig'));
end

%% Plot flux map

% Choose a date to make map:
mapdate = datetime([2000, 12, 31]); % year, month, day

t_ind = find(FLUXES.time == mapdate);
FLUXES = ProcessVICFluxResultsMaps(FLUXES, t_ind);
disp(['Map generated for ' datestr(FLUXES.time(t_ind))]);

% 3D scatter
figure
plot3(FLUXES.lat, FLUXES.lon, FLUXES.maps.prec, '*')
xlabel('lon'), ylabel('lat'), zlabel('value (units)')
title('title'); grid on

% Filled scatter
figure
scatter(FLUXES.lat,FLUXES.lon,50,FLUXES.maps.prec,'filled')
colorbar

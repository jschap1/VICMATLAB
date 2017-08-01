% Generic workflow for processing VIC results. Must be run from the same
% directory where the VIC results are located.

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
precision = 4;
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

plot(FLUXES.time, FLUXES.ts.fluxes_48_1875_120_6875.wdew);

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
disp(['Map generated for ' datestr(FLUXES.time(t))]);

% Scatterplot

figure
plot3(FLUXES.lat, FLUXES.lon, FLUXES.maps.prec, '*')
xlabel('lon'), ylabel('lat'), zlabel('value (units)')
title('title'); grid on



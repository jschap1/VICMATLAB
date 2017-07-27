% Generic workflow for processing VIC results

addpath('/Users/jschapMac/Desktop/VIC/VICMATLAB')

cd('/Users/jschapMac/Desktop/results_VIC/longrun')

rout = 0; % Specify whether or not to process routing results
if rout == 1
    prefix = 'STEHE';
    units = 'mm';
    timestep = 'daily';
    [gridcells, fluxresults, snowresults, routresults] = LoadVICResults(rout, prefix, units, timestep);
else
    [gridcells, fluxresults, snowresults] = LoadVICResults();
end

nlayers = 3;
run_type = 'WATER_BALANCE';
rec_interval = 'daily';

%%
FLUXES = ProcessVICFluxResults(gridcells, fluxresults, nlayers, run_type, rec_interval);

%%
% Plot fluxes for a particular grid cell
saveflag = 0;
figure

plot(FLUXES.time, FLUXES.gridcell_fluxes_48_1875_120_6875.wdew);

titletext = 'Daily wdew for (48.1875, 120.6875)';
xlabel('Time')
ylabel('wdew (mm)')
title(titletext)

if saveflag
    saveas(gcf, fullfile(saveloc, 'PREC.png'));
    savefig(gcf, fullfile(saveloc, 'PREC.fig'));
end

%%

% Make a map

% The strategy here is to avoid looping through the data ntimesteps times.
% Instead, the processing is only done for the time where you want to make
% a map.

precision = 4;
t = 1;
FLUXES = ProcessVICFluxResultsMaps(gridcells, FLUXES, precision, t);
disp(['Map generated for ' datestr(FLUXES.time(t))]);

xyz = [FLUXES.lat, FLUXES.lon, FLUXES.maps.prec];

% xyz should be sufficient to make a georeferenced precipitation map for
% for the whole simulation region at time t. Hopefully, there is some
% Matlab function in the mapping toolbox that will do this.

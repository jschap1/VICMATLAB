function h = PlotVICResults(t, cellnum, fluxresults)

% General function for plotting VIC results. Plots the specified fluxes
% for 
%
% INPUTS
% t = datetime vector for VIC flux results. See GetDateTime.
% cellnum = number of the grid cell whose flux time series is to be plotted
%
% fluxresults = array of all VIC flux results. See LoadVICResults.
% nlayers = number of soil layers, as specified in the VIC global parameter
% file.
% FROZEN_SOILS = 1 if FROZEN_SOILS, 0 otherwise.
% FULL_ENERGY = 1 if FULL_ENERGY, 0 otherwise.
%
% OUTPUTS
% h = figure handle for the VIC flux plot

% What kind of VIC results are contained in the fluxresults input?
% Need to put in some kind of look-up table.

figure, plot(t,fluxresults(:,4,cellnum))

titletext = 'Daily PREC for grid cell 1 (48.1875, 120.6875) (mm)';
xlabel('Time')
ylabel('PREC (mm)')
title(titletext)

if saveflag
    saveas(gcf, fullfile(saveloc, 'PREC.png'));
    savefig(gcf, fullfile(saveloc, 'PREC.fig'));
end

end
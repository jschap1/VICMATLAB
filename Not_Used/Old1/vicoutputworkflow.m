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
precision = 4;
nlayers = 3;
run_type = 'WATER_BALANCE';
rec_interval = 'daily';

saveflag = 1;
saveloc = '/Users/jschapMac/Desktop/Tuolumne/Plots';

%%

% Process VIC flux results, arranging them into a structure with entries 
% for each grid cell
FLUXES = ProcessVICFluxResults(gridcells, fluxresults, nlayers, run_type, rec_interval);

[lat, lon] = GetCoords(gridcells, precision);
FLUXES.lat = lat;
FLUXES.lon = lon;

% Try to consolidate this into one big post-processing file

% Compute time average flux maps
FLUXES = timeavgfluxmaps(FLUXES);

% Compute basin average time series
fluxmaps = NaN(length(FLUXES.time),ncells);
for t_ind =1:length(FLUXES.time)
    fluxmaps(t,:) = ProcessVICFluxResultsMaps(FLUXES, t_ind);
end
FLUXES.avgts.(varnames{p}) = mean(fluxmaps,1);

%% Plot flux time series

figure

% % Choose lat/lon for time series plot (DOES NOT WORK)
% tslat = 48.1875;
% tslon = 120.6875;
% 
% fieldname_flux = ['fluxes_' num2str(tslat) '_' num2str(tslon)];

subplot(2,1,1)
plot(FLUXES.time, FLUXES.ts.fluxes_37_5938_121_0938.prec);
titletext = 'Daily precip. for a grid cell';
xlabel('Time')
ylabel('precipitation (mm)')
title(titletext)

subplot(2,1,2)
plot(FLUXES.time, FLUXES.ts.fluxes_37_5938_121_0938.runoff);
titletext = 'Daily runoff. for a grid cell';
xlabel('Time')
ylabel('runoff (mm)')
title(titletext)

% Save the time series of FLUXES.time
timevector = FLUXES.time;
save(fullfile(saveloc, 'timevector'), 'timevector')

if saveflag
    saveas(gcf, fullfile(saveloc, 'PREC.png'));
    savefig(gcf, fullfile(saveloc, 'RUNOFF.fig'));
end

%% Plot flux map

% Choose a date to make map:
mapdate = datetime([2000, 12, 31]); % year, month, day

t_ind = find(FLUXES.time == mapdate);
FLUXES = ProcessVICFluxResultsMaps(FLUXES, t_ind);
disp(['Map generated for ' datestr(FLUXES.time(t_ind))]);

% 3D scatter
figure
plot3(FLUXES.lon, FLUXES.lat, FLUXES.maps.prec, '*')
xlabel('lon'), ylabel('lat'), zlabel('value (units)')
title('title'); grid on

% Filled scatter
figure
scatter(FLUXES.lon,FLUXES.lat,50,FLUXES.maps.prec,'filled')
colorbar

%% Plot time average flux maps
cellnames = fieldnames(FLUXES.ts);
varnames = FLUXES.ts.(cellnames{1}).Properties.VariableNames;
unitscell = struct2cell(FLUXES.units);

for p=1:length(varnames)
    
    figure
      
    if ~strcmp(varnames{p}, 'moist')   
        scatter(FLUXES.lon,FLUXES.lat,50,FLUXES.avgmaps.(varnames{p}),'filled')
        title([datestr(FLUXES.time(1)) ' to ' datestr(FLUXES.time(end)) ...
        ' average ' varnames{p} ' (' unitscell{p} ')']);
        xlabel('lon (degrees)'); ylabel('lat (degrees)')
        set(gca, 'FontSize', 14)
        colorbar        
    else
        nlayers = size(FLUXES.ts.(cellnames{1}).(varnames{p}),2);
        for n = 1:nlayers
            subplot(nlayers, 1, n)
            scatter(FLUXES.lon,FLUXES.lat,50,FLUXES.avgmaps.(varnames{p})(:,n),'filled')
            title(['layer ' num2str(n) ' soil moisture (' unitscell{p} ')']);
            xlabel('lon (degrees)'); ylabel('lat (degrees)')
            set(gca, 'FontSize', 14)
            colorbar            
        end
    end

   if saveflag
        saveas(gcf, fullfile(saveloc, ['avg_' varnames{p} '.png']));
        savefig(gcf, fullfile(saveloc, ['avg_' varnames{p} '.fig']));
    end

end


% Generic workflow for loading and processing VIC results, and generating
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
% GetCoords

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

%% Post-processing
FLUXES = ProcessVICFluxResults(gridcells, fluxresults, nlayers, run_type, rec_interval);

if saveflag
    timevector = FLUXES.time;
    save(fullfile(saveloc,'timevector.mat'),'timevector')
end

[lat, lon] = GetCoords(gridcells, precision);
FLUXES.lat = lat;
FLUXES.lon = lon;

ncells = length(fieldnames(FLUXES.ts));
cellnames = fieldnames(FLUXES.ts);
varnames = FLUXES.ts.(cellnames{1}).Properties.VariableNames;

% Add basin average time series to FLUXES

for p=1:length(varnames)
    
    if ~strcmp(varnames{p}, 'moist')
        fluxarray = NaN(length(FLUXES.time),ncells);
        for k=1:ncells
            fluxarray(:,k) = FLUXES.ts.(cellnames{k}).(varnames{p});
        end
        FLUXES.avgts.(varnames{p}) = mean(fluxarray,2);
    else
        fluxarray = NaN(length(FLUXES.time),ncells, nlayers);
        for k=1:ncells
            fluxarray(:,k,:) = FLUXES.ts.(cellnames{k}).(varnames{p});
        end
        FLUXES.avgts.(varnames{p}) = squeeze(mean(fluxarray,2));
    end
    
end

% Add time-average flux maps to FLUXES

for p = 1:length(varnames)
    
    FLUXES.avgmaps.(varnames{p}) = NaN(ncells,1);
    if strcmp(varnames{p}, 'moist')
        FLUXES.avgmaps.(varnames{p}) = NaN(ncells,nlayers);
    end
    for k=1:ncells
        FLUXES.avgmaps.(varnames{p})(k,:) = mean(FLUXES.ts.(cellnames{k}).(varnames{p}));
    end
    
end

%% Generate figures

%% Plot basin average time series
unitscell = struct2cell(FLUXES.units);

for p=1:length(varnames)
    
    figure
      
    if ~strcmp(varnames{p}, 'moist')
        plot(FLUXES.time, FLUXES.avgts.(varnames{p}))
        titletext = ['Basin average ' varnames{p} ' (' unitscell{p} ')'];
        title(sprintf('%s_%d',titletext), 'Interpreter', 'none'); 
        xlabel('time'); ylabel(varnames{p})
        set(gca, 'FontSize', 14)  
    else
        for n = 1:nlayers
            subplot(nlayers, 1, n)
            plot(FLUXES.time, FLUXES.avgts.(varnames{p})(:,n))
            titletext = ['layer ' num2str(n) ' soil moisture (' unitscell{p} ')'];
            title(sprintf('%s_%d',titletext), 'Interpreter', 'none'); 
            xlabel('time'); ylabel(varnames{p})
            set(gca, 'FontSize', 14)          
        end
    end

   if saveflag
        saveas(gcf, fullfile(saveloc, ['avg_' varnames{p} 'ts.png']));
        savefig(gcf, fullfile(saveloc, ['avg_' varnames{p} 'ts.fig']));
    end

end

% Note: 'Interpreter','none' disables the LaTeX interpreter which reads
% underbar as "make subscript".

%% Plot time average flux maps
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
        saveas(gcf, fullfile(saveloc, ['avg_' varnames{p} 'map.png']));
        savefig(gcf, fullfile(saveloc, ['avg_' varnames{p} 'map.fig']));
    end

end

% Generic workflow for loading and processing VIC results, and generating
% figures. Should be run from the same directory where the VIC results are
% located. 
%
% VIC4
%
% Loads VIC results (flux, snow) and optionally routing model results and 
% arranges them into a nicely formatted structure. 
%
% Allows easy creation of time series plots and maps for fluxes
%
% Dependencies:
% LoadVICResults
% ProcessVICFluxResults
% ProcessVICSnowResults
% GetCoords

%% Inputs

% Path to VICMATLAB codes
addpath(genpath('/Users/jschap/Documents/Codes/VICMATLAB'))

% Provide info about the VIC model run
precision = 5;
nlayers = 3;
run_type = 'WATER_BALANCE'; % FULL_ENERGY or WATER_BALANCE
rec_interval = 'daily';

invisible = 1; % flag to turn on/off plotting
saveflag = 1;
saveloc = '/Volumes/HD3/SWOTDA/Outputs/VIC_UMRB';

%%
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%                            POST-PROCESSING
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

%% Load results into Matlab

% For large VIC runs, this is very RAM-intensive. Should optimize.

[gridcells, fluxresults, snowresults] = LoadVICResults();

%% Fluxes

FLUXES = ProcessVICFluxResults(gridcells, fluxresults, nlayers, run_type, rec_interval);
timevector = FLUXES.time;

% load('FLUXES.mat')
% load('timevector.mat')
% load('../gridcells.mat')

[lat, lon] = GetCoords(gridcells, precision);
FLUXES.lat = lat;
FLUXES.lon = lon;
scatter(FLUXES.lon,FLUXES.lat,50,ones(14015,1),'filled')

ncells = length(fieldnames(FLUXES.ts));
cellnames = fieldnames(FLUXES.ts);
fluxvarnames = FLUXES.ts.(cellnames{1}).Properties.VariableNames;

% Add basin average time series to FLUXES

for p=1:length(fluxvarnames)
    
    if ~strcmp(fluxvarnames{p}, 'moist')
        fluxarray = NaN(length(FLUXES.time),ncells);
        for k=1:ncells
            fluxarray(:,k) = FLUXES.ts.(cellnames{k}).(fluxvarnames{p});
        end
        FLUXES.avgts.(fluxvarnames{p}) = mean(fluxarray,2);
    else
        fluxarray = NaN(length(FLUXES.time),ncells, nlayers);
        for k=1:ncells
            fluxarray(:,k,:) = FLUXES.ts.(cellnames{k}).(fluxvarnames{p});
        end
        FLUXES.avgts.(fluxvarnames{p}) = squeeze(mean(fluxarray,2));
    end
    
end

% Add time-average flux maps to FLUXES

for p = 1:length(fluxvarnames)
    
    FLUXES.avgmaps.(fluxvarnames{p}) = NaN(ncells,1);
    if strcmp(fluxvarnames{p}, 'moist')
        FLUXES.avgmaps.(fluxvarnames{p}) = NaN(ncells,nlayers);
    end
    for k=1:ncells
        FLUXES.avgmaps.(fluxvarnames{p})(k,:) = mean(FLUXES.ts.(cellnames{k}).(fluxvarnames{p}));
    end
    
end

%% Snow states

SNOW = ProcessVICSnowResults(gridcells, snowresults, run_type, rec_interval);

%%%%%%%%%%%%
% The lat/lon and timevector for the snow states are the same as for
% the flux states (unless the model setup is unusual),so there is no need
% re-run the following lines of code:

% if saveflag 
% recalculate them for snow states
%     timevector = FLUXES.time;
%     save(fullfile(saveloc,'timevector.mat'),'timevector')
% end
% 
% [lat, lon] = GetCoords(gridcells, precision);
% FLUXES.lat = lat;
% FLUXES.lon = lon;
%%%%%%%%%%%%

SNOW.lat = FLUXES.lat;
SNOW.lon = FLUXES.lon;
snowvarnames = SNOW.ts.(cellnames{1}).Properties.VariableNames;

% Add basin average time series to SNOW

for p=1:length(snowvarnames)
    
    snowarray = NaN(length(timevector),ncells);
    for k=1:ncells
        snowarray(:,k) = SNOW.ts.(cellnames{k}).(snowvarnames{p});
    end
    SNOW.avgts.(snowvarnames{p}) = mean(snowarray,2);
    
end

% Add time-average snow state maps to SNOW

for p = 1:length(snowvarnames)
    
    SNOW.avgmaps.(snowvarnames{p}) = NaN(ncells,1);
    
    for k=1:ncells
        SNOW.avgmaps.(snowvarnames{p})(k,:) = mean(SNOW.ts.(cellnames{k}).(snowvarnames{p}));
    end
    
end

if saveflag
    if ~exist(saveloc,'dir')
        mkdir(saveloc)
    end    
    save(fullfile(saveloc,'timevector.mat'),'timevector')
%     save(fullfile(saveloc,'FLUXES.mat'),'FLUXES')
    save(fullfile(saveloc,'FLUXES.mat'),'FLUXES', '-v7.3')
    save(fullfile(saveloc,'SNOW.mat'),'SNOW')
end

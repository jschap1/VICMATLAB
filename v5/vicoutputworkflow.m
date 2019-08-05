% Generic workflow for loading and processing VIC results, and generating
% figures. Should be run from the same directory where the VIC results are
% located. 
%
% Use with VIC5 Classic Driver. 
%
% Use vicoutputworkflow_image.m if results are from the Image Driver.
%
% Loads VIC results (flux, snow) and optionally routing model results and 
% arranges them into a nicely formatted structure. 
%
% Allows easy creation of time series plots and maps for fluxes
%
% Dependencies:
% LoadVIC5Results
% ProcessVICFluxResults
% ProcessVICSnowResults
% GetCoords

%% Inputs

% Path to VICMATLAB codes
addpath(genpath('/Users/jschap/Documents/Codes/VICMATLAB'))

% Provide info about the VIC model run
precision = 5;
nlayers = 3;
run_type = 'WATER_BALANCE'; % FULL_ENERGY or WATER_BALANCE or FROZEN_SOIL
rec_interval = 'daily';

saveflag = 1;
% saveloc = '/Volumes/HD3/SWOTDA/Outputs/VIC_IRB/WB_05012019/Raw';
saveloc = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/MiniDomain/Raw';

%%
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%                            POST-PROCESSING
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

%% Load results into Matlab

% For large VIC runs, this is very RAM-intensive. Should optimize.

% default for VIC5 classic outputs is to have three lines of headers
headers = 3; 

% get number of flux variables in the output files
fluxnames = dir(fullfile(saveloc, 'fluxes*'));
fluxfile = dlmread(fullfile(saveloc, fluxnames(1).name), '\t', headers, 0); % skip headers
nfluxvars = size(fluxfile, 2);

snownames = dir(fullfile(saveloc, 'snow_*'));
snowfile = dlmread(fullfile(saveloc, snownames(1).name), '\t', headers, 0); % skip headers
nsnowvars = size(snowfile, 2);

[gridcells, fluxresults, snowresults] = LoadVIC5Results(headers, ...
    nfluxvars, nsnowvars, saveloc);

%% Fluxes

FLUXES = ProcessVICFluxResults(gridcells, fluxresults, nlayers, run_type, rec_interval);
timevector = FLUXES.time;

% load('FLUXES.mat')
% load('timevector.mat')
% load('../gridcells.mat')

[lat, lon] = GetCoords(gridcells, precision);
FLUXES.lat = lat;
FLUXES.lon = lon;
% scatter(FLUXES.lon,FLUXES.lat,50,ones(14015,1),'filled')

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
        FLUXES.avgts.(fluxvarnames{p}) = nanmean(fluxarray,2);
%         FLUXES.avgts.(fluxvarnames{p}) = mean(fluxarray,2);
    else
        fluxarray = NaN(length(FLUXES.time),ncells, nlayers);
        for k=1:ncells
            fluxarray(:,k,:) = FLUXES.ts.(cellnames{k}).(fluxvarnames{p});
        end
        FLUXES.avgts.(fluxvarnames{p}) = squeeze(nanmean(fluxarray,2));
%         FLUXES.avgts.(fluxvarnames{p}) = squeeze(mean(fluxarray,2));
    end
    
end

% Add time-average flux maps to FLUXES

for p = 1:length(fluxvarnames)
    
    FLUXES.avgmaps.(fluxvarnames{p}) = NaN(ncells,1);
    if strcmp(fluxvarnames{p}, 'moist')
        FLUXES.avgmaps.(fluxvarnames{p}) = NaN(ncells,nlayers);
    end
    for k=1:ncells
        FLUXES.avgmaps.(fluxvarnames{p})(k,:) = nanmean(FLUXES.ts.(cellnames{k}).(fluxvarnames{p}));
    end
    
end

%% Snow states

SNOW = ProcessVICSnowResults(gridcells, snowresults, run_type, rec_interval);

SNOW.lat = FLUXES.lat;
SNOW.lon = FLUXES.lon;

snowvarnames = SNOW.ts.(cellnames{1}).Properties.VariableNames;

% Add basin average time series to SNOW

for p=1:length(snowvarnames)
    
    snowarray = NaN(length(timevector),ncells);
    for k=1:ncells
        snowarray(:,k) = SNOW.ts.(cellnames{k}).(snowvarnames{p});
    end
    SNOW.avgts.(snowvarnames{p}) = nanmean(snowarray,2);
    
end

% Add time-average snow state maps to SNOW

for p = 1:length(snowvarnames)
    
    SNOW.avgmaps.(snowvarnames{p}) = NaN(ncells,1);
    
    for k=1:ncells
        SNOW.avgmaps.(snowvarnames{p})(k,:) = nanmean(SNOW.ts.(cellnames{k}).(snowvarnames{p}));
    end
    
end

%% Save

if saveflag
    if ~exist(saveloc,'dir')
        mkdir(saveloc)
    end    
    save(fullfile(saveloc,'../timevector.mat'),'timevector')
    save(fullfile(saveloc,'../FLUXES.mat'),'FLUXES')
%     save(fullfile(saveloc,'FLUXES.mat'),'FLUXES', '-v7.3')
    save(fullfile(saveloc,'../SNOW.mat'),'SNOW')
end

%     save(fullfile(saveloc,'./timevector.mat'),'timevector')
%     save(fullfile(saveloc,'./FLUXES.mat'),'FLUXES')
%     save(fullfile(saveloc,'./SNOW.mat'),'SNOW')


%%
% %%
% % vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
% %                               PLOTTING
% % vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
% 
% saveloc = '/Volumes/HD3/SWOTDA/Figures/VIC_UMRB/';
% 
% %% Plot basin average flux time series
% fluxunitscell = struct2cell(FLUXES.units);
% 
% for p=1:length(fluxvarnames)
%     
%     h = figure; % means
%     
%     if invisible == 1
%         set(h, 'Visible', 'off');
%     end
%       
%     if ~strcmp(fluxvarnames{p}, 'moist')
%         plot(FLUXES.time, FLUXES.avgts.(fluxvarnames{p}))
%         titletext = ['Basin average ' fluxvarnames{p} ' (' fluxunitscell{p} ')'];
%         title(sprintf('%s_%d',titletext), 'Interpreter', 'none'); 
%         xlabel('time'); ylabel(fluxvarnames{p})
%         set(gca, 'FontSize', 14)  
%     else
%         for n = 1:nlayers
%             subplot(nlayers, 1, n)
%             plot(FLUXES.time, FLUXES.avgts.(fluxvarnames{p})(:,n))
%             titletext = ['layer ' num2str(n) ' soil moisture (' fluxunitscell{p} ')'];
%             title(sprintf('%s_%d',titletext), 'Interpreter', 'none'); 
%             xlabel('time'); ylabel(fluxvarnames{p})
%             set(gca, 'FontSize', 14)          
%         end
%     end
% 
%    if saveflag
%         saveas(gcf, fullfile(saveloc, ['avg_' fluxvarnames{p} '_ts.png']));
%         savefig(gcf, fullfile(saveloc, ['avg_' fluxvarnames{p} '_ts.fig']));
%     end
% 
% end
% 
% % Note: 'Interpreter','none' disables the LaTeX interpreter which reads
% % underbar as "make subscript".
% 
% %% Plot time average flux maps
% fluxunitscell = struct2cell(FLUXES.units);
% 
% for p=1:length(fluxvarnames)
%     
%     h = figure; % means
%     
%     if invisible == 1
%         set(h, 'Visible', 'off');
%     end
%       
%     if ~strcmp(fluxvarnames{p}, 'moist')   
%         scatter(FLUXES.lon,FLUXES.lat,50,FLUXES.avgmaps.(fluxvarnames{p}),'filled')
% %         scatter(FLUXES.lon,FLUXES.lat,50,ones(14015,1),'filled')
%         title([datestr(FLUXES.time(1)) ' to ' datestr(FLUXES.time(end)) ...
%         ' average ' fluxvarnames{p} ' (' fluxunitscell{p} ')']);
%         xlabel('lon (degrees)'); ylabel('lat (degrees)')
%         set(gca, 'FontSize', 14)
%         colorbar        
%     else
%         for n = 1:nlayers
%             subplot(nlayers, 1, n)
%             scatter(FLUXES.lon,FLUXES.lat,50,FLUXES.avgmaps.(fluxvarnames{p})(:,n),'filled')
%             title(['layer ' num2str(n) ' soil moisture (' fluxunitscell{p} ')']);
%             xlabel('lon (degrees)'); ylabel('lat (degrees)')
%             set(gca, 'FontSize', 14)
%             colorbar            
%         end
%     end
% 
%    if saveflag
%         saveas(gcf, fullfile(saveloc, ['avg_' fluxvarnames{p} '_map.png']));
%         savefig(gcf, fullfile(saveloc, ['avg_' fluxvarnames{p} '_map.fig']));
%     end
% 
% end
% 
% %% Plot basin average snow state time series
% 
% % The basin average time series for swe and snow depth look strange.
% % After cell 19, the swe and snow depth time series increase each year.
% % Check the VIC model inputs and parameters for what may be causing this.
% 
% snowunitscell = struct2cell(SNOW.units);
% 
% for p=1:length(snowvarnames)
%     
%     h = figure; % means
%     
%     if invisible == 1
%         set(h, 'Visible', 'off');
%     end
%     
%     plot(timevector, SNOW.avgts.(snowvarnames{p}))
%     titletext = ['Basin average ' snowvarnames{p} ' (' snowunitscell{p} ')'];
%     title(sprintf('%s_%d',titletext), 'Interpreter', 'none'); 
%     xlabel('time'); ylabel(snowvarnames{p})
%     set(gca, 'FontSize', 14)  
% 
%    if saveflag
%         saveas(gcf, fullfile(saveloc, ['avg_' snowvarnames{p} '_ts.png']));
%         savefig(gcf, fullfile(saveloc, ['avg_' snowvarnames{p} '_ts.fig']));
%     end
% 
% end
% 
% %% Plot time average snow state maps
% snowunitscell = struct2cell(SNOW.units);
% 
% for p=1:length(snowvarnames)
%     
%     h = figure; % means
%     
%     if invisible == 1
%         set(h, 'Visible', 'off');
%     end
%       
%     scatter(lon,lat,50,SNOW.avgmaps.(snowvarnames{p}),'filled')
%     title([datestr(timevector(1)) ' to ' datestr(timevector(end)) ...
%     ' average ' snowvarnames{p} ' (' snowunitscell{p} ')']);
%     xlabel('lon (degrees)'); ylabel('lat (degrees)')
%     set(gca, 'FontSize', 14)
%     colorbar        
% 
%    if saveflag
%         saveas(gcf, fullfile(saveloc, ['avg_' snowvarnames{p} '_map.png']));
%         savefig(gcf, fullfile(saveloc, ['avg_' snowvarnames{p} '_map.fig']));
%     end
% 
% end

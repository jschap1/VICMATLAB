% Plots VIC flux outputs

load('./Outputs/VIC_UMRB/FLUXES.mat') % load processed flux data
saveloc = './Figures/VIC_UMRB/Fluxes/';

ncells = length(fieldnames(FLUXES.ts));
cellnames = fieldnames(FLUXES.ts);
fluxvarnames = FLUXES.ts.(cellnames{1}).Properties.VariableNames;
invisible = 1; % flag to turn on/off plotting
saveflag = 1;
nlayers = 3;

%% Plot basin average flux time series
% fluxunitscell = struct2cell(FLUXES.units);

% fluxvarnames = {'prec','evap','runoff','baseflow','wdew','moist','rad_temp', ...
%     'net_short','r_net','latent','evap_canop','evap_veg','evap_bare','sub_canop', ...
%     'sub_snow','sensible','grnd_flux','deltah','fusion','aero_resist','surf_temp', ...
%     'albedo','rel_humid','in_long','air_temp','wind'};

%   1×20 cell array
%   Columns 1 through 13
%     'prec'    'evap'    'runoff'    'baseflow'    'wdew'    'moist'    'net_short'    'r_net'    'evap_canop'    'evap_veg'    'evap_bare'    'sub_canop'    'sub_snow'
%   Columns 14 through 20
%     'aero_resist'    'surf_temp'    'albedo'    'rel_humid'    'in_long'    'air_temp'    'wind'

fluxunitscell = {'mm','mm','mm','mm','mm','mm','W/m^2','W/m^2','mm','mm','mm', ...
    'mm','mm','s/m','deg. C','-','-','W/m^2','deg. C', 'm/s'};

for p=1:length(fluxvarnames)
    
    h = figure; % means
    
    if invisible == 1
        set(h, 'Visible', 'off');
    end
      
    if ~strcmp(fluxvarnames{p}, 'moist')
        plot(FLUXES.time, FLUXES.avgts.(fluxvarnames{p}))
        titletext = ['Basin average ' fluxvarnames{p} ' (' fluxunitscell{p} ')'];
        title(sprintf('%s_%d',titletext), 'Interpreter', 'none'); 
        xlabel('time'); ylabel(fluxvarnames{p})
        set(gca, 'FontSize', 14)  
    else
        for n = 1:nlayers
            subplot(nlayers, 1, n)
            plot(FLUXES.time, FLUXES.avgts.(fluxvarnames{p})(:,n))
            titletext = ['layer ' num2str(n) ' soil moisture (' fluxunitscell{p} ')'];
            title(sprintf('%s_%d',titletext), 'Interpreter', 'none'); 
            xlabel('time'); ylabel(fluxvarnames{p})
            set(gca, 'FontSize', 14)          
        end
    end

   if saveflag
        saveas(gcf, fullfile(saveloc, ['avg_' fluxvarnames{p} '_ts.png']));
%         savefig(gcf, fullfile(saveloc, ['avg_' fluxvarnames{p} '_ts.fig']));
    end

end

% Note: 'Interpreter','none' disables the LaTeX interpreter which reads
% underbar as "make subscript".

%% Plot time average flux maps
% fluxunitscell = struct2cell(FLUXES.units);

for p=1:length(fluxvarnames)
    
    h = figure; % means
    
    if invisible == 1
        set(h, 'Visible', 'off');
    end
      
    if ~strcmp(fluxvarnames{p}, 'moist')   
        scatter(FLUXES.lon,FLUXES.lat,50,FLUXES.avgmaps.(fluxvarnames{p}),'filled')
%         scatter(FLUXES.lon,FLUXES.lat,50,ones(14015,1),'filled')
        title([datestr(FLUXES.time(1)) ' to ' datestr(FLUXES.time(end)) ...
        ' average ' fluxvarnames{p} ' (' fluxunitscell{p} ')']);
        xlabel('lon (degrees)'); ylabel('lat (degrees)')
        set(gca, 'FontSize', 14)
        colorbar        
    else
        for n = 1:nlayers
            subplot(nlayers, 1, n)
            scatter(FLUXES.lon,FLUXES.lat,50,FLUXES.avgmaps.(fluxvarnames{p})(:,n),'filled')
            title(['layer ' num2str(n) ' soil moisture (' fluxunitscell{p} ')']);
            xlabel('lon (degrees)'); ylabel('lat (degrees)')
            set(gca, 'FontSize', 14)
            colorbar            
        end
    end

   if saveflag
        saveas(gcf, fullfile(saveloc, ['avg_' fluxvarnames{p} '_map.png']));
%         savefig(gcf, fullfile(saveloc, ['avg_' fluxvarnames{p} '_map.fig']));
    end

end

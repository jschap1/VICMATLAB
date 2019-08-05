% Plots VIC flux outputs

load('/Volumes/HD3/SWOTDA/Data/IRB/VIC/MiniDomain/FLUXES.mat') % load processed flux data
saveloc = '/Volumes/HD3/SWOTDA/Data/IRB/VIC/MiniDomain/Plots/';

ncells = length(fieldnames(FLUXES.ts));
nsteps = length(FLUXES.time);
cellnames = fieldnames(FLUXES.ts);
fluxvarnames = FLUXES.ts.(cellnames{1}).Properties.VariableNames;
invisible = 1; % flag to turn on/off plotting
saveflag = 1;
nlayers = 3;

VICmode = 'WB';

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

if strcmp(VICmode, 'WB') % Water balance, three soil layers
    fluxunitscell = {'mm','mm','mm','mm','mm','mm','W/m^2','W/m^2','mm','mm','mm', ...
        'mm','mm','s/m','deg. C','-','-','W/m^2','deg. C', 'm/s'};
elseif strcmp(VICmode, 'EB') % Energy balance, three soil layers
    fluxunitscell = {
        'mm', % prec
        'mm', % evap
        'mm', % runoff
        'mm', % baseflow
        'mm', % wdew
        'mm', % moist
        'K', % rad_temp
        'W/m^2', % net_short
        'W/m^2', % r_net
        'W/m^2', % latent
        'mm', % evap_canop
        'mm', % evap_veg
        'mm', % evap_bare
        'mm', % sub_canop
        'mm', % sub_snow
        'W/m^2', % sensible
        'W/m^2', % ground heat flux
        'W/m^2', % delta_h
        'W/m^2', % fusion
        's/m', % aero_resist
        'deg. C', % surf_temp
        '-', % albedo
        '-', % relative_humidity
        'W/m^2', % incoming_longwave
        'deg. C',  % air_temp
        'm/s'}; % wind speed
end

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
   else
       pause;
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

    % save a geotiff version, as well
    A = flipud(xyz2grid(FLUXES.lon, FLUXES.lat, FLUXES.avgmaps.(fluxvarnames{p})(:,1)));
    xres = 1/16;
    yres = 1/16;
    R = makerefmat(min(FLUXES.lon), min(FLUXES.lat), xres, yres);
    geotiffwrite(fullfile(saveloc, ['avg_' fluxvarnames{p} '_map.tif']), A, R)

   else
       pause;
   end

end


%% Make nice plots of water balance components

lwd = 2;
% load('./Outputs/VIC_IRB/WB_corrected_veg/SNOW.mat')

figure

subplot(4,2,1)
plot(FLUXES.time, FLUXES.avgts.prec, 'linewidth', lwd)
title('Precipitation'), xlabel('Time'), ylabel('P (mm)')
set(gca, 'fontsize', 18)

subplot(4,2,2)
plot(FLUXES.time, FLUXES.avgts.evap, 'linewidth', lwd)
title('Evaporation'), xlabel('Time'), ylabel('E (mm)')
set(gca, 'fontsize', 18)

subplot(4,2,3)
plot(FLUXES.time, FLUXES.avgts.runoff, 'linewidth', lwd)
title('Runoff'), xlabel('Time'), ylabel('Q_r (mm)')
set(gca, 'fontsize', 18)

subplot(4,2,4)
plot(FLUXES.time, FLUXES.avgts.baseflow, 'linewidth', lwd)
title('Baseflow'), xlabel('Time'), ylabel('Q_b (mm)')
set(gca, 'fontsize', 18)

subplot(4,2,5)
plot(FLUXES.time, FLUXES.avgts.moist(:,1), 'linewidth', lwd)
title('Soil moisture (layer 1)'), xlabel('Time'), ylabel('SM (mm)')
set(gca, 'fontsize', 18)

subplot(4,2,6)
plot(FLUXES.time, FLUXES.avgts.moist(:,2), 'linewidth', lwd)
title('Soil moisture (layer 2)'), xlabel('Time'), ylabel('SM (mm)')
set(gca, 'fontsize', 18)

subplot(4,2,7)
plot(FLUXES.time, FLUXES.avgts.moist(:,3), 'linewidth', lwd)
title('Soil moisture (layer 3)'), xlabel('Time'), ylabel('SM (mm)')
set(gca, 'fontsize', 18)

subplot(4,2,8)
plot(FLUXES.time, SNOW.avgts.swe, 'linewidth', lwd)
title('Snow-water equivalent'), xlabel('Time'), ylabel('SWE (mm)')
set(gca, 'fontsize', 18)


%% Calculate water balance components (in meaningful units)

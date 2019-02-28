% Plots VIC snow outputs

load('./Outputs/VIC_UMRB/SNOW.mat') % load processed snow state data
saveloc = './Figures/VIC_UMRB/Snow';
invisible = 1;
saveflag = 1;

%% Plot basin average snow state time series

% The basin average time series for swe and snow depth look strange.
% After cell 19, the swe and snow depth time series increase each year.
% Check the VIC model inputs and parameters for what may be causing this.


ncells = length(fieldnames(SNOW.ts));
cellnames = fieldnames(SNOW.ts);
snowvarnames = SNOW.ts.(cellnames{1}).Properties.VariableNames;
snowunitscell = struct2cell(SNOW.units);
timevector = SNOW.time;
lat = SNOW.lat;
lon = SNOW.lon;

for p=1:length(snowvarnames)
    
    h = figure; % means
    
    if invisible == 1
        set(h, 'Visible', 'off');
    end
    
    plot(timevector, SNOW.avgts.(snowvarnames{p}))
    titletext = ['Basin average ' snowvarnames{p} ' (' snowunitscell{p} ')'];
    title(sprintf('%s_%d',titletext), 'Interpreter', 'none'); 
    xlabel('time'); ylabel(snowvarnames{p})
    set(gca, 'FontSize', 14)  

   if saveflag
        saveas(gcf, fullfile(saveloc, ['avg_' snowvarnames{p} '_ts.png']));
%         savefig(gcf, fullfile(saveloc, ['avg_' snowvarnames{p} '_ts.fig']));
    end

end

%% Plot time average snow state maps
snowunitscell = struct2cell(SNOW.units);

for p=1:length(snowvarnames)
    
    h = figure; % means
    
    if invisible == 1
        set(h, 'Visible', 'off');
    end
      
    scatter(lon,lat,50,SNOW.avgmaps.(snowvarnames{p}),'filled')
    title([datestr(timevector(1)) ' to ' datestr(timevector(end)) ...
    ' average ' snowvarnames{p} ' (' snowunitscell{p} ')']);
    xlabel('lon (degrees)'); ylabel('lat (degrees)')
    set(gca, 'FontSize', 14)
    colorbar        

   if saveflag
        saveas(gcf, fullfile(saveloc, ['avg_' snowvarnames{p} '_map.png']));
%         savefig(gcf, fullfile(saveloc, ['avg_' snowvarnames{p} '_map.fig']));
    end

end

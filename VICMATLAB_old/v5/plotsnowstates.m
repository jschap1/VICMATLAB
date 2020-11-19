% Plots VIC snow outputs

load('./Outputs/VIC_IRB/WB_corrected_veg/SNOW.mat') % load processed snow state data
saveloc = './Figures/VIC_IRB/Water_Balance_MERIT/Snow';
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

%% Answering a question from Dennis about snow accumulation

% how many cells are accumulating snow?
% cells with at least seasonal snow (mean max exceeding 25 cm per year)

SNOW.ts.cell_37_03125_74_84375_txt

% Plot maps of mean snow depth in each cell

map1 = flipud(xyz2grid(SNOW.lon, SNOW.lat, SNOW.avgmaps.snow_depth));
figure
imagesc(SNOW.lon, SNOW.lat, map1./10)
set(gca, 'ydir', 'normal')
xlabel('Lon')
ylabel('Lat')
title('Average snow depth (cm)')

% cells with average snow depth higher than 25 cm

% Which cells reach at least 25 cm of snow depth during a typical year?
% Of these cells, how many experience runaway snow accumulation?

%% Plotting all snow variables

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

    % save a geotiff version, as well
    A = flipud(xyz2grid(SNOW.lon, SNOW.lat, SNOW.avgmaps.(snowvarnames{p})));
    xres = 1/16;
    yres = 1/16;
    R = makerefmat(min(SNOW.lon), min(SNOW.lat), xres, yres);
%     figure, imagesc(SNOW.lon, SNOW.lat, A)
%     set(gca, 'ydir', 'normal')
    geotiffwrite(fullfile(saveloc, ['avg_' snowvarnames{p} '_map.tif']), A, R)

    end

end



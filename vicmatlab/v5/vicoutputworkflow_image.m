% Output workflow for VIC Image Driver
% Loads results into Matlab, processes them, and makes plots

clear, clc

%% Inputs

% Path to VICMATLAB codes
addpath('/Users/jschap/Documents/Codes/VICMATLAB')

% Directory of VIC outputs
cd /Volumes/HD4/SWOTDA/Data/Tuolumne/Image_VICGlobal/Results_EB_SB

% Name of fluxes file
fluxname = 'fluxes.1999-01-01.nc';

invisible = 0;

nlayers = 3; % number of soil layers

saveflag = 1;

saveloc = '/Volumes/HD4/SWOTDA/Data/Tuolumne/Image_VICGlobal/Results_EB_SB/Figures';

%%

% Load NetCDF file contents into Matlab
info = ncinfo(fluxname);
numvars = length(info.Variables);
varnames = cell(numvars,1);

for p=1:numvars
    varnames{p} = info.Variables(p).Name;
    var = ncread(fluxname, varnames{p});
    if p==10
        continue;
    end
    FLUXES.(varnames{p}) = permute(var, [2,1,3]);
end

% Convert to datetime:
FLUXES.time = datestr(FLUXES.time+367);
FLUXES.time = datetime(FLUXES.time);
% Adding 367 seems to work for me... I think this is because Matlab's epoch
% date is Jan. 1 0000, and the VIC output is days since Jan. 1 0001. Year
% 0000 was a leap year. Should check this with someone more experienced
% with VIC.

% Make separate entries for each layer of soil moisture
for q = 1:nlayers
    name1 = ['OUT_SOIL_MOIST_' num2str(q)];
    name2 = ['OUT_SOIL_TEMP_' num2str(q)];
%     FLUXES.(name1) = squeeze(FLUXES.OUT_SOIL_MOIST(:,:,q,:));
%     FLUXES.(name2) = squeeze(FLUXES.OUT_SOIL_TEMP(:,:,q,:));
end
% FLUXES = rmfield(FLUXES, 'OUT_SOIL_MOIST');
% FLUXES = rmfield(FLUXES, 'OUT_SOIL_TEMP');

% Update varnames, remove time, lat, lon, etc.
allnames = fieldnames(FLUXES);
varnames = allnames(5:end);
numvars = length(varnames);

% Extract timevector for convenience
timevector = FLUXES.time;

% Compute averages:
for p=1:numvars
    
    % Time average maps
    FLUXES.avgmaps.(varnames{p}) = nanmean(FLUXES.(varnames{p}),3);
    
    % Basin average time series
    for t=1:length(timevector)
        fluxvar = FLUXES.(varnames{p})(:,:,t);
        FLUXES.avgts.(varnames{p})(t) = nanmean(fluxvar(:));
    end
    
end

%% Plot basin average time series

for p=1:numvars
    
    h = figure;
    
    if invisible == 1
        set(h, 'Visible', 'off');
    end
      
    plot(timevector, FLUXES.avgts.(varnames{p}));
    titletext = ['Basin average ' varnames{p}];
    title(sprintf('%s_%d',titletext), 'Interpreter', 'none'); 
    xlabel('time'); ylabel(varnames{p})
    set(gca, 'FontSize', 14)     

   if saveflag
        saveas(gcf, fullfile(saveloc, ['avg_' varnames{p} '_ts.png']));
        savefig(gcf, fullfile(saveloc, ['avg_' varnames{p} '_ts.fig']));
    end

end
        
%% Plot time average maps

for p=1:numvars
    
    h = figure;
    
    if invisible == 1
        set(h, 'Visible', 'off');
    end
      
    if p==6
       continue
    end
    plotraster(FLUXES.lon, FLUXES.lat, FLUXES.avgmaps.(varnames{p}), varnames{p}, 'Lon', 'Lat');
%     imagesc(FLUXES.lon, FLUXES.lat, FLUXES.avgmaps.(varnames{p}));
    set(gca,'YDir','normal')
    titletext = [datestr(timevector(1)) ' to ' datestr(timevector(end)) ...
    ' average ' varnames{p}];
    title(sprintf('%s_%d',titletext), 'Interpreter', 'none');
    xlabel('lon (degrees)'); ylabel('lat (degrees)')
    set(gca, 'FontSize', 14)
    colorbar        

   if saveflag
        saveas(gcf, fullfile(saveloc, ['avg_' varnames{p} '_map.png']));
        savefig(gcf, fullfile(saveloc, ['avg_' varnames{p} '_map.fig']));
    end

end

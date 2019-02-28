% Loads and plots ASCII forcing data
%
% VIC 5 classic
%
% Note: disaggregated forcings are 3-hourly for the nominal VIC5 run
%
% Dependencies:
% GetCoords.m
%
% addpath(genpath('/Users/jschap/Documents/Codes/VICMATLAB/'))

%% Specify inputs

forcingpath = './Data/Livneh_Forcings/Disagg_Forc'; 

precision = 5;

% order of the forcing variables must match the order in the forcing file
varnames = {'prec','air_temp','shortwave','longwave','density','pressure','vp','wind'};
varunits = {'mm','deg C','W/m^2','W/m^2','kg/m^3','kPa','kPa','m/s'};

invisible = 1;
saveflag = 1;
saveloc = './Figures/VIC_UMRB'; 

%% Load forcing data

forcenames = dir(fullfile(forcingpath, 'full_data*'));

ncells = length(forcenames);
addpath(forcingpath)

% dlmread requires ASCII forcing data.
tmpforc = dlmread(forcenames(1).name); 

nsteps = size(tmpforc,1);
nvars = size(tmpforc,2);

%%% Adapted from LoadVICResults
gridcells = cell(ncells, 1);
for k=1:ncells
    tmpstring = forcenames(k).name;
    tmpstring = strrep(tmpstring,'-',''); % remove some characters bc Matlab cannot handle them
    tmpstring = strrep(tmpstring,'.','_');
    gridcells{k} = tmpstring;
end
%%%

[lat, lon] = GetCoords(gridcells, precision);

FORC = NaN(nsteps, nvars, ncells); % this variable is LARGE.
for k=1:ncells
    FORC(:,:,k) = dlmread(forcenames(k).name); 
end

if saveflag
    % save(fullfile(saveloc,'FORC.mat'),'FORC')
    save(fullfile(saveloc,'FORC.mat'),'FORC', 'lat', 'lon', 'varnames', 'varunits', '-v7.3')
end

%% Process data
% The best way to do this might be something like ProcessVICFluxResults
% But for now, this is a quick and dirty way to do it

% Compute basin average time series

AVGFORC = NaN(nsteps,nvars);

for p=1:nvars
    
    forcarray = NaN(nsteps,ncells);
    for k=1:ncells
        forcarray(:,k) = FORC(:,p,k);
    end
    AVGFORC(:,p) = mean(forcarray,2);
    disp(p)
end

%% Plot basin average time series

% Optionally load time series vector from VIC output
timev = 1;
if timev
    timevector = load('./Outputs/VIC_UMRB/timevector.mat');
    timevector = timevector.timevector;
    timevector = linspace(timevector(1), timevector(end), nsteps);
else
    timevector = 1:nsteps;
end

for p = 1:nvars
    
    h = figure; 
    
    if invisible
        set(h, 'Visible', 'off');
    end
    
    plot(timevector, AVGFORC(:,p))
    titletext = ['Basin average ' varnames{p} ' (' varunits{p} ')'];
    title(titletext)
    xlabel('time'); ylabel(varnames{p})
    set(gca, 'FontSize', 14)
    
    if saveflag
        saveas(gcf, fullfile(saveloc, ['avg_' varnames{p} 'ts.png']));
%         savefig(gcf, fullfile(saveloc, ['avg_' varnames{p} 'ts.fig']));
    end    
    
end

%% Plot time average maps

FORC_maps = squeeze(mean(FORC, 1));

for p=1:length(varnames)
    
    h = figure; % means
    
    if invisible == 1
        set(h, 'Visible', 'off');
    end
      
    scatter(lon,lat,50,FORC_maps(p,:),'filled')
%   scatter(FLUXES.lon,FLUXES.lat,50,ones(14015,1),'filled')
    title([datestr(timevector(1)) ' to ' datestr(timevector(end)) ...
    ' average ' varnames{p} ' (' varunits{p} ')']);
    xlabel('lon (degrees)'); ylabel('lat (degrees)')
    set(gca, 'FontSize', 14)
    colorbar        

   if saveflag
        saveas(gcf, fullfile(saveloc, ['avg_' varnames{p} '_map.png']));
%         savefig(gcf, fullfile(saveloc, ['avg_' fluxvarnames{p} '_map.fig']));
    end

end

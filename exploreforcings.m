% Loads and plots ASCII forcing data
%
% Dependencies:
% GetCoords.m

%% Specify inputs

forcingpath = '/Users/jschapMac/Desktop/Tuolumne/ClippedForcings/'; 
% the final slash is needed

precision = 4;

varnames = {'prec','tmin','tmax','wind'};
varunits = {'mm','deg C','deg C','m/s'};
% order of the forcing variables must match the order in the forcing files

saveflag = 1;
saveloc = '/Users/jschapMac/Desktop/Tuolumne/Plots/Forcings'; 

%% Load forcing data

forcenames = dir([forcingpath 'data*']);
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

FORC = NaN(nsteps, nvars, ncells);
for k=1:ncells
    FORC(:,:,k) = dlmread(forcenames(k).name);
end

%% Process data

% Compute basin average time series

AVGFORC = NaN(nsteps,nvars);

for p=1:nvars
    
    forcarray = NaN(nsteps,ncells);
    for k=1:ncells
        forcarray(:,k) = FORC(:,p,k);
    end
    AVGFORC(:,p) = mean(forcarray,2);
    
end

%% Plot basin average time series

% Optionally load time series vector from VIC output
timev = 0;
if timev
    timevector = load('timevector.mat');
    timevector = timevector.timevector;
else
    timevector = 1:nsteps;
end

for p = 1:nvars
    
    figure;
    plot(timevector, AVGFORC(:,p))
    titletext = ['Basin average ' varnames{p} ' (' varunits{p} ')'];
    title(titletext)
    xlabel('time'); ylabel(varnames{p})
    set(gca, 'FontSize', 14)
    
    if saveflag
        saveas(gcf, fullfile(saveloc, ['avg_' varnames{p} 'ts.png']));
        savefig(gcf, fullfile(saveloc, ['avg_' varnames{p} 'ts.fig']));
    end    
    
end

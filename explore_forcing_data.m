% Loads and facilitates viewing of ASCII forcing data

% This code does not work yet; it is a work in progress

%% Load forcing data
forcingpath = '/Users/jschapMac/Desktop/Stehekin_4_2/forcing/'; % the final slash is needed
forcenames = dir([forcingpath 'data*']);
ncells = length(forcenames);
precision = 4;

addpath(forcingpath)
tmpforc = dlmread(forcenames(1).name); % this requires ASCII forcing data. The Stehekin data are binary.
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

forc = NaN(nsteps, nvars, ncells);
for k=1:ncells
    forc(:,:,k) = dlmread(forcenames(k).name);
end

%% Time series

% Choose a grid cell. Eventually, will implement selection by lat/lon
k = 1;

% You need to know which column of forc(:,:,k) corresponds to which met.
% forcing variable.

% A vector of times from the VIC model run is generated in
% vicoutputworkflow.m. Load it here.

timevectorpath = '/Users/jschapMac/Desktop/Stehekin_4_2/results/vic/default';
addpath(timevectorpath)
t = load('timevector.mat');
t = t.timevector;

plot(t, forc(:,1,k))
title(['Precipitation' ' time series for grid cell ' num2str(k)])
xlabel('Time (days)'), ylabel('Precipitation (mm)')

%% Maps



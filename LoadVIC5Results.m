function [gridcells, fluxresults, snowresults] = LoadVIC5Results(headers, nfluxvars, nsnowvars)

% Loads results from VIC simulation
% Must be run from the directory containing the VIC outputs. 
% If the routing model is being run, this directory must contain the 
% routing model outputs, as well.
%
% INPUTS
% headers = number of header rows
% nvars = number of variables in the output files
%
% OUTPUTS
% gridcells = locations of the gridcells (lat/lon)
% fluxresults
% snowresults 

%% Load VIC flux results
fluxnames = dir('fluxes*');
ncells = length(fluxnames);
tmp = dlmread(fluxnames(1).name, ' ', headers, 0); % skip headers
tmp = tmp(:,1:nfluxvars);
nsteps = size(tmp,1);
fluxresults = NaN([nsteps, nfluxvars ,ncells]);

% read headers/variable names
fID = fopen(fluxnames(1).name);
fluxheadernames = textscan(fID, '%s', nfluxvars, 'HeaderLines',headers-1);
fclose(fID);

% Get grid cells locations
gridcells = cell(ncells, 1);
for k=1:ncells
    tmpstring = fluxnames(k).name;
    tmpstring = strrep(tmpstring,'-',''); % remove some characters bc Matlab cannot handle them
    tmpstring = strrep(tmpstring,'.','_');
    tmpstring = strrep(tmpstring,'fluxes_','');
    tmpstring = strcat('cell_', tmpstring);
    gridcells{k} = tmpstring;
end

for i=1:ncells
    fluxresults(:,:,i) = dlmread(fluxnames(i).name, ' ', [headers 0 nsteps+headers-1 nfluxvars-1]);  
end

%% Load VIC snow results

snownames = dir('snow_*');
ncells = length(snownames);
tmp = dlmread(snownames(1).name, ' ', headers, 0); % skip headers
tmp = tmp(:,1:nsnowvars);
nsteps = size(tmp,1);
snowresults = NaN([nsteps, nsnowvars ,ncells]);

% read headers/variable names
fID = fopen(snownames(1).name);
snowheadernames = textscan(fID, '%s', nsnowvars, 'HeaderLines',headers-1);
fclose(fID);

for i=1:ncells
    snowresults(:,:,i) = dlmread(snownames(i).name, ' ', [headers 0 nsteps+headers-1 nsnowvars-1]);  
end

end
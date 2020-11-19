function [gridcells, fluxresults, snowresults] = LoadVIC5Results(headers, nfluxvars, nsnowvars, resultsdir)

% TODO
%
% Modify to handle different (non WB mode) VIC modes
% Only loads 25 columns for EB mode, my EB mode results have 31 columns
% However, I ran with frozen soils, too
% 
%
% Loads results from VIC simulation
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
fluxnames = dir(fullfile(resultsdir, 'fluxes*'));
ncells = length(fluxnames);
tmp = dlmread(fullfile(resultsdir, fluxnames(1).name), '\t', headers, 0); % skip headers
tmp = tmp(:,1:nfluxvars);
nsteps = size(tmp,1);
fluxresults = NaN([nsteps, nfluxvars ,ncells]);

% read headers/variable names
fID = fopen(fullfile(resultsdir, fluxnames(1).name));
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
    fluxresults(:,:,i) = dlmread(fullfile(resultsdir, fluxnames(i).name), '\t', [headers 0 nsteps+headers-1 nfluxvars-1]);  
    
%     for p=4:length(fluxheadernames{1})
%         figure, plot(fluxresults(:,p,i)) 
%         title(fluxheadernames{1}{p})
%     end
    
end



%% Load VIC snow results

snownames = dir(fullfile(resultsdir, 'snow_*'));
ncells = length(snownames);
tmp = dlmread(fullfile(resultsdir, snownames(1).name), '\t', headers, 0); % skip headers
tmp = tmp(:,1:nsnowvars);
nsteps = size(tmp,1);
snowresults = NaN([nsteps, nsnowvars ,ncells]);

% read headers/variable names
fID = fopen(fullfile(resultsdir, snownames(1).name));
snowheadernames = textscan(fID, '%s', nsnowvars, 'HeaderLines',headers-1);
fclose(fID);

for i=1:ncells
    snowresults(:,:,i) = dlmread(fullfile(resultsdir, snownames(i).name), '\t', [headers 0 nsteps+headers-1 nsnowvars-1]);  
end

end
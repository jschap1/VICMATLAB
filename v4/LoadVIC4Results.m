function [varargout] = LoadVIC4Results(results_dir, varargin)

% Loads results from VIC simulation OR the routing model.
%
% INPUTS
% --Need to specify the filename for the routing outputs. 
% --The form is prefix.timestep or prefix.timestep_mm
% prefix = string
% units = string (mm or cfs)
% timestep = string (day, month, or year)
%
% OUTPUTS
% gridcells = locations of the gridcells (lat/lon)
% fluxresults
% snowresults 
% routresults 

%% Load routing results
if nargin>1
    prefix = varargin{1};
    units = varargin{2};
    timestep = varargin{3};
    switch units
        case 'cfs'
           routname = [prefix '.' timestep];
        case 'mm'
           routname = [prefix '.' timestep '_mm'];
        otherwise
            error('Not a valid routname')
    end
    routresults = dlmread(fullfile(results_dir, routname)); 
    varargout{1} = routresults;
else
% OR    
%% Load VIC results
fluxnames = dir(fullfile(results_dir, 'fluxes*'));
ncells = length(fluxnames);
tmp = dlmread(fullfile(results_dir, fluxnames(1).name));  
fluxresults = NaN([size(tmp),ncells]);

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
    fluxresults(:,:,i) = dlmread(fullfile(results_dir, fluxnames(i).name));  
end

snownames = dir(fullfile(results_dir, 'snow*'));
ncells = length(snownames);
tmp = dlmread(fullfile(results_dir, snownames(1).name));  
snowresults = NaN([size(tmp),ncells]);

for i=1:ncells
    snowresults(:,:,i) = dlmread(fullfile(results_dir, snownames(i).name));  
end

varargout{1} = gridcells;
varargout{2} = fluxresults;
varargout{3} = snowresults;

end

end
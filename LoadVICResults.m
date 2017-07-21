function [varargout] = LoadVICResults(varargin)

% Loads results from VIC simulation and the routing model (if specified).
% Must be run from the directory containing the VIC outputs.
%
% INPUTS
% rout = logical (1 if routing model, 0 otherwise)
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

fluxnames = dir('fluxes*');
ncells = length(fluxnames);
tmp = dlmread(fluxnames(1).name);  
fluxresults = NaN([size(tmp),ncells]);

% Get grid cells locations
gridcells = cell(ncells, 1);
for k=1:ncells
    tmpstring = fluxnames(k).name;
    tmpstring = strrep(tmpstring,'-',''); % remove some characters bc Matlab cannot handle them
    tmpstring = strrep(tmpstring,'.','_');
    % tmpstring = strrep(tmpstring,'fluxes_',''); % enable to make variable
    % names a bit shorter
    gridcells{k} = tmpstring;
end

for i=1:ncells
    fluxresults(:,:,i) = dlmread(fluxnames(i).name);  
end

snownames = dir('snow*');
ncells = length(snownames);
tmp = dlmread(snownames(1).name);  
snowresults = NaN([size(tmp),ncells]);

for i=1:ncells
    snowresults(:,:,i) = dlmread(snownames(i).name);  
end

if nargin
    prefix = varargin{1};
    units = varargin{2};
    timestep = varargin{3};
    switch units
        case 'cfs'
           routname = [prefix '.' timestep];
        case 'mm'
           routname = [prefix '.' timestep '_mm'];
    end
    routresults = dlmread(routname); 
    varargout{4} = routresults;
end

varargout{1} = gridcells;
varargout{2} = fluxresults;
varargout{3} = snowresults;

end
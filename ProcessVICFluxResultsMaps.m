function FLUXES = ProcessVICFluxResultsMaps(gridcells, FLUXES, precision, t)

% Adds maps to the FLUXES structure generated in ProcessVICFluxResults.
% Note: this will eventually be combined with ProcessVICFluxResults.
%
% INPUTS
% gridcells = locations of the gridcells (lat/lon). Output from LoadVICResults.
% Struct of VIC flux results, with time series only
% precision = number of decimal points of precision for the VIC flux output file names
% t = time index when you would like to make the map
%
% OUTPUTS
% Struct of VIC flux results, with maps and time series

ncells = length(fieldnames(FLUXES.ts));

% Extract numbers from within the strings in the cell array using regexp
gridcells_numbers =  regexprep(gridcells,'\D',''); % this replaces non-numeric digits with "".

% There is probably a way to implement this without a loop
lon = NaN(ncells,1);
lat = NaN(ncells,1);
for k = 1:ncells
    str1 = gridcells_numbers{k}(1:2);
    str2 = gridcells_numbers{k}(3:2 + precision);
    str3 = gridcells_numbers{k}(3 + precision:1 + 2*precision);
    str4 = gridcells_numbers{k}(2 + 2*precision:1 + 3*precision);
    lat(k) = str2double([str1 '.' str2]);
    lon(k) = str2double([str3 '.' str4]);
end

    % Indexing convention
    % i, x
    % j, y
    % t, time
    % p, varnames
    % k, ncells

varnames = FLUXES.ts.gridcell_fluxes_48_1875_120_6875.Properties.VariableNames;
cellnames = fieldnames(FLUXES.ts);

for p = 1:length(varnames)
    
    FLUXES.maps.(varnames{p}) = NaN(ncells,1);
    
    for k=1:ncells

        FLUXES.maps.(varnames{p})(k) = FLUXES.ts.(cellnames{k}).(varnames{p})(t);
        % Currently only takes the first soil moisture layer. Modify to
        % split up each soil moisture layer into its own variable.
        
    end
    
end

FLUXES.lat = lat;
FLUXES.lon = lon;

end
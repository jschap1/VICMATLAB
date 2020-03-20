function [lat, lon] = GetCoords(gridcells, precision)

% Extracts coordinates (numeric) from flux output file names. 
% Assumes northwestern hemisphere (positive latitudes, negative longitudes)
%
% INPUTS
% gridcells = locations of the gridcells (lat/lon). Output from LoadVICResults.
% precision = number of decimal points of precision for the VIC flux output file names
%
% OUTPUTS
% Coordinates of grid cell centroids, as a numeric array

ncells = length(gridcells);

% Extract numbers from within the strings in the cell array using regexp
gridcells_numbers =  regexprep(gridcells,'\D',''); % this replaces non-numeric digits with "".

% There is probably a way to implement this without a loop
% Also, the switch-case could definitely be rewritten/removed
lon = NaN(ncells,1);
lat = NaN(ncells,1);
for k = 1:ncells
    
    ndigits = size(gridcells_numbers{1},2);
    str1 = gridcells_numbers{k}(1:2);
    str2 = gridcells_numbers{k}(3:2 + precision);   
    switch ndigits
        case 12 % 2 digits for lat and 2 for lon (precision = 4)
            str3 = gridcells_numbers{k}(3 + precision:4 + precision);
            str4 = gridcells_numbers{k}(2*precision:end);
        case 13 % 2 digits for lat and 3 for lon (precision = 4)
            str3 = gridcells_numbers{k}(3 + precision:2*precision + 1);
            str4 = gridcells_numbers{k}(2*precision + 2:end);
        case 14 % 2 digits for lat and 2 for lon (precision = 5)
            str3 = gridcells_numbers{k}(3 + precision:2*precision - 1);
            str4 = gridcells_numbers{k}(2*precision:end);
        case 15 % 2 digits for lat and 3 for lon (precision = 5)
            str3 = gridcells_numbers{k}(3 + precision:2*precision);
            str4 = gridcells_numbers{k}(1 + 2*precision:end);
        otherwise
            error('Incorrect number of characters in gridcells string. Precision must be 4 or 5 digits.');
    end
        
    lat(k) = str2double([str1 '.' str2]);
    lon(k) = str2double([str3 '.' str4]);
    
end

% Convention for western hemisphere
% lon = -lon;
    % Indexing convention
    % i, x
    % j, y
    % k, ncells

end
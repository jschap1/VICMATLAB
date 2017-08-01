function [lat, lon] = GetCoords(gridcells, precision)

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
    % k, ncells

end
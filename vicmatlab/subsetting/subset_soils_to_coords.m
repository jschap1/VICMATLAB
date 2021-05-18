% Subset soils to coords
%
% Subset soil parameter file to a list of coordinates
% 2/28/2021 JRS

% coordinates = [lon, lat];

function soils_sub = subset_soils_to_coords(soils, coordinates, outname)

resolution = 1/16;
grid_decimal = 5;
outformat = '3l';

% Find lines in SPF that are closest to each line of the coordinate list
% In future, make match_coords.m to do this task

slat = soils(:,3);
slon = soils(:,4);

soils_sub = soils;
soils_sub(:,1) = 0;

ncoords = size(coordinates,1);
for k=1:ncoords
    lon1 = coordinates(k,1);
    lat1 = coordinates(k,2);
    sind = find((abs(lat1 - slat) <= resolution/2) & (abs(lon1 - slon) <= resolution/2));
    soils_sub(sind,1) = 1;
end

soils_sub(soils_sub(:,1) == 0,:) = []; % include this line to reduce file size

ncells_in_spf = size(soils_sub,1);
disp(['There are ' num2str(ncells_in_spf) 'cells in the subsetted soil parameter file'])

write_soils(grid_decimal, soils_sub, outname, outformat)

return

% Obtains a list of coordinates from a basin mask
% 
% The coordinates are all the points where the basin mask is not NaN.
%
% Example
% maskname = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif';
% coords = basin_mask2coordinate_list(maskname)

function extent1 = basin_mask2coordinate_list(maskname)

[basin_mask, R] = geotiffread2(maskname);

basin_mask = double(basin_mask);
basin_mask(basin_mask~=1) = NaN;

ncells = length(basin_mask(~isnan(basin_mask)));
[nx, ny] = size(basin_mask);

k = 1;
lon1 = zeros(ncells, 1);
lat1 = zeros(ncells, 1);
for i=1:nx
    for j=1:ny
        if ~isnan(basin_mask(i,j))
            [lat1(k), lon1(k)] = pix2latlon(R, i, j);
            k = k + 1;
        end
    end
end

extent1 = [lon1, lat1];

return
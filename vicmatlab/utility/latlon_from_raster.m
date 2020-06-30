% Lat lon from raster
%
% Get list of latlon coordinates from raster
% 5/23/2020 JRS

function [lon, lat] = latlon_from_raster(raster, R)

ncells = length(raster(~isnan(raster)));
[nx, ny] = size(raster);

k = 1;
lon = zeros(ncells, 1);
lat = zeros(ncells, 1);
for i=1:nx
    for j=1:ny
        if ~isnan(raster(i,j))
            [lat(k), lon(k)] = pix2latlon(R, i, j);
            k = k + 1;
        end
    end
end

end
% Rmat = georefobj2mat(R, 'LL');
% geotiffread2
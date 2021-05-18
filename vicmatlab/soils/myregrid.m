% My Regrid
%
% Regrids georeferenced, gridded data to a finer or coarser resolution
% Used in make_soil_file.m
%
% INPUTS
% original and target coordinates
% original resolution gridded data
% method = 'nearest' or 'linear' or 'cubic', etc.
%
% Example:
% lat = ncread('../HWSD/HWSD_1247/data/T_SAND.nc4', 'lat');
% lon = ncread('../HWSD/HWSD_1247/data/T_SAND.nc4', 'lon');
% target_lon = min(lon):0.0625:max(lon);
% target_lat = min(lat):0.0625:max(lat);
% [lons, lats] = ndgrid(lon, lat);
% [rglons, rglats] = ndgrid(target_lon, target_lat);
% t_sand_rg = myregrid(lons, lats, rglons, rglats, t_sand)';

function rg = myregrid(olons, olats, rglons, rglats, og, method)

    F = griddedInterpolant(olons', olats', og', method);
    rg = F(rglons, rglats); 
    
end
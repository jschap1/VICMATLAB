% Crop DEM
%
% Crops a DEM to a user-defined extent
%
% INPUTS
% demname = name of file containing DEM. Should be a GeoTIFF.
% extent = [minlon, maxlon, minlat, maxlat]
% outname = name of the output file (the cropped DEM)
%
% OUTPUTS
% cropelev = cropped DEM
% Rcrop = georeferencing matrix for the cropped DEM
%
% Sample inputs:
% demname = '/Volumes/HD2/MERIT/DEM/Merged_1_16/merged_merit_dem_1_16.tif';
% extent = [-121.5, -119, 37.5, 38.3];
% outname = '/Volumes/HD4/SWOTDA/Data/Tuolumne/dem.tif';

function [cropelev, Rcrop] = crop_dem(demname, extent, outname)

[dem, R] = geotiffread(demname);
dem = flipud(dem);
Rmat = makerefmat(R.LongitudeLimits(1), R.LatitudeLimits(1), R.CellExtentInLongitude, R.CellExtentInLatitude);
[lon, lat] = pixcenters(Rmat, size(dem));

xres = R.CellExtentInLongitude;
yres = R.CellExtentInLatitude;

disp(['Bounding box is [' num2str(extent, 3) ']'])

xmin = extent(1);
xmax = extent(2);
ymin = extent(3);
ymax = extent(4);

lat_range = [ymin ymax];
lon_range = [xmin xmax]; 

[rect, croplon, croplat] = make_cropping_rectangle(lon, lat, lon_range, lat_range, xres, yres);
cropelev = imcrop(dem, rect);

figure, subplot(1,2,1)
plotraster(lon, lat, dem, 'Original elevation (m)', 'Lon', 'Lat')
subplot(1,2,2)
plotraster(croplon, croplat, cropelev, 'Cropped elevation (m)', 'Lon', 'Lat')

Rcrop = makerefmat(min(croplon), min(croplat), xres, yres);
geotiffwrite(outname, cropelev, Rcrop)
disp(['Saved cropped DEM as ' outname])

return
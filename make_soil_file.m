function make_soil_file()
% Uses HWSD data to make a soil parameter file for VIC

% minlat, minlon, maxlat, maxlon, res

% Do some operations in GDAL to get the domain in a reasonable size and
% format, such as a geotiff file for the study area extent
% gdal_translate -a_srs epsg:4326 ./HWSD_RASTER/hwsd.bil ./HWSD_RASTER/hwsd_geog.tif

% Domain: [66.158333, 24.025, 82.45, 37.083333] % [xll, yll, xur, yur]
% gdal_translate

% load soil parameter data for the region of interest



return
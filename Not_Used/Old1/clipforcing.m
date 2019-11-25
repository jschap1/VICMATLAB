% Read a shapefile and clip the met. forcing file to the shapefile

% Load a shapefile
cd '/Applications/MATLAB_R2014b.app/toolbox/map/mapdata'
concord = shaperead('concord_hydro_area.shp');

% Load a forcing file (netCDF)
ncname = '/Users/jschapMac/Documents/HydrologyData/Livneh/MetNC/prec.1915.nc';
ncdisp(ncname)
prec = ncread(ncname, 'prec');
lat = ncread(ncname, 'lat');
lon = ncread(ncname, 'lon');

% THIS DOES NOT WORK RIGHT NOW, BUT THE PROCESS IS CORRECT. JUST NEED TO UP
% MY UNDERSTANDING OF THE MATLAB IMAGE PROCESSING AND MAPPING TOOLBOXES.
% Create a mask
[X,Y] = meshgrid(lon,lat);
mask = NaN(length(lon), length(lat));
for i=1:length(lon)
    for j=1:length(lat)
        mask(i,j) = inpolygon(X(i,j), Y(i,j), concord.X(i), concord.Y(j));
    end
end

% Clip to the mask
precmask = prec(mask);









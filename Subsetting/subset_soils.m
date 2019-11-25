% Clips soil parameter file to basin extent
%
% INPUTS
% soils = soil parameter file encompassing a region at least as large as the basin
% landmask = land mask used to generate the soil parameter file
% extent = bounding box or shapefile for subsetting
% outname = Name of output
% outformat = format of output ('3l'). See write_soils.m
% grid_decimal = precision used in forcing filenames
%
% OUTPUTS
% Clipped soil parameter file
%
% Updated 10/31/2019 to perform watershed masking

function soils_subset = subset_soils(soils, extent, outname, outformat, grid_decimal)

disp('Assuming resolution is 1/16 degrees');
resolution = 1/16;
maskflag = 0;

if ischar(extent)

    tmp1 = strsplit(extent, '.');
    extension = tmp1{2};

    if strcmp(extension, 'tif')

        disp('Subsetting to DEM');
        [mask1, Rsub] = geotiffread(extent);
        mask1 = flipud(mask1);
        mask1(mask1<0) = NaN;
        Rdemmat = georefobj2mat(Rsub, 'LL');

        maskflag = 1;
        disp('Applying the input raster as a mask')
        [lon_sub, lat_sub] = pixcenters(Rdemmat, size(mask1));

    elseif strcmp(extension, 'shp')

        disp('Subsetting to extent of shapefile');
        extent = shaperead(extent);
        shape_lons = extent.X(1:end-1)';
        shape_lats = extent.Y(1:end-1)';
        lon_sub = min(shape_lons):resolution:max(shape_lons);
        lat_sub = min(shape_lats):resolution:max(shape_lats);
        disp('maskflag not used')

    end

elseif isnumeric(extent)

    disp('Subsetting to bounding box');
    disp('Assuming resolution is 1/16 degrees');
    resolution = 1/16;
    lon_sub = extent(1):resolution:extent(2);
    lat_sub = extent(3):resolution:extent(4);
    disp('maskflag not used')

end

slat = soils(:,3);
slon = soils(:,4);
nlat = length(lat_sub);
nlon = length(lon_sub);

soils_subset = soils;
soils_subset(:,1) = 0;
for i=1:nlat
    for j=1:nlon
    % Get the index of the study area lat/lons that (nearly) matches the basin mask

    sind = find((abs(lat_sub(i) - slat) <= resolution/2) & (abs(lon_sub(j) - slon) <= resolution/2));
    soils_subset(sind,1) = 1;

% Turn off grid cells outside the mask area
    if maskflag
      if mask1(i,j) == 0 || isnan(mask1(i,j))
        soils_subset(sind,1) = 0;
      end
    end

    end
    disp(i)
end

soils_subset(soils_subset(:,1) == 0,:) = []; % include this line to reduce file size

ncells_in_spf = size(soils_subset,1);
disp(['There are ' num2str(ncells_in_spf) 'cells in the subsetted soil parameter file'])

write_soils(grid_decimal, soils_subset, outname, outformat)

return

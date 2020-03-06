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
% Modified 1/22/2020 JRS
% Does a somewhat better job subsetting to a basin mask/matching the extent of the
% basin mask

function soils_subset = subset_soils(soils, extent, outname, outformat, grid_decimal)

disp('Assuming resolution is 1/16 degrees');
resolution = 1/16;

if ischar(extent) 
    tmp1 = strsplit(extent, '.');
    extension = tmp1{2};
    if strcmp(extension, 'tif')
        disp('Subsetting to DEM');
        [dem, Rsub] = geotiffread(extent);
        dem = flipud(dem);
        Rdemmat = georefobj2mat(Rsub, 'LL');
        [lon_sub, lat_sub] = pixcenters(Rdemmat, size(dem));  
    elseif strcmp(extension, 'shp')
        disp('Subsetting to shapefile');
        extent = shaperead(extent);
        shape_lons = extent.X(1:end-1)';
        shape_lats = extent.Y(1:end-1)';
        lon_sub = min(shape_lons):resolution:max(shape_lons);
        lat_sub = min(shape_lats):resolution:max(shape_lats);
        
%         min_lon = min(shape_lons);
%         max_lon = max(shape_lons);
%         min_lat = min(shape_lats);
%         max_lat = max(shape_lats);
        
    end
end

if isnumeric(extent) % in case extent is bbox coordinates
    disp('Subsetting to bounding box');
    disp('Assuming resolution is 1/16 degrees');
    resolution = 1/16;
    lon_sub = extent(1):resolution:extent(2);
    lat_sub = extent(3):resolution:extent(4);
end

slat = soils(:,3);
slon = soils(:,4);
nlat = length(lat_sub);
nlon = length(lon_sub);

% lon = unique(slon);
% lat = unique(slat);
% lat_ind=find(lat>=min_lat & lat<=max_lat);
% lon_ind=find(lon>=min_lon & lon<=max_lon);
% 
% resolution = 1/16; 
% disp('Assuming resolution is 1/16 degrees')
% lat_ind1 = find(abs(lat(:,1) - min_lat) < resolution/2);
% lat_ind2 = find(abs(lat(:,1) - max_lat) < resolution/2);
% lat_ind = lat_ind1:lat_ind2;
% lon_ind1 = find(abs(lon(:,1) - min_lon) < resolution/2);
% lon_ind2 = find(abs(lon(:,1) - max_lon) < resolution/2);
% lon_ind = lon_ind1:lon_ind2;

soils_subset = soils;
soils_subset(:,1) = 0;

% sind = zeros(nlon*nlat,1);
% ind1 =1;
for i=1:nlat
    for j=1:nlon
        
        % Get the index of the study area lat/lons that (nearly) matches the basin mask
%       sind(ind1) = find((abs(lat_sub(i) - slat) <= resolution/2) & (abs(lon_sub(j) - slon) <= resolution/2));
        sind = find((abs(lat_sub(i) - slat) <= resolution/2) & (abs(lon_sub(j) - slon) <= resolution/2));
        
%         figure
%         plotraster(lon_sub, lat_sub, dem, 'Basin Mask', 'Lon', 'Lat')
        
%       Adding this to ensure that cells outside of the domain are not
%       included in the subsetted soil parameter file (JRS 1/22/2020)
        if dem(i,j) == 0
            continue
        end
        
        soils_subset(sind,1) = 1;
%       ind1=ind1+1;
    
    end
    disp(i)
end

soils_subset(soils_subset(:,1) == 0,:) = []; % include this line to reduce file size

ncells_in_spf = size(soils_subset,1);
disp(['There are ' num2str(ncells_in_spf) 'cells in the subsetted soil parameter file'])

write_soils(grid_decimal, soils_subset, outname, outformat)

return
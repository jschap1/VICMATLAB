% Clips soil parameter file to basin extent
%
% Adapted from vicinputworkflow 3/12/2019 JRS
%
% Made into a function, modified to take a DEM as input 8/5/2019 JRS
%
% INPUTS
% Soil parameter file encompassing a region at least as large as the basin
% List of coordinates in basin
% 
% OUTPUTS
% Clipped soil parameter file
%
% EXAMPLE
%
% % Load full soil parameter file
% soilpath = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_2';
% soilname = 'soils_3L_MERIT.txt'; 
% disp('Loading soil parameter file')
% soils = load(fullfile(soilpath, soilname));
% disp('Soil parameter file has been loaded')
% 
% outname = './Data/IRB/VIC/MiniDomain2/soils.SB';
% demname = './Data/IRB/VIC/MiniDomain2/MERIT_small_extent_1_16.tif';
% 
% % MERIT landmask
% [landmask, Rmask] = geotiffread('/Volumes/HD3/VICParametersGlobal/Global_1_16/landmask/merit_mask_1_16.tif');
% Rmask_mat = georefobj2mat(Rmask);
% [slon1, slat1] = pixcenters(Rmask_mat, size(landmask));
% 
% figure
% subplot(1,2,1)
% plotraster(slon1, slat1, landmask, 'MERIT landmask', 'Lon', 'Lat')
% subplot(1,2,2)
% plotraster(lon_dem, lat_dem, dem, 'Cropped DEM', 'Lon', 'Lat')
%
% soils_clipped = clip_soils(soils, './Data/IRB/VIC/MiniDomain2/MERIT_small_extent_1_16.tif' ...
%     , './Data/IRB/VIC/MiniDomain2/soils.SB', '3l')

function soils_clipped = clip_soils(soils, demname, outname, outformat)

grid_decimal = 5; % precision used in forcing filenames

[dem, Rdem] = geotiffread(demname);
dem = flipud(dem);
Rdemmat = georefobj2mat(Rdem);
[lon_dem, lat_dem] = pixcenters(Rdemmat, size(dem));

% xres = Rdem.CellExtentInLongitude;
% yres = Rdem.CellExtentInLatitude;
% if abs(xres - yres) < 1e5
%     resolution = xres;
% else
%     error('x and y resolution are not equal')
% end

[landmask, Rmask] = geotiffread('/Volumes/HD3/VICParametersGlobal/Global_1_16/landmask/merit_mask_1_16.tif');
landmask = flipud(landmask);
% landmask = logical(landmask);

Rmask_mat = georefobj2mat(Rmask);
[slon, slat] = pixcenters(Rmask_mat, size(landmask));

nlon = length(lon_dem);
nlat = length(lat_dem);

% Nlon = length(slon);
Nlat = length(slat);

% landmask_copy = landmask;
for i=1:nlat
    for j=1:nlon
        [row_ind, col_ind] = latlon2pix(Rmask, lat_dem(i), lon_dem(j));
        
        landmaskval = landmask(Nlat - round(row_ind), round(col_ind));
        
        if landmaskval > 0
            landmask(Nlat - round(row_ind), round(col_ind)) = 2;
        end
        
    end
end

% now, the landmask == 2 where the cropped DEM is
% turn off any soil rows where landmask ~=2
landcells = find(landmask);
landvalues = landmask(landcells);

% domaincells = find(landmask == 2);

soils_clipped = soils;
soils_clipped(landvalues == 2, 1) = 1;
soils_clipped(landvalues == 1, :) = []; % include this line to reduce file size

% soils_clipped(landvalues == 1, 1) = 0;
% blank_rows = find(soils_clipped(:,1) == 0);

ncells = sum(dem(:)>0); % number of land cells in the DEM
if size(soils_clipped, 1) < ncells
    disp('Some mask grid cells were not in the soil parameter file')
    disp([num2str(ncells - size(soils_clipped, 1)) ' to be exact'])
end

write_soils(grid_decimal, soils_clipped, outname, outformat)

return

% 
% landcells = find(landmask);
% N = length(landcells); % number of cells in global 1/16 deg. land mask
% 
% Rmask_mat = georefobj2mat(Rmask);
% [slon, slat] = pixcenters(Rmask_mat, size(landmask));
% 
% % slat = soils(:,3);
% % slon = soils(:,4);
% 
% % fulldomain = xyz2grid(slon, slat, ones(length(slon), 1));
% 
% % xrange = [min(slon), max(slon)];
% % yrange = [min(slat), max(slat)];
% 
% % Rfull = georefcells(yrange, xrange, resolution, resolution);
% % Rfull = georefcells(yrange, xrange, size(fulldomain));
% 
% % clip the soil parameter file
% soils_clipped = soils;
% soils_clipped(:,1) = 0; % turn off all cells in the soil parameter file
% 
% 
% 
% % %%%%%%
% % % Experimental
% % for i=1:nlon
% %     for j=1:nlat
% %         
% %         % find the closest entry in the soil parameter file
% %         
% %         [delta_lat, lat_ind] = min(abs(slat - lat_dem(j)));
% %         
% %         lon_ind = find(abs(slon - lon_dem(j)) <= 4*resolution/2);
% %         lat_ind = find(abs(slat - lat_dem(j)) <= 4*resolution/2);
% %         
% %         ismember(lon_ind, lat_ind)
% %         
% %         slon(lon_ind)
% %         
% %         [lon_dem(i), lat_dem(j)]
% %         
% %         ss=soils(lon_ind, 3);
% %         
% %         % There are simply no latitudes less than 25.28 deg. in the soil parameter file corresponding
% %         % to longitude 65 degrees in the cropped DEM. This probably has to
% %         % do with either the land mask or the CRS or both.
% %         
% %     end
% % end
% % %%%%%%%%
% 
% maskxyz = raster2xyz(lon_dem, lat_dem, dem);
% nancells = maskxyz(:,3) == -9999;
% maskxyz(nancells, 3) = 0;
% maskxyz(~nancells, 3) = 1;
% 
% masklon = maskxyz(maskxyz(:,3) == 1,1);
% masklat = maskxyz(maskxyz(:,3) == 1,2);
% ncells = length(masklon);
% 
% % turn on the grid cells in the soils file that match the grid cells in the
% % mask
% 
% for k=1:ncells
%     
%     nearby_lats = find(abs(slat-masklat(k)) < 0.5*resolution);
%     nearby_lons = find(abs(slon-masklon(k)) < 0.5*resolution);
%     
%     ind = find(ismember(nearby_lats, nearby_lons));
%     cellnumber = nearby_lats(ind);
% %     slat(cellnumber)
% %     slon(cellnumber)
%     soils_clipped(cellnumber,1) = 1;
%     
%     % report progres
%     if mod(k, 1e3) == 0
%         disp(k)
%     end
%     
% end
% 
% 
% 
% 
% 
% return
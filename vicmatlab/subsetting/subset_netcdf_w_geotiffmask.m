% Subsets a netCDF file using a geotiff mask
%
% Example: subsetting VIC Image Driver outputs to a particular study basin
% or subsetting NetCDF forcing files to a particular study basin. Based on
% vicinputworkflow.m code
%
% runoff_tuolumne = subset_netcdf_w_geotiffmask(runoff_livneh, livneh_lon, livneh_lat, dem_tuo.tif)
%
% Modified so it actually works. But only for one time step. Multi-temporal
% still needs to be implemented.

function [cropvar, cropped_lon, cropped_lat] = subset_netcdf_w_geotiffmask(nc_data, nc_lon, nc_lat, mask_lon, mask_lat, mask1)

    [~, ~, nt] = size(nc_data);
    [ny, nx] = size(mask1);
    
    % This could be made more efficient vvvvv
%     [mask1, R1, ~, ~] = geotiffread2(geotiff_mask);
%     R1mat = georefobj2mat(R1, 'LL');
%     [masklon1, masklat1] = pixcenters(R1mat, size(mask1));
%     [masklon, masklat, maskvals] = grid2xyz(masklon1', masklat1', mask1);
    
    nanmask = double(mask1);
    nanmask(nanmask==0) = NaN;

    figure, subplot(1,2,1)
    plotraster(nc_lon, nc_lat, nc_data, 'Input data', '', '')
    caxis([0, 300])
    subplot(1,2,2)
    plotraster(mask_lon, mask_lat, nanmask, 'Mask', '', '')

    minx = min(mask_lon);
    miny = min(mask_lat);
    maxx = max(mask_lon);
    maxy = max(mask_lat);
    extent = [minx, maxx, miny, maxy];
    demname = 'temp1.tif';
    outname = 'temp2.tif';
    xres = 1/16;
    yres = 1/16;
    R = makerefmat(min(nc_lon), min(nc_lat), xres, yres);
    
    geotiffwrite(demname, nc_data(:,:,1), R);
    
    % This does not give an exact match to the size of the original mask,
    % so it is not useful.
    [cropvar, Rcrop] = crop_dem(demname, extent, outname);
    
    [cropped_lon, cropped_lat] = pixcenters(Rcrop, size(cropvar));
    
    verbose = 1;
    if verbose
        disp('Ignore title of plot -- it was originally intended for cropping a DEM')
    end
    
    
    
    % demname = '/Volumes/HD2/MERIT/DEM/Merged_1_16/merged_merit_dem_1_16.tif';
% extent = [-121.5, -119, 37.5, 38.3];
% outname = '/Volumes/HD4/SWOTDA/Data/Tuolumne/dem.tif';
    
%     ncells = length(nanmask(:));
%     lat_ind = zeros(ncells,1);
%     lon_ind = zeros(ncells,1);
%     k = 1;
%     for i=1:nx
%         for j=1:ny
%             [~, lat_ind(k)] = min(abs(mask_lat(i) - nc_lat));
%             [~, lon_ind(k)] = min(abs(mask_lon(j) - nc_lon));
%             k = k + 1;
%         end
%     end
%     
%     B = NaN(nt, ncells); % the subsetted data
%     for k=1:ncells      
%         B(:,k) = nc_data(lon_ind(k),lat_ind(k), :);        
%     end
%     
%     x_sub = zeros(ncells, 1);
%     y_sub = zeros(ncells, 1);
%     for k=1:ncells
%         x_sub(k) = nc_lon(lon_ind(k));
%         y_sub(k) = nc_lat(lat_ind(k));
%     end
%     
%     B_map = xyz2grid(x_sub, y_sub, B');
%     figure, plotraster(mask_lon, mask_lat, B_map, '', '', '')
%     
%     % Need to account for mask
%     % Cells that are not in the mask should not be counted toward ncells
%     masklon(maskvals==0) = [];
%     masklat(maskvals==0) = [];
%     ncells = length(masklon);
%     
%     figure
%     plotraster(masklon, flipud(masklat), mask1, '', '','')
%     
%     figure
%     plotraster(nc_lon, nc_lat, flipud(nc_data), '', '','')
% 
% 
%     % The above only has to be done once
%     % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%     
%     B = NaN(nt, ncells); % the subsetted data
%     for k=1:ncells      
%         B(:,k) = nc_data(lon_ind(k),lat_ind(k), :);        
%     end

return
% Subsets a netCDF file using a geotiff mask
%
% Example: subsetting VIC Image Driver outputs to a particular study basin
% or subsetting NetCDF forcing files to a particular study basin. Based on
% vicinputworkflow.m code
%
% runoff_tuolumne = subset_netcdf_w_geotiffmask(runoff_livneh, livneh_lon, livneh_lat, dem_tuo.tif)

function [B, masklon, masklat] = crop_livneh(nc_data, nc_lon, nc_lat, geotiff_mask)

    [~, ~, nt] = size(nc_data);
    
    % This could be made more efficient vvvvv
    [mask1, R1] = geotiffread(geotiff_mask);
    R1mat = georefobj2mat(R1, 'LL');
    [masklon1, masklat1] = pixcenters(R1mat, size(mask1));
    [masklon, masklat, maskvals] = grid2xyz(masklon1', masklat1', mask1);
    
    
    % Need to account for mask
    % Cells that are not in the mask should not be counted toward ncells
    masklon(maskvals==0) = [];
    masklat(maskvals==0) = [];
    ncells = length(masklon);
    
%     figure
%     plotraster(masklon, flipud(masklat), mask1, '', '','')
    
%     figure
%     plotraster(nc_lon, nc_lat, flipud(nc_data), '', '','')

    lat_ind = zeros(ncells,1);
    lon_ind = zeros(ncells,1);
    for k=1:ncells
            [~, lat_ind(k)] = min(abs(masklat(k) - nc_lat));
            [~, lon_ind(k)] = min(abs(masklon(k) - nc_lon));
    end    
    % The above only has to be done once
    % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    B = NaN(nt, ncells); % the subsetted data
    for k=1:ncells      
        B(:,k) = nc_data(lon_ind(k),lat_ind(k), :);        
    end

return
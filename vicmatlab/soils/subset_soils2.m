% Clips soil parameter file to basin extent
%
% INPUTS
% soils: larger soil parameter file to subset from
% maskfile: basin mask to subset to. Must have same resolution as soils.
% outname: savename
% outformat:  livneh, 2l, or 3l
% grid_decimal: e.g. 5
% generate_tif: flag for whether or not to write tif files
% setup: e.g. '3L-no-org-frost-msds'
% resolution: e.g. 1/4 or 1/16
%
% OUTPUTS
% soils_subset
% soil_var_path

function [soils_subset, soil_var_path] = subset_soils2(soils, maskfile, ...
    outname, outformat, grid_decimal, generate_tif, setup, resolution)

[basinmask, Rsub, lon_sub, lat_sub] = geotiffread2(maskfile);

slat = soils(:,3);
slon = soils(:,4);
nlat = length(lat_sub);
nlon = length(lon_sub);

large_mask = xyz2grid(slon, slat, soils(:,1));
figure, plotraster(slon, slat, large_mask, 'mask');

soils_subset = soils;
soils_subset(:,1) = 0;

X = basin_mask2coordinate_list(basinmask, Rsub);
dt = delaunayTriangulation([slon, slat]);
[inds, d] = nearestNeighbor(dt, X); 

soils_subset(inds,1) = 1;
soils_subset(soils_subset(:,1) == 0,:) = []; % include this line to reduce file size
soil_var_path = fileparts(outname);

if length(unique(inds)) ~= length(inds)
    warning('The subset does not perfectly match the basin extent')
    warning('Consider adaptive measures, such as duplicating an adjacent cell')
    
end
% 
%     % there can be duplicates in inds, causing soils_subset to be 
%     % smaller than the mask. essentially, this is roundoff error, 
%     % but we need to do something about it.
% 
%     % Find duplicates
%     [v, w] = unique( inds, 'stable' );
%     duplicate_indices = setdiff( 1:numel(inds), w );
% 
%     % Adds extra lines to soil parameter file to make sure it perfectly
%     % matches the input basinmask
%     numDup = length(duplicate_indices);
% 
%     % Put in grid cell (-87.625, 41.875) Found this by comparing maps.
%     % This cell is not in the domain, so I'll copy the next cell to it.
%     soils_subset(end + numDup,:) = soils(inds(duplicate_indices),:);
%     soils_subset(end,3) = 41.875; % the missing lat 
%     soils_subset(end,4) = -87.625; % the missing lon
%     
% end

if generate_tif
    
    soil_var_path = fullfile(soil_var_path, 'tifs');
    mkdir(soil_var_path)
    varnames = get_soil_var_names(setup); % 3L-no-org-frost-msds
    lat_vect = soils_subset(:,3);
    lon_vect = soils_subset(:,4);
    
    if length(varnames)~=size(soils_subset,2)
        error('Check that the value for setup is correct')
    end
    
    chopped = 0;
    for k=1:length(varnames)
        svar = soils_subset(:,k);
        svar_map = xyz2grid(lon_vect, lat_vect, svar);
%         svar_map = flipud(svar_map);
%         
%         figure,plotraster(lon_vect, lat_vect, svar_map, '')
        
%         svar_map = xyz2grid(lon_vect, lat_vect, flipud(svar));
%         figure,plotraster(lon_vect, lat_vect, svar_map, '')
       
%%%%%%%------------------------------------------------
% unless these figures are exactly the same, this function will not work as
% expected
% figure, imagesc(basinmask), title('DEM')
% figure, imagesc(svar_map), title('Soil parameter')
%%%%%%%------------------------------------------------

        if size(basinmask,1) ~= size(svar_map,1) || size(basinmask,2) ~= size(svar_map,2)
            disp('Matrix sizes do not match. Chopping off empty rows.')
            basinmask = basinmask(:,2:end);
            Rsub = makerefmat(min(lon_vect), min(lat_vect), resolution, resolution);
            chopped = 1;
            disp('Check whether the soil parameter plots are right-side up')
            disp('May need to flip them over for VIC model to run the correct domain')
        end

        soutname = fullfile(soil_var_path, [varnames{k} '.tif']);
        
        if chopped % difference most likely relates to Rsub
            geotiffwrite(soutname, svar_map, Rsub);
        else
            geotiffwrite(soutname, flipud(svar_map), Rsub);
        end
    end
    disp(['Geotiff files saved to ', soil_var_path]);
end

% figure, plotraster(lon_vect, lat_vect, elev_map, 'Elevation', '', '')

ncells_in_spf = size(soils_subset,1);
disp(['There are ' num2str(ncells_in_spf) 'cells in the subsetted soil parameter file'])

write_soils(grid_decimal, soils_subset, outname, outformat)

return
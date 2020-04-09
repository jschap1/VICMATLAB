% Crop soil parameter tif
% 
% Crops a tif generated from a larger soil parameter file to match a
% smaller study area. Based on plot_dsmax.m
%
% Sample inputs
% basin_mask = '/Volumes/HD4/SWOTDA/Data/UpperMiss/umrb_mask.tif';
% outname = '/Volumes/HD4/SWOTDA/Data/UpperMiss/L15/elev.mat';
% tif_in = '/Volumes/HD3/VICParametersCONUS/soil_data/elev.tif';

% elev = crop_soil_parameter_tif(tif_in, basin_mask, 'elev.mat', 1);

function clipped_soil_parameter = crop_soil_parameter_tif(tif_in, basin_mask, outname, plotflag)

[ws, Rvg, lon1, lat1] = geotiffread2(tif_in);

if plotflag
    figure
    plotraster(lon1, lat1, ws, 'Plot 1', '', '')
end

% Mask out the part of dsmax that we want to keep

[mask1, R1] = geotiffread(basin_mask);
mask1 = flipud(mask1);
R1mat = georefobj2mat(R1, 'LL');
[masklon1, masklat1] = pixcenters(R1mat, size(mask1));
[masklon, masklat, ~] = grid2xyz(masklon1', masklat1', mask1);
ncells = length(masklon);

if plotflag
    figure
    plotraster(masklon1, masklat1, mask1, 'Plot 2', '', '')
end

lat_ind = zeros(ncells,1);
lon_ind = zeros(ncells,1);
for k=1:ncells
        [~, lat_ind(k)] = min(abs(masklat(k) - lat1));
        [~, lon_ind(k)] = min(abs(masklon(k) - lon1));
end    

B = NaN(ncells, 1); % the subsetted data
for k=1:ncells      
    B(k) = ws(lat_ind(k),lon_ind(k)); % swapping the indices makes it work...   
end

cropped_soil_parameter = xyz2grid(masklon, masklat, B);

clipmask = double(mask1);
clipmask(clipmask==0) = NaN;

clipped_soil_parameter = cropped_soil_parameter.*clipmask;

if plotflag
    figure
    plotraster(masklon, masklat, clipped_soil_parameter, 'Plot 3', '', '')
end

extension = outname(end-2:end);
if strcmp(extension, 'tif')
    geotiffwrite(outname, flipud(clipped_soil_parameter), R1)
    disp('Save tif')
elseif strcmp(extension, 'mat')
    save(outname, 'cropped_soil_parameter', 'masklon', 'masklat');
    disp(['Saved ' outname])
else
    disp('Please enter a valid file extension')
end

return



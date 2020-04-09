% Read variable from Geotiffs
%
% Returns the variable, loaded into Matlab, cropped to the basin extent,
% and averaged to monthly values. 
%
% INPUTS
% tifdir = directory where tifs are saved
% sz = size of the basin_mask (nlon, nlat)

function var_crop_mon = read_var_from_geotiffs(tifdir, sz)

tifnames = dir([tifdir '/*.tif']);
for tt=1:length(tifnames)
    if tt==1
        var_crop_mon = zeros(sz(1), sz(2), length(tifnames));
    end
    var_crop = geotiffread2(fullfile(tifdir, tifnames(tt).name));
    var_crop_mon(:,:,tt) = var_crop;
end

return
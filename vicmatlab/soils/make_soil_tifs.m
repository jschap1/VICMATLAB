% Make soil tifs
%
% Generates tif files from the soil parameter file
% Taken from subset_soils()
%
% Usage
% soils_subset = load('soil_parameter_file.txt');
% setup = 'livneh';
% soil_var_path = 'output_directory';
% [basin, Rsub, ~, ~] = geotiffread2('basinmask.tif');
% make_soil_tifs(soils_subset, setup, soil_var_path, Rsub);

function make_soil_tifs(soils_subset, setup, soil_var_path, Rsub)

varnames = get_soil_var_names(setup);
lat_vect = soils_subset(:,3);
lon_vect = soils_subset(:,4);

if length(varnames)~=size(soils_subset,2)
    error('Check that the value for setup is correct')
end

for k=1:length(varnames)
    svar = soils_subset(:,k);
    svar_map = xyz2grid(lon_vect, lat_vect, svar);       
    soutname = fullfile(soil_var_path, [varnames{k} '.tif']);
    geotiffwrite(soutname, flipud(svar_map), Rsub);
end
disp(['Geotiff files saved to ', soil_var_path]);

end
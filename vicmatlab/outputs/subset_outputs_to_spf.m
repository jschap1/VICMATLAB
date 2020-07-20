% Subset forcings
%
% Copy over ASCII outputs for a particular basin
%
% Example:
% subset_outputs_to_spf(soilfile, forcdir, 'Forcings_')
% soilfile = '/Volumes/HD3/SWOTDA/Calibration/Bhakra/soils_bhakra.txt';
% prefix = 'wb_';
% forcdir = '/Volumes/HD4/SWOTDA/Data/IRB/Downscaled_Forcings/ascii_1980-2019';
% newforcdir = '/Volumes/HD4/SWOTDA/Data/IRB/ascii_1980-2019';

function [soil_lon, soil_lat] = subset_outputs_to_spf(soilfile, forcdir, prefix, newforcdir)

mkdir(newforcdir)

soils = load(soilfile);
forcnames = dir([fullfile(forcdir, prefix) '*']);

[lon, lat] = get_coordinates_from_VIC_file(forcdir, prefix);

soil_lat = soils(:,3);
soil_lon = soils(:,4);

ncells = length(soil_lat);
T = delaunayn([lat, lon]);
disp(['There are ' num2str(ncells) ' grid cells'])
for j=1:ncells
    k = dsearchn([lat, lon],T,[soil_lat(j) soil_lon(j)]);
    origfile = fullfile(forcdir, forcnames(k).name);
    newfile = fullfile(newforcdir, forcnames(k).name);
    copyfile(origfile, newfile);
    if mod(j,100)==0
        disp(num2str(j))
    end
end

return

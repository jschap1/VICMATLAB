
basin_mask = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif';
% basin_mask = '/Volumes/HD4/SWOTDA/Data/UpperMiss/umrb_mask.tif';
% outname = '/Volumes/HD4/SWOTDA/Data/UpperMiss/L15/elev.mat';

tifnames = dir(['/Volumes/HD3/VICParametersCONUS/soil_data/' '*.tif']);
N = length(tifnames);

for k=1:N
    
    outname = ['/Volumes/HD4/SWOTDA/Data/Colorado/L15/' tifnames(k).name];
    clipped_soil_parameter = crop_soil_parameter_tif(fullfile('/Volumes/HD3/VICParametersCONUS/soil_data/', tifnames(k).name), basin_mask, outname, 0);
    
end
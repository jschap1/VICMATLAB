% Run make_elevband
%
% Executive file to make the elevation band file

coarse_res = 0.0625;
min_delta = 40; % m
min_fract = 0.05; 
numbands = 5;

soils = load('/Users/jschap/Documents/soils.SB');
soils_new = soils(:,1:4);
soils_new(:,5) = soils(:,22);
soils = soils_new;

outfile = '/Users/jschap/Documents/mysnowbands.txt';
dem_name = '/Volumes/HD3/SWOTDA/Data/MERIT/merged_merit_30as.tif';

elevband = make_elevband(soils, dem_name, coarse_res, numbands, min_delta, min_fract, outfile);

% Get just the band-average elevations
% elevs = elevband(:,7:11);
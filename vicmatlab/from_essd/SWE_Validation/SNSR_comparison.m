% Comparing the simulated snowpack in the Tuolumne basin to the Sierra
% Nevada Snow Reanalysis (SNSR) dataset
%
% 10/18/2019 JRS

fname = 'C:/Users/Jacob/Documents/My Documents/Research/Data/SNSR/SN_SWE_WY1999.h5';
info1 = h5info(fname);

lon_snsr = h5read(fname, '/lon');
lat_snsr = h5read(fname, '/lat');

h5disp(fname, '/SWE')
start = [1,1,1];
count = [6601,5701,1];
swe = h5read(fname, '/SWE', start, count);

% Re-scale the SWE to 1/16 degrees for comparison with VIC outputs
% Might be better to do this in gdal
% Also, crop to the Tuolumne basin extent

box_dir = 'C:\Users\Jacob\Box\Margulis_Research_Group\Jacob';
[tuo_dem, R_tuo] = geotiffread(fullfile(box_dir, 'dem.tif'));
tuo_dem = flipud(tuo_dem);

R_tuo_mat = makerefmat(R_tuo.LongitudeLimits(1), R_tuo.LatitudeLimits(1), ...
    R_tuo.CellExtentInLongitude, R_tuo.CellExtentInLatitude);
[lon_tuo, lat_tuo] = pixcenters(R_tuo_mat, size(tuo_dem));

figure, imagesc(lon_tuo, lat_tuo, tuo_dem)
set(gca, 'ydir', 'normal'); 
title('Tuolumne Basin Extent (m)')

% Crop SWE to Tuolumne basin extent and resample to 1/16 resolution
% Definitely a better job for GDAL than for Matlab.
% Modeled Tuolumne for calendar year 1999, so need maps for water years
% 1999 (Oct. 1998 - Sept. 1999) and 2000 (Oct. 1999 - Sept. 2000)
%
% Goal is to compare simulated SWE to reanalysis SWE.

% h5read(fname, '')
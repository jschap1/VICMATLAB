% Wrapper to create vegetation parameter file from MODIS data
%
% 9/6/2019 JRS

clearvars -except soils
addpath('/Users/jschap/Documents/Codes/VICMATLAB')

%% Read MODIS file and resample to 1/16 deg. resolution

modisfile = '/Volumes/HD3/MODIS/MCD12C1/MCD12C1.A2017001.006.2018257171411.hdf';

lc = hdfread(modisfile, 'MOD12C1', 'Fields', 'Land_Cover_Type_1_Percent');

lc = single(lc);

% Interpolate land cover data using nearest neighbor to a 1/16 degree resolution
ores = 0.05;
lon = -180+0.5*ores:ores:180-0.5*ores;
lat = -90+0.5*ores:ores:90-0.5*ores;
[lons, lats] = ndgrid(lon, lat);

tres = 0.0625;
target_lon = -180+tres/2:tres:180-tres/2;
target_lat = -90+tres/2:tres:90-tres/2;

[tlons, tlats] = ndgrid(target_lon, target_lat);

[~, ~, nclasses] = size(lc);
nlon = length(target_lon);
nlat = length(target_lat);
lc_rg = zeros(nlon, nlat, nclasses);

for k=1:nclasses
    lc_rg(:,:,k) = myregrid(lons, lats, tlons, tlats, lc(:,:,k)', 'nearest');
end

lc1 = permute(lc_rg, [2,1,3]);
lc2 = flipud(lc1);
figure, imagesc(target_lon, target_lat, lc2(:,:,17))
set(gca, 'ydir', 'normal')

R = makerefmat(target_lon(1), target_lat(1), tres, tres);
outdir = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Classic';
geotiffwrite(fullfile(outdir, 'open_water_lc.tif'), flipud(lc1), R)

% IGBP land cover classes (key)
fID = fopen('/Volumes/HD3/MODIS/IGBP_classes.txt');
LC.classes = textscan(fID, '%s');
fclose(fID);
LC.classes = LC.classes{1};
LC.class_numbers = 1:nclasses;
LC.IGBP_numbers = LC.class_numbers - 1;

i = 1;
figure
plotraster(target_lon, target_lat, lc2(:,:,i), LC.classes{i}, 'Lon', 'Lat')

% save the resampled vegetation cover fraction data as geotiffs
resampled_modis_dir = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Classic';
for k=1:nclasses
    fname = [LC.classes{k} '_fraction.tif'];
    geotiffwrite(fullfile(resampled_modis_dir, fname), lc2(:,:,k), R); 
end
disp(['Saved resampled vegetation cover fraction data to ', resampled_modis_dir])

%% Crop vegetation cover fraction data to MERIT landmask (in R)

% (Save location = cropped_modis_dir)

%% Load soil parameter file

soilfile = '/Volumes/HD3/VICParametersCONUS/vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt';
% soilfile = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_4/Classic/soils_3L_MERIT.txt';
disp('Loading soil parameter file')
soils = load(soilfile);
disp('Soil parameter file has been loaded')

soils = soils(:,1:5); % only the first five columns are needed

%% Write vegetation parameter file

% There is something going on with the land classification. The vegpar file
% that I wrote out 1/13/20 is only 256 MB, whereas the v1_4 vegpar file was
% 418 MB. Figure this out.

cropped_modis_dir = '/Volumes/HD3/VICParametersGlobal/VICGlobal/vegetation/Vegetation_Fractions/Cropped2MERIT';
% cropped_modis_dir = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Classic';
maskfile = '/Volumes/HD2/MERIT/DEM/Merged_1_16/merit_mask_1_16.tif';
savename = fullfile(outdir, 'global_vegetation_1_16_IGBP_d2.txt');
checkflag = 1;

addpath('/Users/jschap/Documents/Codes/VICMATLAB/Make_Vegetation')

root_data = write_vegparam(soils, cropped_modis_dir, maskfile, savename, checkflag);

save('/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_4/Data/root_data_3.mat', 'root_data', '-v7.3');

% load('/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_4/Data/root_data.mat')

%% Read vegetation parameter file (to check)

indir = '/Volumes/HD3/VICParametersCONUS';
vegfile = fullfile(indir, 'vic.veg.0625.new.cal.adj.can');

% indir = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Classic';
% vegfile = fullfile(indir, '/global_vegetation_1_16_IGBP.txt');
ncells = size(soils, 1);

% name to use to save the outputs of read_vegparam (formatted vegetation
% parameter data)
savename = '/Volumes/HD3/VICParametersCONUS/Vegetation/vegparamtable.mat';
% savename = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Data/Vegetation/vegparamtable_IGBP.mat';

% ncells = 139589;
% Note: read_vegparam saves nvegtable, etc. under savename
% [nvegtable, vegparamtable, latlontable, LC] = read_vegparam(vegfile, soils, ncells, savename);
[nvegtable, vegparamtable, latlontable, LC] = read_vegparam_w_LAI(vegfile, soils, ncells, savename);

% dat = load(savename);
% ind2 = find(dat.vegparamtable.Woody_Savannas(:,2)>0);
% dat.vegparamtable.Woody_Savannas(ind2(1),:)
% 
% % Percent of landmass with each land cover
% length(find(dat.vegparamtable.Deciduous_Broadleaf(:,2)))/ncells*100
% length(find(dat.vegparamtable.Deciduous_Needleleaf(:,2)))/ncells*100
% 
% dat.vegparamtable.Evergreen_Broadleaf(1:10,:)

load(savename);

%% Make Geotiffs

veglib = '/Volumes/HD3/VICParametersCONUS/vic_veglib_nohead.txt';


% Since we're using the northern hemisphere vegetation library, the
% southern hemisphere portions of the maps should be thrown out...

outdir = '/Volumes/HD3/VICParametersCONUS/Vegetation/';
plotflag = 0;
saveflag = 1;

% clearvars -except nvegtable vegparamtable latlontable veglib LC outdir plotflag saveflag

VEGPARAMS = convert_veg_parameters_v3(nvegtable, vegparamtable, latlontable, veglib, LC, outdir, plotflag, saveflag);

VEGPARAMS = load('/Volumes/HD3/VICParametersCONUS/Vegetation/VEGPARAM.mat');
VEGPARAMS = VEGPARAMS.VEGPARAM;

% Remove unnecessary parts of VEGPARAMS
% rmfield(VEGPARAMS, 'sdfdf')
% There are none :(
% It is big.

% savename = '/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_4/Data/VEGPARAMS.mat';
% save(savename, 'VEGPARAMS', '-v7.3')

% This uses a stupid amount of memory.

%% Plot root fractions

everneedle_rf1 = flipud(geotiffread(fullfile(outdir, 'Evergreen_Needleleaf_rf1.tif')));
figure

latmax = max(latlontable(:,2));
latmin = min(latlontable(:,2));
lonmax = max(latlontable(:,3));
lonmin = min(latlontable(:,3));

lonlim = [lonmin, lonmax];
latlim = [latmin, latmax];

plotraster(lonlim, latlim, everneedle_rf1, 'Root fraction (layer 1)', '', '')

%% Plot LAI

% Can run independently, after running the above code

addpath('/Volumes/HD3/VICParametersGlobal/Global_1_16/vegetation/Vegetation_Fractions/Codes')

savename = '/Volumes/HD3/VICParametersGlobal/Global_1_16/vegetation/Vegetation_Fractions/VEGPARAMS.mat';
load(savename)

savename = '/Volumes/HD3/VICParametersGlobal/Global_1_16/vegetation/Vegetation_Fractions/global_vegetation_params_IGBP.mat';
load(savename)

outdir = '/Volumes/HD3/VICParametersGlobal/Global_1_16/vegetation/Vegetation_Fractions/tiffs';
plot_lai(VEGPARAMS.LAI, VEGPARAMS.lon, VEGPARAMS.lat, LC, VEGPARAMS.R, outdir)

%% Check if things are OK

MAPS.Jan.Barren = xyz2grid(VEGPARAMS.lon, VEGPARAMS.lat, VEGPARAMS.LAI.Jan.Barren);
MAPS.Jul.Barren = xyz2grid(VEGPARAMS.lon, VEGPARAMS.lat, VEGPARAMS.LAI.Jul.Barren);

MAPS.Jan.Deciduous_Broadleaf = xyz2grid(VEGPARAMS.lon, VEGPARAMS.lat, VEGPARAMS.LAI.Jan.Deciduous_Broadleaf);
MAPS.Jul.Deciduous_Broadleaf = xyz2grid(VEGPARAMS.lon, VEGPARAMS.lat, VEGPARAMS.LAI.Jul.Deciduous_Broadleaf);

minlat = min(VEGPARAMS.lat);
minlon = min(VEGPARAMS.lon);
maxlat = max(VEGPARAMS.lat);
maxlon = max(VEGPARAMS.lon);
latlim = [minlat, maxlat];
lonlim = [minlon, maxlon];

figure
subplot(2,2,1)
plotraster(lonlim, latlim, MAPS.Jan.Barren, 'Barren LAI, January', 'Lon', 'Lat')
subplot(2,2,2)
plotraster(lonlim, latlim, MAPS.Jul.Barren, 'Barren LAI, July', 'Lon', 'Lat')
subplot(2,2,3)
plotraster(lonlim, latlim, MAPS.Jan.Deciduous_Broadleaf, 'DecidBroad LAI, January', 'Lon', 'Lat')
subplot(2,2,4)
plotraster(lonlim, latlim, MAPS.Jul.Deciduous_Broadleaf, 'DecidBroad LAI, July', 'Lon', 'Lat')



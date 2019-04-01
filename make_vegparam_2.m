% Make Vegetation Parameter File
%
% Script that uses data to make a vegetation parameter file for VIC 4 or
% for VIC 5 classic. Sort of like an inverted version of
% load_veg_parameters.m

function make_vegparam_2(soils, landmask, savename)

%% Load input data

% Load MODIS land cover data

soils = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/soils/second_attempt/global_soils_1_16.txt');

lc = hdfread('/Volumes/HD3/MODIS/MCD12C1/MCD12C1.A2017001.006.2018257171411.hdf', ...
    'MOD12C1', 'Fields', 'Majority_Land_Cover_Type_2');
lc = single(lc);

% Interpolate land cover data using nearest neighbor to a 1/16 degree resolution

ores = 0.05;
lon = -180+0.5*ores:ores:180-0.5*ores;
lat = -90+0.5*ores:ores:90-0.5*ores;
[lons, lats] = ndgrid(lon, lat);

tres = 0.0625;
target_lon = min(lon):tres:max(lon);
target_lat = min(lat):tres:max(lat);
[tlons, tlats] = ndgrid(target_lon, target_lat);
lc_rg = myregrid(lons, lats, tlons, tlats, lc', 'nearest');

lc1 = lc_rg';
lc2 = flipud(lc1);
figure, imagesc(target_lon, target_lat, lc2)
set(gca, 'ydir', 'normal')

R = makerefmat(target_lon(1), target_lat(1), tres, tres);
geotiffwrite('./Data/Global/landcover_umd.tif', lc2, R); % save intermediate result

%% Re-assign land cover types as needed

% land cover types key
% IGBPnames = {'water', 'EvergreenNeedleleafForest', 'EvergreenBroadleafForest', ...
%     'DeciduousNeedleleafForest', 'DeciduousBroadleafForest', 'MixedForest', ...
%     'ClosedShrublands', 'OpenShrublands', 'WoodySavannas', 'Savannas', 'Grasslands', ...
%     'PermanentWetlands', 'Croplands', 'Urban', 'CroplandNaturalVeg', 'SnowIce', 'Barren'};

UMDnames = {'Water','EvergreenNeedleleaf', 'EvergreenBroadleaf', ...
    'DeciduousNeedleleaf', 'DeciduousBroadleaf', 'MixedCover', ...
    'ClosedShrublands', 'OpenShrublands', 'WoodySavannas', 'Savannas', ...
    'Grasslands', 'PermanentWetlands','Croplands', 'Urban', ...
    'Croplandnaturalveg', 'nonvegetated'};
UMDnumbers = 0:15;

lcnames = {'EvergreenNeedleleaf', 'EvergreenBroadleaf', ...
    'DeciduousNeedleleaf', 'DeciduousBroadleaf', 'MixedCover', ...
    'Woodland', 'WoodedGrasslands', 'ClosedShrublands', 'OpenShrublands', 'Grasslands', ...
    'CroplandCorn'};
lcnumbers = 1:11;

landmask = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/landmask/landmask_hwsd_1_16.mat');
landmask = landmask.landmask;

vegtype = lc2(landmask); % land cover for HWSD land cells

% Reclassify from UMD to GLDAS LC numbering system
vegcopy = vegtype;

vegcopy(vegtype == 0) = 0; % re-assign water to bare soil
vegcopy(vegtype == 11) = 0; % re-assign wetlands to bare soil
vegcopy(vegtype == 13) = 0; % re-assign urban to bare soil
vegcopy(vegtype == 15) = 0; % re-assign nonvegetated to bare soil

vegcopy(vegtype == 6) = 8; % ClosedShrublands
vegcopy(vegtype == 7) = 9; % OpenShrublands

vegcopy(vegtype == 8) = 6; % WoodySavannas -> Woodland
vegcopy(vegtype == 9) = 7; % Savannas -> WoodedGrasslands

vegcopy(vegtype == 12) = 11; % Croplands
vegcopy(vegtype == 14) = 11; % CroplandNaturalVeg -> Croplands

vegtype = vegcopy;

ntypes = length(unique(vegtype)); % number of land cover types

%%

cellID = soils(:,2);
ncells = length(cellID);
nveg = ones(ncells,1); % always 1 for the current setup

% Initialize variables
cv = ones(ncells,1);
rootdepth = zeros(ncells, 2);
rootfract = zeros(ncells, 2);
lai = zeros(ncells, 12);

% Assign parameter values from lookup table
lookuptable = dlmread('/Volumes/HD3/VICParametersGlobal/Global_1_16/vegetation/vegpars_lut.txt', '\t', 1, 1);

% Note: vegtypes other than the 11 listed here are classified implicitly at bare soil (barren, urban, and water)

for cl=1:(ntypes-1)
    n_of_type = sum(vegtype == cl);
    row = find(lookuptable(:,1)==cl); % get row of lookup table
    rootdepth(vegtype == cl,:) = repmat(lookuptable(row, 2:3), n_of_type, 1);
    rootfract(vegtype == cl,:) = repmat(lookuptable(row, 4:5), n_of_type, 1);
    lai(vegtype == cl,:) = repmat(lookuptable(row, 6:17), n_of_type, 1);
end

%% Write the vegetation parameter file

fID = fopen(savename, 'w');

current_cellID = 1;

while current_cellID<=ncells

    current_nveg = nveg(current_cellID);

    if vegtype(current_cellID) == 0
        fmt = '%d %d\n';
        fprintf(fID, fmt, [current_cellID, 0]); % "bare soil"
        current_cellID = current_cellID + 1;
        continue
    else
        fmt = '%d %d\n';
        fprintf(fID, fmt, [current_cellID, current_nveg]);
    end
            
    % For each vegetation class, write parameters, then monthly LAI

    current_cv = cv(current_cellID, :);
    current_rd = rootdepth(current_cellID, :);
    current_rf = rootfract(current_cellID, :);

    fmt = '%d %4.3f %0.2f %0.2f %0.2f %0.2f\n';
    fprintf(fID, fmt, [vegtype(current_cellID) current_cv, current_rd, current_rf]);

    current_lai = lai(current_cellID, :);

    fmt = '%4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f\n';
    fprintf(fID, fmt, current_lai);
                
    current_cellID = current_cellID + 1;
        
    % show progress
    if mod(current_cellID, 1e5)==0
        disp(round(current_cellID/ncells*100))
    end
    
end

fclose(fID);

return


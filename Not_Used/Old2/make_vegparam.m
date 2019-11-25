% Make Vegetation Parameter File
%
% Script that uses data to make a vegetation parameter file for VIC 4 or
% for VIC 5 classic. Sort of like an inverted version of
% load_veg_parameters.m

function make_vegparam(soils, landmask, savename)

%% Load input data

% now, I just need to create a structure like VEGPAR from a high-resolution
% data source, and I can write the output to a VIC vegetation parameter
% file and plot it :)

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

% Make VEGPAR object
% VEGPAR_coarse = VEGPAR; % saving the old one (remove this line later)

% IGBP land cover types key
% IGBPnames = {'water', 'EvergreenNeedleleafForest', 'EvergreenBroadleafForest', ...
%     'DeciduousNeedleleafForest', 'DeciduousBroadleafForest', 'MixedForest', ...
%     'ClosedShrublands', 'OpenShrublands', 'WoodySavannas', 'Savannas', 'Grasslands', ...
%     'PermanentWetlands', 'Croplands', 'Urban', 'CroplandNaturalVeg', 'SnowIce', 'Barren'};

UMDnames = {'EvergreenNeedleleaf', 'EvergreenBroadleaf', ...
    'DeciduousNeedleleaf', 'DeciduousBroadleaf', 'MixedCover', ...
    'Woodland', 'WoodedGrasslands', 'ClosedShrublands', 'OpenShrublands', 'Grasslands', ...
    'CroplandCorn'};

ntypes = length(UMDnames); % number of land cover types

landmask = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/landmask/landmask_hwsd_1_16.mat');
landmask = landmask.landmask;
landcells = find(landmask);
ncells = length(landcells);

vegtype = lc2(landmask); % land cover for HWSD land cells

cellID = soils(:,2);
nveg = ones(ncells,1); % always 1 for the current setup

% Initialize variables
cv = ones(ncells,1);
rootdepth = zeros(ncells, 2);
rootfract = zeros(ncells, 2);
lai = zeros(ncells, 12);

% Assign parameter values from lookup table
lookuptable = dlmread('/Volumes/HD3/VICParametersGlobal/Global_1_16/vegetation/vegpars_lut.txt', '\t', 1, 1);

for cl=1:ntypes
    rootdepth(vegtype == cl,:) = lookuptable(cl, 2:3);
    rootfract(vegtype == cl,:) = lookuptable(cl, 4:5);
    lai(vegtype == cl,:) = lookuptable(cl, 6:17);
end

% % Create VEGPAR structure for writing to file (or better yet, don't use
% % VEGPAR structure)
% ncells = size(soils, 1);
% VEGPAR(ncells) = struct();
% for k=1:ncells % loop is kind of slow (~15 minutes). Also, the VEGPAR file is HUGE.
%     
%     VEGPAR(k).cellID = soils(k,2);
%     VEGPAR(k).nveg = 1; % for now, just using the majority value for each pixel
%     
% %     lc2(k) % vegetation class number for this cell
%     
%     VEGPAR(k).(lcnames{lc2(k)}).cv = 1;
%     VEGPAR(k).(lcnames{lc2(k)}).rootdepth = [0.3 0.3]; % should have lookup tables for these parameters...
%     VEGPAR(k).(lcnames{lc2(k)}).rootfract = [0.7 0.7];
%     VEGPAR(k).(lcnames{lc2(k)}).LAI = [0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7];
%         
%     if mod(k, 1e4)==0
%         disp(num2str(k))
%     end
%     
% end


%% Scrap 

% Assign vegetation parameter values to each land cover tile
% Use lookup table based on the UMD land cover classification
% Re-process the land cover data, but use UMD classifications, not IGBP

% only want the data that covers the study area...
% other pixels can be set to NULL values
% it would be nice to start off with a global dataset and then subset it

% download the HMA static data to PC
% fill in NA values where HMA processing is not done yet
% do the rest of the HMA processing

% VEGPAR % using the same structure as load_veg_parameters.m

% HMA project data example ---------------------------------------------
% A = load('/Users/jschap/N48_0E91_0/DATA_STATIC/static_agg_16_N48_0E91_0_SRTM.mat');
% A.lon
% A.lat
% A.agg_data.landcover
% figure, imagesc(A.agg_data.lon, A.agg_data.lat, A.agg_data.landcover)
% set(gca, 'ydir', 'normal')
% ----------------------------------------------------------------------

%% Write the vegetation parameter file

fID = fopen(savename, 'w');

current_cellID = 1;

while current_cellID<=ncells

    current_nveg = VEGPAR(current_cellID).nveg;

    fmt = '%d %d\n';
    fprintf(fID, fmt, [current_cellID, current_nveg]);
    
    if current_nveg==0
        current_cellID = current_cellID + 1;
        continue
    end
    
    % get the vegetation classes for this grid cell by checking where the
    % empty values are
    fn = fieldnames(VEGPAR);
    vegnames = fn(3:end);
    veg_in_gridcell = zeros(length(vegnames),1);
    for k=1:length(vegnames)
        veg_in_gridcell(k) = ~isempty(VEGPAR(current_cellID).(vegnames{k}).cv);
    end
    classnums = find(veg_in_gridcell);
    
    % For each vegetation class, write parameters, then monthly LAI
    for cl=[classnums']
        current_cv = VEGPAR(current_cellID).(vegnames{cl}).cv;
        current_rd = VEGPAR(current_cellID).(vegnames{cl}).rootdepth;
        current_rf = VEGPAR(current_cellID).(vegnames{cl}).rootfract;
        
        fmt = '%d %4.3f %0.2f %0.2f %0.2f %0.2f\n';
        fprintf(fID, fmt, [cl current_cv, current_rd, current_rf]);
        
        current_lai = VEGPAR(current_cellID).(vegnames{cl}).LAI;
        
        fmt = '%4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f\n';
        fprintf(fID, fmt, current_lai);
        
    end
        
    current_cellID = current_cellID + 1;
        
end

fclose(fID);

return


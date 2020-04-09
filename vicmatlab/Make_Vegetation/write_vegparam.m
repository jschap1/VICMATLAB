% Make Vegetation Parameter File
%
% Script that uses data to make a vegetation parameter file for VIC 4 or
% for VIC 5 classic. Sort of like an inverted version of
% load_veg_parameters.m
%
% Update: Version 4 allows partial vegetation cover fractions from
% different land covers in the same tile. Source dataset is MOD12Q1, the
% 500 m MODIS land cover dataset.
%
% INPUTS
% soil parameter file
% modis_dir = directory containing vegetation cover data from MODIS
% MOD12C1, resampled to 1/16 degree resolution and cropped to the extent of the MERIT mask
% landmask from MERIT, corresponding to soil parameter file
% filename to save vegetation parameter file
%
% OUTPUTS
% Vegetation parameter file
%
% Example:
% root_data = write_vegparam(soils, cropped_modis_dir, maskfile, savename, checkflag);

function LC = write_vegparam(soils, modis_dir, meritmaskfile, savename, checkflag)

outdir = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Classic';

%% Prepare to write vegetation parameter file

% Get cell IDs and coordinates from soil parameter file
soils = soils(:,1:5);
ncells = size(soils, 1);
cellID = soils(:,2);

% Load landmask from MERIT
[landmask, R_merit] = geotiffread(meritmaskfile);
landmask = flipud(logical(landmask));

% Load MODIS land cover data
lc_names = dir(fullfile(modis_dir, '*fraction.tif'));
lc_names_full = lc_names;
nclasses = length(lc_names);
for i=1:nclasses
    lc_names(i).name = lc_names(i).name(1:end-13);
end

% Make vegetation cover fraction array
lc_all = struct();
vegtype = zeros(ncells, nclasses);
for i=1:nclasses
    lc_all.(lc_names(i).name) = geotiffread(fullfile(modis_dir, lc_names_full(i).name));
    lc_all.(lc_names(i).name) = flipud(lc_all.(lc_names(i).name));
    vegtype(:, i) = lc_all.(lc_names(i).name)(landmask);
    disp(i)
end
vegtype = vegtype/100; % convert from percentage to fraction

alphabetical_order = [16, 8, 7, 6, 5, 10, 2, 11, 17, 13, 9, 12, 4, 15, 3, 14, 1];
vegtype = vegtype(:, alphabetical_order); % this should fix the issue with IGBP numbers/order
% note that this assumes the lc_names are in alphabetical order

% Get number of vegetation covers in each grid cell
nveg = zeros(ncells, 1);
for k=1:ncells
    nveg(k) = sum(vegtype(k,:) > 0);
    if mod(k, 5e5) == 0
        disp(k)
    end
end

% Make cell arrays to store vegetation covers (names and percentages) for
% each grid cell
cover_fraction = cell(ncells, 1);
vgn = cell(ncells, 1); % (vegtype number)
for k=1:ncells 
    vgn{k} = find(vegtype(k,:));
    cover_fraction{k} = vegtype(k, vgn{k});
    if mod(k, 5e5) == 0
        disp(k)
    end
end

% save('/Volumes/HD3/VICParametersGlobal/Global_1_16/vegetation/Vegetation_Fractions/mats/vgn.mat', 'vgn');
% outdir = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Classic';
save(fullfile(outdir, 'vgn.mat'), 'vgn');

%% Calculate root depth and root fraction values

% Set root depth and root fraction values based on the land cover type
% Using the method of Xubin Zeng (2001, JHM)

% Define land cover names and IGBP codes in alphabetical order
LC.land_cover = lc_all;
LC.names = cell(nclasses, 1);
for i=1:nclasses
    LC.names{i} = lc_names(i).name;
end
LC.names = LC.names(alphabetical_order); 
% now the names are in order of the IGBP numbers, so LC.names(1) is water 
% bodies, and LC.names(17) is barren

LC.alphabetical_order = alphabetical_order; % for accounting purposes
LC.IGBP_Number = [17, 7, 13, 15, 5, 4, 3, 2, 11, 6, 8, 1, 12, 9, 16, 14, 10];

% Load parameters from Zeng (2001)
zeng_table = load('/Users/jschap/Documents/Codes/VICMATLAB/Make_Vegetation/zeng_table.txt');
zeng_table = table(zeng_table(:,1), LC.names, zeng_table(:,2), zeng_table(:,3), zeng_table(:,4));
zeng_table.Properties.VariableNames = {'IGBP_Number','IGBP_Name', 'a','b','dr'};

% Define cumulative root fraction formula from Zeng (2001)
cum_root_fract = @(a, b, d) 1 - 0.5.*(exp(-a.*d) + exp(-b.*d));

% Calculate root fraction at several depths: 0.1, 0.7, and maximum (dr)

% Test
% i = 1; % land cover 1, aka Barren
% d = 0:0.1:zeng_table.dr(LC.IGBP_Number(i));
% Y = cum_root_fract(zeng_table.a(LC.IGBP_Number(i)), zeng_table.b(LC.IGBP_Number(i)), d);
% figure, plot(d,Y), xlabel('Depth'), ylabel('Cumulative root fraction')

for i=1:nclasses
    LC.root_depth.(LC.names{i}) = zeros(1, 3);
    LC.root_fract.(LC.names{i}) = zeros(1, 3);
end

for i=1:nclasses
    disp(['Calculating root parameters for ' LC.names{i} ' ' zeng_table.IGBP_Name{i}])
    
    d = [0.1, 0.7 - 0.1, zeng_table.dr(i) - 0.7];
%     d = [0.1, 0.7, zeng_table.dr(i)];
    
    Y = cum_root_fract(zeng_table.a(i), zeng_table.b(i), d);
    
    % Fix the cumulative root fraction at dr to be 1
    Y(3) = 1;
    
    LC.root_fract.(LC.names{i}) = [Y(1), Y(2) - Y(1), Y(3) - Y(2)];
    LC.root_depth.(LC.names{i}) = d;
end

%% Perform optional checks

if checkflag

    % Plot a sample vegetation cover map
    ii = 1;
    lc = geotiffread(fullfile(modis_dir, lc_names_full(ii).name));
    lc = flipud(lc);    
    minlat = R_merit.LatitudeLimits(1);
    maxlat = R_merit.LatitudeLimits(2);
    minlon = R_merit.LongitudeLimits(1);
    maxlon = R_merit.LongitudeLimits(2);  
    figure
    plotraster([minlon maxlon], [minlat maxlat], lc, lc_names(ii).name, 'Lon', 'Lat')

    % Plot land and water    
    figure, 
    subplot(2,1,1)
    imagesc([minlon maxlon], [minlat maxlat], lc_all.Water_Bodies) 
    title('MODIS Water Bodies')
    set(gca, 'ydir', 'normal')

    subplot(2,1,2)
    imagesc([minlon maxlon], [minlat maxlat], landmask)
    title('MERIT Land Mask')
    set(gca, 'ydir', 'normal')
    
end

%% Write vegetation parameter file

fID = fopen(savename, 'w');

for k=1:ncells

    fmt = '%d %d\n';
    fprintf(fID, fmt, [cellID(k), nveg(k)]);
  
    % For each vegetation class, write parameters, then (optionally) monthly LAI

    for ii=1:nveg(k)
        current_vgn = vgn{k}(ii);
        current_cf = cover_fraction{k}(ii);
        current_rd = LC.root_depth.(LC.names{current_vgn});
        current_rf = LC.root_fract.(LC.names{current_vgn});
        fmt = '%d %4.3f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f\n';
        
%         vegpars = [current_vgn current_cf, current_rd, current_rf];
        vegpars = [current_vgn current_cf, current_rd(1), current_rf(1), current_rd(2), current_rf(2), current_rd(3), current_rf(3)];
        
        fprintf(fID, fmt, vegpars);
    end
           
    % show progress
    if mod(k, 1e4)==0
        disp(round(k/ncells*100))
    end
    
end

fclose(fID);

return
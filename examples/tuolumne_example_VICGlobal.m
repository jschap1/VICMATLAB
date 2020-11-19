% Using VIC to model the Tuolumne River Basin (with VICGlobal v1.6)
%
% Written 8/4/2020 JRS
% Write this up in Markdown and post to vicmatlab/docs
%
% This is a tutorial to demonstrate how to use VICMATLAB to prepare inputs
% for and analyze results from the VIC model. The test basin is the Upper
% Tuolumne River Basin, near Yosemite National Park, California.
%
% Before running this code, make sure that all files and paths are set
% properly according to your system. You will need VIC-4 and VIC-5, as
% well as the VICMATLAB toolbox.
%
% This script is set up for the VIC parameters from Livneh et al. (2013). 
% It could easily be adjusted for another set of VIC parameters, such as
% VICGlobal (actually do this).
% 
% You should have:
% -shapefile of the study area
% -classic mode meteorological forcing data
% -basin mask
% -VICGlobal parameters (classic and image mode)

addpath(genpath('~/Documents/Codes/VICMATLAB/vicmatlab'))

%% Subset the VICGlobal soil parameter file

soilfile = '/hdd/Data/VICParametersGlobal/VICGlobal/v1.6/Classic/soils_3L_MERIT.txt';
disp('Loading soil parameter file')
soils = load(soilfile);
disp('Soil parameter file has been loaded')

% extent of the study area; can be specified multiple ways
% extent = '/Users/jschap/Documents/Codes/VICMATLAB/data/revised_basin/upptuo_mask.tif'; 
extent = '/home/jschap/Documents/Codes/VICMATLAB/data/upptuo_mask.tif';
grid_decimal = 5; % number of decimals used in forcing filenames

outformat = '3l'; % format of input soil parameter file (number of soil layers)
outname = '/home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/tuo_soils_VG.txt';

generate_tif = 1; % generates tifs in same directory as subsetted soil file
setup = '3L-no-org-frost-msds'; % options are 2L, 2L-no-org-fs-july_tavg, 3L, 3L-no-org-frost-msds, and livneh
[soils_tuo, soilvarpath] = subset_soils(soils, extent, outname, outformat, grid_decimal, generate_tif, setup);

%% Run classic simulation

% Do on command line

%% Make forcing files for the image mode simulation

forcdir = '/home/jschap/Documents/Codes/VICMATLAB/data/Forcings/disagg_forc_2009-2011';
outname = '/home/jschap/Documents/Codes/VICMATLAB/data/Forcings/netcdf_forcings/tuo_forc';
start_date = datetime(2009, 1, 1, 0, 0, 0); % need to specify hours
end_date = datetime(2011, 12, 31, 23, 0, 0);
nt_per_day = 24;
prefix = 'full_data_';
precision = 5;
convert_forcing(forcdir, prefix, outname, precision, start_date, end_date, nt_per_day)

%% Subset parameters for image mode simulation

compname = 'atlantic';
switch compname
    case 'atlantic'
        param_name = '/home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/tuolumne_domain.nc';
        domain_name = '/home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/tuolumne_params.nc';
        basinmaskname = '/home/jschap/Documents/Codes/VICMATLAB/data/upptuo_mask.tif';
        global_domain = '/hdd/Data/VICParametersGlobal/VICGlobal/v1.6/output_latest_aug14/output_latest/VICGlobal_domain.nc';
        global_params = '/hdd/Data/VICParametersGlobal/VICGlobal/v1.6/output_latest_aug14/output_latest/VICGlobal_params.nc';
    case 'pacific'
        param_name = '/home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/tuolumne_domain.nc';
        domain_name = '/home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/tuolumne_params.nc';
        basinmaskname = '/home/jschap/Documents/Codes/VICMATLAB/data/upptuo_mask.tif';
        global_domain = '/hdd/Data/VICParametersGlobal/VICGlobal/v1.6/output_latest_aug14/output_latest/VICGlobal_domain.nc';
        global_params = '/hdd/Data/VICParametersGlobal/VICGlobal/v1.6/output_latest_aug14/output_latest/VICGlobal_params.nc';        
end
subset_domain(basinmaskname, global_domain, domain_name)
subset_parameter(basinmaskname, global_params, param_name)

%% Run image simulation

% Do on command line

%% Process classic mode outputs

% Borrowing code from /ESSD/codes/wrapper.m to load the classic mode
% outputs

timestep_out = 'daily';
prefix = 'fluxes_'; % the underscore is needed
classic.outdir = '/home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/out_2009-2011';
classic.metadata = get_vic_run_metadata2(classic.outdir, timestep_out, prefix);
classic.processed_dir = '/home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/processed_classic';
basin_mask = geotiffread2('/home/jschap/Documents/Codes/VICMATLAB/data/upptuo_mask.tif');

% starting index at 4 because the first three are year, month, and day
for i=4:length(classic.metadata.vars)
    varname = classic.metadata.vars{i};
    savename = fullfile(classic.processed_dir, [varname '.mat']);
    try
        [classic.timevector, ~, classic.info] = load_vic_output2(classic.outdir, varname, prefix, 'daily', 1, savename);
        masksavename = fullfile(classic.processed_dir, [varname '_masked.mat']);
        mask_results_to_bb(savename, basin_mask, masksavename);
    catch
        disp(savename)
    end
end

%% Compare image and classic mode outputs

% Load classic outputs
for i=4:length(classic.metadata.vars)
    varname = classic.metadata.vars{i};
    tempname = strsplit(varname, 'OUT_');
    masksavename = fullfile(classic.processed_dir, [varname '_masked.mat']);
    temp = load(masksavename);
    classic.out.(tempname{2}) = temp.output_map_masked;
end

% Load image outputs
vic_out_name_image = '/home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/out_image_2009-2011/fluxes.2009-01-01.nc';
image = load_vic_output_image(vic_out_name_image);

%% SWE plots

% Plot SWE time-average maps
classic.avgmap.swe = nanmean(classic.out.SWE,3);
image.avgmap.swe = nanmean(image.OUT_SWE,3);

figure, subplot(2,1,1)
plotraster(classic.metadata.lon, classic.metadata.lat, classic.avgmap.swe, 'Classic Mode SWE (mm)')
subplot(2,1,2)
plotraster(image.lon, image.lat, image.avgmap.swe, 'Image Mode SWE (mm)')

% Plot SWE basin-average time series
classic.avgts.swe = calc_basin_avg_ts(classic.out.SWE);
image.avgts.swe = calc_basin_avg_ts(image.OUT_SWE);

figure
plot(classic.metadata.time, classic.avgts.swe)
hold on
plot(image.time, image.avgts.swe)
xlabel('Time')
ylabel('SWE (mm)')
legend('Classic Mode','Image Mode')

%% Baseflow plots

% Baseflow
classic.avgmap.qb = nanmean(classic.out.BASEFLOW,3);
image.avgmap.qb = nanmean(image.OUT_BASEFLOW,3);
classic.avgts.qb = calc_basin_avg_ts(classic.out.BASEFLOW);
image.avgts.qb = calc_basin_avg_ts(image.OUT_BASEFLOW);

figure, subplot(2,1,1)
plotraster(classic.metadata.lon, classic.metadata.lat, classic.avgmap.qb, 'Classic Mode Baseflow (mm)')
subplot(2,1,2)
plotraster(image.lon, image.lat, image.avgmap.qb, 'Image Mode Baseflow (mm)')

figure
plot(classic.metadata.time, classic.avgts.qb)
hold on
plot(image.time, image.avgts.qb)
xlabel('Time')
ylabel('Baseflow (mm)')
legend('Classic Mode','Image Mode')

%% Runoff plots

% Runoff
classic.avgmap.qd = nanmean(classic.out.RUNOFF,3);
image.avgmap.qd = nanmean(image.OUT_RUNOFF,3);
classic.avgts.qd = calc_basin_avg_ts(classic.out.RUNOFF);
image.avgts.qd = calc_basin_avg_ts(image.OUT_RUNOFF);

figure, subplot(2,1,1)
plotraster(classic.metadata.lon, classic.metadata.lat, classic.avgmap.qd, 'Classic Mode Runoff (mm)')
subplot(2,1,2)
plotraster(image.lon, image.lat, image.avgmap.qd, 'Image Mode Runoff (mm)')

figure
plot(classic.metadata.time, classic.avgts.qd)
hold on
plot(image.time, image.avgts.qd)
xlabel('Time')
ylabel('Runoff (mm)')
legend('Classic Mode','Image Mode')

%% Look at an individual grid cell

ts_classic = squeeze(classic.out.BASEFLOW(1,1,:));
ts_image = squeeze(image.OUT_BASEFLOW(1,1,:));
figure
plot(classic.metadata.time, ts_classic, '--');
hold on
plot(image.time, ts_image)
legend('Classic Mode','Image Mode')

%% Compare soil parameters

dsmax_classic = geotiffread2('/home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/soiltifs_classic/dsmax.tif');

image_params = '/home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/tuolumne_params.nc';
image_params_info = ncinfo(image_params);
dsmax_image = ncread(image_params, 'Dsmax');
dsmax_image = dsmax_image';

figure, subplot(2,1,1)
plotraster(classic.metadata.lon, classic.metadata.lat, dsmax_classic, 'Classic Mode Dsmax (mm/day)')
subplot(2,1,2)
plotraster(image.lon, image.lat, basin_mask.*dsmax_image, 'Image Mode Dsmax (mm/day)')

% They check out. Exactly the same.


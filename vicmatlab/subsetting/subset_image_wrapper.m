%% Image mode

% Subset Image
%
% Wrapper for subsetting the global VIC parameter file
% to a user-specified domain
%
% There are different codes to use for different parameter sets
% For VICGlobal, use subset_domain and subset_parameter
% For L2013, use classic2image
% For BV2019, use subset_domain and subset_parameter

addpath(genpath('/home/jschap/Documents/Codes/VICMATLAB/vicmatlab'))
outdir = '/home/jschap/Documents/Codes/VICMATLAB/vicmatlab/subsetting/subset_image/';
basinname = 'ucrb';

%% Specify extent

% Input may be a shapefile or a list of coordinates
% The coordinate system must have the same resolution, origin, and CRS as
% the image mode parameter files that you are subsetting.

ucrb_shp = '/home/jschap/Documents/ESSD/data/bb.shp';
ucrb_tif = '/home/jschap/Documents/ESSD/data/colo_mask.tif';
tuolumne_tif = '/home/jschap/Documents/Codes/VICMATLAB/vicmatlab/subsetting/subset_image/upper_tuolumne_basin_wgs.tif';

% If you want to use a basin mask (geotiff), you must get a list of coords
% basin_coords = basin_mask2coordinate_list(ucrb_tif);

%% Subset BV2019

% BV2019
bv2019domain = '/media/jschap/HD_ExFAT/vicglobal-image-mode-inputs/VICGlobal_domain.nc';
bv2019params = '/home/jschap/Documents/Data/BV2019/params.CONUS_MX.MOD_IGBP.mode.2000_2016.nc';
subset_domain(ucrb_tif, bv2019domain, fullfile(outdir, [basinname, '_bv2019_domain.nc']));
subset_parameter(ucrb_tif, bv2019params, fullfile(outdir, [basinname, '_bv2019_params.nc']));

%% Subset L2013

% L2013
% Convert Classic mode inputs to image mode and do subsetting in one step
% The subsetting is controlled by the forcings, so they need to be subset
% to the domain first. 
%
% NEEDS WORK

parpath = '/home/jschap/Documents/Data/VICParametersCONUS/';
inputs.veglib = fullfile(parpath, 'vic_veglib_nohead.txt');
inputs.soilparfile = fullfile(parpath, 'vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt');
inputs.snowband = fullfile(parpath, 'vic.snow.0625.new.cal.adj.can.5bands');
inputs.vegparam = fullfile(parpath, 'vic.veg.0625.new.cal.adj.can');

% Forcings
inputs.forcdir = '/media/jschap/HD_ExFAT/ucrb/disagg_forcings_L13_clip/full_data*';
inputs.domainfile_name = fullfile(outdir, [basin_name, '_domain.nc']);
inputs.params_name = fullfile(outdir, [basin_name, '_params.nc']);
subset_forcings
subset_forcings_to_spf
classic2image()

%% Subset VICGlobal

% VICGlobal
vicglobaldomain = '/media/jschap/HD_ExFAT/vicglobal-image-mode-inputs/VICGlobal_domain.nc';
vicglobalparams = '/media/jschap/HD_ExFAT/vicglobal-image-mode-inputs/VICGlobal_params.nc';
subset_domain(basin_coords, l2013domain, fullfile(outdir,' upper_tuolumne_domain.nc'))

%% Aggregate forcings (if desired)

forcdir = '/media/jschap/HD_ExFAT/ucrb/disagg_forcings_L13_clip';
delt_in = 1; % hours
delt_out = 6; % hours
starttime = datetime(1980, 1, 1, 0, 0, 0);
outdir = '/media/jschap/HD_ExFAT/ucrb/L13_forc_6hr_ucrb';
precip_col = 1; % need to specify because precip is summed, not averaged,
% to keep units consistant (mm/timestep)
hpc = 1;
[aggf, aggt] = aggregate_forcing(forcdir, delt_in, delt_out, starttime, outdir, precip_col, hpc);
% ^^^ takes about 2 hours to run, makes about 20 GB of output data

%% Convert forcings

forcdir = '/media/jschap/HD_ExFAT/ucrb/L13_forc_6hr_ucrb';
outname_domain = fullfile('/media/jschap/HD_ExFAT/ucrb/L13_forc_6hr_ucrb_nc', basinname);
start_date = datetime(1980, 1, 1, 0, 0, 0); % need to specify hours
end_date = date
time(2011, 12, 31, 23, 0, 0);
nt_per_day = 4;
prefix = 'full_data_';
precision = 5;
convert_forcing(forcdir, prefix, outname_domain, precision, start_date, end_date, nt_per_day)
% ^^ takes 23 minutes per year for the UCRB (7833 grid cells)

%%


% Define extent
% extent = horzcat([75; 76], [34; 35]);
extent = '/Users/jschap/Documents/Research/SWOTDA_temp/Tuolumne/Tuolumne2/Shapefiles/upper_tuolumne.shp';
% extent = horzcat([-85.5625; -97.6875], [36.625; 48.0625]);
% extent = horzcat([65; 85], [23; 38]);

% Subset domain
global_domain = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Image/VICGlobal_domain.nc';
outname_domain = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_7/Image_VG/tuo_domain.nc';
subset_domain(extent, global_domain, outname_domain)

% Subset parameters
global_params = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Image/VICGlobal_params.nc';
% outname_params = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_4/Image_VICGlobal/params_subL.nc';
outname_params = '/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_7/Image_VG/tuo_params.nc';
subset_parameter(extent, global_params, outname_params)

% takes a long time to run (took a couple hours for the Tuolumne basin test case)
% also uses a lot of RAM, like 20 GB
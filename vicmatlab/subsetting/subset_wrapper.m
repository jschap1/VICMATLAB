% Wrapper
%
% Wrapper for subsetting the VICGlobal input files to a particular area of interest
% The VICGlobal subsetting functions require the image processing toolbox.

% Not necessary, but we will run this from within its own directory
addpath('/home/jschap/Documents/Codes/VICMATLAB/vicmatlab/subsetting/')
% cd('/Users/jschap/Documents/Codes/VICMATLAB/Subsetting/Subsetting')
clearvars -except soils_vg soils

%% Classic mode: Subsets the soil parameter file.

% Load soil parameter file
soilfile = '/Volumes/HD3/VICParametersCONUS/vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt';
% soilfile = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_5/Classic/soils_3L_MERIT.txt';
disp('Loading soil parameter file')
% soils_vg = load(soilfile);
soils = load(soilfile);
disp('Soil parameter file has been loaded')

% Define extent
% extent = [75, 76, 34, 35];
extent = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif';
% extent = '/Users/jschap/Documents/Research/Glaciers/Skagit/skagit_mask.tif';
% extent = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/dem.tif';
% use full name; do not use . in file name
% extent = './Data/IRB/VIC/34N_75E/dem.tif';
% extent = '/Volumes/HD4/SWOTDA/Data/UMRB/dem.tif'; % UMRB

grid_decimal = 5; % number of decimals used in forcing filenames
% outformat = '3l';
outformat = 'livneh'; % format of input soil parameter file (number of soil layers)
outname = '/Volumes/HD4/SWOTDA/Data/Colorado/colo_soils_L15.txt';
% outname = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/Classic/soils_.txt';
% outname = '/Volumes/HD4/SWOTDA/Data/UMRB/soils_umrb_vg.txt';
% outname = '/Users/jschap/Documents/Research/Glaciers/Skagit/skagit_soils.txt';

soils_subset = subset_soils(soils, extent, outname, outformat, grid_decimal, generate_tif, setup);

%% Image mode

% Subset Image
%
% Wrapper for subsetting the global VIC parameter file
% to a user-specified domain

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

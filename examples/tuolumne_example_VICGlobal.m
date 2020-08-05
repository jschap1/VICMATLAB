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

classic.out.swe = load(fullfile(classic.processed_dir, 'OUT_SWE.mat'))































%% Subset meteorological forcing data 

% Create ASCII input files for the meteorological forcing data. Subset it to the basin.
% force_in = '/Volumes/HD3/Livneh_2013/MetNC'; % directory where forcing data are located
force_in = '/media/jschap/HD_ExFAT/Livneh_2013/MetNC/';
numforcings = 4; % number of forcings in daily CONUS daily forcing file

% Beginning and ending years of simulation
beginyear = 2009;
endyear = 2011;

% Directory where clipped forcing files should be saved
force_out = ['./data/forc_' num2str(beginyear) '-' num2str(endyear)];

grid_decimal = 5; % number of decimals used in forcing filenames
maskname = './data/upptuo_mask.tif'; % basin mask

temp = subset_forcings(force_in, force_out, beginyear, endyear, maskname);

% temp = subset_forcings(force_in, force_out, beginyear, endyear, grid_decimal, numforcings, maskname);

%% Plot the forcing data

forcingpath = './data/forc_2009-2011';
precision = 5; 
varnames = {'PRECIP','TMIN','TMAX','WIND'};
prefix = 'data_';
forc = load_ascii_forcings(forcingpath, prefix, precision, varnames);

figure, subplot(2,1,1)
tmax_map = xyz2grid(forc.lon, forc.lat, mean(forc.TMAX)');
plotraster(forc.lon, forc.lat, tmax_map, 'TMAX (deg. C)')

subplot(2,1,2)
prec_map = xyz2grid(forc.lon, forc.lat, mean(forc.PRECIP)');
plotraster(forc.lon, forc.lat, prec_map, 'PREC (mm/day)')

% Save as Geotiff
geotiffwrite('./data/livneh_precipitation_2009-2011_average.tif', ...
    flipud(prec_map), R)

%% Disaggregate meteorological forcing data

% Run the VIC model as a meteorological forcing disaggregator.

% Create a global parameter file and run the following code to disaggregate 
% the met. forcing data with MT-CLIM

disagg_force_out = ['./data/disagg_forc_' num2str(beginyear) '-' num2str(endyear)];
mkdir(disagg_force_out)
disp(['Created directory ' disagg_force_out ' for disaggregated forcings'])

disp('Running met. forcing disaggregation')
tic
system('/home/jschap/Documents/Software/VIC/src/vicNl -g ./data/global_param_disagg.txt')
toc
% Takes 54 sec on my laptop for 175 grid cells, 3 years
% Takes 27 seconds on my desktop! Woohoo for computing power!

%% Plot the disaggregated forcings

forcingpath = './data/disagg_forc_2009-2011/';
precision = 5; 
varnames = {'PRECIP','AIR_TEMP','SHORTWAVE','LONGWAVE','DENSITY','PRESSURE','VP','WIND'};
prefix = 'full_data_';
forc = load_ascii_forcings(forcingpath, prefix, precision, varnames);

nvars =length(varnames);

avg_maps = struct();
avg_maps.names = varnames;
for i=1:nvars
    avg_maps.(varnames{i}) = xyz2grid(forc.lon, forc.lat, mean(forc.(varnames{i}),1)');
%     avg_maps.(varnames{i}) = fliplr(xyz2grid(forc.lon, forc.lat, mean(forc.(varnames{i}),1)'));
end

labels = {'Precip (mm/hr)','T (deg. C)','SW (W/m^2)','LW (W/m^2)','\rho (kg/m^3)','P (kPa)','VP (kPa)','u (m/s)'};

figure
for i=1:nvars
    subplot(4,2,i)
    plotraster(forc.lon, forc.lat, avg_maps.(varnames{i}), labels{i})
%     plotraster(forc.lon, forc.lat, avg_maps.(varnames{i}), varnames{i})
end

% Save as Geotiff
geotiffwrite('./data/livneh_precipitation_downscaled_2009-2011_average.tif', ...
    flipud(avg_maps.PRECIP), R)

%% Run VIC for water years 2010-2011. This takes about 15 minutes on my computer.

% Set up the global parameter file for this model run, first. 

% This setup is for the Upper Tuolumne Basin, in energy balance mode,
% with the frozen soils module disabled to reduce computational requirements. 
% Use the template global parameter file to generate the specific list of outputs 
% that we are using in this example.

% system('/Volumes/HD3/SWOTDA/Software/VIC-VIC.5.1.0.rc1/vic/drivers/classic/vic_classic.exe -g ./data/global_param.txt')

system('/home/jschap/Documents/Software/VIC/vic/drivers/classic/vic_classic.exe -g ./data/global_param.txt')

% Analyze the outputs from the VIC simulation. 
% First, re-organize the VIC outputs by entering the following commands on
% the command line
% cd /home/jschap/Documents/Codes/VICMATLAB/data/out_2009-2011/
% mkdir eb; mv eb*.txt eb
% mkdir wb; mv wb*.txt wb

%% Process VIC output data with VICMATLAB

% Get VIC run metadata
vic_out_dir = './data/out_2009-2011/';
timestep_out = 'daily';
info = get_vic_run_metadata(vic_out_dir, timestep_out);

% Create directories to store processed outputs
results_dir = fullfile(vic_out_dir, 'processed');
figdir = fullfile(vic_out_dir, 'figures');
mkdir(results_dir)
disp(['Created directory for results: ' results_dir]);
mkdir(figdir)
disp(['Created directory for figures: ' figdir]);

% Save the metadata from the VIC run
save(fullfile(results_dir, 'vic_run_metadata.mat'), 'info');
disp(['Saved VIC run metadata as ' fullfile(results_dir, 'vic_run_metadata.mat')])

% Read in and plot the VIC results
[~, swe, ~] = load_vic_output(vic_out_dir, 'swe');

swe_map = xyz2grid(info.lon, info.lat, mean(swe,2));
figure, plotraster(info.lon, info.lat, swe_map, 'SWE (mm)')

% mkdir wb; mv wb*.txt wb
[~, precip, ~] = load_vic_output(vic_out_dir, 'precipitation');
% only works with water balance variables right now.
% need to generalize the code to handle energy balance variables
% also should find a less clunky way to specify the variable name, aside
% from the column number

precip_map = xyz2grid(info.lon, info.lat, mean(precip,2));
figure, plotraster(info.lon, info.lat, precip_map, 'Precip (mm/day)')

% Save as Geotiff
geotiffwrite('./data/precipitation_output_2009-2011_average.tif', ...
    flipud(precip_map), R)

%% Convert inputs from ASCII to NetCDF

% addpath(genpath('/home/jschap/Documents/Research/Codes/VICMATLAB/vicmatlab_dev'))

% must use entire CONUS domain file for the current setup

wkpath = '/home/jschap/Documents/Codes/VICMATLAB/';
parpath = '/home/jschap/Documents/Data/VICParametersCONUS/';

inputs.veglib = fullfile(parpath, 'vic_veglib_nohead.txt');
inputs.soilparfile = fullfile(parpath, 'vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt');
inputs.snowband = fullfile(parpath, 'vic.snow.0625.new.cal.adj.can.5bands');
inputs.vegparam = fullfile(parpath, 'vic.veg.0625.new.cal.adj.can');

inputs.forcdir = fullfile(wkpath, '/data/disagg_forc_2009-2011/full_data*');
inputs.domainfile_name = fullfile(wkpath, '/data/tuolumne_domain2.nc');
inputs.params_name = fullfile(wkpath, '/data/tuolumne_params2.nc');

classic2image(inputs);

%% Convert forcings from ASCII to NetCDF

cd ~/Documents/Codes/VICMATLAB/
addpath(genpath('./vicmatlab'))

% forcingpath = './data/disagg_forc_2009-2011/';
forcingpath = '/home/jschap/Documents/Backups/data/disagg_forc_2009-2011';
precision = 5; 
start_date = datetime(2009,1,1,0,0,0);
end_date = datetime(2011,12,31,23,0,0);
nt_per_day = 24;
prefix = 'full_data_';
outname = './data/netcdf_forcings/forc_tuo';

convert_forcing(forcingpath, prefix, outname, precision, start_date, end_date, nt_per_day)


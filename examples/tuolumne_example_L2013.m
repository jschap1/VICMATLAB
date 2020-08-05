% Using VIC to model the Tuolumne River Basin
% A demo for new users of the VICMATLAB Toolbox
% Written 3/6/2020 JRS
% Last updated 4/13/2020 JRS

cd ~/Documents/Codes/VICMATLAB/

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

%% Obtain a shapefile of the study area. 

% You can download this from the USGS National Map Viewer website here: 
% https://viewer.nationalmap.gov/basic/

% The first step is to extract the Tuolumne River Basin in R. 

% Here is some sample code.
% huc8 <- readOGR(".../Shape/WBDHU8.shp")
% upptuo <- huc8[huc8$Name == "Upper Tuolumne",]
% writeOGR(tuolumne, "./data/upptuo.shp", driver = "ESRI Shapefile", layer = "bdy)
% Load the RoutR package, and use it to make an elevation map masked to the basin (also in R).
% merit_mask <- raster(".../MERIT/DEM/Merged_1_16/merit_mask_1_16.tif")
% merit_dem <- raster(".../MERIT/DEM/Merged_1_16/merged_merit_dem_1_16.tif")
% upptuo_dem <- clip_raster(merit_dem, upptuo) # or see alternative method in R script.
% writeRaster(upptuo_dem, "./data/revised_basin/upptuo_dem.tif", NAflag = 0, datatype = "INT2S")

%% Download forcing data

% Download Livneh et al. (2013) or Livneh et al. (2015) NetCDF meteorological forcing data.
% We will use Livneh et al. (2013) data for this tutorial. 
% The data can be downloaded from https://www.esrl.noaa.gov/psd/data/gridded/data.livneh.html.

%% Make a basin mask from the DEM.

addpath(genpath('./vicmatlab'));
% addpath(genpath('./vicmatlab1.0'));

[upptuo_dem, R, lon, lat] = geotiffread2('./data/upptuo_dem.tif');
basin_mask = double(upptuo_dem);
basin_mask(basin_mask == 0) = NaN;
basin_mask(~isnan(basin_mask)) = 1;

figure, plotraster(lon, lat, basin_mask,'Basin Mask');

geotiffwrite('./data/upptuo_mask.tif', flipud(basin_mask), R);

%% Create soil parameter file

% Create the soil parameter file for the basin by subsetting the 
% soil parameter file for the CONUS, from Livneh et al. (2013). 
% This can be downloaded from the same website as the meteorological forcing data.

soilfile = '/home/jschap/Documents/Data/VICParametersCONUS/vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt';
% soilfile = '/Volumes/HD3/VICParametersCONUS/vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt';
disp('Loading soil parameter file')
soils = load(soilfile);
disp('Soil parameter file has been loaded')

% extent = '/Users/jschap/Documents/Codes/VICMATLAB/data/revised_basin/upptuo_mask.tif'; % extent of the study area; can be specified multiple ways
extent = fullfile(pwd, './data/upptuo_mask.tif');
grid_decimal = 5; % number of decimals used in forcing filenames

outformat = 'livneh'; % format of input soil parameter file (number of soil layers)
outname = fullfile(pwd, '/data/upptuo_soils_livneh.txt');

generate_tif = 1;
setup = 'livneh';
[soils_tuo, soilvarpath] = subset_soils(soils, extent, outname, outformat, grid_decimal, generate_tif, setup);

%% Plot the soil parameters

figure, subplot(2,1,1)
[elev, R, lon, lat] = geotiffread2(fullfile(soilvarpath, 'elev.tif'));
plotraster(lon, lat, elev, 'Elevation (m)')

subplot(2,1,2)
[dsmax, R, lon, lat] = geotiffread2(fullfile(soilvarpath, 'dsmax.tif'));
plotraster(lon, lat, dsmax, 'Dsmax (mm/day)')

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


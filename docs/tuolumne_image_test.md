# Modeling the Tuolumne Basin in Image Mode

This tutorial demonstrates the use of the VICMATLAB subsetting functions for subsetting the VICGlobal dataset to a particular domain. Also, we test the image mode parameters and make sure they produce the same model results as the classic mode parameters.

```matlab
addpath(genpath('~/Documents/Codes/VICMATLAB/vicmatlab'))
outdir = '/hdd/ESSD/image_mode/';
basinname = 'tuolumne';

% Subset domain file to Tuolumne Basin
basinmask = '/home/jschap/Documents/Codes/VICMATLAB/data/tifs/upptuo_mask.tif';

fulldomain = '/hdd/Data/VICParametersGlobal/VICGlobal/v1.6/Image/VICGlobal_domain.nc';
fullparams = '/hdd/Data/VICParametersGlobal/VICGlobal/v1.6/Image/VICGlobal_params.nc';
subset_domain(basinmask, fulldomain, fullfile(outdir, [basinname, '_domain.nc']));

% Need 20 GB of RAM for this to run
subset_parameter2(basinmask, fullparams, fullfile(outdir, [basinname, '_params.nc']));
```

The last line requires 20 GB of RAM to run, so I exported my input files and code to Atlantic, where I have 48 GB of RAM, and I ran the subsetting there. (It would be convenient to split the image mode parameters into several files so that people without access to this much RAM can still work with the files.)

## Debugging
For some reason, the image mode Tuolumne simulation does not run. There is an error 

		[ERROR] ../shared_image/src/get_nc_field.c:49: errno: NetCDF: Index exceeds dimension bound: Error getting values for tas in /hdd/ESSD/tuo_image_example/v1.5/forc_tuo_2009.nc 

I changed the domain and parameter subsetting codes to use the correct basin mask for the mask and run_cell variables, but this did not help. Perhaps I need to change more/all of the parameters to match the mask?

### Comparing format of forcing files

The Stehekin simulation runs, but the VICGlobal simulation does not run; the difference is due to the different input files. 

Stehekin:
Dimensions:
lat = 8
lon = 10
time = 240
float: pres, wind, tas, dlwrf, prcp, dswrf, vp
double: lon, lat
int: time (hours since 1949/1/1 00:00:00)
10 days, hourly forcing data
Uses nodata values for places outside the basin (though this probably isn't important, given that run_cell controls which cells are used). 

VICGlobal (Tuo):
Dimensions:
lat = 11
lon = 32
time = 8760
float: pres, wind, tas, dlwrf, prcp, dswrf, vp
double: lon, lat
int: mask, time (hours since 2009/1/1)
10 days, hourly forcing data
Uses nodata values for places outside the basin (though this probably isn't important, given that run_cell controls which cells are used). 

Potential solution: just make forcing files that cover the whole basin extent. Modify the forcing creation method to do this, as an option. Unless the time is the issue, in which case I might need to add hours, etc. to the description for the time variable. Does the image driver use this info?

### Comparing format of domain files

		scp jschap@american.seas.ucla.edu:~/Documents/Research/VICGlobal/Data/tuolumne_image_subsetting/stehekin_domain.nc .

Stehekin: 
Dimensions:
lat = 8
lon = 10
month = 12
Variables: 
float: t_pk, dur, elev, frac, area
int: month, mask
double: lon, lat

t_pk and dur are special options for precipitation.

month = 0:11
mask = ones(8,10)
elev = 8 x 10 matrix of elevs
frac, area = similar to elev

VICGlobal:
Dimensions
lat = 8
lon = 10
Variables:
float: frac, area
int: mask
double: lat, lon

The VICGlobal and Stehekin domain files are basically identical.

VICGlobal (full, from Dongyue):

Dimensions
lat = 2267
lon = 5760
Variables:
double: frac, area, lat, lon
int: mask (either 1 or nodata value)

### Comparing format of parameter files

Stehekin:

Dimensions:
veg_class = 12
string20 = 20
month = 12
lat = 8
lon = 10
root_zone = 3
snow_band = 5
nlayer = 3

Variables:
char: veg_descr, 
double: veg_rough, displacement, LAI, albedo, wind_atten, wind_h, RGL, rmin, trunk_ratio, rarc, Cv, rad_atten, overstory, root_depth, root_fract, lon, snow_rough, elev, off_gmt, avg_T, lons, annual_prec, Ds, lats, infilt, Ws, c, Dsmax, dp, rough, cellnum, lat, AreaFract, elevation, Pfactor, quartz, depth, Ksat, bulk_density, bubble, resid_moist, Wcr_FRACT, soil_density, phi_s, init_moist, expt, Wpwp_FRACT, 
int: month, veg_class, fs_active, run_cell, mask, gridcell, Nveg, snow_band, root_zone, layer, 

Months go from 1:12
Mask has all cells in it, so does gridcell
Run_cell is only the domain you wish to run (the basin) (1, 0)

VICGlobal (subset)

Dimensions: 
lat = 8
lon = 10
nlayer = 3
snow_band = 5
veg_class = 17
root_zone = 3
month = 12

Variables:
double: lat, lon, mask, lons, lats, infilt, Ds, Dsmax, Ws, c, elev, avg_T, dp, off_gmt, rough, snow_rough, annual_prec, July_Tavg, cellnum, expt, Ksat, phi_s, init_moist, depth, bubble, quartz, bulk_density, soil_density, Wcr_FRACT, Wpwp_FRACT, resid_moist, AreaFract, elevation, Pfactor, Cv, rarc, rmin, wind_h, RGL, rad_atten, wind_atten, trunk_ratio, root_depth, root_fract, LAI, albedo, veg_rough, displacement, fcanopy
int: run_cell, gridcell, fs_active, Nveg, overstory

Months go from 1:12
Mask has all cells in it
Run_cell has all cells in it

VICGlobal (full, from Dongyue)

Dimensions: 
lat = 2267
lon = 5760
nlayer = 3
snow_band = 5
veg_class = 17
root_zone = 3
month = 12

Variables:
char: veg_descr
double: lat, lon, lats, lons, infilt, Ds, Dsmax, Ws, c, expt, Ksat, phi_s, init_moist, elev, depth, avg_T, dp, bubble, quartz, bulk_density, soil_density, off_gmt, Wcr_FRACT, Wpwp_FRACT, rough, snow_rough, annual_prec, resid_moist, July_Tavg, cellnum, AreaFract, elevatoin, Pfactor, Cv, root_depth, root_fract, LAI, rarc, rmin, wind_h, RGL, rad_atten, wind_atten, trunk_ratio, albedo, veg_rough, displacement, fcanopy
int: mask, layer, run_cell, gridcell, fs_active, snow_band, veg_class, root_zone, month, Nveg, overstory, 

Months go from 1:12
Mask has all cells in it
Run_cell has all cells in it

## Splitting VICGlobal parameters into pieces

Subsetting the whole parameter file (140 GB) is too much for most computers. Therefore, I have split the parameter file into smaller, more manageable pieces. 

## Running the Stehekin example

./vic_image.exe -g /home/jschap/Documents/Software/VIC_sample_data/image/Stehekin/parameters/global_param.Stehekin.L2015.txt

The Stehekin simulation worked fine, although I was not able to write to the /VIC_sample_data/ directory for some reason. Instead, I wrote the result to /home/jschap/Documents/.

Inspecting the Stehekin parameters, the masks cover the full domain, but the run_cell parameter matches the forcing files.

Try subsetting VICGlobal to the Stehekin domain and see if it runs. If so, then the issue is with my forcing files.

### Procedure

Create basin mask for Stehekin basin.

```r
library(raster)
library(ncdf4)
r <- raster("file.tif")
```

Send mask to remote computer

```bash
LWKDIR=/home/jschap/Documents/Software/VIC_sample_data/image/Stehekin/parameters
RWKDIR=/Users/jschap/Documents/Research/VICGlobal/Data/tuolumne_image_subsetting
scp ${LWKDIR}/stehekin_mask.tif jschap@american.seas.ucla.edu:${RWKDIR}/ 
```

Push the latest version of VICMATLAB

```bash
cd /home/jschap/Documents/Codes/VICMATLAB
git add .
git commit -m "subsetting update"
git push origin master
```

Connect to remote computer
```bash
ssh jschap@american.seas.ucla.edu
RWKDIR=/Users/jschap/Documents/Research/VICGlobal/Data/tuolumne_image_subsetting
cd $RWKDIR
matlab -nodisplay -nosplash
```

Subset VICGlobal domain and parameter files to Stehekin basin mask
```matlab
addpath(genpath('~/Documents/Codes/VICMATLAB/vicmatlab'))
outdir = '/Users/jschap/Documents/Research/VICGlobal/Data/tuolumne_image_subsetting';
basinmask = '/Users/jschap/Documents/Research/VICGlobal/Data/tuolumne_image_subsetting/stehekin_mask.tif';
fulldomain = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/Image/VICGlobal_domain.nc';
fullparams = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/Image/VICGlobal_params.nc';
basinname = 'stehekin';
subset_domain(basinmask, fulldomain, fullfile(outdir, [basinname, '_domain.nc']));
subset_parameter2(basinmask, fullparams, fullfile(outdir, [basinname, '_params.nc']));
```

Copy the subsetted domain and parameter files back to local
```bash
scp jschap@american.seas.ucla.edu:${RWKDIR}/stehekin_* ${LWKDIR}/
```

Set up global parameter file

Run VIC
```bash 
# Run VIC
/Documents/Software/VIC/vic/drivers/image/vic_image.exe -g ${LWKDIR}/global_param.txt
```
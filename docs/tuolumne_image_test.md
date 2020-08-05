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

Subsetting the whole parameter file (140 GB) is too much for most computers. Therefore, I will split the parameter file into smaller, more manageable pieces. This is one of those things that I will do "in the future."

## Running the Stehekin example

./vic_image.exe -g /home/jschap/Documents/Software/VIC_sample_data/image/Stehekin/parameters/global_param.Stehekin.L2015.txt

The Stehekin simulation worked fine, although I was not able to write to the /VIC_sample_data/ directory for some reason. Instead, I wrote the result to /home/jschap/Documents/.

Inspecting the Stehekin parameters, the masks cover the full domain, but the run_cell parameter matches the forcing files.

Try subsetting VICGlobal to the Stehekin domain and see if it runs. If so, then the issue is with my forcing files.

### Procedure

Create basin mask for Stehekin basin in R.

```r
library(raster)
library(ncdf4)
LWKDIR <- "/home/jschap/Documents/Software/VIC_sample_data/image/Stehekin/parameters"
stehekin <- raster(file.path(LWKDIR, "params.Stehekin.L2015.nc"))
stehe <- nc_open(file.path(LWKDIR, "params.Stehekin.L2015.nc")
mask1 <- ncvar_get(stehe, "run_cell")
lat <- ncvar_get(stehe, "lat")
lon <- ncvar_get(stehe, "lon")
values(stehekin) <- as.vector(mask1)
writeRaster(stehekin, file.path(LWKDIR, "stehekin_mask.tif"))
```

Switch to command line. Send mask to remote computer

```bash
LWKDIR=/home/jschap/Documents/Software/VIC_sample_data/image/Stehekin/parameters
RWKDIR=/Users/jschap/Documents/Research/VICGlobal/Data/tuolumne_image_subsetting
scp ${LWKDIR}/stehekin_mask.tif jschap@american.seas.ucla.edu:${RWKDIR}/ 
```

Push the latest version of VICMATLAB.

```bash
cd /home/jschap/Documents/Codes/VICMATLAB
git add .
git commit -m "subsetting update"
git push origin master
```

Connect to remote computer, start MATLAB
```bash
ssh jschap@american.seas.ucla.edu
cd ~/Documents/Codes/VICMATLAB/VICMATLAB/
git pull
RWKDIR=/Users/jschap/Documents/Research/VICGlobal/Data/tuolumne_image_subsetting
cd $RWKDIR
# remove any existing .nc files
rm stehekin_domain.nc
rm stehekin_params.nc
matlab -nodisplay -nosplash
```

Subset VICGlobal domain and parameter files to Stehekin basin mask using the VICMATLAB subsetting functions.
```matlab
addpath(genpath('~/Documents/Codes/VICMATLAB/VICMATLAB/vicmatlab'))
rwkdir = '/Users/jschap/Documents/Research/VICGlobal/Data/tuolumne_image_subsetting';
basinmask = fullfile(rwkdir, 'stehekin_mask.tif');
fulldomain = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/Image/VICGlobal_domain.nc';
fullparams = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/Image/VICGlobal_params.nc';
basinname = 'stehekin';
subset_domain(basinmask, fulldomain, fullfile(rwkdir, [basinname, '_domain.nc']));
subset_parameter2(basinmask, fullparams, fullfile(rwkdir, [basinname, '_params.nc']));
```

Copy the subsetted domain and parameter files back to local
```bash
scp jschap@american.seas.ucla.edu:${RWKDIR}/stehekin_domain.nc ${LWKDIR}/
scp jschap@american.seas.ucla.edu:${RWKDIR}/stehekin_params.nc ${LWKDIR}/
```

Set up global parameter file
```bash
cp global_param.Stehekin.L2015.txt global_param_VICGlobal.txt
xdg-open global_param_VICGlobal.txt
```

Run VIC
```bash 
# Run VIC
./Documents/Software/VIC/vic/drivers/image/vic_image.exe -g ${LWKDIR}/global_param_VICGlobal.txt
```

The above analysis confirms that the model runs using the subsetted VICGlobal parameters. But do the results make sense? How do they compare to the Stehekin results?

 <!-- ----------------------------------------------- -->

## Repeat for Tuolumne Basin

<!-- Steps:
1. Create basin mask for Stehekin basin in R.
1. Switch to command line. Send mask to remote computer
1. Push the latest version of VICMATLAB.
1. Connect to remote computer, start MATLAB
1. Subset VICGlobal domain and parameter files to Stehekin basin mask using the VICMATLAB subsetting functions.
1. Copy the subsetted domain and parameter files back to local
1. Set up global parameter file
1. Run VIC -->

1. Check that the forcing maps exactly match the basin mask.

I checked this out, and in the end, the basin mask from Tuolumne used NaNs where the Stehekin basin mask used 0s. I created a 1's and 0's version of the mask called `upptuo_mask_10.tif.`

The previous analysis seems to imply that the issue has do with the forcing files, since the domain and parameter files worked OK with the Stehekin forcing files. 

1. Send Tuolumne basin mask to remote computer

```bash
LWKDIR=/hdd/ESSD/tuo_image_example
RWKDIR=/Users/jschap/Documents/Research/VICGlobal/Data/tuolumne_image_subsetting
scp /home/jschap/Documents/Codes/VICMATLAB/data/tifs/upptuo_mask_10.tif jschap@american.seas.ucla.edu:${RWKDIR}/ 
```

1. Connect to remote computer, start MATLAB
```bash
ssh jschap@american.seas.ucla.edu
RWKDIR=/Users/jschap/Documents/Research/VICGlobal/Data/tuolumne_image_subsetting
cd $RWKDIR
# remove any existing .nc files
rm tuolumne_domain.nc
rm tuolumne_params.nc
matlab -nodisplay -nosplash
```

1. Subset VICGlobal domain and parameter files to Tuolumne basin mask using the VICMATLAB subsetting functions.
```matlab
addpath(genpath('~/Documents/Codes/VICMATLAB/VICMATLAB/vicmatlab'))
rwkdir = '/Users/jschap/Documents/Research/VICGlobal/Data/tuolumne_image_subsetting';
basinmask = fullfile(rwkdir, 'upptuo_mask_10.tif');
fulldomain = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/Image/VICGlobal_domain.nc';
fullparams = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/Image/VICGlobal_params.nc';
basinname = 'tuolumne';
subset_domain(basinmask, fulldomain, fullfile(rwkdir, [basinname, '_domain.nc']));
subset_parameter2(basinmask, fullparams, fullfile(rwkdir, [basinname, '_params.nc']));
```

Copy the subsetted domain and parameter files back to local
```bash
scp jschap@american.seas.ucla.edu:${RWKDIR}/tuolumne_domain.nc ${LWKDIR}/
scp jschap@american.seas.ucla.edu:${RWKDIR}/tuolumne_params.nc ${LWKDIR}/
```

Set up global parameter file
```bash
xdg-open ${LWKDIR}/tuolumne_image_global.txt
```

Run VIC
```bash 
# Run VIC
LWKDIR=/hdd/ESSD/tuo_image_example
~/Documents/Software/VIC/vic/drivers/image/vic_image.exe -g ${LWKDIR}/tuolumne_image_global.txt
```

Still getting the same error: 
		[ERROR] ../shared_image/src/get_nc_field.c:49: errno: NetCDF: Index exceeds dimension bound: Error getting values for tas in /home/jschap/Documents/Codes/VICMATLAB/data/netcdf_forcings/forc_tuo.2009.nc

OK, now it is time to dig into the VIC source code.

## Looking at VIC source code

The issue comes up in `get_nc_field.c`, in the function `get_nc_field_double`. 

```c
int
get_nc_field_double(nameid_struct *nc_nameid,
                    char          *var_name,
                    size_t        *start,
                    size_t        *count,
                    double        *var)
{
    int status;
    int var_id;

    /* get NetCDF variable */
    status = nc_inq_varid(nc_nameid->nc_id, var_name, &var_id);
    check_nc_status(status, "Error getting variable id for %s in %s", var_name,
                    nc_nameid->nc_filename);

    status = nc_get_vara_double(nc_nameid->nc_id, var_id, start, count, var);
    check_nc_status(status, "Error getting values for %s in %s", var_name,
                    nc_nameid->nc_filename);

    return status;
}
```

This function has five arguments.

`nameid_struct *nc_nameid`, all of which are declared as pointers.

`nameid_struct` is the datatype, and `*nc_nameid` is a pointer. A pointer allocates memory dynamically at runtime. The * is to show that the variable is a pointer variable, not a normal variable. Whereas variables store values, pointer variables store the address of the variable. Pointers are initialized to null (value 0). To get the address of the variable, call it with `&`. To get the value of the variable, call it with `*`. 

First, the function calls `nc_inq_varid`, which takes three arguments. The function call uses the arrow operator `->`, which is used to access elements in structures and unions. `nc_nameid` is a pointer of type structure. The `->` is used to access the variable `nc_id` in `nc_nameid`. 

`get_nc_status` (line 49) is where the VIC code crashed. So `nc_inq_varid` executed successfully. This is a function from the NetCDF library that returns the ID of a netCDF variable, given its name. varname is `tas`, `status` should equal the ID of the `tas` variable in the netCDF forcing file. 

Then, the function `check_nc_status` is called. This is a VIC function that checks whether the mask variable in the input domain is integer type. 

```c
#define check_nc_status(A, M, ...) 
if (A != NC_NOERR)
{
	log_ncerr(A, M, ## __VA_ARGS__);  
	errno = 0; 
	exit(EXIT_FAILURE); 
}
```

Printing status reports in the VIC image driver shows me that the start and count numbers are totally messed up, and the status of the NetCDF file from `nc_get_vara_double` is -40. Start and count get passed into `get_nc_field_double`, so they are wrong from the beginning. Check how they are calculated. It looks like they are calculated in the file `vic_force.c`, in the function `get_forcing_file_info`.

See where `get_nc_field_double` is called. This should give a clue as to why the values for start, count, and var are so strange. It's in line 213 of `vic_driver_shared_image.h`. OK, so the Stehekin forcings have a `_FillValue` field, whereas mine do not. I will re-generate my forcing files to be a better match to the Stehekin forcing files and see if that causes the simulation to run/to not have messed up start/count/var values. I used the forcing preparation code from `wrapper.m` to do this with the VICMATLAB functions `convert_forcing` and `write_netcdf_forcing`.

```matlab
forcdir = '/home/jschap/Documents/Codes/VICMATLAB/data/disagg_forc_2009-2011/';
outname_domain = '/hdd/ESSD/tuo_image_example/forc/upptuo';
start_date = datetime(2009, 1, 1, 0, 0, 0); % need to specify hours
end_date = datetime(2011, 12, 31, 23, 0, 0);
nt_per_day = 24;
prefix = 'full_data_';
precision = 5;
convert_forcing(forcdir, prefix, outname_domain, precision, start_date, end_date, nt_per_day)
```

Need to modify the `write_netcdf_forcing` file to include the proper fill values. Previously, I had

		nccreate(outname,'prcp',...
		    'Datatype','single',...
		    'Dimensions',{'lon',length(lon),'lat',length(lat),'time',ndays*nt_per_day},...
		          'Format','netcdf4_classic')   

But I want to add another field `_FillValue`.

## This bug happened before

OK, so it turns out I had the same bug in September 2019. Looking back through Slack messages to find anything useful.

Error, Sept. 26, 2019. 

"Here are the changes that should be made to the image mode parameter files in order to run VIC with them:
The mask variable in the domain file needs to be of type int. Use 'Datatype', 'int32' when calling nccreate in Matlab to do this.
The variables run_cell, gridcell, fs_active, and Nveg in the parameter file also need to be of type int.
Also the overstory parameter.
Finally, the parameter "fcan" should be named "fcanopy" instead.
Let's hold off on splitting the parameter file into 5 degree tiles. We'll be providing subsetting code, anyway."

But all of these changes have already been applied. And I am still getting the error.

The FillValue did not matter when Dongyue did his streamflow forecasting work. He says that one reason it shouldn't matter is because there is a 0-1 mask (presumably run_cell) telling VIC which cells to use. 

Possibly ask on Gitter? "Hi, is anyone knowledgeable about input file formats for the VIC image driver? I am able to run the image driver with the forcing data provided for the Stehekin example, but not with a different forcing dataset that I made myself." --> Really should come up with a better/more specific question.

## Back to Stehekin

Try creating forcing files for Stehekin. Use the original Stehekin parameter and domain files. This will tell me if there is a problem with my forcing files, specifically. Otherwise, the issue is a combination of the forcing files and parameter/domain files, aka each is OK on their own, but the combination is no good. 

Requires creating forcing files from Livneh, from scratch.

* Subset Livneh NetCDF files to basin

* Run VIC as disaggregator to get forcing files

* Convert from Classic Driver forcing files to Image Driver forcing files. 

```matlab
forcdir = '/home/jschap/Documents/Software/VIC_sample_data/classic/Stehekin/forcings/';
outname_domain = '/hdd/ESSD/tuo_image_example/forc/stehe';
start_date = datetime(1949, 1, 1, 0, 0, 0); % need to specify hours
end_date = datetime(1949, 1, 10, 23, 0, 0);
nt_per_day = 24;
prefix = 'full_data_';
precision = 5;
convert_forcing(forcdir, prefix, outname_domain, precision, start_date, end_date, nt_per_day)
```

I generated NetCDF forcings for Stehekin, and the same error occurred, leading me to strongly believe the error is due to the forcing files, and not due to the parameter/domain files or a combination. 

OK, small breakthrough. Changing the timestamp in the NetCDF forcing file changes the force_skip value slightly. This lends evidence to my hypothesis that the issue is with the format of the timestamp in the forcing file. Perhaps force_skip is very large because VIC has the wrong time origin? Or perhaps it is using seconds instead of hours? Will need to investigate the documentation a bit and try out different setups to find a solution.

I compared my forcing file to the VICGlobal forcing file, and they have identical entries for the time variable. I am using print statements in the source code to see if there is any difference in the way VIC reads in the start dates from these two forcing files. Added `printf("Time is %s", nc_time_origin);` to `vic_force.c` after it gets the date/time of the first entry in the forcing file. Recompiled VIC with `make clean` and `make`. 

OK, the problem is solved! I had the time variable starting from 1, but VIC requires the time variable to start from 0! So instead of 1:365, it should be 0:364.






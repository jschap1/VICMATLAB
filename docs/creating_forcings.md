# Creating meteorological forcing input files

The VIC model requires meteorological forcing input files containing time series of precipitation, temperature, and other meteorological variables at each grid cell. This article describes preparation of meteorological forcing files for the Stehekin domain, using the Livneh et al. (2013) meteorological input data. 

1. Subset precipitation, wind speed, and minimum and maximum temperature data to the basin. 
1. Run VIC-4 as a meteorological disaggregator to calculate temperature, air pressure, vapor pressure, and downwelling shortwave and longwave radiation. 
1. Convert forcing files from classic mode to image mode.

## Subset forcings to Stehekin basin

(This will be done in the ESSD/codes directory, since part of this project is to figure out what is wrong with my image mode forcing files, if anything.)

Get Stehekin domain
Run VICMATLAB subsetting code

See stehekin_wrapper.m.

## Run VIC as disaggregator

Set up and run VIC-4 model for Stehekin. 

Data are in `/hdd/ESSD/data/stehekin`. The Stehekin soil, vegetation, and elevation band parameter files are in the `~/Documents/Software/vic_sample_data` directory. However, these will not work because they are not at the same 1/16 degree resolution as the forcing files. Instead, I am subsetting the VICGlobal soil parameter file.

## Convert from classic to image mode

Run VICMATLAB conversion code

## Run VIC in classic and image mode

Run VIC in classic mode
Run VIC in image mode
Confirm that outputs are the same from both modes
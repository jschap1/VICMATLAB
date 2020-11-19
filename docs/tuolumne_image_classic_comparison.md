# Tuolumne basin image mode/classic mode comparison

Files are in `/home/jschap/Documents/Codes/VICMATLAB/examples/`. Demo the subsetting codes and publish the instructions in a Markdown file, which is available for reference from the VICGlobal Zenodo respository. Compare Classic/Image mode outputs for Tuolumne basin.

Simulation uses the Livneh et al. (2013) meteorological forcing data, and is for the calendar years 2009-2011, with no accounting for soil moisture spin-up. Energy balance mode with five snowbands. 

See MATLAB wrapper for the comparison. Will write this up here as a Markdown file, once it is running. 

<!-- Run the classic mode simulation -->
```bash
./vic_classic.exe -g /home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/global_param.txt
```

<!-- Run the image mode simulation -->
```bash
./vic_classic.exe -g /home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/tuolumne_image_global.txt
```

The basin average and time average values look very similar. However, I checked baseflow for a particular grid cell, and it was pretty different. Maybe there is something to do with the global parameter files? No, it doesn't look like there are any major differences there.

## Setup for classic mode

* Energy balance mode, hourly time steps.
* Frozen soils set to false
* Implicit false
* exp_trans false
* 2 m wind height
* 5 snow bands
* July avgT supplied
* ARNO baseflow formulation

## Setup for image mode

* Energy balance mode, hourly time steps
* Frozen soils set to false
* Wind measurement height read from NetCDF file. Might be 2 m or 10 m
* 5 snow bands
* July avgT supplied
* ARNO baseflow formulation

This makes me think that maybe the soil parameters are different between the image mode and classic mode formats. Look at plots of soil parameters from each source and confirm that they are in fact the same.

OK, I figured out the problem. The tuolumne_image NetCDF parameters use the old rmin and r0 values, and likely use the same albedo, fractional canopy cover, and LAI (non-snow-free), so they don't match the v1.6 vegetation parameters, which I updated this summer. 

## Revisiting with updated image mode parameters

### Subset image mode parameters to the Tuolumne basin

Transfer input files to Atlantic.

```bash

scp /hdd/Data/VICParametersGlobal/VICGlobal/v1.6/Older_Versions/output_latest_aug14.tar.gz jschap@american.seas.ucla.edu:/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/

scp /home/jschap/Documents/Codes/VICMATLAB/data/upptuo_mask.tif jschap@american.seas.ucla.edu:/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/sub/

scp -r /hdd/Data/VICParametersGlobal/VICGlobal/continent_masks jschap@american.seas.ucla.edu:/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/sub/

```

Run subsetting on Atlantic.

```matlab

addpath(genpath('~/Documents/Codes/VICMATLAB/vicmatlab'))

domain_name = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/sub/tuolumne_domain.nc';

param_name = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/sub/tuolumne_params.nc';

basinmaskname = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/sub/upptuo_mask.tif';

global_domain = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/output_latest/VICGlobal_domain.nc';

global_params = '/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/output_latest/VICGlobal_params.nc';    

subset_domain(basinmaskname, global_domain, domain_name)
subset_parameter(basinmaskname, global_params, param_name)

```

Also subset to continents.

```matlab

cd('/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/sub/continent_masks/')
basinmasknames = dir('*.tif');
n = length(basinmasknames)
for k=1:n
	basinmaskname = basinmasknames(k).name;
	param_name = [basinmaskname '_param.nc'];
	domain_name = [basinmaskname '_domain.nc'];
	subset_domain(basinmaskname, global_domain, domain_name)
	subset_parameter(basinmaskname, global_params, param_name)
end

```

Copy the Tuolumne parameters back to Atlantic

```bash
scp jschap@american.seas.ucla.edu:/Volumes/HD3/VICParametersGlobal/VICGlobal/v1.6/sub/*.nc .
```

Run the Tuolumne simulation in image mode.

<!-- Run the image mode simulation -->
```bash
./vic_image.exe -g /home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/tuolumne_image_global.txt > /home/jschap/Documents/Codes/VICMATLAB/data/VICGlobal/vic_image_log.txt
```

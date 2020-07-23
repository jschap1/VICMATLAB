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

## Running the Stehekin example

./vic_image.exe -g /home/jschap/Documents/Software/VIC_sample_data/image/Stehekin/parameters/global_param.Stehekin.L2015.txt

The Stehekin simulation worked fine, although I was not able to write to the /VIC_sample_data/ directory for some reason. Instead, I wrote the result to /home/jschap/Documents/.

Inspecting the Stehekin parameters, the masks need to cover the full domain, but the run_cell parameter should match the forcing files.
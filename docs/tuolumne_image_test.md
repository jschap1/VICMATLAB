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

The last line requires 20 GB of RAM to run, so I exported my input files and code to Atlantic, where I have 48 GB of RAM, and I ran the subsetting there.
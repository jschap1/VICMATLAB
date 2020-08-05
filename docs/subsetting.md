# Subsetting VICGlobal

This documents explains how to use the subsetting codes in VICMATLAB to subset VICGlobal to your study area. The subsetting codes are also included separately with the VICGlobal dataset. 

Written: August 4, 2020 JRS
Dependencies: MATLAB Mapping Toolbox.

1. Download VICGlobal from its [link][vicglobal] on Zenodo
2. Extract the VICGlobal files from the relevant archive (classic or image). Warning: the image mode parameter file is very large after extracting (about 140 GB).  
3. Choose a study area. You should have either a shapefile, a raster, or a list of coordinates defining your study area.

## Classic mode

Subset the soil parameter file to your study basin. Also generates Geotiffs of the soil parameters for easy visualization. 

```matlab
% add path to the subsetting codes
addpath('${your_path}/subsetting/') 

% Load full soil parameter file
soilfile = '${your_path}/Classic/soils_3L_MERIT.txt';
disp('Loading soil parameter file')
soils = load(soilfile);
disp('Soil parameter file has been loaded')

% Subset the soil parameter file
extent = '${your_path}/colo_mask.tif';
grid_decimal = 5; % precision of VIC grid cells
outformat = '3l'; % may be '2l', 3l' or 'livneh'
outname = '${your_path}/colo_soils_VG_v1.6.txt';

generate_tif = 1; % option to generate tifs of soil parameters
setup = '3L-no-org-frost-msds'; % can be 2L, 2L-no-org-fs-july_tavg, 3L, 3L-no-org-frost-msds, or livneh

soils_subset = subset_soils(soils, extent, outname, outformat, grid_decimal, generate_tif, setup);
```

## Image mode

Subset the image mode domain and parameter files to your study basin. Note that this requires a lot of RAM because of the large size of the VICGlobal image mode parameter file. For the Upper Colorado River Basin (7833 1/16 degree grid cells), I found that 16 GB was not enough RAM, but 48 GB is enough. 

```matlab
% add path to the subsetting codes
addpath('${your_path}/subsetting/') 

% Specify output directory and basin name
outdir = '${your_path}/subset_image/';
basinname = 'ucrb';
extent = '${your_path}/colo_mask.tif';

vicglobaldomain = '${your_path}/Image/VICGlobal_domain.nc';
vicglobalparams = '${your_path}/Image/VICGlobal_params.nc';
subset_domain(ucrb_tif, vicglobaldomain, 'ucrb_domain.nc'))
subset_parameter(ucrb_tif, vicglobalparams, 'ucrb_params.nc'));
```

# List of all functions

* subset_domain
* subset_parameter
* geotiffread2
* georefobj2mat
* basin_mask2coordinate_list
* xyz2grid
* get_soil_var_names
* subset_soils
* write_soils

[vicglobal]:https://zenodo.org/record/3475602
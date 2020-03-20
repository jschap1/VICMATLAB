# VICMATLAB
Matlab scripts and functions for processing output files from the Variable Infiltration Capacity (VIC) model for version 4.2.

## Scripts
* vicinputworkflow - creates basin-specific VIC input files from met. forcing data and soil data over a larger domain (CONUS)
* vicoutputworkflow - loads VIC results into Matlab, formats them so they are easy to work with, and generates plots of fluxes and snow states
* correct_flowdir - plots flow direction map, VIC grid cells, basin boundary, stream network, and gauge locations to make it easy to manually correct the flow direction input file
* find_stnloc - finds row and column indices of stream gauges for the station location file
* plotforcings - plots met. forcings (VIC input)
* routoutputworkflow - loads, processes, and makes plots for VIC (Lohmann) routing model outputs

## Functions
* GetCoords - extracts lat/lon coordinates from VIC output file names. Called by vicoutputworkflow.
* LoadVICResults - loads VIC results into Matlab. Called by vicoutputworkflow.
* ProcessVICFluxResults - formats VIC flux outputs. Called by vicoutputworkflow.
* ProcessVICSnowResults - formats VIC snow state outputs. Called by vicoutputworkflow.

## R scripts
* routinputworkflow - creates fraction file, watershed mask, and clipped DEM at DEM and/or VIC modeling resolution. Does not create flow direction file (this is best done with ArcMap).
* modfdir - script for modifying the flow direction file

## MATLAB functions provided in VICMATLAB 1.0
* georefobj2mat
* geotiffread2
* plotraster
* geotiffwrite
* imagescnan
* subset_soils
* write_soils
* subset_forcings
* grid2xyz
* xyz2grid
* check_outputs_wrapper
* get_vic_run_metadata
* get_coordinates_from_VIC_file
* get_vic_header
* check_latlon
* readVIC_ds
* make_outputs_struct
* plot_spatial_avg_ts
* plot_time_avg_maps
* jsplot
* convert_soil_parameters
* get_soil_var_names


## Requirements
* MATLAB R2016a or later
* Parallel Computing Toolbox
* Mapping Toolbox
* Image Processing Toolbox

## Recommended
* R, RStudio
* raster, rgeos, and rgdal packages for R


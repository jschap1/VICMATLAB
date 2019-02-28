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

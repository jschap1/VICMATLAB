## VICMATLAB
Matlab functions for preparing inputs for - and processing outputs from - the Variable Infiltration Capacity (VIC) model. Has functionality for the following VIC versions: VIC-4.2.d, VIC-5 Classic Driver, and VIC-5 Image Driver.

## Tutorial

See the documentation page --- [https://jschap1.github.io/VICMATLAB/](https://jschap1.github.io/VICMATLAB/) --- for a tutorial.

### MATLAB functions provided in VICMATLAB

#### Main functions

* **subset_forcings**: Subsets forcings to a particular basin.

* **load_ascii_forcings**: Loads ASCII meteorological forcing files into MATLAB.

* **get_vic_run_metadata**: Gets metadata, such as variable names, length of simulation, etc. from a VIC simulation. 

* **load_vic_output**: Loads a VIC output variable, such as SWE or net radiation, into MATLAB.

* **subset_soils**: Clips soil parameter file to basin extent.

* **classic2image** (credit to Dongyue Li): Converts VIC-4 or VIC-5 Classic Driver inputs to NetCDF format for use with the VIC-5 Image Driver.

* **plotraster**: Main plotting function to make maps of meteorological variables, VIC outputs, soil parameters, etc.

#### Utility functions

* save_forcings
* georefobj2mat
* geotiffread2
* get_soil_var_names
* write_soils
* get_coordinates_from_VIC_file
* get_vic_header
* GetCoords
* grid2xyz (from MATLAB File Exchange)
* xyz2grid (from MATLAB File Exchange)
* imagescnan (from MATLAB File Exchange)

### Requirements
MATLAB R2016a or later  
Parallel Computing Toolbox  
Mapping Toolbox  
Image Processing Toolbox

### Recommended
R, RStudio  
raster, rgeos, and rgdal packages for R

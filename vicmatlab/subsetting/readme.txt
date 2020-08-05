This directory brings together all the VIC image mode subsetting files to make a unified approach to preparing inputs for the VIC image driver. 

/home/jschap/Documents/Codes/VICMATLAB/vicmatlab/subsetting/
1) subset_domain.m -- subsets the domain file
2) subset_parameter.m -- subsets the parameter file, requires a ton of RAM and doesn't always work
3) subset_netcdf_w_geotiffmask.m -- subsets a NetCDF file to a (geotiff) mask 

/home/jschap/Documents/Codes/VICMATLAB/vicmatlab/utility/
4) classic2image.m -- converts the L2013/L2015 VIC classic parameter files to NetCDF format for the image driver. Does not do subsetting.

/home/jschap/Documents/ESSD/codes/
5) subset_vicglobal_to_extent.m -- wrapper for subsetting the image mode VICGlobal parameters (converted from Classic parameters using Dongyue's classic to image converter) to the basin mask. Would be nice to generalize so it would work with the BV2019 parameters (currently, it gives an error). Calls subset_parameter.m

/home/jschap/Documents/Codes/VICMATLAB/vicmatlab/forcings/
6) write_netcdf_forcing.m -- converts ASCII forcing data to NetCDF forcings
7) make_netcdf_forcing.m -- converts ASCII forcing data to NetCDF forcings. Only handles one year. Other versions are probably better.
8) convert_forcing.m -- converts ASCII forcing data to NetCDF forcings

/home/jschap/Documents/ESSD/codes/Old/
9) subset_netcdf_to_bb.R -- subsets a NetCDF file to a basin boundary. Intended for use with VIC image driver output files.

Desired functionality:
Given classic mode parameter files and a geotiff or shapefile with the basin mask, convert the classic mode inputs to image mode and subset them to the domain. 
Subset ASCII forcings to the domain and convert them to NetCDF format
Ideally works with VICGlobal, L2013, and BV2019 inputs. 

Minimum functionality: 
Subset image mode parameters to the domain
Convert image mode forcing files to NetCDF format


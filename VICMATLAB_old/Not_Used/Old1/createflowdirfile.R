# (Preparing VIC Routing Model Input Files)
#
# DOES NOT WORK (mostly). ATTEMPTING SIMILAR IN QGIS
#
# INPUTS
# Basin mask at the VIC model resolution
# DEM
#
# OUTPUTS (in ASCII GRID format)
# Flow direction file

basinmask.coarse <- raster(file.path("/Users/jschapMac/Desktop/UpperKaweah", 
                                     "basinmask_coarse.asc"))

# Path and filename of DEM
dempath <- "/Users/jschapMac/Documents/HydrologyData/GTOPO30/gt30w140n40_dem"
dem <- raster(file.path(dempath, "gt30w140n40.dem"))

saveloc <- "/Users/jschapMac/Desktop/UpperKaweah/"

# No need to modify below here.
# ----------------------------------------------------------------------

library(rgdal)
library(raster)

### THis bit of code is important, at least, as far as following the instructions 
# on the VIC website goes.
# Clip the DEM using the coarse resolution basin mask
dem.crop <- crop(dem, basinmask.coarse)
basinmask.boostres <- resample(basinmask.coarse, dem.crop, progress="text")
dem.clipped <- dem.crop*basinmask.boostres
writeRaster(dem.clipped, file.path(saveloc, "clipped_dem_coarse"), 
            format = "ascii", overwrite = T)
###

# Use either RSAGA or TauDEM to perform fill, flow direction, and 
# flow accumulation calculations

# OK, so it may be more appropriate to do this in a GIS software, like
# ArcMap or QGIS. I will try it in QGIS.

# The workflow from here goes like this:
# 1. Fill, flowdir, flowacc in ArcMap or QGIS
# 2. Re-map the flow direction numbering so it matches the Stehekin sample basin

# I used the fill tool from SAGA, accessible through QGIS. 
# I set minimum slope to 0.01 degree.





library(RSAGA) # an R library for geospatial terrain analysis. TauDEM is another option.
rsaga.fill.sinks(in.dem = "GTOPO30_clipped", out.dem = "dem.filled.sgrd", 
                 method = "planchon.darboux.2001", minslope = 0.1)

# Based on code from http://ibis.geog.ubc.ca/~rdmoore/Rcode.htm
# http://ibis.geog.ubc.ca/~rdmoore/rcode/ShannonFallsMap.r


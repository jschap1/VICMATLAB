# (Preparing VIC Routing Model Input Files)
#
# INPUTS
# Basin boundary shapefile
# DEM
#
# OUTPUTS (in ASCII GRID format)
# Fraction File at VIC model resolution
# Watershed mask at DEM resolution
# Watershed mask at VIC model resolution
# DEM clipped to basin boundaries

# Note: all coordinate systems should be geographic. I've been using WGS84.

# Basin boundary shapefile; can use ogrListLayers() to get layer name
# ogrListLayers("/Users/jschapMac/Desktop/UpperKaweah")
bb <- readOGR(dsn = "/Users/jschapMac/Desktop/UpperKaweah", 
              layer = "upper_kaweah")

# Path and filename of DEM
dempath <- "/Users/jschapMac/Documents/HydrologyData/GTOPO30/gt30w140n40_dem"
dem <- raster(file.path(dempath, "gt30w140n40.dem"))

saveloc <- "/Users/jschapMac/Desktop/UpperKaweah"

# Name to save clipped raster
rastername <- "GTOPO30_clipped"

fineres <- 30 # DEM resolution (arc-seconds - "fine")
fineres <- fineres/3600

coarseres <- 1/16 # VIC model resolution (degrees - "coarse")

# No need to modify below here.
# ----------------------------------------------------------------------

library(raster)
library(rgdal)

# A user-written function from Carsten Neumann:
# https://stat.ethz.ch/pipermail/r-sig-geo/2013-July/018912.html
clip<-function(raster,shape) {
  a1_crop<-crop(raster,shape)
  step1<-rasterize(shape,a1_crop)
  a1_crop*step1}

# Create a clipped DEM using the basin boundary shapefile
dem.clipped <- clip(dem, bb)
writeRaster(dem.clipped, file.path(saveloc, rastername), 
            format = "ascii", overwrite = T)

# Create a watershed mask at the DEM resolution
dem.clipped.copy <- dem.clipped
dem.clipped.copy[!is.na(dem.clipped.copy)] <- 1

# Use a slightly larger extent than that of the clipped DEM
e <- extent(dem.clipped)
e@xmax <- ceiling(e@xmax)
e@xmin <- floor(e@xmin)
e@ymax <- ceiling(e@ymax)
e@ymin <- floor(e@ymin)

finemask <- raster(resolution = fineres, ext = e)
values(finemask) <- 0

basinmask.fine <- merge(dem.clipped.copy, finemask)
writeRaster(basinmask.fine, file.path(saveloc, "basinmask_fine"), 
            format = "ascii", overwrite = T)

# Make a grid with the VIC model resolution.
# Extent must be larger than the real basin boundaries.

coarsegrid <- raster(resolution = coarseres , ext = e)
for (i in 1:ncell(coarsegrid)) {
  coarsegrid[i] <- i
}
# Must have a unique value in each cell (defines "zones").

# Resample the coarse grid to the DEM resolution, without changing the values
finegrid <- resample(coarsegrid, basinmask.fine, method = "ngb")

# Compute zonal mean
zmean <- zonal(basinmask.fine, z = finegrid, progress = "text")

# Output fraction file
fract <- raster(coarsegrid)
values(fract) <- zmean[,2]
writeRaster(fract, file.path(saveloc, "fract"), 
            format = "ascii", overwrite = T)

# Create a basin mask at the VIC model resolution
basinmask.coarse <- fract
basinmask.coarse[basinmask.coarse != 0] <- 1
writeRaster(basinmask.coarse, file.path(saveloc, "basinmask_coarse"), 
            format = "ascii", overwrite = T)

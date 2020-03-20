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
# If the pixels are not square (in geographic coordinates), this becomes a pain.

library(raster)
library(rgdal)

setwd("/Volumes/HD3/SWOTDA")

# make a shapefile for the IRB basin boundary
irb_mask <- raster("/Volumes/HD4/SWOTDA/Data/IRB/irb_mask.tif")

fract1 <- irb_mask
fract1[is.na(fract1)] <- 0
plot(fract1)
writeRaster(fract, file.path(saveloc, "irb_fract"), 
            format = "ascii", overwrite = T)

irb_bb <- rasterToPolygons(irb_mask)
writeOGR(irb_bb, dsn = "/Volumes/HD4/SWOTDA/Data/IRB/irb_bb_1_16", layer = "bb", driver = "ESRI Shapefile")

fd_masked <- fd*irb_mask
fd_clip <- clip(fd, irb_bb)
fd_crop <- crop(fd, irb_bb)

# fd <- raster("/Volumes/HD4/SWOTDA/Data/IRB/flow_directions_vic_d4_masked.tif")

fd_latest <- raster("/Volumes/HD3/SWOTDA/FDT/v12032019/fd_grass_d5_masked.tif")

writeRaster(fd_latest, file = "/Volumes/HD3/SWOTDA/FDT/v12032019/fd_grass_d5_masked_v2.tif", 
            NAflag = 0, datatype = "INT2S", overwrite=TRUE)


writeRaster(fd_masked, filename = "/Volumes/HD4/SWOTDA/Data/IRB/flow_directions_vic_d4_masked.tif")

# Basin boundary shapefile; can use ogrListLayers() to get layer name
# ogrListLayers("/Users/jschap/Desktop/UMRB/GIS")
# bb <- readOGR(dsn = "/Volumes/HD3/SWOTDA/Data/UMRB/Basin", 
#               layer = "bb_merc")
# wgs84 <- "+init=epsg:4326"
# bb <- spTransform(bb, crs(wgs84))

# bb <- readOGR("./Delineation/basin_3as.shp")
bb <- readOGR("./FDT/v10282019/pandoh.shp")
bb <- readOGR("./FDT/v10282019/bhakra.shp")
bb <- readOGR("./FDT/v10282019/bhakra.shp")
# bb <- readOGR("./FDT/v10282019/tarbela.shp")
# bb <- readOGR("./FDT/v10282019/mangla.shp")

# Path and filename of DEM
# dempath <- "/Volumes/HD3/SWOTDA/Data/UMRB/DEM/HydroSHEDS"
# dem <- raster(file.path(dempath, "clipped_DEM.tif"))

# dem <- raster("/Volumes/HD2/MERIT/Merged_IRB/merged_merit_30as.tif")

# doesn't have to be a DEM, any fine resolution raster will do.

dem <- raster("./FDT/merged_merit_30as.tif")

saveloc <- "/Volumes/HD4/SWOTDA/Data/IRB/"

fineres <- res(dem)[1] # DEM resolution (arc-seconds - "fine")
# fineres <- fineres/3600

coarseres <- 1/16 # VIC model resolution (degrees - "coarse")

# Make a shapefile from the text file

gages <- read.table("./FDT/Rout/smaller/gage_coordinates_revised.txt")
names(gages) <- c("dam","river","lon","lat","V5")
coordinates(gages) <- ~lon+lat
crs(gages) <- crs(dem)
writeOGR(gages, dsn = "./FDT/Rout/smaller/gage_coordinates_revised.shp", layer = "gages", driver = "ESRI Shapefile")

# No need to modify below here.
# ----------------------------------------------------------------------

# A user-written function from Carsten Neumann:
# https://stat.ethz.ch/pipermail/r-sig-geo/2013-July/018912.html
clip<-function(raster,shape) {
  a1_crop<-crop(raster,shape)
  step1<-rasterize(shape,a1_crop)
  a1_crop*step1}

# Create a clipped DEM using the basin boundary shapefile
dem.clipped <- clip(dem, bb)
writeRaster(dem.clipped, file.path(saveloc, "MERIT_30as_clipped"), 
            format = "ascii", overwrite = T)

# writeRaster(dem.clipped, file.path(saveloc, "bhakra_dem_30as"), 
#             format = "ascii", overwrite = T)

# Also save a version that is cropped to the basin extent
# dem.cropped <- crop(dem, bb)
# writeRaster(dem.cropped, file.path(saveloc, "GTOPO30_cropped"), 
#             format = "ascii", overwrite = T)

# Create a watershed mask at the DEM resolution
dem.clipped.copy <- dem.clipped
dem.clipped.copy[!is.na(dem.clipped.copy)] <- 1

# Use a slightly larger extent than that of the clipped DEM
e <- extent(dem.clipped)
e@xmax <- ceiling(e@xmax)
e@xmin <- floor(e@xmin)
e@ymax <- ceiling(e@ymax)
e@ymin <- floor(e@ymin)

finemask <- raster(resolution = fineres, ext = e, crs = crs(bb))
values(finemask) <- 0

basinmask.fine <- merge(dem.clipped.copy, finemask)
writeRaster(basinmask.fine, file.path(saveloc, "basinmask_fine"), 
            format = "ascii", overwrite = T)

# Make a grid with the VIC model resolution.
# Extent must be larger than the real basin boundaries.

# Must have a unique value in each cell (defines "zones").
coarsegrid <- raster(resolution = coarseres , ext = e, crs = crs(bb))
values(coarsegrid) <- 1:ncell(coarsegrid)

# Resample the coarse grid to the DEM resolution, without changing the values
finegrid <- resample(coarsegrid, basinmask.fine, method = "ngb")

# Compute zonal mean
zmean <- zonal(basinmask.fine, z = finegrid, progress = "text")

# Output fraction file
fract <- raster(coarsegrid)
values(fract) <- zmean[,2]
writeRaster(fract, file.path(saveloc, "irb_fract"), 
            format = "ascii", overwrite = T)

# Create a basin mask at the VIC model resolution
basinmask.coarse <- fract
basinmask.coarse[basinmask.coarse != 0] <- 1
writeRaster(basinmask.coarse, file.path(saveloc, "mangla_basinmask_coarse"), 
            format = "ascii", overwrite = T)
writeRaster(basinmask.coarse, file.path(saveloc, "mangla_basinmask_coarse.tif"), overwrite = T)

# Check outputs
plot(basinmask.coarse)
plot(fract)

# Load Wu et al. (2012) flow direction data

flowdir <- raster("../DRT/upscaled_global_hydrography/by_HYDRO1K/by_HYDRO1K/WGS84/flow_direction/DRT_FDR_globe_ARCGIS_qd.asc")
fd.crop <- crop(flowdir, basinmask.coarse)
origin(fd.crop) <- c(0,0) # sketchy...
fd <- fd.crop*basinmask.coarse
# plot(fd)
# fd
# basinmask.coarse
# origin(basinmask.coarse)
crs(fd) <- crs(basinmask.coarse)
writeRaster(fd, "./FDT/Rout/flowdir_wu.tif")


## Change flow directions from VIC convention to GRASS convention

fd_copy <- flowdir
fd_copy[flowdir == 1] <- 2
fd_copy[flowdir == 2] <- 1
fd_copy[flowdir == 3] <- 8
fd_copy[flowdir == 4] <- 7
fd_copy[flowdir == 5] <- 6
fd_copy[flowdir == 6] <- 5
fd_copy[flowdir == 7] <- 4
fd_copy[flowdir == 8] <- 3
writeRaster(fd_copy, file = "./Data/IRB/Experimental/flow_directions_grass.tif", 
            NAflag = 0, datatype = "INT2S", overwrite=TRUE)

# Crop MERIT DEM to exact extent used for topographical forcing downscaling
elev_fine <- raster('/Volumes/HD2/MERIT/Merged_1_16/merged_merit_dem_1_16.tif')
t0_example <- raster("./Data/IRB/Experimental/temperature_example.tif")
crs(t0_example) <- crs(elev_fine)
elev_fine_crop <- crop(elev_fine, t0_example)
plot(elev_fine_crop)
writeRaster(elev_fine_crop, file = "./Data/IRB/Experimental/cropped_DEM_for_downscaling_1_16.tif")

# Crop the coarse basin mask to the exact extent used for topographical forcing downscaling
basinmask <- raster("/Volumes/HD3/SWOTDA/Data/IRB/Experimental/fract.asc")
basinmask[basinmask>=0] <- 1
crs(basinmask) <- crs(elev_fine)
basinmask_crop <- crop(basinmask, t0_example)
plot(basinmask_crop)
writeRaster(basinmask_crop, file = "./Data/IRB/Experimental/cropped_basinmask_for_downscaling_1_16.tif")

# Slope, aspect, and curvature calculations for wind speed downscaling
slope1 <- raster("./Data/IRB/Experimental/slope1.tif")
aspect1 <- raster("./Data/IRB/Experimental/aspect1.tif")
curv1 <- raster("./Data/IRB/Experimental/pcurv1.tif")
slope1_crop <- crop(slope1, elev_fine_crop)
aspect1_crop <- crop(aspect1, elev_fine_crop)
curv1_crop <- crop(curv1, elev_fine_crop) # curvature is close to 0 everywhere
writeRaster(slope1_crop, file = "./Data/IRB/Experimental/cropped_slope_for_downscaling_1_16.tif", 
            overwrite = TRUE, NAflag = -9999)
writeRaster(aspect1_crop, file = "./Data/IRB/Experimental/cropped_aspect_for_downscaling_1_16.tif", 
            overwrite = TRUE, NAflag = -9999)
writeRaster(curv1_crop, file = "./Data/IRB/Experimental/cropped_curvature_for_downscaling_1_16.tif", 
            overwrite = TRUE, NAflag = -9999)
plot(aspect1_crop)
plot(basin, add=TRUE)



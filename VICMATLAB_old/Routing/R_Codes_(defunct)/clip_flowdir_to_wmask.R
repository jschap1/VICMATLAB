library(rgdal)
library(raster)

# Clips Wu et al. flow direction file to basin boundary
#
# INPUTS
# Basin boundary shapefile
# Flow direction raster
# Watershed mask
#
# OUTPUTS (in ASCII GRID format)
# DEM clipped to basin boundaries

# Note: all coordinate systems should be geographic. I've been using WGS84.

# Watershed mask at 1/16 resolution
wmask <- raster("/Volumes/HD3/SWOTDA/Data/IRB/basinmask_coarse.asc")
crs(wmask) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# Path and filename of DEM
flowdirpath <- "/Volumes/HD3/DRT/upscaled_global_hydrography/by_HYDRO1K/by_HYDRO1K/WGS84/flow_direction"
flowdir <- raster(file.path(flowdirpath, "DRT_FDR_globe_ARCGIS_16th.asc"))
crs(flowdir) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

saveloc <- "/Volumes/HD3/SWOTDA/Data/IRB"

flowdir.cropped <- crop(flowdir, wmask)

# Forcing the origin to be the same so I can multiply by the wmask raster
fd <- resample(flowdir.cropped, wmask, method="ngb", progress = "text")

flowdir.clipped <- fd*wmask

# Convert flow direction from ArcGIS to VIC convention
arcmap2vic <- function(fd)
{
  fd[fd==64] <- 1
  fd[fd==128] <- 2
  fd[fd==1] <- 3
  fd[fd==2] <- 4
  fd[fd==4] <- 5
  fd[fd==8] <- 6
  fd[fd==16] <- 7
  fd[fd==32] <- 8
  return(fd)
}
flowdir.clipped.vic <- arcmap2vic(flowdir.clipped)

# Convert flow direction from ArcGIS to GRASS convention
arcmap2grass <- function(fd)
{
  fd[fd==128] <- 1
  fd[fd==64] <- 2
  fd[fd==32] <- 3
  fd[fd==16] <- 4
  fd[fd==8] <- 5
  fd[fd==4] <- 6
  fd[fd==2] <- 7
  fd[fd==1] <- 8
  return(fd)
}
flowdir.clipped.vic <- arcmap2vic(flowdir.clipped)

writeRaster(flowdir.clipped.vic, file.path(saveloc, "drt_flowdir_clipped_vic.asc"), 
            format = "ascii", overwrite = T, NAflag = 0)


# Remove river polylines where the width is less than 100 m
keep.ind <- which(rivs$a_WIDTH >= 100)
rivs2 <- rivs[keep.ind,]
plot(rivs2)
writeOGR(rivs2, dsn = saveloc, layer = "rivs", driver = "ESRI Shapefile")

# Keep river polylines with a sufficiently large upstream drainage area
keep.ind <- which(rivs$a_AREA >= 5000)
rivs3 <- rivs[keep.ind,]
plot(rivs3)

# Convert format of flow direction file
fdir_drt <- raster("./Data/IRB/ROUT/irb.flowdir.asc")
crs(fdir_drt) <- "+init=epsg:4326"
fdir_mod <- raster("./Data/IRB/ROUT/mod_fdir.asc")
fdir_out <- resample(fdir_mod, fdir_drt, method = "ngb")
# fdir2 <- 
writeRaster(fdir, file = "./Data/IRB/ROUT/irb.flowdir_corrected.asc", overwrite = TRUE, NAflag = 0)


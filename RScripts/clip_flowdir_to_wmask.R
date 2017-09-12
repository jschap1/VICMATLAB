library(rgdal)
library(raster)

# Clips Wu et al. flow direction file to basin boundary
#
# INPUTS
# Basin boundary shapefile
# Flow direction raster
#
# OUTPUTS (in ASCII GRID format)
# DEM clipped to basin boundaries

# Note: all coordinate systems should be geographic. I've been using WGS84.

# Watershed mask at 1/16 resolution
wmask <- raster("/Users/jschapMac/Desktop/Tuolumne/RoutingInputs/basinmask_coarse.asc")
crs(wmask) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# Path and filename of DEM
flowdirpath <- "/Users/jschapMac/Documents/HydrologyData/DRT/data/DRT/upscaled_global_hydrography/by_HYDRO1K/by_HYDRO1K/WGS84/flow_direction"
flowdir <- raster(file.path(flowdirpath, "DRT_FDR_globe_ARCGIS_16th.asc"))
crs(flowdir) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

saveloc <- "/Users/jschapMac/Desktop/Tuolumne5/LFP"

flowdir.cropped <- crop(flowdir, wmask)

# Forcing the origin to be the same so I can multiply by the wmask raster
fd <- resample(flowdir.cropped, wmask, method="ngb")

flowdir.clipped <- flowdir.cropped*wmask
writeRaster(flowdir.clipped, file.path(saveloc, "flowdir_clipped.asc"), 
            format = "ascii", overwrite = T)

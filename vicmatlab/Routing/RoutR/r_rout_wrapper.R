# R_Rout wrapper
# 
# Executive file for the R_Rout R package

# Maintaining the RoutR package ----------------------------------------------

library(devtools)
library(roxygen2)

parent_dir <- "/Users/jschap/Documents/Codes/VICMATLAB/Routing/"
# create("RoutR")

setwd(file.path(parent_dir, "RoutR"))
document()
setwd("..")
install("RoutR")

library(RoutR)
library(raster)
library(rgdal)

# Load an existing flow direction file ---------------------

fd <- raster("/Volumes/HD3/SWOTDA/FDT/v12032019/fd_grass_d5_masked.tif")
plot(fd, main = "Flow directions (GRASS convention)")

# Write out a flow direction file for input into GRASS GIS ---------------------

fd[fd==0] <- NA
writeRaster(fd, file = "/Volumes/HD3/SWOTDA/FDT/v12032019/fd_grass_d5_masked_v2.tif", 
            NAflag = 0, datatype = "INT2S", overwrite=TRUE)

# Create a shapefile with gage locations -----------------------------------------

gages <- read.table("/Volumes/HD3/SWOTDA/FDT/Rout/smaller/gage_coordinates_revised.txt")
names(gages) <- c("dam","river","lon","lat","V5")
outname <- "/Volumes/HD3/SWOTDA/FDT/v12032019/gage_coordinates_revised.shp"
gages_spdf <- coordlist2shape(gages, outname, fd)

# Create a fraction file ----------------------------------------------------------

dem <- raster("/Volumes/HD4/SWOTDA/Data/IRB/irb_dem.tif")
bb <- readOGR("/Volumes/HD4/SWOTDA/Data/IRB/irb_bb_1_16/bb.shp")
saveloc <- "/Volumes/HD3/SWOTDA/FDT/v12032019/Routing_Inputs"
basename <- "irb"
fract <- make_fract(dem, bb, saveloc, basename, fineres = 1/16, coarseres = 1/16)

# Create a trivial fraction file -------------------------------------------------

basinmask <- raster("/Volumes/HD3/SWOTDA/FDT/v12032019/fd_grass_d6_masked.tif")
basinmask[basinmask==0] <- NA
basinmask[basinmask>0] <- 1
saveloc <- "/Volumes/HD3/SWOTDA/FDT/v12032019/Routing_Inputs"
basename <- "irb"
fract1 <- make_fract_trivial(basinmask, saveloc, basename)

# Make a basin boundary shapefile based on a basin mask --------------------------

irb_bb <- rasterToPolygons(irb_mask)
writeOGR(irb_bb, dsn = "/Volumes/HD4/SWOTDA/Data/IRB/irb_bb_1_16", layer = "bb", driver = "ESRI Shapefile")

# Use a basin mask to subset a flow direction file ------------------------------
flowdir <- raster("/Volumes/HD3/SWOTDA/FDT/v10282019/flow_direction_vic_d4_00.asc")
basinmask <- raster("/Volumes/HD3/SWOTDA/FDT/v10282019/mangla_basinmask_coarse.tif")
savename <- "/Volumes/HD3/SWOTDA/Calibration/Mangla/mangla_flowdir_vic"
fd <- clip_flowdir(flowdir, basinmask, savename, out_format = "ascii")

# Write out x and y coordinates of the domain ----------------------------------

basinmask <- raster("/Volumes/HD4/SWOTDA/Data/IRB/irb_mask.tif")
saveloc <- "/Volumes/HD4/SWOTDA/Data/IRB/"
domain_points <- get_domain_points(basinmask, saveloc)

# Make station location file ----------------------------------------------------

r = raster("/Volumes/HD3/SWOTDA/FDT/v12032019/Routing_Inputs/irb_fract_trivial.asc")
gage_coords = "/Volumes/HD3/SWOTDA/FDT/v12032019/gage_coords_in.txt"
saveloc <- "/Volumes/HD3/SWOTDA/FDT/v12032019/Routing_Inputs"
basename <- "irb5"
make_stnloc(r, gage_coords, basename, saveloc)

# Create a flow direction file for the UW routing model --------------------------




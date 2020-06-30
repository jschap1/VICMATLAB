# File conversions before running subset_image_wrapper.m
#
# 6/5/2020 JRS

library(raster)

ut <- raster("/home/jschap/Documents/Codes/VICMATLAB/vicmatlab/subsetting/subset_image/upper_tuolumne_basin_merc.tif")
ucrb <- raster("/home/jschap/Documents/ESSD/data/colo_mask.tif")

ut
ucrb[is.na(ucrb)] <- 0
plot(ucrb)

ut.wgs <- projectRaster(ut, crs = "+init=epsg:4326", method = "ngb", res = 1/16)
origin(ut.wgs) <- c(0,0) # not great practice...
writeRaster(ut.wgs, filename = "/home/jschap/Documents/Codes/VICMATLAB/vicmatlab/subsetting/subset_image/upper_tuolumne_basin_wgs.tif",
            overwrite = TRUE)

# ut.wgs <- projectRaster(ut, crs = "+init=epsg:4326", method = "ngb", alignOnly = TRUE)
# ut.wgs <- projectRaster(ut, ucrb, method = "ngb")
# ut_extent <- projectExtent(ut, crs = "+init=epsg:4326")
# ut.wgs <- crop(ut.wgs, ut_extent)
#
# plot(ut.wgs)

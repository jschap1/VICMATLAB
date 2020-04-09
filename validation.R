# Comparing inputs and outputs from VIC to ensure that VICMATLAB is working properly
#
# 4/4/2020 JRS

library(raster)
library(rgdal)
setwd("~/Documents/Research/Codes/VICMATLAB/")

# We are looking at precipitation data from Livneh (2013)
# Values are averaged over the calendar years 2009-2011

# Original Livneh data
precip1 <- raster("./data/validation/livneh_precipitation_2009-2011_average.tif")
plot(precip1)

  # Downscaled Livneh data
precip2 <- raster("./data/validation/livneh_precipitation_downscaled_2009-2011_average.tif")
plot(precip2)

# Output data from VIC
precip3 <- raster("./data/validation/precipitation_output_2009-2011_average.tif")
plot(precip3)

# Everything is in agreement! Fixed the issue where the VIC outputs were getting all jumbled up.
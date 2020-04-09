# Compare Livneh et al. (2013) soil depths with other soil depth data

library(raster)
library(rgdal)

d1 <- raster("/Volumes/HD3/VICParametersCONUS/soil_data/depth1.tif")
d2 <- raster("/Volumes/HD3/VICParametersCONUS/soil_data/depth2.tif")
d3 <- raster("/Volumes/HD3/VICParametersCONUS/soil_data/depth3.tif")
d_tot <- d1+d2+d3

e <- extent(d_tot)

p2016 <- raster("/Volumes/HD3/Soil_Depth/Global_Soil_Regolith_Sediment_1304/data/average_soil_and_sedimentary-deposit_thickness.tif")
p2016_crop <- crop(p2016, e)

# Handle no-data values
# units are meters
# Values of 50 m are actually greater than or equal to 50 m
vals <- values(p2016_crop)
vals[vals == 255] <- NA
values(p2016_crop) <- vals

# Resample to 1/16 degrees

p2016_res <- resample(p2016_crop, d_tot, progress = "text")

# Mask out the area we're not interested in 

p2016_mask <- mask(p2016_res, d_tot)

plot(p2016_mask, main = "Pelletier et al. (2016) soil thickness (m)")
plot(d_tot, main = "Livneh et al. (2013) soil thickness (m)")


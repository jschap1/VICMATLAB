area.rast <- dsmax
values(area.rast) <- 1
plot(area.rast)

area.vals <- area(area.rast)
writeRaster(area.vals, filename = "/Volumes/HD3/VICParametersGlobal/Global_1_16/landmask/area_raster_km2.tif")

library(rgdal)
library(raster)
merit <- raster("/Volumes/HD2/MERIT/DEM/Merged_1_16/merged_merit_dem_1_16.tif")
plot(merit)
plot(irb, add=TRUE)
e <- drawExtent()
cropped.dem <- crop(merit, e)
plot(cropped.dem)
irb <- readOGR("/Volumes/HD3/SWOTDA/Data/IRB/irb_bb.shp")
plot(irb, add = TRUE)
writeRaster(cropped.dem, "/Volumes/HD4/SWOTDA/Data/IRB/VIC/dem_really_extended.tif")

# Find a streamgauge near the basin outlet

library(raster)
library(rgdal)
elev.tuo <- raster("/Volumes/HD4/SWOTDA/Data/Tuolumne/dem.tif")
plot(elev.tuo, main = "Tuolumne River Basin")
bb <- readOGR("/Users/jschap/Documents/Research/SWOTDA_temp/Tuolumne/Tuolumne1/Shapefiles/upper_tuolumne_wgs.shp")
plot(bb, add = TRUE)
gages <- readOGR("/Users/jschap/Documents/Classes/GEOG207 Seminar/Data/QMIN/Stream Gage QC/basin_gage_shapefiles/Geog207_station_list.shp")
wgs84 <- crs(bb)
gages.wgs <- spTransform(gages, wgs84)
plot(gages.wgs, add = TRUE, pch = 19, col = "darkgreen", cex = 1.3)

# gages.in.tuo <- extract(elev.tuo, gages.wgs)
gages.in.tuo <- intersect(gages.wgs, elev.tuo)
relevant.gage.data <- "/Users/jschap/Documents/Classes/GEOG207 Seminar/Data/QMIN/streamflow_data/archived_data/data/11274790.txt"

library(sp)
gageloc <- data.frame(lat = 37.91659, lon = -119.6599)
coordinates(gageloc) <- ~lon+lat
# gageloc <- c(-119.6599, 37.91659)
# coords(gageloc)
gageloc.sp <- SpatialPoints(gageloc, proj4string = wgs84)
plot(bb)
points(gageloc.sp, pch = 19, col = "red")


A <- area(elev.tuo)
rect.area <- sum(values(A))


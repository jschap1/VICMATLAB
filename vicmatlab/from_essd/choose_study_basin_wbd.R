library(rgdal)
library(raster)
library(RoutR)

# Upper Tuolumne River Basin ------------------------------

dat <- readOGR("/Volumes/HD3/WBD/wbdhu8_a_us_september2016.gdb")
plot(dat)
ca_ind <- which(dat$STATES == "CA")
ca <- dat[ca_ind,]

ca$Shape_Area

tuo_ind <- which(ca$NAME == "Upper Tuolumne")  
tuo <- ca[tuo_ind,]
area(tuo)/1e6

tuo <- tuolumne

#tuo_dem <- clip_raster(merit_dem, tuo)
#a1_crop <- crop(merit_dem, tuo)
#step1 <- rasterize(shape, a1_crop)
#a1_crop * step1

cls <- cellFromPolygon(merit_dem, tuo, weights = TRUE)[[1]][,"cell"]
merit_dem_copy <- merit_dem
merit_dem[][-cls] <- NA
tuo_dem <- trim(merit_dem)
plot(tuo_dem)
plot(tuo, add=TRUE)

# compare area of true and model study area
area(tuo)/1e6 # 4850 km^2
ind <- !is.na(values(tuo_dem)) # find cells in domain
area_vect <- values(area(tuo_dem))
sum(area_vect[ind]) # 6650 km^2
# The discrepancy will be less important for larger basins

# Upper Mississippi River Basin ------------------------------

huc2 <- readOGR("/Volumes/HD3/WBD/HUC2/Upp_Miss_07/Shape/WBDHU2.shp")
plot(huc2)

umrb <- huc2

area(umrb)/1e6

# Use WBD polygon to make a basin mask
merit_mask <- raster("/Volumes/HD2/MERIT/DEM/Merged_1_16/merit_mask_1_16.tif")
merit_dem <- raster("/Volumes/HD2/MERIT/DEM/Merged_1_16/merged_merit_dem_1_16.tif")

skagit <- readOGR("/Users/jschap/Documents/Research/Glaciers/Skagit/boundary/SkagitRiver_BasinBoundary.shp")
target_mask <- clip_raster(merit_dem, tuo)


# umrb_mask <- clip_raster(merit_mask, huc2)
plot(target_mask)
# plot(skagit, add = TRUE)
writeRaster(target_mask, "/Users/jschap/Documents/Research/Glaciers/Skagit/skagit_dem.tif",
            NAflag = 0, datatype = "INT2S")

# Upper Colorado River Basin ------------------------------

huc2 <- readOGR("/Volumes/HD3/WBD/wbdhu2_a_us_september2019.gdb")
upp_colo <- huc2[4,]
colo_mask <- clip_raster(merit_mask, upp_colo)
plot(colo_mask)
writeRaster(colo_mask, "/Users/jschap/Documents/Research/VICGlobal/Data/Colorado/colo_mask.tif",
            NAflag = 0, datatype = "INT2S")

# Make figure for presentation -----------------------------

library(maps)

map("state", boundary = TRUE)
plot(umrb, add = TRUE, border = "darkblue")
plot(upp_colo, add = TRUE, border = "darkblue")

# Find stream gages ---------------------------------------------------

gages2 <- readOGR("/Volumes/HD3/GAGES-II/gagesII_9322_point_shapefile/gagesII_9322_sept30_2011_wgs.shp")

plot(umrb)
plot(rivers, add = TRUE, col = "blue")
e <- drawExtent()
umrb_sub <- crop(umrb, e)
plot(umrb_sub)
plot(gages2, add = TRUE, pch = 19, cex = 0.5)

rivers <- readOGR("/Volumes/HD3/NaturalEarthData/ne_50m_rivers_lake_centerlines/ne_50m_rivers_lake_centerlines.shp")
plot(rivers, add = TRUE, col = "blue")



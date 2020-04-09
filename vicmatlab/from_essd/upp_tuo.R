library(rgdal)

dat <- readOGR("/Volumes/HD3/WBD/wbdhu8_a_us_september2016.gdb")
plot(dat)
ca_ind <- which(dat$STATES == "CA")
ca <- dat[ca_ind,]

ca$Shape_Area

tuo_ind <- which(ca$NAME == "Upper Tuolumne")  
tuo <- ca[tuo_ind,]

# Upper Mississippi River Basin

huc2 <- readOGR("/Volumes/HD3/WBD/HUC2/Upp_Miss_07/Shape/WBDHU2.shp")
plot(huc2)

library(raster)
area(huc2)/1e6

# Upper Colorado River Basin
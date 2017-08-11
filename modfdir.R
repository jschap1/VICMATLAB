library(raster)
library(rgdal)

setwd("C:\\Users\\Jacob\\Documents\\My Documents\\Research")

fdir <- raster("flowdir")

# Coordinates of the flow direction cells to change
lon <- c(-120.952956, -120.838742, -120.837473)
lat <- c(37.595848, 37.593309, 37.650417)
coords <- cbind(lon, lat)

ind <- extract(fdir,SpatialPoints(coords), cellnumbers=TRUE) 

modfdir <- fdir

# Specify new values
modfdir[ind[1,1]] <- 16
modfdir[ind[2,1]] <- 16
modfdir[ind[3,1]] <- 16

writeRaster(modfdir, "modfdir", format = "ascii")

## Convert numbering system from ArcMap to that 
## used by the Lohmann routing model

modfdir[modfdir == 1] <- 3
modfdir[modfdir == 2] <- 4
modfdir[modfdir == 4] <- 5
modfdir[modfdir == 8] <- 6
modfdir[modfdir == 16] <- 7
modfdir[modfdir == 32] <- 8
modfdir[modfdir == 64] <- 1
modfdir[modfdir == 128] <- 2

fdir_in <- modfdir
writeRaster(fdir_in, "fdir_in", format = "ascii")




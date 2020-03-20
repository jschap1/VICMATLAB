# Convert soil parameter file to geotiffs

## Do it with the whole Livneh soil parameter file ----------------------------------------------------------------------

soils <- read.table("/Volumes/HD3/VICParametersCONUS/vic.soil.0625.new.cal.adj.conus.plus.crb.can_no_July_T_avg.txt", stringsAsFactors = FALSE)

basin <- raster("/Users/jschap/Documents/Codes/VICMATLAB/data/upptuo_mask.tif")
dem <- raster("/Users/jschap/Documents/Codes/VICMATLAB/data/upptuo_dem.tif")

lat <- soils[,3]
lon <- soils[,4]

xyz <- cbind(lon, lat, soils[,22])
elev <- rasterFromXYZ(xyz, crs = "+init=epsg:4326")
plot(elev, main = "Elevation")

xyz <- cbind(lon, lat, soils[,5])
b <- rasterFromXYZ(xyz, crs = "+init=epsg:4326")
plot(b, main = "VIC parameter")
 
## OK, now try with the irregular data from the clipped soil parameter file ---------------------------------------------

soils_tuo <- read.table("/Users/jschap/Documents/Codes/VICMATLAB/data/upptuo_soils_livneh.txt", stringsAsFactors = FALSE)
lat <- soils_tuo[,3]
lon <- soils_tuo[,4]
xyz <- cbind(lon, lat, soils_tuo[,22])
elev <- rasterFromXYZ(xyz, crs = "+init=epsg:4326")
plot(elev, main = "Elevation") # elevation has been flipped upside-down. This is likely a consequence of the soil-parameter subsetting function.


plot(dem)


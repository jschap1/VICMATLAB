
fd <- raster("tuolumne_flowdir.tif")
plot(fd)

felev <- raster("tuolumne_elev.tif")
dem <- raster("tuolumne_dem.tif")

# Output the X and Y coordinates of the domain:

spts <- rasterToPoints(fd, spatial = TRUE)
coords <- as.data.frame(spts)

x <- unique(coords[,2])
y <- unique(coords[,3])

write.csv(x, file="x.csv")
write.csv(y, file="y.csv")

#' Clip raster
#' 
#' Clips a raster to a shapefile
#' A user-written function from Carsten Neumann:
#' https://stat.ethz.ch/pipermail/r-sig-geo/2013-July/018912.html
#' @export

clip_raster <- function(raster,shape) {
  a1_crop<-crop(raster,shape)
  step1<-rasterize(shape,a1_crop)
  a1_crop*step1}

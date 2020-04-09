#' coordlist2shape
#' 
#' Make shapefile from list of coordinates
#' @param coordlist data frame containing lat and lon values for each station
#' @param outname savename for the output shapefile
#' @param r raster or vector object with the desired output coordinate system
#' @export

coordlist2shape <- function(gages, outname, r)
{
  coordinates(gages) <- ~lon+lat
  crs(gages) <- crs(r)
  writeOGR(gages, dsn = outname, layer = "gages", driver = "ESRI Shapefile", overwrite = TRUE)
  print(paste("Saved", outname))
  return(gages)
}



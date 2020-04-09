#' Get raster coordinates
#' 
#' Gets coordinates of a raster and optionally saves them as csv files
#' @param r raster
#' @param saveloc location to (optionally) save raster pixel coordinates
#' @value
#' x and y coordinates of each point in the raster
#' @export
#' @examples 
#' basinmask <- raster("/Volumes/HD4/SWOTDA/Data/IRB/irb_mask.tif")
#' saveloc <- "/Volumes/HD4/SWOTDA/Data/IRB/"
#' domain_points <- get_domain_points(basinmask, saveloc)

get_raster_coords <- function(r, saveloc = FALSE)
{
  spts <- rasterToPoints(fd, spatial = TRUE)
  coords <- as.data.frame(spts)
  
  x <- unique(coords[,2])
  y <- unique(coords[,3])
  
  if (saveloc != FALSE)
  {
    write.csv(x, file=file.path(saveloc, "x.csv"))
    write.csv(y, file=file.path(saveloc, "y.csv"))
  }
  
  cc <- data.frame(x = x, y = y)
  return(cc)
  
}

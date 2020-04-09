#' Subset Study Area from GeoTiff
#'
#' Subsets soil or vegetation parameters from a larger GeoTiff file to the basin mask of some other study area
#' @param r soil or vegetation parameter raster
#' @param basin_mask basin mask
#' @examples 
#' basin_mask <- raster("/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif")
#' r <- raster("/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_4/Figures/tifs/jan_lai_nh.tif")
#' subset_geotiff(r, basin_mask)
#' @export

subset_geotiff <- function(r, basin_mask, savename = NULL)
{

  r2 <- resample(r, basin_mask)
  
  dat <- values(r2)
  
  dat[is.na(values(basin_mask))] <- NA
  
  r_clip <- r2
  values(r_clip) <- dat
  
  if (!is.null(savename))
  {
    writeRaster(r_clip, savename)
  }
  
  return(r_clip)
  
}

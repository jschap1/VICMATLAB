#' Make stnloc
#' 
#' Makes a station location file for the UW routing model
#' @param r raster with the same dimensions and coordinate system as the flow
#' direction file (could just use flow direction direction file)
#' @param gage_coords text file containing lon/lat data for gage coordinaates
#' @param basename base name for the save file
#' @param saveloc location to save the station location file
#' @details Creates a station location file and saves it to the specified directory. 
#' Adapted from an earlier Matlab version with the same name
#' @export
#' @examples 
#' r = raster("/Volumes/HD3/SWOTDA/Calibration/TuoSub/fdir.tif")
#' gage_coords = "/Volumes/HD3/SWOTDA/Calibration/TuoSub/gage_coords.txt"
#' saveloc <- "/Volumes/HD3/SWOTDA/Calibration/TuoSub"
#' basename <- "tuosub"
#' make_stnloc(r, gage_coords, basename, saveloc)

make_stnloc <- function(r, gage_coords, basename, saveloc)
{
  latlon <- read.table(gage_coords, col.names = c("lon", "lat"), sep=",")
  
  # Find row and column locations of gage_coords in r
  cols <- colFromX(r, latlon$lon)
  rows <- rowFromY(r, latlon$lat)
  
  rows = nrow(r) - round(rows) + 1;
  # rows = round(rows);
  
  cols = round(cols);
  
  # Write out the stnloc file
  savename <- file.path(saveloc, paste0(basename, "_stnloc.txt"))
  nstations = length(latlon$lon);
  for (k in 1:nstations)
  {
    pracma::fprintf("%d \t %s \t %d \t %d \t %d\n", 
                    1, 
                    paste0(substr(basename, 1, 3), k), 
                    cols[k], 
                    rows[k], 
                    -9999, 
                    file = savename, 
                    append = TRUE
                    )
    pracma::fprintf("%s\n", "NONE", file = savename, append = TRUE)
  }
  print(paste("Saved station location file as ", savename))
  
}

#' Make fract (trivial)
#' 
#' Makes a trivial fraction file for the UW routing model, where all the fraction values are 1
#' @param basinmask at the modeling resolution
#' @param saveloc save location
#' @param basename basename
#' @export

make_fract_trivial <- function(basinmask, saveloc, basename)
{
  fract1 <- basinmask
  fract1[is.na(fract1)] <- 0
  plot(fract1)
  writeRaster(fract1, file.path(saveloc, paste0(basename, "_fract_trivial")), 
              format = "ascii", overwrite = T)
  print(paste("Fraction file saved to", file.path(saveloc, paste0(basename, "_fract_trivial"))))
  return(fract1)
}
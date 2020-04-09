#' Make fraction file
#' 
#' Creates a fraction file for the UW routing model
#' An alternative function is provided for creating a fraction file of all 1's
#'
#' @param dem fine resolution DEM at least as large as the study area
#' @param bb basin boundary shapefile defining the study area
#' @param saveloc directory to save outputs
#' @param basename basename for the outputs
#' @param target_res target (coarse) resolution
#' @export
#' @return Returns the following rasters, in ASCII grid format:
#' Fraction File at VIC model resolution
#' Watershed mask at DEM resolution
#' Watershed mask at VIC model resolution
#' DEM clipped to basin boundaries
#' @details All coordinate systems should be geographic. I've been using WGS84.
#' If the pixels are not square (in geographic coordinates), this becomes a pain.
#' @examples
#' dem <- raster("/Volumes/HD4/SWOTDA/Data/IRB/irb_dem.tif")
#' bb <- readOGR("/Volumes/HD4/SWOTDA/Data/IRB/irb_bb_1_16/bb.shp")
#' saveloc <- "/Volumes/HD3/SWOTDA/FDT/v12032019/Routing_Inputs"
#' basename <- "irb"
#' fract <- make_fract(dem, bb, saveloc, basename, fineres = 1/16, coarseres = 1/16)

make_fract <- function(dem, bb, saveloc, basename, target_res = 1/16)
{
  
  print("Creating fraction file")
  
  fineres <- res(dem)[1] # DEM resolution (arc-seconds - "fine")
  coarseres <- target_res
  
  # Create a clipped DEM using the basin boundary shapefile
  dem.clipped <- clip_raster(dem, bb)
  writeRaster(dem.clipped, file.path(saveloc, paste0(basename, "_clipped_DEM_fine")),
              format = "ascii", overwrite = T)
  
  # Create a watershed mask at the DEM resolution
  dem.clipped.copy <- dem.clipped
  dem.clipped.copy[!is.na(dem.clipped.copy)] <- 1
  
  # Use a slightly larger extent than that of the clipped DEM
  e <- extent(dem.clipped)
  e@xmax <- ceiling(e@xmax)
  e@xmin <- floor(e@xmin)
  e@ymax <- ceiling(e@ymax)
  e@ymin <- floor(e@ymin)
  
  finemask <- raster(resolution = fineres, ext = e, crs = crs(bb))
  values(finemask) <- 0
  
  basinmask.fine <- raster::merge(dem.clipped.copy, finemask)
  
  writeRaster(basinmask.fine, file.path(saveloc, paste0(basename, "_basinmask_fine")), 
              format = "ascii", overwrite = T)
  
  
  # Make a grid with the VIC model resolution.
  # Extent must be larger than the real basin boundaries.
  
  # Must have a unique value in each cell (defines "zones").
  coarsegrid <- raster(resolution = coarseres , ext = e, crs = crs(bb))
  values(coarsegrid) <- 1:ncell(coarsegrid)
  
  # Resample the coarse grid to the DEM resolution, without changing the values
  finegrid <- resample(coarsegrid, basinmask.fine, method = "ngb")
  
  # Compute zonal mean
  zmean <- zonal(basinmask.fine, z = finegrid, progress = "text")
  
  # Output fraction file
  fract <- raster(coarsegrid)
  values(fract) <- zmean[,2]
  writeRaster(fract, file.path(saveloc, paste0(basename, "_fract")), 
              format = "ascii", overwrite = T)
  
  # Create a basin mask at the VIC model resolution
  basinmask.coarse <- fract
  basinmask.coarse[basinmask.coarse != 0] <- 1
  writeRaster(basinmask.coarse, file.path(saveloc, paste0(basename, "_basinmask_coarse")), 
              format = "ascii", overwrite = T)
  writeRaster(basinmask.coarse, file.path(saveloc, paste0(basename, "_basinmask_coarse")), 
                                          format = "GTiff", overwrite = T)
  
  print(paste("Outputs written to", savedir))
  return(fract)
}
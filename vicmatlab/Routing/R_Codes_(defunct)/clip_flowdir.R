 #' Clip FlowDir
 #' 
 #' Clips the flow direction file to the basinmask
 #' 
 #' @include raster
 #' @examples 
 flowdir <- "/Volumes/HD3/SWOTDA/FDT/v10282019/flow_direction_vic_d4_00.asc"
 basinmask <- "/Volumes/HD3/SWOTDA/FDT/v10282019/mangla_basinmask_coarse.tif"
 savename <- "/Volumes/HD3/SWOTDA/Calibration/Mangla/mangla_flowdir_vic"
 fd <- clip_flowdir(flowdir, basinmask, savename, input_format = "VIC", out_format = "ascii")

clip_flowdir <- function(flowdir, basinmask, savename, input_format = "GRASS", out_format = "ascii")
{
  
  fd <- raster(flowdir)
  
  # Convert format if necessary
  if (input_format == "GRASS")
  {
    fd_copy <- fd
    fd_copy[fd == 2] <- 1
    fd_copy[fd == 1] <- 2
    fd_copy[fd == 8] <- 3
    fd_copy[fd == 7] <- 4
    fd_copy[fd == 6] <- 5
    fd_copy[fd == 5] <- 6
    fd_copy[fd == 4] <- 7
    fd_copy[fd == 3] <- 8
    fd <- fd_copy
    print('Converted flow direction convention from GRASS to VIC')
  }
  
 mask1 <- raster(basinmask)
 fd_crop <- crop(fd, mask1)
 fd_clip <- overlay(fd_crop, mask1, fun = function(x,y) x*y, progress = "text")
 writeRaster(fd_clip, filename = savename, format = out_format, overwrite = TRUE)  
 
 return(fd_clip)
  
}


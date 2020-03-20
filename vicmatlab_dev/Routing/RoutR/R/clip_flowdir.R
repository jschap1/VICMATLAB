#' Clip FlowDir
#' 
#' Clips the flow direction file to the basinmask
#' @param fd flow direction map
#' @param basinmask mask to which to clip the fd file
#' @param savename output name to save the clipped flow direction file
#' @param out_format format to output the flow direction file. Can be GTiff or ascii.
#' @value Clipped flow direction file. Also saves the files as savename
#' @examples 
#' flowdir <- raster("/Volumes/HD3/SWOTDA/FDT/v10282019/flow_direction_vic_d4_00.asc")
#' basinmask <- raster("/Volumes/HD3/SWOTDA/FDT/v10282019/mangla_basinmask_coarse.tif")
#' savename <- "/Volumes/HD3/SWOTDA/Calibration/Mangla/mangla_flowdir_vic"
#' fd <- clip_flowdir(flowdir, basinmask, savename, out_format = "ascii")

clip_flowdir <- function(fd, basinmask, savename, out_format = "ascii")
{

 fd_crop <- crop(fd, basinmask)
 fd_clip <- overlay(fd_crop, basinmask, fun = function(x,y) x*y, progress = "text")
 writeRaster(fd_clip, filename = savename, format = out_format, overwrite = TRUE)  
 
 return(fd_clip)
  
}


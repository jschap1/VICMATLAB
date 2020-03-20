#' Convert flow direction convention
#' 
#' Converts flow direction convention among three different conventions
#' @param flowdir flow direction raster
#' @param conv1 name of input flow direction convention
#' @param conv2 name of desired output flow direction convention
#' @export

convert_flowdir <- function(flowdir, conv1, conv2)
{
  
  fd_copy <- flowdir
  
  if (conv1 == "VIC" & conv2 == "GRASS")
  {
    print("Converting from VIC to GRASS")
    fd_copy[flowdir == 1] <- 2
    fd_copy[flowdir == 2] <- 1
    fd_copy[flowdir == 3] <- 8
    fd_copy[flowdir == 4] <- 7
    fd_copy[flowdir == 5] <- 6
    fd_copy[flowdir == 6] <- 5
    fd_copy[flowdir == 7] <- 4
    fd_copy[flowdir == 8] <- 3    
  } else if (conv1 == "GRASS" & conv2 == "VIC")
  {
    print("Converting from GRASS to VIC")
    fd_copy(flowdir==2) = 1;
    fd_copy(flowdir==1) = 2;
    fd_copy(flowdir==8) = 3;
    fd_copy(flowdir==7) = 4;
    fd_copy(flowdir==6) = 5;
    fd_copy(flowdir==5) = 6;
    fd_copy(flowdir==4) = 7;
    fd_copy(flowdir==3) = 8;    
  } else if (conv1 == "ARCGIS" & conv2 == "GRASS")
  {
    print("Converting from ARCGIS to GRASS")
    fd_copy[flowdir==128] <- 1
    fd_copy[flowdir==64] <- 2
    fd_copy[flowdir==32] <- 3
    fd_copy[flowdir==16] <- 4
    fd_copy[flowdir==8] <- 5
    fd_copy[flowdir==4] <- 6
    fd_copy[flowdir==2] <- 7
    fd_copy[flowdir==1] <- 8    
  } else if (conv1 == "GRASS" & conv2 == "ARCGIS")
  {
    print("GRASS to ARCGIS not yet implemented")
  } else if (conv1 == "ARCGIS" & conv2 == "VIC")
  {
    print("Converting from ARCGIS to VIC")
    fd_copy(flowdir==64) = 1;
    fd_copy(flowdir==128) = 2;
    fd_copy(flowdir==1) = 3;
    fd_copy(flowdir==2) = 4;
    fd_copy(flowdir==4) = 5;
    fd_copy(flowdir==8) = 6;
    fd_copy(flowdir==16) = 7;
    fd_copy(flowdir==32) = 8;    
  } else if (conv1 == "VIC" & conv2 == "ARCGIS")
  {
    print("VIC to ARCGIS not yet implemented")
  } else
  {
    print("Please enter valid values for conv1 and conv2")
  }
  
  return(fd_copy)
  
}
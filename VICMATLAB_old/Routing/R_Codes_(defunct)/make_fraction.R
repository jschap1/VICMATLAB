# Make Fraction File
#
# INPUTS
# bb = basin boundary shapefile
# dem = fine resolution DEM
# resolution = target (coarse) resolution
# fraction = name of output fraction file
# finemask = name of output fine resolution mask file
# coarsemask = name of output coarse resolution mask file
# clipped_dem = name of clipped DEM to output
#
# OUTPUTS (in ASCII GRID format)
# Fraction File at VIC model resolution
# Watershed mask at DEM resolution
# Watershed mask at VIC model resolution
# DEM clipped to basin boundaries

# dem <- "/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_3/dem.tif"
# fract <- make_fraction_file(dem = dem, 
#                             bb = bb,
#                             resolution = 1/16,
#                             fraction = fractname,
#                             )

# Note: all coordinate systems should be geographic. I've been using WGS84.
# If the pixels are not square (in geographic coordinates), this becomes a pain.
#
# Requires raster, rgdal
#
# Example
# fract  <-  make_fraction_file(dem = "./FDT/merged_merit_30as.tif",
#                    bb = "./FDT/Rout/smaller/basin_1_4.shp",
#                    resolution = 1/4,
#                    fraction = "./Routing_Input_Prep/fract.asc",
#                    clipped_dem = "./Routing_Input_Prep/clipped_dem.tif",
#                    finemaskname = "./Routing_Input_Prep/finemask.tif",
#                    coarsemaskname = "./Routing_Input_Prep/coarsemask.tif")
# 
# bb.true <- readOGR("./Delineation/basin_3as.shp")
# par(mar = c(5,5,5,5))
# plot(fract, xlab = "Lon", ylab = "Lat", cex.lab = 2, cex.axis = 2)
# plot(bb.true, add = TRUE, lwd = 1)

make_fraction_file <- function(dem, bb, resolution, fraction, 
                               clipped_dem = "none", 
                               finemaskname = "none", 
                               coarsemaskname = "none")
{
  
  dem <- raster(dem)
  bb <- readOGR(bb)
  
  fineres <- res(dem)[1] # DEM resolution (arc-seconds - "fine")
  coarseres <- resolution # VIC model resolution (degrees - "coarse")
  
  # A user-written function from Carsten Neumann:
  # https://stat.ethz.ch/pipermail/r-sig-geo/2013-July/018912.html
  clip <- function(raster, shape) 
    {
      a1_crop <- crop(raster, shape)
      step1 <- rasterize(shape, a1_crop)
      return(a1_crop*step1)
    }
  
  # Create a clipped DEM using the basin boundary shapefile
  print("Clipping DEM to match basin boundary")
  dem.clipped <- clip(dem, bb)
  # dem.clipped <- dem
  
  if(clipped_dem!="none")
  {
    writeRaster(dem.clipped, clipped_dem, overwrite = T)
    print(paste("Saved clipped DEM:", clipped_dem))
  }

  # Create a watershed mask at the DEM resolution
  dem.clipped.copy <- dem.clipped
  dem.clipped.copy[!is.na(dem.clipped.copy)] <- 1
  
  # Make a grid with the VIC model resolution.
  # Extent must be larger than the real basin boundaries.
  
  # Use a slightly larger extent than that of the clipped DEM
  e <- extent(dem.clipped)
  e@xmax <- ceiling(e@xmax)
  e@xmin <- floor(e@xmin)
  e@ymax <- ceiling(e@ymax)
  e@ymin <- floor(e@ymin)
  
  finemask <- raster(resolution = fineres, ext = e, crs = crs(bb))
  values(finemask) <- 0
  basinmask.fine <- merge(dem.clipped.copy, finemask)
  
  if(finemaskname!="none")
  {
    writeRaster(basinmask.fine, finemaskname, overwrite = T)
    print(paste("Saved fine resolution basin mask:", finemaskname))
  }
  
  # Must have a unique value in each cell (defines "zones").
  coarsegrid <- raster(resolution = coarseres , ext = e, crs = crs(bb))
  values(coarsegrid) <- 1:ncell(coarsegrid)  
  
  # Resample the coarse grid to the DEM resolution, without changing the values
  finegrid <- resample(coarsegrid, basinmask.fine, method = "ngb")
  
  # Compute zonal mean
  print("Calculating fractions")
  zmean <- zonal(basinmask.fine, z = finegrid, progress = "text")
  
  # Output fraction file
  fract <- raster(coarsegrid)
  values(fract) <- zmean[,2]
  writeRaster(fract, fraction, format = "ascii", overwrite = T)
  print(paste("Saved fraction file:", fraction))
  
  # Create a basin mask at the VIC model resolution
  basinmask.coarse <- fract
  basinmask.coarse[basinmask.coarse != 0] <- 1
  
  if(coarsemaskname!="none")
  {
    writeRaster(basinmask.coarse, coarsemaskname, overwrite = T)
    print(paste("Saved coarse resolution basin mask:", coarsemaskname))
  }
 
  return(fract)
   
}

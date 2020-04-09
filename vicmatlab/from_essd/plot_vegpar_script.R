# Script for plotting UCRB vegetation parameters
#
# 2/18/2020 JRS

library(raster)
library(rgdal)
source("/Users/jschap/Documents/Research/VICGlobal/Codes/subset_geotiff.R")

basin_mask <- raster("/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif")

tifdir <- "/Volumes/HD3/VICParametersGlobal/VICGlobal/v1_4/Figures/tifs"
# lai.names <- list.files(file.path(tifdir), pattern = glob2rx("*lai*.tif"))

savedir <- "/Users/jschap/Documents/Research/VICGlobal/Data/uppcolo_vg_vegpars"

# LAI ---------------------------------------------------------------------------

r <- raster(file.path(tifdir, "jan_lai_nh.tif"))
january_lai <- subset_geotiff(r, basin_mask, savename = file.path(savedir, "jan_lai.tif"))

r <- raster(file.path(tifdir, "jul_lai_nh.tif"))
july_lai <- subset_geotiff(r, basin_mask, savename = file.path(savedir, "jul_lai.tif"))

opar <- par()
par(mfrow = c(1,2))
plot(january_lai, main = "January LAI")
plot(july_lai, main = "July LAI")

# FCAN ---------------------------------------------------------------------------

r <- raster(file.path(tifdir, "jan_fcan_nh.tif"))
january_fcan <- subset_geotiff(r, basin_mask, savename = file.path(savedir, "jan_fcan.tif"))

r <- raster(file.path(tifdir, "jul_fcan_nh.tif"))
july_fcan <- subset_geotiff(r, basin_mask, savename = file.path(savedir, "jul_fcan.tif"))

par(mfrow = c(1,2))
plot(january_fcan, main = "January FCAN")
plot(july_fcan, main = "July FCAN")

# VICGLobal vegetation cover -------------------------------------------------------------

tifdir <- "/Volumes/HD3/VICParametersGlobal/VICGlobal/vegetation/Vegetation_Fractions/Resampled_1_16"
vegnames <- list.files(tifdir, pattern = glob2rx("*.tif"), full.names = TRUE)

par(mfrow = c(3,6))
for (nam in vegnames)
{
  r <- raster(nam)
  # r_clip <- subset_geotiff(r, basin_mask, savename = file.path(savedir, basename(nam)))
  r_clip <- subset_geotiff(r, basin_mask)
  plot(r_clip, main = file_path_sans_ext(basename(nam)))
}

# Livneh vegetation cover and LAI -------------------------------------------------------------
# Note that the Livneh vegetation cover categories are somewhat different than VICGlobal categories
# Note that the Livneh LAI is defined on a gridcell-to-gridcell basis

vegnames <- list.files("/Volumes/HD3/VICParametersCONUS/Vegetation", pattern = glob2rx("*vegfract.tif"), full.names = TRUE)

par(mfrow = c(3,4))
for (nam in vegnames)
{
  r <- raster(nam)
  # r_clip <- subset_geotiff(r, basin_mask, savename = file.path(savedir, basename(nam)))
  r_clip <- subset_geotiff(r, basin_mask)
  plot(r_clip, main = file_path_sans_ext(basename(nam)))
}

# jan_alb <- raster("/Volumes/HD3/VICParametersCONUS/Vegetation/average_lai_alb/january_albedo.tif")
# jul_alb <- raster("/Volumes/HD3/VICParametersCONUS/Vegetation/average_lai_alb/july_albedo.tif")
jan_lai <- raster("/Volumes/HD3/VICParametersCONUS/Vegetation/average_lai_alb/january_LAI.tif")
jul_lai <- raster("/Volumes/HD3/VICParametersCONUS/Vegetation/average_lai_alb/july_LAI.tif")

savedir <- "/Users/jschap/Documents/Research/VICGlobal/Data/uppcolo_l15_vegpars"
# dir.create(savedir)
jan_lai_clip <- subset_geotiff(jan_lai, basin_mask, savename = file.path(savedir, "jan_lai.tif"))
jul_lai_clip <- subset_geotiff(jul_lai, basin_mask, savename = file.path(savedir, "jul_lai.tif"))

par(mfrow = c(1,2))
plot(jan_lai_clip, main = "January LAI")
plot(jul_lai_clip, main = "July LAI")

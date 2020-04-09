# Subsets the SWE variable from the Livneh NetCDF output files to the extent of a specified basin mask
#
# Written 1/31/2020 JRS
#
# Could be generalized to any VIC image mode output variable

library(raster)
library(lubridate)
library(ncdf4)
source("/Users/jschap/Documents/Research/VICGlobal/Codes/subset_livneh_func.R")

r <- raster('/Volumes/HD3/SWOTDA/temp1.tif')
plot(r)

livneh_dir <- "/Volumes/HD3/Livneh_2013/Fluxes"
basin_mask <- "/Volumes/HD4/SWOTDA/Data/Colorado/colo_mask.tif"


mask1 <- raster(basin_mask)
plot(mask1)

start_date <- ymd("2000-10-01")
end_date <- ymd("2011-09-30")

# Subset SWE ---------------------------------------------------------------------

outdir <- "/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/SWE"
dir.create(outdir)
swe_crop <- subset_livneh_func(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name = "SWE")
 
tifnames <- list.files(outdir, "*.tif", full.names = TRUE)
swe_all <- stack(tifnames)
average_swe <- calc(swe_all, fun = mean, progress = "text")

# Subset ET -------------------------------------------------------------

outdir <- "/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/ET"
dir.create(outdir)
evap_crop <- subset_livneh_func(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name = "ET")

# # Subset potential ET ------------------------------------------------------
# outdir <- "/Users/jschap/Documents/Research/VICGlobal/Data/uppcololivneh_pet_long/"
# dir.create(outdir)
# pet_crop <- subset_livneh_func(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name = "PETLong")

# Subset soil moisture -------------------------------------------------------------

outdir <- "/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/SM1/"
dir.create(outdir)
sm_level <- 1
sm_crop <- subset_livneh_soil_moisture(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name = "SoilMoist", sm_level)

outdir <- "/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/SM2/"
dir.create(outdir)
sm_level <- 2
sm_crop <- subset_livneh_soil_moisture(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name = "SoilMoist", sm_level)

outdir <- "/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/SM3/"
dir.create(outdir)
sm_level <- 3
sm_crop <- subset_livneh_soil_moisture(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name = "SoilMoist", sm_level)


# Subset runoff  ------------------------------------------------------
outdir <- "/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/runoff/"
dir.create(outdir)
runoff_crop <- subset_livneh_func(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name = "Q")

# Subset baseflow  ------------------------------------------------------

outdir <- "/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/baseflow/"
dir.create(outdir)
baseflow_crop <- subset_livneh_func(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name = "Qsb")

# Subset canopy moisture  ------------------------------------------------------

outdir <- "/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/wdew/"
dir.create(outdir)
wdew_crop <- subset_livneh_func(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name = "WDew")

# Subset precipitation  ------------------------------------------------------

outdir <- "/Users/jschap/Documents/Research/VICGlobal/Data/UCRB_L13/PREC/"
dir.create(outdir)
precip_crop <- subset_livneh_func(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name = "Prec")

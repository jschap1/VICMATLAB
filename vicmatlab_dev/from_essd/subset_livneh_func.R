#' Subset Livneh Function
#'
#' Subsets the Livneh et al. (2013) VIC flux outputs to a particular basin and time period.
#' Was written with SWE in mind, so some of the variable names still reflect this.
#' @param variable_name SWE, for example
#'

subset_livneh_func <- function(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name)
{
  
  first_iter <- 1
  current_date <- start_date
  ndays <- end_date - start_date + 1
  
  while (current_date < end_date)
  {
    
    current_year <- year(current_date)
    current_month <- month(current_date)
    if (current_month < 10)
    {
      month1 <- paste0("0", current_month)
    } else
    {
      month1 <- current_month
    }
    
    fname <- paste0("VIC_fluxes_Livneh_CONUSExt_v.1.2_2013.", current_year, month1, ".nc")
    # fname <- paste0("Fluxes_Livneh_NAmerExt_15Oct2014.", current_year, month1, ".nc")
    
    filename <- file.path(livneh_dir, fname)
    nc_obj <- nc_open(filename)
    
    if (first_iter == 1)
    {
      lon <- ncvar_get(nc_obj, "lon")
      lat <- ncvar_get(nc_obj, "lat")
      swe <- ncvar_get(nc_obj, variable_name)
      nx <- dim(swe)[1]
      ny <- dim(swe)[2]
      r1 <- raster(filename) # getting a template
      first_iter <- 0
    }
    
    swe_raster <- raster(r1)
    swe <- ncvar_get(nc_obj, variable_name)
    nt <- dim(swe)[3]
    
    for (t in 1:nt) # create a SWE stack for the current month
    {
      swe_vect <- as.vector(swe[,,t]) # for one time step
      if (t > 1)
      {
        swe_raster <- addLayer(swe_raster, swe_raster[[t-1]])  
        values(swe_raster[[t]]) <- swe_vect
      } else
      {
        values(swe_raster) <- swe_vect  
      }
      # if (t%%10==0)
      # {
      #   print(t)
      # }
    }
    
    # Probably quickest/simplest to work in monthly averages for now
    mean_swe <- mean(swe_raster)
    mean_swe <- flip(mean_swe, direction = "y")
    swe_crop <- crop(mean_swe, mask1)
    
    fname_out <- paste0("monthly_", variable_name, "_UCRB_", current_date, ".tif")
    writeRaster(swe_crop, file = file.path(outdir, fname_out), overwrite = TRUE)
    print(paste("Saved cropped, monthly", variable_name, "data:", fname_out))
    
    current_date <- current_date %m+% months(1) # advance to next month
    
  }
  
  return(swe_crop)
  
}

#' Subset Livneh Soil Moisture
#' 
#' If the variable is soil moisture, then a separate function is needed because soil moisture has multiple levels, usually three.

subset_livneh_soil_moisture <- function(livneh_dir, basin_mask, start_date, end_date, outdir, variable_name = "SoilMoist", sm_level)
{
  
  first_iter <- 1
  current_date <- start_date
  ndays <- end_date - start_date + 1
  
  while (current_date < end_date)
  {
    
    current_year <- year(current_date)
    current_month <- month(current_date)
    if (current_month < 10)
    {
      month1 <- paste0("0", current_month)
    } else
    {
      month1 <- current_month
    }
    
    fname <- paste0("VIC_fluxes_Livneh_CONUSExt_v.1.2_2013.", current_year, month1, ".nc")
    # fname <- paste0("Fluxes_Livneh_NAmerExt_15Oct2014.", current_year, month1, ".nc")    
    filename <- file.path(livneh_dir, fname)
    nc_obj <- nc_open(filename)
    
    if (first_iter == 1)
    {
      lon <- ncvar_get(nc_obj, "lon")
      lat <- ncvar_get(nc_obj, "lat")
      soil_moist <- ncvar_get(nc_obj, variable_name)
      nx <- dim(soil_moist)[1]
      ny <- dim(soil_moist)[2]
      r1 <- raster(filename) # getting a template
      first_iter <- 0
    }
    
    sm_raster <- raster(r1)
    soil_moist <- ncvar_get(nc_obj, variable_name)
    nt <- dim(soil_moist)[4]
    
    for (t in 1:nt) # create a SWE stack for the current month
    {
      sm_vect <- as.vector(soil_moist[,,sm_level,t]) # for one time step
      if (t > 1)
      {
        sm_raster <- addLayer(sm_raster, sm_raster[[t-1]])  
        values(sm_raster[[t]]) <- sm_vect
      } else
      {
        values(sm_raster) <- sm_vect  
      }
      # if (t%%10==0)
      # {
      #   print(t)
      # }
    }
    
    # Probably quickest/simplest to work in monthly averages for now
    mean_sm <- mean(sm_raster)
    mean_sm <- flip(mean_sm, direction = "y")
    soil_moisture_crop <- crop(mean_sm, mask1)
    
    fname_out <- paste0("monthly_", variable_name, "_UCRB_", current_date, ".tif")
    writeRaster(soil_moisture_crop, file = file.path(outdir, fname_out), overwrite = TRUE)
    print(paste("Saved cropped, monthly", variable_name, "data:", fname_out))
    
    current_date <- current_date %m+% months(1) # advance to next month
    
  }
  
  return(swe_crop)
  
}

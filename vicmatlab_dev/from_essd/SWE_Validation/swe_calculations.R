library(raster)

saveloc <- "/Volumes/HD4/SWOTDA/Data/Tuolumne/v1_6/SWE_Validation/overlapping_times"

# swe_map_vg <- readRDS(file.path(saveloc, "mean_swe_map_vg.rds"))
# swe_map_livneh <- readRDS(file.path(saveloc, "mean_swe_map_livneh.rds"))
swe_ts_vg <- readRDS(file.path(saveloc, "mean_swe_ts_vg.rds"))
swe_ts_livneh <- readRDS(file.path(saveloc, "mean_swe_ts_livneh.rds"))

# ------------------------------------------------------------------------------------

swe_map_vg <- raster(file.path(saveloc, "swe_vicglobal.tif"))
swe_map_livneh <- raster(file.path(saveloc, "swe_livneh.tif"))

mean(values(swe_map_livneh), na.rm = TRUE)
mean(values(swe_map_vg), na.rm = TRUE)

swe_all_livneh <- readRDS(file.path(saveloc, "swe_data_livneh_all.rds"))

# I have everything except the VG time series
# Create swe_ts_vg time series for the overlapping and subset it to the overlapping dates

s_brick
s_brick_copy
s_brick_sub

swe_ts_vg <- vector(length = nlayers(s_brick_sub))
for (k in 1:length(swe_ts_vg))
{
  swe_ts_vg[k] <- mean(values(s_brick_sub[[k]]), na.rm=T)
  if (k%%1000==0)
  {
    print(k)
  }
}

saveRDS(swe_ts_vg, file.path(saveloc, "mean_swe_ts_vg.rds"))
write.table(swe_ts_vg, file.path(saveloc, "swe_ts_vicglobal.txt"), row.names = FALSE)

# Experiment

sample_swe <- brick(nrows = 2, ncols = 2, nl = 2)
values(sample_swe[[1]]) <- c(4,3,2,NA)
values(sample_swe[[2]]) <- c(0,7,0,NA)

sample_swe_ts <- 1
sample_swe_map <- 1

# Try calculating average time series
sample_swe_ts <- vector(length = nlayers(sample_swe))
for (k in 1:nlayers(sample_swe))
{
  sample_swe_ts[k] <- mean(values(sample_swe[[k]]), na.rm=T)
}
# Average value as expected (16/6)

# Try calculating average map
sample_swe_map <- mean(sample_swe)
mean(values(sample_swe_map), na.rm = TRUE)


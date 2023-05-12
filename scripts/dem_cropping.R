pacman::p_load("terra", "raster", "sf", "sp", "tidyverse")

java_crop <- vect("large/large_file_crops/java_crop.kml")

files <- list.files(path = "large/srtm_30/", full.names = T)

vert_file <- vrt(c(files), "vert_dem.tif", overwrite = T)

vert_agr <- terra::aggregate(vert_file, fact = 10, fun = 'mean')

dem_crop <- crop(vert_agr, java_crop)

plot(dem_crop)

writeRaster(dem_crop, "large/java_300_dem.tif")
writeRaster(vert_agr, "large/indo_dem_resample_300.tif")



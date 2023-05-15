pacman::p_load("terra", "raster", "sf", "sp", "tidyverse")

## Resamplin srtm 30 m to match Jung habitat map ----
# loading in additional files for crop and resample
java_crop <- vect("large/large_file_crops/java_crop.kml")
java_hab <- rast("large/jung_hab_java.tif")
esa_java <- rast("large/esa_java.tif")

# loading DEM files and creating vrt raster to merge
files <- list.files(path = "large/srtm_30/", full.names = T)
vert_file <- vrt(c(files), "vert_dem.tif", overwrite = T)

# resampling DEM raster to match jung hab
vert_resamp <- resample(rast("large/vert_dem.tif"), java_hab, method = 'bilinear') # having to load in the vert file as a rast as it didnt like if using the vert file created earlier

# masking dem
dem_masked <- mask(vert_resamp, java_crop)

# comparing rasters 
compareGeom(java_hab, dem_masked)

# writing file
writeRaster(vert_resamp, "large/indo_dem_resample_300.tif", overwrite = T)
writeRaster(dem_masked, "large/java_jung_dem.tif", overwrite = T)

## Resampling Jung hab map to 1 km ----
java_hab_1km <- resample(java_hab, esa_java, method = 'near')

# comparing java hab and esa_java
compareGeom(java_hab_1km, esa_java)

# writing file 
writeRaster(java_hab_1km, "large/java_jung_hab_1km.tiff", overwrite = T)

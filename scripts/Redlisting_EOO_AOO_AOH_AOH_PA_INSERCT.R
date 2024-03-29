##### EOO, AOO, AOH, PA intersect script with smoothr. #####
## Seth Ratcliffe - s.ratcliffe@kew.org
## 20/06/2023

### Packages, libraries, functions and occurrence upload ----
## 1. Packages and Libraries
install.packages(c("sf", "sp", "raster", "rCAT", "pacman", "tidyverse", 
                   "stars", "terra", "smoothr")) 

remotes::install_github('rstudio/leaflet') # this script requires the development version of to create the AOH maps

pacman::p_load(sf, leaflet, raster, rCAT, tidyverse, stars, terra, smoothr, sp) 

## 2. Functions 
source("scripts/seth_functions.R")

## 3. Data Import and cleaning
occs_raw <- read.csv("IUCN_point_files/senecio_fresenii_iucn_pointfile.csv")
occ_points <-  occs_raw %>% 
    filter.occurences(T) 

### EOO, AOO and Mapping ----
## 4. Calculate EOO and AOO
eoo_aoo <- cal.eoo.aoo(occ_points)

## 5. Produce Map 
make.map(occ_points)

### AOH -----
## 2. Making boundary - if buffer wants to be added define the distance as b where 0.1 = 10km 
boundary <- make.boundary(occ_points, eoo = T, buffer = F) 

## 3. AOH Parameters
# setting min and max elevation
elevmin <- 1500 
elevmax <- 4300

# elevation raster
DEMrast <- rast("large/dem.tif")

# habitat raster
habstack <- rast("large/jung_hab_1km_global.tiff")

# Defining habitat codes
ESA_codes <- data.frame(ESA_codes = c(105, 307, 407))  

## 4. Generate the AOH
theAOH <- calc.aoh.sing(DEMrast, habstack, ESA_codes, boundary, elevmin, elevmax)

## 6. Generate the Red List stats from AOH 
cal.aoh.stats(theAOH)

## 7. Smoothing AoH
# making aoh a polygon 
aoh_polygon <- theAOH %>%
    terra::as.polygons() %>%
    st_as_sf() %>%
    st_union() %>%
    st_make_valid()

# dropping crumbs 
aoh_no_crumbs <- drop_crumbs(aoh_polygon, units::set_units(3, km^2))

# smoothing aoh
aoh_smooth <- smooth(aoh_no_crumbs, method = 'ksmooth', smoothness = 3) %>%
    st_union() %>%
    st_cast('POLYGON')

## 8. Map View 
make.aoh.map(occ_points, aoh_smooth, boundary_aoh = T, aoh_raster = F)

### Protected Area Occurrence and Percentage of AoH covered by protected areas ----
## 9. loading in WDPA data and wrangling AoH data
# setting path
files_path <- list.files(path = "large/wdpa_data/", pattern = "*shp", full.names = T )

# reading in and combing into single spatvect
wdpa_comb <- lapply(files_path, read_sf) %>%
    do.call(rbind, .) %>%
    vect()

# Occurrence in PA's 
pa.occurrence(wdpa_comb, occ_points)

# Converting AoH to polygon 
aoh_polygon <- as.polygons(theAOH) %>%
    project(., crs(wdpa_comb))

## 10. Masking wpda by aoh and calculating area 
wpda_masked <- terra::intersect(aoh_polygon, wdpa_comb)
area_of_aoh_in_pa <- print(sum(expanse(wpda_masked))/sum(expanse(aoh_polygon))*100)

### Data Export ----
## 11. Export raw AOH 
# exporting aoh raw raster as shp file
writeVector(as.polygons(theAOH), "aoh_outs/senecio_fresenii_raw_aoh.shp",
            filetype = 'ESRI Shapefile', overwrite = T) 

## 12. Exporting smoothed AOH with SIS datatable
# adding required dataframe for shp. file 
sis_dataframe <- as.data.frame(occs_raw[1,]) %>%
    select(-dec_lat, -dec_long, -spatialref, -event_year, -basisofrec, -catalog_no, 
           -recordedby, -recordno)
sis_dataframe$source <- NA

aoh_with_sis <- sp::merge(aoh_smooth, sis_dataframe)

# writing aoh.shp file  
st_write(aoh_with_sis, "aoh_outs/senecio_fresenii_smoothed_aoh.shp", overwrite = T)

## 13. Exporting for external data tools
# exporting boundary as shape file
st_write(boundary, "aoh_outs/boundaries/.kml", 
         driver = 'kml')

# exporting occurrences as shape file with boundary
points_spat <- occ_points %>% 
    select(long, lat) %>% 
    st_as_sf(coords = c("long", "lat"), crs = "4236") %>%
    st_buffer(dist = 0.01) #1 km buffer

st_write(points_spat, "aoh_outs/boundaries/.kml",
         driver = 'kml')

# exporting AOH PA intersection 
writeVector(wpda_masked, "aoh_outs/pa_intersects/echinops_ellenbeckii_pa_intersect.shp",
            filetype = 'ESRI Shapefile', overwrite = T)
    

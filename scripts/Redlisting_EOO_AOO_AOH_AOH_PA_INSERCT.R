##### EOO, AOO, AOH, PA intersect script with smoothr. #####
## Seth Ratcliffe - s.ratcliffe@kew.org
## 20/06/2023

### Packages, libraries, functions and occurrence upload ----
## 1. Packages and Libraries
install.packages(c("sf", "sp", "leaflet", "raster", "rCAT", "pacman", "tidyverse", 
                   "stars", "terra", "smoothr")) # requires development version of leaflet for mapping spatrast

pacman::p_load(sf, leaflet, raster, rCAT, tidyverse, stars, terra, smoothr, sp) 

## 2. Functions 
source("scripts/seth_functions.R")

## 3. Data Import and cleaning
occs_raw <- read.csv("IUCN_point_files/helichrysum_harennense_iucn_pointfile.csv")
occ_points <-  occs_raw %>% 
    filter.occurences(T) 

### EOO, AOO and Mapping ----
## 4. Calculate EOO and AOO
eoo_aoo <- cal.eoo.aoo(occ_points)

## 5. Produce Map 
make.map(occ_points)

### AOH -----
## 2. Making boundary - if buffer wants to be added define the distance as b where 0.1 = 10km 
boundary <- make.boundary(occ_points, 0.01, eoo = F, buffer = T) 

## 3. AOH Parameters
# setting min and max elevation
elevmin <- 3200  
elevmax <- 3500 

# elevation raster
DEMrast <- rast("large/eth_DEM_100.tif")

# habitat raster
habstack <- rast("large/eth_jung.tif")

# Defining habitat codes
ESA_codes <- data.frame(ESA_codes = 307)  

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
make.aoh.map(occ_points, theAOH, boundary_aoh = F, aoh_raster = T)

### Percentage of AoH covered by protected areas ----
## 9. loading in WDPA data and wrangling AoH data
# setting path
files_path <- list.files(path = "large/wdpa_data/", pattern = "*shp", full.names = T )

# reading in and combing into single spatvect
wdpa_comb <- lapply(files_path, read_sf) %>%
    do.call(rbind, .) %>%
    vect()

# Converting raw AoH to polygon 
aoh_polygon <- as.polygons(terra::rast(theAOH)) %>%
    project(., crs(wdpa_comb))

## 10. Masking wpda by aoh and calculating area 
wpda_masked <- terra::intersect(aoh_polygon, wdpa_comb)
area_of_aoh_in_pa <- print(sum(expanse(wpda_masked))/sum(expanse(aoh_polygon))*100)

### Data Export ----
## 11. Export raw AOH 
# converting and exporting aoh raw raster to sf
aoh_sf <- st_as_stars(theAOH) %>% # converting to stars object for sf transformation
    st_as_sf(as_points = F, merge = T) # converting to sf and merging points

st_write(aoh_sf, "aoh_outs/.shp")

## 12. Exporting smoothed AOH with SIS datatable
# adding required dataframe for shp. file 
sis_dataframe <- as.data.frame(occs_raw[1,]) %>%
    select(-dec_lat, -dec_long, -spatialref, -event_year, -basisofrec, -catalog_no, 
           -recordedby, -recordno)
sis_dataframe$source <- NA

aoh_with_sis <- sp::merge(aoh_smooth, sis_dataframe)

# writing aoh.shp file  
st_write(aoh_with_sis, "aoh_outs/.shp", overwrite = T)

## 13. Exporting for external data tools
# exporting boundary as shape file
st_write(boundary, "aoh_outs/boundaries/.kml", 
         driver = 'kml')

# exporting occurences as shape file with boundary
points_spat <- occ_points %>% 
    select(long, lat) %>% 
    st_as_sf(coords = c("long", "lat"), crs = "4236") %>%
    st_buffer(dist = 0.01) #1 km buffer

st_write(points_spat, "aoh_outs/boundaries/.kml",
         driver = 'kml')

# exporting AOH PA intersection 
writeVector(wpda_masked, "aoh_outs/pa_intersects/.shp",
            filetype = 'ESRI Shapefile', overwrite = T)

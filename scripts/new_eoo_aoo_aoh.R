##### EOO, AOO, AOH, PA intersect script with smoothr. #####
## Seth Ratcliffe - s.ratcliffe@kew.org
## 20/06/2023

### Packages, libraries, functions and occurrence upload ----
## 1. Packages and Libraries
install.packages(c("sf", "sp", "leaflet", "raster", "rCAT", "pacman", "tidyverse", 
                   "stars", "terra", "smoothr")) 

pacman::p_load(sf, leaflet, raster, rCAT, tidyverse, stars, terra, smoothr, sp) 

## 2. Functions 
source("scripts/functions.R")

## 3. Data Import and cleaning
occs_raw <- read.csv("IUCN_point_files/Onobrychis_richardi_IUCN_pointfile.csv")
occ_points <-  occs_raw %>% 
    filter.occurences(T) 

### EOO, AOO and Mapping ----
## 4. Calculate EOO and AOO
eoo_aoo <- cal.eoo.aoo(occ_points)

## 5. Produce Map 
make.map(occ_points)

### AOH -----
## 2. Data Import
# Google earth KMLimport
est_range <- st_read("", type = 3) 

# Multi polygon second kml import
poly_1 <- st_read('', type = 3)

# joining the two kmls
poly_union <- st_union(est_range, poly_1)

## 3. Making boundary
boundary <- make.boundary(occ_points, 0.1, eoo = T, buffer = T) # if buffer wants to be added define the distance as b where 0.1 = 1km

## 4. Parameters
# setting min and max elevation
elevmin <- 2100 # Varies dependent on species 
elevmax <- 3500

# creating mask 
mask <- boundary # This can either be the est_range KML or the boundary object from the EOO calculation                             

# elevation raster
DEMrast <- raster::raster("large/dem.tif") # elevation data

# habitat raster
habstack <- raster::raster("large/ethio_jung_hab_1km.tiff")

# Defining habitat codes
ESA_codes <- data.frame(ESA_codes = c(105, 305, 307, 405, 407)) # This will vary dependent on the habitat type 

## 4. Generate the AOH
theDEM <- dem(DEMrast, mask, elevmin, elevmax)
theHAB <- habitat(habstack, mask,  ESA_codes)
theAOH <- aoh(theHAB, theDEM, mask)

## 5. Generate the Red List stats from AOH 
cal.aoh.stats(theAOH)

## 6. Smoothing AoH
# making aoh a polygon 
aoh_polygon <- theAOH %>%
    rasterToPolygons() %>%
    st_as_sf() %>%
    st_union() %>%
    st_make_valid()

# dropping crumbs 
aoh_no_crumbs <- drop_crumbs(aoh_polygon, units::set_units(3, km^2))

# smoothing aoh
aoh_smooth <- smooth(aoh_no_crumbs, method = 'ksmooth', smoothness = 3) %>%
    st_union() %>%
    st_cast('POLYGON')

## 7. Map View 
make.aoh.map(occ_points, aoh_smooth, boundary_aoh = F, aoh_raster = F)
cal.aoh.stats(theAOH)

### Percentage of AoH covered by protected areas ----
## 8. loading in WDPA data and wrangling AoH data
# setting path
files_path <- list.files(path = "large/wdpa_data/", pattern = "*shp", full.names = T )

# reading in and combing into single spatvect
wdpa_comb <- lapply(files_path, read_sf) %>%
    do.call(rbind, .) %>%
    vect()

# Converting raw AoH to polygon 
aoh_polygon <- as.polygons(terra::rast(theAOH)) %>%
    project(., crs(wdpa_comb))

## 9. Masking wpda by aoh and calculating area 
wpda_masked <- terra::intersect(aoh_polygon, wdpa_comb)
area_of_aoh_in_pa <- print(sum(expanse(wpda_masked))/sum(expanse(aoh_polygon))*100)

### Data Export ----
## 10. Export raw AOH 
# converting and exporting aoh raw raster to sf
aoh_sf <- st_as_stars(theAOH) %>% # converting to stars object for sf transformation
    st_as_sf(as_points = F, merge = T) # converting to sf and merging points

st_write(aoh_sf, "aoh_outs/argyrolobium_schimperianum_aoh_raw.shp")

## 11. Exporting smoothed AOH with SIS datatable
# adding required dataframe for shp. file 
sis_dataframe <- as.data.frame(occs_raw[1,]) %>%
    select(-dec_lat, -dec_long, -spatialref, -event_year, -basisofrec, -catalog_no, 
           -recordedby, -recordno)
sis_dataframe$source <- NA

aoh_with_sis <- sp::merge(aoh_smooth, sis_dataframe)

# writing aoh.shp file  
st_write(aoh_with_sis, "aoh_outs/argyrolobium_schimperianum_distribution_polygon.shp", overwrite = T)

## 12. Exporting for external data tools
# exporting boundary as shape file
st_write(boundary, "aoh_outs/boundaries/Rhipidoglossum_candidum.kml", 
         driver = 'kml')

# exporting occurences as shape file with boundary
points_spat <- occ_points %>% 
    select(long, lat) %>% 
    st_as_sf(coords = c("long", "lat"), crs = "4236") %>%
    st_buffer(dist = 0.01) #1 km buffer

st_write(points_spat, "aoh_outs/boundaries/vaccinium_cuneifoliums_points_1km.kml",
         driver = 'kml')

# exporting AOH PA intersection 
writeVector(wpda_masked, "aoh_outs/pa_intersects/argyrolobium_schimperianum_pa_intersect.shp",
            filetype = 'ESRI Shapefile', overwrite = T)

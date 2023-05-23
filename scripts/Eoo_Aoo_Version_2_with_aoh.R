##### EOO and AOO for Use by the CAA Assessment Team. #####
## Seth Ratcliffe - s.ratcliffe@kew.org
## 25/11/2020


### EOO, AOO and Mapping ----
## 1. Packages and Libraries
install.packages(c("sf", "leaflet", "raster", "rCAT", "pacman", "tidyverse", 
                   "stars", "terra")) # here we are installing the required packages. 

# loading packages
pacman::p_load(sf, leaflet, raster, rCAT, tidyverse, stars, terra) # here we are loading the packages we need for this session. 

## 2. Functions 
source("scripts/functions.R")

## 3. Data Import and cleaning
occ_points <- read.csv("IUCN_point_files/Shirakiopsis_virgata_IUCN_Pointfile.csv") %>% 
  filter.occurences(T) # If you want to use all occurrences change the T to an F

## 4. Calculate EOO and AOO
eoo_aoo <- cal.eoo.aoo(occ_points)

## 5. Produce Map 
make.map(occ_points)

### AOH -----
## 2. Data Import
# use this if loading in from a Google earth KML
est_range <- st_read("polys/festuca_simien.kml", type = 3) 

# if you want to use a multi polygon the individual components have to be loaded 
# in separately then unioned 
poly_1 <- st_read('polys/festuca_arsi.kml', type = 3)

poly_union <- st_union(est_range, poly_1)

# use this if using the spp EOO as the boundary
boundary <- make.boundary(occ_points)

# add buffer to EOO mcp if wanted
buffer <- st_buffer(boundary, dist = 0.1)
boundary_buffer <- st_union(buffer, boundary) # joining mcp with buffer
    
## 3. Parameters
# setting min and max elevation
elevmin <- 0 # Varies dependent on species 
elevmax <- 500

# creating mask 
mask <- boundary_buffer # This can either be the est_range KML or the boundary object from the EOO calculation                             

# elevation raster
DEMrast <- raster::raster("large/dem.tif") # elevation data

# habitat raster
habstack <- raster::raster("large/java_jung_hab_1km.tiff")

# Defining habitat codes
ESA_codes <- data.frame(ESA_codes = c(106, 107, 108)) # This will vary dependent on the habitat type 

## 4. Generate the AOH
theDEM <- dem(DEMrast, mask, elevmin, elevmax)
theHAB <- habitat(habstack, mask,  ESA_codes)
theAOH <- aoh(theHAB, theDEM, mask)

## 5. Generate the Red List stats from AOH 
cal.aoh.stats(theAOH)

## 6. Map View 
make.aoh.map(occ_points, theAOH, F)

### Percentage of AoH covered by protected areas ----
## 7. loading in WDPA data and wrangling AoH data
# setting path
files_path <- list.files(path = "large/wpda_data/indo/", pattern = "*shp", full.names = T )

# reading in and combing into single spatvect
wdpa_comb <- lapply(files_path, read_sf) %>%
    do.call(rbind, .) %>%
    vect()

# Converting AoH to polygon 
aoh_polygon <- as.polygons(terra::rast(theAOH)) %>%
    project(., crs(wdpa_comb))

## 8. Masking wpda by aoh and calculating area 
wpda_masked <- terra::intersect(aoh_polygon, wdpa_comb)
area_of_aoh_in_pa <- print(sum(expanse(wpda_masked))/sum(expanse(aoh_polygon))*100)

### Population estimate ----
## 10. estimating population density from AoH
# setting mean population density (from external data)
mean_pop_dens_km <- 679.4

# calculating aoh area
aoh_area <- terra::expanse(aoh_terra, unit = 'km', transform = F)

# estimating population 
mean_pop_dens_km * aoh_area

### Data Export ----
## 11. Export AOH .shp for upload to SIS
# converting aoh raster to sf
aoh_sf <- st_as_stars(theAOH) %>% # converting to stars object for sf transformation
    st_as_sf(as_points = F, merge = T) # converting to sf and merging points

# exporting boundary as shape file
st_write(boundary, "aoh_outs/boundaries/vaccinium_cuneifolium.kml", 
         driver = 'kml')

# exporting occurences as shape file with boundary
points_spat <- occ_points %>% 
    select(long, lat) %>% 
    st_as_sf(coords = c("long", "lat"), crs = "4236") %>%
    st_buffer(dist = 0.01) #1 km buffer

st_write(points_spat, "aoh_outs/boundaries/vaccinium_cuneifoliums_points_1km.kml",
         driver = 'kml')

# writing .shp file  
st_write(aoh_sf, "aoh_outs/alangium_villosum_esa_1km.shp")


# exporting AOH PA intersection 
writeVector(wpda_mask, "aoh_outs/vaccinium_cuneifolium_pa_intersect.kml",
            filetype = 'kml', overwrite = T)

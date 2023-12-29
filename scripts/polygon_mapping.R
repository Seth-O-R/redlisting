## Polygon mapping ##
## Seth Ratcliffe - s.ratcliffe@kew.org
## 07/11/2023

### Packages, libraries, functions and occurrence upload ----
## 1. Packages and Libraries
install.packages(c("sf", "sp", "pacman", "tidyverse", 
                   "stars", "terra", "rCAT")) 

remotes::install_github('rstudio/leaflet') 

pacman::p_load(sf, tidyverse, stars, terra, sp, leaflet, rCAT) 

## 2. Functions 
source("scripts/seth_functions.R")

## 3. Data Import and cleaning
# occurrence points
occs_raw <- read.csv("IUCN_point_files/Lasjia_griseifolia_IUCN_pointfile.csv")
occ_points <-  occs_raw %>% 
    filter.occurences(T) 

# Polygon
polygon_dist <- vect("large/Re_ Lasjia Full Assessment/SUL_ultramafic_above_2000_singlepart.shp", crs = '+proj=longlat')

## 4. Intersect Points and polygons 
# vectorising occs for intersect 
occ_vect <- vect(occ_points, geom = c('long', 'lat'), crs = '+proj=longlat')

# intersect
poly_intersect <- polygon_dist[occ_vect]

## 5. Merging SIS data with occupied polygons 
# extracting required data for SIS
sis_dataframe <- as.data.frame(occs_raw[1,]) %>%
    select(-dec_lat, -dec_long, -spatialref, -event_year, -basisofrec, -catalog_no, 
           -recordedby, -recordno)
sis_dataframe$source <- NA

# converting vect to sf selecting required data and merging
poly_with_sis <- st_as_sf(poly_intersect) %>%
    select(ID, geometry) %>%
    sp::merge(sis_dataframe)

## 6. Adding Possibly Extant to SIS polygpon map
# selecting only polygons where the species is not recorded
poly_no_intersect <- erase(polygon_dist, occ_vect) 

# Altering presence coding and merging with polygns
sis_df_poss_extant <- mutate(sis_dataframe, 'presence' = 3)

poss_poly_with_sis <- st_as_sf(poly_no_intersect) %>%
    select(ID, geometry) %>%
    sp::merge(sis_df_poss_extant)

## 7. Merging all polygons 
polys_all_coding <- rbind(poly_with_sis, poss_poly_with_sis)

## 8. Calculating EOO and AOO 
point_stats <- cal.eoo.aoo(occ_points)
polygon_extant_stats <- polygon.eoo.aoo(poly_with_sis)
polygon_all_stats <- polygon.eoo.aoo(polys_all_coding)

## 9. Map
make.poly.map(occ_points, polys_all_coding)

## 10. exporting polygons 
st_write(polys_all_coding, "aoh_outs/polygon_maps/lasjia_griseifolia_polygon_distribution.shp",
         overwrite = T)

st_write(poly_with_sis, "aoh_outs/boundaries/lasjia_greseifolia_extant_poly.kml", 
         driver = 'kml')

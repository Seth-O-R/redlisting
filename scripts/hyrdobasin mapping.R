## Hydrobasin mapping ##
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
occs_raw <- read.csv("IUCN_point_files/haplocarpha_hastata_IUCN_pointfile.csv")
occ_points <-  occs_raw %>% 
    filter.occurences(T) 

# hydrobasins
hydro_bas <- vect("large/hybas_af_lev01-12_v1c/hybas_af_lev08_v1c.shp", crs = '+proj=longlat')

## 4. Intersect Points and hydrobasins 
# vectorising occs for intersect 
occ_vect <- vect(geom = c('long', 'lat'), crs = '+proj=longlat')

# intersect
hydro_intersect <- hydro_bas[occ_vect]

## 5. Merging SIS data with hydrobasins 
# extracting required data for SIS
sis_dataframe <- as.data.frame(occs_raw[1,]) %>%
    select(-dec_lat, -dec_long, -spatialref, -event_year, -basisofrec, -catalog_no, 
           -recordedby, -recordno)
sis_dataframe$source <- NA

# converting vect to sf selecting required data and merging
hydro_with_sis <- st_as_sf(hydro_intersect) %>%
    select(HYBAS_ID, geometry) %>%
    sp::merge(sis_dataframe) %>%
    vect()

## 6. Calculating EOO and AOO for points
point_stats <- cal.eoo.aoo(occ_points)

polygon_stats <- polygon.eoo.aoo(hydro_with_sis)

## 7. Map
make.poly.map(occ_points, hydro_with_sis)



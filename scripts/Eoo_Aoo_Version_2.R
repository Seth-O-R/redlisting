##### EOO, AOO and AOH for Use by the CAA Assessment Team. #####
## Seth Ratcliffe - s.ratcliffe@kew.org
## 25/11/2020

################ EOO, AOO and Mapping ######################
## 1. Packages and Libraries ----
install.packages(c("sf", "leaflet", "raster", "rCAT", "pacman", "tidyverse"))  
pacman::p_load(sf, leaflet, raster, rCAT, tidyverse) 

## 2. Functions ----
filter.occurences <- function(occurences, T){
  if(T){
    occurences %>%
      dplyr::select("presence", "dec_lat", "dec_long") %>%
      filter(! presence %in% c(4,6)) %>%
      dplyr::select(-"presence") %>%
      rename(lat = dec_lat, long = dec_long)
  } else {
    occurences %>%
      dplyr::select("dec_lat", "dec_long") %>%
      rename(lat = dec_lat, long = dec_long)
  }
}

## 2. Data Import and cleaning ----
occ_points <- read.csv("IUCN point files/Merxmuellera_grandiflora_IUCN.csv") %>% 
  filter.occurences(T) # If you want to use all occurrences change the T to an F

## 3. EOO and AOO ----
# EOO
center <- trueCOGll(occ_points) 
eoo_proj <- rCAT::simProjWiz(occ_points, center)
eoom2 <- EOOarea(eoo_proj)
eookm2 <- abs(eoom2/1000000) #calculates km2

eoo_cat <- EOORating(eookm2, T) # calculating the category

# AOO
cell_size_m <- 2000 # sets to 2x2 km grid
aoo_no_cells <- AOOsimp (eoo_proj, cell_size_m)
aookm2 <- aoo_no_cells * (cell_size_m/1000)^2

aoo_cat <- AOORating(aookm2, T) # calculating the catagory

print(paste0(" EOO = ", ceiling(eookm2)," ", eoo_cat, "; AOO = ", aookm2," ", aoo_cat)) 
# The above code prints numbers and categories in the Script tab below.

## 4. Mapping ----
# swapping lat long columns round for mapping and making spatial
points_spat <- occ_points %>% 
  select(long, lat) %>% 
  st_as_sf(coords = c("long", "lat"), crs = "4236")

# making boundary
boundary <- st_convex_hull(st_union(points_spat)) %>% 
  st_as_sf(type = 3)

# Mapping 
leaflet() %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "esri") %>%
  addScaleBar(position = "bottomright") %>%
  addCircleMarkers(data = points_spat, 
                   color = "blue", 
                   stroke = F, 
                   fillOpacity = 0.8) %>%
  addPolygons(data = boundary,
              color = "red", weight = 2, fill = F) 

################ Area of Habitat ###########################
## 1. Sourcing AOH functions ----
source("aoh_functions.R")
hab_crosswalk <- read.csv("hab_crosswalk.csv")

## 2. Data Import ----
est_range <- st_read("polys/merx_poly_small.kml", type = 3) # use this if loading in from a Google earth KML

## 3. Parameters ----
elevmin <- 3760 # Varies dependent on species 
elevmax <- 4620
mask <- est_range # This can either be the est_range KML or the boundary object from the EOO calculation                             
DEMrast <- raster::raster("eth_DEM_100.tif") # elevation data
habstack <- raster::raster("eth_jung.tif")
ESA_codes <- data.frame(ESA_codes = 407) # This will vary dependent on the habitat type 

## 4. Generate the AOH ----
theDEM <- dem(DEMrast, mask, elevmin, elevmax)
theHAB <- habitat(habstack, mask,  ESA_codes)
theAOH <- aoh(theHAB, theDEM, mask)

## 5. Generate the Red List stats from AOH ---- 
# EOO
aoh_eoo_df <- as.data.frame(rasterToPoints(theAOH))
aoh_eoo_ll <- aoh_eoo_df %>%
  dplyr::rename(lat = y, long = x) %>%
  dplyr::select(-layer)
aoh_center <- trueCOGll(aoh_eoo_ll)
aoh_eoo_proj <- rCAT::simProjWiz(aoh_eoo_ll, aoh_center)
aoh_eoom2 <- EOOarea(aoh_eoo_proj)
aoh_eookm2 <- abs(aoh_eoom2/1000000) #calculates km2

# AOO
aoh_cell_size_m <- 2000 # sets to 2x2 km grid
aoh_aoo_no_cells <- AOOsimp (aoh_eoo_proj, aoh_cell_size_m)
aoh_aookm2 <- aoh_aoo_no_cells * (aoh_cell_size_m/1000)^2

print(paste0(" EOO = ", ceiling(eookm2)," ", eoo_cat, "; AOO = ", aookm2," ", aoo_cat, "; AOH = ", aoh_aookm2))
print(paste0(" EOO = ", ceiling(aoh_eookm2), "; AOO = ", aoh_aookm2))


## 6. Map View ---- 
leaflet() %>%
  addTiles() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "ersi") %>%
  addScaleBar(position = "bottomright") %>%
  addRasterImage(theAOH, group = "aoh", color = "red", opacity = 0.5) %>%
  addPolygons(data = est_range, # change this depending on range source
              color = "black",
              weight = 1,
              fillColor = "yellow",
              group = "shape") %>% 
  addPolygons(data = est_range, # change this depending on range source
              color = "black", weight = 2, fill = F) %>%
  addCircleMarkers(data = points_spat, 
                   color = "blue", 
                   stroke = F, 
                   fillOpacity = 0.8) 


print(paste0(" EOO = ", ceiling(eookm2)," ", eoo_cat, "; AOO = ", aookm2," ", aoo_cat, "; AOH = ", aoh_aookm2))
print(paste0(" EOO = ", ceiling(aoh_eookm2), "; AOO = ", aoh_aookm2))


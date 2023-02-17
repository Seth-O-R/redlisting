##### EOO and AOO for Use by the CAA Assessment Team. #####
## Seth Ratcliffe - s.ratcliffe@kew.org
## 25/11/2020


################ EOO, AOO and Mapping ######################
## 1. Packages and Libraries
install.packages(c("sf", "leaflet", "raster", "rCAT", "pacman", "tidyverse", 
                   "stars")) # here we are installing the required packages. 

# loading packages
pacman::p_load(sf, leaflet, raster, rCAT, tidyverse, stars) # here we are loading the packages we need for this session. 

## 2. Functions 
source("scripts/functions.R")

## 3. Data Import and cleaning
occ_points <- read.csv("IUCN_point_files/phagnalon_phagnaloides_iucn_point.csv") %>% 
  filter.occurences(T) # If you want to use all occurrences change the T to an F

## 4. Calculate EOO and AOO
eoo_aoo <- cal.eoo.aoo(occ_points)

## 5. Produce Map 
make.map(occ_points)

############### AOH ###########################
## 2. Data Import ----
est_range <- st_read("polys/echinops_inf_range.kml", type = 3) # use this if loading in from a Google earth KML
boundary <- make.boundary(occ_points)

## 3. Parameters ----
# setting min and max elevation
elevmin <- 1830 # Varies dependent on species 
elevmax <- 3200

# creating mask 
mask <- boundary # This can either be the est_range KML or the boundary object from the EOO calculation                             

# elevation raster
DEMrast <- raster::raster("large/eth_DEM_100.tif") # elevation data

# habitat raster
habstack <- raster::raster("large/eth_jung.tif")

# Defining habitat codes
ESA_codes <- data.frame(ESA_codes = c(105, 307)) # This will vary dependent on the habitat type 

## 4. Generate the AOH ----
theDEM <- dem(DEMrast, mask, elevmin, elevmax)
theHAB <- habitat(habstack, mask,  ESA_codes)
theAOH <- aoh(theHAB, theDEM, mask)

## 5. Generate the Red List stats from AOH ---- 
cal.aoh.stats(theAOH)

## 6. Map View ---- 
make.aoh.map(occ_points, theAOH, F)

## 7. Export AOH .shp for upload to SIS ----
# converting aoh raster to sf
aoh_sf <- theAOH %>%
    st_as_stars() %>% # converting to stars object for sf transformation
    st_as_sf(as_points = F, merge = T) # converting to sf and merging points

# writing .shp file  
st_write(aoh_sf, "aoh_outs/phagnalon_phagnaloides_aoh.shp")

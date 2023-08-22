# External Data Cleaning and Standardizing 
# Seth Ratcliffe - s.ratcliffe@kew.org
# 31/03/2023

# Packages and libraries ----
install.packages('CoordinateCleaner')
pacman::p_load('tidyverse', 'CoordinateCleaner', 'sf', 'leaflet', 'rWCVP', 'terra')

# Shiny GeoCAT ---- 
# Data import 
raw_occs <- read.csv("shiny_geo_raw/trifolium_mattirolianum_shiny_raw.csv") %>%
    mutate(sci_name = paste(genus, specificEpithet))

# Cliping to native range
# getting native range
native_range <- wcvp_distribution(paste(raw_occs$sci_name[1]), taxon_rank = 'species',
                                  introduced = F, extinct = F, location_doubtful = F)

# making points spatial
occs_spatial <- raw_occs %>%
    st_as_sf(coords = c('longitude', 'latitude'), crs = st_crs(4326))

# ploting points against native range
(p_nat <- wcvp_distribution_map(native_range, crop_map = T) +
        theme(legend.position = 'none') +
        geom_sf(data = occs_spatial, fill = 'black', col = 'black', shape = 21))

# merging range polygons and adding buffer
buff_dist <- native_range %>%
    st_transform(3857) %>% # transforming to mercator to enable buffer
    st_union() %>%
    st_buffer(0.009) %>% # this is roughly 1km at equator and will expand and need to be altered at high/low latitudes
    st_transform(4326) # transforming back to wgs84

# Tagging occs to native range
occs_spatial$native_buffer <- st_intersects(occs_spatial, buff_dist, 
                                            sparse = F)[,1]

# filtering non-native 
occs_native <- occs_spatial %>%
    filter(native_buffer == T)

(p_clean <- wcvp_distribution_map(native_range, crop_map = T) +
        theme(legend.position = 'none') +
        geom_sf(data = occs_native, fill = 'black', col = 'black', shape = 21))    

# converting native range to spatvect for extraction
range_vect <- vect(native_range)

# extracting country coords and extracting lat/longs from geometery
occs_with_cc <- extract(range_vect, vect(occs_native)) %>%
    cbind.data.frame(occs_native, .) %>%
    mutate(longitude = unlist(map(.$geometry, 1)), 
           latitude = unlist(map(.$geometry, 2))) %>%
    select(-geometry) 

# Cleaning coords
# creating flags
flags <- occs_with_cc %>%
    clean_coordinates(lon = 'longitude', 
                      lat = 'latitude', 
                      species = 'sci_name', 
                      countries = 'LEVEL3_COD')

# Plotting for investigation 
leaflet() %>%
    addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "esri") %>%
    addScaleBar(position = "bottomright") %>%
    addCircleMarkers(data = filter(flags, .summary == 'TRUE'),
                     ~longitude, ~latitude,
                     color = 'blue', 
                     stroke = F, 
                     fillOpacity = 0.8) %>%
    addCircleMarkers(data = filter(flags, .summary == 'FALSE'),
                     ~longitude, ~latitude, 
                     color = 'red',
                     stroke = F, 
                     fillOpacity = 0.8)

# summary stats 
summary(flags)

# filter flag records if necesserary 
occs_flag_rm <- flags %>%
    filter(.summary == T)

# Data wrangling for SIS and merge with external data
data_cleaned <- occs_with_cc %>% # this will vary based on if you flitered occs with coordinate cleaner
    rename(dec_lat = latitude, dec_long = longitude) %>%
    add_column(subspecies = NA, subpop = NA, sens_comm = NA, dist_comm = NA, tax_comm = NA) %>%
    select(sci_name, presence, origin, seasonal, compiler, yrcompiled, citation, dec_lat, dec_long, 
           spatialref, subspecies, subpop, data_sens, sens_comm, event_year, source, basisOfRecord, catalogNumber,
           recordedBy, recordNumber, dist_comm, tax_comm) %>%
    mutate(presence = 1, origin = 1, compiler = 'Seth Ratcliffe', yrcompiled = 2023, citation = 'Royal Botanic Garden, Kew',
           data_sens = 0)

# Genesys ----
# loading in data 
gen_cor <- read.csv("genesys/trifolium_mattirolanum_genesys/core.csv")
gen_geo <- read.csv("genesys/trifolium_mattirolanum_genesys/geo.csv", row.names = NULL) %>%
    select(-uncertainty, -method) %>%
    rename(genesysId1 = row.names, latitude1 = genesysId, longitude1 = latitude, elevation1 = longitude) %>%
    select(-elevation, -datum) %>%
    rename(genesysId = genesysId1, latitude = latitude1, longitude = longitude1, elevation = elevation1) # for some reason the geo file is loading weird and using the genesysId coloumn as row names

# merging by genesys ID
gen_joined <- gen_geo %>% 
    merge(gen_cor, by = 'genesysId', all = T)

# filtering occs with no lat long
gen_filt <- gen_joined %>%
    filter(!is.na(latitude) & !is.na(longitude))

# data wrangling for merging and SIS
gen_clean <- gen_filt %>%
    unite('sci_name', genus:species) %>%
    rename(dec_lat = latitude, dec_long = longitude, event_year = acqDate,
           source = doi, recordedBy = instCode, recordNumber = genesysId, catalogNumber = acceNumb) %>%
    add_column(presence = 1, origin = 1, seasonal = NA, compiler = 'Seth Ratcliffe', 
               yrcompiled = 2023, citation = 'Royal Botanic Garden, Kew', spatialref = 'WGS84', basisOfRecord = 'Living_Specimen', 
               subspecies = NA, subpop = NA, data_sens = 0, sens_comm = NA, dist_comm = NA, tax_comm = NA) %>%
    select(sci_name, presence, origin, seasonal, compiler, yrcompiled, citation, dec_lat, dec_long, 
           spatialref, subspecies, subpop, data_sens, sens_comm, event_year, source, basisOfRecord, catalogNumber,
           recordedBy, recordNumber, dist_comm, tax_comm)
    
# Brahms ---- 
# loading data
ext_data <- read.csv("IUCN_point_files/trifolium_mattirolianum_iucn_pointfile.csv") %>%
    rename(basisOfRecord = basisofrec, catalogNumber = catalog_no, recordedBy = recordedby, 
           recordNumber = recordno)
ext_data$source <- as.character(ext_data$source) 
ext_data$catalogNumber <- as.character(ext_data$catalogNumber)

summary(data_cleaned)
summary(ext_data)

# Merging data ----    
# joining to cleaned data 
data_all <- list(data_cleaned, gen_clean, ext_data) %>% # creating list for join
    bind_rows()

# removing duplicated records
duplicate_rm <- data_all %>%
    cc_dupl(lon = 'dec_long', 
            lat = 'dec_lat', 
            species = 'sci_name')

# converting NAs so blanks for later
duplicate_rm[is.na.data.frame(duplicate_rm)] <- ""
data_all[is.na.data.frame(data_all)] <- ""
    
# Data Export ---- 
write.csv(data_all, file = "shiny_geo_clean/trifolium_mattirolianum_standrarised_all.csv")
write.csv(duplicate_rm, file = "shiny_geo_clean/trifolium_mattirolianum_IUCN_pointfile_no_dups.csv")

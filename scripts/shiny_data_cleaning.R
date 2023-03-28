# ShinyGeoCAT Data Cleaning 
# Seth Ratcliffe - s.ratcliffe@kew.org
# 31/03/2023

# Packages and libraries ----
install.packages('CoordinateCleaner')
pacman::p_load('tidyverse', 'CoordinateCleaner', 'sf', 'leaflet', 'rWCVP', 'terra')

# Data import ----
raw_occs <- read.csv("shiny_geo_raw/justicia_heterocarpa_shiny_raw.csv") %>%
    mutate(sci_name = paste(genus, specificEpithet))

# Cliping to native range ----
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

# extracting country coords and extracting latg longs from geometery
occs_with_cc <- extract(range_vect, vect(occs_native)) %>%
    cbind.data.frame(occs_native, .) %>%
    mutate(longitude = unlist(map(.$geometry, 1)), 
           latitude = unlist(map(.$geometry, 2))) %>%
    select(-geometry) 

# Cleaning coords ----
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

# Data Export ---- 
# Data wrangling for SIS 
data_cleaned <- occs_with_cc %>% # this will vary based on if you flitered occs with coordinate cleaner
    rename(dec_lat = latitude, dec_long = longitude) %>%
    add_column(subspecies = NA, subpop = NA, sens_comm = NA, dist_comm = NA, tax_comm = NA) %>%
    select(sci_name, presence, origin, seasonal, compiler, yrcompiled, citation, dec_lat, dec_long, 
           spatialref, subspecies, subpop, data_sens, sens_comm, event_year, source, basisOfRecord, catalogNumber,
           recordedBy, recordNumber, dist_comm, tax_comm)

write.csv(data_cleaned, file = "shiny_geo_clean/justicia_heterocarpa_iucn_pointfile.csv")

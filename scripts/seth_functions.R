# Seth redlisting functions - 28/09/23

# filter.occurences - filter occurences coded as 4, 6 on upload 
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

# cal.eoo.aoo - almagamation of rCAT functions to reduce code and calculate EOO and AOO
cal.eoo.aoo <- function(occurences) {
    
    # calculate eoo
    centre <- trueCOGll(occurences)
    eoo_projection <- rCAT::simProjWiz(occurences, centre)
    eoom2 <- EOOarea(eoo_projection)
    eookm2 <- abs(eoom2/1000000) # calculates km2
    
    eoo_cat <- EOORating(eookm2, T)
    
    # calculate aoo
    cell_size_m <- 2000 # sets to 2x2 km grid
    aoo_no_cells <- AOOsimp (eoo_projection, cell_size_m)
    aookm2 <- aoo_no_cells * (cell_size_m/1000)^2
    
    aoo_cat <- AOORating(aookm2, T) # calculating the catagory
    
    # print eoo and aoo with categories
    return(print(paste0(" EOO = ", ceiling(eookm2)," ", eoo_cat, "; AOO = ", aookm2," ", aoo_cat))) 
    
}

# make.boundary - makes boundary object from occurrence points
make.boundary <- function(a, b, eoo = T, buffer = F){
    
    if(eoo == T & buffer == F) {
        
        # make points spatial
        points_spat <- a %>% 
            select(long, lat) %>% 
            st_as_sf(coords = c("long", "lat"), crs = "4236")
        
        # make mcp
        st_convex_hull(st_union(points_spat)) %>% 
            st_as_sf(type = 3)
        
    } else if(eoo == T & buffer == T){ 
        
        # make points spatial
        points_spat <- a %>% 
            select(long, lat) %>% 
            st_as_sf(coords = c("long", "lat"), crs = "4236")
        
        # make mcp
        st_convex_hull(st_union(points_spat)) %>% 
            st_as_sf(type = 3) 
        
        # add buffer
        buf <- st_buffer(points_spat, dist = b)
        buf_bound <- st_union(buf, points_spat) %>%
            st_as_sf(type = 3)
        
    } else if(eoo == F & buffer == T){
        
        # make points spatial
        points_spat <- a %>% 
            select(long, lat) %>% 
            st_as_sf(coords = c("long", "lat"), crs = "4236")
        
        # add buffer
        buf <- st_buffer(points_spat, dist = b)
        st_union(buf, points_spat)
    }
}

# calc.aoh.sing - New AOH function for single species as the old one was not clipping by elevation
calc.aoh.sing <- function(dem, habitat, hab_codes, boundary, elv_min, elv_max) {
    suppressWarnings({  
        # cropping elevation by boundary
        dem_crop <- crop(dem, boundary)
        
        # reducing values down to elevation boundary, everything elsa NA
        dem_crop[dem_crop < elv_min] <- NA 
        dem_crop[dem_crop > elv_max] <- NA 
        
        # cropping hab
        hab_crop <- crop(habitat, boundary)
        
        # converting hab to df to filter by hab codes
        hab_vals <- terra::as.data.frame(hab_crop, xy = T)
        
        colnames(hab_vals) <- c('lat', 'long', 'aoh')
        
        reduced_hab <- hab_vals %>%
            filter(.[,3 ] %in% hab_codes)
        
        red_hab_rast <- rast(reduced_hab, type = "xyz", crs = 'WGS84')
        
        # hab masked by elv 
        aoh_rast <- mask(red_hab_rast, dem_crop)
        aoh_rast[!is.na(aoh_rast)] <- 1
        
        aoh_to_eoo <- mask(aoh_rast, boundary)
        
        # calculating aoh EOO and AOH re-scalled to 2x2 
        # EOO
        aoh_eoo_df <- terra::as.data.frame(aoh_to_eoo, xy = T)
        
        aoh_eoo_ll <- aoh_eoo_df %>%
            dplyr::rename(lat = y, long = x) %>%
            select(lat, long)
        
        aoh_center <- trueCOGll(aoh_eoo_ll)
        aoh_eoo_proj <- rCAT::simProjWiz(aoh_eoo_ll, aoh_center)
        aoh_eoom2 <- EOOarea(aoh_eoo_proj)
        aoh_eookm2 <- abs(aoh_eoom2/1000000) #calculates km2
        
        # AOO
        aoh_cell_size_m <- 2000 # sets to 2x2 km grid
        aoh_aoo_no_cells <- AOOsimp (aoh_eoo_proj, aoh_cell_size_m)
        aoh_aookm2 <- aoh_aoo_no_cells * (aoh_cell_size_m/1000)^2
        
        print(paste0("AOH EOO = ", ceiling(aoh_eookm2), "; AOH AOO = ", ceiling(aoh_aookm2)))
        
        return(aoh_to_eoo)
    })
}

# cal.aoh.stats -  Calls AOH stats 
cal.aoh.stats <- function(the_aoh){
    # EOO
    aoh_eoo_df <- terra::as.data.frame(theAOH, xy = T)
    
    aoh_eoo_ll <- aoh_eoo_df %>%
        dplyr::rename(lat = y, long = x) %>%
        select(lat, long)
    
    aoh_center <- trueCOGll(aoh_eoo_ll)
    aoh_eoo_proj <- rCAT::simProjWiz(aoh_eoo_ll, aoh_center)
    aoh_eoom2 <- EOOarea(aoh_eoo_proj)
    aoh_eookm2 <- abs(aoh_eoom2/1000000) #calculates km2
    
    # AOO
    aoh_cell_size_m <- 2000 # sets to 2x2 km grid
    aoh_aoo_no_cells <- AOOsimp (aoh_eoo_proj, aoh_cell_size_m)
    aoh_aookm2 <- aoh_aoo_no_cells * (aoh_cell_size_m/1000)^2
    
    print(paste0(eoo_aoo, "; AOH EOO = ", ceiling(aoh_eookm2), "; AOH AOO = ", ceiling(aoh_aookm2)))
}

# make.map - produces map using leaflet 
make.map <- function(occurences){
    
    # swapping lat long columns round for mapping and making spatial
    points_spat <- occurences %>% 
        select(long, lat) %>% 
        st_as_sf(coords = c("long", "lat"), crs = "4236")
    
    # creating boundary
    mcp <- make.boundary(occurences)
    
    # Mapping 
    leaflet() %>%
        addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
        addProviderTiles(providers$Esri.WorldImagery, group = "esri") %>%
        addScaleBar(position = "bottomright") %>%
        addCircleMarkers(data = points_spat, 
                         color = "blue", 
                         stroke = F, 
                         fillOpacity = 0.8) %>%
        addPolygons(data = mcp,
                    color = "red", weight = 2, fill = F)  # Need to figure out how to add AOO squares onto the map
}

# make.aoh.map - make map for aoh using leaflet
make.aoh.map <- function(occurences, the_aoh, boundary_aoh = F, aoh_raster = F){
    
    if(aoh_raster == T & boundary_aoh == F){
        
        # making occurrences spatial
        points_spat <- occurences %>% 
            select(long, lat) %>% 
            st_as_sf(coords = c("long", "lat")) %>%
            st_set_crs(4326)
        
        boundary_object <- st_convex_hull(st_union(points_spat)) %>% 
            st_as_sf() %>%
            st_set_crs(4326)
        
        # producing AOH map
        leaflet() %>%
            addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
            addProviderTiles(providers$Esri.WorldImagery, group = "ersi") %>%
            addScaleBar(position = "bottomright") %>%
            addRasterImage(the_aoh,
                           group = "aoh",
                           color = "red", 
                           opacity = 0.5) %>%
            addCircleMarkers(data = points_spat, 
                             color = "blue", 
                             stroke = F, 
                             fillOpacity = 0.8) %>%
            addPolygons(data = boundary_object, 
                        color = "black",
                        weight = 1,
                        fillColor = "yellow",
                        group = "shape") %>% 
            addPolygons(data = boundary_object, 
                        color = "black", weight = 2, fill = F)
        
    } else if(aoh_raster == T & boundary_aoh == T){
        
        # making occurences spatial
        points_spat <- occurences %>% 
            select(long, lat) %>% 
            st_as_sf(coords = c("long", "lat")) %>%
            st_set_crs(4326)
        
        # making boundary on aoh
        boundary_object <- as.data.frame(rasterToPoints(the_aoh)) %>%
            rename(lat = x, long = y) %>%
            st_as_sf(coords = c("lat", "long")) %>%
            st_set_crs(4326)
        
        boundary_object <- st_convex_hull(st_union(boundary_object)) %>% 
            st_as_sf() %>%
            st_set_crs(4326)
        
        # producing AOH map
        leaflet() %>%
            addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
            addProviderTiles(providers$Esri.WorldImagery, group = "ersi") %>%
            addScaleBar(position = "bottomright") %>%
            addRasterImage(the_aoh,
                           group = "aoh",
                           color = "red",
                           opacity = 0.5) %>%
            addCircleMarkers(data = points_spat, 
                             color = "blue", 
                             stroke = F, 
                             fillOpacity = 0.8) %>%
            addPolygons(data = boundary_object, 
                        color = "black",
                        weight = 1,
                        fillColor = "yellow",
                        group = "shape") %>% 
            addPolygons(data = boundary_object, 
                        color = "black", weight = 2, fill = F)
        
    }else if(aoh_raster == F & boundary_aoh == F){
        
        # making occurrences spatial
        points_spat <- occurences %>% 
            select(long, lat) %>% 
            st_as_sf(coords = c("long", "lat")) %>%
            st_set_crs(4326)
        
        # making boundary
        boundary_object <- st_convex_hull(st_union(points_spat)) %>% 
            st_as_sf() %>%
            st_set_crs(4326)
        
        # producing AOH map
        leaflet() %>%
            addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
            addProviderTiles(providers$Esri.WorldImagery, group = "ersi") %>%
            addScaleBar(position = "bottomright") %>%
            addPolygons(data = the_aoh,
                        color = "red",
                        fillColor = "red",
                        opacity = 0.5) %>%
            addCircleMarkers(data = points_spat, 
                             color = "blue", 
                             stroke = F, 
                             fillOpacity = 0.8) %>%
            addPolygons(data = boundary_object, 
                        color = "black",
                        weight = 1,
                        fillColor = "yellow",
                        group = "shape") %>% 
            addPolygons(data = boundary_object, 
                        color = "black", weight = 2, fill = F)
        
    } else if(aoh_raster == F & boundary_aoh == T){
        
        # making occurrences spatial
        points_spat <- occurences %>% 
            select(long, lat) %>% 
            st_as_sf(coords = c("long", "lat")) %>%
            st_set_crs(4326)
        
        # making boundary
        boundary_object <- st_convex_hull(st_union(the_aoh)) %>% 
            st_as_sf() %>%
            st_set_crs(4326)
        
        # producing AOH map
        leaflet() %>%
            addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
            addProviderTiles(providers$Esri.WorldImagery, group = "ersi") %>%
            addScaleBar(position = "bottomright") %>%
            addPolygons(data = the_aoh,
                        color = "red",
                        fillColor = "red",
                        opacity = 0.5) %>%
            addCircleMarkers(data = points_spat, 
                             color = "blue", 
                             stroke = F, 
                             fillOpacity = 0.8) %>%
            addPolygons(data = boundary_object, 
                        color = "black",
                        weight = 1,
                        fillColor = "yellow",
                        group = "shape") %>% 
            addPolygons(data = boundary_object, 
                        color = "black", weight = 2, fill = F)
        
    }
}


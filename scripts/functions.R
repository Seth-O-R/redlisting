# AOH functions from Steve


hab_crosswalk <- read.csv("hab_crosswalk.csv")

# aoh_simple - combined AOH but no integrated analysis
aoh_simple = function(spnam,
                      powoID,
                      rast_mask,
                      ESA_codes,
                      raster_years,
                      DEMrast,
                      hab_stack,
                      min_elev,
                      max_elev){
  print(paste0("running ", spnam))
  
  #### 2 - make the DEM
  myDEM = dem(DEMrast, rast_mask, min_elev, max_elev)
  print("(1/3 - DEM ready)")
  #plot(myDEM)
  
  #### 3 - generate habitat raster(s)
  myhab = habitat(hab_stack, rast_mask, raster_years, ESA_codes)
  print("(2/3 - Habitat ready)")
  
  #### 4 - make AOH 
  myaoh = aoh(myhab, myDEM, rast_mask, raster_years)
  print("(3/3 - Area of habitat ready)")
  
  return(myaoh)
}

# dem - crop DEM to mask
dem = function(DEMrast, mask, elevmin, elevmax){
  
  sp_DEM <- raster::crop(DEMrast, mask) 
  
  # reduce to elevation range - everything else NA
  sp_DEM[sp_DEM < elevmin] <- NA 
  sp_DEM[sp_DEM > elevmax] <- NA 
  
  vals <- unique(values(sp_DEM)) %>% 
    as.data.frame(.)
  
  # some format wrangling to get ready to join 
  colnames(vals)[which(names(vals) == ".")] <- "orig"
  
  vals$new <- 1
  vals <- vals[-1, ]
  
  # reclassify the habitat raster
  sp_DEM <- raster::reclassify(sp_DEM, rcl = vals) 
  
  return(sp_DEM)
  
}



# habitat - reclassify the habitat layers or stack to suitable habitat
habitat = function(habstack, mask, ESA_codes){
  
  cropped <- raster::crop(habstack, mask)
  
  reclassed_rast = recl_hab(cropped, ESA_codes)
  
  return(reclassed_rast)
  
}

recl_hab = function(singl_hab, ESA_codes){
  
  # get the raster vals
  vals <- unique(values(singl_hab)) %>%
      as.data.frame()
  
  colnames(vals)[which(names(vals) == ".")] <- "code"
  
  # # get the ESA codes and reformat to do the reclassification
  # ESA_codes <- read.table(text = ESA_codes, sep =",", 
  #                         header = FALSE, 
  #                         stringsAsFactors = FALSE)
  # 
  # ESA_codes = ESA_codes %>%
  #   tidyr::pivot_longer(cols = starts_with("V"), 
  #                       #names_to = "val", 
  #                       values_to = "code")
  
  ESA_codes$val <- 1
  ESA_codes$code <- ESA_codes$ESA_codes
  ESA_codes = ESA_codes %>% dplyr::select(-1) 
  
  # now join to make the reclassification table
  recl <- dplyr::full_join(vals, ESA_codes, by = "code")
  
  recl = as.matrix(recl)
  
  #recl <- matrix(recl, ncol=2, byrow=TRUE)
  
  #run the reclassifier
  rast =raster::reclassify(singl_hab, recl)
  
}

# aoh - make area of habitat
aoh = function(hab, dem, mask ){
  
  # combine habitat and dem 
  AOH = hab * dem
  
  # clip to mask
  AOH <- mask(AOH, mask) 
  
  return(AOH)
  
}


# SETH FUNCTIONS

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

# make.boundary 
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


# Seth AOH functions
# cal.aoh.stats -  amalgamation of steves AOH functions to reduce code and calculate AOH
cal.aoh.stats <- function(the_aoh){
  
  # EOO
  aoh_eoo_df <- as.data.frame(rasterToPoints(the_aoh))
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
  
  print(paste0(eoo_aoo, "; AOH EOO = ", ceiling(aoh_eookm2), "; AOH AOO = ", ceiling(aoh_aookm2)))
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


# CSV point data to kml conversion script 
# Seth Ratcliffe (23/08/2023)

# Packages and Libraries
install.packages(c("sf", "sp", "pacman")) 

pacman::p_load(sf, sp) # exporting occurences as shape file with boundary

# Importing csv file
occs_raw <- read.csv("............") # input file path 
 
# converting csv to spatial sf class and adding WGS84 crs
points_spat <- occ_points %>% 
    select(long, lat) %>% # you will need to make sure the latitude/longitude columns are named accordingly
    st_as_sf(coords = c("long", "lat"), crs = "4236") %>%
    st_buffer(dist = 0.01) #1 km buffer

# exporting kml file
st_write(points_spat, "................", 
         driver = 'kml') # input output path where the dots are

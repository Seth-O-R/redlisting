##### GBIF and Brahms occ data pull and summary #####
## Seth Ratcliffe - s.ratcliffe@kew.org
## 15/06/2023

## 1. Packages and Libraries ----
install.packages(c("rgbif", "tidyverse", "readr"))
pacman::p_load("tidyverse", "readr", "rgbif", 'vctrs')

## 2. Data input ---- 
spp_names <- read.csv("c:/Users/sra10kg/OneDrive - The Royal Botanic Gardens, Kew/GCBC/Rapid_lc/eth_p2_names.csv")
brahms_raw <- read.csv("c:/Users/sra10kg/OneDrive - The Royal Botanic Gardens, Kew/GCBC/brahms_data_base_pull.CSV")

## 3. GBIF Pull ----
# getting taxon keys 
taxon_keys <- spp_names %>%
    name_backbone_checklist() %>%
    filter(!matchType == 'None') %>%
    pull(usageKey)

# searching gbif database 
gbif_data <- occ_search(taxonKey = taxon_keys, fields = 'all')

# extracting data section of gbif data and converting to a list 
gbif_list <- vector('list', length(gbif_data))
names(gbif_list) <- as.vector(spp_names$name_in) # setting names 

for (s in as.vector(spp_names$name_in)) {
    data <- gbif_data[[s]]$data
    gbif_list[[s]] <- data.frame(data)
    
}

# cutting empty list sectioin (species with no records in GBIF)
gbif_no_empt <- list_drop_empty(gbif_list) %>%
    lapply(data_frame)

# filtering out uncessarry columns in prep for collapsing into single df
gbif_standard <- vector('list', length(gbif_no_empt))
names(gbif_standard) <- names(gbif_no_empt)
    
for (i in as.vector(names(gbif_standard))) {
    data <- gbif_no_empt[[i]]
    
        if(is.null(data$decimalLatitude) || is.null(data$decimalLongitude)) {
             
            data <- data %>%
                mutate(decimalLatitude = 0, decimalLongitude = 0)
            
            data <- data %>%
                select('scientificName', 'acceptedScientificName', 'family', 'genus', 'specificEpithet', 'basisOfRecord', 
                       'catalogNumber', 'institutionCode', 'collectionCode', 'decimalLatitude', 'decimalLongitude',
                       'gbifID')
        } else {
            data_stand <- data %>%
                select('scientificName', 'acceptedScientificName', 'family', 'genus', 'specificEpithet', 'basisOfRecord', 
                       'catalogNumber', 'institutionCode', 'collectionCode', 'decimalLatitude', 'decimalLongitude',
                       'gbifID')
        }
    
    gbif_standard[[i]] <- data.frame(data_stand)
    
}

# collapsing into dataframe 
gbif_df <- as.data.frame(do.call(rbind, gbif_standard))

## 4. Data summaries ----
# filtering and summarizing brahms data 
brahms_clean <- brahms_raw %>%
    unite('name', GENUS:SP1, sep = " ") %>%
    filter(name %in% spp_names$name_in)

brahms_summary <- brahms_clean %>% 
    group_by(name) %>%
    summarise(no_specimes_brahms = n())

# summarizing gbif data
gbif_summary <- gbif_df %>%
    group_by(acceptedScientificName) %>%
    summarise(no_specimens_gbif = n())


## 5. Data Output ----
write.csv(brahms_clean, 'brahms_pull_phase_2.csv')
write.csv(gbif_summary, 'gbif_summary_phase_2.csv')
write.csv(brahms_summary, 'brahms_summary_phase_2.csv')

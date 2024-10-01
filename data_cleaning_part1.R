################################################################################
### Pre-processing of raw data: Part 1
###
### Script for cleaning occurrence records
### Part of the methods for the manuscript:
### How well do we understand geographic range size?: 
### A case study of Australia’s frogs and citizen science projects
###
###
################################################################################
###
###
###
################################################################################
### Load required libraries and data objects into R

### Required libraries
library(tidyverse)
library(galah)
library(CoordinateCleaner)
library(countrycode)
library(sf)

### Load raw data

# Load ALA data
load("data/ala_frog_occurrences_raw_20231017.Rda") #ala_all_rawdata
ala_occurrences_raw <- ala_all_rawdata #change data object name

# Load FrogID data (FrogID datset 4.0).
# Note that this study obtained FrogID data through licence agreement.
# The dataset includes sensitive records and cannot be shared publicly.
# FrogID dataset can be downloaded from https://www.frogid.net.au/explore
# Sensitive records obtained directly from the website
# have been desensitised with reduced geolocation accuracy.

####
# IMPORTANT: f you do not have the desensitised FrogID data,
# skip rows 41 - 80 since desensitised records of FrogID are included in ALA download 
###
frogid_occurrences_raw <- read_csv("data/FrogID4_final_dataset_SENSITIVE_INCLUDED_AUSTRALIA.csv")
nrow(frogid_occurrences_raw)
#[1] 484665

## Check unique species names from FrogID 
frogid_species <-sort(unique(frogid_occurrences_raw$scientificName))
checked_frogid_species <-search_taxa(frogid_species) 
# Three names with fuzzyMatch
# "Crinia nimbus", "Cyclorana vagitus", "Litoria nudidigitus"

### Standardise scientific names to match Australian Faunal Directory
frogid_occurrences_raw$scientificName<- gsub("Crinia nimbus", "Crinia nimba",
                                        gsub("Cyclorana vagitus", "Cyclorana vagita",
                                        gsub("Litoria nudidigitus", "Litoria nudidigita",
                                        gsub("Lechriodus fletcheri","Platyplectrum fletcheri", #synonym
                                        frogid_occurrences_raw$scientificName))))

save(ala_occurrences_raw, file ="data/ala_occurrences_raw_20240617_namestandardised.Rda")
save(frogid_occurrences_raw, file = "data/frogid_occurrences_raw_20240617_namestandardised.Rda")


### Load rawdata with standardised names and select relevant columns for data filtering
load("data/ala_occurrences_raw_20240617_namestandardised.Rda")
load("data/frogid_occurrences_raw_20240617_namestandardised.Rda")

### Check the total number of species
length(unique(frogid_occurrences_raw$scientificName)) #208 species
length(unique(ala_occurrences_raw$species)) #242 species

### Group ALA data by dataResourceName which is the data resource that supplies the record.

occurrences_by_insitution <- ala_occurrences_raw %>%
  group_by(dataResourceName) %>%
  summarise(Records=n())  
# FrogID is the largest contributor ofoccurrence records 
# (N = 394,140 records which is ~ 1/3 (32.2%) of raw records)

### Exclude records in ALA which were supplied by FrogID to remove duplicates
ala_occurrences_raw <- subset(ala_occurrences_raw, dataResourceName != "FrogID")
nrow(ala_occ_raw) 
# 826,181 records retained


### Plot distribution of frog species occurrences for visualisation ####
ala_occurrences <- ala_occurrences_raw %>% 
  dplyr::select("species", "decimalLatitude", "decimalLongitude")%>%
  mutate(datatype = "ALA")


frogid_occurrences <- frogid_occurrences_raw %>% 
  dplyr::select("scientificName", "decimalLatitude", "decimalLongitude") %>% 
  mutate(datatype = "FrogID")

colnames(frogid_occurrences)[colnames(frogid_occurrences)=="scientificName"] <-"species"
frog_alloccurrences_raw <-bind_rows(ala_occurrences,
                                    frogid_occurrences)

################################################################################
### Data cleaning processes for ALA dataset

### Select and keep relevant fields
ala_rawdata <- ala_occ_raw %>%
  dplyr::select("scientificName","species","taxonRank","decimalLongitude",
                "decimalLatitude", "dataResourceName", "basisOfRecord",
                "eventDate","day","month","year","sensitive","country","cl1048", 
                "coordinatePrecision","coordinateUncertaintyInMeters","stateProvince",
                "recordedBy", "locality","identifiedBy","institutionCode","matchType")


### Flag records for basis of records, taxonomy and time
ala_preprocess <- ala_rawdata %>%
  mutate(
  # Basis of record test 
    basisOfRecord_test = ifelse(basisOfRecord != "MATERIAL_SAMPLE" & basisOfRecord != "MATERIAL_CITATION", TRUE, FALSE), 
  # Taxon rank and match type test       
    taxon_matchtype_test = ifelse((taxonRank == "species" | taxonRank == "subspecies") & matchType != "higherMatch", TRUE, FALSE), 
  # Old record test
    time_test = ifelse(year >= 1950 & !is.na(year), TRUE, FALSE))

### Flag duplicated records (i.e., additional records of the same taxon at the same location and collection date)

ala_preprocess$duplicate_test <-cc_dupl(
  x=ala_preprocess, 
  lon = "decimalLongitude", 
  lat = "decimalLatitude",
  species = "species",
  value= "flagged", 
  additions = c("month","year","day")) 

# Flagged 254,297 records as duplicates

### Flag records with geospatial issues

ala_preprocess <- ala_preprocess %>% 
  mutate(
  # Records without coordinates test
    coord_test = ifelse(!is.na(decimalLongitude) & !is.na(decimalLatitude), TRUE, FALSE), 
  # Flag records outside Australia
    country_test = ifelse(country == "Australia" & !is.na(country), TRUE, FALSE), 
  # Flag records with just zeros in decimals
    dec.zeros_test =ifelse(!(grepl('^[^\\.]+$|\\.0*$', decimalLatitude)) & !(grepl('^[^\\.]+$|\\.[0*]$',decimalLongitude)), TRUE, FALSE)) 


# Flag records with more geographic coordinates issues 
# (https://ropensci.github.io/CoordinateCleaner/articles/Cleaning_GBIF_data_with_CoordinateCleaner.html)
# Flag potential erroneous coordinates that fall in centriods, equal, institutions, zeros, and seas

# Convert country name to ISO3c 3-letter-code
ala_preprocess$countryCode <- ala_preprocess$country
ala_preprocess$countryCode <- countrycode(ala_preprocess$countryCode,
                                         origin = "country.name",
                                         destination = "iso3c")

ala_preprocess <- {
  
  # Split dataframe into two subsets
  with_coords <- ala_preprocess%>% filter(coord_test == TRUE)
  without_coords <- ala_preprocess %>% filter(coord_test == FALSE)
  
  # Apply coordinate cleaning to a subset with coordinates
  cleaned_coords <- clean_coordinates(x=with_coords,
                                      lon = "decimalLongitude", 
                                      lat= "decimalLatitude",
                                      countries = "countryCode",
                                      species ="species",
                                      tests = c("captials", 
                                                "centroids",
                                                "equal",
                                                "institutions",
                                                "gbif",
                                                "seas",
                                                "zeros"), 
                                      value = "spatialvalid")
  
  # Combine the cleaned subset with the original subset without coordinates
  ala_preprocess <- bind_rows(cleaned_coords, 
                             without_coords)    
  return(ala_preprocess)
  
}


summary(cleaned_coords)
#.val     .equ     .zer     .cen     .sea     .gbf    .inst .summary 
#  0        0        0       52    21012        0      270    21321 
#.val tests for coordinate validity, .equ for equal lat/long, .zer for zero coordinates, .cen = country capitals


## Flag records with coordinate uncertainty higher than 2 km
ala_preprocess <- ala_preprocess %>% 
  # test for records with coordinate uncertainty above 2 km
  mutate(coord_pres_test = ifelse(coordinateUncertaintyInMeters <= 2000 & !is.na(coordinateUncertaintyInMeters), TRUE, FALSE))

save(ala_preprocess, file="data/ala_preprocess_withflaggedrecords_20240618.Rda")

### Remove records that have been flagged
ala_initial_clean <- ala_preprocess %>% 
  filter(basisOfRecord_test == TRUE) %>%
  filter (taxon_matchtype_test == TRUE) %>%
  filter (time_test == TRUE) %>%
  filter (duplicate_test == TRUE) %>%
  filter (coord_test == TRUE) %>%
  filter (country_test == TRUE) %>% # false 228+ na 4 #check this later
  filter (dec.zeros_test == TRUE)  %>%
  filter (.summary == TRUE) %>%
  filter (coord_pres_test == TRUE)

# Create a column for data aggregator 
ala_initial_clean$dataprovider <- "ALA"
length(unique(ala_initial_clean$species)) #220 species remain

### Save ALA clean data
save(ala_initial_clean, file = "data/ala_initial_clean_20240618.Rda")

################################################################################
### Data cleaning processes for FrogID dataset
### IMPORTANT: If using FrogID data aggregated through ALA, skip the rest of the code

### Data harmonisation and integration of fields for FrogID data
frogid_occ_raw <- frogid_occurrences_raw 
frogid_occ_raw$species <- frogid_occ_raw$scientificName
frogid_occ_raw$taxonRank <- "species" 
frogid_occ_raw$dataResourceName <-frogid_occ_raw$datasetName
frogid_occ_raw$eventDate_original <-frogid_occ_raw$eventDate
frogid_occ_raw <-frogid_occ_raw %>% separate(eventDate_original, sep = "/", into = c("day","month","year"))
frogid_occ_raw$sensitive <-frogid_occ_raw$dataGeneralizations

save(frogid_occ_raw, file = "Data/frogid_occurrences_raw_20240618_namestandardised_addcol.Rda")

nrow(frogid_occ_raw) 
# 484665 records

### Select and keep relevant fields
frogid_rawdata <- frogid_occ_raw %>%
  dplyr::select("scientificName","species","taxonRank","decimalLongitude",
                "decimalLatitude", "dataResourceName", "basisOfRecord",
                "eventDate","day","month","year","sensitive","country",
                "coordinateUncertaintyInMeters","stateProvince","recordedBy")                                               

### Flag records for basis of records, taxonomy and time
frogid_preprocess <- frogid_rawdata %>%
  mutate(
  # Basis of record test  
    basisOfRecord_test = ifelse(basisOfRecord != "MATERIAL_SAMPLE" & basisOfRecord != "MATERIAL_CITATION", TRUE, FALSE),
  # Taxon rank and match type test    
    taxon_matchtype_test = ifelse(taxonRank == "species" | taxonRank == "subspecies", TRUE, FALSE), 
  # Old record test
    time_test = ifelse(year >= 1950 & !is.na(year), TRUE, FALSE)) 

### Flag duplicated records (i.e. additional records of the same taxon at the same location and collection date)

frogid_preprocess$duplicate_test <-cc_dupl(
  x=frogid_preprocess, 
  lon = "decimalLongitude", 
  lat = "decimalLatitude",
  species = "species",
  value= "flagged", 
  additions = c("month","year","day")) #Flagged 37,810 records

### Flag records with geospatial issues

frogid_preprocess <- frogid_preprocess %>% 
  mutate(
  # Records without coordinates tes
    coord_test = ifelse(!is.na(decimalLongitude) & !is.na(decimalLatitude), TRUE, FALSE), 
  # Flag records outside Australia
    country_test = ifelse(country == "Australia" & !is.na(country), TRUE, FALSE), 
  # Flag records with just zeros in decimals
    dec.zeros_test =ifelse(!(grepl('^[^\\.]+$|\\.0*$', decimalLatitude)) & !(grepl('^[^\\.]+$|\\.[0*]$',decimalLongitude)), TRUE, FALSE)) 


# Flag more  geographic coordinates issues 
# (https://ropensci.github.io/CoordinateCleaner/articles/Cleaning_GBIF_data_with_CoordinateCleaner.html)
# Flag potential erroneous coordinates that fall in centriods, equal, institutions, zeros, and seas

# Convert country name to ISO3c 3-letter-code
frogid_preprocess$countryCode <-frogid_preprocess$country
frogid_preprocess$countryCode <-countrycode(frogid_preprocess$countryCode,
                                            origin = "country.name",
                                            destination = "iso3c")

frogid_preprocess <- {
  # Split dataframe into two subsets
  with_coords <- frogid_preprocess %>% filter(coord_test == TRUE)
  without_coords <- frogid_preprocess %>% filter(coord_test == FALSE)
  
  # Apply coordinate cleaning to a subset with coordinates
  cleaned_coords <- clean_coordinates(with_coords,
                                      lon = "decimalLongitude", 
                                      lat= "decimalLatitude",
                                      countries = "countryCode",
                                      species ="species",
                                      tests = c("captials", 
                                                "centroids",
                                                "equal",
                                                "institutions",
                                                "gbif",
                                                "seas",
                                                "zeros"), 
                                      value = "spatialvalid")
  
  # Combine the cleaned subset with the original subset without coordinates
  frogid_preprocess <-bind_rows(cleaned_coords, 
                                without_coords)    
  return(frogid_preprocess)
  
}

summary(cleaned_coords)
#.val     .equ     .zer     .cen     .sea     .gbf    .inst .summary 
#   0        0        0        1    13125        0      263    13387 
#.val tests for coordinate validity, .equ for equal lat/long, .zer for zero coordinates, .cen = country capitals


### Flag records with coordinate uncertainty higher than 2 km
frogid_preprocess <- frogid_preprocess %>% 
  # Test for records with coordinate uncertainty above 2 km
  mutate(coord_pres_test = ifelse(coordinateUncertaintyInMeters <= 2000 & !is.na(coordinateUncertaintyInMeters), TRUE, FALSE))
# 22,077 records flagged

save(frogid_preprocess, file="Data/frogid_preprocess_withflaggedrecords_20240618.Rda")


## remove records which have been flagged
frogid_initial_clean <- frogid_preprocess %>% 
  filter(basisOfRecord_test == TRUE) %>%
  filter (taxon_matchtype_test == TRUE) %>%
  filter (time_test == TRUE) %>%
  filter (duplicate_test == TRUE) %>%
  filter (coord_test == TRUE) %>%
  filter (country_test == TRUE) %>% 
  filter (dec.zeros_test == TRUE)  %>%
  filter (.summary == TRUE) %>%
  filter (coord_pres_test == TRUE)

# Create a column for data aggregator name
frogid_initial_clean$dataprovider <- "FrogID"
length(unique(frogid_initial_clean$species)) #207 species remain

# save FrogID cleaned data

save(frogid_initial_clean, file = "data/frogid_initial_clean_20240618.Rda")

################################################################################ 
### Combining ala_initial_clean with frogid_initial_clean two one data object
# load("data/ala_initial_clean_20240618.Rda")
# load("data/frogid_initial_clean_20240618.Rda")

### Data harmonisation and integration of fields
# Modify FrogID data frame to match ALA data
frogid_initial_clean$eventDate <-as.Date(frogid_initial_clean$eventDate, 	
                                         format="%d/%m/%Y")

frogid_initial_clean <- frogid_initial_clean %>% 
                        mutate_at(c("day","month","year"), as.numeric)
frogid_initial_clean <- frogid_initial_clean %>% 
                        mutate_at(c("recordedBy"), as.character)

save(frogid_initial_clean, file = "data/frogid_initial_clean_20240618.Rda") 

# Replace old oject 
# Columns and structures standardised
## Bind ALA data cleaned with FrogID data
allfrog_occ_intial_clean <- bind_rows(ala_initial_clean,
                                      frogid_initial_clean)

allfrog_occ_intial_clean$stateProvince <- gsub("TAS", "Tasmania",
                                          gsub("VIC", "Victoria", 
                                          gsub("SA","South Australia",
                                          gsub("NSW","New South Wales",  
                                          gsub("ACT", "Australian Capital Territory",
                                          gsub("WA","Western Australia", 
                                          gsub("QLD","Queensland",    
                                          gsub("NT","Northern Territory",
                                          allfrog_occ_intial_clean$stateProvince))))))))

save(allfrog_occ_intial_clean , file = "data/allfrog_occ_intial_clean_20240618.Rda") 

################# End of Data Processing: Part 1 ###############################
################################################################################

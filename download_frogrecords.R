################################################################################
### Downloading data from Atlas of Living Australia (ALA) with galah package  
###
### Script for download data from ALA using galah package
### Part of the methods for the manuscript:
### How well do we understand geographic range size?: 
### A case study of Australia’s frogs and citizen science projects
###
### 
################################################################################
### Exact same dataset can be directly downloaded from ALA using the links provided
### See rows 146-152 
################################################################################ 
### Install packages and load required libraries
# Install galah package - the development version from GitHub:
#install.packages("remotes")
#remotes::install_github("AtlasOfLivingAustralia/galah")
  
# Load library and configure email

library(galah)
library(tidyverse)

galah_config(email = "youremailaddress", atlas="Australia", verbose=TRUE) 
# Change this to own email
preserve = TRUE


################################################################################
### Load species list and download occurrence records
# Load species list obtained from Australian Faunal Directory (AFD)
spp_list <- read_csv("data/Australian-frogs-species-list-from-AFD-SpeciesLevel.csv")
frog_list <- as.character(unique(spp_list$SPECIES_NAME)) #248 unique species names

# Search species name in ALA
name_check <-search_taxa(frog_list)
save(name_check, file = "data/frogspp_namecheck_ALA.Rda")

species <-subset(name_check, match_type == "exactMatch") 
# 237 species names with exact match at the time this was checked
# Note that ALA had updated their Taxonomic Backbone in April 2024

# Species as character
spp_download_list <-as.character(species$scientific_name)

### Download species occurrence records where taxonomic names had exactMatch in ALA
raw_data <-galah_call() %>%
  galah_identify(spp_download_list) %>% 
  # Subset species for download if spp_download_list cannot be downloaded all at once 
  # e.g. raw_data <- galah_identify(spp_download_list[1:60])
  # e.g. raw_data2 <- galah_identify (spp_download_list[61:120])
  galah_filter(profile="ALA")%>%
  galah_select(group = "basic",# Return basic fields; https://rdrr.io/cran/galah/man/galah_select.html
               vernacularName, # The common name the ALA has matched this record to
               phylum,
               class,
               order, 
               family,
               genus,
               speciesGroup, # Higher level group for this record e.g. Birds
               species,
               specificEpithet,
               taxonRank,
               taxonomicIssues,
               sensitive, # This fields is populated if the record is sensitive.
               basisOfRecord, 
               day,
               month,
               year,
               country,
               stateConservation, # The state conservation status for the taxon.
               countryConservation,
               geodeticDatum, 
               establishmentMeans, 
               coordinatePrecision,
               coordinateUncertaintyInMeters, 
               stateProvince, 
               spatiallyValid, # Whether this record is suitable for use in species distribution
               cl22, # Australian States and Territories Australian States and Territories
               cl1048, # IBRA 7
               cl10000, # Forest of Australia 2018B
               cl617, # Vegetation types- native Pre-European major vegetation groups
               cl620, # Vegetation types - present current major vegetation groups
               recordedBy,
               recordedByID,
               identifiedBy,
               identifiedByID,	
               institutionName,
               institutionCode,
               isDuplicateOf,
               license,
               locality,
               matchType, # An indication of how a match to a concept was achieved.
               outlierLayer)%>%
  atlas_occurrences(mint_doi = TRUE) 
# Download occurrences with the all of the specified field names

### atlas_occurrences() could not download all records at once
# so we subset the species list and download multiple times 
#
# For species where match_type != "exactMatch"
# search_taxa() could not find information for "Platyplectrum fletcheri",
# "Anstisia alba", "Anstisia lutea","Anstisia rosea","Anstisia vitellina"
# So we searched https://amphibiansoftheworld.amnh.org/ for synonyms and found
# "Lechriodus fletcheri", "Geocrinia alba","Geocrinia lutea","Geocrinia rosea","Geocrinia vitellina"
#
# Create spp_download_list2 to donwload additional species
spp_download_list2 <- c("Lechriodus fletcheri", "Geocrinia alba","Geocrinia lutea","Geocrinia rosea","Geocrinia vitellina")
#
# Replace spp_download_list on line 49 with spp_download_list2
# to download additional spp. in the raw data
# then download species occurrence data for five species 
# (replacing raw_data with raw_data_synonyms)

### Replace synonyms and standardise with AFD taxonmic names
# Use nested gsub to replace multiple values for species and genus names to match AFD
raw_data_synonyms_standardised <- raw_data_synonyms
raw_data_synonyms_standardised$scientificName <- gsub("Lechriodus fletcheri", "Platyplectrum fletcheri",
                                                 gsub("Geocrinia alba", "Anstisia alba",
                                                 gsub("Geocrinia lutea", "Anstisia lutea",
                                                 gsub("Geocrinia rosea","Anstisia rosea",
                                                 gsub("Geocrinia vitellina", "Anstisia vitellina", 
                                                 raw_data_synonyms_standardised$scientificName)))))


raw_data_synonyms_standardised$genus <- gsub("Geocrinia", "Anstisia",
                                        gsub("Lechriodus", "Platyplectrum",
                                        raw_data_synonyms_standardised$genus))

### combine downloaded data into one object and save
ala_all_rawdata <- bind_rows(raw_data,
                         raw_data_synonyms_standardised)
nrow(ala_all_rawdata)
#[1] 1220321 
# total of 1220321 records downloaded

save(ala_all_rawdata, file = "data/ala_frog_occurrences_raw_20231017.Rda" )  


################################################################################ 
### DOI for all download
# a)  first download - spp_download_list[1:60], occurrences = 305,254 records (https://doi.ala.org.au/doi/fb7754b9-f084-4ee3-b4ab-a87dd313d467; https://doi.org/10.26197/ala.fb7754b9-f084-4ee3-b4ab-a87dd313d467).
# b)	second download - spp_download_list[61:120], occurrences = 523,471 records (https://doi.ala.org.au/doi/1e0ded65-7fa5-439f-8298-081537b871a7; https://doi.org/10.26197/ala.1e0ded65-7fa5-439f-8298-081537b871a7).
# c)	third download - spp_download_list[121:180], occurrences = 268,351 records (https://doi.ala.org.au/doi/abd60ea2-f837-4f21-afdc-a5d41340c820; https://doi.org/10.26197/ala.abd60ea2-f837-4f21-afdc-a5d41340c820).
# d)	fourth download - spp_download_list[181:237], occurrences = 121,255 records (https://doi.ala.org.au/doi/d14ae6e3-c8b6-4d7a-bb75-3dbe25e9a753; https://doi.org/10.26197/ala.d14ae6e3-c8b6-4d7a-bb75-3dbe25e9a753).
# e)	spp_download_list2, 5 synonyms, occurrences = 1990 records (https://doi.ala.org.au/doi/42348355-ed0e-40ac-853f-98fa30a2d5ac; https://doi.org/10.26197/ala.42348355-ed0e-40ac-853f-98fa30a2d5ac).
### copy and paste links to download
################### End of code for downloading data ###########################
################################################################################
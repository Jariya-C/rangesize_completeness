################################################################################
### Pre-processing of raw data: Part 2
###
### Script for cleaning occurrence records
### Part of the methods for the manuscript:
### How well do we understand species’ geographic range size?: 
### A case study of Australia’s frogs
###
### Jariya Chanachai (jariya.chanachai@hdr.mq.edu.au)
################################################################################

### Install packages, load required libraries and data objects
# install.packages("foreach")
# install.packages("doParallel")

library(foreach)
library(doParallel)
library(spdep)
library(tidyverse)
library(sf)

### Load all frog occurrences 
load("data/allfrog_occ_intial_clean_20240618.Rda")

allrecords <- as.data.frame(allfrog_occ_intial_clean)

### Create NeighbourDataFrame

# Create a unique ID for each occurrence record from 1:length of dataframe
allrecords$recordID <-as.numeric(rownames(allrecords))

# Load Interim Biogeographic Regionalisation for Australia (IBRA)
ibra_shp <- st_read("data/IBRA7_regions/ibra7_regions.shp")

# Check coordinate system
st_crs(ibra_shp) 

# Transform to EPSG:4326
ibra_shp <- st_transform(ibra_shp, 4326) 

# Transform records into sf object
allrecords_sf <- st_as_sf(allrecords, 
                          coords = c("decimalLongitude","decimalLatitude"), crs=st_crs(ibra_shp))

ibra_shp <- st_make_valid(ibra_shp)

# Perform spatial join to obtain IBRA attributes for each occurrence
ibra_occ_join <- st_join(allrecords_sf, ibra_shp)


# Extract x and y coordinates for each record
ibra_occ_join_df <- (ibra_occ_join) %>% 
  dplyr::mutate(decimalLongitude = sf::st_coordinates(.)[,1],
                decimalLatitude = sf::st_coordinates(.)[,2]) %>% as.data.frame()


# Retain only records which intersect with an IBRA
# REG_CODE_7 is not NA
ibra_occ_joined_df <- ibra_occ_join_df %>% 
  filter(!is.na(REG_CODE_7))
# Records removed = 639

# One species, Uperoleia inundata, has recently been merged into Uperoleia crassa
# Replace species name with accepted name

ibra_occ_joined_df$scientificName <- gsub("Uperoleia inundata", "Uperoleia crassa", 
                                          ibra_occ_joined_df$scientificName)

ibra_occ_joined_df$species <- gsub("Uperoleia inundata", "Uperoleia crassa", 
                                   ibra_occ_joined_df$species)


################################################################################

### Extract IBRA and neighbouring IBRA values which intersect with species distribution.

# Set row names of IBRA to 3-letter-code REG_CODE_7
row.names(ibra_shp) <- as.character(ibra_shp$REG_CODE_7)

# Construct neighbours list from polygon list (spdep package)
queen_nb <- poly2nb(ibra_shp, queen = TRUE, row.names = as.character(ibra_shp$REG_CODE_7))

#plot(queen_nb, coords= st_coordinates(thPoly_AM))
#source("1_codes/functions_.R")

#1:29
# Function for converting neighbours list into dataframe
neighborsDataFrame <- function(nb) {
  stopifnot( inherits(nb, 'nb'))
  ks = data.frame(k = unlist( mapply(rep, 1:length(nb), sapply(nb, length), SIMPLIFY = FALSE) ), k_nb = unlist(nb) )
  nams = data.frame(id = attributes(nb)$region.id, k = 1:length(nb) )
  o = merge(ks, nams, by.x = 'k', by.y = 'k')
  o = merge(o, nams, by.x = 'k_nb', by.y = 'k', suffixes = c("","_neigh"))
  o[, c("id", "id_neigh")]
}

# Convert neighbour list to dataframe
df1 <- neighborsDataFrame(queen_nb)

# Rename columns as "src" for source IBRA and "nbr" for adjacent IBRA

colnames(df1) <- c("src", "nbr");nrow(df1)


# Find the intersect between IBRA and frog distribution map 
# And return IBRA that intersects with species range 
# Use parallel processing to speed up the process

unique_spp <- sort(unique(ibra_occ_joined_df$species))


# Register a parallel backend using doParallel
cores <- detectCores()
cl <- makeCluster(cores-7) # use just 5 cores
registerDoParallel(cl)

# Initialize the list to store results
spp_ibra_intersect_list <- list()

### WARNING - The below loop will take some time to run ###

# Parallel processing with foreach
spp_ibra_intersect_list <- foreach(i = 1:length(unique_spp), .packages = c("sf")) %dopar% {
  
  species <- gsub(" ","_", unique_spp[i])
  
  # Load species map for species i - Data source = Australian Frog Atlas (AFA)
  species_map <- st_read(paste0("data/australian_frog_atlas_frogid/species_map_v2_shp/", species,".shp"))
  
  # Transform species map
  species_map <- st_transform(species_map,4326) # EPSG:4326
  
  species_map <- st_make_valid(species_map)
  
  # Find the intersection between AFA and bioregion (IBRA)
  afa_ibra_intersect <- st_intersection(species_map, ibra_shp)
  
  # Create a dataframe to share data and obtain unique IBRA ID which intersect with AFA
  afa_ibra_intersect_df <-cbind.data.frame("species" = unique_spp[i], "IBRA_ID" = afa_ibra_intersect$REG_CODE_7)
  
  # Return the result for this iteration
  return(afa_ibra_intersect_df)
  # Store each data frame in a list for each species.
  #spp_ibra_intersect[[i]] <- afa_ibra_intersect_df
  
  
}

# Stop the parallel backend
stopCluster(cl)


# Combine the list of data frames into a single data frame
spp_ibra_intersect_df <- do.call(rbind, spp_ibra_intersect_list)

save(spp_ibra_intersect_df, file = "data/frog_distribution_ibra_intersect_df.Rda")

################################################################################

### Run the range test to remove records which fall outside defined range

# Steps for obtaining the src and nrc IBRA for each species
unique_spp <- sort(unique(ibra_occ_joined_df$species))
#unique_spp <- unique_spp[1:3]

occ_range_test <- ibra_occ_joined_df

# Create a new column called range_test
occ_range_test$range_test <- NA

occ_cleaned_list <- list()
for (i in 1:length(unique_spp)) {
  
  # Subset IBRA for species [i]
  spp_ibra <- subset(spp_ibra_intersect_df, species == unique_spp[i])
  
  # Filter for the source IBRA for species [i]
  nbr_ <- filter(df1, src %in% unique(spp_ibra$IBRA_ID))
  
  # Obtain a unique list of IBRA for source and neighbouring IBRAs for species [i]
  nbr_ <- unique(c(nbr_$src, nbr_$nbr))
  
  # Filter records for species [i]
  frog_occ_current <- occ_range_test %>% filter(species == unique_spp[i])
  
  # Test whether the occurrences fall within the source or neighbouring IBRA
  frog_occ_current$range_test <- ifelse(frog_occ_current$REG_CODE_7 %in% nbr_, TRUE, FALSE)
  
  # Return the result in a list
  occ_cleaned_list[[i]] <- frog_occ_current
  
  
}

# Combine the list into a single data frame of frog occurrence for all species
occ_range_test_df <- do.call(rbind, occ_cleaned_list)
summary(factor(occ_range_test_df$range_test))

save(occ_range_test_df, file = "data/occ_range_test_20240626.Rda")

################################################################################

### Plot distribution of records for each frog species species
### To identify records outside of range

records_by_spp <- occ_range_test_df %>% 
                  group_by(species) %>% 
                  summarise(records = n())

rangetest_true <- occ_range_test_df %>% 
                  filter(range_test == "TRUE") %>% 
                  group_by(species) %>% 
                  summarise(cleaned_records = n())

records_by_spp <- left_join(records_by_spp, rangetest_true, by = "species")

unique_spp <- sort(unique(ibra_occ_joined_df$species))

### WARNING - This loop generates distribution map for each spp ###
# Takes a while to run

for(i in 1:length(unique_spp)) {
  
  species <- gsub(" ","_", unique_spp[i])
  spp_records <- subset(records_by_spp, species == unique_spp[i])
  
  # Load species map for species i. Data source = Australian Frog Atlas (AFA)
  species_map <- st_read(paste0("data/australian_frog_atlas_frogid/species_map_v2_shp/", species,".shp"))
  
  # Subset the data for the current species
  occurrences <- occ_range_test_df[occ_range_test_df$species == unique_spp[i], ]
  
  # Plot distribution of records for each species
  occurrence_map <- ggplot() + 
    geom_sf(data = ibra_shp, fill = "#FBFBEF", linewidth=0.1, colour = "black")+ 
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = (occurrences), 
               aes(x =decimalLongitude, 
                   y=decimalLatitude, 
                   colour = factor(dataprovider),
                   shape = ifelse(occurrences$range_test =="TRUE",  "within range", "outliers")),
               size = 1) +
    scale_shape_manual(values = c(4, 1)) +
    geom_sf(data = species_map, alpha =0.2, fill = NA, colour = "red") +
    coord_sf(expand=TRUE, xlim = c(110,175), ylim = c(-45,-10)) + #limit extent to Oceania
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    guides(shape = guide_legend(override.aes = list(size = 3))) +
    annotate("text", label = unique_spp[i], x = 125, y = -10, size =6,  fontface = 'italic') +
    annotate("text", label = bquote("total records = " ~.(spp_records$records)), x = 125, y = -45, size =6) +
    annotate("text", label = bquote("records within range = " ~.(spp_records$cleaned_records)), x = 165, y = -45, size =6) +
    theme_classic() +
    theme(legend.position = c(0.9,0.8),
          legend.title = element_blank(),
          legend.text = element_text(colour="black", size =15),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.key.size = unit(1,"cm"),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15),
          panel.background = element_blank())
  ggsave(plot = occurrence_map, paste0("result/species_distribution_with_afa_frogidmap/",unique_spp[i],"_distribution.png"), dpi = 1200, width = 12, height = 6)
}

################################################################################

### Final cleaning
# Removing records which fail the range_test or have been translocated
# Remove species which have been split into two or more species
# This info is inferred from FrogID website & literature

# load("data/occ_range_test_20240626.Rda")

frog_occ_cleaned_final <- occ_range_test_df %>% 
  filter(!(species %in% c("Limnodynastes terraereginae",
                          "Litoria rothii",
                          "Litoria lesueuri",
                          "Litoria dentata",
                          "Geocrinia victoriana",
                          "Limnodynastes dumerilii",
                          "Litoria ewingii",
                          "Mixophyes balbus",
                          "Litoria watjulumensis"))) %>% #remove L. watjulumensis 
  # The above step will remove 87959 records belong to the above species  
  filter(range_test == "TRUE") %>%
  # The above step removes additional 111 records found to be outside range   
  #
  # Remove ClimateWatch records for L. infrafrenata
  filter(!(species == "Litoria infrafrenata" & dataResourceName == "ClimateWatch")) %>% 
  filter(recordID != 281996) %>% # Remove one record for L. infrafrenata from NSW_BioNet Altas 
  filter(!(recordID %in% c(297179,297265, 295992))) %>% # Remove 3 records for Limnodynastes convexiusculus 
  filter(!(recordID %in% c(309927,475952))) %>% #Remove two records for L. rubella with known translocation 
  filter(recordID != 377017) %>% # Remove 1 record for Litoria verreauxii in MUL
  filter(!(recordID %in% c(362344,363657, 363766, 371065))) %>% # Remove 3 records for Neobatrachus pictus in SYB and 1 in BBS
  filter(recordID != 415446) %>% # Remove one record for Pseudophryne raveni in SYB  
  filter(recordID != 431357) #remove 1 record for Pseudophryne coriacea in MUL

# Remaining records = 777472

### Divide records into citizen science and non-citizen science datasets
# Identify citizen science projects from dataResource Name
unique(frog_occ_cleaned_final$dataResourceName)

# Create vector of citizen science data
cs_data <- c("FrogID", "iNaturalist Australia", "ClimateWatch", 
             "Melbourne Water Frog Census", "ALA species sightings and OzAtlas", 
             "Gaia Guide","Earth Guardians Weekly Feed")

# Create a column for classifying records as citizen/non-citizen science data

frog_occ_cleaned_final <- frog_occ_cleaned_final %>%
  mutate(datatype = ifelse(dataResourceName %in% cs_data, "citizenscience", "non_citizenscience"))

summary(factor(frog_occ_cleaned_final$datatype))
# citizenscience     non_citizenscience  
#         394986                382486

save(frog_occ_cleaned_final, file = "data/frog_occ_cleaned_final_with_datatype_20240702.Rda")

# clean dataset = frog_occ_cleaned_final
# data cleaning completed
################################################################################

### Plot species distribution for clean data

unique_spp <- unique(frog_occ_cleaned_final$species)

records_by_spp <- frog_occ_cleaned_final %>% group_by(species) %>% summarise(records = n())

### This plots distribution map for each species - with clean dataset
### WARNING - This loop will take time to run ###

for(i in 1:length(unique_spp)) {
  
  species <- gsub(" ","_", unique_spp[i])
  spp_records <- subset(records_by_spp, species == unique_spp[i])
  
  # load species map for species i. Data source = Australian Frog Atlas (AFA)
  species_map <- st_read(paste0("Data/australian_frog_atlas_frogid/species_map_v2_shp/", species,".shp"))
  # Subset the data for the current species
  occurrences <- frog_occ_cleaned_final[frog_occ_cleaned_final$species == unique_spp[i], ]
  
  occurrence_map <- ggplot() + 
    geom_sf(data = ibra_shp, fill = "#FBFBEF", linewidth=0.1, colour = "black")+ 
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = (occurrences), 
               aes(x =decimalLongitude, y=decimalLatitude, 
                   colour = factor(dataprovider)),
               #shape = ifelse(occurrences$range_test =="TRUE",  "within range", "outliers")),  #color = ifelse(abs(x) < 1, "within range", "outliers")))
               size = 1) +
    scale_shape_manual(values = c(4, 1)) +
    geom_sf(data = species_map, alpha =0.2, fill = NA, colour = "red") +
    coord_sf(expand=TRUE, xlim = c(110,175), ylim = c(-45,-10)) + #limit extent to Oceania
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    guides(shape = guide_legend(override.aes = list(size = 3)))+
    annotate("text", label = unique_spp[i], x = 125, y = -10, size =6,  fontface = 'italic') +
    annotate("text", label = bquote("total records = " ~.(spp_records$records)), x = 125, y = -45, size =6) +
    #annotate("text", label = bquote("records within range = " ~.(spp_records$cleaned_records)), x = 165, y = -45, size =6) +
    theme_classic() +
    theme(legend.position = c(0.9,0.8),
          legend.title = element_blank(),
          legend.text = element_text(colour="black", size =15),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.key.size = unit(1,"cm"),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15),
          panel.background = element_blank())
  ggsave(plot = occurrence_map, paste0("Results/frog_species_distribution_cleaned_with_afa_frogidmap_20240702/",unique_spp[i],"_distribution.png"), dpi = 1200, width = 12, height = 6)
}
#
#
################# End of Data Processing: Part 1 ###############################
###################### End of data cleaning ####################################
################################################################################
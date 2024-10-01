################################################################################
### Pre-processing of raw data: Part 2
###
### Script for cleaning occurrence records
### Part of the methods for the manuscript:
### How well do we understand geographic range size?: 
### A case study of Australia’s frogs and citizen science
###
### 
################################################################################

### Install packages, load required libraries and data objects
# install.packages("foreach")
# install.packages("doParallel")

library(foreach)
library(doParallel)
library(spdep)
library(tidyverse)
library(sf)

# Load all frog occurrences 
load("data/allfrog_occ_intial_clean_20240618.Rda")
allrecords <- as.data.frame(allfrog_occ_intial_clean)

# neighbourDataFrame

# Create a unique ID for each occurrence record from 1:length of dataframe
allrecords$recordID <-as.numeric(rownames(allrecords))

# Interim Biogeographic Regionalisation for Australia (IBRA)
ibra_shp <- st_read("data/IBRA7_regions/ibra7_regions.shp")
st_crs(ibra_shp) 
ibra_shp <- st_transform(ibra_shp, 4326) #Transform to EPSG:4326

# Transform occurrences into sf object
allrecords_sf <- st_as_sf(allrecords, 
                          coords = c("decimalLongitude","decimalLatitude"), crs=st_crs(ibra_shp))

ibra_shp <- st_make_valid(ibra_shp)

# Perform spatial join to obtain IBRA attributes for each occurrence
ibra_occ_join <- st_join(allrecords_sf, ibra_shp)


#ibra_occ_join_df <- cbind(st_drop_geometry(ibra_occ_join), st_coordinates(ibra_occ_join))
ibra_occ_join_df <- (ibra_occ_join) %>% dplyr::mutate(decimalLongitude = sf::st_coordinates(.)[,1],
                                                      decimalLatitude = sf::st_coordinates(.)[,2]) %>% as.data.frame()

# Retain only occurrences which intersect with an IBRA (REG_CODE_7 is not NA)
ibra_occ_joined_df <- ibra_occ_join_df %>% 
  filter(!is.na(REG_CODE_7))
# records removed = 639


# One species, Uperoleia inundata, has been merged into Uperoleia crassa
# Replace species name with accepted name

ibra_occ_joined_df$scientificName <- gsub("Uperoleia inundata", "Uperoleia crassa", ibra_occ_joined_df$scientificName)

ibra_occ_joined_df$species <- gsub("Uperoleia inundata", "Uperoleia crassa", ibra_occ_joined_df$species)


################################################################################
### Extract IBRA and neighbouring IBRA values which intersect with species distribution map.

# Set row names of IBRA to 3-letter-code REG_CODE_7
row.names(ibra_shp) <- as.character(ibra_shp$REG_CODE_7)

## Construct neighbours list from polygon list
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

# Convert neighbour list to datafram
df1 <- neighborsDataFrame(queen_nb)

# Rename columns as "src" for source IBRA and "nbr" for adjacent IBRA

colnames(df1) <- c("src", "nbr");nrow(df1)



# Find the intersect between IBRA and frog distribution map and return IBRA that intersects with species range 

unique_spp <- sort(unique(ibra_occ_joined_df$species))


# Register a parallel backend using doParallel
cores <- detectCores()
cl <- makeCluster(cores-7) # choose to use just 5 cores
registerDoParallel(cl)

# Initialize the list to store results
spp_ibra_intersect_list <- list()

start_time <- Sys.time()
### WARNING: This loops takes about 20 minutes to run
# Parallel processing with foreach
spp_ibra_intersect_list <- foreach(i = 1:length(unique_spp), .packages = c("sf")) %dopar% {
  
  species <- gsub(" ","_", unique_spp[i])
  
  # Load species map for species i. Data source = Australian Frog Atlas (AFA)
  species_map <- st_read(paste0("data/australian_frog_atlas_frogid/species_map_v2_shp/", species,".shp"))
  
  # Transform species map
  species_map <- st_transform(species_map,4326) #EPSG:4326
  
  species_map <- st_make_valid(species_map)
  
  # Find the intersection between AFA and bioregion (IBRA)
  afa_ibra_intersect <- st_intersection(species_map, ibra_shp)
  
  # Create a dataframe to share data and obtain the unique IBRA ID which intersect with AFA
  afa_ibra_intersect_df <-cbind.data.frame("species" = unique_spp[i], "IBRA_ID" = afa_ibra_intersect$REG_CODE_7)
  
  
  # Return the result for this iteration
  return(afa_ibra_intersect_df)
  # Store each data frame in a list for each species.
  #spp_ibra_intersect[[i]] <- afa_ibra_intersect_df
  
  
}

# Stop the parallel backend
stopCluster(cl)


end_time <- Sys.time()
taken_time <- end_time - start_time


# Combine the list of data frames into a single data frame
spp_ibra_intersect_df <- do.call(rbind, spp_ibra_intersect_list)

save(spp_ibra_intersect_df, file = "data/frog_distribution_ibra_intersect_df.Rda")

################################################################################
### Run the data through the range test to remove records which fall outside defined range

### Steps for obtaining the src and nrc IBRA for each species
unique_spp <- sort(unique(ibra_occ_joined_df$species))
#unique_spp <- unique_spp[1:3]

frog_occ_initial_cleaned <- ibra_occ_joined_df

# Create a new column called range_test
frog_occ_initial_cleaned$range_test <- NA

occ_cleaned_list <- list()
for (i in 1:length(unique_spp)) {
  
  # Subset IBRA for species [i]
  spp_ibra <- subset(spp_ibra_intersect_df, species == unique_spp[i])
  
  # Filter for the source IBRA for species [i]
  nbr_ <- filter(df1, src %in% unique(spp_ibra$IBRA_ID))
  
  # Obtain a unique list of IBRA for source and neighbouring IBRAs for species [i]
  nbr_ <- unique(c(nbr_$src, nbr_$nbr))
  
  
  # Filter records for species [i]
  frog_occ_current <- frog_occ_initial_cleaned %>% filter(species == unique_spp[i])
  
  # Test whether the occurrences fall within the source or neighbouring IBRA
  frog_occ_current$range_test <- ifelse(frog_occ_current$REG_CODE_7 %in% nbr_, TRUE, FALSE)
  #frog_occ_initial_cleaned$range_test <- ifelse(frog_occ_initial_cleaned$REG_CODE_7 %in% nbr_, TRUE, FALSE)
  
  # Put the result in a list
  occ_cleaned_list[[i]] <- frog_occ_current
  
  
}

# Combine the list of data frames for each frog species into a single data frame of frog occurrence for all species
frog_occ_cleaned_df <- do.call(rbind, occ_cleaned_list)
summary(factor(frog_occ_cleaned_df$range_test))

################################################################################

### Plot distribution of records for each frog species species identitified with records outside of range


records_by_spp <- frog_occ_cleaned_df %>% group_by(species) %>% summarise(records = n())

records_by_spp_rangetest_true <- frog_occ_cleaned_df %>% filter(range_test == "TRUE") %>% 
  group_by(species) %>% summarise(cleaned_records = n())

records_by_spp <- left_join(records_by_spp, records_by_spp_rangetest_true, by = "species")

unique_spp <- sort(unique(ibra_occ_joined_df$species))

start_time <- Sys.time()


for(i in 1:length(unique_spp)) {
  
  species <- gsub(" ","_", unique_spp[i])
  spp_records <- subset(records_by_spp, species == unique_spp[i])
  
  # Load species map for species i. Data source = Australian Frog Atlas (AFA)
  species_map <- st_read(paste0("data/australian_frog_atlas_frogid/species_map_v2_shp/", species,".shp"))
  
  # Subset the data for the current species
  occurrences <- frog_occ_cleaned_df[frog_occ_cleaned_df$species == unique_spp[i], ]
  
  # Plot distribution of records for each species
  occurrence_map <- ggplot() + 
    geom_sf(data = ibra_shp, fill = "#FBFBEF", linewidth=0.1, colour = "black")+ 
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = (occurrences), 
               aes(x =decimalLongitude, y=decimalLatitude, 
                   colour = factor(dataprovider),
                   shape = ifelse(occurrences$range_test =="TRUE",  "within range", "outliers")),  #color = ifelse(abs(x) < 1, "within range", "outliers")))
               size = 1) +
    scale_shape_manual(values = c(4, 1)) +
    geom_sf(data = species_map, alpha =0.2, fill = NA, colour = "red") +
    coord_sf(expand=TRUE, xlim = c(110,175), ylim = c(-45,-10)) + #limit extent to Oceania
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    guides(shape = guide_legend(override.aes = list(size = 3)))+
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
  ggsave(plot = occurrence_map, paste0("Results/species_distribution_with_afa_frogidmap_20240626/",unique_spp[i],"_distribution.png"), dpi = 1200, width = 12, height = 6)
}


end_time <-Sys.time()

taken_time <- end_time - start_time

################################################################################

### Final cleaning
### Removing records which fail the range_test or have been translocated
### Remove species which have been split into two or more species
### IMPORTANT: Note that unique recordID are assigned to the dataset and will vary
### Depending on the dataset used

frog_occ_cleaned_final <- frog_occ_cleaned_df %>% 
  filter(!(species %in% c("Limnodynastes terraereginae",
                          "Litoria rothii",
                          "Litoria lesueuri",
                          "Litoria dentata",
                          "Geocrinia victoriana",
                          "Limnodynastes dumerilii",
                          "Litoria ewingii",
                          "Mixophyes balbus",
                          "Litoria watjulumensis"))) %>% # Remove L. watjulumensis 
  # The above step will remove 87959 records belong to the above species  
  filter(range_test == "TRUE") %>%
  # The above step removes additional 111 records found to be outside range   
  filter(!(species == "Litoria infrafrenata" & dataResourceName == "ClimateWatch")) #%>% # remove ClimateWatch records for L. infrafrenata
  # The below code is only relevant to our dataset
#  filter(recordID != 281996) %>% #record one record for L. infrafrenata from NSW_BioNet Altas 
#  filter(!(recordID %in% c(297179,297265, 295992))) %>% # remove 3 records for Limnodynastes convexiusculus 
#  filter(!(recordID %in% c(309927,475952))) %>% #remove two records for L. rubella with known translocation 
#  filter(recordID != 377017) %>% # remove 1 record for Litoria verreauxii in MUL
#  filter(!(recordID %in% c(362344,363657, 363766, 371065))) %>% # remove 3 records for Neobatrachus pictus in SYB and 1 in BBS
#  filter(recordID != 415446) %>% # remove one record for Pseudophryne raveni in SYB  
#  filter(recordID != 431357) #remove 1 record for Pseudophryne coriacea in MUL

# remaining records = 777472


save(frog_occ_cleaned_final, file ="data/frog_occ_cleaned_final_20240701.Rda")                       

################################################################################
##########Divide data into citizen and non-citizen science data#################


load("data/frog_occ_cleaned_final_20240701.Rda")


### Group data by data-providing institution 
records_by_institution <- frog_occ_cleaned_final %>% 
  group_by(dataResourceName) %>%
  summarise(records = n())

write.csv(records_by_institution, file = "Results/clean_records_by_institution_final_20240702.csv", row.names = FALSE)


### Identify citizen science projects from dataResource Name
unique(frog_occ_cleaned_final$dataResourceName)

# A vector of citizen science data
cs_data <- c("FrogID", "iNaturalist Australia", "ClimateWatch", 
             "Melbourne Water Frog Census", "ALA species sightings and OzAtlas", 
             "Gaia Guide","Earth Guardians Weekly Feed")

## Create a column for classifying records as citizen/non-citizen science data

frog_occ_cleaned_final <- frog_occ_cleaned_final %>%
  mutate(datatype = ifelse(dataResourceName %in% cs_data, "citizenscience", "non_citizenscience"))

save(frog_occ_cleaned_final, file = "frog_occ_cleaned_final_with_datatype_20240702.Rda")

summary(factor(frog_occ_cleaned_final$datatype))
# citizenscience     non_citizenscience  
#         394986                382486
# Group data by citizen vs non-citizen science data
records_by_species <- frog_occ_cleaned_final %>%
  group_by(species) %>%
  summarise(total_records = n(),
            citizen_science_records = sum(datatype == "citizenscience"),
            non_citizen_science_records = sum(datatype != "citizenscience"))

write.csv(records_by_species,  file = "Results/clean_records_by_species_csvsnoncs_final_20240702.csv", row.names = FALSE )


################################################################################
############Plot distribution of records for cleaned data#######################

### Load data and packages if needed

load("data/frog_occ_cleaned_final_with_datatype_20240702.Rda")

library(ggplot2)
library(ozmaps) 

### Retrieve Australia map data
aus <- ozmap_data()


### Plot distribution of occurrences for cleaned data

# Distribution of ALA vs FrogID
frog_distribution_cleaned <-ggplot() + 
  geom_sf(data = aus, fill = "#FBFBEF", size=NA) +
  geom_point(data = frog_occ_cleaned_final, 
             aes(x =decimalLongitude , y=decimalLatitude, colour = factor(dataprovider)),
             size =0.0001, shape =20) +
  coord_sf(expand=TRUE, xlim = c(110,175), ylim = c(-45,-10)) +
  labs(x="Longitude", y="Latitude") + 
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  theme(legend.position = c(0.9,0.9),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =14),
        legend.key = element_rect(colour = NA, fill = NA),
        #legend.key.size = unit(5,"cm"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.background = element_blank())

ggsave(plot = frog_distribution_cleaned, "Results/frog_distribution_cleaneddata_map_20240702.png", width = 12, height = 6, dpi = 1200)


# Distribution of Noncitizen vs citizen science data
frog_distribution_cleaned_2 <-ggplot() + 
  geom_sf(data = aus, fill = "#FBFBEF", size=NA) +
  geom_point(data = subset(frog_occ_cleaned_final, datatype != "citizenscience"), 
             aes(x =decimalLongitude, y=decimalLatitude, colour = "Non-citizen science data"), size = 0.00001) +
  geom_point(data = subset(frog_occ_cleaned_final, datatype == "citizenscience") , 
             aes(x =decimalLongitude, y=decimalLatitude, colour="Citizen science data"), size = 0.00001) +
  coord_sf(expand=TRUE, xlim = c(110,175), ylim = c(-45,-10)) +
  labs(x="Longitude", y="Latitude") + 
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  theme(legend.position = c(0.9,0.9),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =14),
        legend.key = element_rect(colour = NA, fill = NA),
        #legend.key.size = unit(5,"cm"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        panel.background = element_blank())

ggsave(plot = frog_distribution_cleaned_2, "Results/frog_distribution_cleaneddata_map_cs_noncs_20240702.png", width = 12, height = 6, dpi = 1200)

############## End of Data Processing: Part 2 ##################################
###################### End of data cleaning ####################################
################################################################################
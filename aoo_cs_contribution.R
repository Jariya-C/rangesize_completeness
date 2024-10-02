################################################################################
### AOO Estimates and Completeness without citizen science data
###
### Script for estimating observed AOO and AOO completeness without citizen science data
### Part of the methods for the manuscript:
### How well do we understand geographic range size?: 
### A case study of Australia’s frogs and citizen science project
###
### 
################################################################################

####  Load required packages and data objects
library(raster)
library(tidyverse)
library(sf)
library(iNEXT)
library(sp)
library(reshape2)
library(ozmaps)
library(viridis)
library(cowplot)
library(rstatix)

load("data/frog_occ_cleaned_final_with_datatype_20240702.Rda")

iucn_redlist <- read_csv("data/iucn_redlist_namesstandardised.csv")

noncs_data <- frog_occ_cleaned_final %>% 
  filter(datatype == "non_citizenscience")

noncs_data_sf <- st_as_sf(noncs_data, coords = c("decimalLongitude", "decimalLatitude"))

st_crs(noncs_data_sf) <- 4326


noncs_records <- noncs_data %>% dplyr::select(c(species,decimalLongitude,decimalLatitude))
noncs_records <- as.data.frame(noncs_records)

noncs_records_grouped <- noncs_records %>% group_by(species) %>% summarise(records = n())

### species list
unique_species <- unique(records$species) #203 spp

################################################################################
### Estimating AOO completeness with non-CS data using Chao 1 method ####

## covert frog_sp to spatial object
frogs_sp <- noncs_records

## Convert data frame to sp object
coordinates(frogs_sp) <-c("decimalLongitude","decimalLatitude")

## Set the coordinate reference system to EPSG:4326 and convert to EPGS:3577 for Australian Albers coordinate system
proj4string(frogs_sp) <- '+init=epsg:4326'
frogs_sp <- spTransform(frogs_sp, CRS('+init=epsg:3577'))
noncs_data_sf <- st_as_sf(frogs_sp)

ausmap <- st_transform(ozmaps::ozmap_country, 3577)

## Create an empty list to store results
chao_noncs_results_list <- list()
cell_sizes <- c(1, 2, 4, 8, 16, 32, 64, 128)


## Loop through each cell size
pb_aoo_comp <- txtProgressBar(max = length(cell_sizes), style = 3)
for (i in 1:length(cell_sizes)) {
  # List of cell sizes in kilometers
  cell_size_km <- cell_sizes[i]
  # Convert cell size to meters
  cell_size_m <- cell_size_km * 1000
  
  extent_ausmap <- extent(ausmap)
  
  # Calculate the number of rows and columns based on the cell size
  num_rows <- ceiling((extent_ausmap@ymax - extent_ausmap@ymin) / cell_size_m)
  num_cols <- ceiling((extent_ausmap@xmax - extent_ausmap@xmin) / cell_size_m)
  
  # Create a raster with the specified cell size
  raster_ausmap <- raster(extent(extent_ausmap), nrow = num_rows, ncol = num_cols, 
                          crs = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs", res = cell_size_m)
  
  ausmap_raster <- rasterize(ausmap, raster_ausmap)
  
  # Stack the original ausmap_raster with a unique values raster
  unique_values_raster <- raster(ausmap_raster)
  unique_values_raster[] <- 1:ncell(unique_values_raster) # create second layer which represents values of cell ID
  aus_stacked_raster <- stack(ausmap_raster, unique_values_raster) 
  
  # Extract raster values to the spatial data
  frogs_joined <- raster::extract(aus_stacked_raster, noncs_data_sf, df=TRUE, sp=TRUE)
  allfrogs_data <- as.data.frame(frogs_joined)
  allfrogs_data <- allfrogs_data[, c("layer.2", "species", "coords.x1", "coords.x2")]
  allfrogs_data$var <- 1
  
  # Reshape the data
  frog_aoo_completeness <- dcast(allfrogs_data, layer.2 ~ species, value.var = "var")
  frog_aoo_completeness[is.na(frog_aoo_completeness)] <- 0
  
  # Calculate ChaoRichness
  aoo_completeness <- ChaoRichness(frog_aoo_completeness[, -1], datatype = "abundance")
  
  aoo_completeness$species <-row.names(aoo_completeness) 
  
  aoo_completeness$completeness <- aoo_completeness$Observed/aoo_completeness$Estimator
  aoo_completeness$aoo_area <- aoo_completeness$Observed*cell_size_km^2
  aoo_completeness <- left_join(aoo_completeness, noncs_records_grouped, by = c("species"))
  aoo_completeness <-left_join(aoo_completeness, iucn_redlist, by = c("species"))
  
  # Add results to the list
  chao_noncs_results_list[[paste("CellSize", cell_size_km, "km")]] <- aoo_completeness
  
  
  setTxtProgressBar(pb_aoo_comp, value = i)
  
}

close(pb_aoo_comp)

save(chao_noncs_results_list, file = "result/aoo_chaoestimate_results_list_noncitizenscience_final.Rda")


## Convert non_citizen_chao_results_list to data frame
non_citizen_chao_results_df <- do.call(rbind, chao_noncs_results_list)
non_citizen_chao_results_df$Cellsize <- row.names(non_citizen_chao_results_df)

non_citizen_chao_results_df$Cellsize[1:203] <- 1
non_citizen_chao_results_df$Cellsize[204:406] <- 2
non_citizen_chao_results_df$Cellsize[407:609] <- 4
non_citizen_chao_results_df$Cellsize[610:812] <- 8
non_citizen_chao_results_df$Cellsize[813:1015] <- 16
non_citizen_chao_results_df$Cellsize[1016:1218] <- 32
non_citizen_chao_results_df$Cellsize[1219:1421] <- 64
non_citizen_chao_results_df$Cellsize[1422:1624] <- 128

non_citizen_chao_results_df$Cellsize <- as.character(non_citizen_chao_results_df$Cellsize)

chao_results_df_noncs <- non_citizen_chao_results_df

save(chao_results_df_noncs, file = "result/aoo_chaoestimate_results_df_noncitizenscience_final.Rda")

################################################################################
### Contribution of Citizen science to AOO estimates and completeness
### @ 2km spatial resolution

# Load the data if needed
#load("result/aoo_chaoestimate_results_df_noncitizenscience_final.Rda")
#load("result/aoo_estimates_completeness_allcellsizes_withredlist.Rda")


# Bind two dataframes together and save 

noncs_chao2km <- chao_results_df_noncs %>% 
  filter(Cellsize==2)

noncs_chao2km <- noncs_chao2km[,c(1:9)]
colnames(noncs_chao2km) <- c("noncs_observed", "noncs_estimated", "noncs_est.se", "noncs_95%lower", 
                             "noncs_95%upper", "species", "noncs_completeness", "noncs_aoo", "noncs_recs")

alldata_chao2km <- chao_results_df_withredlist%>% filter(Cellsize ==2)

chao2km_cs_contribution_allspp <- left_join(alldata_chao2km, 
                                     noncs_chao2km, by = "species")

# For species where there is no non-citizen science records, replace aoo_area with 0

chao2km_cs_contribution_allspp$noncs_aoo[is.na(chao2km_cs_contribution_allspp$noncs_aoo)] <- 0

# Calculate citizen science contribution to AOO estimate
## AOO C = (AOO_alldata - AOO_ncsdata)/Aoo_alldata

chao2km_cs_contribution_allspp$cs_to_AOOestimate <- 
  (chao2km_cs_contribution_allspp$aoo_area - chao2km_cs_contribution_allspp$noncs_aoo)/
  chao2km_cs_contribution_allspp$aoo_area

# Calculate citizen science contribution to AOO completeness
## AOO C = (AOO_alldata - AOO_ncsdata)/Aoo_predicted_with_alldata

chao2km_cs_contribution_allspp$cs_to_AOOcompleteness <- 
  (chao2km_cs_contribution_allspp$aoo_area - chao2km_cs_contribution_allspp$noncs_aoo)/
  (chao2km_cs_contribution_allspp$Estimator * 4)



save(chao2km_cs_contribution_allspp, file = "result/aoo_estimates_and_completeness_with_cs_contribution_2km_allspp.Rda")


# Filter for species with at least 100 records
chao2km_cs_contribution_filtered <- subset(chao2km_cs_contribution_allspp, records >= 100)

mean(chao2km_cs_contribution_filtered$cs_to_AOOestimate)
#[1]  0.2829664
sd(chao2km_cs_contribution_filtered$cs_to_AOOestimate)
#[1]  0.2973092
mean(chao2km_cs_contribution_filtered$cs_to_AOOcompleteness)
# 0.1587121
sd(chao2km_cs_contribution_filtered$cs_to_AOOcompleteness)
# 0.1871617

save(chao2km_cs_contribution_filtered, file = "result/contribution_of_cs_to_aoo_estimates_and_completeness_138spp_final.Rda")



hist_cs_aooest <- ggplot(chao2km_cs_contribution_filtered, aes(x=cs_to_AOOestimate)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill ="grey", lwd =0.2) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(0,35),  breaks = seq(0, 30, by = 5)) +
  geom_vline(xintercept = 0.283, linetype = "longdash") +
  annotate("text", x = 0.33, y = 30, label = "0.283", size = 5) +
  labs(x = "Proportion of citizen science contribution to AOO estimate", y = "Number of species")+
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"))


hist_cs_aoocomp <- ggplot(chao2km_cs_contribution_filtered, aes(x=cs_to_AOOcompleteness)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill ="grey", lwd =0.2) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(0,35),  breaks = seq(0, 30, by = 5)) +
  geom_vline(xintercept = 0.159, linetype = "longdash") +
  annotate("text", x = 0.20, y = 30, label = "0.159", size = 5) +
  labs(x = "Proportion of citizen science contribution to AOO completeness", y = "Number of species")+
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"))

hist_cs_aoo_combined <- plot_grid(hist_cs_aooest, hist_cs_aoocomp, nrow = 2, ncol = 1,labels=c("(a)", "(b)"), label_size = 16)

ggsave(plot = hist_cs_aoo_combined, "result/figS12_cs_contribution_to_aoo_estimate_completeness_histograms_new.png", width = 12, height = 15, dpi = 1200)

##### End of script for  CS contribution to AOO estimates and completeness #####
################################################################################
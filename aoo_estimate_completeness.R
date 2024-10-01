################################################################################
### AOO Estimates and Completeness
###
### Script for estimating area of occupancy (AOO) and completeness of AOO
### Part of the methods for the manuscript:
### How well do we understand geographic range size?: 
### A case study of Australia’s frogs and citizen science
### Note: This script estimates AOO and completeness AOO for all data
### 
################################################################################

### Load required packages and data objects
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
frogdata <- frog_occ_cleaned_final

frogdata_sf <- st_as_sf(frogdata, coords = c("decimalLongitude", "decimalLatitude"))

st_crs(frogdata_sf) <- 4326


iucn_redlist <- read_csv("data/IUCN_Redlist_20240521_final/assessments.csv")

records <- frogdata %>% dplyr::select(c(species,decimalLongitude,decimalLatitude))
records <- as.data.frame(records)

records_grouped <- records %>% group_by(species) %>% summarise(records = n())

################################################################################
### 1. Calculating range area:area of occupancy (AOO)
### using different cell sizes
### The code in this section was adopted from:
### Marsh, C.J., Syfert, M.M., Aletrari, E. et al (2023). 
### The effect of sampling effort and methodology on range size estimates of poorly-recorded species for IUCN Red List assessments. 
### Biodiversity and Conservation 32, 1105–1123. https://doi.org/10.1007/s10531-023-02543-9

unique_species <- unique(records$species) 

### Data frame for storing results at different cell widths
aoo <- data.frame(species = unique_species,
                  aoo1 = NA,
                  aoo2 = NA,
                  aoo4 = NA,
                  aoo8 = NA,
                  aoo16 = NA,
                  aoo32 = NA,
                  aoo64 = NA,
                  aoo128 = NA)

### Loop through species

pb_aoo <- txtProgressBar(max = length(unique_species), style = 3)
for(sp in 1:length(unique_species)) {
  # Subset points for selected species and remove duplicates
  sp_coords <- records[records$species == as.character(unique_species[sp]), c("decimalLongitude", "decimalLatitude")]
  sp_coords <- sp_coords[!duplicated(sp_coords), ]
  
  if(nrow(sp_coords) > 0) {
    sp_coords <- SpatialPoints(sp_coords, proj4string = CRS("+proj=longlat +ellps=WGS84 +degrees=TRUE"))
    sp_coords <- spTransform(sp_coords, CRS("+proj=cea +datum=WGS84 +units=km"))
    
    
    # Loop through generating grids at the different cell widths (km)
    cellWidths <- c(1, 2, 4, 8, 16, 32, 64, 128)
    
    for(cellWidth in cellWidths) {
      # If only a single point then AOO is area of single cell
      if(length(sp_coords) == 1) {
        aoo[sp, paste0("aoo", cellWidth)] <- cellWidth ^ 2
      }
      
      # If more than one point then rasterise and sum area of occupied cells
      if(length(sp_coords)) {
        r <- raster(xmn =   floor(extent(sp_coords)[1] / cellWidth) * cellWidth,
                    xmx = ceiling(extent(sp_coords)[2] / cellWidth) * cellWidth,
                    ymn =   floor(extent(sp_coords)[3] / cellWidth) * cellWidth,
                    ymx = ceiling(extent(sp_coords)[4] / cellWidth) * cellWidth,
                    resolution = cellWidth)
        aooRaster <- rasterize(sp_coords, r) # rasterise sp_coords into r raster so each occurrence point is assigned to a raster cell
        vals <- values(aooRaster)
        vals[vals > 1] <- 1 # set cell values greater than one to be one 
        vals[is.na(vals)] <- 0 # set NA cell values to be zero
        aoo[sp, paste0("aoo", cellWidth)] <- sum(vals) * (cellWidth ^ 2) # calculate AOO by summing the grid cells x area of each cell
      }
    }
  }
  setTxtProgressBar(pb_aoo, value = sp)
}
close(pb_aoo)

aoo_with_records <- left_join(aoo, records_grouped, by = c("species"))

### Save results as CSV
write.csv(aoo_with_records, "result/aoo_gridoverlay_various_cellsizes.csv", quote = FALSE, row.names = FALSE)

################################################################################
### Obtain IUCN Redlist status for each frog species
redlist <-iucn_redlist[,c("scientificName","redlistCategory", "redlistCriteria", "yearPublished")]

names_check<- setdiff(aoo$species, redlist$scientificName)

### Standardise names in redlist to match current species list
redlist$scientificName <- gsub("Crinia nimbus", "Crinia nimba",
                          gsub("Cyclorana vagitus","Cyclorana vagita",
                          gsub("Litoria burrowsi", "Litoria burrowsae",
                          gsub("Litoria lesueurii","Litoria lesueuri",
                          gsub("Neobatrachus sudelli","Neobatrachus sudellae",
                          gsub("Lechriodus fletcheri", "Platyplectrum fletcheri",
                          gsub("Geocrinia vitellina", "Anstisia vitellina",
                          gsub("Geocrinia alba", "Anstisia alba",
                          gsub("Geocrinia rosea", "Anstisia rosea",
                          gsub("Geocrinia lutea", "Anstisia lutea", 
                         redlist$scientificName))))))))))

names_check<- setdiff(aoo$species, redlist$scientificName)
redlist <- as.data.frame(redlist)
colnames(redlist)[colnames(redlist) == "scientificName"] <- "species"

write.csv(redlist, "data/iucn_redlist_namesstandardised.csv", quote = TRUE, row.names = FALSE)

aoo_with_records_and_redlist <-left_join(aoo_with_records, redlist, by = c("species"))

save(aoo_with_records_and_redlist, file = "result/aoo_various_cellsizes_with_records_and_redlist.Rda")

write.csv(aoo_with_records_and_redlist, "result/aoo_various_cellsize_with_records_and_redlist.csv", quote = TRUE, row.names = FALSE)

################################################################################
### Further analysis

aoo_with_records_and_redlist$threats_by_aoo <- 
  ifelse(aoo_with_records_and_redlist$aoo2 < 10, "Critically Endangered",
  ifelse(aoo_with_records_and_redlist$aoo2 >= 10 & 
  aoo_with_records_and_redlist$aoo2 < 500, "Endangered",
  ifelse(aoo_with_records_and_redlist$aoo2 >= 500 & 
  aoo_with_records_and_redlist$aoo2 < 2000, "Vulnerable",
                                          "Not threatened")))

save(aoo_with_records_and_redlist, file = "result/aoo_estimate_with_inferred_redlist_category.Rda")

### Sumnmarise observed threat status based on AOO calculation at 2 x 2 km cell size
observed_threatstatus_subset_atleast100records <- aoo_with_records_and_redlist %>% 
  filter(records >= 100) %>%
  group_by(threats_by_aoo) %>% 
  summarise(total.spp = n())

# A tibble: 3 × 2
# threats_by_aoo total.spp
#<chr>                 <int>
#1 Endangered            36
#2 Not threatened        52
#3 Vulnerable            50

# Actual IUCN threat status
actual_threatstatus_subset_atleast100records <- aoo_with_records_and_redlist %>% 
  filter(records >= 100) %>%
  group_by(redlistCategory) %>% 
  summarise(total.spp = n())

#actual_threatstatus_subset_atleast100records
# A tibble: 5 × 2
#redlistCategory       total.spp
#<chr>                     <int>
#1 Critically Endangered         1
#2 Endangered                   11
#3 Least Concern               112
#4 Near Threatened               6
#5 Vulnerable                    8

################################################################################
#### 3. Estimating range size completeness for AOO with Chao 1 ####


### Load data if required
#load("data/frog_occ_cleaned_final_with_datatype_20240702.Rda")
#frogdata <- frog_occ_cleaned_final
#records <- frogdata %>% dplyr::select(c(species,decimalLongitude,decimalLatitude))


### Covert frog_sp to spatial object
frogs_sp <- records

coordinates(frogs_sp) <-c("decimalLongitude","decimalLatitude")

# Set the coordinate reference system to EPSG:4326 and convert to EPGS:3577 for Australian Albers coordinate system
proj4string(frogs_sp) <- '+init=epsg:4326'
frogs_sp <- spTransform(frogs_sp, CRS('+init=epsg:3577'))
frogs_sf <- st_as_sf(frogs_sp)

ausmap <- st_transform(ozmaps::ozmap_country, 3577)

## Create an empty list to store results
chao_results_list <- list()
cell_sizes <- c(1, 2, 4, 8, 16, 32, 64, 128)


### WARNING - This loop takes time to run
## Because it loops through each cell size
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
 
  # Rasterize ausmap_transformed
  ausmap_raster <- rasterize(ausmap, raster_ausmap)
  
  # Stack the original ausmap_raster with a unique values raster
  unique_values_raster <- raster(ausmap_raster)
  # Create second layer which represents values of cell ID
  unique_values_raster[] <- 1:ncell(unique_values_raster)
  aus_stacked_raster <- stack(ausmap_raster, unique_values_raster) 
  
  # Extract raster values to the spatial data
  frogs_joined <- raster::extract(aus_stacked_raster, frogs_sf, df=TRUE, sp=TRUE)
  allfrogs_data <- as.data.frame(frogs_joined)
  
  allfrogs_data <- allfrogs_data[, c("layer.2", "species", "coords.x1", "coords.x2")]
  allfrogs_data$var <- 1
  
  # Reshape the data
  frog_aoo_completeness <- dcast(allfrogs_data, layer.2 ~ species, value.var = "var")
  frog_aoo_completeness[is.na(frog_aoo_completeness)] <- 0
  
  # Calculate expected cellsize with ChaoRichness
  aoo_completeness <- ChaoRichness(frog_aoo_completeness[, -1], datatype = "abundance")
  
  aoo_completeness$species <-row.names(aoo_completeness) 
  
  # Calculate AOO completeness as no. of grid cells observed/ no. of predicted grid cells
  aoo_completeness$completeness <- aoo_completeness$Observed/aoo_completeness$Estimator
  aoo_completeness$aoo_area <- aoo_completeness$Observed*cell_size_km^2
  aoo_completeness <- left_join(aoo_completeness, records_grouped, by = c("species"))
  aoo_completeness <-left_join(aoo_completeness, redlist, by = c("species"))
  
  # Add results to the list
  chao_results_list[[paste("CellSize", cell_size_km, "km")]] <- aoo_completeness
  
  setTxtProgressBar(pb_aoo_comp, value = i)
  
}
close(pb_aoo_comp)

chao_results_list_all_final <- chao_results_list

save(chao_results_list_all_final, file = "result/aoo_chaoestimate_results_list_final.Rda")


## Convert chao_results_list to data frame

chao_results_df <- do.call(rbind, chao_results_list_all_final)
chao_results_df$Cellsize <- row.names(chao_results_df)

chao_results_df$Cellsize[1:227] <- 1
chao_results_df$Cellsize[228:454] <- 2
chao_results_df$Cellsize[455:681] <- 4
chao_results_df$Cellsize[682:908] <- 8
chao_results_df$Cellsize[909:1135] <- 16
chao_results_df$Cellsize[1136:1362] <- 32
chao_results_df$Cellsize[1363:1589] <- 64
chao_results_df$Cellsize[1590:1816] <- 128

chao_results_df$Cellsize <- as.character(chao_results_df$Cellsize)

save(chao_results_df, file = "result/chao_rangesize_estimates_and_completeness_final.Rda")

### Box plot
#install.packages("hrbrthemes")

### Convert Cellsize to a factor with ordered levels
chao_results_df$Cellsize <- factor(chao_results_df$Cellsize, levels = unique(chao_results_df$Cellsize))

chao_aoo_boxplot <-ggplot(chao_results_df, aes(x = Cellsize, y = completeness, fill = Cellsize)) +
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  stat_summary(fun=mean, geom="point", shape = 17, size = 3, colour = "black")+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  xlab("Cell resolution (km)") + ylab("Completeness of AOO estimates") +
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(),
        axis.ticks.y = element_line(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size =20))

ggsave(plot = chao_aoo_boxplot, "result/chao_aoo_boxplot_allresults.png", dpi = 1200, width = 12, height = 6)



### Plot number of records vs AOO completeness

df <- na.omit(cbind(chao_results_df$completeness,chao_results_df$records))
y <- df[,1]
x <- log10(df[,2])
lo <- loess(y~x)
xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
out = predict(lo,xl)
infl <- c(FALSE, diff(diff(out)>0)!=0)
plot(x,y,type="p")
lines(xl, out, col='red', lwd=2)
points(xl[infl ], out[infl ], col="blue", lwd = 5)

### Determine the start of monotonic relationship
10^xl[infl ]
#106.0092

comp_vs_n <- ggplot(data = chao_results_df, aes(x = records, y = completeness)) + 
  geom_point(aes(colour = factor(Cellsize))) +
  scale_colour_viridis_d() + 
  scale_shape_manual (values = c(3,4,5,6,7,8,9,10)) +
  scale_x_continuous(trans = "log10", breaks =c(1,10,100,1000,10000, 100000), labels = c("1","10","100","1000", "10000", "100000")) + 
  geom_smooth(method = "loess") + 
  annotate("point", x =107.4324, y = 0.6604759, shape = 16, size = 4, colour = "black") +
  geom_vline(xintercept = 107, linetype = "dashed" ) +
  labs(x = "Number of records", y = "Completeness of AOO estimates", colour = "Cell size in km")+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(size=18, colour = "black")) + 
  theme(legend.position = "bottom",
        legend.text = element_text(size = 16))

ggsave(plot=comp_vs_n, "result/completeness_minimumcutoff_100record.png", width = 8, height = 8, dpi =1200)  

# Combine boxplot and point graph
aoo_completeness_plots <- plot_grid(chao_aoo_boxplot, comp_vs_n, nrow = 2, ncol = 1,labels=c("(a)", "(b)"), label_size = 16)

ggsave(plot = aoo_completeness_plots, "result/figS1_aoo_and_completeness_boxxplot_and_pointgraph.png", width = 12, height = 15, dpi = 1200)

### Monotonic relationship between AOO completeness and no of records was observed for species with >= 100 records
### AOO completeness increases with cell size
### Filter and plot species with at least 100 records for 2 x 2 km cell size

chao_results_df_withredlist <- chao_results_df


chao_results_df_withredlist$redlistID <-chao_results_df_withredlist$redlistCategory 
chao_results_df_withredlist$redlistID <- gsub(" ", "_", chao_results_df_withredlist$redlistID)

unique(chao_results_df_withredlist$redlistID)

chao_results_df_withredlist$redlistID <- gsub("Critically_Endangered","1",
                                   gsub("Endangered", "2",
                                   gsub("Vulnerable", "3",
                                   gsub("Near_Threatened", "4",
                                   gsub("Least_Concern", "5", 
                                   gsub("Data_Deficient", "6",
                                   chao_results_df_withredlist$redlistID)))))) 


chao_results_df_withredlist$redlistID <- gsub("Critically_2", "1",chao_results_df_withredlist$redlistID)

chao_results_df_withredlist$threatstatus <-chao_results_df_withredlist$redlistID
chao_results_df_withredlist$threatstatus <-  gsub("1", "Threatened",
                                       gsub("2","Threatened",
                                       gsub("3","Threatened",
                                       gsub("4", "Not threatened",
                                       gsub("5", "Not threatened", 
                                       gsub("6","Not threatened",
                                       chao_results_df_withredlist$threatstatus)))))) 


save(chao_results_df_withredlist, file = "result/aoo_estimates_completeness_allcellsizes_withredlist.Rda")

filtered_chao_results_2km <- chao_results_df_withredlist %>% filter(Cellsize == 2)
save(filtered_chao_results_2km, file= "result/filtered_chao_results_2km_alldatatype_final.Rda")

### Boxplot aoo and threatstatus

chao_results_2km_subset <-  chao_results_df_withredlist %>% filter(Cellsize == 2 & records >= 100)

# Summary statistics

chao_results_2km_subset %>%  
  group_by(threatstatus) %>%
  get_summary_stats(completeness)

# A tibble: 2 × 14
#threatstatus     variable         n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
#<chr>             <fct>        <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#1 Not threatened completeness   118 0.317 0.841  0.507 0.442 0.599 0.157 0.118 0.523 0.106 0.01  0.019
#2 Threatened     completeness    20 0.284 0.91   0.594 0.52  0.673 0.152 0.117 0.604 0.159 0.036 0.075


chao_results_2km_subset %>% 
  wilcox_test(completeness ~ threatstatus) %>%
  add_significance()

# A tibble: 1 × 8
#.y.                group1         group2        n1    n2 statistic      p p.signif
#<chr>              <chr>          <chr>      <int> <int>     <dbl>  <dbl> <chr>   
#  1 completeness Not threatened Threatened   118    20       794 0.0197 *    


chao_results_2km_subset %>% rstatix::wilcox_effsize(completeness~threatstatus)
# A tibble: 1 × 7
#   .y.          group1         group2     effsize    n1    n2 magnitude
# * <chr>        <chr>          <chr>        <dbl> <int> <int> <ord>    
#  1 completeness Not threatened Threatened   0.199   118    20 small   


##Plot historgram of AOO completeness for all 138 spp @ 2km 

# Get mean completeness value 
summary(chao_results_2km_subset$completeness)
# Mean = 0.54

plothist <- ggplot(chao_results_2km_subset, aes(x=completeness)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill ="grey", lwd =0.2) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(0,35),  breaks = seq(0, 30, by = 5)) +
  geom_vline(xintercept = 0.535, linetype = "longdash") +
  annotate("text", x = 0.57, y = 30, label = "0.54", size = 5) +
  labs(x = "AOO completeness", y = "Number of species")+
  theme_classic(base_size = 16) +
  theme(legend.text = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 14, colour = "black")) 

ggsave(plot = plothist, "result/aoo_completeness_hist_138spp.png", width = 5.5, height = 4, dpi = 1200)


### histogram of AOO for threatened and non-threatened spp

plothist2 <- ggplot(chao_results_2km_subset, aes(x = completeness, fill = threatstatus)) +
  geom_histogram(binwidth = 0.05, 
                 lwd = 0.2, 
                 colour = "black",
                 alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5)) +
  labs(x = "AOO completeness", y = "Number of species") + 
  scale_fill_manual(values = c("Threatened" = "#F0E442", "Not threatened" = "#009E73"), 
                    name = "") +
  theme_classic(base_size = 16) +
  theme(legend.position = "top",
        legend.text = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"))


hist_combined <- plot_grid(plothist, plothist2, nrow = 2, ncol = 1,labels=c("(a)", "(b)"), label_size = 16)

ggsave(plot = hist_combined, "result/figS7_aoo_completeness_histograms.png", width = 12, height = 15, dpi = 1200)

###### End of script for AOO estimates and completeness with all data ##########
################################################################################
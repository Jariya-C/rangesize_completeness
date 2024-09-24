################################################################################
### Data trends and biases
###
### Script for assessing data trends and biases
### Part of the methods for the manuscript:
### How well do we understand species’ geographic range size?: 
### A case study of Australia’s frogs
###
### Jariya Chanachai (jariya.chanachai@hdr.mq.edu.au)
################################################################################
###
###
###
################################################################################
### Load required packages and data objects

#install.packages("remotes")
#remotes::install_github("robboyd/occAssess")

library(tidyverse)
library(ozmaps)
library(ggspatial)
library(sf)
library(terra)
library(viridis)
library(spatialEco) #for NNI test
library(cowplot)
library(rstatix)
library(occAssess)
library(raster)

load("data/frog_occ_cleaned_final_with_datatype_20240702.Rda")

# Australia
aus <- ozmap_data(data = "states")

# Australia bioregion shapefiles
ibra_shp <- st_read("data/IBRA7_regions/ibra7_regions.shp")

# Select relevant columns 
allfrogdata <- as.data.frame(frog_occ_cleaned_final[,c(1:15,39,53:56)])

###############################################################################
### Distribution of records and distribution of sampling efforts across Australia


state_centroids <- aus %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  as.data.frame() %>%
  mutate(state_name = aus$NAME) 


state_centroids <- state_centroids[c(1:8),]  # Remove Other Territories
# Use acronyms for the states
state_centroids$state_name <- gsub("New South Wales", "NSW", 
                              gsub("Victoria", "VIC", 
                              gsub("Queensland", "QLD",
                              gsub("South Australia", "SA",
                              gsub("Western Australia", "WA",
                              gsub("Tasmania", "TAS", 
                              gsub("Northern Territory", "NT", 
                              gsub("Australian Capital Territory", "ACT", 
                                   state_centroids$state_name))))))))


plot_distribution <- ggplot() + 
  geom_sf(data = aus, fill = "#FBFBEF", size=NA) +
  geom_point(data = subset(allfrogdata, datatype == "non_citizenscience"), 
             aes(x =decimalLongitude, y=decimalLatitude, 
                 colour = "Non-citizen science data"), 
             size = 0.00001) +
  geom_point(data = subset(allfrogdata, datatype == "citizenscience"), 
             aes(x =decimalLongitude, y=decimalLatitude, 
                 colour="Citizen science data"),
             size = 0.00001) +
  geom_text(data = state_centroids, 
            aes(x = X, y = Y, label = state_name), 
            size = 4, color = "black", fontface = "bold") +  # Add state names with geom_text
  annotation_scale()+
  coord_sf(expand = TRUE,
           xlim =  c(110, 155), 
           ylim = c(-45, -10)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic()+
  theme(legend.position = c(1,0.8),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =15),
        legend.background = element_blank(),
        legend.key.size = unit(0.6,"cm")
  )

ggsave(plot = plot_distribution, "result/fig4_distribution_map_withstatenames.png", dpi = 1200, width = 10, height = 6)

###
### Spatial join occurrences with bioregion and calculate record density

ibra_dataframe <- as.data.frame(ibra_shp)

ibra_valid <- st_make_valid(ibra_shp)

# Covert frogdata to sf object
frog_records_sf <- st_as_sf(allfrogdata, coords = c("decimalLongitude","decimalLatitude"), crs=st_crs(ibra_shp))

# Spatial join records and bioregions
ibra_occ_join <- st_join(frog_records_sf, ibra_valid)

ibra_rec_df <- as.data.frame(ibra_occ_join)

ibra_density <- ibra_rec_df %>% 
  group_by(REG_NAME_7) %>%
  summarise(records =n(), richness = length(unique(species)))

# Calculate sampling effort as number of records/km2
ibra_density <- merge(ibra_dataframe, ibra_density, by = "REG_NAME_7", all.x =TRUE) %>%
  mutate(samplingeffort = records/SQ_KM) %>% 
  mutate(samplingeffort_log10 = log10(records/SQ_KM))

ibra_density_sf<- st_as_sf(ibra_density)

### Plot sampling effort map
plot_sampling_effort <- ggplot(ibra_density_sf, fill= NA, size =0.5) + 
  geom_sf(data = subset(ibra_density_sf, records >=1), 
          aes(fill=(samplingeffort_log10)), 
          colour ="grey50", lwd =0.001) +
  scale_fill_viridis(256, option = "D",
                     name = bquote("Sampling effort (records/km"^2*")"),
                     na.value = NA,
                     #type = "seq",
                     #palette = "YlOrBr",
                     direction = 1,
                     breaks = c(-5, -4, -3, -2, -1, 0),  # Specifying breaks matching the labels
                     labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1","1"),
                     guide = guide_colorsteps(direction = "horizontal",
                                              label.position = "bottom",
                                              title.position = "left")) +
  annotation_scale() +
  #annotation_scale(location = "bl",  pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm")) +
  coord_sf(xlim = c(110, 155), ylim = c(-45, -10)) +
  labs(x = "decimalLongitude", y = "decimalLatitude") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.key.height = unit(2, 'mm'),
        legend.key.width = unit(20, 'mm'))

ggsave(plot = plot_sampling_effort, "result/figS3_samplingeffort_map.png", width = 10, height = 6, dpi = 1200)

################################################################################
### Test for spatial autocorrelation with Nearest Neighbour
### Assess spatial bias

cs_records_df <- as.data.frame(frog_occ_cleaned_final) %>% filter(datatype == "citizenscience")
noncs_records_df <- as.data.frame(frog_occ_cleaned_final) %>% filter(datatype != "citizenscience")

cs_records_sf <- st_as_sf(cs_records_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 3577)
noncs_records_sf <- st_as_sf(noncs_records_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 3577)

frogs_sf <- st_as_sf(allfrogdata, coords = c("decimalLongitude", "decimalLatitude"), crs = 3577)

records_by_species <- allfrogdata %>% group_by(species) %>% 
  summarise(totalrecords =n(), 
            cs_records = sum(datatype == "citizenscience"),
            noncs_records = sum(datatype == "non_citizenscience"))

### Calculate nearest-neighbour index

# Create list to store results
alldata_nni_list <- list()
csdata_nni_list <- list()
noncsdata_nni_list <- list()

species <- sort(unique(frogs_sf$species))

# Set progression bar to track progress
pb <- txtProgressBar(max = length(species), style = 3)

### Create a loop to calculate NNI for each species
for(i in 1:length(species)) {
  # Subset the data for the current species
  species_data_cs <- cs_records_sf[cs_records_sf$species == species[i], ]
  species_data_noncs <- noncs_records_sf[noncs_records_sf$species == species[i],]
  species_data_allrecords <- frogs_sf[frogs_sf$species == species[i], ]
  
  if (nrow(species_data_cs) > 3) {
    # Calculate NNI index for just citizen science records for species [i]
    cs_nni <- nni(species_data_cs, win = "extent")
    csdata_nni_list[[species[i]]] <-cs_nni
  }
  else {
    # If there is no data for the species, print a message
    cat("No data found for citizen science data of species:", species[i], "\n")
    csdata_nni_list[[species[i]]] <-NA 
  }
  if(nrow(species_data_noncs) > 3) {
    # Calculate NNI index for just non-citizen science records for species [i]
    noncs_nni <- nni(species_data_noncs, win = "extent")
    noncsdata_nni_list[[species[i]]] <- noncs_nni
  }
  else {
    # If there is no data for the species, print a message
    cat("No data found for non-citizen science data of species:", species[i], "\n")
    noncsdata_nni_list[[species[i]]] <- NA
  } 
  if(nrow(species_data_allrecords) > 3) {
    # Calculate NNI index for all records for species [i]
    allrecords_nni <- nni(species_data_allrecords, win = "extent")
    alldata_nni_list[[species[i]]] <- allrecords_nni
  }
  else {
    # If there is no data for the species, print a message
    cat("No data found for species:", species[i], "\n")
    alldata_nni_list[[species[i]]] <- NA
  }  
  
  setTxtProgressBar(pb, value = i) 
}

cs_nni_df <- do.call(rbind.data.frame, csdata_nni_list)
cs_nni_df$species <- rownames(cs_nni_df)
cs_nni_df$datatype <- "citizenscience"
noncs_nni_df <- do.call(rbind.data.frame, noncsdata_nni_list)
noncs_nni_df$species <- rownames(noncs_nni_df)
noncs_nni_df$datatype <- "non_citizenscience"
alldata_nni_df <- do.call(rbind.data.frame, alldata_nni_list)
alldata_nni_df$species <- rownames(alldata_nni_df)
alldata_nni_df$datatype <- "alldatatype"


combined_data_nni_df <- rbind(cs_nni_df, noncs_nni_df, alldata_nni_df)
save(combined_data_nni_df, file = "results/combined_data_nni_df_final_20240702.Rda")

combined_data_nni_df_with_records <- left_join(records_by_species, combined_data_nni_df, by = "species")

save(combined_data_nni_df_with_records, file = "Results/combined_data_nni_df_with_records_final_20240702.Rda")
###
### Create boxplot and
### Test for significant difference between citizen and non-citizen science data


boxplot_nni <- ggplot(data = combined_data_nni_df_with_records, 
                      aes(x =  reorder(datatype,NNI), 
                          y = NNI, fill = datatype)) +
              geom_boxplot() +
              scale_fill_viridis(discrete = TRUE, alpha=0.6) +
              scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
              xlab("") + ylab("Nearest Neighbour Index (NNI)") +
              theme(panel.background = element_blank(),
                    legend.position = "none",
                    axis.line = element_line(),
                    axis.ticks.y = element_line(),
                    axis.text = element_text(size = 12),
                    axis.title = element_text(size =12))

ggsave(plot=boxplot_nni, "result/figS2_nni_boxplot.png", dpi = 1200,  width = 12, height = 8)


### Plot histogram of frequency of NNI to visualise the distribution.

histogram_plot_cs <-ggplot(subset(combined_data_nni_df_with_records, datatype == "citizenscience"), 
                           aes(x = NNI)) +
  geom_histogram(binwidth = 0.1, color = "black", fill ="grey", lwd =0.2) +
  labs(x = "Nearest neighbour index (NNI)", y = "Number of species")+
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(size=12, colour = "black")) +
  (annotate("text", label= "Citizen science data", x = 1, y= 60) )

ggsave(plot = histogram_plot_cs, "result/figS2_nni_histogram_csdata.png", width = 8, height = 5.5, dpi = 600)


histogram_plot_noncs <-ggplot(subset(combined_data_nni_df_with_records, datatype == "non_citizenscience"), 
                              aes(x = NNI)) +
  geom_histogram(binwidth = 0.1, color = "black", fill ="grey", lwd =0.2) +
  labs(x = "Nearest neighbour index (NNI)", y = "Number of species")+
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(size=12, colour = "black")) +
  (annotate("text", label= "Non citizen science data", x = 1, y= 60) )

ggsave(plot = histogram_plot_noncs, "result/figS2_nni_histogram_noncsdata.png", width = 8, height = 5.5, dpi = 600)

### Combine boxplot and histrograms into one figure
nni_hist <- plot_grid(histogram_plot_cs, histogram_plot_noncs, labels = c("(a)","(b)"), label_size = 12)
nni_boxplot_and_hist <- plot_grid(nni_hist, boxplot_nni, labels = c("", "(c)"), label_size = 12, ncol = 1 )

ggsave(plot = nni_boxplot_and_hist, "result/figS2_combined_nni_histogram_and_boxplot.png", width = 8, height = 8, dpi = 1200)

### Statistical analysis of NNI

library(rstatix)


# Obtain summary statistics by group

filtered_combined_nni <- combined_data_nni_df_with_records %>%
  filter(is.finite(NNI))

filtered_combined_nni %>%  
  group_by(datatype) %>%
  get_summary_stats(NNI)

#datatype           variable     n   min    max  median    q1    q3   iqr   mad  mean    sd    se    ci
#<chr>              <fct>    <dbl> <dbl>   <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#1 alldatatype        NNI        213 0.007  4.64  0.177 0.113 0.305 0.193 0.113 0.31  0.448 0.031 0.06 
#2 citizenscience     NNI        187 0.007  2.26  0.138 0.08  0.28  0.201 0.111 0.294 0.412 0.03  0.059
#3 non_citizenscience NNI        178 0      2.27  0.224 0.16  0.324 0.164 0.107 0.349 0.368 0.028 0.054
 
# Filter dataset for just citizen and non-citizen science data
cs_vs_noncs_nni <- subset(filtered_combined_nni, datatype != "alldatatype")

# Test for significant differences in the distribution of NNI
# Between CS and non-CS data (NNI for both datasets are not normally distributed)
# Using non-parametric Wilcoxon rank sum test

stat.test <- cs_vs_noncs_nni %>% 
  wilcox_test(NNI ~ datatype) %>%
  add_significance()

stat.test
# A tibble: 1 × 8
#.y.        group1         group2                n1    n2 statistic            p p.signif
#<chr>      <chr>          <chr>              <int> <int>     <dbl>        <dbl> <chr>   
#  1 NNI   citizenscience non_citizenscience   187   178     11145 0.0000000487 ****    

### There is a significant difference in the NNI distribution between groups ###

cs_vs_noncs_nni %>% rstatix::wilcox_effsize(NNI ~ datatype)
# A tibble: 1 × 7
#.y.      group1          group2              effsize    n1    n2 magnitude
#* <chr>  <chr>           <chr>                <dbl>   <int> <int> <ord>    
#  1 NNI   citizenscience non_citizenscience   0.286    187   178 small    

### With a small effect size ###
################################################################################
### Taxonomic and temporal trends and biases

### Plot taxonomic trends with histogram
# (log10 transform the number of records per spp.)
plothist_records <- ggplot(records_by_species, aes(x=log10(totalrecords))) + 
  geom_histogram(binwidth = 0.3, color = "black", fill ="grey", lwd =0.2) +
  labs(x = "Number of records", y = "Number of species")+
  theme_classic(base_size = 10) +
  scale_x_continuous(breaks =c(0,1,2, 3, 4, 5),
                     labels=c("1","10","100","1000", "10000","100000")) +
  theme(axis.text = element_text(size=10, colour = "black"))

ggsave(plot = plothist_records, "result/figS4_histogram_records_frequency.png", width = 5.5, height = 4, dpi = 1200)


### Plot number of records per year and trend in record accumulative

# Group records by year
dfhx <- allfrogdata %>% 
  group_by(factor(year))%>% 
  summarise(N = n())

colnames(dfhx) <- c("Year", "N")
dfhx$Year <- as.double(as.character(dfhx$Year))
dfhx <- filter(dfhx, between(Year, 1950,2023))
# Calculate cumulative sum of records
dfhx$cum.sum <- cumsum(dfhx$N)

cs_by_year <- cs_records_df %>% 
  group_by(factor(year)) %>% 
  summarise(N = n())
colnames(cs_by_year) <- c("Year", "N")

### Plot number of records per year and trend in record accumulation

plot_records_accum <-ggplot(data=dfhx)+
  geom_col(aes(x = as.numeric(as.character(Year)), y= N/1e3),width = .8, fill = "dodgerblue4") +
  geom_col(data = cs_by_year, aes(x = as.numeric(as.character(Year)), y= N/1e3),width = .8, fill = "darkred") +
  geom_line(aes(x = as.numeric(as.character(Year)), y= (cum.sum/1e3)/5), colour = "grey30", linewidth=1) +
  scale_y_continuous(name = 'Number of records (x thousand)',
                     expand = c(0,0),
                     sec.axis = sec_axis(~.*5,name="Number of cumulated records (x thousand)")) +
  scale_x_continuous(name = "Year", expand = c(0,0), 
                     breaks = seq(1950,2023,10),limits= c(1950,2023)) +
  theme_classic(base_size = 14)+
  theme(axis.line = element_line(linewidth = 0.4),
        axis.line.y.left = element_line(colour = "dodgerblue4"),
        axis.text.y.left = element_text(colour = "dodgerblue4"),
        axis.title.y.left = element_text(colour = "dodgerblue4"),
        axis.ticks.y.left = element_line(colour = "dodgerblue4"),
        axis.line.y.right = element_line(colour = "grey30"),
        axis.text.y.right = element_text(colour = "grey30"),
        axis.title.y.right = element_text(colour = "grey30"),
        axis.ticks.y.right = element_line(colour = "grey30"))

ggsave(plot=plot_records_accum, "result/fig5_temporaltrend_record_accumulation.png", dpi = 1200, width = 8, height = 5)

################################################################################
### Assessing biases in occurrence records with occAssess package 
### Assess changes in data coverage over time

# Define periods of time (years)
periods <- list(1950:1959, 1960:1969, 1970:1979, 1980:1989, 1990:1999, 2000:2009, 2010:2019, 2020:2023)

brks <- c(2,4,6,8)
lb <- c("1960", "1980","2000", "2020")

### Identifies temporal variation in sampling intensity
# Subset data and only include citizen science collected from 2000 onwards
frogdata_subset <-  frog_occ_cleaned_final[!(frog_occ_cleaned_final$datatype == "citizenscience" & frog_occ_cleaned_final$year < 2000), ]

record_bias <- assessRecordNumber(dat = frogdata_subset,
                                  periods = periods,
                                  species = "species",
                                  x = "decimalLongitude",
                                  y = "decimalLatitude",
                                  year = "year", 
                                  spatialUncertainty = "coordinateUncertaintyInMeters",
                                  identifier = "datatype",
                                  normalize = FALSE)

plot_record_bias <- record_bias$plot +  
  scale_x_continuous(name ="Period",breaks = brks,labels =lb) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))


### Assess changes in the number of species collected over time

species_bias <- assessSpeciesNumber(dat = frogdata_subset,
                                    periods = periods,
                                    species = "species",
                                    x = "decimalLongitude",
                                    y = "decimalLatitude",
                                    year = "year", 
                                    spatialUncertainty = "coordinateUncertaintyInMeters",
                                    identifier = "datatype",
                                    normalize = FALSE)

plot_species_bias <- species_bias$plot +  
  scale_x_continuous(name ="Period",breaks = brks,labels =lb) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))

### Assess taxonomic biases over time
# Indicates whether rare species are over-represented in the data
# and whether the degree to which they are overrepresented changes over time

rarity_bias <- assessRarityBias(dat = frogdata_subset,
                                periods = periods,
                                res = 0.1, # 0.1 degree or about 10km
                                prevPerPeriod = FALSE,
                                species = "species",
                                x = "decimalLongitude",
                                y = "decimalLatitude",
                                year = "year", 
                                spatialUncertainty = "coordinateUncertaintyInMeters",
                                identifier = "datatype")

plot_rarity_bias <-rarity_bias$plot  +
  scale_x_continuous(name ="Period",breaks = brks,labels =lb) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))

### Assess spatial bias over time with NNI 

# Define the study region, a raster object, and resolution of the raster
aus_extent <- extent(aus)

aus_raster <- raster(crs = crs(aus), 
                     vals = 0, resolution = c(0.1, 0.1), ext = aus_extent) 


spatial_bias <- assessSpatialBias(dat = frogdata_subset,
                                  periods = periods,
                                  mask = aus_raster,
                                  nSamps = 10, # number of iterations of random samples
                                  degrade = TRUE, # remove duplicated coordinates
                                  species = "species",
                                  x = "decimalLongitude",
                                  y = "decimalLatitude",
                                  year = "year", 
                                  spatialUncertainty = "coordinateUncertaintyInMeters",
                                  identifier = "datatype")

plot_spatial_bias <- spatial_bias$plot + 
  scale_x_continuous(name ="Period",breaks = brks,labels =lb) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

# Combine plots into one figure
plot_biases <- plot_grid(plot_record_bias,
                         plot_species_bias,
                         plot_rarity_bias,
                         plot_spatial_bias, nrow = 2, ncol = 2,
                         labels=c("(a)", "(b)", "(c)", "(d)"), 
                         hjust = 0,
                         label_size = 14)

ggsave(plot = plot_biases, "result/figS5_occAssess_datacover_over_time.png", width = 15, height = 8, dpi = 1200)

### Assess the extent to which the data are spatiao-temporally biased
## which is the extent to which the same area has been sampled over time

spatialcov_bias <- assessSpatialCov(dat = frogdata_subset,
                                    periods = periods,
                                    res = 0.5 ,
                                    logCount = TRUE,
                                    countries = "Australia",
                                    species = "species",
                                    x = "decimalLongitude",
                                    y = "decimalLatitude",
                                    year = "year", 
                                    spatialUncertainty = "coordinateUncertaintyInMeters",
                                    identifier = "datatype",
                                    output = "nPeriods" )


ncs_spatialcov <- na.omit(spatialcov_bias$non_citizenscience) + 
  scale_fill_viridis_d(name = "Number of periods sampled") +
  coord_sf(xlim = c(110, 155), ylim = c(-45, -10)) +
  labs(x = "Longitude", y = "Latitude") 

cs_spatialcov <- na.omit(spatialcov_bias$citizenscience) + 
  scale_fill_viridis_d(name = "Number of periods sampled") +
  coord_sf(xlim = c(110, 155), ylim = c(-45, -10)) +
  labs(x = "Longitude", y = "Latitude") 


### Assess spatial-temporal bias for the entire dataset

alldata_spatialcov_bias <- assessSpatialCov(dat = frogdata_subset,
                                            periods = periods,
                                            res = 0.5 , 
                                            logCount = TRUE,
                                            countries = "Australia",
                                            species = "species",
                                            x = "decimalLongitude",
                                            y = "decimalLatitude",
                                            year = "year", 
                                            spatialUncertainty = "coordinateUncertaintyInMeters",
                                            identifier = "country" ,
                                            output = "nPeriods" )


alldata_spatialcov <-  na.omit(alldata_spatialcov_bias$Australia) + 
  scale_fill_viridis_d(name = "Number of periods sampled") +
  coord_sf(xlim = c(110, 155), ylim = c(-45, -10)) +
  labs(x = "Longitude", y = "Latitude") 



### Combine spatial-temporal bias plots into one figure

plot_spatialcov <- plot_grid(alldata_spatialcov,
                             ncs_spatialcov,
                             cs_spatialcov,
                             plot_spatial_bias, nrow = 3, ncol = 1,
                             labels=c("(a)", "(b)", "(c)"), 
                             hjust = 0,
                             label_size = 14)

ggsave(plot = plot_spatialcov, "result/figS6_sampling_frequency_bias_0.5res.png", width = 8, height = 10, dpi = 1200)

############### End of analysis for data trends and biases #####################
################################################################################
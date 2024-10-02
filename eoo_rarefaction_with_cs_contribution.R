################################################################################
### EOO Estimates and Completeness and Citizen science contribution to EOO estimate
###
### Script for estimating the extent of occurrence (EOO) and the completeness of EOO
### EOO completeness is calculated from the slope of the last 20% of EOO rarefaction curve
### Which is 1 - the slope of the last 20% of the EOO rarefaction curve
### This script also estimates EOO without citizen science data and 
### Quantify the contribution of CS to EOO estimate as:
### (EOO_alldata - EOO_noncs_data)/EOO_alldata
### Part of the methods for the manuscript:
### How well do we understand geographic range size?: 
### A case study of Australia’s frogs and citizen science project
###
###
################################################################################
# Load packages and data objects
library(adehabitatHR)
library(sp)
library(sf)
library(tidyverse)
library(ozmaps)
library(cowplot)

load("data/frog_occ_cleaned_final_with_datatype_20240702.Rda")


allrecords <-frog_occ_cleaned_final %>% 
  dplyr::select(c(species,decimalLongitude,decimalLatitude,datatype))

## Check to ensure there aren't any missing values
allrecords <- allrecords[!is.na(records$decimalLatitude) & !is.na(records$decimalLongitude),]

records_by_species <- allrecords %>% group_by(species) %>% 
  summarise(totalrecords =n(), 
            cs_records = sum(datatype == "citizenscience"),
            noncs_records = sum(datatype == "non_citizenscience"))

check <- allrecords %>% group_by(species) %>% 
  summarise(relocations=length(species)) %>% 
  dplyr::filter(relocations <6)

# Retain species with at least 6 occurrence records
records_filtered <- allrecords %>% anti_join(check)

frogs_sp <- records_filtered

frogs_mcp <- records_filtered[,-4]

## Convert data frame to sp object
coordinates(frogs_sp) <-c("decimalLongitude","decimalLatitude")
coordinates(frogs_mcp) <-c("decimalLongitude","decimalLatitude")

## Set the coordinate reference system to EPSG:4326 and convert to EPGS:3577 for Australian Albers coordinate system
proj4string(frogs_sp) <- '+init=epsg:4326'
proj4string(frogs_mcp) <- '+init=epsg:4326'

frogs_sp <- spTransform(frogs_sp, CRS('+init=epsg:3577'))
frogs_mcp <- spTransform(frogs_mcp,  CRS('+init=epsg:3577'))

### Calculate EOO for all species with >=6 records using Minimum convex polygon (MCP)
frogs_eoo <- mcp(frogs_mcp, percent = 100, unout = "km2")
frogs_eoo_df <- data.frame(frogs_eoo)
colnames(frogs_eoo_df) <- c("species", "eoo_area")

################################################################################
### Estimate EOO and EOO rarefaction curve

### Create list to store results
mcp_allrarefaction <- list() # results for all datatypes
mcp_cs_rarefaction <- list() # results for citizen science data
mcp_noncs_rarefaction <- list() # results for non-citizen science data

raretab <- matrix(nrow = 100, ncol = 20)
cs_raretab <- matrix(nrow = 100, ncol = 20)
noncs_raretab <- matrix(nrow = 100, ncol = 20)

### Sequence from 5% to 100% of the data
rareseq <- seq(from = 5, to = 100, by=5)/100

# Number of iterations = 100
runs <- 100

species <- sort(unique(frogs_sp$species))
#species <- "Litoria rheocola"
#species <- "Litoria watjulumensis"
species <- c("Crinia signifera", "Adelotus brevis")
#species <- 	"Mixophyes balbus"

### WARNING - This loops takes a very long time to run
# This loops calculate EOO at some proportion of the data rareseq j)
# For 100 times 
# For each species

for (k in 1:length(species)) {
  for (i in 1:runs) {
    for (j in 1:length(rareseq)){
      species_tab <-frogs_sp[frogs_sp$species == species[k],]
      curve_points <- round(rareseq*nrow(species_tab))
      rando <- sample( 1: nrow(species_tab), size = curve_points[j])
      species_tab <- species_tab[rando,]
      if(nrow(species_tab) > 5){
        species_mcp <- mcp(species_tab,percent = 100, unout = "km2" )
        raretab[i,j] <-species_mcp$area
        mcp_allrarefaction[[species[k]]] <- raretab}
      else {raretab [i,j] <- NA
      mcp_allrarefaction[[species[k]]] <- raretab}
      
      #for citizen science data
      cs_species_tab <- subset(species_tab, datatype == "citizenscience")
      if(nrow(cs_species_tab) > 5) {
        cs_species_mcp <- mcp(cs_species_tab,percent = 100, unout = "km2" )
        cs_raretab[i,j] <-cs_species_mcp$area
        mcp_cs_rarefaction[[species[k]]] <- cs_raretab
      } else { cs_raretab[i, j] <- NA
      mcp_cs_rarefaction[[species[k]]] <- cs_raretab}
      
      #for non_citizen science data
      noncs_species_tab <- subset(species_tab, datatype == "non_citizenscience")
      if(nrow(noncs_species_tab) > 5) {
        noncs_species_mcp <- mcp(noncs_species_tab,percent = 100, unout = "km2" )
        noncs_raretab[i,j] <-noncs_species_mcp$area
        mcp_noncs_rarefaction[[species[k]]] <- noncs_raretab
      } else {noncs_raretab[i, j] <- NA
      mcp_noncs_rarefaction[[species[k]]] <- noncs_raretab}
      
      
    }
  }
  setTxtProgressBar(pb, value = k)
}

save(mcp_allrarefaction, file = "result/eoo_rarefaction_results_final.Rda")
save(mcp_cs_rarefaction,  file = "result/eoo_rarefaction_results_final_cs_records.Rda")
save(mcp_noncs_rarefaction, file = "result/eoo_rarefaction_results_final_noncs_records.Rda")

################################################################################
### Plotting EOO rarefraction curve

mean_mcp_tab <- list()
cs_mean_mcp_tab <- list()
noncs_mean_mcp_tab <- list()
cs.contribution_tab <- list()

aus_states <- ozmap_states

species <- sort(unique(frogs_sp$species))

### WARNING: This loops takes about half an hour to run
### This loop estimates the proportion of EOO calculated with certain % of records
### And plots EOO rarefaction curve as well as species distribution

start.time <- Sys.time()

pb <- txtProgressBar(max = length(species), style = 3)
for (k in 1:length(species)) {
  mcp_tab <- as.data.frame(mcp_allrarefaction[(species[k])])
  mcp_tab[is.na(mcp_tab)] <- 0
  mcp_mean <- colMeans(mcp_tab)
  mcp_sd <- apply(mcp_tab, 2,sd)
  percentages <- seq(from = 5, to = 100, by=5)
  average_mcp <- data.frame(mcp_mean, mcp_sd, percentages)
  mean_mcp_tab[[species[k]]] <- average_mcp
  
  # For citizen science data
  cs_mcp_tab <- as.data.frame(mcp_cs_rarefaction[(species[k])])
  cs_mcp_tab[is.na(cs_mcp_tab)] <-0
  cs_mcp_mean <- colMeans(cs_mcp_tab)
  cs_mcp_sd <- apply(cs_mcp_tab, 2,sd)
  percentages <- seq(from = 5, to = 100, by=5)
  cs_average_mcp <- data.frame(cs_mcp_mean, cs_mcp_sd, percentages)
  cs_mean_mcp_tab[[species[k]]] <- cs_average_mcp
  
  # For non citizen science data
  noncs_mcp_tab <- as.data.frame(mcp_noncs_rarefaction[(species[k])])
  noncs_mcp_tab[is.na(noncs_mcp_tab)] <- 0
  noncs_mcp_mean <- colMeans(noncs_mcp_tab)
  noncs_mcp_sd <- apply(noncs_mcp_tab, 2,sd)
  percentages <- seq(from = 5, to = 100, by=5)
  noncs_average_mcp <- data.frame(noncs_mcp_mean, noncs_mcp_sd, percentages)
  noncs_mean_mcp_tab[[species[k]]] <- noncs_average_mcp
  
  # EOO Completeness
  # Percentage of EOO calculated with 80% and 100% of records (%EOO80 - %EOO100)
  # Percentage of EOO calculated with 100% of records will be 100, thus it's %EOO80 - 100
  eoo_prop_at_80percent <- (average_mcp[16,"mcp_mean"]/average_mcp[20,"mcp_mean"])*100 - 100 
  # Completeness is calculated as 1 - the slope of the last 20% of EOO rarefaction:
  # 1 - (%EOO80 - %EOO100)/(N80-N100) where  N80 and N100 is the 80 and 100 percent of records, respectively.
  mean_mcp_completeness  <- 1 - (eoo_prop_at_80percent/-20) 
  mcp_completeness_sigfig <- signif(mean_mcp_completeness, digits = 3)
  species_records <- records_by_species[records_by_species$species == species[k],]
  #allrds_position <-  average_mcp$mcp_mean[1] - average_mcp$mcp_sd[1]
  #cs_rds_position <- mcp_tab[18,19]
  #noncs_rds_position <- mcp_tab[16,19]
  #mcp_tab$slope <- ((mcp_tab[,20] - mcp_tab[,16])/mcp_tab[,20]*100)/(100-80)
  #mean_slope <- mean(mcp_tab$slope)
  #mean_slope_3sigfig <- signif(mean_slope, digits = 3)
  
  # Citizen science contribution
  cs_mcp_tab_na_rm <- cs_mcp_tab
  cs_mcp_tab_na_rm[is.na(cs_mcp_tab_na_rm)] <- 0
  noncs_mcp_tab_na_rm <- noncs_mcp_tab
  noncs_mcp_tab_na_rm[is.na(noncs_mcp_tab_na_rm)] <-0
  cs.contribution <- (mcp_tab[,20] - noncs_mcp_tab_na_rm[,20])/mcp_tab[,20]*100
  mean_cs.contribution <- mean(cs.contribution)
  mean_cs.contribution_3sigfig <- signif((mean_cs.contribution), digits =3)
  eoo_rarefraction_summary_tab <- data.frame(species_records,mean_mcp_completeness, mean_slope, mean_cs.contribution)
  cs.contribution_tab[[species[k]]] <- eoo_rarefraction_summary_tab
  
  
  # Create a map of species distribution
  spp <- gsub(" ","_", species[k])
  
  # Load species map for species i. Data source = Australian Frog Atlas (AFA)
  species_map <- st_read(paste0("data/australian_frog_atlas_frogid/species_map_v2_shp/", spp,".shp"))
  # Subset the data for the current species
  occurrences <- frog_occ_cleaned_final[frog_occ_cleaned_final$species == species[k], ]
  
  # Estimate the arrow length for citizen science contribution
  y_noncs <- noncs_mcp_tab_na_rm[1,20]
  y_allrecs <- mcp_tab[1,20]
  y_dist <- ((mcp_tab[1,20] - noncs_mcp_tab_na_rm[1,20])/2) + noncs_mcp_tab_na_rm[1,20]
  
  # Plot distribution of occurrence records
  occurrence_map <- ggplot() + 
    geom_sf(data = aus_states, fill = "#FBFBEF", linewidth=0.3, colour = "black")+ 
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = (occurrences), 
               aes(x =decimalLongitude, y=decimalLatitude, 
                   colour = factor(datatype)),
               size = 0.8) +
    scale_colour_manual(values = c("citizenscience" = "#D81B60", "non_citizenscience" = "#1E88E5"),
                        labels=c("Citizen science", "Non-citizen science"))+
    scale_shape_manual(values = c(4, 1)) +
    geom_sf(data = species_map, fill = NA, colour = "#FFC107",linewidth=1) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    guides(shape = guide_legend(override.aes = list(size = 3))) +
    annotate("text", label = species[k], x = 125, y = -10, size =6,  fontface = 'italic') +
    annotate("text", label = bquote("cs records = " ~.(species_records$cs_records)), x = 125, y = -45, size = 6) +
    annotate("text", label = bquote("non-cs records = " ~.(species_records$noncs_records)), x = 155, y = -45, size = 6)+
    theme(legend.position = c(0.85,0.8),
          legend.title = element_blank(),
          legend.text = element_text(colour="black", size =15),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.key.size = unit(1,"cm"),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          panel.background = element_blank())

  # Plot EOO rarefaction curve
  eoo_plot <- ggplot(data = average_mcp, aes(x = percentages, y = mcp_mean)) +
    geom_smooth() +
    geom_smooth(data = noncs_average_mcp, aes(x = percentages, y = noncs_mcp_mean ), size = 0.8, colour = "black") +
    scale_x_continuous(breaks = seq(from = 0, to = 100, by = 10)) +
    labs(x = "Percentage of total records", y = bquote ("Mean EOO in " ~ km^2)) + 
    geom_segment(aes(x=100, y = y_noncs , xend = 100, yend= y_allrecs ),
                 arrow = arrow(length=unit(0.60,"cm")), size = 2, colour = "red") +
    annotate("text", label=bquote("Citizen science contribution (%) ="~.(paste0(mean_cs.contribution_3sigfig))), x = 72, y = y_dist, size = 5) +
    annotate("text", label = paste0(species[k]), x = 30, y = mcp_tab[20,20], size = 5, fontface = 'italic') +
    annotate("text", label = bquote("Proportion of EOO completeness ="~.(paste0(mcp_completeness_sigfig))), x = 60, y = y_allrecs*1.1) +
    theme_classic() +
    theme(axis.title = element_text(size = 16, colour = "black"),
          axis.text = element_text(size = 14, colour = "black"),
          panel.background = element_blank())
  
  # Combine plots into 1 figure
  combined_plots <- plot_grid(occurrence_map, eoo_plot, nrow = 2, ncol = 1,labels=c("(a)", "(b)", label_size = 14))
  ggsave(plot = combined_plots, paste0("result/occ_distribution_and_eoo_rarefaction/",species[k],"_map_and_eoo_rarefraction.png"), width = 8, height = 12, dpi = 1200)
  
  setTxtProgressBar(pb, value = k) 
  
}

end.time <- Sys.time()
time_diff2 <-  end.time - start.time

save(cs.contribution_tab, file = "result/cs_contribution_to_eoo_rarefaction_list.Rda")
save(mean_mcp_tab, file = "result/mean_eoo_rarefaction_list.Rda")
save(cs_mean_mcp_tab, file = "result/citizenscience_mean_eoo_rarefaction_list.Rda")
save(noncs_mean_mcp_tab, file = "result/noncitizenscience_mean_eoo_rarefaction_list.Rda")

### Covert list to dataframe
df_cs_contribution <- do.call(rbind.data.frame, cs.contribution_tab)

# Add estimated EOO column "eoo_area", to df_cs_contribution dataframe
df_cs_contribution <- left_join(frogs_eoo_df, df_cs_contribution, by = "species")
save(df_cs_contribution, file = "result/eoo_completeness_and_cs_contribution_table_final.Rda")
write.csv(df_cs_contribution, file = "result/eoo_completeness_and_cs_contribution_final.csv", row.names = FALSE)

mean_cs_contribtuion <- df_cs_contribution %>% 
  subset(totalrecords >= 100) %>%
  get_summary_stats(mean_cs.contribution)

# A tibble: 1 × 13
#variable                    n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
#<fct>                     <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#  1 mean_cs.contribution   138     0   100   2.99 0.061  15.9  15.8  4.44  18.0  30.9  2.63   5.2

hist_cs_eooest <- ggplot(subset(df_cs_contribution, totalrecords >=100), aes(x=mean_cs.contribution/100)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill ="grey", lwd =0.2) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(0,20),  breaks = seq(0, 20, by = 4)) +
  geom_vline(xintercept = 0.18, linetype = "longdash") +
  annotate("text", x = 0.2, y = 8, label = "0.18", size = 5) +
  labs(x = "Proportion of citizen science contribution to EOO estimate", y = "Number of species")+
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"))


ggsave(plot = hist_cs_eooest, "result/figS12_cs_contribution_to_eoo_estimate_histogram.png", width = 12, height = 8, dpi = 1200)


######## End of script for estimating EOO estimates and completeness ###########
########### And for calcualting CS contribution to EOO estimates ###############
################################################################################
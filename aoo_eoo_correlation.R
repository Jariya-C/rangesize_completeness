################################################################################
### AOO and EOO Correlation
###
### Script for analysing the relationship between
### Area of Occupancy (AOO) and Extent of Occurrence (EOO)
### Part of the methods for the manuscript:
### How well do we understand species’ geographic range size?: 
### A case study of Australia’s frogs
###
### Jariya Chanachai (jariya.chanachai@hdr.mq.edu.au)
################################################################################

### Load required packages and data objects
library(tidyverse)
library(rstatix)
library(viridis)
library(cowplot)

load("result/aoo_estimates_and_completeness_with_cs_contribution_2km_allspp.Rda")
load("result/eoo_completeness_and_cs_contribution_table_final.Rda")
load("data/frog_occ_cleaned_final_with_datatype_20240702.Rda")


records_by_species <- frog_occ_cleaned_final %>% group_by(species) %>% 
  summarise(totalrecords =n(), 
            cs_records = sum(datatype == "citizenscience"),
            noncs_records = sum(datatype == "non_citizenscience"))

records_by_species$cs_proportion <- records_by_species$cs_records/records_by_species$totalrecords

chao2km_cs_contr_with_datatype <- left_join(records_by_species, chao2km_cs_contribution_allspp, by = "species")

df_cs_contribution <- df_cs_contribution[,c(1,2,6,7,8)]

aoo_eoo_joined <- left_join(chao2km_cs_contr_with_datatype,df_cs_contribution, by ="species")


aoo_eoo_joined_df <- aoo_eoo_joined %>% select("species", "totalrecords","cs_records",
                                               "noncs_records","cs_proportion", "aoo_area",
                                               "eoo_area", "completeness", "mean_mcp_completeness",
                                               "cs_to_AOOestimate", "cs_to_AOOcompleteness", "Cellsize",
                                               "mean_cs.contribution", "redlistCategory", "threatstatus")
  
  
save(aoo_eoo_joined_df, file = "result/aoo_eoo_estimates_and_completeness_df_with_cs_contribution_allspp.Rda")

aoo_eoo_subset <- aoo_eoo_joined_df %>% 
  subset(totalrecords >=100)

aoo_eoo_subset %>%  
  get_summary_stats(eoo_area) #eoo area summary

#variable        n   min      max  median     q1      q3     iqr     mad    mean       sd      se      ci
#    <fct>    <dbl> <dbl>    <dbl>   <dbl>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>   <dbl>   <dbl>
#  1 eoo_area   138  136. 7931538. 267466. 45707. 984850. 939143. 370035. 826733. 1287485. 109598. 216723.


aoo_eoo_subset %>%  
  get_summary_stats(mean_mcp_completeness)

#   variable                  n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
#<fct>                 <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#  1 mean_mcp_completeness   138 0.104 0.996  0.858 0.753 0.931 0.178 0.125 0.812 0.166 0.014 0.028


aoo_eoo_subset %>% 
  group_by(threatstatus) %>%
  get_summary_stats(mean_mcp_completeness)
# A tibble: 2 × 14
#threatstatus.x variable                    n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
#<chr>            <fct>                 <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#1 Not threatened mean_mcp_completeness   118 0.205 0.996  0.869 0.789 0.938 0.149 0.111 0.838 0.141 0.013 0.026
#2 Threatened     mean_mcp_completeness    20 0.104 0.931  0.684 0.543 0.859 0.317 0.257 0.66  0.218 0.049 0.102

### Test for significant difference in distribtuion of eoo completeness 
#   between threatened and non-threatened groups
stat.test <- aoo_eoo_subset %>% 
  wilcox_test(mean_mcp_completeness ~ threatstatus) %>%
  add_significance()

stat.test

# A tibble: 1 × 8
#   .y.                   group1         group2        n1    n2 statistic        p p.signif
#   <chr>                 <chr>          <chr>      <int> <int>     <dbl>    <dbl> <chr>   
#  1 mean_mcp_completeness Not threatened Threatened   118    20      1805 0.000159 ***     

aoo_eoo_subset %>% rstatix::wilcox_effsize(mean_mcp_completeness ~ threatstatus)
# A tibble: 1 × 7
#    .y.                   group1         group2     effsize    n1    n2 magnitude
#   * <chr>                 <chr>          <chr>        <dbl> <int> <int> <ord>    
#  1 mean_mcp_completeness Not threatened Threatened   0.322   118    20 moderate 


## Test for monotonic relationship between AOO and EOO completeness
aoo_eoo_cor <-cor.test(aoo_eoo_subset$completeness, aoo_eoo_subset$mean_mcp_completeness, method = c("spearman"))
aoo_eoo_cor
##Spearman's rank correlation rho

#data:  aoo_eoo_joined$completeness and aoo_eoo_joined$mean_mcp_completeness
#S = 458602, p-value = 0.5832
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.04706283 


cor_plot <- ggplot(data = aoo_eoo_subset, aes(x = mean_mcp_completeness, y =completeness))+ 
  geom_point(size = 1.5)+
  geom_smooth()+
  labs(x = "EOO completeness", y = "AOO completeness")+
  theme_classic()+
  theme(axis.title = element_text(size = 19),
      axis.text = element_text(size = 17),
      panel.background = element_blank())

ggsave(plot = cor_plot, "result/fig6_aoo_eoo_correlation_plot_v2.png", dpi = 600, width = 12, height = 8)

################################################################################
##### Plotting

### Histogram of the distribution of EOO completeness

plothist <- ggplot(aoo_eoo_subset, aes(x=mean_mcp_completeness)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill ="grey", lwd =0.2) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(0,35),  breaks = seq(0, 30, by = 5)) +
  ### plot dash line at the mean EOO completeness value (mean = 0.812)
  geom_vline(xintercept = 0.812, linetype = "longdash") + 
  annotate("text", x = 0.85, y = 30, label = "0.81", size = 5) +
  labs(x = "EOO completeness", y = "Number of species")+
  theme_classic(base_size = 16) +
  theme(legend.text = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 14, colour = "black")) 

### Histogram of EOO completeness for threatened and non-threatened spp

plothist2 <- ggplot(aoo_eoo_subset, 
                    aes(x = mean_mcp_completeness, fill = threatstatus)) +
  geom_histogram(binwidth = 0.05, 
                 lwd = 0.2, 
                 colour = "black",
                 alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5)) +
  labs(x = "EOO completeness", y = "Number of species") + 
  scale_fill_manual(values = c("Threatened" = "#F0E442", "Not threatened" = "#009E73"), 
                    name = "") +
  theme_classic(base_size = 16) +
  theme(legend.position = "top",
        legend.text = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"))

hist_combined <- plot_grid(plothist, plothist2, nrow = 2, ncol = 1,labels=c("(a)", "(b)"), label_size = 16)

ggsave(plot = hist_combined, "result/figS8_eoo_completeness_histograms.png", width = 12, height = 15, dpi = 1200)

### Boxplot of AOO and EOO completeness
### Compared between threatened and non-theatened groups

aoo_threat_boxplot <-ggplot(aoo_eoo_subset, 
                            aes(x = reorder(threatstatus,-completeness), 
                                y = completeness, fill = threatstatus)) +
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  stat_summary(fun.y=mean, geom="point", shape = 17, size = 3, colour = "black")+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ### p-value obtained from wilcoxon test (analysed in aoo_estimate_completeness.R)
  annotate("text", x = 1.5, y = 0.4, label ="p-value = 0.0197", size = 5) +
  xlab("") + ylab("AOO completeness") +
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(),
        axis.ticks.y = element_line(),
        axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size =20, colour = "black"))

eoo_threat_boxplot <-ggplot(aoo_eoo_subset, 
                     aes(x = reorder(threatstatus,mean_mcp_completeness), 
                   y = mean_mcp_completeness, fill = threatstatus)) +
  geom_boxplot()+
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  stat_summary(fun.y=mean, geom="point", shape = 17, size = 3, colour = "black")+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  annotate("text", x = 1.5, y = 0.4, label ="p-value = 0.00016", size = 5) +
  xlab("") + ylab("EOO completeness") +
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(),
        axis.ticks.y = element_line(),
        axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size =20, colour = "black"))

aoo_eoo_threatstatus_bp <- plot_grid(aoo_threat_boxplot, eoo_threat_boxplot, labels = c("(a)","(b)"), nrow = 2, ncol = 1,label_size = 14)

ggsave(plot = aoo_eoo_threatstatus_bp, "result/figS9_aoo_eoo_threatstatus_boxplots.png", width = 10, height = 12, dpi = 1200)


################################################################################
############# Proportion of cs contribution to AOO completeness ################

# Proportion of cs contribution to AOO completeness vs proportion of cs records
# For 138 species with >= 100 records

cs_aooc_vs_cs_proportion <- ggplot(data = aoo_eoo_subset, aes(x = cs_proportion, y =cs_to_AOOcompleteness))+ 
  geom_point()+
  geom_smooth()+
  #geom_smooth(method = "loess")+
  #scale_x_log10() +
  labs(x = "Proportion of citizen science records", y = "Proportion of citizen science contribution to AOO completeness")+
  theme_classic()+
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        panel.background = element_blank())

ggsave(plot = cs_aooc_vs_cs_proportion, "result/figS11_cs_contribution_to_aooc_vs_proportion_of_csrecords_loessfitted_v2.png", dpi = 600, width = 12, height = 8)


aoocomp_cs_contr_and_cspro_spearman <-cor.test(aoo_eoo_subset$cs_proportion, 
                                               aoo_eoo_subset$cs_to_AOOcompleteness, method = c("spearman"))
aoocomp_cs_contr_and_cspro_spearman

#Spearman's rank correlation rho

#data:  aoo_eoo_subset$cs_proportion and aoo_eoo_subset$cs_to_AOOcompleteness
#S = 56484, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.8710383 


## Proportion of CS contribution to AOO and EOO estimates
load("C:/Users/45790779/OneDrive - Macquarie University/PhD Project/2.-Experiment/SpeciesRangeSize/Results/records_aoo_eoo_20240823.Rda")

combined_df$cs_proportion <- combined_df$citizen_science_records/combined_df$total_records 

cs_aooest <- ggplot(data = aoo_eoo_subset, aes(x =cs_proportion, y =cs_to_AOOestimate))+ 
  geom_point() + 
  geom_smooth() +
  labs(x = "Proportion of citizen science records", y = "Proportion of citizen science contribution to AOO estimate")+
  theme_classic()+
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        panel.background = element_blank())


ggsave(plot = cs_aooest, "result/figS10a_cscontribution_aooest_loessfitted.png", dpi = 600, width = 12, height = 8)


cs_aooest_cortest <-cor.test(aoo_eoo_subset$cs_proportion, 
                             aoo_eoo_subset$cs_to_AOOestimate, method = c("spearman"))


#Spearman's rank correlation rho

#data:  aoo_eoo_subset$cs_proportion and aoo_eoo_subset$cs_to_AOOestimate
#S = 70150, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.8398357 


cs_eooest <- ggplot(data = aoo_eoo_subset, aes(x = cs_proportion, y =mean_cs.contribution/100))+ 
  geom_point()+
  geom_smooth()+
  labs(x = "Proportion of citizen science records", y = "Proportion of citizen science contribution to EOO estimate")+
  theme_classic()+
  theme(axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        panel.background = element_blank())



ggsave(plot = cs_eooest, "result/figS10b_cscontribution_eooest_loessfitted.png", dpi = 600, width = 12, height = 8)


cs_eooest_cortest <-cor.test(aoo_eoo_subset$cs_proportion, 
                             aoo_eoo_subset$mean_cs.contribution, method = c("spearman"))
cs_eooest_cortest

#Spearman's rank correlation rho

#data:  combined_df$cs_proportion and combined_df$mean_cs.contribution
#S = 180801, p-value = 3.74e-14
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.5872022

cscontribution_eoo_aoo <- plot_grid(cs_aooest, cs_eooest, nrow = 2, ncol = 1,labels=c("(a)", "(b)"), label_size = 16)

ggsave(plot = cscontribution_eoo_aoo, "result/figS10_cscontribution_eoo_aooestimates_loessfitted.png", width = 12, height = 15, dpi = 1200)

################################################################################
########## End of script for AOO and EOO correlation ###########################
################################################################################
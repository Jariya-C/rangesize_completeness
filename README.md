# Estimates of geographic range size completeness

## Overview

This study evaluates the completeness of geographic range size for Australia's frog species. It assesses trends and biases in the occurrence records, estimates how well these records capture the likely geographic range of Australia's frogs and quantifies the additional value of citizen science to this knowledge. Specifically, it estimates the completeness of the area of occupancy (AOO) and the extent of occurrence (EOO) based on the statistical distribution of rarely sampled AOO grid cells (2 x 2 km) and EOO rarefaction curves. It also compares whether the completeness of AOO and EOO differences between threatened and non-threatened species and examines the relationship between AOO and EOO completeness.

A set of R scripts used for data collection, data cleaning and analyses and the required dataset are provided.

## Manuscript title

How well do we understand geographic range size?: A case study of Australia's frogs and citizen science projects (under review)

## Data download and initial set up

1.  Install R (<https://www.r-project.org/>) and R studio for desktop (<https://posit.co/download/rstudio-desktop/>)
2.  Download all data files and R scripts and set project directory (name it `rangesize_completeness`). Store all data files in this project folder
3.  Download additional data that are publicly accessible: (i) [Australian Frog Atlas](https://zenodo.org/doi/10.5281/zenodo.5493513) and (ii) [Australia's bioregions (IBRA)](https://datasets.seed.nsw.gov.au/dataset/interim-biogeographic-regionalisation-for-australia-ibra-version-7-subregions)
4.  Rename the `Australian Frog Atlas` folder as `australian_frog_atlas_frogid`
5.  Rename the `Australia's bioregions (IBRA)` folder as `IBRA7_regions` and rename all files within this folder as `ibra_regions`
6.  Create two sub folders within the project: `data` and `result`
7.  Move all datasets to the `data` sub folder and keep R code in the main project directory, these include all file types, except for the **.R files**.
8.  Create another sub folder within the data sub folder called `iucn_redlist_20240521_final`
9.  Install the following are packages:

-   `galah`
-   `CoordinateCleaner`
-   `countrycode`
-   `tidyverse`
-   `reshape2`
-   `sf`
-   `spdep`
-   `terra`
-   `raster`
-   `sp`
-   `ggspatial`
-   `viridis`
-   `cowplot`
-   `ozmaps`
-   `ggrepel`
-   `foreach`
-   `doParallel`
-   `occAssess`
-   `iNext`
-   `adehabitatHR`
-   `spatialEco`
-   `rstatix`

### *Data formats*

-   **.R files**: Contain R script for data download, cleaning and analysis.
-   **.Rda files**: Contains the raw and clean datasets of species occurrence data.
-   **.csv files**: Contains the frog species list, and the IUCN Redlist assessment for the frog species.
-   **.shp files**: Contains shapefiles for the Australian Frog Altas and the IBRA bioregions.

### *Dataset*

1.  `Australian-frogs-species-list-from-AFD-SpeciesLevel-20230614.csv` is the frog species list from the [Australian Faunal Directory](https://biodiversity.org.au/afd/home)
2.  `ala_frog_occurrences_raw_20231017.Rda` is the raw species occurrence record downloaded from the [Atlas of Living Australia](https://www.ala.org.au/)
3.  [`australian_frog_atlas_frogid`](https://zenodo.org/doi/10.5281/zenodo.5493513) is the [Australian Frog Atlas](https://zenodo.org/doi/10.5281/zenodo.5493513)
4.  `IBRA7_regions` is the [Australia's bioregions (IBRA)](https://datasets.seed.nsw.gov.au/dataset/interim-biogeographic-regionalisation-for-australia-ibra-version-7-subregions)
5.  `assessments_20240521.csv` is the IUCN Redlist assessment for Australia's frog species downloaded from [IUCN](https://www.iucnredlist.org/)
6.  `iucn_redlist_namesstandardised.csv` is the IUCN Redlist assessment with species name standardised
7.  `frog_occ_cleaned_final_with_datatype_20240702.Rda`is the clean data after data processing

## Data privacy statement

This study involves threatened frog species of Australia. We obtained species occurrence data from (i) the publicly accessible biodiversity data portal: [Atlas of Living Australia](https://www.ala.org.au/) and (ii) the Australian Museum under a data licence agreement (FrogID data). [Atlas of Living Australia](https://www.ala.org.au/) handles sensitive records by (i) withholding the geographic locations or (ii) generalising the coordindates. FrogID data obtained under a licence agreement includes sensitive species records and cannot be shared publicly and is not provided in this dataset. However, users can download publicly available FrogID data (sensitive species records have either been generalised or excluded) from [FrogID](https://www.frogid.net.au/explore).

## Code running

This workflow provides a clear guideline on data collection, processing and analysis and the following steps can be followed to run the code.

**1.** **Download occurrence data (optional):**

-   Open the `download_frogrecords.R` script to download frog species occurrence data from the [Atlas of Living Australia](https://www.ala.org.au/).
-   This code will help to generate a dataset called `ala_frog_occurrences_raw_20231017.Rda`.
-   *(Note: This dataset is already provided, so downloading is optional.)*
-   **Imporant:** The ALA dataset includes FrogID data with sensitive records generalised/excluded. You do not need to download additional data from [FrogID](https://www.frogid.net.au/explore). However, this dataset is slightly different from the actual dataset used for the analyse which includes sensitive species data from FrogID.

**2.** **Data cleaning (optional):**

-   Open `data_cleaning_part1.R` for initial data cleaning.
-   Then, open `data_cleaning_part2.R` to complete data cleaning.
-   The clean dataset is called `frog_occ_cleaned_final_with_datatype_20240702.Rda`
-   *(Note: The clean dataset is already provided, so these steps are optional.)*

**3.** **Assess data trends and biases:**

-   Open and run `assess_cleandata.R` to evaluate trends and biases in the clean dataset

**4.** **AOO estimates and completeness:**

-   Open and run `aoo_estimate_completeness.R` to estimate AOO and the completeness within AOO estimates for all records in clean dataset.

**5.** **Citizen science contribution to AOO estimates and completeness:**

-   Open and run `aoo_cs_contribution.R` to quantify citizen science contribution to AOO estimates and completeness.

**6.** **EOO estimate and completeness and citizen science contribution:**

-   Open and run `eoo_rarefaction_with_cs_contribution.R` to estimate EOO,completeness in EOO, and assess citizen science contribution to EOO estimates.

**7.** **Conduct further analysis:**

-   Open and run `aoo_eoo_further_analysis.R` for additional analysis.
-   This includes comparing AOO and EOO completeness between threatened and non-threatened species, assessing the relationship between AOO and EOO completeness and assessing the relationship between the proportion of citizen science records and the proportion of citizen science contribution.

## Reproducibility

Due to the exclusion or generalisation of sensitive species data from FrogID, results presented in the manuscript may not be exactly reproducible with the data provided. However, results produced should not vary greatly from those presented in the paper.

## References to methods and datasets

-   Atlas of Living Australia. (2023). Occurrence records download on 2023-10-09. Retrieved 9 October 2023 from <https://doi.org/10.26197/ala.fb7754b9-f084-4ee3-b4ab-a87dd313d467>; <https://doi.org/10.26197/ala.1e0ded65-7fa5-439f-8298-081537b871a7>; <https://doi.org/10.26197/ala.abd60ea2-f837-4f21-afdc-a5d41340c820>; <https://doi.org/10.26197/ala.d14ae6e3-c8b6-4d7a-bb75-3dbe25e9a753>; <https://doi.org/10.26197/ala.42348355-ed0e-40ac-853f-98fa30a2d5ac>
-   Atlas of Living Australia. (2024). Retrieved 10 October 2024 from <https://www.ala.org.au/>
-   Australian Faunal Directory. (2023). Retrieved 14 June 2023 from <https://biodiversity.org.au/afd/home>
-   Boyd, R. J., Powney, G. D., Carvell, C., & Pescott, O. L. (2021). occAssess: An R package for assessing potential biases in species occurrence data. *Ecology and Evolution*, *11*(22), 16177-16187. <https://doi.org/10.1002/ece3.8299>
-   Chao, A. (1984). Nonparametric estimation of the number of classes in a population. *Scandinavian Journal of Statistics*, *11*(4), 265-270. <http://www.jstor.org/stable/4615964>
-   Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. *Biometrics*, *43*(4), 783-791. <https://doi.org/10.2307/2531532>
-   Clark, P. J., & Evans, F. C. (1954). Distance to nearest neighbor as a measure of spatial relationships in populations. *Ecology*, *35*(4), 445-453. <https://doi.org/10.2307/1931034>
-   Commonwealth of Australia and Department of Climate Change, Energy, the Environment and Water. (2024). Interim Biogeographic Regionalisation for Australia (IBRA), Version 7 (Subregions). Retrieved from <https://datasets.seed.nsw.gov.au/dataset/b1284c2c-f3bd-4b50-ab07-e1b593b8ee67>
-   Cutajar, T. P., Portway, C. D., Gillard, G. L., & Rowley, J. J. L. (2022). Australian Frog Atlas: species’ distribution maps informed by the FrogID dataset. *Technical Reports of the Australian Museum*, *36*, 1–48. <https://doi.org/10.3853/j.1835-4211.36.2022.1789>
-   Gotelli, N. J., & Colwell, R. K. (2011). Estimating species richness. In *Biological diversity: Frontier in measurement and assessment* (pp. 39-54). Oxford University Press.
-   Gupta, G., Dunn, J., Sanderson, R., Fuller, R., & McGowan, P. J. K. (2020). A simple method for assessing the completeness of a geographic range size estimate. *Global Ecology and Conservation*, *21*, 7. <https://doi.org/10.1016/j.gecco.2019.e00788>
-   IUCN. (2023). The IUCN Red List of Threatened Species. Version 2024-1. <https://www.iucnredlist.org>
-   IUCN Standards and Petitions Committee. (2024). Guidelines for Using the IUCN Red List Categories and Criteria. Version 16. <https://www.iucnredlist.org/documents/RedListGuidelines.pdf>
-   Nipperess, D. A., Faith, D. P., Williams, K. J., King, D., Manion, G., Ware, C., Schmidt, B., Love, J., Drielsma, M., Allen, S., & Gallagher, R. V. (2020). *Using representative sets of known species and habitat condition to inform change in biodiversity status: a case example for vascular plants.* Sydney, NSW, Australia: NSW Department of Planning Industry and Environment. Retrieved from <https://www.environment.nsw.gov.au/research-and-publications/publications-search/using-representative-sets-known-species-habitat-condition-inform-change-biodiversity-status>
-   Ronquillo, C., Stropp, J., & Hortal, J. (2024). OCCUR Shiny application: A user-friendly guide for curating species occurrence records. Methods in Ecology and Evolution, 15(5), 816-823. <https://doi.org/10.1111/2041-210X.14271>
-   Rowley, J. J. L., & Callaghan, C. T. (2020). The FrogID dataset: expert-validated occurrence records of Australia's frogs collected by citizen scientists. *ZooKeys*, *912*(8), 139-151. <https://doi.org/10.3897/zookeys.912.38253>
-   Rowley, J. J. L., Callaghan, C. T., Cutajar, T. P., Portway, C. D., Potter, K., Mahony, S., Trembath, D. F., Flemons, P., & Woods, A. (2019). FrogID: Citizen scientists provide validated biodiversity data on frogs of Australia. *Herpetological Conservation and Biology*, *14*(1), 155 -170.

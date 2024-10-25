# muscidae_phenology

## Content
This repository contains the code and data necessary to replicate data analyses, figures and tables in the manuscript:

***Species' traits modulate rapid changes in flight time in high-Arctic muscid flies under climate change***

## Contact
Anonymized

## Data usage guidelines

### Data

- Phenological observations at Zackenberg were provided by the Greenland Ecosystem Monitoring Programme. Family data for Muscidae is available at: https://data.g-e-m.dk/. Species data is available in this repository. 

- Environmental predictors: Temperature and snowmelt observations for Zackenberg were also provided by the Greenland Ecosystem Monitoring Programme. Data available at: https://data.g-e-m.dk/ The raw data downloaded from the database is included in this repository (downloaded 13. January 2022) along with a formatted version including all estimated temperature and timing of snowmelt estimates.

### Code 
All code provided for data preparation and analysis is licensed under a MIT License. In accordance with the license the code is available to be shared and adapted, but we would appreciate attribution to the authors, e.g. through citation of the above manuscript, and indications where changes were made. Although not mandatory, we additionally suggest that code users contact and collaborate with contributors should the code form a substantial proportion of a particular publication or analysis.

# Data preparation and clean up

Raw Muscidae species abundance data can be found here: 

```
/Data/Zackenberg_Muscidae_SL.xlsx

```

The data preparation, cleaning and assembly scripts can be found here:

```
/Dry_seg_data_cleanup.R
/Dry_gro_data_cleanup.R
/Lim_gro_data_cleanup.R
/Pha_bid_data_cleanup.R
/Spi_alm_data_cleanup.R
/Spi_def_data_cleanup.R
/Spi_dor_data_cleanup.R
/Spi_mal_data_cleanup.R
/Spi_mel_data_cleanup.R
/Spi_mic_data_cleanup.R
/Spi_nov_data_cleanup.R
/Spi_san_data_cleanup.R
/Spi_tun_data_cleanup.R
/Spi_zai_data_cleanup.R

/Data_clean_combine_species_data.R

```

*Important:* Please note that this summarised data is for archival purposes only. If you intend to use the phenological observations in this dataset please refer to the data usage guidance for the raw data sets described above. 

The following path leads directly to the cleaned abundance data (zero capture dates included etc.) used to conduct the analysis:

```
/Data/Species_data_clean/All_species_clean.csv

```

The following path leads directly to the climate variables used to conduct the analysis:


```
/Data/Climate_data/Rolling_mean_temp/Air_temp_30_days_rolling.csv
/Data/Climate_data/Snowmelt/Snowmelt_Climatestation.xlsx

```

*Please note:* All_species_clean.csv includes phenology estimates for analysis of temporal trends. The following datasets includes temperature data related to arthropod phenology estimates.

```
/Data/Species_data_clean/All_species_cleanup_onset.csv
/Data/Species_data_clean/All_species_cleanup_peak.csv
/Data/Species_data_clean/All_species_cleanup_end.csv

```  

**Note:** The scripts used to calculate the rolling mean temperature 30 days before an average phenological event is found here:

```
/Temperature_rolling_mean.R

```

# Generalized Additive Modelling analysis scripts
Scripts to derive phenology estimates through GAM:

```
/peak_week_final.R
/onset_end_week_final.R
/All_models_gam.R

```

# Analysis scripts

```
/peak_week_final.R
/onset_end_week_final.R
/All_models_gam.R

```


# Supplementary material
The following quarto documents with all relevant code and detailed descriptions of each step of GAM, trends in climate variables and linear mixed modelling analysis performed can be used to reproduce the analysis:


```
/Supplementary_info.qmd

```

The following R scripts contains relevant code to reproduce figures 2 and 3

```
/figures_manuscript.R

```

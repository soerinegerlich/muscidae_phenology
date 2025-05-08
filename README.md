# muscidae_phenology

## Content
This repository contains the code and data necessary to replicate data analyses, figures and tables in:

***Species’ traits modulate rapid changes in flight time in high-Arctic muscid flies under climate change***

Authors: Hannah Sørine Gerlich1, Sarah Loboda2, Gavin L. Simpson3, Jade Savage4, Niels M. Schmidt5, Martin Holmstrup1 and Toke T. Høye1

1 Department of Ecoscience and Arctic Research Centre, Aarhus University, C.F. Møllers Allé 4-8, DK-8000 Aarhus C, Denmark 

2 Fisheries and Oceans Canada, Maurice Lamontagne Institute, 850 Rte de la Mer, Mont-Joli, QC G5H 3Z4, Canada

3 Department of Veterinary and Animal Sciences, Aarhus University, Tjele, Denmark

4 Department of Biological Sciences, Bishop’s University, 2600 College Street, Sherbrooke, Qc, Canada, J1M 1Z7

5 Department of Ecoscience and Arctic Research Centre, Aarhus University, Frederiksborgvej 399, 4000 Roskilde, Denmark

Corresponding author: Hannah Sørine Gerlich; E-mail: soger@ecos.au.dk



- Scripts for phenology modeling and trait-based analysis
- Data files derived from raw environmental and phenological observations
- Figure-generation code for the main manuscript and supplementary materials

## Data sources and usage

# Phenology Observations
Phenological data from Zackenberg, Northeast Greenland, were provided by the Greenland Ecosystem Monitoring Programme (GEM).
Source: https://data.g-e-m.dk
Both the raw data (as downloaded) and a cleaned/processed version with estimated phenological events are included in this repository.

# Environmental Predictors
Temperature and snowmelt timing data were also obtained from GEM.
Source: https://data.g-e-m.dk
This repository includes both the original data and a processed version with estimated environmental metrics.

# Species Data
Species trait and phenology data used in the analysis are publicly available on Dryad:
🔗 https://doi.org/10.5061/dryad.3r2280gtm

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
The following path leads directly to the cleaned abundance data (zero capture dates included etc.) used to conduct the analysis:

```
/Data/Species_data_clean/All_species_clean.csv

```

# Data preparation
Cleaned phenological metrics and environmental inputs are located in the following directories:

```
/Data/phenology_data
/Data/Climate_data

```

DNA sequences for phylogenetic analysis:

```
/sequences.fasta

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


# Supplementary material
The following quarto documents with all relevant code and detailed descriptions of each step of GAM, trends in climate variables and linear mixed modelling analysis performed can be used to reproduce the analysis:


```
/Supplementary_info_proceedings.qmd

```

# Trait-based Linear Models

Scripts for linear regression analyses exploring relationships between climate trends and species traits:

```
/temporal_trend_traits.R
/snowmelt_temperature_trends_traits.R

```

# Figures

Scripts to recreate Figures 1–3 in the main manuscript:

```
/figures_manuscript.R
/fig3_new_24022025.R
/figure_peak_with_zeros_removed.R
/figure_onset_with_zeros_removed.R
/figure_end_with_zeros_removed.R

```


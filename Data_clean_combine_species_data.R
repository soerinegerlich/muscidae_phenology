#Combine clean dataset from all species in one data frame

library(tidyverse)


df1 <- read.csv("Data/Species_data_clean\\Dry_gro.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df2 <- read.csv("Data/Species_data_clean\\Dry_seg.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df3 <- read.csv("Data/Species_data_clean\\Pha_bid.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df4 <- read.csv("Data/Species_data_clean\\Spi_alm.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df6 <- read.csv("Data/Species_data_clean\\Spi_dor.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df7 <- read.csv("Data/Species_data_clean\\Spi_mal.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df8 <- read.csv("Data/Species_data_clean\\Spi_nov.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df9 <- read.csv("Data/Species_data_clean\\Spi_san.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df10 <- read.csv("Data/Species_data_clean\\Spi_tun.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df11 <- read.csv("Data/Species_data_clean\\Spi_zai.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df12 <- read.csv("Data/Species_data_clean\\Lim_gro.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df13 <- read.csv("Data/Species_data_clean\\Spi_mel.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df14 <- read.csv("Data/Species_data_clean\\Spi_mic.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)
df15 <- read.csv("Data/Species_data_clean\\Spi_def.csv",sep=",",stringsAsFactors = FALSE, header = TRUE)


df_all <- bind_rows(df1,df2,df3,df4,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15) 

df_all$Abundance <- round(df_all$Abundance, digits = 0) # Need to remove decimals from abundance values

write.csv(df_all, file = "Data/Species_data_clean\\All_species_clean.csv", row.names=FALSE)



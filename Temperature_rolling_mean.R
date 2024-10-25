# Calculating temperature variable using rolling mean function

library(zoo)
library(tidyverse)

#Get the temperature data
dfair1 <- read.csv("Data/Climate_data/Temp_for_sliding_win/Temp_air_sliding_win.csv",sep=",",
                   stringsAsFactors = FALSE, header = TRUE)

#Calculating the rolling mean (air temp)
#dfair1$doymean<-rollmean(dfair1$DOYTemp,k=30,fill=NA,align="right")
#Problem with NA values in 1996 which can be fixed with rollapply
dfair1$doymean<-rollapply(dfair1$DOYTemp,width=30,FUN=mean,na.rm = TRUE,fill=NA,align="right")

ggplot(dfair1, aes(x=DOY, y=doymean)) + 
  geom_line() +
  ylab("Mean rolling 30 day temperature (degrees C)") +
  facet_wrap(~Year)

#Get phenology data for all species

df_peak_all <- read.csv("Data/phenology_estimates/peak_pheno_muscidae.csv",sep=",",
                   stringsAsFactors = FALSE, header = TRUE)

df_onset_all <- read.csv("Data/phenology_estimates/onset_pheno_muscidae.csv",sep=",",
                        stringsAsFactors = FALSE, header = TRUE)

df_end_all <- read.csv("Data/phenology_estimates/end_pheno_muscidae.csv",sep=",",
                        stringsAsFactors = FALSE, header = TRUE)


df_peak_all%>%
  group_by(Species_name)%>% 
  summarise(Peak_meanDOY=mean(.peak,na.rm=T),
            Peak_SD=sd(.peak,na.rm=T))-> df_phen_event_peak_mean

df_onset_all%>%
  group_by(Species_name)%>% 
  summarise(Onset_meanDOY=mean(.onset,na.rm=T),
            Onset_SD=sd(.onset,na.rm=T))-> df_phen_event_onset_mean

df_end_all%>%
  group_by(Species_name)%>% 
  summarise(End_meanDOY=mean(.end,na.rm=T),
            End_SD=sd(.end,na.rm=T))-> df_phen_event_end_mean


df_phen_event_peak_mean%>%
  group_by(Species_name)%>% 
  summarise(Peak_DOY=Peak_meanDOY-Peak_SD)-> df_phen_event_peak_min

df_phen_event_onset_mean%>%
  group_by(Species_name)%>% 
  summarise(Onset_DOY=Onset_meanDOY-Onset_SD)-> df_phen_event_onset_min

df_phen_event_end_mean%>%
  group_by(Species_name)%>% 
  summarise(End_DOY=End_meanDOY-End_SD)-> df_phen_event_end_min

#df_phen_event_min[,-1:-2] <- round(df_phen_event_min[,-1:-2], 0)
df_phen_event_peak_min$Peak_DOY <- round(df_phen_event_peak_min$Peak_DOY, 0)

df_phen_event_onset_min$Onset_DOY <- round(df_phen_event_onset_min$Onset_DOY, 0)

df_phen_event_end_min$End_DOY <- round(df_phen_event_end_min$End_DOY, 0)

#Match with original dataframe

df_peak_all$Peak_DOY <- df_phen_event_peak_min$Peak_DOY[match(paste0(df_peak_all$Species_name),
                                                    paste0(df_phen_event_peak_min$Species_name))]

df_onset_all$Onset_DOY <- df_phen_event_onset_min$Onset_DOY[match(paste0(df_onset_all$Species_name),
                                                    paste0(df_phen_event_onset_min$Species_name))]


df_end_all$End_DOY <- df_phen_event_end_min$End_DOY[match(paste0(df_end_all$Species_name),
                                                    paste0(df_phen_event_end_min$Species_name))]


#Match rollmean temperature 30 days with doy for phen event

df_peak_all$Peak_Temp <- dfair1$doymean[match(paste0(df_peak_all$Year,df_peak_all$Peak_DOY),
                                         paste0(dfair1$Year,dfair1$DOY))]

df_onset_all$Onset_Temp <- dfair1$doymean[match(paste0(df_onset_all$Year,df_onset_all$Onset_DOY),
                                              paste0(dfair1$Year,dfair1$DOY))]

df_end_all$End_Temp <- dfair1$doymean[match(paste0(df_end_all$Year,df_end_all$End_DOY),
                                              paste0(dfair1$Year,dfair1$DOY))]

df_onset_all <- subset(df_onset_all, !df_onset_all$Species_name == "Lim.gro")
df_onset_all <- subset(df_onset_all, !df_onset_all$Species_name == "Spi.mel")
df_onset_all <- subset(df_onset_all, !df_onset_all$Species_name == "Spi.def")

df_peak_all <- subset(df_peak_all, !df_peak_all$Species_name == "Lim.gro")
df_peak_all <- subset(df_peak_all, !df_peak_all$Species_name == "Spi.mel")
df_peak_all <- subset(df_peak_all, !df_peak_all$Species_name == "Spi.def")

df_end_all <- subset(df_end_all, !df_end_all$Species_name == "Lim.gro")
df_end_all <- subset(df_end_all, !df_end_all$Species_name == "Spi.mel")
df_end_all <- subset(df_end_all, !df_end_all$Species_name == "Spi.def")

#Save dataframes for figure
write_xlsx(df_onset_all, "Data/Climate_data/Figure\\df_onset_all.xlsx", col_names = TRUE)
write_xlsx(df_peak_all, "Data/Climate_data/Figure\\df_peak_all.xlsx", col_names = TRUE)
write_xlsx(df_end_all, "Data/Climate_data/Figure\\df_end_all.xlsx", col_names = TRUE)

#Figure of temperature trends

ggplot() +
  geom_point(df_onset_all, mapping = aes(Year, Onset_Temp, color = Species_name), shape = 1) +
  geom_line(df_onset_all, mapping = aes(Year, Onset_Temp, color = Species_name)) +
  geom_point(df_peak_all, mapping = aes(Year, Peak_Temp, color = Species_name), shape = 6) +
  geom_line(df_peak_all, mapping = aes(Year, Peak_Temp, color = Species_name)) +
  geom_point(df_end_all, mapping = aes(Year, End_Temp, color = Species_name), shape = 5) +
  geom_line(df_end_all, mapping = aes(Year, End_Temp, color = Species_name)) +
  ylim(0,9) +
  ylab("Temperature (ºC)") +
  xlab("") +
  scale_color_viridis_d(option = "plasma", name = "Species") +
  facet_wrap( ~ Species_name) +
  theme_bw() +
  theme(strip.background.x = element_rect(fill = "white"),
        strip.text = element_text(face = "italic"),
        plot.margin = margin(6,1,6,1, unit = "cm"),
        axis.title.y = element_text(size = 12, vjust = 4),
        legend.position = "top",
        legend.text = element_text(face = "italic"))

df_peak_all |>
  ggplot(mapping = aes(Year, Peak_Temp, color = Species_name)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d(option = "plasma") +
  facet_wrap(~ Species_name)

df_end_all |>
  ggplot(mapping = aes(Year, End_Temp, color = Species_name)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d(option = "plasma") +
  facet_wrap(~ Species_name)

###################Linear regression of temperature##########################
df_summary_temp <-
  data.frame(
    Species_name = character(),
    pheno_event = character(),
    slope1 = numeric(),
    SE1 = numeric(),
    Tvalue1 = numeric(),
    Pvalue1 = numeric(),
    Rsquare = numeric()
  )

for (i in unique(df_onset_all$Species_name)) {
  print(i)
  df8a <- subset(df_onset_all, Species_name == i)
    
    if (sum(!is.na(df8a$.onset)) < 6) {
      df_other <- data.frame(
        Species_name = df8a$Species_name[1],
        pheno_event = "Onset",
        slope1 = NA,
        SE1 = NA,
        Tvalue1 = NA,
        Pvalue1 = NA,
        Rsquare = NA
      )
    }
    else{
      mod1 <- lm(Onset_Temp ~ Year, data = df8a)
      df_other <- data.frame(
        Species_name = df8a$Species_name[1],
        pheno_event = "Onset",
        slope1 = summary(mod1)$coefficients[2],
        SE1 = summary(mod1)$coefficients[4],
        Tvalue1 = summary(mod1)$coefficients[6],
        Pvalue1 = summary(mod1)$coefficients[8],
        Rsquare = summary(mod1)$r.squared
      )
      df_summary_temp <- bind_rows(df_summary_temp, df_other)
    }
    #plot(mod1)
  
}


df_summary_temp_peak <-
  data.frame(
    Species_name = character(),
    pheno_event = character(),
    slope1 = numeric(),
    SE1 = numeric(),
    Tvalue1 = numeric(),
    Pvalue1 = numeric(),
    Rsquare = numeric()
  )

for (i in unique(df_peak_all$Species_name)) {
  print(i)
  df8a <- subset(df_peak_all, Species_name == i)
  
  if (sum(!is.na(df8a$.peak)) < 6) {
    df_other <- data.frame(
      Species_name = df8a$Species_name[1],
      pheno_event = "Peak",
      slope1 = NA,
      SE1 = NA,
      Tvalue1 = NA,
      Pvalue1 = NA,
      Rsquare = NA
    )
  }
  else{
    mod1 <- lm(Peak_Temp ~ Year, data = df8a)
    df_other <- data.frame(
      Species_name = df8a$Species_name[1],
      pheno_event = "Peak",
      slope1 = summary(mod1)$coefficients[2],
      SE1 = summary(mod1)$coefficients[4],
      Tvalue1 = summary(mod1)$coefficients[6],
      Pvalue1 = summary(mod1)$coefficients[8],
      Rsquare = summary(mod1)$r.squared
    )
    df_summary_temp_peak <- bind_rows(df_summary_temp_peak, df_other)
  }
  #plot(mod1)
  
}


df_summary_temp_end <-
  data.frame(
    Species_name = character(),
    pheno_event = character(),
    slope1 = numeric(),
    SE1 = numeric(),
    Tvalue1 = numeric(),
    Pvalue1 = numeric(),
    Rsquare = numeric()
  )

for (i in unique(df_end_all$Species_name)) {
  print(i)
  df8a <- subset(df_end_all, Species_name == i)
  
  if (sum(!is.na(df8a$.end)) < 6) {
    df_other <- data.frame(
      Species_name = df8a$Species_name[1],
      pheno_event = "End",
      slope1 = NA,
      SE1 = NA,
      Tvalue1 = NA,
      Pvalue1 = NA,
      Rsquare = NA
    )
  }
  else{
    mod1 <- lm(End_Temp ~ Year, data = df8a)
    df_other <- data.frame(
      Species_name = df8a$Species_name[1],
      pheno_event = "End",
      slope1 = summary(mod1)$coefficients[2],
      SE1 = summary(mod1)$coefficients[4],
      Tvalue1 = summary(mod1)$coefficients[6],
      Pvalue1 = summary(mod1)$coefficients[8],
      Rsquare = summary(mod1)$r.squared
    )
    df_summary_temp_end <- bind_rows(df_summary_temp_end, df_other)
  }
  #plot(mod1)
  
}

#Save excel files
#require(writexl)
#
#write_xlsx(df_summary_temp, "Data/model_summaries\\df_summary_temp.xlsx", col_names = TRUE)
#write_xlsx(df_summary_temp_peak, "Data/model_summaries\\df_summary_temp_peak.xlsx", col_names = TRUE)
#write_xlsx(df_summary_temp_end, "Data/model_summaries\\df_summary_temp_end.xlsx", col_names = TRUE)

############INCLUDE TEMP VARIABLE IN PHENOLOGY DATA###########################



df_peak <- read.csv("Data/Species_data_Clean/All_species_clean.csv",sep=",",
               stringsAsFactors = FALSE, header = TRUE)

df_onset <- read.csv("Data/Species_data_Clean/All_species_clean.csv",sep=",",
                    stringsAsFactors = FALSE, header = TRUE)

df_end <- read.csv("Data/Species_data_Clean/All_species_clean.csv",sep=",",
                    stringsAsFactors = FALSE, header = TRUE)

df_peak$Temperature <- df_peak_all$Peak_Temp[match(paste0(df_peak$Year,df_peak$Species_name,df_peak$Plot),
                                              paste0(df_peak_all$Year,df_peak_all$Species_name, df_peak_all$Plot))]

df_onset$Temperature <- df_onset_all$Onset_Temp[match(paste0(df_onset$Year,df_onset$Species_name,df_onset$Plot),
                                                   paste0(df_onset_all$Year,df_onset_all$Species_name, df_onset_all$Plot))]

df_end$Temperature <- df_end_all$End_Temp[match(paste0(df_end$Year,df_end$Species_name,df_end$Plot),
                                                      paste0(df_end_all$Year,df_end_all$Species_name, df_end_all$Plot))]

# Include snowmelt date 

df_snow <- readxl::read_xlsx("Data/Climate_data/Snowmelt\\Snowmelt_Climatestation.xlsx")

#Temporal changes in snowmelt timing figure

df_snow <- subset(df_snow, !df_snow$Year > 2014)

ggplot(df_snow, aes(Year, SnowmeltDOY)) +
  geom_point(size = 3) +
  geom_line() +
  ylab("Snowmelt timing (Day of Year)") +
  xlab("") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14, vjust = 4),
        axis.text = element_text(size = 12),
        plot.margin = margin(6,2,6,2, unit = "cm"))

m.1 <- lm(SnowmeltDOY ~ Year, df_snow)  
summary(m.1)

df_peak <- merge(df_snow, df_peak, by = "Year")

df_onset <- merge(df_snow, df_onset, by = "Year")

df_end <- merge(df_snow, df_end, by = "Year")

#Save new dataframe with climate variables

#write.csv(df_peak, file = "Data/Species_data_clean\\All_species_cleanup_peak.csv", row.names=FALSE)
#write.csv(df_onset, file = "Data/Species_data_clean\\All_species_cleanup_onset.csv", row.names=FALSE)
#write.csv(df_end, file = "Data/Species_data_clean\\All_species_cleanup_end.csv", row.names=FALSE)

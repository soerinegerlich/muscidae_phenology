######### Trends in climate responses according to traits ##############

library(tidyverse)

########ONSET#############

df_onset <- read.csv("Data/phenology_estimates/onset_temp_trend_nas_removed.csv",sep=",",
                     stringsAsFactors = FALSE, header = TRUE)

# Get climate data

df_temp <- read.csv("Data/Species_data_clean/All_species_cleanup_onset.csv",sep=",",
                    stringsAsFactors = FALSE, header = TRUE)

df_temp$Species_name[df_temp$Species_name == "Dry.gro"] <- "D.groenlandica"
df_temp$Species_name[df_temp$Species_name == "Dry.seg"] <- "D.segnis"
df_temp$Species_name[df_temp$Species_name == "Lim.gro"] <- "L.groenlandica"
df_temp$Species_name[df_temp$Species_name == "Pha.bid"] <- "P.bidentata"
df_temp$Species_name[df_temp$Species_name == "Spi.alm"] <- "S.almqvistii"
df_temp$Species_name[df_temp$Species_name == "Spi.def"] <- "S.deflorata"
df_temp$Species_name[df_temp$Species_name == "Spi.dor"] <- "S.dorsata"
df_temp$Species_name[df_temp$Species_name == "Spi.mal"] <- "S.malaisei"
df_temp$Species_name[df_temp$Species_name == "Spi.mel"] <- "S.melanosoma"
df_temp$Species_name[df_temp$Species_name == "Spi.nov"] <- "S.novaesibiriae"
df_temp$Species_name[df_temp$Species_name == "Spi.san"] <- "S.sanctipauli"
df_temp$Species_name[df_temp$Species_name == "Spi.tun"] <- "S.tundrae"
df_temp$Species_name[df_temp$Species_name == "Spi.zai"] <- "S.zaitzevi"



df_onset$Temperature <- df_temp$Temperature[match(paste0(df_onset$year,df_onset$species,df_onset$plot),
                                                  paste0(df_temp$Year,df_temp$Species_name,df_temp$Plot))]

df_onset$Snowmelt <- df_temp$SnowmeltDOY[match(paste0(df_onset$year,df_onset$species,df_onset$plot),
                                               paste0(df_temp$Year,df_temp$Species_name,df_temp$Plot))]

df_summary_onset <-
  data.frame(
    species = character(),
    plot = character(),
    pheno_event = character(),
    slope1 = numeric(),
    slope2 = numeric(),
    SE1 = numeric(),
    SE2 = numeric(),
    Tvalue1 = numeric(),
    Tvalue2 = numeric(),
    Pvalue1 = numeric(),
    Pvalue2 = numeric(),
    Rsquare = numeric(),
    AdjRsquare = numeric(),
    #Count = numeric(),
    n = numeric(),
    Residual = numeric()
  )

for (i in unique(df_onset$species)) {
  print(i)
  df8b <- subset(df_onset, species == i)
  for (j in unique(df8b$plot)) {
    df8a <- subset(df8b, plot == j)
    
    if (sum(!is.na(df8a$trend)) < 6) {
      df_other <- data.frame(
        species = df8a$species[1],
        plot = df8a$plot[1],
        pheno_event = "Onset",
        slope1 = NA,
        slope2 = NA,
        SE1 = NA,
        SE2 = NA,
        Tvalue1 = NA,
        Tvalue2 = NA,
        Pvalue1 = NA,
        Pvalue2 = NA,
        Rsquare = NA,
        AdjRsquare = NA,
        #Count = NA,
        n = NA,
        Residual = NA
      )
    }
    else{
      mod1 <- lm(trend ~ Snowmelt + Temperature, data = df8a, weights = inv_std_err / mean(inv_std_err))
      Residual1 <- sqrt(deviance(mod1) / df.residual(mod1))
      df_other <- data.frame(
        species = df8a$species[1],
        plot = df8a$plot[1],
        pheno_event = "Onset",
        slope1 = summary(mod1)$coefficients[2],
        slope2 = summary(mod1)$coefficients[3],
        SE1 = summary(mod1)$coefficients[5],
        SE2 = summary(mod1)$coefficients[6],
        Tvalue1 = summary(mod1)$coefficients[8],
        Tvalue2 = summary(mod1)$coefficients[9],
        Pvalue1 = summary(mod1)$coefficients[11],
        Pvalue2 = summary(mod1)$coefficients[12],
        Rsquare = summary(mod1)$r.squared,
        AdjRsquare = summary(mod1)$adj.r.squared,
        #Count = sum(df8a$TotalAbundance),
        n = sum(!is.na(df8a$trend)),
        Residual = Residual1
      )
      df_summary_onset <- bind_rows(df_summary_onset, df_other)
    }
    #plot(mod1)
    
  }
}


df_summary_onset$body_size <- c(5.7, 6, 6, 6, 5.6, 6.05, 6.05, 6.05, 5.1, 3.85, 4.15, 4.25, 4.25, 4.25, 5.25, 5.05)
df_summary_onset$flower_visiting <- c("Irregular", "Frequent", "Frequent", "Frequent", "Frequent","Frequent", "Frequent","Frequent","Frequent", "Irregular", "NA", "Frequent", "Frequent", "Frequent", "Irregular", "Irregular")


df_mean_onset <- read.csv("Data/phenology_estimates/df_mean_onset.csv",sep=",",
                          stringsAsFactors = FALSE, header = TRUE)


df_summary_onset$average_doy <- df_mean_onset$average_doy[match(paste0(df_summary_onset$species),
                                                                paste0(df_mean_onset$species))]


########PEAK#############

df_peak <- read.csv("Data/phenology_estimates/peak_temp_trend_nas_removed.csv",sep=",",
                     stringsAsFactors = FALSE, header = TRUE)

# Get climate data

df_temp <- read.csv("Data/Species_data_clean/All_species_cleanup_peak.csv",sep=",",
                    stringsAsFactors = FALSE, header = TRUE)

df_temp$Species_name[df_temp$Species_name == "Dry.gro"] <- "D.groenlandica"
df_temp$Species_name[df_temp$Species_name == "Dry.seg"] <- "D.segnis"
df_temp$Species_name[df_temp$Species_name == "Lim.gro"] <- "L.groenlandica"
df_temp$Species_name[df_temp$Species_name == "Pha.bid"] <- "P.bidentata"
df_temp$Species_name[df_temp$Species_name == "Spi.alm"] <- "S.almqvistii"
df_temp$Species_name[df_temp$Species_name == "Spi.def"] <- "S.deflorata"
df_temp$Species_name[df_temp$Species_name == "Spi.dor"] <- "S.dorsata"
df_temp$Species_name[df_temp$Species_name == "Spi.mal"] <- "S.malaisei"
df_temp$Species_name[df_temp$Species_name == "Spi.mel"] <- "S.melanosoma"
df_temp$Species_name[df_temp$Species_name == "Spi.nov"] <- "S.novaesibiriae"
df_temp$Species_name[df_temp$Species_name == "Spi.san"] <- "S.sanctipauli"
df_temp$Species_name[df_temp$Species_name == "Spi.tun"] <- "S.tundrae"
df_temp$Species_name[df_temp$Species_name == "Spi.zai"] <- "S.zaitzevi"



df_peak$Temperature <- df_temp$Temperature[match(paste0(df_peak$year,df_peak$species,df_peak$plot),
                                                 paste0(df_temp$Year,df_temp$Species_name,df_temp$Plot))]

df_peak$Snowmelt <- df_temp$SnowmeltDOY[match(paste0(df_peak$year,df_peak$species,df_peak$plot),
                                              paste0(df_temp$Year,df_temp$Species_name,df_temp$Plot))]

df_peak <- df_peak |>
  mutate(inv_std_err = 1/se)


df_summary_peak <-
  data.frame(
    species = character(),
    plot = character(),
    pheno_event = character(),
    slope1 = numeric(),
    slope2 = numeric(),
    SE1 = numeric(),
    SE2 = numeric(),
    Tvalue1 = numeric(),
    Tvalue2 = numeric(),
    Pvalue1 = numeric(),
    Pvalue2 = numeric(),
    Rsquare = numeric(),
    AdjRsquare = numeric(),
    #Count = numeric(),
    n = numeric(),
    Residual = numeric()
  )

for (i in unique(df_peak$species)) {
  print(i)
  df8b <- subset(df_peak, species == i)
  for (j in unique(df8b$plot)) {
    df8a <- subset(df8b, plot == j)
    
    if (sum(!is.na(df8a$trend)) < 6) {
      df_other <- data.frame(
        species = df8a$species[1],
        plot = df8a$plot[1],
        pheno_event = "Peak",
        slope1 = NA,
        slope2 = NA,
        SE1 = NA,
        SE2 = NA,
        Tvalue1 = NA,
        Tvalue2 = NA,
        Pvalue1 = NA,
        Pvalue2 = NA,
        Rsquare = NA,
        AdjRsquare = NA,
        #Count = NA,
        n = NA,
        Residual = NA
      )
    }
    else{
      mod1 <- lm(trend ~ Snowmelt + Temperature, data = df8a, weights = inv_std_err / mean(inv_std_err))
      Residual1 <- sqrt(deviance(mod1) / df.residual(mod1))
      df_other <- data.frame(
        species = df8a$species[1],
        plot = df8a$plot[1],
        pheno_event = "Peak",
        slope1 = summary(mod1)$coefficients[2],
        slope2 = summary(mod1)$coefficients[3],
        SE1 = summary(mod1)$coefficients[5],
        SE2 = summary(mod1)$coefficients[6],
        Tvalue1 = summary(mod1)$coefficients[8],
        Tvalue2 = summary(mod1)$coefficients[9],
        Pvalue1 = summary(mod1)$coefficients[11],
        Pvalue2 = summary(mod1)$coefficients[12],
        Rsquare = summary(mod1)$r.squared,
        AdjRsquare = summary(mod1)$adj.r.squared,
        #Count = sum(df8a$TotalAbundance),
        n = sum(!is.na(df8a$trend)),
        Residual = Residual1
      )
      df_summary_peak <- bind_rows(df_summary_peak, df_other)
    }
    #plot(mod1)
    
  }
}

df_summary_peak$body_size <- c(5.7, 6, 6, 6, 5.6, 6.05, 6.05, 6.05, 5.1, 3.85, 4.15, 4.25, 4.25, 4.25, 5.25, 5.05)
df_summary_peak$flower_visiting <- c("Irregular", "Frequent", "Frequent", "Frequent", "Frequent","Frequent", "Frequent","Frequent","Frequent", "Irregular", "NA", "Frequent", "Frequent", "Frequent", "Irregular", "Irregular")

df_mean_peak <- read.csv("Data/phenology_estimates/df_mean_peak.csv",sep=",",
                         stringsAsFactors = FALSE, header = TRUE)

df_summary_peak$average_doy <- df_mean_peak$average_doy[match(paste0(df_summary_peak$species),
                                                              paste0(df_mean_peak$species))]


##############END###################

df_end <- read.csv("Data/phenology_estimates/end_temp_trend_nas_removed.csv",sep=",",
                     stringsAsFactors = FALSE, header = TRUE)

df_end$species[df_end$species == "Dry.gro"] <- "D.groenlandica"
df_end$species[df_end$species == "Dry.seg"] <- "D.segnis"
df_end$species[df_end$species == "Lim.gro"] <- "L.groenlandica"
df_end$species[df_end$species == "Pha.bid"] <- "P.bidentata"
df_end$species[df_end$species == "Spi.alm"] <- "S.almqvistii"
df_end$species[df_end$species == "Spi.def"] <- "S.deflorata"
df_end$species[df_end$species == "Spi.dor"] <- "S.dorsata"
df_end$species[df_end$species == "Spi.mal"] <- "S.malaisei"
df_end$species[df_end$species == "Spi.mel"] <- "S.melanosoma"
df_end$species[df_end$species == "Spi.nov"] <- "S.novaesibiriae"
df_end$species[df_end$species == "Spi.san"] <- "S.sanctipauli"
df_end$species[df_end$species == "Spi.tun"] <- "S.tundrae"
df_end$species[df_end$species == "Spi.zai"] <- "S.zaitzevi"

# Get climate data

df_temp <- read.csv("Data/Species_data_clean/All_species_cleanup_end.csv",sep=",",
                    stringsAsFactors = FALSE, header = TRUE)

df_temp$Species_name[df_temp$Species_name == "Dry.gro"] <- "D.groenlandica"
df_temp$Species_name[df_temp$Species_name == "Dry.seg"] <- "D.segnis"
df_temp$Species_name[df_temp$Species_name == "Lim.gro"] <- "L.groenlandica"
df_temp$Species_name[df_temp$Species_name == "Pha.bid"] <- "P.bidentata"
df_temp$Species_name[df_temp$Species_name == "Spi.alm"] <- "S.almqvistii"
df_temp$Species_name[df_temp$Species_name == "Spi.def"] <- "S.deflorata"
df_temp$Species_name[df_temp$Species_name == "Spi.dor"] <- "S.dorsata"
df_temp$Species_name[df_temp$Species_name == "Spi.mal"] <- "S.malaisei"
df_temp$Species_name[df_temp$Species_name == "Spi.mel"] <- "S.melanosoma"
df_temp$Species_name[df_temp$Species_name == "Spi.nov"] <- "S.novaesibiriae"
df_temp$Species_name[df_temp$Species_name == "Spi.san"] <- "S.sanctipauli"
df_temp$Species_name[df_temp$Species_name == "Spi.tun"] <- "S.tundrae"
df_temp$Species_name[df_temp$Species_name == "Spi.zai"] <- "S.zaitzevi"



df_end$Temperature <- df_temp$Temperature[match(paste0(df_end$year,df_end$species,df_end$plot),
                                                paste0(df_temp$Year,df_temp$Species_name,df_temp$Plot))]

df_end$Snowmelt <- df_temp$SnowmeltDOY[match(paste0(df_end$year,df_end$species,df_end$plot),
                                             paste0(df_temp$Year,df_temp$Species_name,df_temp$Plot))]

df_end <- df_end |>
  mutate(inv_std_err = 1/se)

df_summary_end <-
  data.frame(
    species = character(),
    plot = character(),
    pheno_event = character(),
    slope1 = numeric(),
    slope2 = numeric(),
    SE1 = numeric(),
    SE2 = numeric(),
    Tvalue1 = numeric(),
    Tvalue2 = numeric(),
    Pvalue1 = numeric(),
    Pvalue2 = numeric(),
    Rsquare = numeric(),
    AdjRsquare = numeric(),
    #Count = numeric(),
    n = numeric(),
    Residual = numeric()
  )

for (i in unique(df_end$species)) {
  print(i)
  df8b <- subset(df_end, species == i)
  for (j in unique(df8b$plot)) {
    df8a <- subset(df8b, plot == j)
    
    if (sum(!is.na(df8a$trend)) < 6) {
      df_other <- data.frame(
        species = df8a$species[1],
        plot = df8a$plot[1],
        pheno_event = "End",
        slope1 = NA,
        slope2 = NA,
        SE1 = NA,
        SE2 = NA,
        Tvalue1 = NA,
        Tvalue2 = NA,
        Pvalue1 = NA,
        Pvalue2 = NA,
        Rsquare = NA,
        AdjRsquare = NA,
        #Count = NA,
        n = NA,
        Residual = NA
      )
    }
    else{
      mod1 <- lm(trend ~ Snowmelt + Temperature, data = df8a, weights = inv_std_err / mean(inv_std_err))
      Residual1 <- sqrt(deviance(mod1) / df.residual(mod1))
      df_other <- data.frame(
        species = df8a$species[1],
        plot = df8a$plot[1],
        pheno_event = "End",
        slope1 = summary(mod1)$coefficients[2],
        slope2 = summary(mod1)$coefficients[3],
        SE1 = summary(mod1)$coefficients[5],
        SE2 = summary(mod1)$coefficients[6],
        Tvalue1 = summary(mod1)$coefficients[8],
        Tvalue2 = summary(mod1)$coefficients[9],
        Pvalue1 = summary(mod1)$coefficients[11],
        Pvalue2 = summary(mod1)$coefficients[12],
        Rsquare = summary(mod1)$r.squared,
        AdjRsquare = summary(mod1)$adj.r.squared,
        #Count = sum(df8a$TotalAbundance),
        n = sum(!is.na(df8a$trend)),
        Residual = Residual1
      )
      df_summary_end <- bind_rows(df_summary_end, df_other)
    }
    #plot(mod1)
    
  }
}


df_summary_end$body_size <- c(5.7, 6, 6, 6, 5.6, 6.05, 6.05, 6.05, 5.1, 3.85, 4.15, 4.25, 4.25, 4.25, 5.25, 5.05)

df_summary_end$flower_visiting <- c("Irregular", "Frequent", "Frequent", "Frequent", "Frequent","Frequent", "Frequent","Frequent","Frequent", "Irregular", "NA", "Frequent", "Frequent", "Frequent", "Irregular", "Irregular")


df_mean_end <- read.csv("Data/phenology_estimates/df_mean_end.csv",sep=",",
                          stringsAsFactors = FALSE, header = TRUE)


df_summary_end$average_doy <- df_mean_end$average_doy_end[match(paste0(df_summary_end$species),
                                                                paste0(df_mean_end$species))]


############SNOWMELT##################

onset_size <- ggplot(df_summary_onset, aes(body_size, slope1, fill = species, shape = plot), color = "black")+
  xlim(3.5,6.5) +
  ylim(-1,1)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(1, 0.1, 0.5, 0), "cm"))

lm.size <- lm(slope1 ~ body_size, df_summary_onset)

summary(lm.size)

peak_size <- ggplot(df_summary_peak, aes(body_size, slope1, fill = species, shape = plot), color = "black")+
  xlim(3.5,6.5) +
  ylim(-1,1)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  scale_shape_manual(values = c(21, 22, 24)) +
  xlab("Wing length (mm)")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(1, 0.1, 0.5, 0), "cm"))

lm.size <- lm(slope1 ~ body_size, df_summary_peak)

summary(lm.size)


end_size <- ggplot(df_summary_end, aes(body_size, slope1, fill = species, shape = plot), color = "black")+
  xlim(3.5,6.5) +
  ylim(-1,1)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(1, 0.1, 0.5, 0), "cm"))

lm.size <- lm(slope1 ~ body_size, df_summary_end)

summary(lm.size)

body_size <- ggarrange(onset_size, peak_size, end_size, labels = c("A", "B", "C"), hjust = -3, vjust = 2,  nrow = 1) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


#phenological niche

onset_niche <- ggplot(df_summary_onset, aes(average_doy, slope1, fill = species, shape = plot), color = "black")+
  ylim(-1,1)+
  xlim(179,192)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #scale_x_continuous(breaks = c(176, 178, 180, 182, 184, 186, 188, 190, 192, 194)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0, 0.1, 1.5, 0), "cm"))

lm.niche <- lm(slope1 ~ average_doy, df_summary_onset)

summary(lm.niche)

peak_niche <- ggplot(df_summary_peak, aes(average_doy, slope1, fill = species, shape = plot), color = "black")+
  ylim(-1,1)+
  xlim(188,212)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #scale_x_continuous(breaks = c(186, 190, 194, 198, 202, 206, 210, 214)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0, 0.1, 1.5, 0), "cm"))

lm.niche <- lm(slope1 ~ average_doy, df_summary_peak)

summary(lm.niche)


end_niche <- ggplot(df_summary_end, aes(average_doy, slope1, fill = species, shape = plot), color = "black")+
  ylim(-1,1)+
  xlim(200,232)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #scale_x_continuous(breaks = c(200, 208, 216, 224, 232, 240)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0, 0.1, 1.5, 0), "cm"))

lm.niche <- lm(slope1 ~ average_doy, df_summary_end)

summary(lm.niche)

niche <- ggarrange(onset_niche, peak_niche, end_niche, labels = c("G", "H", "I"), hjust = c(-3, -3, -7), vjust = -1.5, nrow = 1) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

#flower visitors

df_summary_onset <- subset(df_summary_onset, !df_summary_onset$species == "S.novaesibiriae")
df_summary_peak <- subset(df_summary_peak, !df_summary_peak$species == "S.novaesibiriae")
df_summary_end <- subset(df_summary_end, !df_summary_end$species == "S.novaesibiriae")


onset_visit <- ggplot(df_summary_onset, aes(flower_visiting, slope1, fill = species, shape = plot), color = "black")+
  ylim(-1,1)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("Phenological shifts (days/shifted snowmelt day)")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0.5, 0.1, 1, 0), "cm"))

lm.visit <- lm(slope1 ~ flower_visiting, df_summary_onset)

summary(lm.visit)

peak_visit <- ggplot(df_summary_peak, aes(flower_visiting, slope1, fill = species, shape = plot), color = "black")+
  ylim(-1,1)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0.5, 0.1, 1, 0), "cm"))

lm.visit <- lm(slope1 ~ flower_visiting, df_summary_peak)

summary(lm.visit)


end_visit <- ggplot(df_summary_end, aes(flower_visiting, slope1, fill = species, shape = plot), color = "black")+
  ylim(-1,1)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0.5, 0.1, 1, 0), "cm"))

lm.visit <- lm(slope1 ~ flower_visiting, df_summary_end)

summary(lm.visit)

visit <- ggarrange(onset_visit, peak_visit, end_visit, labels = c("D", "E", "F"), hjust = -3, vjust = 0.3,  nrow = 1) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))



ggarrange(body_size, visit, niche, nrow = 3) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))


############TEMPERATURE##################

onset_size <- ggplot(df_summary_onset, aes(body_size, slope2, fill = species, shape = plot), color = "black")+
  xlim(3.5,6.5) +
  ylim(-13,3)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(1, 0.1, 0.5, 0), "cm"))

lm.size <- lm(slope2 ~ body_size, df_summary_onset)

summary(lm.size)

peak_size <- ggplot(df_summary_peak, aes(body_size, slope2, fill = species, shape = plot), color = "black")+
  xlim(3.5,6.5) +
  ylim(-13,3)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  scale_shape_manual(values = c(21, 22, 24)) +
  xlab("Wing length (mm)")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(1, 0.1, 0.5, 0), "cm"))

lm.size <- lm(slope2 ~ body_size, df_summary_peak)

summary(lm.size)


end_size <- ggplot(df_summary_end, aes(body_size, slope2, fill = species, shape = plot), color = "black")+
  xlim(3.5,6.5) +
  ylim(-13,3)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(1, 0.1, 0.5, 0), "cm"))

lm.size <- lm(slope2 ~ body_size, df_summary_end)

summary(lm.size)

body_size <- ggarrange(onset_size, peak_size, end_size, labels = c("A", "B", "C"), hjust = -3, vjust = 2,  nrow = 1) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


#phenological niche

onset_niche <- ggplot(df_summary_onset, aes(average_doy, slope2, fill = species, shape = plot), color = "black")+
  ylim(-13,3)+
  xlim(179,192)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #scale_x_continuous(breaks = c(176, 178, 180, 182, 184, 186, 188, 190, 192, 194)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0, 0.1, 1.5, 0), "cm"))

lm.niche <- lm(slope2 ~ average_doy, df_summary_onset)

summary(lm.niche)

peak_niche <- ggplot(df_summary_peak, aes(average_doy, slope2, fill = species, shape = plot), color = "black")+
  ylim(-13,3)+
  xlim(188,212)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #scale_x_continuous(breaks = c(186, 190, 194, 198, 202, 206, 210, 214)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0, 0.1, 1.5, 0), "cm"))

lm.niche <- lm(slope2 ~ average_doy, df_summary_peak)

summary(lm.niche)


end_niche <- ggplot(df_summary_end, aes(average_doy, slope2, fill = species, shape = plot), color = "black")+
  ylim(-13,3)+
  xlim(200,232)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #scale_x_continuous(breaks = c(200, 208, 216, 224, 232, 240)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0, 0.1, 1.5, 0), "cm"))

lm.niche <- lm(slope2 ~ average_doy, df_summary_end)

summary(lm.niche)

niche <- ggarrange(onset_niche, peak_niche, end_niche, labels = c("G", "H", "I"), hjust = c(-3, -3, -7), vjust = -1.5, nrow = 1) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

#flower visitors

df_summary_onset <- subset(df_summary_onset, !df_summary_onset$species == "S.novaesibiriae")
df_summary_peak <- subset(df_summary_peak, !df_summary_peak$species == "S.novaesibiriae")
df_summary_end <- subset(df_summary_end, !df_summary_end$species == "S.novaesibiriae")


onset_visit <- ggplot(df_summary_onset, aes(flower_visiting, slope2, fill = species, shape = plot), color = "black")+
  ylim(-13,3)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("Phenological shifts (days/shifted snowmelt day)")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0.5, 0.1, 1, 0), "cm"))

lm.visit <- lm(slope2 ~ flower_visiting, df_summary_onset)

summary(lm.visit)

peak_visit <- ggplot(df_summary_peak, aes(flower_visiting, slope2, fill = species, shape = plot), color = "black")+
  ylim(-13,3)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0.5, 0.1, 1, 0), "cm"))

lm.visit <- lm(slope2 ~ flower_visiting, df_summary_peak)

summary(lm.visit)


end_visit <- ggplot(df_summary_end, aes(flower_visiting, slope2, fill = species, shape = plot), color = "black")+
  ylim(-13,3)+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = "dashed") +
  geom_point(size = 4)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  theme_test()+
  theme(axis.title.x = element_text(size = 12, color = "black", vjust = -1),
        axis.title.y = element_text(size = 12, color = "black", vjust = 2),
        axis.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0.5, 0.1, 1, 0), "cm"))

lm.visit <- lm(slope2 ~ flower_visiting, df_summary_end)

summary(lm.visit)

visit <- ggarrange(onset_visit, peak_visit, end_visit, labels = c("D", "E", "F"), hjust = -3, vjust = 0.3,  nrow = 1) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))



ggarrange(body_size, visit, niche, nrow = 3) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))




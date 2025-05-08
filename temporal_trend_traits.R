# Temporal trends based on traits

library(readxl)
library(ggpubr)

df_summary_onset <- read_xlsx("Data/phenology_estimates/df_summary_onset.xlsx")
df_summary_peak <- read_xlsx("Data/phenology_estimates/df_summary_peak.xlsx")
df_summary_end <- read_xlsx("Data/phenology_estimates/df_summary_end.xlsx")

#Remove years with less than 5 years

df_summary_onset <- subset(df_summary_onset, df_summary_onset$plot!= "Art5" | df_summary_onset$species!= "S.dorsata")
df_summary_peak <- subset(df_summary_peak, df_summary_peak$plot!= "Art5" | df_summary_peak$species!= "S.dorsata")
df_summary_end <- subset(df_summary_end, df_summary_end$plot!= "Art5" | df_summary_end$species!= "S.dorsata")

df_summary_onset <- subset(df_summary_onset, df_summary_onset$plot!= "Art3" | df_summary_onset$species!= "D.groenlandica")
df_summary_peak <- subset(df_summary_peak, df_summary_peak$plot!= "Art3" | df_summary_peak$species!= "D.groenlandica")
df_summary_end <- subset(df_summary_end, df_summary_end$plot!= "Art3" | df_summary_end$species!= "D.groenlandica")


df_summary_onset$slope_decade <- df_summary_onset$slope*10
df_summary_peak$slope_decade <- df_summary_peak$slope*10
df_summary_end$slope_decade <- df_summary_end$slope*10

df_summary_onset$se_decade <- df_summary_onset$SE*10
df_summary_peak$se_decade <- df_summary_peak$SE*10
df_summary_end$se_decade <- df_summary_end$SE*10

df_mean_onset <- read.csv("Data/phenology_estimates/df_mean_onset.csv",sep=",",
                        stringsAsFactors = FALSE, header = TRUE)

df_mean_end <- read.csv("Data/phenology_estimates/df_mean_end.csv",sep=",",
                        stringsAsFactors = FALSE, header = TRUE)

df_mean_peak <- read.csv("Data/phenology_estimates/df_mean_peak.csv",sep=",",
                         stringsAsFactors = FALSE, header = TRUE)

df_summary_onset$average_doy <- df_mean_onset$average_doy[match(paste0(df_summary_onset$species),
                                                paste0(df_mean_onset$species))]

df_summary_peak$average_doy <- df_mean_peak$average_doy[match(paste0(df_summary_peak$species),
                                                                paste0(df_mean_peak$species))]

df_summary_end$average_doy <- df_mean_end$average_doy_end[match(paste0(df_summary_end$species),
                                                                paste0(df_mean_end$species))]


onset_size <- ggplot(df_summary_onset, aes(body_size, slope_decade, fill = species, shape = plot, 
                                           color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlim(3.5,6.5) +
  #ylim(-40,25)+
  coord_cartesian(ylim=c(-40, 25))+
  geom_errorbar(mapping = aes(ymin = slope_decade - se_decade, ymax = slope_decade + se_decade, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))), width = 0, linewidth = 0.5) + # Adding error bars
  geom_point(size = 3, stroke = 0.6)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("Phenological shifts (days/decade)")+
  scale_fill_viridis_d(option = "plasma")+
  scale_color_manual(values = c('grey', 'grey48', 'black'),
                     limits = c('(0,0.05]', '(0.05,0.099]', '(0.099,1]'))+
  theme_test()+
  theme(axis.title.x = element_text(size = 10, color = "black", vjust = -1),
        axis.title.y = element_text(size = 10, color = "black", vjust = 2),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(1, 0.1, 0.5, 0), "cm"))

lm.size <- lm(slope ~ body_size, df_summary_onset)

summary(lm.size)

peak_size <- ggplot(df_summary_peak, aes(body_size, slope_decade, fill = species, shape = plot, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlim(3.5,6.5) +
  coord_cartesian(ylim=c(-40, 25)) +
  geom_errorbar(mapping = aes(ymin = slope_decade - se_decade, ymax = slope_decade + se_decade, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))), width = 0, linewidth = 0.5) + # Adding error bars
  geom_point(size = 3, stroke = 0.6)+
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  scale_shape_manual(values = c(21, 22, 24)) +
  xlab("Wing length (mm)")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  scale_color_manual(values = c('grey', 'grey48', 'black'),
                     limits = c('(0,0.05]', '(0.05,0.099]', '(0.099,1]'))+
  theme_test()+
  theme(axis.title.x = element_text(size = 10, color = "black", vjust = -1),
        axis.title.y = element_text(size = 10, color = "black", vjust = 2),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(1, 0.1, 0.5, 0), "cm"))

lm.size <- lm(slope ~ body_size, df_summary_peak)

summary(lm.size)


end_size <- ggplot(df_summary_end, aes(body_size, slope_decade))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_shape_manual(values = c(21, 22, 24)) +
  xlim(3.5,6.5) +
  coord_cartesian(ylim=c(-40, 25))+
  geom_errorbar(mapping = aes(ymin = slope_decade - se_decade, ymax = slope_decade + se_decade, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))), width = 0, linewidth = 0.5) + # Adding error bars
  geom_point(df_summary_end, mapping = aes(body_size, slope_decade, fill = species, shape = plot, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))), size = 3, stroke = 0.6)+
  geom_smooth(df_summary_end, mapping = aes(body_size, slope_decade), method = "lm", se = FALSE, color = "black")+
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  scale_color_manual(values = c('grey', 'grey48', 'black'),
                     limits = c('(0,0.05]', '(0.05,0.099]', '(0.099,1]'))+
  theme_test()+
  theme(axis.title.x = element_text(size = 10, color = "black", vjust = -1),
        axis.title.y = element_text(size = 10, color = "black", vjust = 2),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(1, 0.1, 0.5, 0), "cm"))



#df_summary_end <- df_summary_end %>% 
  #filter(!(species == "S.novaesibiriae" & plot == "Art2"))

lm.size <- lm(slope ~ body_size, df_summary_end)

summary(lm.size)

body_size <- ggarrange(onset_size, peak_size, end_size, labels = c("A", "B", "C"), hjust = -3, vjust = 2,  nrow = 1) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


#phenological niche

onset_niche <- ggplot(df_summary_onset, aes(average_doy, slope_decade, fill = species, shape = plot, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim=c(-40, 25))+
  geom_errorbar(mapping = aes(ymin = slope_decade - se_decade, ymax = slope_decade + se_decade, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))), width = 0, linewidth = 0.5) + # Adding error bars
  xlim(179,192)+
  geom_point(size = 3, stroke = 0.6)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #scale_x_continuous(breaks = c(176, 178, 180, 182, 184, 186, 188, 190, 192, 194)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("Phenological shifts (days/decade)")+
  scale_fill_viridis_d(option = "plasma")+
  scale_color_manual(values = c('grey', 'grey48', 'black'),
                     limits = c('(0,0.05]', '(0.05,0.099]', '(0.099,1]'))+
  theme_test()+
  theme(axis.title.x = element_text(size = 10, color = "black", vjust = -1),
        axis.title.y = element_text(size = 10, color = "black", vjust = 2),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0, 0.1, 1.5, 0), "cm"))

lm.niche <- lm(slope ~ average_doy, df_summary_onset)

summary(lm.niche)

peak_niche <- ggplot(df_summary_peak, aes(average_doy, slope_decade, fill = species, shape = plot, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim=c(-40, 25))+
  geom_errorbar(mapping = aes(ymin = slope_decade - se_decade, ymax = slope_decade + se_decade, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))), width = 0, linewidth = 0.5) + # Adding error bars
  xlim(188,212)+
  geom_point(size = 3, stroke = 0.6)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #scale_x_continuous(breaks = c(186, 190, 194, 198, 202, 206, 210, 214)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  scale_color_manual(values = c('grey', 'grey48', 'black'),
                     limits = c('(0,0.05]', '(0.05,0.099]', '(0.099,1]'))+
  theme_test()+
  theme(axis.title.x = element_text(size = 10, color = "black", vjust = -1),
        axis.title.y = element_text(size = 10, color = "black", vjust = 2),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0, 0.1, 1.5, 0), "cm"))

lm.niche <- lm(slope ~ average_doy, df_summary_peak)

summary(lm.niche)


end_niche <- ggplot(df_summary_end, aes(average_doy, slope_decade, fill = species, shape = plot, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim=c(-40, 25))+
  geom_errorbar(mapping = aes(ymin = slope_decade - se_decade, ymax = slope_decade + se_decade, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))), width = 0, linewidth = 0.5) + # Adding error bars
  xlim(200,232)+
  geom_point(size = 3, stroke = 0.6)+
  scale_shape_manual(values = c(21, 22, 24)) +
  #scale_x_continuous(breaks = c(200, 208, 216, 224, 232, 240)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  scale_color_manual(values = c('grey', 'grey48', 'black'),
                     limits = c('(0,0.05]', '(0.05,0.099]', '(0.099,1]'))+
  theme_test()+
  theme(axis.title.x = element_text(size = 10, color = "black", vjust = -1),
        axis.title.y = element_text(size = 10, color = "black", vjust = 2),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0, 0.1, 1.5, 0), "cm"))

lm.niche <- lm(slope ~ pheno_niche, df_summary_end)

summary(lm.niche)

niche <- ggarrange(onset_niche, peak_niche, end_niche, labels = c("G", "H", "I"), hjust = c(-3, -3, -7), vjust = -1.5, nrow = 1) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

#flower visitors

df_summary_onset <- subset(df_summary_onset, !df_summary_onset$species == "S.novaesibiriae")
df_summary_peak <- subset(df_summary_peak, !df_summary_peak$species == "S.novaesibiriae")
df_summary_end <- subset(df_summary_end, !df_summary_end$species == "S.novaesibiriae")


onset_visit <- ggplot(df_summary_onset, aes(flower_visiting, slope_decade, fill = species, shape = plot, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-40, 25)) +
  geom_errorbar(mapping = aes(ymin = slope_decade - se_decade, ymax = slope_decade + se_decade, 
                              color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1)), group = interaction(species, plot, Pvalue)), 
                position = position_dodge(width = 0.7), width = 0, linewidth = 0.5) + # Adding error bars
  geom_point(position = position_dodge(width = 0.7), size = 3, stroke = 0.6, aes(group = interaction(species, plot, Pvalue))) +
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("Phenological shifts (days/decade)")+
  scale_fill_viridis_d(option = "plasma")+
  scale_color_manual(values = c('grey', 'grey48', 'black'),
                     limits = c('(0,0.05]', '(0.05,0.099]', '(0.099,1]'))+
  theme_test()+
  theme(axis.title.x = element_text(size = 10, color = "black", vjust = -1),
        axis.title.y = element_text(size = 10, color = "black", vjust = 2),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0.5, 0.1, 1, 0), "cm"))

lm.visit <- lm(slope ~ flower_visiting, df_summary_onset)

summary(lm.visit)

peak_visit <- ggplot(df_summary_peak, aes(flower_visiting, slope_decade, fill = species, shape = plot, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-40, 25)) +
  geom_errorbar(mapping = aes(ymin = slope_decade - se_decade, ymax = slope_decade + se_decade, 
                              color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1)), group = interaction(species, plot, Pvalue)), 
                position = position_dodge(width = 0.7), width = 0, linewidth = 0.5) + # Adding error bars
  geom_point(position = position_dodge(width = 0.7), size = 3, stroke = 0.6, aes(group = interaction(species, plot, Pvalue))) +
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  scale_color_manual(values = c('grey', 'grey48', 'black'),
                     limits = c('(0,0.05]', '(0.05,0.099]', '(0.099,1]'))+
  theme_test()+
  theme(axis.title.x = element_text(size = 10, color = "black", vjust = -1),
        axis.title.y = element_text(size = 10, color = "black", vjust = 2),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0.5, 0.1, 1, 0), "cm"))

lm.visit <- lm(slope ~ flower_visiting, df_summary_peak)

summary(lm.visit)


end_visit <- ggplot(df_summary_end, aes(flower_visiting, slope_decade, fill = species, shape = plot, color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1))))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-40, 25)) +
  geom_errorbar(mapping = aes(ymin = slope_decade - se_decade, ymax = slope_decade + se_decade, 
                              color = cut(Pvalue, breaks = c(0, 0.05, 0.099, 1)), group = interaction(species, plot, Pvalue)), 
                position = position_dodge(width = 0.7), width = 0, linewidth = 0.5) + # Adding error bars
  geom_point(position = position_dodge(width = 0.7), size = 3, stroke = 0.6, aes(group = interaction(species, plot, Pvalue))) +
  scale_shape_manual(values = c(21, 22, 24)) +
  #geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("")+
  ylab("")+
  scale_fill_viridis_d(option = "plasma")+
  scale_color_manual(values = c('grey', 'grey48', 'black'),
                     limits = c('(0,0.05]', '(0.05,0.099]', '(0.099,1]'))+
  theme_test()+
  theme(axis.title.x = element_text(size = 10, color = "black", vjust = -1),
        axis.title.y = element_text(size = 10, color = "black", vjust = 2),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "none",
        plot.margin=unit(c(0.5, 0.1, 1, 0), "cm"))

lm.visit <- lm(slope ~ flower_visiting, df_summary_end)

summary(lm.visit)

visit <- ggarrange(onset_visit, peak_visit, end_visit, labels = c("D", "E", "F"), hjust = -3, vjust = 0.3,  nrow = 1) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))



ggarrange(body_size, visit, niche, nrow = 3) +
  theme(plot.margin = unit(c(6.5, 4, 6.5, 4), "cm"))

#ggarrange(body_size, nrow = 1) +
 # theme(plot.margin = unit(c(11, 1.5, 11, 1.5), "cm"))








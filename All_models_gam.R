# GAM full model with species and plot
# This model is close to the final model and takes a long time to run

pkgs <- c("tibble", "readr", "dplyr", "here", "tidyr", "ggplot2", "mgcv",
          "gratia", "purrr", "furrr")
vapply(pkgs, library, FUN.VALUE = logical(1L), character.only = TRUE,
       logical.return = TRUE)

df <- read_csv("Data/Species_data_Clean/All_species_clean.csv")

#In the above, data from all species are combined in one dataframe. Initially,
# data that have passed the following criteria: 1) Abundance > 3 in a season,
# 2) The species was present for at least 2 weeks.
# Include = 1 fulfill these criteria


#Visualize abundance data across seasons for each years

df %>%
  ggplot(aes(x = DOY, y = Abundance, group = Year, colour = Year)) +
  geom_line(alpha = 0.5)+
  scale_color_viridis_c(option = "plasma") +
  facet_grid(Species_name ~ Plot, scales = "free", space = "free_x")

df %>%
  ggplot(aes(x = DOY, y = Abundance)) +
  geom_point()

# Simplified model

df <- df |>
  mutate(Plot = as.factor(Plot),
         Species_name = as.factor(Species_name),
         fyear = as.factor(df$Year),
         exposure = rep(7, length(Plot)))

m_gls <- bam(Abundance ~ s(DOY) + s(Year) +
               s(Year, Species_name, bs = "fs") +
               ti(DOY, Year) +
               #ti(DOY, fyear, Species_name, bs = c("cr", "re", "re"),
               #ti excludes main effect
               #k = c(10, 18, 13)) +
               s(Plot, bs = "re"),
             data = df, family = nb(), method = "fREML",
             discrete = TRUE, nthreads = 4, drop.unused.levels = FALSE)

#Two things to be aware of:
# 1. We still need to include fyear in the model
# 2. offset(log(exposure)) is included to create weekly doys

m_gls2 <- bam(Abundance ~ s(DOY) + s(Year) + ti(DOY, Year) +
                s(Plot, bs = "re") + s(Species_name, bs = "re") +
                s(DOY, fyear, bs = "fs", xt = "cr", k = 15) +
                #s(DOY, Species_name, bs = "fs") +
                #s(Year, Species_name, bs = "fs") +
                #s(DOY, Plot, bs = "fs") +
                #s(Year, Plot, bs = "fs") +
                ti(DOY, Year, Plot, Species_name,
                   bs = c("cr", "cr", "re", "re"), k = c(10, 18, 3, 13)) +
                offset(log(exposure)),
              data = df, family = nb(), method = "fREML", discrete = TRUE,
              nthreads = 4, control = gam.control(trace = TRUE),
drop.unused.levels = FALSE)


#new model 27.10.2023
m_gls3 <- bam(Abundance ~ s(DOY) + s(DOY, Species_name, bs = "sz") +
                s(DOY, Plot, bs = "sz") + 
                s(Year, k = 15) + s(Year, Species_name, bs = "sz") +
                s(Year, Plot, bs = "sz") +
                #s(fyear, bs = "re") +
                offset(log(exposure)),
                #s(DOY, fyear, bs = "fs", xt = "cc", k = 15) +
                #s(DOY, Species_name, bs = "fs") +
                #s(Year, Species_name, bs = "fs") +
                #s(DOY, Plot, bs = "fs") +
                #s(Year, Plot, bs = "fs") +
                #ti(DOY, Year, Plot, Species_name,
                   #bs = c("cr", "cr", "re", "re"), k = c(10, 18, 3, 13)) +
              data = df, family = nb(), method = "fREML", discrete = TRUE,
              nthreads = 3, control = gam.control(trace = TRUE),
drop.unused.levels = FALSE)

draw(m_gls3, residuals = FALSE, rug = FALSE, overall_uncertainty = FALSE)
summary(m_gls3)
appraise(m_gls2, method = "simulate")

sz_bs <- list(bs = "cr")
m_gls4 <- bam(Abundance ~ s(DOY, bs = "cr") +
  s(DOY, Species_name, bs = "sz", xt = sz_bs) +
  s(DOY, Plot, bs = "sz", xt = sz_bs) +
  s(Year, k = 10, bs = "cr") +
  s(Year, Species_name, bs = "sz", xt = sz_bs) +
  s(Year, Plot, bs = "sz", xt = sz_bs) +
  # ti(DOY, Year, bs = c("cr", "cr")) +
  s(DOY, fyear, bs = "fs", xt = "cc", k = 15) +
  offset(log(exposure)),
# s(DOY, fyear, bs = "fs", xt = "cc", k = 15) +
# s(DOY, Species_name, bs = "fs") +
# s(Year, Species_name, bs = "fs") +
# s(DOY, Plot, bs = "fs") +
# s(Year, Plot, bs = "fs") +
# ti(DOY, Year, Plot, Species_name,
# bs = c("cr", "cr", "re", "re"), k = c(10, 18, 3, 13)) +
data = df, family = nb(), method = "fREML", discrete = TRUE,
nthreads = 3, control = gam.control(trace = TRUE),
drop.unused.levels = FALSE)

m_gls5 <- bam(Abundance ~ #s(DOY, bs = "cr") +
  #s(DOY, Species_name, bs = "sz", xt = sz_bs) +
  #s(DOY, Plot, bs = "sz", xt = sz_bs) +
  #s(DOY, Plot, Species_name, bs = "sz", xt = sz_bs) +
  #s(Year, k = 10, bs = "cr") +
  #s(Year, Species_name, bs = "sz", xt = sz_bs) +
  #s(Year, Plot, bs = "sz", xt = sz_bs) +
  #s(Year, Plot, Species_name, bs = "sz", xt = sz_bs) +
  # ti(DOY, Year, bs = c("cr", "cr")) +
  # s(DOY, fyear, bs = "fs", xt = "cc", k = 15) +
  offset(log(exposure)) +
# s(DOY, fyear, bs = "fs", xt = "cc", k = 15) +
# s(DOY, Species_name, bs = "fs") +
# s(Year, Species_name, bs = "fs") +
# s(DOY, Plot, bs = "fs") +
# s(Year, Plot, bs = "fs") +
t2(DOY, fyear, Plot, Species_name, bs = c("cr", "re", "re", "re"), k = c(10, 18, 3, 13)),
data = df, family = nb(), method = "fREML", discrete = TRUE,
nthreads = 3, control = gam.control(trace = TRUE),
drop.unused.levels = FALSE)

# save all the models
write_rds(m_gls,  "m_gls.rds",  compress = "xz")
write_rds(m_gls2, "m_gls2.rds", compress = "xz")
write_rds(m_gls3, "m_gls3.rds", compress = "xz")
write_rds(m_gls4, "m_gls4.rds", compress = "xz")
write_rds(m_gls5, "m_gls5.rds", compress = "xz")

# read in all models
m_gls  <- read_rds("m_gls.rds")
m_gls2 <- read_rds("m_gls2.rds")
m_gls3 <- read_rds("m_gls3.rds")
m_gls4 <- read_rds("m_gls4.rds")
m_gls5 <- read_rds("m_gls5.rds")


AIC(m_gls, m_gls2, m_gls3, m_gls4, m_gls5)

summary(m_gls2, re.test = FALSE)
draw(m_gls2)

par( mfrow= c(2,2) )
gam.check(m_gls2)
dev.off()

####calculating phenologial events####
levs <- with(df, levels(fyear))

m_use <- m_gls2

#Prepare a data slice by chosen covariates
ds <- data_slice(m_use, DOY = evenly(DOY),
                 fyear = evenly(fyear, by = 1),
                 Plot = evenly(Plot))%>%
  #mutate(fyear = factor(Year, levels = levs)) %>%
  filter(DOY <= 260) # We filter doys less than 260 as model extrapolates
# an increase in abundance in later doys
# This is only for the alphabetically first species.
# To get all species include Species = evenly(Species_name)
# Include all species, the output is too big to be retrieved

fv <- fitted_values(m_use, data = ds, scale = "response", n.threads = 4)
# Fitted values from model but only from dataslice created (ds)

fv |>
  mutate(year = as.numeric(as.character(fyear))) |>
  ggplot(aes(x = DOY, y = .fitted, group = year, colour = year)) +
  geom_line(linewidth = 1) +
  scale_colour_viridis_c(option = "cividis") +
  facet_grid(Species_name ~ Plot, scales = "free", space = "free_x") +
  labs(y = "Abundance", x = "DOY",
       title = "Annual seasonal trend in abundance of Muscidae species",
       colour = NULL)

#Check with all species and plots
ds_check <- data_slice(m_use, DOY = evenly(DOY, upper = 260),
                       Year = evenly(Year, by = 1),
                       Plot = evenly(Plot),
                       Species_name = evenly(Species_name)) %>%
  mutate(fyear = factor(Year, levels = levs)) %>%
  filter(DOY <= 260)

fv_check <- fitted_values(m_use, data = ds_check, scale = "response",
  n.threads = 4)

d_segn <- fv_check |>
  filter(Species_name == "Dry.seg") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line() +
  ylab("Fitted values") +
  xlab("Day of Year") +
  labs(title = "D. segnis") +
  geom_point(data = df |> filter(Species_name == "Dry.seg"),
             aes(x = DOY, y = Abundance), size = 1) +
  facet_grid(Year ~ Plot) +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(1,2,1,1, unit = "cm"),
        legend.position = "none")


d_gro <- fv_check |>
  filter(Species_name == "Dry.gro") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line() +
  coord_cartesian(ylim=c(0, 20))+
  ylab("Fitted values") +
  xlab("Day of Year") +
  labs(title = "D. groenlandica") +
  geom_point(data = df |> filter(Species_name == "Dry.gro"),
             aes(x = DOY, y = Abundance), size = 1) +
  facet_grid(Year ~ Plot) +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(1,1,1,0, unit = "cm"))

ggarrange(d_segn,d_gro)

p_bid <- fv_check |>
  filter(Species_name == "Pha.bid") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line() +
  ylab("Fitted values") +
  xlab("Day of Year") +
  labs(title = "P. bidentata") +
  coord_cartesian(ylim=c(0, 20))+
  geom_point(data = df |> filter(Species_name == "Pha.bid"),
             aes(x = DOY, y = Abundance), size = 1) +
  facet_grid(Year ~ Plot) +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(1,2,1,1, unit = "cm"),
        legend.position = "none")


s_alm <- fv_check |>
  filter(Species_name == "Spi.alm") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line() +
  coord_cartesian(ylim=c(0, 150))+
  ylab("Fitted values") +
  xlab("Day of Year") +
  labs(title = "S. almqvistii") +
  geom_point(data = df |> filter(Species_name == "Spi.alm"),
             aes(x = DOY, y = Abundance), size = 1) +
  facet_grid(Year ~ Plot) +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(1,1,1,0, unit = "cm"))

ggarrange(p_bid,s_alm)

s_def <- fv_check |>
  filter(Species_name == "Spi.def") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line() +
  ylab("Fitted values") +
  xlab("Day of Year") +
  labs(title = "S. deflorata") +
  coord_cartesian(ylim=c(0, 20))+
  geom_point(data = df |> filter(Species_name == "Spi.def"),
             aes(x = DOY, y = Abundance), size = 1) +
  facet_grid(Year ~ Plot) +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(1,2,1,1, unit = "cm"),
        legend.position = "none")


s_dor <- fv_check |>
  filter(Species_name == "Spi.dor") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line() +
  coord_cartesian(ylim=c(0, 150))+
  ylab("Fitted values") +
  xlab("Day of Year") +
  labs(title = "S. dorsata") +
  geom_point(data = df |> filter(Species_name == "Spi.dor"),
             aes(x = DOY, y = Abundance), size = 1) +
  facet_grid(Year ~ Plot) +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(1,1,1,0, unit = "cm"))

ggarrange(s_def,s_dor)

s_mal <- fv_check |>
  filter(Species_name == "Spi.mal") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line() +
  ylab("Fitted values") +
  xlab("Day of Year") +
  labs(title = "S. malaisei") +
  coord_cartesian(ylim=c(0, 20))+
  geom_point(data = df |> filter(Species_name == "Spi.mal"),
             aes(x = DOY, y = Abundance), size = 1) +
  facet_grid(Year ~ Plot) +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(1,2,1,1, unit = "cm"),
        legend.position = "none")


s_nov <- fv_check |>
  filter(Species_name == "Spi.nov") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line() +
  coord_cartesian(ylim=c(0, 100))+
  ylab("Fitted values") +
  xlab("Day of Year") +
  labs(title = "S. novaesibiriae") +
  geom_point(data = df |> filter(Species_name == "Spi.nov"),
             aes(x = DOY, y = Abundance), size = 1) +
  facet_grid(Year ~ Plot) +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(1,1,1,0, unit = "cm"))

ggarrange(s_mal,s_nov)

s_san <- fv_check |>
  filter(Species_name == "Spi.san") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line() +
  ylab("Fitted values") +
  xlab("Day of Year") +
  labs(title = "S. sanctipauli") +
  coord_cartesian(ylim=c(0, 200))+
  geom_point(data = df |> filter(Species_name == "Spi.san"),
             aes(x = DOY, y = Abundance), size = 1) +
  facet_grid(Year ~ Plot) +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(1,2,1,1, unit = "cm"),
        legend.position = "none")


s_tun <- fv_check |>
  filter(Species_name == "Spi.tun") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line() +
  coord_cartesian(ylim=c(0, 130))+
  ylab("Fitted values") +
  xlab("Day of Year") +
  labs(title = "S. tundrae") +
  geom_point(data = df |> filter(Species_name == "Spi.tun"),
             aes(x = DOY, y = Abundance), size = 1) +
  facet_grid(Year ~ Plot) +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(1,1,1,0, unit = "cm"))

ggarrange(s_san,s_tun)

fv_check |>
  filter(Species_name == "Spi.zai") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2,
              inherit.aes = FALSE) +
  geom_line() +
  coord_cartesian(ylim=c(0, 200))+
  ylab("Fitted values") +
  xlab("Day of Year") +
  labs(title = "S. zaitzevi") +
  geom_point(data = df |> filter(Species_name == "Spi.zai"),
             aes(x = DOY, y = Abundance), size = 1) +
  facet_grid(Year ~ Plot) +
  theme(plot.title = element_text(face = "italic"),
        plot.margin = margin(1,14,1,1, unit = "cm"))
# Geom_ribbon includes uncertainty bands

# Calculation of peak pheno event

ds_peak <- data_slice(m_use, DOY = evenly(DOY, by = 1),
                       Year = evenly(Year, by = 1),
                       Plot = evenly(Plot),
                       Species_name = evenly(Species_name)) %>%
  mutate(fyear = factor(Year, levels = levs)) %>%
  filter(DOY <= 260)

fv_dry_seg <- fitted_values(m_use,
  data = ds_peak |> filter(Species_name == "Dry.seg"),
  scale = "response", n.threads = 4)

fv_dry_seg |>
  ggplot(aes(x = DOY, y = .fitted, colour = Year, group = Year)) +
  geom_line() +
  scale_colour_viridis_c(option = "cividis") +
  facet_wrap(~ Plot)

fv2 <- fitted_values(m_use, data = ds_check, scale = "response", n.threads = 4)

fv2 |>
  filter(!is.na(.fitted)) |>
  group_by(Species_name, Plot, fyear) |>
  slice_max(order_by = .fitted) |>
  ggplot(aes(x = Year, y = DOY)) +
  geom_point() +
  facet_grid(Species_name ~ Plot, scales = "free", space = "free_x")


# Draw fitted values from the posterior distribution

# generate draws from model of choice
drws <- generate_draws(m_use, n = 100, seed = 42, n_cores = 4)

ds_peak_dry_seg <- ds_peak |> filter(Species_name == "Dry.seg")
fs <- fitted_samples(m_use, data = ds_peak_dry_seg,
  scale = "response", method = "user", draws = drws)

fs <- fs |>
  left_join(ds_peak_dry_seg |>
    mutate(.row = row_number()), by = join_by(.row))

# visulaise a few draws for a selected year
fs |>
  filter(Species_name == "Dry.seg", Year == 1996, .draw %in% 1:10) |>
  ggplot(aes(x = DOY, y = .fitted, group = .draw)) + geom_line() +
  facet_wrap(~Plot) +
  labs(title = "Posterior draws for Dry seg in 1996")

# work through a function that will compute the peak week as an example
peak_wk <- fs |>
  group_by(Species_name, Plot, Year, .draw) |>
  slice_max(order_by = .fitted) |>
  group_by(Species_name, Plot, Year) |>
  summarise(.peak = quantile(DOY, prob = 0.5),
            .lower_ci = quantile(DOY, prob = 0.025),
            .upper_ci = quantile(DOY, prob = 0.975),
            .se = sd(DOY) / sqrt(nrow(fs)),
            .groups = "drop")

# visualise this
peak_wk |>
  mutate(Year = as.numeric(as.character(fyear))) |>
  ggplot(aes(x = Year, y = .peak, group = Plot)) +
  geom_pointrange(aes(ymin = .lower_ci, ymax = .upper_ci)) +
  geom_line() +
  facet_wrap(~ Plot)

ds_peak_row <- ds_peak |>
  mutate(.row = row_number()) |>
  relocate(.row, .before = 1L)

peak_week <- function(data, model, draws) {
  fit <- fitted_samples(model, data = data, scale = "response", method = "user",
    draws = draws)
  data <- data |>
    mutate(.row = row_number())
  fit <- fit |>
    left_join(data, by = join_by(.row))
  fit |>
    group_by(Species_name, Plot, Year, .draw) |>
    slice_max(order_by = .fitted) |>
    group_by(Species_name, Plot, Year) |>
    summarise(.peak = quantile(DOY, prob = 0.5, na.rm = TRUE),
      .lower_ci = quantile(DOY, prob = 0.025, na.rm = TRUE),
      .upper_ci = quantile(DOY, prob = 0.975, na.rm = TRUE),
      .se = sd(DOY) / sqrt(nrow(fit)),
      .groups = "drop")
}

plan(multisession, workers = 3)
op <- options(future.globals.maxSize = 1300 * 1024^2)

pheno <- ds_peak |>
  # filter(Species_name %in% levels(Species_name)[1:3]) |>
  # (\(x) {x <- droplevels(x); split(x, f = x$Species_name)})() |>
  group_split(Species_name) |>
  future_map(.f = \(x) peak_week(x, model = m_use, draws = drws))

options(op)

peak_week |>
  ggplot(aes(x = Year, y = peak)) +
  geom_point() +
  #geom_line()+
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci)) +
  #facet_grid(Species_name ~ Plot, scales = "free", space = "free_x") +
  labs(y = "DOY of peak", x = NULL,
       title = "Annual peak in abundance of Muscidae species",
       caption = "95% credible interval") +
  theme_bw()


#Onset
df |>
  group_by(Year) |>
  distinct()

ds_pheno <- data_slice(m_gls, DOY = evenly(DOY, by = 1),
                       Year = evenly(Year, by = 1)) %>%
  mutate(fyear = factor(Year, levels = levs)) %>%
  filter(DOY <= 260) #Rethink when we filter

fv_pheno <- fitted_values(m_gls, data = ds_pheno, scale = "response")


#Onset of activity - 10% abundance
onset_week_fv <- fv_pheno |>
  group_by(Species_name,Year) |>
  reframe(csum = cumsum(fitted/7)/sum(fitted/7)) |>
  select(csum) |>
  bind_cols(arrange(ds_pheno, Species_name, Year, DOY))|>
  group_by(Species_name, Year) |>
  filter(csum <= 0.10) |>
  slice_max(order_by = DOY)

end_week <- fv_pheno |>
  group_by(Species_name,Year) |>
  reframe(csum = cumsum(fitted/7)/sum(fitted/7)) |>
  select(csum) |>
  bind_cols(arrange(ds_pheno, Species_name, Year, DOY))|>
  group_by(Species_name, Year) |>
  filter(csum <= 0.90) |>
  slice_max(order_by = DOY)

#Attempt with fitted_samples

onset_week <- fs |>
  group_by(Species_name, Year, draw) |>
  slice_max(order_by = fitted) |>
  group_by(Species_name, Year, DOY) |>
  reframe(csum = cumsum(fitted/7)/sum(fitted/7)) |>
  select(Species_name, Year, DOY, csum) |>
  #bind_cols(arrange(ds_pheno, Species_name, Year, DOY))|>
  group_by(Species_name, Year) |>
  filter(csum <= 0.10) |>
  slice_max(order_by = DOY)


#Figures

onset_week |>
  ggplot(aes(x = Year, y = DOY)) +
  geom_point()+
  #geom_line()+
  #geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci)) +
  #facet_grid(Species_name ~ Plot, scales = "free", space = "free_x") +
  labs(y = "DOY of peak", x = NULL,
       title = "Annual peak in abundance of Muscidae species",
       caption = "95% credible interval") +
  theme_bw()

end_week |>
  ggplot(aes(x = Year, y = DOY)) +
  geom_point()+
  #geom_line()+
  #geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci)) +
  #facet_grid(Species_name ~ Plot, scales = "free", space = "free_x") +
  labs(y = "DOY of peak", x = NULL,
       title = "Annual peak in abundance of Muscidae species",
       caption = "95% credible interval") +
  theme_bw()



#Gavin additions from 27-10-2023

ds_fv <- data_slice(m_gls2, DOY = evenly(DOY, by = 1), 
                       Year = evenly(Year, by = 1), 
                       Plot = evenly(Plot),
                       Species_name = evenly(Species_name)) %>%
  mutate(fyear = factor(Year, levels = levs))

fv_check <- fitted_values(m_gls2, data = ds_fv, scale = "response")

fv_check |>
  filter(Species_name == "Dry.seg") |>
  ggplot(aes(x = DOY, y = .fitted, group = Year, colour = Year)) +
  geom_ribbon(aes(x = DOY, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, 
              inherit.aes = FALSE)+
  geom_line() +
  geom_point(data = df |> filter(Species_name == "Dry.seg"), 
             aes( x = DOY, y = Abundance)) +
  facet_grid(Year ~ Plot)

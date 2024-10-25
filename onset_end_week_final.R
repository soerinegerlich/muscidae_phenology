# GAM full model with species and plot
# This model is close to the final model and takes a long time to run

# Note that this use the very latest version of gratia, that I'll be pushing
# shortly.

# packages - you might need to install some of these
pkgs <- c("tibble", "readr", "dplyr", "here", "tidyr", "ggplot2", "mgcv",
          "gratia", "purrr", "furrr")
vapply(pkgs, library, FUN.VALUE = logical(1L), character.only = TRUE,
       logical.return = TRUE)

df <- read.csv("Data/Species_data_Clean/All_species_clean.csv",sep=",",
               stringsAsFactors = FALSE, header = TRUE)

# Prep data
df <- df |>
  mutate(Plot = as.factor(Plot),
    Species_name = as.factor(Species_name),
    fyear = as.factor(df$Year),
    exposure = rep(7, length(Plot))) # response is  weekly sums

# fit model - slow!
# m_gls2 <- bam(Abundance ~ s(DOY) + s(Year) + ti(DOY, Year) +
#         s(Plot, bs = "re") + s(Species_name, bs = "re") +
#         s(DOY, fyear, bs = "fs", xt = "cr", k = 15) +
#         ti(DOY, Year, Plot, Species_name,
#             bs = c("cr", "cr", "re", "re"), k = c(10, 18, 3, 13)) +
#         offset(log(exposure)),
#     data = df, family = nb(), method = "fREML", discrete = TRUE,
#     nthreads = 4, control = gam.control(trace = TRUE),
#     drop.unused.levels = FALSE)

# Instead, try reading in this saved model
m_gls2 <- read_rds("m_gls2.rds")

####calculating phenologial events####
levs <- with(df, levels(fyear))

# I wa scopying the model into m_use as I was trying different models
m_use <- m_gls2

ds_peak <- data_slice(m_use, DOY = evenly(DOY, by = 1),
  Year = evenly(Year, by = 1),
  Plot = evenly(Plot),
  Species_name = evenly(Species_name)) %>%
  mutate(fyear = factor(Year, levels = levs)) %>%
  filter(DOY <= 260)

# Draw fitted values from the posterior distribution

# generate draws from model of choice - this is where we generate the needed
# set of random draws from the posterior distribution of the model parameters
# this way we use the same set of draws for each species when we run in
# parallel / serial later
#
# You might want to set n to be lower say 100 just to try
drws <- generate_draws(m_use, n = 10000, seed = 42, n_cores = 3)

# wrapper function that does all the steps needed to run the peak week calc
# on any subset of the data you want to pass it.
onset_day <- function(data, model, draws) {
  fit <- fitted_samples(model = model, data = data, scale = "response",
    method = "user", draws = draws)
  data <- data |>
    mutate(.row = row_number())
  fit <- fit |>
    left_join(data, by = join_by(.row))
  fit |>
    group_by(Species_name, Plot, Year, .draw) |>
    arrange(DOY, .by_group = TRUE) |>
    mutate(cum_sum = cumsum(.fitted),
      cum_sum = cum_sum / sum(.fitted)) |>
    dplyr::filter(cum_sum >= 0.1) |>
    slice_min(order_by = cum_sum) |>
    ungroup() |>
    group_by(Species_name, Plot, Year) |>
    summarise(.onset = quantile(DOY, prob = 0.5, na.rm = TRUE),
      .lower_ci = quantile(DOY, prob = 0.025, na.rm = TRUE),
      .upper_ci = quantile(DOY, prob = 0.975, na.rm = TRUE),
      .se = sd(DOY) / sqrt(n()),
      .groups = "drop")
}

# wrapper function that does all the steps needed to run the peak week calc
# on any subset of the data you want to pass it.
end_day <- function(data, model, draws) {
  fit <- fitted_samples(model = model, data = data, scale = "response",
    method = "user", draws = draws)
  data <- data |>
    mutate(.row = row_number())
  fit <- fit |>
    left_join(data, by = join_by(.row))
  fit |>
    group_by(Species_name, Plot, Year, .draw) |>
    arrange(DOY, .by_group = TRUE) |>
    mutate(cum_sum = cumsum(.fitted),
        cum_sum = cum_sum / sum(.fitted)) |>
    dplyr::filter(cum_sum >= 0.9) |>
    slice_min(order_by = cum_sum) |>
    ungroup() |>
    group_by(Species_name, Plot, Year) |>
    summarise(.end = quantile(DOY, prob = 0.5, na.rm = TRUE),
      .lower_ci = quantile(DOY, prob = 0.025, na.rm = TRUE),
      .upper_ci = quantile(DOY, prob = 0.975, na.rm = TRUE),
      .se = sd(DOY) / sqrt(n()),
      .groups = "drop")
}

# Option 1:
# Slowest as it processes each species in turn, one after another, but you
# can run this overnight for example

# this does the work
#Onset
onset <- ds_peak |> # pass our data slice and
    group_split(Species_name) |> # split it by species and then
    purrr::map(.f = \(x) onset_day(x, model = m_use, draws = drws))

#End
end <- ds_peak |> # pass our data slice and
  group_split(Species_name) |> # split it by species and then
  purrr::map(.f = \(x) end_day(x, model = m_use, draws = drws))

# pheno now contains a list, one per species with the peak week per species,
# each element of the list contains the peak DOY info for all Plots and years

# Option 2:
# Using parallel processing - idea is you split this over a few CPUs, processing
# a single species at a time
# This requires the furrr package, which uses the future package for the
# parallel processing

# set up to run using 3 cpu cores
plan(multisession, workers = 3)

# we need to pass some large objects to each of the forks (workers), so we need
# to increase the limit. This works for me but if you get warning about a
# needing to increase this, the 1300 is equivalent of 1.3 Gb so the 1300 is the
# value to increase
op <- options(future.globals.maxSize = 1300 * 1024^2)

# this does the work - you can ignore the warnings
onset <- ds_peak |> # pass our data slice and
    group_split(Species_name) |> # split it by species and then
    future_map(.f = \(x) onset_day(x, model = m_use, draws = drws))
# apply our function to each species, running in parallel with future_map

# The reason for the warnings is that fitted_samples still calls set.seed()
# even when method = "user". I could fix that but I also want to look into the
# options with future more generally for the other methods so right now, just
# ignore the warnings, they are harmless

# reset the options to defaults
options(op)

# pheno now contains a list, one per species with the peak week per species,
# each element of the list contains the peak DOY info for all Plots and years

# Back to common code

# bind all the list elements together
pheno_onset_all <- onset |>
    bind_rows()

pheno_end_all <- end |>
  bind_rows()

# visualise
pheno_onset_all |>
    ggplot(aes(x = Year, y = .onset, colour = Plot)) +
    geom_pointrange(aes(ymin = .lower_ci, ymax = .upper_ci),
        position = position_dodge(width = 0.4), size = 0.3) +
    facet_wrap(~ Species_name)

pheno_end_all |>
  ggplot(aes(x = Year, y = .end, colour = Plot)) +
  geom_pointrange(aes(ymin = .lower_ci, ymax = .upper_ci),
                  position = position_dodge(width = 0.4), size = 0.3) +
  facet_wrap(~ Species_name)

# add something for the weights
pheno_onset_all <- pheno_onset_all |>
    mutate(.inv_std_err = 1/.se)

# Investigating the effect of phylogenetics on phenological responses 
# to snowmelt timing and temperature

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

df_peak$temp_scaled <- as.numeric(scale(df_peak$Temperature, center = F))
df_peak$snow_scaled <- as.numeric(scale(df_peak$Snowmelt, center = F))
df_peak$trend_scaled <- as.numeric(scale(df_peak$trend, center = F))
df_peak$inv_std_error_scaled <- as.numeric(scale(df_peak$inv_std_err, center = F))

# Assign body size to species-specific rows
species_body_size <- c("D.groenlandica" = 5.7, "D.segnis" = 6, "P.bidentata" = 5.6, "S.almqvistii" = 6.05, 
                       "S.dorsata" = 5.10, "S.malaisei" = 3.85, "S.novaesibiriae" = 4.15, "S.sanctipauli" = 4.25, 
                       "S.tundrae" = 5.25, "S.zaitzevi" = 5.05)

df_peak$body_size <- species_body_size[df_peak$species]

species_flower_visitor <- c("D.groenlandica" = "Infrequent", "D.segnis" = "Frequent", "P.bidentata" = "Frequent", "S.almqvistii" = "Frequent", 
                            "S.dorsata" = "Frequent", "S.malaisei" = "Infrequent", "S.novaesibiriae" = "Infrequent", "S.sanctipauli" = "Frequent", 
                            "S.tundrae" = "Infrequent", "S.zaitzevi" = "Infrequent")

df_peak$flower_visiting <- species_flower_visitor[df_peak$species]


df_mean_peak <- read.csv("Data/phenology_estimates/df_mean_peak.csv",sep=",",
                         stringsAsFactors = FALSE, header = TRUE)

df_peak$average_doy <- df_mean_peak$average_doy[match(paste0(df_peak$species),
                                                        paste0(df_mean_peak$species))]


###############################################################################

# Read the CSV file
sequences_df <- read.csv("Data/dna_sequence_muscidae.csv", sep = ";", stringsAsFactors = FALSE)

# Open a text file to write the FASTA format
# Function to write FASTA file
write_fasta <- function(df, output_file) {
  if (!nrow(df)) {
    stop("Data frame is empty! Check the input file.")
  }
  
  # Open a text file to write the FASTA format
  fileConn <- file(output_file, "w")
  
  for (i in 1:nrow(df)) {
    # Write the header and sequence
    writeLines(paste0(">", as.character(df$id[i])), fileConn)  # FASTA header
    writeLines(as.character(df$sequence[i]), fileConn)  # DNA sequence
  }
  
  close(fileConn)
  message("FASTA file written successfully: ", output_file)
}

# Run the function to create a FASTA file
output_file <- "sequences.fasta"
write_fasta(sequences_df, output_file)

# Check if the file exists
if (file.exists(output_file)) {
  message("Output file exists: ", output_file)
} else {
  stop("Output file was not created!")
}

fasta_content <- readLines("sequences.fasta")
print(fasta_content)

library(DECIPHER)
library(Biostrings)

fasta_content <- readDNAStringSet("sequences.fasta")


# Align the sequences
aligned_sequences <- AlignSeqs(fasta_content)

# Check the alignment
print(aligned_sequences)

library(ape)

# Convert to DNAbin
aligned_sequences_dnabin <- as.DNAbin(aligned_sequences)

# Create a distance matrix (using Jukes-Cantor model or similar)
dist_matrix <- dist.dna(aligned_sequences_dnabin, model = "JC")

# Generate the UPGMA tree
upgma_tree <- hclust(dist_matrix, method = "average")  # UPGMA uses the "average" method

# Convert to a phylogenetic tree
phylo_tree <- as.phylo(upgma_tree)

phylo_tree$tip.label <- gsub("^>MW ", "", phylo_tree$tip.label) # Remove ">MW "


# Construct the neighbor-joining tree
nj_tree <- nj(dist_matrix)

nj_tree <- chronos(nj_tree)  # Make tree ultrametric

nj_tree$tip.label <- gsub("^>MW ", "", nj_tree$tip.label) # Remove ">MW "

# Create a phylogenetic variance-covariance matrix
phylo_cor <- corBrownian(phy = nj_tree)

################################################################################

library(nlme)

name_corrections <- c(
  "D.groenlandica" = "Drymeia groenlandica",
  "D.segnis" = "Drymeia segnis",
  "P.bidentata" = "Phaonia bidentata",
  "S.almqvistii" = "Spilogona almqvistii",
  "S.deflorata" = "Spilogona deflorata",
  "S.dorsata" = "Spilogona dorsata",
  "S.malaisei" = "Spilogona malaisei",
  "S.novaesibiriae" = "Spilogona novaesibiriae",
  "S.sanctipauli" = "Spilogona sanctipauli",
  "S.tundrae" = "Spilogona tundrae",
  "S.zaitzevi" = "Spilogona zaitzevi"
)

library(dplyr)  # For recoding

df_peak <- df_peak %>%
  mutate(species = recode(species, !!!name_corrections))

df_peak <- na.omit(df_peak)

species_to_keep <- df_peak$species  # Only keep species in dataset
phylo_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, species_to_keep))
nj_tree <- drop.tip(nj_tree, setdiff(nj_tree$tip.label, species_to_keep))


setdiff(df_peak$species, phylo_tree$tip.label) # Check for species in df_peak missing from tree
setdiff(phylo_tree$tip.label, df_peak$species) # Check for species in tree missing from df_peak
setdiff(nj_tree$tip.label, df_peak$species) # Check for species in tree missing from df_peak


#df_peak$species <- factor(df_peak$species, levels = phylo_tree$tip.label)

################################################################################

# Decided not to use gls model because there are problem with multiple observations
# of the same species. 

# Define the phylogenetic correlation structure
phylo_cor <- corBrownian(phy = phylo_tree, form = ~ species)


pgls_model <- nlme::gls(trend ~ snow_scaled + temp_scaled, 
                        data = df_peak, 
                        correlation = phylo_cor, 
                        weights = varFixed(~ inv_std_error_scaled / mean(inv_std_error_scaled)), 
                        method = "REML")

summary(pgls_model)

library(caper)

# Ensure no duplicates in row names before proceeding
#df_peak <- df_peak[!duplicated(df_peak$species), ]

# Now, convert your data to a 'comparative.data' object
comp_data <- comparative.data(phy = nj_tree, data = df_peak, 
                              names.col = "species", na.omit = TRUE)


# Perform PGLS with the phylogenetic correlation structure
pgls_model_caper <- pgls(trend ~ snow_scaled + temp_scaled, 
                         data = comp_data, 
                         lambda = "ML")  # "ML" estimates the optimal lambda

summary(pgls_model_caper)

###############################################################################

# Using MCMCglmm to model the effects of phylogenetics and species traits on 
# phenological responses to snowmelt timing and temperature

#This guide is great: 
# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://devillemereuil.legtux.org/wp-content/uploads/2021/09/tuto_en.pdf

library(MCMCglmm)
library(lme4)

# Fit a non-phylogenetic mixed model using lme4
lm_model <- lmer(trend ~ snow_scaled + temp_scaled + (1 | species), data = df_peak)

summary(lm_model)


#DIC_mcmc <- mcmc_model$DIC
DIC_mcmc_expand <- mcmc_model_expand$DIC
DIC_lm <- AIC(lm_model)  # AIC for non-phylogenetic model

#print(DIC_mcmc)
print(DIC_mcmc_expand)
print(DIC_lm)

#prior <- list(
#  R = list(V = 1, nu = 0.002),
#  G = list(G1 = list(V = 1, nu = 0.002),
#           G2 = list(V = 1, nu = 0.002))
#)
#
#priorA <- list(
#  R = list(V = 1, fix = 1), 
#  G = list(G1 = list(V = 1, nu = 50), 
#           G2 = list(V = 1, nu = 0.002)))
#
priorB <- list(
  R = list(V = 1, nu = 0.002), 
  G = list(G1 = list(V = 1, nu = 50), 
           G2 = list(V = 1, nu = 50)))

# MCMC expanded
#mcmc_model_expand <- MCMCglmm(trend ~ snow_scaled + temp_scaled + body_size + flower_visiting + average_doy, 
#                       random = ~ species, 
#                       family = "gaussian", 
#                      ginverse = list(species = inverseA(nj_tree, nodes = "ALL")$Ainv),
#                      data = df_peak, 
#                      prior = prior, 
#                      nitt = 500000, burnin = 100000, thin = 100)

df_peak$species_plot <- interaction(df_peak$species, df_peak$plot)

mcmc_model_expand <- MCMCglmm(trend ~ snow_scaled + temp_scaled + year + body_size + flower_visiting + average_doy, 
                              random = ~ species + species_plot,  # Use species_plot as random effect
                              family = "gaussian", 
                              ginverse = list(species = inverseA(nj_tree, nodes = "ALL")$Ainv), 
                              data = df_peak, 
                              prior = priorB, 
                              nitt = 800000, burnin = 200000, thin = 100)

summary(mcmc_model_expand)

# Trace plots

# Each of the three couples of graphs shows us the trace (left), 
# i.e. the evolution of the sampled values along the iterations. 
# It allows us to check the convergence (we should not see any trend in the trace) 
# and that autocorrelation is weak (the values are widely spread, not following 
# a traceable path). On the right of these graphs, we have an estimation of the 
# posterior distribution for each component (Intercept, animal and units).

plot(mcmc_model_expand$Sol)  # Fixed effects trace plot
plot(mcmc_model_expand$VCV)  # Random effects & residual variance trace plot

# Effect size

# the aim of MCMC is to estimate the posterior distribution of the parameter of 
# interest by sampling from it. To do so, we need the largest number of independent 
# values as possible, which means a large effective sample size. 
# In practice, an effective size above 1,000 is recommended.

effectiveSize(mcmc_model_expand[["Sol"]])
effectiveSize(mcmc_model_expand[["VCV"]])

# Can the differences in effective sample size be explained by differences in 
# auto-correlation (as they should)?

autocorr(mcmc_model_expand$Sol)  # Fixed effects
autocorr(mcmc_model_expand$VCV)  # Random effects & residual variance

# Here Lag 10 stands for ‘autocorrelation every 10 iteration values’. 
# Since our thin parameter was 10, this refers actually to the correlation of every 
# sampled value with the following one. We can see that there is little 
# autocorrelation on the mean (Intercept). On the contrary, the autocorrelation 
# on variance componentsbecomes negligible rather with a lag of 50. 
# The best way to decide whether this level of auto-correlation is problematic 
# or not is rather to look at the effective sample size, whether it is acceptable 
# or not (as draws samples from longer chains with higher autocorrelation can 
# have a higher effective sample size than from a shorter chains with lesser 
# autocorrelation).

# Diagnostic test of convergence - Heidelberg stationarity test

# The“HalfwidthMean”partof heidel.diagisnotrelatedtoconvergence,buttrytoassess 
# whether the chain was run long enough to get a given precision level. 
# It can be quite sensitive to departure from normality and is not as informative 
# as the computation of effective sample size in the end.

heidel.diag(mcmc_model_expand[["VCV"]]) # The test needs to pass for all fixed
# and random variables. If they fail, try to change priors or no. of iterations 
# and burnins.

#lambda_est <- mcmc_model_expand$VCV[,"species"] / 
#  (mcmc_model_expand$VCV[,"species_plot"] + mcmc_model_expand$VCV[,"units"] + mcmc_model_expand$VCV[,"species"])

#mean(lambda_est)  # Average lambda

# Posterior distribution of heritability

# To answer the question: Is the phylogenetic signal in phenological responses strong?

herit <- mcmc_model_expand[["VCV"]][ , "species"] / 
  (mcmc_model_expand$VCV[,"species_plot"] + mcmc_model_expand$VCV[,"units"] + 
     mcmc_model_expand$VCV[,"species"])

effectiveSize(herit)

mean(herit)

HPDinterval(herit)

plot(herit)


lm.test <- lm(trend ~ snow_scaled + temp_scaled, data = df_peak, weights = inv_std_error_scaled / mean(inv_std_error_scaled))

summary(lm.test)

# The phylogenetic signal is poor: h = 0.05, CI = 0.03, 0.07


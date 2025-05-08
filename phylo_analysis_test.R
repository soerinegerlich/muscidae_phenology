# Investigating the effect of phylogenetics on phenological responses 
# to snowmelt timing and temperature

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

df_onset$temp_scaled <- as.numeric(scale(df_onset$Temperature, center = F))
df_onset$snow_scaled <- as.numeric(scale(df_onset$Snowmelt, center = F))
df_onset$trend_scaled <- as.numeric(scale(df_onset$trend, center = F))
df_onset$inv_std_error_scaled <- as.numeric(scale(df_onset$inv_std_err, center = F))

# Assign body size to species-specific rows
species_body_size <- c("D.groenlandica" = 5.7, "D.segnis" = 6, "P.bidentata" = 5.6, "S.almqvistii" = 6.05, 
                       "S.dorsata" = 5.10, "S.malaisei" = 3.85, "S.novaesibiriae" = 4.15, "S.sanctipauli" = 4.25, 
                       "S.tundrae" = 5.25, "S.zaitzevi" = 5.05)

df_onset$body_size <- species_body_size[df_onset$species]


species_flower_visitor <- c("D.groenlandica" = "Infrequent", "D.segnis" = "Frequent", "P.bidentata" = "Frequent", "S.almqvistii" = "Frequent", 
                       "S.dorsata" = "Frequent", "S.malaisei" = "Infrequent", "S.novaesibiriae" = "Infrequent", "S.sanctipauli" = "Frequent", 
                       "S.tundrae" = "Infrequent", "S.zaitzevi" = "Infrequent")

df_onset$flower_visiting <- species_flower_visitor[df_onset$species]


df_mean_onset <- read.csv("Data/phenology_estimates/df_mean_onset.csv",sep=",",
                          stringsAsFactors = FALSE, header = TRUE)


df_onset$average_doy <- df_mean_onset$average_doy[match(paste0(df_onset$species),
                                                                paste0(df_mean_onset$species))]


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
df_onset <- df_onset %>%
  mutate(species = recode(species, !!!name_corrections))

df_onset <- na.omit(df_onset)

species_to_keep <- df_onset$species  # Only keep species in dataset
phylo_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, species_to_keep))
nj_tree <- drop.tip(nj_tree, setdiff(nj_tree$tip.label, species_to_keep))



setdiff(df_onset$species, phylo_tree$tip.label) # Check for species in df_onset missing from tree
setdiff(phylo_tree$tip.label, df_onset$species) # Check for species in tree missing from df_onset
setdiff(nj_tree$tip.label, df_onset$species) # Check for species in tree missing from df_onset


#df_onset$species <- factor(df_onset$species, levels = phylo_tree$tip.label)

# Define the phylogenetic correlation structure
phylo_cor <- corBrownian(phy = phylo_tree, form = ~ species)


pgls_model <- nlme::gls(trend ~ snow_scaled + temp_scaled, 
                  data = df_onset, 
                  correlation = phylo_cor, 
                  weights = varFixed(~ inv_std_error_scaled / mean(inv_std_error_scaled)), 
                  method = "REML")

summary(pgls_model)

library(caper)

# Ensure no duplicates in row names before proceeding
#df_onset <- df_onset[!duplicated(df_onset$species), ]

# Now, convert your data to a 'comparative.data' object
comp_data <- comparative.data(phy = nj_tree, data = df_onset, names.col = "species", na.omit = TRUE)


# Perform PGLS with the phylogenetic correlation structure
pgls_model_caper <- pgls(trend ~ snow_scaled + temp_scaled, 
                         data = comp_data, 
                         lambda = "ML")  # "ML" estimates the optimal lambda

summary(pgls_model_caper)


library(MCMCglmm)

pmm_model <- MCMCglmm(trend ~ snow_scaled + temp_scaled, 
                      random = ~species,  # Allows for repeated measures
                      family = "gaussian",
                      ginverse = list(species = inverseA(nj_tree, nodes = "TIPS")$Ainv), 
                      data = df_onset)

summary(pmm_model)

plot(pmm_model$Sol)  # Fixed effects trace plot
plot(pmm_model$VCV)  # Random effects & residual variance trace plot

autocorr(pmm_model$Sol)  # Fixed effects
autocorr(pmm_model$VCV)  # Random effects & residual variance

print(pmm_model$DIC)

lambda_est <- pmm_model$VCV[,"species"] / 
  (pmm_model$VCV[,"species"] + pmm_model$VCV[,"units"])

mean(lambda_est)  # Average lambda


library(MCMCglmm)

# Define prior
prior <- list(G = list(G1 = list(V = 1, nu = 0.002)), 
              R = list(V = 1, nu = 0.002))

# Fit the phylogenetic model
mcmc_model <- MCMCglmm(trend ~ snow_scaled + temp_scaled, 
                       random = ~ species, 
                       family = "gaussian", 
                       ginverse = list(species = inverseA(nj_tree, nodes = "ALL")$Ainv),
                       data = df_onset, 
                       prior = prior, 
                       nitt = 200000, burnin = 50000, thin = 100)

summary(mcmc_model)

plot(mcmc_model$Sol)  # Fixed effects trace plot
plot(mcmc_model$VCV)  # Random effects & residual variance trace plot

effectiveSize(mcmc_model[["Sol"]])
effectiveSize(mcmc_model[["VCV"]])

pp_check(mcmc_model)

autocorr(mcmc_model$Sol)  # Fixed effects
autocorr(mcmc_model$VCV)  # Random effects & residual variance

heidel.diag(mcmc_model[["VCV"]])

lambda_est <- mcmc_model$VCV[,"species"] / 
  (mcmc_model$VCV[,"species"] + mcmc_model$VCV[,"units"])

mean(lambda_est)  # Average lambda

# Extract variance components from the MCMCglmm model
vc <-mcmc_model$VCV  # Extract posterior samples of variance components

vc <- posterior.mode(mcmc_model$VCV)
Vphy <- vc["species"]  # Phylogenetic variance
Vres <- vc["units"]  # Residual variance
heritability <- Vphy / (Vphy + Vres)
print(heritability)

herit <- mcmc_model[["VCV"]][ , "species"] / (mcmc_model[["VCV"]][ , "species"] + mcmc_model[["VCV"]][ , "units"])

effectiveSize(herit)

mean(herit)

HPDinterval(herit)

plot(herit)



library(lme4)

# Fit a non-phylogenetic mixed model using lme4
lm_model <- lmer(trend ~ snow_scaled + temp_scaled + (1 | species), data = df_onset)

summary(lm_model)


DIC_mcmc <- mcmc_model$DIC
DIC_mcmc_expand <- mcmc_model_expand$DIC
DIC_lm <- AIC(lm_model)  # AIC for non-phylogenetic model

print(DIC_mcmc)
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
#                      data = df_onset, 
#                      prior = prior, 
#                      nitt = 500000, burnin = 100000, thin = 100)

df_onset$species_plot <- interaction(df_onset$species, df_onset$plot)

mcmc_model_expand <- MCMCglmm(trend ~ snow_scaled + temp_scaled + year + body_size + flower_visiting + average_doy, 
                              random = ~ species + species_plot,  # Use species_plot as random effect
                              family = "gaussian", 
                              ginverse = list(species = inverseA(nj_tree, nodes = "ALL")$Ainv), 
                              data = df_onset, 
                              prior = priorB, 
                              nitt = 800000, burnin = 200000, thin = 100)

summary(mcmc_model_expand)


plot(mcmc_model_expand$Sol)  # Fixed effects trace plot
plot(mcmc_model_expand$VCV)  # Random effects & residual variance trace plot

effectiveSize(mcmc_model_expand[["Sol"]])
effectiveSize(mcmc_model_expand[["VCV"]])

autocorr(mcmc_model_expand$Sol)  # Fixed effects
autocorr(mcmc_model_expand$VCV)  # Random effects & residual variance

heidel.diag(mcmc_model_expand[["VCV"]])

#lambda_est <- mcmc_model_expand$VCV[,"species"] / 
#  (mcmc_model_expand$VCV[,"species_plot"] + mcmc_model_expand$VCV[,"units"] + mcmc_model_expand$VCV[,"species"])

#mean(lambda_est)  # Average lambda


herit <- mcmc_model_expand[["VCV"]][ , "species"] / (mcmc_model_expand$VCV[,"species_plot"] + mcmc_model_expand$VCV[,"units"] + mcmc_model_expand$VCV[,"species"])

effectiveSize(herit)


mean(herit)

HPDinterval(herit)

plot(herit)


lm.test <- lm(trend ~ snow_scaled + temp_scaled, data = df_onset, weights = inv_std_error_scaled / mean(inv_std_error_scaled))

summary(lm.test)




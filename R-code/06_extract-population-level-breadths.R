#### PROJECT: Performance curves of Mimulus species
#### PURPOSE: Calculation of species-level breadth for tpcs
#### AUTHOR/DATE: RCW/2020-07-15


### load packages 
library(ggplot2)
library(dplyr)


#############################################
### read data 
# this dataset contains family-level estimates of thermal performance curve
# parameters across 10 mimulus species
#############################################

tpc <- read.csv("Analysis-output/family-mean-estimates-sxp-af.csv") 


#############################################
### for thermal performance curves, select minimum lower limit and 
### maximum upper limit across families within each of the 10 species
### note: we are considering t.min and t.max that are limited to fall
### within the measurement intervel (15 or above, and 50 or below)
#############################################

# thermal maximum
tpc_t.max <- tpc %>%
  group_by(species) %>%
  dplyr::select(species, range, pair, t.max2) %>%
  slice(which.max(t.max2))
# thermal minimum
tpc_t.min <- tpc %>%
  group_by(species) %>%
  dplyr::select(species, t.min2) %>%
  slice(which.min(t.min2))



#############################################
### bind columns together in one dataset, then clean up
#############################################

tpc_species <- dplyr::bind_cols(tpc_t.max, tpc_t.min)

tpc_species2 <- tpc_species %>% dplyr::rename(species = species...1)


tpc_species3 <- tpc_species2 %>%
  dplyr::select(species, range, pair, t.min2, t.max2)

# calculate species breadth
tpc_species4 <- tpc_species3 %>%
  mutate(t.breadth = t.max2 - t.min2)

# these are the temperature limits/breadths for each species
write.csv(tpc_species4, "Analysis-output/Species/species-limits.csv")


# now calculate how widespread and restricted species differ in thermal breadth

tpc_species4$pair <- as.character(tpc_species4$pair)
pairs <- unique(tpc_species4$pair)
species_diffs <- data.frame(pair = c(pairs),
                            t.breadth = c(-diff(tpc_species4$t.breadth[which(tpc_species4$pair==pairs[1])]),
                                          -diff(tpc_species4$t.breadth[which(tpc_species4$pair==pairs[2])]),
                                          -diff(tpc_species4$t.breadth[which(tpc_species4$pair==pairs[3])]),
                                          -diff(tpc_species4$t.breadth[which(tpc_species4$pair==pairs[4])]),
                                          -diff(tpc_species4$t.breadth[which(tpc_species4$pair==pairs[5])])))
# positive values mean that the widespread species has a higher breadth, 
# and negative values mean the restricted species has a higher breadth
#write.csv(species_diffs, "Analysis-output/Species/species-breadth-differences-sxp-af.csv")

# create description of values
species_diffs_metadata <- "positive values mean that the widespread species has a higher breadth, and negative values mean the restricted species has a higher breadth"
#write.csv(species_diffs_metadata, "Analysis-output/Species/species-breadth-differences_metadata-sxp-af.csv")






#### PROJECT: Performance curves of Mimulus species
#### PURPOSE: Take the parameter estimates from the models and join 
#### them together in one dataset
#### AUTHOR/DATE: RCW/2020-07-07

## So far, we have used hierarchical models to estimate thermal performance curves
## and moisture performance curves of each family, using separate models for each
## species. Performance curves have been estimated based on stem RGR for some
## species (bic, fil, gut, lac) and leaf RGR for others (car, eas, flo, nor, 
## par, ver). In this R script, we will import the mean parameter estimates 
## from each model and merge them into one dataset to use in mixed effects models.

# load appropriate packages
library(ggplot2)
library(dplyr)



#############################################
### read in group means for each MPC parameter
#############################################
# mean_mpc_car <- read.csv("Analysis-output/MPC/Mean-dfs/mean_df_sxp_car.csv")[,-1] # cardinalis (widespread)
# mean_mpc_par <- read.csv("Analysis-output/MPC/Mean-dfs/mean_df_sxp_par.csv")[,-1] # parishii (restricted)
# mean_mpc_ver <- read.csv("Analysis-output/MPC/Mean-dfs/mean_df_sxp_ver.csv")[,-1] # verbenaceus (widespread)
# mean_mpc_eas <- read.csv("Analysis-output/MPC/Mean-dfs/mean_df_sxp_eas.csv")[,-1] # eastwoodiae (restricted)
# mean_mpc_flo <- read.csv("Analysis-output/MPC/Mean-dfs/mean_df_sxp_flo.csv")[,-1] # floribundus (widespread)
# mean_mpc_nor <- read.csv("Analysis-output/MPC/Mean-dfs/mean_df_sxp_nor.csv")[,-1] # norrisii (restricted)
# mean_mpc_bic <- read.csv("Analysis-output/MPC/Mean-dfs/mean_df_sxp_bic.csv")[,-1] # bicolor (widespread)
# mean_mpc_fil <- read.csv("Analysis-output/MPC/Mean-dfs/mean_df_sxp_fil.csv")[,-1] # filicaulis (restricted)
# mean_mpc_gut <- read.csv("Analysis-output/MPC/Mean-dfs/mean_df_sxp_gut.csv")[,-1] # guttatus (widespread)
# mean_mpc_lac <- read.csv("Analysis-output/MPC/Mean-dfs/mean_df_sxp_lac.csv")[,-1] # laciniatus (restricted)

#############################################
### read in group means for each TPC parameter
#############################################
mean_tpc_car <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_car.csv")[,-1] # cardinalis (widespread)
mean_tpc_par <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_par.csv")[,-1] # parishii (restricted)
mean_tpc_ver <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ver.csv")[,-1] # verbenaceus (widespread)
mean_tpc_eas <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_eas.csv")[,-1] # eastwoodiae (restricted)
mean_tpc_flo <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_flo.csv")[,-1] # floribundus (widespread)
mean_tpc_nor <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_nor.csv")[,-1] # norrisii (restricted)
mean_tpc_bic <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_bic.csv")[,-1] # bicolor (widespread)
mean_tpc_fil <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_fil.csv")[,-1] # filicaulis (restricted)
mean_tpc_gut <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_gut.csv")[,-1] # guttatus (widespread)
mean_tpc_lac <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_lac.csv")[,-1] # laciniatus (restricted)


#############################################
### read in reference dataframes for original families
#############################################
refdat <- read.csv("Processed-data/temp.refdat_all-fams.csv")


#############################################
### join mpc/tpc datasets together 
#############################################

# first, delineate families and species in each of the MPC datasets
# mean_mpc_car$family <- mean_mpc_car$species # the 'species' from the analysis actually represent different families
# mean_mpc_car$species <- rep("car", dim(mean_mpc_car)[1]) # add species character
# 
# mean_mpc_par$family <- mean_mpc_par$species 
# mean_mpc_par$species <- rep("par", dim(mean_mpc_par)[1]) 
# 
# mean_mpc_ver$family <- mean_mpc_ver$species 
# mean_mpc_ver$species <- rep("ver", dim(mean_mpc_ver)[1])
# 
# mean_mpc_eas$family <- mean_mpc_eas$species 
# mean_mpc_eas$species <- rep("eas", dim(mean_mpc_eas)[1]) 
# 
# mean_mpc_flo$family <- mean_mpc_flo$species 
# mean_mpc_flo$species <- rep("flo", dim(mean_mpc_flo)[1])
# 
# mean_mpc_nor$family <- mean_mpc_nor$species 
# mean_mpc_nor$species <- rep("nor", dim(mean_mpc_nor)[1]) 
# 
# mean_mpc_bic$family <- mean_mpc_bic$species 
# mean_mpc_bic$species <- rep("bic", dim(mean_mpc_bic)[1])
# 
# mean_mpc_fil$family <- mean_mpc_fil$species 
# mean_mpc_fil$species <- rep("fil", dim(mean_mpc_fil)[1]) 
# 
# mean_mpc_gut$family <- mean_mpc_gut$species 
# mean_mpc_gut$species <- rep("gut", dim(mean_mpc_gut)[1])
# 
# mean_mpc_lac$family <- mean_mpc_lac$species 
# mean_mpc_lac$species <- rep("lac", dim(mean_mpc_lac)[1]) 


# repeat for each of the TPC datasets
mean_tpc_car$family <- mean_tpc_car$species # the 'species' from the analysis actually represent different families
mean_tpc_car$species <- rep("car", dim(mean_tpc_car)[1]) # add species character

mean_tpc_par$family <- mean_tpc_par$species 
mean_tpc_par$species <- rep("par", dim(mean_tpc_par)[1]) 

mean_tpc_ver$family <- mean_tpc_ver$species 
mean_tpc_ver$species <- rep("ver", dim(mean_tpc_ver)[1])

mean_tpc_eas$family <- mean_tpc_eas$species 
mean_tpc_eas$species <- rep("eas", dim(mean_tpc_eas)[1]) 

mean_tpc_flo$family <- mean_tpc_flo$species 
mean_tpc_flo$species <- rep("flo", dim(mean_tpc_flo)[1])

mean_tpc_nor$family <- mean_tpc_nor$species 
mean_tpc_nor$species <- rep("nor", dim(mean_tpc_nor)[1]) 

mean_tpc_bic$family <- mean_tpc_bic$species 
mean_tpc_bic$species <- rep("bic", dim(mean_tpc_bic)[1])

mean_tpc_fil$family <- mean_tpc_fil$species 
mean_tpc_fil$species <- rep("fil", dim(mean_tpc_fil)[1]) 

mean_tpc_gut$family <- mean_tpc_gut$species 
mean_tpc_gut$species <- rep("gut", dim(mean_tpc_gut)[1])

mean_tpc_lac$family <- mean_tpc_lac$species 
mean_tpc_lac$species <- rep("lac", dim(mean_tpc_lac)[1]) 




# second, bind rows for mpc data and then for tpc data
# mean_mpc <- bind_rows(mean_mpc_car, mean_mpc_par,
#                       mean_mpc_ver, mean_mpc_eas,
#                       mean_mpc_flo, mean_mpc_nor,
#                       mean_mpc_bic, mean_mpc_fil,
#                       mean_mpc_gut, mean_mpc_lac)
mean_tpc <- bind_rows(mean_tpc_car, mean_tpc_par,
                      mean_tpc_ver, mean_tpc_eas,
                      mean_tpc_flo, mean_tpc_nor,
                      mean_tpc_bic, mean_tpc_fil,
                      mean_tpc_gut, mean_tpc_lac)

# third, bind columns for tpc + mpc data
#colnames(mean_mpc) <- paste0("m.", colnames(mean_mpc)) 
colnames(mean_tpc) <- paste0("t.", colnames(mean_tpc))

mean_df <- mean_tpc


### retain important columns
mean_df_1 <- mean_df %>% select(t.species, t.family,
                                t.maximaBT, t.max_RGRBT,  # temperature optimum, pmax
                                t.x_minBT, t.x_maxBT, t.breadthBT, # temperature critical minimum, critical maximum, and tolerance range
                                t.B50, t.B50_low, t.B50_high) # temperature B50, lower limit, upper limit 
# rename columns
mean_df_2 <- mean_df_1 %>% rename(species = t.species, 
                                  family = t.family,
                                  t.opt = t.maximaBT, # temperature optimum
                                  t.pmax = t.max_RGRBT, # temperature performance maximum
                                  t.cmin = t.x_minBT, # temperature critical minimum
                                  t.cmax = t.x_maxBT, # temperature critical maximum
                                  t.cbreadth = t.breadthBT, # temperature critical breadth / tolerance range
                                  t.breadth = t.B50, # temperature breadth
                                  t.min = t.B50_low, # temperature minimum
                                  t.max = t.B50_high) # temperature maximum

# limit t.min and t.max to fall at or above 15 (t.min) and at or below 50 (t.max)
# because measurement interval is 15 to 50
mean_df_2 <- mean_df_2 %>%
  mutate(t.min2=replace(t.min, t.min<15, 15), t.max2=replace(t.max, t.max>50, 50))


#############################################
### create a column for original family numbers
#############################################
mean_df_2$family.orig <- rep(NA, dim(mean_df_2)[1])
for(i in 1:dim(mean_df_2)[1]){
  mean_df_2$family.orig[i] <- refdat$family[which(refdat$species==mean_df_2$species[i] & refdat$familyI==mean_df_2$family[i])]
}
length(unique(mean_df_2$family.orig)) # check



#############################################
### create a column for widespread vs. restricted
#############################################
mean_df_2$range <- rep("restricted", dim(mean_df_2)[1])

mean_df_2$range[which(mean_df_2$species %in% c("car", "ver", "flo", "bic", "gut"))] <- "widespread"

table(mean_df_2$species, mean_df_2$range) # check


#############################################
### create a column for species pair
#############################################
mean_df_2$pair <- rep("car-par", dim(mean_df_2)[1])

mean_df_2$pair[which(mean_df_2$species %in% c("ver", "eas"))] <- "ver-eas"
mean_df_2$pair[which(mean_df_2$species %in% c("flo", "nor"))] <- "flo-nor"
mean_df_2$pair[which(mean_df_2$species %in% c("bic", "fil"))] <- "bic-fil"
mean_df_2$pair[which(mean_df_2$species %in% c("gut", "lac"))] <- "gut-lac"

table(mean_df_2$species, mean_df_2$pair) # check





#############################################
### export dataset
#############################################

write.csv(mean_df_2, "Analysis-output/family-mean-estimates-sxp-af.csv")

#############################################
### export medadata
#############################################

metadat <- data.frame(colname = c(colnames(mean_df_2)),
                      description = c("three-letter code for one of ten Mimulus species: cardinalis, parishii, verbenaceus, eastwoodiae, floribundus, norrisii, bicolor, filicaulis, guttatus, or laciniatus",
                                      "integer representing family in the performance curve models: 1 to number of families in each species",
                                      "thermal optimum, or temperature at which performance is maximized", 
                                      "performance maximum along the temperature gradient",
                                      "critical thermal minimum, or the lower temperature where performance reaches zero", 
                                      "critical thermal maximum, or the upper temperature where performance reaches zero",
                                      "critical breadth for temperature, or tolerance range where performance is above zero across the temperature gradient",
                                      "thermal breadth, or span of temperature across which at least 50% of t.pmax is achieved",
                                      "lower thermal limit, or the lower bound of t.breadth (50% threshold)",
                                      "upper thermal limit, or the upper bound of t.breadth (50% threshold)",
                                      "lower thermal limit, or the lower bound of t.breadth (50% threshold) after being limited to measurement interval",
                                      "upper thermal limit, or the upper bound of t.breadth (50% threshold) after being limited to measurement interval",
                                      "family number assigned in original datasets, which is unique across all families: 1 to 306",
                                      "indicates whether each family is from a widespread or restricted species",
                                      "indicates the species pair that each family falls into"))

write.csv(metadat, "Analysis-output/family-mean-estimates-metadata-sxp-af.csv")







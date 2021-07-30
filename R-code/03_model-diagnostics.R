#### PROJECT: Performance curves of Mimulus species
#### PURPOSE: Explore the fits of performance curve models
#### AUTHOR/DATE: RCW/2020-06-24

library(dplyr)
library(ggplot2)

# we ran tpc models for each species with family as the grouping variable
# we set max_tree_depth to 15, adapt_delta to 0.96, & number of iterations to 6000 
# for models of all species to get proper diagnostic statistics across species
# statistics we care aboue:
# rhat: should be 1 (indicating convergence)
# neff: should be close to the number of posterior samples



##################################################################
# TPCs
##################################################################

#################################
# read in the model summaries
#################################

ms_bic <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_bic.csv")
ms_car <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_car.csv")
ms_eas <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_eas.csv")
ms_fil <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_fil.csv")
ms_flo <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_flo.csv")
ms_gut <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_gut.csv")
ms_lac <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_lac.csv")
ms_nor <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_nor.csv")
ms_par <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_par.csv")
ms_ver <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_ver.csv")


# merge datasets so it's easier to explore
ms_all <- dplyr::bind_rows(ms_bic, ms_car, ms_eas, ms_fil, ms_flo, 
                           ms_gut, ms_lac, ms_nor, ms_par, ms_ver)
ms_all <- ms_all %>% mutate(species = rep(c("bic", "car", "eas", "fil", "flo",
                                            "gut", "nor", "lac", "par", "ver"), 
                                          times=c(dim(ms_bic)[1], dim(ms_car)[1], dim(ms_eas)[1], dim(ms_fil)[1], dim(ms_flo)[1], 
                                                  dim(ms_gut)[1], dim(ms_lac)[1], dim(ms_nor)[1], dim(ms_par)[1], dim(ms_ver)[1])))

#################################
# examine at rhats 
#################################
table(ms_all$species, ms_all$Rhat)
# rhats are 1 for all parameters of each species



#################################
# examine at n_eff 
#################################
ggplot(ms_all, aes(x=n_eff)) +
  geom_histogram() +
  facet_wrap(species~.)
ms_all %>% 
  group_by(species) %>%
  slice(which.min(n_eff))
# all n_eff are > 700,  
# the lowest n_eff are for parameters involved in critical limits and zero-inflation

ms_all %>% 
  group_by(species) %>%
  summarize(mean_neff = mean(n_eff, na.rm = TRUE))







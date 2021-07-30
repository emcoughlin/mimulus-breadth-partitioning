#### PROJECT: Performance curves of Mimulus species
#### PURPOSE: Model temperature performance curves of 
#### families within each species
#### AUTHOR/DATE: RCW/2020-02-10

###### Warnings:
###### 1) This script requires a lot of computational power.
###### 2) The stan_performance output files are very large and thus not located in this git repo. Please contact Seema Sheth (ssheth3@ncsu.edu) if you are interested in obtaining these files.



#load performr
devtools::install_github("silastittes/performr", local = FALSE, ref="zin", force=FALSE) 

#load other libraries 
library(dplyr)
library(devtools)
library(performr)
library(tidyverse)
library(ggridges)
library(ggpubr)


# Other settings
options(mc.cores = parallel::detectCores())
extract <- rstan::extract



############ 
### The data:
############ 
# This is the raw thermal performance data from Sheth & Angert 2014 (Evolution). 
# Stem RGR was used for some species (bic, fil, gut, and lac), 
# and leaf RGR was used for the remaining species (car, eas, flo, nor, par, ver), 
# so we need to choose a different RGR column depending on the species. 
dat <- read.csv("Raw-data/thermal_data.csv")
## Column descriptions:
# species: species of Mimulus
# family: unique ID corresponding to full-sibling seed family
# dayTemp: daytime temperature (in degrees C) of growth chamber
# nightTemp: nighttime temperature (in degrees C) of growth chamber
# stemLen1: stem length (in cm) before going into growth chamber
# leafNum1: leaf # before going into growth chamber
# stemLen2: stem length (in cm) after 7 days in growth chamber
# leafNum2: leaf # after 7 days in growth chamber
# rgrStem: relative growth rate in stem length (in units of cm/cm/day) 
# rgrLeaf: relative growth rate in leaf number (in units of leaf #/leaf #/day) 
                                              
          



############ 
############ SCALING DATA
############ 



############ 
# reference dataframe for family:
# to create a column that groups the families as integers from 
# 1 to the number of families within each species
refdat <- as.data.frame(dat %>% dplyr::group_by(species,family) %>% dplyr::summarise(count=n()))
famcounts <- as.data.frame(refdat %>% dplyr::group_by(species)%>% dplyr::summarise(count=n()))
write.csv(famcounts, "Processed-data/temp.famcounts_all-fams.csv")
refdat$familyI <- c(1:famcounts$count[1],1:famcounts$count[2],1:famcounts$count[3],1:famcounts$count[4],
                    1:famcounts$count[5],1:famcounts$count[6],1:famcounts$count[7],1:famcounts$count[8],
                    1:famcounts$count[9],1:famcounts$count[10])
write.csv(refdat[,c("species","family","familyI")], "Processed-data/temp.refdat_all-fams.csv")

# apply new fam ord to dat and check dat
dat$familyI <- rep(NA, dim(dat)[1])
for(i in 1:dim(dat)[1]){
  dat$familyI[i] <- refdat$familyI[which(refdat$family==dat$family[i])]
}
tab <- dat[,c("species", "family", "familyI")]



############ 
# make a new column of the specific RGR for each species
# Stem RGR was used for some species (bic, fil, gut, and lac), 
# and leaf RGR was used for the remaining species (car, eas, flo, nor, par, ver), 
dat$RGR <- dat$rgrStem
dat$RGR[which(dat$species %in% c("car", "eas", "flo", "nor", "par", "ver"))] <- dat$rgrLeaf[which(dat$species %in% c("car", "eas", "flo", "nor", "par", "ver"))]


############ 
# Get raw overall mean of RGR and Temp from dataset of each species
# We will use this to scale RGR and center Temp data for the bayesian models
meansTot <- dat %>% 
  dplyr::group_by(species) %>%
  dplyr::select(dayTemp, RGR) %>% 
  dplyr::summarize(RGR = mean(RGR, na.rm=TRUE),
            dayTemp = mean(dayTemp, na.rm=TRUE))
### TO STANDARDIZE X-AXIS PARAMETERS ACROSS SPECIES, 
### REPLACE DAYTEMP WITH 33, WHICH IS CLOSE TO 
### THE MEAN TREATMENT ACROSS ALL SPECIES 
meansTot$dayTemp <- rep(33, dim(meansTot)[1])
meansTot <- as.data.frame(meansTot)
write.csv(meansTot, "Processed-data/temp.meansTot_all-fams.csv")


############ 
# Center temp around zero and scale RGR by species means
dat$dayTempc <- rep(NA, dim(dat)[1]) # centered temperature
dat$RGRs <- rep(NA, dim(dat)[1]) # scaled rgr
for(i in 1:dim(dat)[1]){
  dat$dayTempc[i] <- dat$dayTemp[i] - meansTot$dayTemp[which(meansTot$species==dat$species[i])]
  dat$RGRs[i] = dat$RGR[i] / meansTot$RGR[which(meansTot$species==dat$species[i])]
}
write.csv(dat, "Processed-data/temp.mod.dat_all-fams.csv")




# plot raw data for floribundus
ggplot(dat[which(dat$species=="flo"),], aes(x=dayTemp, y=RGR)) +
  geom_jitter(width=0.1, height=0.01, shape=1, size=0.5) +
  facet_wrap(family~.) +
  ylim(0,max(dat$RGR[which(dat$species=="flo")])) + xlim(5,55)

# plot transformed data for floribundus
ggplot(dat[which(dat$species=="flo"),], aes(x=dayTempc, y=RGRs)) +
  geom_jitter(width=0.05, height=0.1, shape=1, size=0.5) +
  facet_wrap(familyI~.)





############ 
############ BAYESIAN TPCS FOR EACH SPECIES, WITH FAMILY AS THE GROUPING FACTOR
############ !CAUTION! THESE MODELS TAKE A LONG TIME TO RUN
############
### SETTING ALL X_MIN AND X_MAX PRIORS TO -18 AND 17, RESPECTIVELY

stan_performance(df = filter(dat, species=="bic"),
                 response = RGRs,
                 treatment = dayTempc,
                 group_ids = familyI,
                 min_pr_mu = -18,
                 max_pr_mu = 17,
                 seed = 1234,
                 max_treedepth = 12,
                 adapt_delta = 0.95,
                 file_id =  "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_bic",
                 iter = 6000)
stan_performance(df = filter(dat, species=="car"),
                 response = RGRs,
                 treatment = dayTempc,
                 group_ids = familyI,
                 min_pr_mu = -18,
                 max_pr_mu = 17,
                 seed = 1234,
                 max_treedepth = 12,
                 adapt_delta = 0.95,
                 file_id =  "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_car",
                 iter = 6000)
stan_performance(df = filter(dat, species=="eas"),
                 response = RGRs,
                 treatment = dayTempc,
                 group_ids = familyI,
                 min_pr_mu = -18,
                 max_pr_mu = 17,
                 seed = 1234,
                 max_treedepth = 12,
                 adapt_delta = 0.95,
                 file_id =  "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_eas",
                 iter = 6000)
stan_performance(df = filter(dat, species=="fil"),
                 response = RGRs,
                 treatment = dayTempc,
                 group_ids = familyI,
                 min_pr_mu = -18,
                 max_pr_mu = 17,
                 seed = 1234,
                 max_treedepth = 12,
                 adapt_delta = 0.95,
                 file_id =  "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_fil",
                 iter = 6000)
stan_performance(df = filter(dat, species=="flo"),
                 response = RGRs,
                 treatment = dayTempc,
                 group_ids = familyI,
                 min_pr_mu = -18,
                 max_pr_mu = 17,
                 seed = 1234,
                 max_treedepth = 12,
                 adapt_delta = 0.95,
                 file_id =  "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_flo",
                 iter = 6000)
stan_performance(df = filter(dat, species=="gut"),
                 response = RGRs,
                 treatment = dayTempc,
                 group_ids = familyI,
                 min_pr_mu = -18,
                 max_pr_mu = 17,
                 seed = 1234,
                 max_treedepth = 12,
                 adapt_delta = 0.95,
                 file_id =  "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_gut",
                 iter = 6000)
stan_performance(df = filter(dat, species=="lac"),
                 response = RGRs,
                 treatment = dayTempc,
                 group_ids = familyI,
                 min_pr_mu = -18,
                 max_pr_mu = 17,
                 seed = 1234,
                 max_treedepth = 12,
                 adapt_delta = 0.95,
                 file_id =  "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_lac",
                 iter = 6000)
stan_performance(df = filter(dat, species=="nor"),
                 response = RGRs,
                 treatment = dayTempc,
                 group_ids = familyI,
                 min_pr_mu = -18,
                 max_pr_mu = 17,
                 seed = 1234,
                 max_treedepth = 12,
                 adapt_delta = 0.95,
                 file_id =  "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_nor",
                 iter = 6000)
stan_performance(df = filter(dat, species=="par"),
                 response = RGRs,
                 treatment = dayTempc,
                 group_ids = familyI,
                 min_pr_mu = -18,
                 max_pr_mu = 17,
                 seed = 1234,
                 max_treedepth = 12,
                 adapt_delta = 0.95,
                 file_id =  "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_par",
                 iter = 6000)
stan_performance(df = filter(dat, species=="ver"),
                 response = RGRs,
                 treatment = dayTempc,
                 group_ids = familyI,
                 min_pr_mu = -18,
                 max_pr_mu = 17,
                 seed = 1234,
                 max_treedepth = 12,
                 adapt_delta = 0.95,
                 file_id =  "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_ver",
                 iter = 6000)


############ 
############ THE FOLLOWING CODE IS FOR PROCESSING MODEL OUTPUT  
############ 




############
# write model summaries
############
model_fits_bic <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_bic.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_bic.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_bic.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_bic.samples_4.csv"))
write.csv(round(rstan::summary(model_fits_bic)$summary, digits=2), file="Analysis-output/TPC/Model-summaries/model-sxp-af-summary_bic.csv")
model_fits_car <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_car.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_car.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_car.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_car.samples_4.csv"))
write.csv(round(rstan::summary(model_fits_car)$summary, digits=2), file="Analysis-output/TPC/Model-summaries/model-sxp-af-summary_car.csv")
model_fits_eas <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_eas.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_eas.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_eas.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_eas.samples_4.csv"))
write.csv(round(rstan::summary(model_fits_eas)$summary, digits=2), file="Analysis-output/TPC/Model-summaries/model-sxp-af-summary_eas.csv")
model_fits_fil <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_fil.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_fil.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_fil.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_fil.samples_4.csv"))
write.csv(round(rstan::summary(model_fits_fil)$summary, digits=2), file="Analysis-output/TPC/Model-summaries/model-sxp-af-summary_fil.csv")
model_fits_flo <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_flo.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_flo.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_flo.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_flo.samples_4.csv"))
write.csv(round(rstan::summary(model_fits_flo)$summary, digits=2), file="Analysis-output/TPC/Model-summaries/model-sxp-af-summary_flo.csv")
model_fits_gut <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_gut.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_gut.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_gut.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_gut.samples_4.csv"))
write.csv(round(rstan::summary(model_fits_gut)$summary, digits=2), file="Analysis-output/TPC/Model-summaries/model-sxp-af-summary_gut.csv")
model_fits_lac <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_lac.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_lac.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_lac.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_lac.samples_4.csv"))
write.csv(round(rstan::summary(model_fits_lac)$summary, digits=2), file="Analysis-output/TPC/Model-summaries/model-sxp-af-summary_lac.csv")
model_fits_nor <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_nor.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_nor.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_nor.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_nor.samples_4.csv"))
write.csv(round(rstan::summary(model_fits_nor)$summary, digits=2), file="Analysis-output/TPC/Model-summaries/model-sxp-af-summary_nor.csv")
model_fits_par <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_par.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_par.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_par.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_par.samples_4.csv"))
write.csv(round(rstan::summary(model_fits_par)$summary, digits=2), file="Analysis-output/TPC/Model-summaries/model-sxp-af-summary_par.csv")
model_fits_ver <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_ver.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_ver.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_ver.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_ver.samples_4.csv"))
write.csv(round(rstan::summary(model_fits_ver)$summary, digits=2), file="Analysis-output/TPC/Model-summaries/model-sxp-af-summary_ver.csv")







                                    

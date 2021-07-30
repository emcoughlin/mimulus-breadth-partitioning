#### PROJECT: Performance curves of Mimulus species
#### PURPOSE: Calculate pairwise comparisons between species
#### and between families within each species
#### AUTHOR/DATE: RCW/2020-07-06

###### Warnings:
###### 1) This script requires a lot of computational power.
###### 2) The tidy_perf files required to run this code (lines 42-51) are not located in this git repo because they are too large. Please contact Seema Sheth (ssheth3@ncsu.edu) if you are interested in obtaining these files.



#load performr
devtools::install_github("silastittes/performr", local = FALSE, ref="zin", force=FALSE) 
# note: tried running regular models that do NOT account for zero inflation 
# (ref="zin" argument in the install_github function above), but these models fit 
# horribly for the MPC models. so we're sticking with zero-inflation models

#load other libraries 
library(dplyr)
library(devtools)
library(performr)
library(tidyverse)
library(ggridges)
library(ggpubr)
library(bayesplot)


# Other settings
options(mc.cores = parallel::detectCores())
extract <- rstan::extract





### species:
sp <- c("car","par","ver","eas","flo","nor","bic","fil","gut","lac")



### read in model draws
tidy_perf_bic <- read.csv("Analysis-output/TPC/Models/tidy_perf_sxp_af_bic.csv") 
tidy_perf_fil <- read.csv("Analysis-output/TPC/Models/tidy_perf_sxp_af_fil.csv") 
tidy_perf_car <- read.csv("Analysis-output/TPC/Models/tidy_perf_sxp_af_car.csv")
tidy_perf_par <- read.csv("Analysis-output/TPC/Models/tidy_perf_sxp_af_par.csv") 
tidy_perf_eas <- read.csv("Analysis-output/TPC/Models/tidy_perf_sxp_af_eas.csv") 
tidy_perf_ver <- read.csv("Analysis-output/TPC/Models/tidy_perf_sxp_af_ver.csv") 
tidy_perf_flo <- read.csv("Analysis-output/TPC/Models/tidy_perf_sxp_af_flo.csv") 
tidy_perf_nor <- read.csv("Analysis-output/TPC/Models/tidy_perf_sxp_af_nor.csv") 
tidy_perf_gut <- read.csv("Analysis-output/TPC/Models/tidy_perf_sxp_af_gut.csv") 
tidy_perf_lac <- read.csv("Analysis-output/TPC/Models/tidy_perf_sxp_af_lac.csv") 


# make new columns for B50 thresholds such than none are >50 or <15
tidy_perf_bic$tmax2 <- tidy_perf_bic$B50_high; tidy_perf_bic$tmax2[which(tidy_perf_bic$tmax2>50)] <- 50
tidy_perf_bic$tmin2 <- tidy_perf_bic$B50_low; tidy_perf_bic$tmin2[which(tidy_perf_bic$tmin2<15)] <- 15
tidy_perf_car$tmax2 <- tidy_perf_car$B50_high; tidy_perf_car$tmax2[which(tidy_perf_car$tmax2>50)] <- 50
tidy_perf_car$tmin2 <- tidy_perf_car$B50_low; tidy_perf_car$tmin2[which(tidy_perf_car$tmin2<15)] <- 15
tidy_perf_par$tmax2 <- tidy_perf_par$B50_high; tidy_perf_par$tmax2[which(tidy_perf_par$tmax2>50)] <- 50
tidy_perf_par$tmin2 <- tidy_perf_par$B50_low; tidy_perf_par$tmin2[which(tidy_perf_par$tmin2<15)] <- 15
tidy_perf_ver$tmax2 <- tidy_perf_ver$B50_high; tidy_perf_ver$tmax2[which(tidy_perf_ver$tmax2>50)] <- 50
tidy_perf_ver$tmin2 <- tidy_perf_ver$B50_low; tidy_perf_ver$tmin2[which(tidy_perf_ver$tmin2<15)] <- 15
tidy_perf_eas$tmax2 <- tidy_perf_eas$B50_high; tidy_perf_eas$tmax2[which(tidy_perf_eas$tmax2>50)] <- 50
tidy_perf_eas$tmin2 <- tidy_perf_eas$B50_low; tidy_perf_eas$tmin2[which(tidy_perf_eas$tmin2<15)] <- 15
tidy_perf_flo$tmax2 <- tidy_perf_flo$B50_high; tidy_perf_flo$tmax2[which(tidy_perf_flo$tmax2>50)] <- 50
tidy_perf_flo$tmin2 <- tidy_perf_flo$B50_low; tidy_perf_flo$tmin2[which(tidy_perf_flo$tmin2<15)] <- 15
tidy_perf_nor$tmax2 <- tidy_perf_nor$B50_high; tidy_perf_nor$tmax2[which(tidy_perf_nor$tmax2>50)] <- 50
tidy_perf_nor$tmin2 <- tidy_perf_nor$B50_low; tidy_perf_nor$tmin2[which(tidy_perf_nor$tmin2<15)] <- 15
tidy_perf_fil$tmax2 <- tidy_perf_fil$B50_high; tidy_perf_fil$tmax2[which(tidy_perf_fil$tmax2>50)] <- 50
tidy_perf_fil$tmin2 <- tidy_perf_fil$B50_low; tidy_perf_fil$tmin2[which(tidy_perf_fil$tmin2<15)] <- 15
tidy_perf_gut$tmax2 <- tidy_perf_gut$B50_high; tidy_perf_gut$tmax2[which(tidy_perf_gut$tmax2>50)] <- 50
tidy_perf_gut$tmin2 <- tidy_perf_gut$B50_low; tidy_perf_gut$tmin2[which(tidy_perf_gut$tmin2<15)] <- 15
tidy_perf_lac$tmax2 <- tidy_perf_lac$B50_high; tidy_perf_lac$tmax2[which(tidy_perf_lac$tmax2>50)] <- 50
tidy_perf_lac$tmin2 <- tidy_perf_lac$B50_low; tidy_perf_lac$tmin2[which(tidy_perf_lac$tmin2<15)] <- 15

# check
range(tidy_perf_car$tmin2); range(tidy_perf_car$tmax2)


# make new column for B50 calculated from tmin2 and tmax2
tidy_perf_bic$tbreadth <- tidy_perf_bic$tmax2 - tidy_perf_bic$tmin2
tidy_perf_car$tbreadth <- tidy_perf_car$tmax2 - tidy_perf_car$tmin2
tidy_perf_par$tbreadth <- tidy_perf_par$tmax2 - tidy_perf_par$tmin2
tidy_perf_ver$tbreadth <- tidy_perf_ver$tmax2 - tidy_perf_ver$tmin2
tidy_perf_eas$tbreadth <- tidy_perf_eas$tmax2 - tidy_perf_eas$tmin2
tidy_perf_flo$tbreadth <- tidy_perf_flo$tmax2 - tidy_perf_flo$tmin2
tidy_perf_nor$tbreadth <- tidy_perf_nor$tmax2 - tidy_perf_nor$tmin2
tidy_perf_fil$tbreadth <- tidy_perf_fil$tmax2 - tidy_perf_fil$tmin2
tidy_perf_gut$tbreadth <- tidy_perf_gut$tmax2 - tidy_perf_gut$tmin2
tidy_perf_lac$tbreadth <- tidy_perf_lac$tmax2 - tidy_perf_lac$tmin2


# check
range(tidy_perf_car$tbreadth)



# set specs for calculations
divby <- 1 # sample every nth iteration for credible intervals
ndraws <- max(tidy_perf_bic$draw)/divby



### data frame that we will put step-wise calculations into
tidy_sp <- data.frame(draw=rep(1:ndraws,10),
                      sp=rep(sp, each=ndraws),
                      popBreadth=rep(NA, ndraws*10),
                      meanFamBreadth=rep(NA, ndraws*10),
                      varTopt=rep(NA, ndraws*10),
                      varFamBreadth=rep(NA, ndraws*10)) 


### fill out step-wise calcuations for each species
for(i in 1:ndraws){ 
  # bic:
  maxnum <- max(tidy_perf_bic$tmax2[which(tidy_perf_bic$draw==i*divby)])
  minnum <- min(tidy_perf_bic$tmin2[which(tidy_perf_bic$draw==i*divby)])
  tidy_sp$popBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="bic")] <- maxnum-minnum
  tidy_sp$meanFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="bic")] <- mean(tidy_perf_bic$tbreadth[which(tidy_perf_bic$draw==i*divby)])
  tidy_sp$varTopt[which(tidy_sp$draw==i & tidy_sp$sp=="bic")] <- var(tidy_perf_bic$maximaBT[which(tidy_perf_bic$draw==i*divby)])
  tidy_sp$varFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="bic")] <- var(tidy_perf_bic$tbreadth[which(tidy_perf_bic$draw==i*divby)])

  # car:
  maxnum <- max(tidy_perf_car$tmax2[which(tidy_perf_car$draw==i*divby)])
  minnum <- min(tidy_perf_car$tmin2[which(tidy_perf_car$draw==i*divby)])
  tidy_sp$popBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="car")] <- maxnum-minnum
  tidy_sp$meanFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="car")] <- mean(tidy_perf_car$tbreadth[which(tidy_perf_car$draw==i*divby)])
  tidy_sp$varTopt[which(tidy_sp$draw==i & tidy_sp$sp=="car")] <- var(tidy_perf_car$maximaBT[which(tidy_perf_car$draw==i*divby)])
  tidy_sp$varFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="car")] <- var(tidy_perf_car$tbreadth[which(tidy_perf_car$draw==i*divby)])
  
  # par:
  maxnum <- max(tidy_perf_par$tmax2[which(tidy_perf_par$draw==i*divby)])
  minnum <- min(tidy_perf_par$tmin2[which(tidy_perf_par$draw==i*divby)])
  tidy_sp$popBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="par")] <- maxnum-minnum
  tidy_sp$meanFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="par")] <- mean(tidy_perf_par$tbreadth[which(tidy_perf_par$draw==i*divby)])
  tidy_sp$varTopt[which(tidy_sp$draw==i & tidy_sp$sp=="par")] <- var(tidy_perf_par$maximaBT[which(tidy_perf_par$draw==i*divby)])
  tidy_sp$varFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="par")] <- var(tidy_perf_par$tbreadth[which(tidy_perf_par$draw==i*divby)])
  
  # fil:
  maxnum <- max(tidy_perf_fil$tmax2[which(tidy_perf_fil$draw==i*divby)])
  minnum <- min(tidy_perf_fil$tmin2[which(tidy_perf_fil$draw==i*divby)])
  tidy_sp$popBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="fil")] <- maxnum-minnum
  tidy_sp$meanFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="fil")] <- mean(tidy_perf_fil$tbreadth[which(tidy_perf_fil$draw==i*divby)])
  tidy_sp$varTopt[which(tidy_sp$draw==i & tidy_sp$sp=="fil")] <- var(tidy_perf_fil$maximaBT[which(tidy_perf_fil$draw==i*divby)])
  tidy_sp$varFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="fil")] <- var(tidy_perf_fil$tbreadth[which(tidy_perf_fil$draw==i*divby)])
  
  # flo:
  maxnum <- max(tidy_perf_flo$tmax2[which(tidy_perf_flo$draw==i*divby)])
  minnum <- min(tidy_perf_flo$tmin2[which(tidy_perf_flo$draw==i*divby)])
  tidy_sp$popBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="flo")] <- maxnum-minnum
  tidy_sp$meanFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="flo")] <- mean(tidy_perf_flo$tbreadth[which(tidy_perf_flo$draw==i*divby)])
  tidy_sp$varTopt[which(tidy_sp$draw==i & tidy_sp$sp=="flo")] <- var(tidy_perf_flo$maximaBT[which(tidy_perf_flo$draw==i*divby)])
  tidy_sp$varFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="flo")] <- var(tidy_perf_flo$tbreadth[which(tidy_perf_flo$draw==i*divby)])
  
  # gut:
  maxnum <- max(tidy_perf_gut$tmax2[which(tidy_perf_gut$draw==i*divby)])
  minnum <- min(tidy_perf_gut$tmin2[which(tidy_perf_gut$draw==i*divby)])
  tidy_sp$popBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="gut")] <- maxnum-minnum
  tidy_sp$meanFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="gut")] <- mean(tidy_perf_gut$tbreadth[which(tidy_perf_gut$draw==i*divby)])
  tidy_sp$varTopt[which(tidy_sp$draw==i & tidy_sp$sp=="gut")] <- var(tidy_perf_gut$maximaBT[which(tidy_perf_gut$draw==i*divby)])
  tidy_sp$varFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="gut")] <- var(tidy_perf_gut$tbreadth[which(tidy_perf_gut$draw==i*divby)])
  
  # lac:
  maxnum <- max(tidy_perf_lac$tmax2[which(tidy_perf_lac$draw==i*divby)])
  minnum <- min(tidy_perf_lac$tmin2[which(tidy_perf_lac$draw==i*divby)])
  tidy_sp$popBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="lac")] <- maxnum-minnum
  tidy_sp$meanFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="lac")] <- mean(tidy_perf_lac$tbreadth[which(tidy_perf_lac$draw==i*divby)])
  tidy_sp$varTopt[which(tidy_sp$draw==i & tidy_sp$sp=="lac")] <- var(tidy_perf_lac$maximaBT[which(tidy_perf_lac$draw==i*divby)])
  tidy_sp$varFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="lac")] <- var(tidy_perf_lac$tbreadth[which(tidy_perf_lac$draw==i*divby)])
  
  # nor:
  maxnum <- max(tidy_perf_nor$tmax2[which(tidy_perf_nor$draw==i*divby)])
  minnum <- min(tidy_perf_nor$tmin2[which(tidy_perf_nor$draw==i*divby)])
  tidy_sp$popBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="nor")] <- maxnum-minnum
  tidy_sp$meanFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="nor")] <- mean(tidy_perf_nor$tbreadth[which(tidy_perf_nor$draw==i*divby)])
  tidy_sp$varTopt[which(tidy_sp$draw==i & tidy_sp$sp=="nor")] <- var(tidy_perf_nor$maximaBT[which(tidy_perf_nor$draw==i*divby)])
  tidy_sp$varFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="nor")] <- var(tidy_perf_nor$tbreadth[which(tidy_perf_nor$draw==i*divby)])
  
  # eas:
  maxnum <- max(tidy_perf_eas$tmax2[which(tidy_perf_eas$draw==i*divby)])
  minnum <- min(tidy_perf_eas$tmin2[which(tidy_perf_eas$draw==i*divby)])
  tidy_sp$popBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="eas")] <- maxnum-minnum
  tidy_sp$meanFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="eas")] <- mean(tidy_perf_eas$tbreadth[which(tidy_perf_eas$draw==i*divby)])
  tidy_sp$varTopt[which(tidy_sp$draw==i & tidy_sp$sp=="eas")] <- var(tidy_perf_eas$maximaBT[which(tidy_perf_eas$draw==i*divby)])
  tidy_sp$varFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="eas")] <- var(tidy_perf_eas$tbreadth[which(tidy_perf_eas$draw==i*divby)])
  
  # ver:
  maxnum <- max(tidy_perf_ver$tmax2[which(tidy_perf_ver$draw==i*divby)])
  minnum <- min(tidy_perf_ver$tmin2[which(tidy_perf_ver$draw==i*divby)])
  tidy_sp$popBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="ver")] <- maxnum-minnum
  tidy_sp$meanFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="ver")] <- mean(tidy_perf_ver$tbreadth[which(tidy_perf_ver$draw==i*divby)])
  tidy_sp$varTopt[which(tidy_sp$draw==i & tidy_sp$sp=="ver")] <- var(tidy_perf_ver$maximaBT[which(tidy_perf_ver$draw==i*divby)])
  tidy_sp$varFamBreadth[which(tidy_sp$draw==i & tidy_sp$sp=="ver")] <- var(tidy_perf_ver$tbreadth[which(tidy_perf_ver$draw==i*divby)])
  
  if(i %in% seq(1,ndraws, by=ndraws/10)){
    cat(round(i/ndraws,2)*100,"% done.", '\n')
  }
  
  }



# get means
tidy_sp_means <- tidy_sp %>%
  dplyr::group_by(sp) %>%
  dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))

# get lower confidence threshold
tidy_sp_lq <- tidy_sp %>%
  dplyr::group_by(sp) %>%
  dplyr::summarise(across(where(is.numeric), ~quantile(.x, 0.025)))

# get upper confidence threshold
tidy_sp_hq <- tidy_sp %>%
  dplyr::group_by(sp) %>%
  dplyr::summarise(across(where(is.numeric), ~quantile(.x, 0.975)))

# round
tidy_sp_means2 <- tidy_sp_means %>%
  dplyr::summarise(across(where(is.numeric), ~round(.x, 2)))
tidy_sp_lq2 <- tidy_sp_lq %>%
  dplyr::summarise(across(where(is.numeric), ~round(.x, 2)))
tidy_sp_hq2 <- tidy_sp_hq %>%
  dplyr::summarise(across(where(is.numeric), ~round(.x, 2)))

# transform to data frames
tidy_sp_means2 <- as.data.frame(tidy_sp_means2)
tidy_sp_lq2 <- as.data.frame(tidy_sp_lq2)
tidy_sp_hq2 <- as.data.frame(tidy_sp_hq2)

# apply species names to row names
rownames(tidy_sp_means2) <- tidy_sp_means$sp
rownames(tidy_sp_lq2) <- tidy_sp_lq$sp
rownames(tidy_sp_hq2) <- tidy_sp_hq$sp

### order by species:
sp <- c("car","par","ver","eas","flo","nor","bic","fil","gut","lac")
tidy_sp_means3 <- tidy_sp_means2[sp,]
tidy_sp_lq3 <- tidy_sp_lq2[sp,]
tidy_sp_hq3 <- tidy_sp_hq2[sp,]

# create sp column
tidy_sp_means3$sp <- as.factor(rownames(tidy_sp_means3))
tidy_sp_lq3$sp <- as.factor(rownames(tidy_sp_lq3))
tidy_sp_hq3$sp <- as.factor(rownames(tidy_sp_hq3))

# paste together
d <- data.frame(sp=sp,
                popBreadth=rep(NA, length(sp)),
                meanFamBreadth=rep(NA, length(sp)),
                varTopt=rep(NA, length(sp)),
                varFamBreadth=rep(NA, length(sp)))

# fill in mean & credible interval
for(i in 1:dim(d)[1]){
  d$popBreadth[i] <- paste0(tidy_sp_means3$popBreadth[which(tidy_sp_means3$sp==d$sp[i])],
                          " [",  tidy_sp_lq3$popBreadth[which(tidy_sp_lq3$sp==d$sp[i])],
                          ", ",  tidy_sp_hq3$popBreadth[which(tidy_sp_hq3$sp==d$sp[i])], "]")
  d$meanFamBreadth[i] <- paste0(tidy_sp_means3$meanFamBreadth[which(tidy_sp_means3$sp==d$sp[i])],
                          " [",  tidy_sp_lq3$meanFamBreadth[which(tidy_sp_lq3$sp==d$sp[i])],
                          ", ",  tidy_sp_hq3$meanFamBreadth[which(tidy_sp_hq3$sp==d$sp[i])], "]")
  d$varTopt[i] <- paste0(tidy_sp_means3$varTopt[which(tidy_sp_means3$sp==d$sp[i])],
                              " [",  tidy_sp_lq3$varTopt[which(tidy_sp_lq3$sp==d$sp[i])],
                              ", ",  tidy_sp_hq3$varTopt[which(tidy_sp_hq3$sp==d$sp[i])], "]")
  d$varFamBreadth[i] <- paste0(tidy_sp_means3$varFamBreadth[which(tidy_sp_means3$sp==d$sp[i])],
                            " [",  tidy_sp_lq3$varFamBreadth[which(tidy_sp_lq3$sp==d$sp[i])],
                            ", ",  tidy_sp_hq3$varFamBreadth[which(tidy_sp_hq3$sp==d$sp[i])], "]")
}


write.csv(d, "Analysis-output/TPC/Creds/Table-1_parameter-variance-calculations.csv")


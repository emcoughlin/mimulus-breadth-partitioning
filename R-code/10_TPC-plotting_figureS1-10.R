#### PROJECT: Performance curves of Mimulus species
#### PURPOSE: Process temperature performance models of 
#### families within each species
#### AUTHOR/DATE: RCW/2020-02-10


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
############ 
# read everything back in so that we can plot tpcs 
# [, -1] means that we read in without the first column, which we don't need
############
############ 

# load data for the model and create subsetted datasets
dat <- read.csv("Processed-data/temp.mod.dat_all-fams.csv")
bic_dat <- dat[which(dat$species=="bic"),]
car_dat <- dat[which(dat$species=="car"),]
eas_dat <- dat[which(dat$species=="eas"),]
fil_dat <- dat[which(dat$species=="fil"),]
flo_dat <- dat[which(dat$species=="flo"),]
gut_dat <- dat[which(dat$species=="gut"),]
lac_dat <- dat[which(dat$species=="lac"),]
nor_dat <- dat[which(dat$species=="nor"),]
par_dat <- dat[which(dat$species=="par"),]
ver_dat <- dat[which(dat$species=="ver"),]
meansTot <- read.csv("Processed-data/temp.meansTot_all-fams.csv")

creds_bic <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_bic.csv")[,-1]
mean_df_bic <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_bic.csv")[,-1] 
mean_df_ci_bic <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_bic.csv")[,-1] 

creds_car <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_car.csv")[,-1]
mean_df_car <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_car.csv")[,-1] 
mean_df_ci_car <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_car.csv")[,-1] 

creds_eas <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_eas.csv")[,-1]
mean_df_eas <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_eas.csv")[,-1] 
mean_df_ci_eas <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_eas.csv")[,-1] 

creds_fil <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_fil.csv")[,-1]
mean_df_fil <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_fil.csv")[,-1] 
mean_df_ci_fil <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_fil.csv")[,-1] 

creds_flo <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_flo.csv")[,-1]
mean_df_flo <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_flo.csv")[,-1] 
mean_df_ci_flo <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_flo.csv")[,-1] 

creds_gut <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_gut.csv")[,-1]
mean_df_gut <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_gut.csv")[,-1] 
mean_df_ci_gut <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_gut.csv")[,-1] 

creds_lac <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_lac.csv")[,-1]
mean_df_lac <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_lac.csv")[,-1] 
mean_df_ci_lac <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_lac.csv")[,-1] 

creds_nor <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_nor.csv")[,-1]
mean_df_nor <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_nor.csv")[,-1] 
mean_df_ci_nor <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_nor.csv")[,-1] 

creds_par <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_par.csv")[,-1]
mean_df_par <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_par.csv")[,-1] 
mean_df_ci_par <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_par.csv")[,-1] 

creds_ver <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_ver.csv")[,-1]
mean_df_ver <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ver.csv")[,-1] 
mean_df_ci_ver <- read.csv("Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_ver.csv")[,-1] 


####################
##### PLOT CURVES
####################

# This plot shows the tpc's for each family and 95% confidence intervals
# Repeat for each species
# Note: this is script used to make final plots for supporting info of manuscript as of 20210728

plwd <- 0.25 # parameter lwd
pla <- 1 # parameter line alpha
tpclwd <- 1.5 # tpc lwd
tpca <- 1 # tpc line alpha
blwd <- 0.5 # background line width
bla <- 0.5 # background line alpha

flo_dat$speciesF <- as.factor(flo_dat$familyI)
creds_flo$speciesF <- as.factor(creds_flo$species)
mean_df_flo$speciesF <- as.factor(mean_df_flo$species)
bayesFit_flo <- ggplot(data = filter(creds_flo, level == 95)) +
  #geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  geom_ribbon(aes(x=x, ymin=upper, ymax=lower, fill=speciesF),alpha=0.3, inherit.aes=F) + 
  geom_point(data=flo_dat, position=position_jitter(w = 0.1, h=0), shape=21, alpha = 0.35, aes(dayTemp, RGR, group=speciesF), inherit.aes=F, size=1.25) +
  labs(x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + 
  geom_line(aes(x, mu, group=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  facet_wrap(~speciesF, ncol=6) + guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_segment(data=mean_df_flo, aes(x=maximaBT, y=max_RGRBT, xend=maximaBT, yend=0,group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # vertical lines for thermal optima
  geom_segment(data=mean_df_flo, aes(x=B50_low, y=max_RGRBT*0.5, xend=B50_high, yend=max_RGRBT*0.5, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for B50 
  geom_segment(data=mean_df_flo, aes(x=maximaBT, y=max_RGRBT, xend=5, yend=max_RGRBT, group=speciesF), alpha=1, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for performance maximum
  geom_segment(data=mean_df_flo, aes(x=x_minBT, y=max(creds_flo$mu)/10, xend=x_minBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df_flo, aes(x=x_maxBT, y=max(creds_flo$mu)/10, xend=x_maxBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1)+ # upper thermal limit
  xlim(5,58) +
  ylim(-0.007,NA) +
  theme_set(theme_minimal()) +
  labs(title="M. floribundus") +
  theme(plot.title = element_text(face = "italic"))
bayesFit_flo

# export plot
pdf("Figures/TPC/flo-fams-sxp-af.pdf", height=6, width=10); figure <- ggarrange(bayesFit_flo, ncol=1, nrow=1); figure; dev.off()



nor_dat$speciesF <- as.factor(nor_dat$familyI)
creds_nor$speciesF <- as.factor(creds_nor$species)
mean_df_nor$speciesF <- as.factor(mean_df_nor$species)
bayesFit_nor <- ggplot(data = filter(creds_nor, level == 95)) +
  #geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  geom_ribbon(aes(x=x, ymin=upper, ymax=lower, fill=speciesF),alpha=0.3, inherit.aes=F) + 
  geom_point(data=nor_dat, position=position_jitter(w = 0.1, h=0), shape=21, alpha = 0.35, aes(dayTemp, RGR, group=speciesF), inherit.aes=F, size=1.25) +
  labs(x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + 
  geom_line(aes(x, mu, group=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  facet_wrap(~speciesF, ncol=6) + guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_segment(data=mean_df_nor, aes(x=maximaBT, y=max_RGRBT, xend=maximaBT, yend=0,group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # vertical lines for thermal optima
  geom_segment(data=mean_df_nor, aes(x=B50_low, y=max_RGRBT*0.5, xend=B50_high, yend=max_RGRBT*0.5, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for B50 
  geom_segment(data=mean_df_nor, aes(x=maximaBT, y=max_RGRBT, xend=5, yend=max_RGRBT, group=speciesF), alpha=1, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for performance maximum
  geom_segment(data=mean_df_nor, aes(x=x_minBT, y=max(creds_nor$mu)/10, xend=x_minBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df_nor, aes(x=x_maxBT, y=max(creds_nor$mu)/10, xend=x_maxBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1)+ # upper thermal limit
  xlim(5,58) +
  ylim(-0.007,NA) +
  theme_set(theme_minimal()) +
  labs(title="M. norrisii") +
  theme(plot.title = element_text(face = "italic"))
bayesFit_nor

# export plot
pdf("Figures/TPC/nor-fams-sxp-af.pdf", height=6, width=10); figure <- ggarrange(bayesFit_nor, ncol=1, nrow=1); figure; dev.off()


car_dat$speciesF <- as.factor(car_dat$familyI)
creds_car$speciesF <- as.factor(creds_car$species)
mean_df_car$speciesF <- as.factor(mean_df_car$species)
bayesFit_car <- ggplot(data = filter(creds_car, level == 95)) +
  #geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  geom_ribbon(aes(x=x, ymin=upper, ymax=lower, fill=speciesF),alpha=0.3, inherit.aes=F) + 
  geom_point(data=car_dat, position=position_jitter(w = 0.1, h=0), shape=21, alpha = 0.35, aes(dayTemp, RGR, group=speciesF), inherit.aes=F, size=1.25) +
  labs(x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + 
  geom_line(aes(x, mu, group=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  facet_wrap(~speciesF, ncol=6) + guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_segment(data=mean_df_car, aes(x=maximaBT, y=max_RGRBT, xend=maximaBT, yend=0,group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # vertical lines for thermal optima
  geom_segment(data=mean_df_car, aes(x=B50_low, y=max_RGRBT*0.5, xend=B50_high, yend=max_RGRBT*0.5, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for B50 
  geom_segment(data=mean_df_car, aes(x=maximaBT, y=max_RGRBT, xend=5, yend=max_RGRBT, group=speciesF), alpha=1, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for performance maximum
  geom_segment(data=mean_df_car, aes(x=x_minBT, y=max(creds_car$mu)/10, xend=x_minBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df_car, aes(x=x_maxBT, y=max(creds_car$mu)/10, xend=x_maxBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1)+ # upper thermal limit
  xlim(5,58) +
  ylim(-0.007,NA) +
  theme_set(theme_minimal()) +
  labs(title="M. cardinalis") +
  theme(plot.title = element_text(face = "italic"))
bayesFit_car

# export plot
pdf("Figures/TPC/car-fams-sxp-af.pdf", height=8, width=10); figure <- ggarrange(bayesFit_car, ncol=1, nrow=1); figure; dev.off()



par_dat$speciesF <- as.factor(par_dat$familyI)
creds_par$speciesF <- as.factor(creds_par$species)
mean_df_par$speciesF <- as.factor(mean_df_par$species)
bayesFit_par <- ggplot(data = filter(creds_par, level == 95)) +
  #geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  geom_ribbon(aes(x=x, ymin=upper, ymax=lower, fill=speciesF),alpha=0.3, inherit.aes=F) + 
  geom_point(data=par_dat, position=position_jitter(w = 0.1, h=0), shape=21, alpha = 0.35, aes(dayTemp, RGR, group=speciesF), inherit.aes=F, size=1.25) +
  labs(x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + 
  geom_line(aes(x, mu, group=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  facet_wrap(~speciesF, ncol=6) + guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_segment(data=mean_df_par, aes(x=maximaBT, y=max_RGRBT, xend=maximaBT, yend=0,group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # vertical lines for thermal optima
  geom_segment(data=mean_df_par, aes(x=B50_low, y=max_RGRBT*0.5, xend=B50_high, yend=max_RGRBT*0.5, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for B50 
  geom_segment(data=mean_df_par, aes(x=maximaBT, y=max_RGRBT, xend=5, yend=max_RGRBT, group=speciesF), alpha=1, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for performance maximum
  geom_segment(data=mean_df_par, aes(x=x_minBT, y=max(creds_par$mu)/10, xend=x_minBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df_par, aes(x=x_maxBT, y=max(creds_par$mu)/10, xend=x_maxBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1)+ # upper thermal limit
  xlim(5,58) +
  ylim(-0.007,NA) +
  theme_set(theme_minimal()) +
  labs(title="M. parishii") +
  theme(plot.title = element_text(face = "italic"))
bayesFit_par

# export plot
pdf("Figures/TPC/par-fams-sxp-af.pdf", height=18, width=10); figure <- ggarrange(bayesFit_par, ncol=1, nrow=1); figure; dev.off()



ver_dat$speciesF <- as.factor(ver_dat$familyI)
creds_ver$speciesF <- as.factor(creds_ver$species)
mean_df_ver$speciesF <- as.factor(mean_df_ver$species)
bayesFit_ver <- ggplot(data = filter(creds_ver, level == 95)) +
  #geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  geom_ribbon(aes(x=x, ymin=upper, ymax=lower, fill=speciesF),alpha=0.3, inherit.aes=F) + 
  geom_point(data=ver_dat, position=position_jitter(w = 0.1, h=0), shape=21, alpha = 0.35, aes(dayTemp, RGR, group=speciesF), inherit.aes=F, size=1.25) +
  labs(x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + 
  geom_line(aes(x, mu, group=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  facet_wrap(~speciesF, ncol=6) + guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_segment(data=mean_df_ver, aes(x=maximaBT, y=max_RGRBT, xend=maximaBT, yend=0,group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # vertical lines for thermal optima
  geom_segment(data=mean_df_ver, aes(x=B50_low, y=max_RGRBT*0.5, xend=B50_high, yend=max_RGRBT*0.5, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for B50 
  geom_segment(data=mean_df_ver, aes(x=maximaBT, y=max_RGRBT, xend=5, yend=max_RGRBT, group=speciesF), alpha=1, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for performance maximum
  geom_segment(data=mean_df_ver, aes(x=x_minBT, y=max(creds_ver$mu)/10, xend=x_minBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df_ver, aes(x=x_maxBT, y=max(creds_ver$mu)/10, xend=x_maxBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1)+ # upper thermal limit
  xlim(5,58) +
  ylim(-0.007,NA) +
  theme_set(theme_minimal()) +
  labs(title="M. verbenaceus") +
  theme(plot.title = element_text(face = "italic"))
bayesFit_ver

# export plot
pdf("Figures/TPC/ver-fams-sxp-af.pdf", height=8, width=10); figure <- ggarrange(bayesFit_ver, ncol=1, nrow=1); figure; dev.off()



eas_dat$speciesF <- as.factor(eas_dat$familyI)
creds_eas$speciesF <- as.factor(creds_eas$species)
mean_df_eas$speciesF <- as.factor(mean_df_eas$species)
bayesFit_eas <- ggplot(data = filter(creds_eas, level == 95)) +
  #geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  geom_ribbon(aes(x=x, ymin=upper, ymax=lower, fill=speciesF),alpha=0.3, inherit.aes=F) + 
  geom_point(data=eas_dat, position=position_jitter(w = 0.1, h=0), shape=21, alpha = 0.35, aes(dayTemp, RGR, group=speciesF), inherit.aes=F, size=1.25) +
  labs(x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + 
  geom_line(aes(x, mu, group=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  facet_wrap(~speciesF, ncol=6) + guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_segment(data=mean_df_eas, aes(x=maximaBT, y=max_RGRBT, xend=maximaBT, yend=0,group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # vertical lines for thermal optima
  geom_segment(data=mean_df_eas, aes(x=B50_low, y=max_RGRBT*0.5, xend=B50_high, yend=max_RGRBT*0.5, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for B50 
  geom_segment(data=mean_df_eas, aes(x=maximaBT, y=max_RGRBT, xend=5, yend=max_RGRBT, group=speciesF), alpha=1, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for performance maximum
  geom_segment(data=mean_df_eas, aes(x=x_minBT, y=max(creds_eas$mu)/10, xend=x_minBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df_eas, aes(x=x_maxBT, y=max(creds_eas$mu)/10, xend=x_maxBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1)+ # upper thermal limit
  xlim(5,58) +
  ylim(-0.007,NA) +
  theme_set(theme_minimal()) +
  labs(title="M. eastwoodiae") +
  theme(plot.title = element_text(face = "italic"))
bayesFit_eas

# export plot
pdf("Figures/TPC/eas-fams-sxp-af.pdf", height=14, width=10); figure <- ggarrange(bayesFit_eas, ncol=1, nrow=1); figure; dev.off()


gut_dat$speciesF <- as.factor(gut_dat$familyI)
creds_gut$speciesF <- as.factor(creds_gut$species)
mean_df_gut$speciesF <- as.factor(mean_df_gut$species)
bayesFit_gut <- ggplot(data = filter(creds_gut, level == 95)) +
  #geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  geom_ribbon(aes(x=x, ymin=upper, ymax=lower, fill=speciesF),alpha=0.3, inherit.aes=F) + 
  geom_point(data=gut_dat, position=position_jitter(w = 0.1, h=0), shape=21, alpha = 0.35, aes(dayTemp, RGR, group=speciesF), inherit.aes=F, size=1.25) +
  labs(x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") + 
  geom_line(aes(x, mu, group=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  facet_wrap(~speciesF, ncol=6) + guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_segment(data=mean_df_gut, aes(x=maximaBT, y=max_RGRBT, xend=maximaBT, yend=0,group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # guttical lines for thermal optima
  geom_segment(data=mean_df_gut, aes(x=B50_low, y=max_RGRBT*0.5, xend=B50_high, yend=max_RGRBT*0.5, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for B50 
  geom_segment(data=mean_df_gut, aes(x=maximaBT, y=max_RGRBT, xend=5, yend=max_RGRBT, group=speciesF), alpha=1, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for performance maximum
  geom_segment(data=mean_df_gut, aes(x=x_minBT, y=max(creds_gut$mu)/10, xend=x_minBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df_gut, aes(x=x_maxBT, y=max(creds_gut$mu)/10, xend=x_maxBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1)+ # upper thermal limit
  xlim(5,58) +
  ylim(-0.007,NA) +
  theme_set(theme_minimal()) +
  labs(title="M. guttatus") +
  theme(plot.title = element_text(face = "italic"))
bayesFit_gut

# export plot
pdf("Figures/TPC/gut-fams-sxp-af.pdf", height=4, width=10); figure <- ggarrange(bayesFit_gut, ncol=1, nrow=1); figure; dev.off()



lac_dat$speciesF <- as.factor(lac_dat$familyI)
creds_lac$speciesF <- as.factor(creds_lac$species)
mean_df_lac$speciesF <- as.factor(mean_df_lac$species)
bayesFit_lac <- ggplot(data = filter(creds_lac, level == 95)) +
  #geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  geom_ribbon(aes(x=x, ymin=upper, ymax=lower, fill=speciesF),alpha=0.3, inherit.aes=F) + 
  geom_point(data=lac_dat, position=position_jitter(w = 0.1, h=0), shape=21, alpha = 0.35, aes(dayTemp, RGR, group=speciesF), inherit.aes=F, size=1.25) +
  labs(x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") + 
  geom_line(aes(x, mu, group=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  facet_wrap(~speciesF, ncol=6) + guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_segment(data=mean_df_lac, aes(x=maximaBT, y=max_RGRBT, xend=maximaBT, yend=0,group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # vertical lines for thermal optima
  geom_segment(data=mean_df_lac, aes(x=B50_low, y=max_RGRBT*0.5, xend=B50_high, yend=max_RGRBT*0.5, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for B50 
  geom_segment(data=mean_df_lac, aes(x=maximaBT, y=max_RGRBT, xend=5, yend=max_RGRBT, group=speciesF), alpha=1, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for performance maximum
  geom_segment(data=mean_df_lac, aes(x=x_minBT, y=max(creds_lac$mu)/10, xend=x_minBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df_lac, aes(x=x_maxBT, y=max(creds_lac$mu)/10, xend=x_maxBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1)+ # upper thermal limit
  xlim(5,58) +
  ylim(-0.007,NA) +
  theme_set(theme_minimal()) +
  labs(title="M. laciniatus") +
  theme(plot.title = element_text(face = "italic"))
bayesFit_lac

# export plot
pdf("Figures/TPC/lac-fams-sxp-af.pdf", height=6, width=10); figure <- ggarrange(bayesFit_lac, ncol=1, nrow=1); figure; dev.off()


bic_dat$speciesF <- as.factor(bic_dat$familyI)
creds_bic$speciesF <- as.factor(creds_bic$species)
mean_df_bic$speciesF <- as.factor(mean_df_bic$species)
bayesFit_bic <- ggplot(data = filter(creds_bic, level == 95)) +
  #geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  geom_ribbon(aes(x=x, ymin=upper, ymax=lower, fill=speciesF),alpha=0.3, inherit.aes=F) + 
  geom_point(data=bic_dat, position=position_jitter(w = 0.1, h=0), shape=21, alpha = 0.35, aes(dayTemp, RGR, group=speciesF), inherit.aes=F, size=1.25) +
  labs(x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") + 
  geom_line(aes(x, mu, group=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  facet_wrap(~speciesF, ncol=6) + guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_segment(data=mean_df_bic, aes(x=maximaBT, y=max_RGRBT, xend=maximaBT, yend=0,group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # bictical lines for thermal optima
  geom_segment(data=mean_df_bic, aes(x=B50_low, y=max_RGRBT*0.5, xend=B50_high, yend=max_RGRBT*0.5, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for B50 
  geom_segment(data=mean_df_bic, aes(x=maximaBT, y=max_RGRBT, xend=5, yend=max_RGRBT, group=speciesF), alpha=1, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for performance maximum
  geom_segment(data=mean_df_bic, aes(x=x_minBT, y=max(creds_bic$mu)/10, xend=x_minBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df_bic, aes(x=x_maxBT, y=max(creds_bic$mu)/10, xend=x_maxBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1)+ # upper thermal limit
  xlim(5,58) +
  ylim(-0.007,NA) +
  theme_set(theme_minimal()) +
  labs(title="M. bicolor") +
  theme(plot.title = element_text(face = "italic"))
bayesFit_bic

# export plot
pdf("Figures/TPC/bic-fams-sxp-af.pdf", height=8, width=10); figure <- ggarrange(bayesFit_bic, ncol=1, nrow=1); figure; dev.off()



fil_dat$speciesF <- as.factor(fil_dat$familyI)
creds_fil$speciesF <- as.factor(creds_fil$species)
mean_df_fil$speciesF <- as.factor(mean_df_fil$species)
bayesFit_fil <- ggplot(data = filter(creds_fil, level == 95)) +
  #geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  geom_ribbon(aes(x=x, ymin=upper, ymax=lower, fill=speciesF),alpha=0.3, inherit.aes=F) + 
  geom_point(data=fil_dat, position=position_jitter(w = 0.1, h=0), shape=21, alpha = 0.35, aes(dayTemp, RGR, group=speciesF), inherit.aes=F, size=1.25) +
  labs(x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") + 
  geom_line(aes(x, mu, group=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  facet_wrap(~speciesF, ncol=6) + guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_segment(data=mean_df_fil, aes(x=maximaBT, y=max_RGRBT, xend=maximaBT, yend=0,group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # vertical lines for thermal optima
  geom_segment(data=mean_df_fil, aes(x=B50_low, y=max_RGRBT*0.5, xend=B50_high, yend=max_RGRBT*0.5, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for B50 
  geom_segment(data=mean_df_fil, aes(x=maximaBT, y=max_RGRBT, xend=5, yend=max_RGRBT, group=speciesF), alpha=1, lwd=plwd, inherit.aes=FALSE, lty=1) + # horizontal line for performance maximum
  geom_segment(data=mean_df_fil, aes(x=x_minBT, y=max(creds_fil$mu)/10, xend=x_minBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1) + # lower thermal limit
  geom_segment(data=mean_df_fil, aes(x=x_maxBT, y=max(creds_fil$mu)/10, xend=x_maxBT, yend=0, group=speciesF), alpha=pla, lwd=plwd, inherit.aes=FALSE, lty=1)+ # upper thermal limit
  xlim(5,58) +
  ylim(-0.007,NA) +
  theme_set(theme_minimal()) +
  labs(title="M. filicaulis") +
  theme(plot.title = element_text(face = "italic"))
bayesFit_fil

# export plot
pdf("Figures/TPC/fil-fams-sxp-af.pdf", height=6, width=10); figure <- ggarrange(bayesFit_fil, ncol=1, nrow=1); figure; dev.off()





### The following plots show the tpc's for each family with global curve

#######
# bicolor
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_bic.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_tpc <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_tpcBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="bic")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="bic")])
bic_dat$speciesF <- as.factor(bic_dat$familyI)
creds_bic$speciesF <- as.factor(creds_bic$species)
mean_df_bic$speciesF <- as.factor(mean_df_bic$species)
bayesFit_bic2 <- ggplot(data = filter(creds_bic, level == 95)) +
  geom_point(data=bic_dat, position=position_jitter(w = 0.5), shape=21, alpha = 0.35, aes(dayTemp, RGR, color=speciesF), inherit.aes=F, size=1.25) +
  labs(title="M. bicolor", x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") + xlim(5,55) + ylim(0,max(na.omit(c(fil_dat$RGR,bic_dat$RGR)))) +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  geom_line(data=mu_tpcBT, aes(x, y), inherit.aes = F, lwd = 2) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(data=mean_df_bic, aes(x=maximaBT, y=0,color=speciesF), inherit.aes=FALSE) + # thermal optima
  geom_point(data=mean_df_bic, aes(x=5, y=max_RGRBT, color=speciesF), inherit.aes=FALSE)  # performance maximum

#######
# filicaulis
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_fil.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_tpc <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_tpcBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="fil")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="fil")])
fil_dat$speciesF <- as.factor(fil_dat$familyI)
creds_fil$speciesF <- as.factor(creds_fil$species)
mean_df_fil$speciesF <- as.factor(mean_df_fil$species)
bayesFit_fil2 <- ggplot(data = filter(creds_fil, level == 95)) +
  geom_point(data=fil_dat, position=position_jitter(w = 0.5), shape=21, alpha = 0.35, aes(dayTemp, RGR, color=speciesF), inherit.aes=F, size=1.25) +
  labs(title="M. filicaulis", x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") + xlim(5,55) + ylim(0,max(na.omit(c(fil_dat$RGR,bic_dat$RGR)))) +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  geom_line(data=mu_tpcBT, aes(x, y), inherit.aes = F, lwd = 2) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(data=mean_df_fil, aes(x=maximaBT, y=0,color=speciesF), inherit.aes=FALSE) + # thermal optima
  geom_point(data=mean_df_fil, aes(x=5, y=max_RGRBT, color=speciesF), inherit.aes=FALSE)  # performance maximum

#######
# guttatus
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_gut.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_tpc <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_tpcBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="gut")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="gut")])
gut_dat$speciesF <- as.factor(gut_dat$familyI)
creds_gut$speciesF <- as.factor(creds_gut$species)
mean_df_gut$speciesF <- as.factor(mean_df_gut$species)
bayesFit_gut2 <- ggplot(data = filter(creds_gut, level == 95)) +
  geom_point(data=gut_dat, position=position_jitter(w = 0.5), shape=21, alpha = 0.35, aes(dayTemp, RGR, color=speciesF), inherit.aes=F, size=1.25) +
  labs(title="M. guttatus", x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") + xlim(5,55) + ylim(0,max(na.omit(c(lac_dat$RGR,gut_dat$RGR)))) +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  geom_line(data=mu_tpcBT, aes(x, y), inherit.aes = F, lwd = 2) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(data=mean_df_gut, aes(x=maximaBT, y=0,color=speciesF), inherit.aes=FALSE) + # thermal optima
  geom_point(data=mean_df_gut, aes(x=5, y=max_RGRBT, color=speciesF), inherit.aes=FALSE)  # performance maximum

#######
# laciniatus
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_lac.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_tpc <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_tpcBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="lac")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="lac")])
lac_dat$speciesF <- as.factor(lac_dat$familyI)
creds_lac$speciesF <- as.factor(creds_lac$species)
mean_df_lac$speciesF <- as.factor(mean_df_lac$species)
bayesFit_lac2 <- ggplot(data = filter(creds_lac, level == 95)) +
  geom_point(data=lac_dat, position=position_jitter(w = 0.5), shape=21, alpha = 0.35, aes(dayTemp, RGR, color=speciesF), inherit.aes=F, size=1.25) +
  labs(title="M. laciniatus", x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") + xlim(5,55) + ylim(0,max(na.omit(c(lac_dat$RGR,gut_dat$RGR)))) +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  geom_line(data=mu_tpcBT, aes(x, y), inherit.aes = F, lwd = 2) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(data=mean_df_lac, aes(x=maximaBT, y=0,color=speciesF), inherit.aes=FALSE) + # thermal optima
  geom_point(data=mean_df_lac, aes(x=5, y=max_RGRBT, color=speciesF), inherit.aes=FALSE)  # performance maximum

#######
# cardinalis
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_car.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_tpc <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_tpcBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="car")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="car")])
car_dat$speciesF <- as.factor(car_dat$familyI)
creds_car$speciesF <- as.factor(creds_car$species)
mean_df_car$speciesF <- as.factor(mean_df_car$species)
bayesFit_car2 <- ggplot(data = filter(creds_car, level == 95)) +
  geom_point(data=car_dat, position=position_jitter(w = 0.5), shape=21, alpha = 0.35, aes(dayTemp, RGR, color=speciesF), inherit.aes=F, size=1.25) +
  labs(title="M. cardinalis", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + xlim(5,55) + ylim(0,max(na.omit(c(car_dat$RGR,par_dat$RGR)))) +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  geom_line(data=mu_tpcBT, aes(x, y), inherit.aes = F, lwd = 2) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(data=mean_df_car, aes(x=maximaBT, y=0,color=speciesF), inherit.aes=FALSE) + # thermal optima
  geom_point(data=mean_df_car, aes(x=5, y=max_RGRBT, color=speciesF), inherit.aes=FALSE)  # performance maximum

#######
# parishii
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_par.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_tpc <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_tpcBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="par")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="par")])
par_dat$speciesF <- as.factor(par_dat$familyI)
creds_par$speciesF <- as.factor(creds_par$species)
mean_df_par$speciesF <- as.factor(mean_df_par$species)
bayesFit_par2 <- ggplot(data = filter(creds_par, level == 95)) +
  geom_point(data=par_dat, position=position_jitter(w = 0.5), shape=21, alpha = 0.35, aes(dayTemp, RGR, color=speciesF), inherit.aes=F, size=1.25) +
  labs(title="M. parishii", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + xlim(5,55) + ylim(0,max(na.omit(c(car_dat$RGR,par_dat$RGR)))) +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  geom_line(data=mu_tpcBT, aes(x, y), inherit.aes = F, lwd = 2) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(data=mean_df_par, aes(x=maximaBT, y=0,color=speciesF), inherit.aes=FALSE) + # thermal optima
  geom_point(data=mean_df_par, aes(x=5, y=max_RGRBT, color=speciesF), inherit.aes=FALSE)  # performance maximum

#######
# verbenaceus
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_ver.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_tpc <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_tpcBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="ver")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="ver")])
ver_dat$speciesF <- as.factor(ver_dat$familyI)
creds_ver$speciesF <- as.factor(creds_ver$species)
mean_df_ver$speciesF <- as.factor(mean_df_ver$species)
bayesFit_ver2 <- ggplot(data = filter(creds_ver, level == 95)) +
  geom_point(data=ver_dat, position=position_jitter(w = 0.5), shape=21, alpha = 0.35, aes(dayTemp, RGR, color=speciesF), inherit.aes=F, size=1.25) +
  labs(title="M. verbenaceus", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + xlim(5,55) + ylim(0,max(na.omit(c(ver_dat$RGR,eas_dat$RGR)))) +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  geom_line(data=mu_tpcBT, aes(x, y), inherit.aes = F, lwd = 2) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(data=mean_df_ver, aes(x=maximaBT, y=0,color=speciesF), inherit.aes=FALSE) + # thermal optima
  geom_point(data=mean_df_ver, aes(x=5, y=max_RGRBT, color=speciesF), inherit.aes=FALSE)  # performance maximum

#######
# eastwoodiae
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_eas.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_tpc <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_tpcBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="eas")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="eas")])
eas_dat$speciesF <- as.factor(eas_dat$familyI)
creds_eas$speciesF <- as.factor(creds_eas$species)
mean_df_eas$speciesF <- as.factor(mean_df_eas$species)
bayesFit_eas2 <- ggplot(data = filter(creds_eas, level == 95)) +
  geom_point(data=eas_dat, position=position_jitter(w = 0.5), shape=21, alpha = 0.35, aes(dayTemp, RGR, color=speciesF), inherit.aes=F, size=1.25) +
  labs(title="M. eastwoodiae", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + xlim(5,55) + ylim(0,max(na.omit(c(ver_dat$RGR,eas_dat$RGR)))) +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  geom_line(data=mu_tpcBT, aes(x, y), inherit.aes = F, lwd = 2) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(data=mean_df_eas, aes(x=maximaBT, y=0,color=speciesF), inherit.aes=FALSE) + # thermal optima
  geom_point(data=mean_df_eas, aes(x=5, y=max_RGRBT, color=speciesF), inherit.aes=FALSE)  # performance maximum

#######
# floribundus
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_flo.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_tpc <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_tpcBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="flo")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="flo")])
flo_dat$speciesF <- as.factor(flo_dat$familyI)
creds_flo$speciesF <- as.factor(creds_flo$species)
mean_df_flo$speciesF <- as.factor(mean_df_flo$species)
bayesFit_flo2 <- ggplot(data = filter(creds_flo, level == 95)) +
  geom_point(data=flo_dat, position=position_jitter(w = 0.5), shape=21, alpha = 0.35, aes(dayTemp, RGR, color=speciesF), inherit.aes=F, size=1.25) +
  labs(title="M. floribundus", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  geom_line(data=mu_tpcBT, aes(x, y), inherit.aes = F, lwd = 2) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(data=mean_df_flo, aes(x=maximaBT, y=0,color=speciesF), inherit.aes=FALSE) + # thermal optima
  geom_point(data=mean_df_flo, aes(x=5, y=max_RGRBT, color=speciesF), inherit.aes=FALSE) +  # performance maximum
  lims(x=c(5,55), y=range(dat$RGR[which(dat$species %in% c("flo", "nor"))]))

#######
# norrisii
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_nor.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_tpc <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_tpcBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="nor")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="nor")])
nor_dat$speciesF <- as.factor(nor_dat$familyI)
creds_nor$speciesF <- as.factor(creds_nor$species)
mean_df_nor$speciesF <- as.factor(mean_df_nor$species)
bayesFit_nor2 <- ggplot(data = filter(creds_nor, level == 95)) +
  geom_point(data=nor_dat, position=position_jitter(w = 0.5), shape=21, alpha = 0.35, aes(dayTemp, RGR, color=speciesF), inherit.aes=F, size=1.25) +
  labs(title="M. norrisii", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + 
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 1, alpha=0.7) +
  geom_line(data=mu_tpcBT, aes(x, y), inherit.aes = F, lwd = 2) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(data=mean_df_nor, aes(x=maximaBT, y=0,color=speciesF), inherit.aes=FALSE) + # thermal optima
  geom_point(data=mean_df_nor, aes(x=5, y=max_RGRBT, color=speciesF), inherit.aes=FALSE) +  # performance maximum
  lims(x=c(5,55), y=range(dat$RGR[which(dat$species %in% c("flo", "nor"))]))







# export plot
sp_pairs <- c("bic-fil", "car-par", "eas-ver", "flo-nor", "gut-lac")
pdf("Figures/TPC/bic-fil-sxp-af.pdf", height=5, width=9); figure <- ggarrange(bayesFit_bic2, bayesFit_fil2, ncol=2, nrow=1); figure; dev.off()
pdf("Figures/TPC/car-par-sxp-af.pdf", height=5, width=9); figure <- ggarrange(bayesFit_car2, bayesFit_par2, ncol=2, nrow=1); figure; dev.off()
pdf("Figures/TPC/eas-ver-sxp-af.pdf", height=5, width=9); figure <- ggarrange(bayesFit_eas2, bayesFit_ver2, ncol=2, nrow=1); figure; dev.off()
pdf("Figures/TPC/flo-nor-sxp-af.pdf", height=5, width=9); figure <- ggarrange(bayesFit_flo2, bayesFit_nor2, ncol=2, nrow=1); figure; dev.off()
pdf("Figures/TPC/gut-lac-sxp-af.pdf", height=5, width=9); figure <- ggarrange(bayesFit_gut2, bayesFit_lac2, ncol=2, nrow=1); figure; dev.off()
pdf("Figures/TPC/all-species-sxp-af.pdf", height=15, width=6); figure <- ggarrange(bayesFit_bic2, bayesFit_fil2,
                                                                            bayesFit_car2, bayesFit_par2,
                                                                            bayesFit_ver2, bayesFit_eas2, 
                                                                            bayesFit_flo2, bayesFit_nor2,
                                                                            bayesFit_gut2, bayesFit_lac2, 
                                                                            ncol=2, nrow=5); figure; dev.off()











#######################
# figures for lab group presentation July 10, 2020


#######
# floribundus
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_flo.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_TPC <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_TPCBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="flo")], y=mu_TPC*meansTot$RGR[which(meansTot$species=="flo")])
flo_dat$speciesF <- as.factor(flo_dat$familyI)
creds_flo$speciesF <- as.factor(creds_flo$species)
mean_df_flo$speciesF <- as.factor(mean_df_flo$species)
bayesFit_flo3 <- ggplot(data = filter(creds_flo, level == 95)) +
  labs(title="M. floribundus", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 0.75, alpha=0.5) +
  scale_color_manual(values = rep("red", length(unique(creds_flo$speciesF)))) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  lims(x=c(5,60), y=c(0,0.3))

#######
# norrisii
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-afsummary_nor.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_TPC <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_TPCBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="nor")], y=mu_TPC*meansTot$RGR[which(meansTot$species=="nor")])
nor_dat$speciesF <- as.factor(nor_dat$familyI)
creds_nor$speciesF <- as.factor(creds_nor$species)
mean_df_nor$speciesF <- as.factor(mean_df_nor$species)
bayesFit_nor3 <- ggplot(data = filter(creds_nor, level == 95)) +
  labs(title="M. norrisii", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + 
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 0.75, alpha=0.5) +
  scale_color_manual(values = rep("blue", length(unique(creds_nor$speciesF)))) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  lims(x=c(5,60), y=c(0,0.3))







#######
# cardinalis
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-afsummary_car.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_TPC <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_TPCBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="car")], y=mu_TPC*meansTot$RGR[which(meansTot$species=="car")])
car_dat$speciesF <- as.factor(car_dat$familyI)
creds_car$speciesF <- as.factor(creds_car$species)
mean_df_car$speciesF <- as.factor(mean_df_car$species)
bayesFit_car3 <- ggplot(data = filter(creds_car, level == 95)) +
  labs(title="M. cardinalis", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 0.75, alpha=0.5) +
  scale_color_manual(values = rep("red", length(unique(creds_car$speciesF)))) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  lims(x=c(5,60), y=c(0,0.15))

#######
# parishii
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-afsummary_par.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_TPC <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_TPCBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="par")], y=mu_TPC*meansTot$RGR[which(meansTot$species=="par")])
par_dat$speciesF <- as.factor(par_dat$familyI)
creds_par$speciesF <- as.factor(creds_par$species)
mean_df_par$speciesF <- as.factor(mean_df_par$species)
bayesFit_par3 <- ggplot(data = filter(creds_par, level == 95)) +
  labs(title="M. parishii", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + 
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 0.75, alpha=0.5) +
  scale_color_manual(values = rep("blue", length(unique(creds_par$speciesF)))) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  lims(x=c(5,60), y=c(0,0.15))









#######
# verbenaceus
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-afsummary_ver.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_TPC <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_TPCBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="ver")], y=mu_TPC*meansTot$RGR[which(meansTot$species=="ver")])
ver_dat$speciesF <- as.factor(ver_dat$familyI)
creds_ver$speciesF <- as.factor(creds_ver$species)
mean_df_ver$speciesF <- as.factor(mean_df_ver$species)
bayesFit_ver3 <- ggplot(data = filter(creds_ver, level == 95)) +
  labs(title="M. verbenaceus", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 0.75, alpha=0.5) +
  scale_color_manual(values = rep("red", length(unique(creds_ver$speciesF)))) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  lims(x=c(5,60), y=c(0,0.15))

#######
# eastwoodiae
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-afsummary_eas.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_TPC <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_TPCBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="eas")], y=mu_TPC*meansTot$RGR[which(meansTot$species=="eas")])
eas_dat$speciesF <- as.factor(eas_dat$familyI)
creds_eas$speciesF <- as.factor(creds_eas$species)
mean_df_eas$speciesF <- as.factor(mean_df_eas$species)
bayesFit_eas3 <- ggplot(data = filter(creds_eas, level == 95)) +
  labs(title="M. eastwoodiae", x=expression(paste("Temperature (°C)")), y="RGR (number/number/day)") + 
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 0.75, alpha=0.5) +
  scale_color_manual(values = rep("blue", length(unique(creds_eas$speciesF)))) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  lims(x=c(5,60), y=c(0,0.15))








#######
# bicolor
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-afsummary_bic.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_TPC <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_TPCBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="bic")], y=mu_TPC*meansTot$RGR[which(meansTot$species=="bic")])
bic_dat$speciesF <- as.factor(bic_dat$familyI)
creds_bic$speciesF <- as.factor(creds_bic$species)
mean_df_bic$speciesF <- as.factor(mean_df_bic$species)
bayesFit_bic3 <- ggplot(data = filter(creds_bic, level == 95)) +
  labs(title="M. bicolor", x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 0.75, alpha=0.5) +
  scale_color_manual(values = rep("red", length(unique(creds_bic$speciesF)))) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  lims(x=c(5,60), y=c(0,0.45))

#######
# filicaulis
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-afsummary_fil.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_TPC <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_TPCBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="fil")], y=mu_TPC*meansTot$RGR[which(meansTot$species=="fil")])
fil_dat$speciesF <- as.factor(fil_dat$familyI)
creds_fil$speciesF <- as.factor(creds_fil$species)
mean_df_fil$speciesF <- as.factor(mean_df_fil$species)
bayesFit_fil3 <- ggplot(data = filter(creds_fil, level == 95)) +
  labs(title="M. filicaulis", x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") + 
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 0.75, alpha=0.5) +
  scale_color_manual(values = rep("blue", length(unique(creds_fil$speciesF)))) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  lims(x=c(5,60), y=c(0,0.45))








#######
# guttatus
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-afsummary_gut.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_TPC <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_TPCBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="gut")], y=mu_TPC*meansTot$RGR[which(meansTot$species=="gut")])
gut_dat$speciesF <- as.factor(gut_dat$familyI)
creds_gut$speciesF <- as.factor(creds_gut$species)
mean_df_gut$speciesF <- as.factor(mean_df_gut$species)
bayesFit_gut3 <- ggplot(data = filter(creds_gut, level == 95)) +
  labs(title="M. guttatus", x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") +
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 0.75, alpha=0.5) +
  scale_color_manual(values = rep("red", length(unique(creds_gut$speciesF)))) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  lims(x=c(5,60), y=c(0,1.3))

#######
# laciniatus
modsum <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-afsummary_lac.csv")
xs <- seq(modsum$mean[which(modsum$X=="mu_min")], modsum$mean[which(modsum$X=="mu_max")], length.out=100)
mu_TPC <- performance_mu(xs = xs, 
                         shape1 = modsum$mean[which(modsum$X=="mu_shape1")], 
                         shape2 = modsum$mean[which(modsum$X=="mu_shape2")], 
                         stretch = modsum$mean[which(modsum$X=="mu_stretch")], 
                         x_min = modsum$mean[which(modsum$X=="mu_min")], 
                         x_max = modsum$mean[which(modsum$X=="mu_max")])
mu_TPCBT <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="lac")], y=mu_TPC*meansTot$RGR[which(meansTot$species=="lac")])
lac_dat$speciesF <- as.factor(lac_dat$familyI)
creds_lac$speciesF <- as.factor(creds_lac$species)
mean_df_lac$speciesF <- as.factor(mean_df_lac$species)
bayesFit_lac3 <- ggplot(data = filter(creds_lac, level == 95)) +
  labs(title="M. laciniatus", x=expression(paste("Temperature (°C)")), y="RGR (cm/cm/day)") + 
  geom_line(aes(x, mu, color=speciesF), inherit.aes = F, lwd = 0.75, alpha=0.5) +
  scale_color_manual(values = rep("blue", length(unique(creds_lac$speciesF)))) +
  geom_hline(yintercept=0, lwd=1, lty=1, alpha=1) + 
  guides(fill=FALSE, color=FALSE) + theme_bw() + 
  theme(plot.title = element_text(size=14, face="italic", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=13), strip.background=element_rect(colour=NA,fill=NA), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  lims(x=c(5,60), y=c(0,1.3))








pdf("Figures/TPC/all-spp_3-sxp-af.pdf", height=12, width=8)
figure <- ggarrange(bayesFit_car3, bayesFit_par3, 
                    bayesFit_ver3, bayesFit_eas3, 
                    bayesFit_flo3, bayesFit_nor3,
                    bayesFit_bic3, bayesFit_fil3,
                    bayesFit_gut3, bayesFit_lac3,
                    ncol=2, nrow=5); figure; dev.off()











#### PROJECT: Mimulus niche breadth partitioning
#### PURPOSE: Plot thermal performance curves of families within each species
#### AUTHOR/LAST UPDATED: RCW & EMC/2021-07-29


#load performr
devtools::install_github("silastittes/performr", local = FALSE, ref="zin", force=FALSE) 

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



############ 
# load in family temperature curves
creds_bic_t <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_bic.csv")[,-1]
creds_car_t <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_car.csv")[,-1]
creds_eas_t <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_eas.csv")[,-1]
creds_fil_t <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_fil.csv")[,-1]
creds_flo_t <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_flo.csv")[,-1]
creds_gut_t <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_gut.csv")[,-1]
creds_lac_t <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_lac.csv")[,-1]
creds_nor_t <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_nor.csv")[,-1]
creds_par_t <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_par.csv")[,-1]
creds_ver_t <- read.csv("Analysis-output/TPC/Creds/creds_sxp_af_ver.csv")[,-1]

############ 
# load in family temperature curve summaries (this is for plotting population-level curve)
sum_bic_t <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_bic.csv")[,-1]
sum_car_t <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_car.csv")[,-1]
sum_eas_t <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_eas.csv")[,-1]
sum_fil_t <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_fil.csv")[,-1]
sum_flo_t <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_flo.csv")[,-1]
sum_gut_t <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_gut.csv")[,-1]
sum_lac_t <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_lac.csv")[,-1]
sum_nor_t <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_nor.csv")[,-1]
sum_par_t <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_par.csv")[,-1]
sum_ver_t <- read.csv("Analysis-output/TPC/Model-summaries/model-sxp-af-summary_ver.csv")[,-1]


############ 
# backtransform family temperature curve summaries 
meansTot <- read.csv("Processed-data/temp.meansTot_all-fams.csv")

numfams <- length(unique(creds_bic_t$species))
mus <- sum_bic_t$mean[1:5+(numfams*6)]
xs <- seq(mus[4], mus[5], length.out=100)
mu_tpc <- performance_mu(xs = xs, shape1 = mus[1], shape2 = mus[2], stretch = mus[3], x_min = mus[4], x_max = mus[5])
mu_tpcBT_bic <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="bic")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="bic")])

numfams <- length(unique(creds_car_t$species))
mus <- sum_car_t$mean[1:5+(numfams*6)]
xs <- seq(mus[4], mus[5], length.out=100)
mu_tpc <- performance_mu(xs = xs, shape1 = mus[1], shape2 = mus[2], stretch = mus[3], x_min = mus[4], x_max = mus[5])
mu_tpcBT_car <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="car")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="car")])

numfams <- length(unique(creds_eas_t$species))
mus <- sum_eas_t$mean[1:5+(numfams*6)]
xs <- seq(mus[4], mus[5], length.out=100)
mu_tpc <- performance_mu(xs = xs, shape1 = mus[1], shape2 = mus[2], stretch = mus[3], x_min = mus[4], x_max = mus[5])
mu_tpcBT_eas <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="eas")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="eas")])

numfams <- length(unique(creds_fil_t$species))
mus <- sum_fil_t$mean[1:5+(numfams*6)]
xs <- seq(mus[4], mus[5], length.out=100)
mu_tpc <- performance_mu(xs = xs, shape1 = mus[1], shape2 = mus[2], stretch = mus[3], x_min = mus[4], x_max = mus[5])
mu_tpcBT_fil <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="fil")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="fil")])

numfams <- length(unique(creds_flo_t$species))
mus <- sum_flo_t$mean[1:5+(numfams*6)]
xs <- seq(mus[4], mus[5], length.out=100)
mu_tpc <- performance_mu(xs = xs, shape1 = mus[1], shape2 = mus[2], stretch = mus[3], x_min = mus[4], x_max = mus[5])
mu_tpcBT_flo <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="flo")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="flo")])

numfams <- length(unique(creds_gut_t$species))
mus <- sum_gut_t$mean[1:5+(numfams*6)]
xs <- seq(mus[4], mus[5], length.out=100)
mu_tpc <- performance_mu(xs = xs, shape1 = mus[1], shape2 = mus[2], stretch = mus[3], x_min = mus[4], x_max = mus[5])
mu_tpcBT_gut <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="gut")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="gut")])

numfams <- length(unique(creds_lac_t$species))
mus <- sum_lac_t$mean[1:5+(numfams*6)]
xs <- seq(mus[4], mus[5], length.out=100)
mu_tpc <- performance_mu(xs = xs, shape1 = mus[1], shape2 = mus[2], stretch = mus[3], x_min = mus[4], x_max = mus[5])
mu_tpcBT_lac <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="lac")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="lac")])

numfams <- length(unique(creds_nor_t$species))
mus <- sum_nor_t$mean[1:5+(numfams*6)]
xs <- seq(mus[4], mus[5], length.out=100)
mu_tpc <- performance_mu(xs = xs, shape1 = mus[1], shape2 = mus[2], stretch = mus[3], x_min = mus[4], x_max = mus[5])
mu_tpcBT_nor <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="nor")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="nor")])

numfams <- length(unique(creds_par_t$species))
mus <- sum_par_t$mean[1:5+(numfams*6)]
xs <- seq(mus[4], mus[5], length.out=100)
mu_tpc <- performance_mu(xs = xs, shape1 = mus[1], shape2 = mus[2], stretch = mus[3], x_min = mus[4], x_max = mus[5])
mu_tpcBT_par <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="par")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="par")])

numfams <- length(unique(creds_ver_t$species))
mus <- sum_ver_t$mean[1:5+(numfams*6)]
xs <- seq(mus[4], mus[5], length.out=100)
mu_tpc <- performance_mu(xs = xs, shape1 = mus[1], shape2 = mus[2], stretch = mus[3], x_min = mus[4], x_max = mus[5])
mu_tpcBT_ver <- data.frame(x= xs+meansTot$dayTemp[which(meansTot$species=="ver")], y=mu_tpc*meansTot$RGR[which(meansTot$species=="ver")])



############ 
# load in mbreadth and tbreadth limits calculated for each species
tpc_species4 <- read.csv("Analysis-output/Species/species-limits.csv")



#######################
# plot temperature curves

clwd <- 0.5
plwd <- 3
plwd2 <- 2
calp <- 0.5
xlims <- c(-2,62)

ax1 <- 10
ax2 <- 10
tlen <- 0.25
bys <- 0
len <- 0.25
ang <- 90
lsz <- 2
lalph <- 0.8

#######
# floribundus-norrisii

creds_flo_t$speciesF <- as.factor(creds_flo_t$species)
creds_nor_t$speciesF <- as.factor(creds_nor_t$species)
t_flonor <- ggplot() +
  labs(title="", x="", 
       y=expression(paste("RGR (number number"^-1, "day"^-1, ")"))) +
  geom_line(data = creds_flo_t, aes(x, mu, group=speciesF), color="royalblue3", inherit.aes = F, lwd = clwd) +
  geom_line(data = creds_nor_t, aes(x, mu, group=speciesF), color="lightskyblue1", inherit.aes = F, lwd = clwd) + 
  geom_line(data=mu_tpcBT_flo, aes(x, y), inherit.aes = F, lwd = plwd, color="black") +
  geom_line(data=mu_tpcBT_flo, aes(x, y), inherit.aes = F, lwd = plwd2, color="royalblue3") +
  geom_line(data=mu_tpcBT_nor, aes(x, y), inherit.aes = F, lwd = plwd, color="black") +
  geom_line(data=mu_tpcBT_nor, aes(x, y), inherit.aes = F, lwd = plwd2, color="lightskyblue1") +
  theme_bw() + 
  theme(axis.text=element_text(size=ax2), axis.title=element_text(size=ax1), 
        strip.background=element_rect(colour=NA,fill=NA), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border = element_rect(size=1), 
        plot.margin = margin(0,0.2,0,0, "cm"),
        axis.ticks.length=unit(tlen, "cm")) +
  scale_y_continuous(limits=c(0,0.235),expand=c(0,0)) + 
  scale_x_continuous(limits=xlims, expand=c(0,0))+
  annotate("segment", x=tpc_species4$t.min[which(tpc_species4$species=="flo")], 
           y=bys, xend=tpc_species4$t.max[which(tpc_species4$species=="flo")], yend=bys,
           col="royalblue3", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +
  annotate("segment", x=tpc_species4$t.min[which(tpc_species4$species=="nor")], 
           y=bys, xend=tpc_species4$t.max[which(tpc_species4$species=="nor")], yend=bys,
           col="lightskyblue1", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +  
  annotate("segment", x=tpc_species4$t.max[which(tpc_species4$species=="flo")], 
           y=bys, xend=tpc_species4$t.min[which(tpc_species4$species=="flo")], yend=bys,
           col="royalblue3", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +
  annotate("segment", x=tpc_species4$t.max[which(tpc_species4$species=="nor")], 
           y=bys, xend=tpc_species4$t.min[which(tpc_species4$species=="nor")], yend=bys,
           col="lightskyblue1", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +  
  coord_cartesian(clip="off")



#######
# cardinalis-parishii

creds_car_t$speciesF <- as.factor(creds_car_t$species)
creds_par_t$speciesF <- as.factor(creds_par_t$species)
t_carpar <- ggplot() +
  labs(title="", x="", 
       y=expression(paste("RGR (number number"^-1, "day"^-1, ")"))) +
  geom_line(data = creds_par_t, aes(x, mu, group=speciesF), color="thistle", inherit.aes = F, lwd = clwd) + 
  geom_line(data = creds_car_t, aes(x, mu, group=speciesF), color="mediumpurple3", inherit.aes = F, lwd = clwd) +
  geom_line(data=mu_tpcBT_par, aes(x, y), inherit.aes = F, lwd = plwd, color="black") +
  geom_line(data=mu_tpcBT_par, aes(x, y), inherit.aes = F, lwd = plwd2, color="thistle") +
  geom_line(data=mu_tpcBT_car, aes(x, y), inherit.aes = F, lwd = plwd, color="black") +
  geom_line(data=mu_tpcBT_car, aes(x, y), inherit.aes = F, lwd = plwd2, color="mediumpurple3") +
  theme_bw() + 
  theme(axis.text=element_text(size=ax2), axis.title=element_text(size=ax1), 
        strip.background=element_rect(colour=NA,fill=NA), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border = element_rect(size=1), 
        plot.margin = margin(0,0.2,0,0, "cm"),
        axis.ticks.length=unit(tlen, "cm")) +
  scale_y_continuous(limits=c(0,0.15),expand=c(0,0)) + 
  scale_x_continuous(limits=xlims, expand=c(0,0))+
  annotate("segment", x=tpc_species4$t.min[which(tpc_species4$species=="car")], 
           y=bys, xend=tpc_species4$t.max[which(tpc_species4$species=="car")], yend=bys,
           col="mediumpurple3", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +
  annotate("segment", x=tpc_species4$t.min[which(tpc_species4$species=="par")], 
           y=bys, xend=tpc_species4$t.max[which(tpc_species4$species=="par")], yend=bys,
           col="thistle", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +  
  annotate("segment", x=tpc_species4$t.max[which(tpc_species4$species=="car")], 
           y=bys, xend=tpc_species4$t.min[which(tpc_species4$species=="car")], yend=bys,
           col="mediumpurple3", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +
  annotate("segment", x=tpc_species4$t.max[which(tpc_species4$species=="par")], 
           y=bys, xend=tpc_species4$t.min[which(tpc_species4$species=="par")], yend=bys,
           col="thistle", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +  
  coord_cartesian(clip="off")



#######
# verbenaceus-eastwoodiae

creds_ver_t$speciesF <- as.factor(creds_ver_t$species)
creds_eas_t$speciesF <- as.factor(creds_eas_t$species)
t_vereas <- ggplot() +
  labs(title="", x="", 
       y=expression(paste("RGR (number number"^-1, "day"^-1, ")"))) +
  geom_line(data = creds_eas_t, aes(x, mu, group=speciesF), color="yellow", inherit.aes = F, lwd = clwd) + 
  geom_line(data = creds_ver_t, aes(x, mu, group=speciesF), color="goldenrod2", inherit.aes = F, lwd = clwd) +
  geom_line(data=mu_tpcBT_eas, aes(x, y), inherit.aes = F, lwd = plwd, color="black") +
  geom_line(data=mu_tpcBT_eas, aes(x, y), inherit.aes = F, lwd = plwd2, color="yellow") +
  geom_line(data=mu_tpcBT_ver, aes(x, y), inherit.aes = F, lwd = plwd, color="black") +
  geom_line(data=mu_tpcBT_ver, aes(x, y), inherit.aes = F, lwd = plwd2, color="goldenrod2") +
  theme_bw() + 
  theme(axis.text=element_text(size=ax2), axis.title=element_text(size=ax1), 
        strip.background=element_rect(colour=NA,fill=NA), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border = element_rect(size=1), 
        plot.margin = margin(0,0.2,0,0, "cm"),
        axis.ticks.length=unit(tlen, "cm")) +
  scale_y_continuous(limits=c(0,0.1225),expand=c(0,0)) + 
  scale_x_continuous(limits=xlims, expand=c(0,0))+
  annotate("segment", x=tpc_species4$t.min[which(tpc_species4$species=="ver")], 
           y=bys, xend=tpc_species4$t.max[which(tpc_species4$species=="ver")], yend=bys,
           col="goldenrod3", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +
  annotate("segment", x=tpc_species4$t.min[which(tpc_species4$species=="eas")], 
           y=bys, xend=tpc_species4$t.max[which(tpc_species4$species=="eas")], yend=bys,
           col="lightgoldenrod1", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +  
  annotate("segment", x=tpc_species4$t.max[which(tpc_species4$species=="ver")], 
           y=bys, xend=tpc_species4$t.min[which(tpc_species4$species=="ver")], yend=bys,
           col="goldenrod3", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +
  annotate("segment", x=tpc_species4$t.max[which(tpc_species4$species=="eas")], 
           y=bys, xend=tpc_species4$t.min[which(tpc_species4$species=="eas")], yend=bys,
           col="lightgoldenrod1", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +  
  coord_cartesian(clip="off")




#######
# bicolor-filicaulis

creds_bic_t$speciesF <- as.factor(creds_bic_t$species)
creds_fil_t$speciesF <- as.factor(creds_fil_t$species)
t_bicfil <- ggplot() +
  labs(title="", x="", 
       y=expression(paste("RGR (cm cm"^-1, "day"^-1, ")"))) +
  geom_line(data = creds_fil_t, aes(x, mu, group=speciesF), color="palegreen1", inherit.aes = F, lwd = clwd) + 
  geom_line(data = creds_bic_t, aes(x, mu, group=speciesF), color="green4", inherit.aes = F, lwd = clwd) +
  geom_line(data=mu_tpcBT_fil, aes(x, y), inherit.aes = F, lwd = plwd, color="black") +
  geom_line(data=mu_tpcBT_fil, aes(x, y), inherit.aes = F, lwd = plwd2, color="palegreen1") +
  geom_line(data=mu_tpcBT_bic, aes(x, y), inherit.aes = F, lwd = plwd, color="black") +
  geom_line(data=mu_tpcBT_bic, aes(x, y), inherit.aes = F, lwd = plwd2, color="green4") +
  theme_bw() + 
  theme(axis.text=element_text(size=ax2), axis.title=element_text(size=ax1), 
        strip.background=element_rect(colour=NA,fill=NA), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border = element_rect(size=1), 
        plot.margin = margin(0,0.2,0,0, "cm"),
        axis.ticks.length=unit(tlen, "cm")) +
  scale_y_continuous(limits=c(0,0.45),expand=c(0,0)) + 
  scale_x_continuous(limits=xlims, expand=c(0,0))+
  annotate("segment", x=tpc_species4$t.min[which(tpc_species4$species=="bic")], 
           y=bys, xend=tpc_species4$t.max[which(tpc_species4$species=="bic")], yend=bys,
           col="green4", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +
  annotate("segment", x=tpc_species4$t.min[which(tpc_species4$species=="fil")], 
           y=bys, xend=tpc_species4$t.max[which(tpc_species4$species=="fil")], yend=bys,
           col="palegreen1", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +  
  annotate("segment", x=tpc_species4$t.max[which(tpc_species4$species=="bic")], 
           y=bys, xend=tpc_species4$t.min[which(tpc_species4$species=="bic")], yend=bys,
           col="green4", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +
  annotate("segment", x=tpc_species4$t.max[which(tpc_species4$species=="fil")], 
           y=bys, xend=tpc_species4$t.min[which(tpc_species4$species=="fil")], yend=bys,
           col="palegreen1", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +  
  coord_cartesian(clip="off")


#######
# guttatus-laciniatus

creds_gut_t$speciesF <- as.factor(creds_gut_t$species)
creds_lac_t$speciesF <- as.factor(creds_lac_t$species)
t_gutlac <- ggplot() +
  labs(title="", x="", 
       y=expression(paste("RGR (cm cm"^-1, "day"^-1, ")"))) +
  geom_line(data = creds_lac_t, aes(x, mu, group=speciesF), color="violetred4", inherit.aes = F, lwd = clwd) + 
  geom_line(data = creds_gut_t, aes(x, mu, group=speciesF), color="palevioletred1", inherit.aes = F, lwd = clwd) +
  geom_line(data=mu_tpcBT_lac, aes(x, y), inherit.aes = F, lwd = plwd, color="black") +
  geom_line(data=mu_tpcBT_lac, aes(x, y), inherit.aes = F, lwd = plwd2, color="violetred4") +
  geom_line(data=mu_tpcBT_gut, aes(x, y), inherit.aes = F, lwd = plwd, color="black") +
  geom_line(data=mu_tpcBT_gut, aes(x, y), inherit.aes = F, lwd = plwd2, color="palevioletred1") +
  theme_bw() + 
  theme(axis.text=element_text(size=ax2), axis.title=element_text(size=ax1), 
        strip.background=element_rect(colour=NA,fill=NA), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.border = element_rect(size=1), 
        plot.margin = margin(0,0.2,0,0, "cm"),
        axis.ticks.length=unit(tlen, "cm")) +
  scale_y_continuous(limits=c(0,1.25),expand=c(0,0)) + 
  scale_x_continuous(limits=xlims, expand=c(0,0))+
  annotate("segment", x=tpc_species4$t.min[which(tpc_species4$species=="gut")], 
           y=bys, xend=tpc_species4$t.max[which(tpc_species4$species=="gut")], yend=bys,
           col="palevioletred1", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +
  annotate("segment", x=tpc_species4$t.min[which(tpc_species4$species=="lac")], 
           y=bys, xend=tpc_species4$t.max[which(tpc_species4$species=="lac")], yend=bys,
           col="violetred3", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +  
  annotate("segment", x=tpc_species4$t.max[which(tpc_species4$species=="gut")], 
           y=bys, xend=tpc_species4$t.min[which(tpc_species4$species=="gut")], yend=bys,
           col="palevioletred1", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +
  annotate("segment", x=tpc_species4$t.max[which(tpc_species4$species=="lac")], 
           y=bys, xend=tpc_species4$t.min[which(tpc_species4$species=="lac")], yend=bys,
           col="violetred3", arrow=arrow(length=unit(len, "cm"), angle=ang), size=lsz, alpha=lalph) +  
  coord_cartesian(clip="off")




#######
# make labels for species pairs
sizes <- 5
  
text = expression(paste(italic("M. cardinalis"), " & "))
car <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text, angle=90) + 
  theme_void()

text = expression(paste(italic("M. verbenaceus"), " & "))
ver <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text, angle=90) + 
  theme_void()

text = expression(paste(italic("M. floribundus"), " & "))
flo <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text, angle=90) + 
  theme_void()

text = expression(paste(italic("M. bicolor"), " & "))
bic <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text, angle=90) + 
  theme_void()

text = expression(paste(italic("M. guttatus"), " & "))
gut <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text, angle=90) + 
  theme_void()

text = expression(paste(italic("M. parishii")))
par <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text, angle=90) + 
  theme_void()

text = expression(paste(italic("M. eastwoodiae")))
eas <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text, angle=90) + 
  theme_void()

text = expression(paste(italic("M. norrisii")))
nor <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text, angle=90) + 
  theme_void()

text = expression(paste(italic("M. filicaulis")))
fil <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text, angle=90) + 
  theme_void()

text = expression(paste(italic(" M. laciniatus")))
lac <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text, angle=90) + 
  theme_void()


#######
# make labels for x axes
sizes <- 6

text = ""
notext <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text) + 
  theme_void()

text = expression(paste("Temperature (°C)"))
temp <- ggplot() + 
  annotate("text", x = 4, y = 25, size=sizes, label = text) + 
  theme_void() +
  scale_x_continuous(limits=c(2,6), expand = c(0, 0), breaks = NULL) +
  scale_y_continuous(limits=c(20,30), expand = c(0, 0), breaks = NULL) 
  




#######
# create legend
#######

splist <- c(expression(paste(italic("M. parishii"))),
            expression(paste(italic("M. cardinalis"))),
            expression(paste(italic("M. verbenaceus"))),
            expression(paste(italic("M. eastwoodiae"))),
            expression(paste(italic("M. floribundus"))),
            expression(paste(italic("M. norrisii"))),
            expression(paste(italic("M. bicolor"))),
            expression(paste(italic("M. fillicaulis"))),
            expression(paste(italic("M. laciniatus"))),
            expression(paste(italic("M. guttatus"))))
spcol <- c("mediumpurple3",
           "thistle",
           "goldenrod2",
           "yellow",
           "royalblue3",
           "lightskyblue1",
           "green4",
           "palegreen1",
           "violetred4",
           "palevioletred1")

leg <- ggplot() +
  annotate(geom="text", x=0.55, y=c(10:1), color="black", label=c(splist), size=6, hjust = 0) + 
  annotate(geom="segment", x=0.3, y=c(10:1), xend=0.45, yend=c(10:1), color=spcol, size=2, lineend="round", linejoin="mitre") + 
  scale_x_continuous(limits=c(0,1.5), expand = c(0, 0), breaks = NULL) +
  scale_y_continuous(limits=c(0,11), expand = c(0, 0), breaks = NULL) +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA), text=element_text(size=0), 
        strip.background=element_rect(colour=NA,fill=NA), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())


# add panel labels
library(grid)
library(gridExtra)
t_carpar2 <- arrangeGrob(t_carpar, top = textGrob("A", x = unit(0.125, "npc"), y = unit(0.1, "npc"), just=c("left","top"), gp=gpar(fontsize=16)))
t_vereas2 <- arrangeGrob(t_vereas, top = textGrob("B", x = unit(0.125, "npc"), y = unit(0.1, "npc"), just=c("left","top"), gp=gpar(fontsize=16)))
t_flonor2 <- arrangeGrob(t_flonor, top = textGrob("C", x = unit(0.125, "npc"), y = unit(0.1, "npc"), just=c("left","top"), gp=gpar(fontsize=16)))
t_bicfil2 <- arrangeGrob(t_bicfil, top = textGrob("D", x = unit(0.125, "npc"), y = unit(0.1, "npc"), just=c("left","top"), gp=gpar(fontsize=16)))
t_gutlac2 <- arrangeGrob(t_gutlac, top = textGrob("E", x = unit(0.125, "npc"), y = unit(0.1, "npc"), just=c("left","top"), gp=gpar(fontsize=16)))


# now export the plot
# note that this is script to make Fig. 2 in manuscript as of 20210728
pdf("Figures/2_Figure-2_family-level-curves.pdf", height=12, width=10)
figure <- ggarrange(t_carpar2, t_bicfil2,
                    notext, notext, 
                    t_vereas2, t_gutlac2,
                    notext, temp, 
                    t_flonor2, leg,
                    temp, notext, 
                    ncol=2, nrow=6, 
                    widths=c(1,1), 
                    heights=c(1, 0.06, 1, 0.06, 1, 0.06))
figure
dev.off()



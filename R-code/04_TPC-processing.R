#### PROJECT: Performance curves of Mimulus species
#### PURPOSE: Process temperature performance models of 
#### families within each species
#### AUTHOR/DATE: RCW/2020-02-10

###### Warnings:
###### 1) This script requires a lot of computational power.
###### 2) The model output files required to run this code (lines 36-45) and tidy_perf files exported from this code (e.g., line 166) are not located in this git repo because they are too large. Please contact Seema Sheth (ssheth3@ncsu.edu) if you are interested in obtaining these files.



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
# read in model output
############
model_fits_bic <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_bic.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_bic.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_bic.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_bic.samples_4.csv"))
model_fits_car <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_car.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_car.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_car.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_car.samples_4.csv"))
model_fits_eas <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_eas.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_eas.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_eas.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_eas.samples_4.csv"))
model_fits_fil <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_fil.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_fil.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_fil.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_fil.samples_4.csv"))
model_fits_flo <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_flo.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_flo.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_flo.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_flo.samples_4.csv"))
model_fits_gut <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_gut.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_gut.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_gut.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_gut.samples_4.csv"))
model_fits_lac <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_lac.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_lac.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_lac.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_lac.samples_4.csv"))
model_fits_nor <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_nor.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_nor.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_nor.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_nor.samples_4.csv"))
model_fits_par <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_par.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_par.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_par.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_par.samples_4.csv"))
model_fits_ver <- rstan::read_stan_csv(c("Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_ver.samples_1.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_ver.samples_2.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_ver.samples_3.csv", "Analysis-output/TPC/Models-stand-x-priors_all-fams/stan_ver.samples_4.csv"))







# settings and functions
color_scheme_set("blue")
bayes_p <- function(stan_df, raw_df, raw_group, raw_treatment, raw_response, ndraw = nrow(stan_df)){
  
  1:ndraw %>% 
    map_df(~{
      draw_i <- stan_df[.x, ]
      tidy_spp <- stan_df$species[.x]
      df_i <- raw_df[raw_df[[raw_group]] == tidy_spp,]
      
      
      mus <- performance_mu(xs = df_i[[raw_treatment]], 
                            shape1 = draw_i$shape1, 
                            shape2 = draw_i$shape2, 
                            stretch = draw_i$stretch, 
                            x_min = draw_i$x_min,
                            x_max = draw_i$x_max
      )
      
      pseudo <- posterior_predict(x = df_i[[raw_treatment]], draw_i)  
      
      data.frame(
        ssq_obs = sum((mus - df_i[[raw_response]])^2),
        ssq_pseudo = sum((mus - pseudo$trait)^2),
        species = tidy_spp
      )
      
    })
}
optimum_breadth <- 
  function(par_df, prop_max = 0.5, n_grid = 100){
    1:nrow(par_df) %>% map_df(function(i){
      xs = seq(par_df$x_min[i],
               par_df$x_max[i],
               length.out = n_grid)
      
      sim_mu <- performr::performance_mu(
        xs = xs, 
        shape1 = par_df$shape1[i],
        shape2 = par_df$shape2[i],
        stretch = par_df$stretch[i],
        x_min = par_df$x_min[i],
        x_max = par_df$x_max[i]
      )
      
      max_y <- which.max(sim_mu)
      prop_max_y <- sim_mu[max_y] * prop_max
      idx_max_low <- which.min((sim_mu[1:(max_y-1)] - prop_max_y)^2)
      idx_max_high <- which.min((sim_mu[max_y:length(sim_mu)] - prop_max_y)^2) + max_y
      tibble("opt_breadth_low" = xs[idx_max_low],
             "opt_breadth_high" = xs[idx_max_high], 
             "opt_breadth" = xs[idx_max_high] - xs[idx_max_low])
    })
  }

predict_interval <- function (x, spp, par_df, x_draws, p){
  if (missing(x)) {
    x <- seq(min(par_df$x_min), max(par_df$x_max), length.out = 100)
  }
  sub_df <- par_df %>% filter(draw %in% x_draws)
  
  p %>% map_df(~{
    posterior_quantile(x = x, par_df = sub_df, p = .x) %>%
      group_by(species, x) %>%
      summarise_all(.funs = mean) %>%
      mutate(level = .x) %>%
      arrange(x) %>%
      dplyr::select(-draw) %>%
      mutate(level = round(.x * 100, 0))
  }) %>%
    arrange(species, x)
}


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


############ 
# load model, check diagnostics, visualize sampling, generate tidy dataframe, 
# calculate pmax, bayesian p value, and B50, generate TPC prediction,
# back-transform posterior draws, calculate parameter averages and credible intervals
# be sure to add code for each species from here on


tidy_perf_bic <- performr::perform_df(model_fits_bic, species_order = c(1:length(unique(bic_dat$familyI))))
tidy_perf_bic %<>% dplyr::ungroup() %>% mutate(idx = 1:n()) %>% group_by(idx) %>% do({ mutate(.data = ., max_RGR = performance_mu(xs = .$maxima, .$shape1, .$shape2, .$stretch, .$x_min, .$x_max))})
#ssq_df_bic <- bayes_p(stan_df = tidy_perf_bic, raw_df = bic_dat, raw_treatment = "dayTempc", raw_response = "RGRs", raw_group = "familyI")
tidy_perf_bic %<>% ungroup() %>% bind_cols(optimum_breadth(par_df = tidy_perf_bic,  prop_max = 0.5,  n_grid = 100))
tidy_perf_bic$B50_low <- tidy_perf_bic$opt_breadth_low + meansTot$dayTemp[which(meansTot$species=="bic")] 
tidy_perf_bic$B50_high <- tidy_perf_bic$opt_breadth_high + meansTot$dayTemp[which(meansTot$species=="bic")]  
tidy_perf_bic$B50 <- tidy_perf_bic$B50_high - tidy_perf_bic$B50_low 
draws_bic <- rstan::extract(model_fits_bic); ndraws_bic <- length(draws_bic$lp__)
x_seq_bic = seq(min(draws_bic$x_min),max(draws_bic$x_max),length.out = 100)
poly_draws_bic <- sample(1:ndraws_bic, 100)
creds_bic <-  predict_interval(x = x_seq_bic, spp = species, par_df = tidy_perf_bic, x_draws = poly_draws_bic, p = c(0.95, 0.5))
tidy_perf_bic$x_minBT <- tidy_perf_bic$x_min + meansTot$dayTemp[which(meansTot$species=="bic")] 
tidy_perf_bic$x_maxBT <- tidy_perf_bic$x_max + meansTot$dayTemp[which(meansTot$species=="bic")] 
tidy_perf_bic$maximaBT <- tidy_perf_bic$maxima + meansTot$dayTemp[which(meansTot$species=="bic")] 
tidy_perf_bic$breadthBT <- tidy_perf_bic$x_maxBT-tidy_perf_bic$x_minBT
tidy_perf_bic$max_RGRBT <- tidy_perf_bic$max_RGR * meansTot$RGR[which(meansTot$species=="bic")] 
write.csv(tidy_perf_bic, "Analysis-output/TPC/Models/tidy_perf_sxp_af_bic.csv")
creds_bic$x <- creds_bic$x + meansTot$dayTemp[which(meansTot$species=="bic")] ; creds_bic$mu <- creds_bic$mu * meansTot$RGR[which(meansTot$species=="bic")] 
creds_bic$upper <- creds_bic$upper * meansTot$RGR[which(meansTot$species=="bic")] ; creds_bic$lower <- creds_bic$lower * meansTot$RGR[which(meansTot$species=="bic")] 
write.csv(creds_bic, "Analysis-output/TPC/Creds/creds_sxp_af_bic.csv")
mean_df_bic <- tidy_perf_bic %>% group_by(species) %>% summarise_if(is.numeric, .funs = c(mean))
mean_df_bic$speciesn <- as.numeric(mean_df_bic$species); mean_df_bic <- mean_df_bic[order(mean_df_bic$speciesn),]
write.csv(mean_df_bic, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_bic.csv")
mean_df_ci_bic <- tidy_perf_bic %>%  group_by(species) %>% summarise(maximaBT_lci = quantile(maximaBT, probs=c(0.025)), maximaBT_uci = quantile(maximaBT, probs=c(0.975)), 
                                                                     B50_lci = quantile(B50, probs=c(0.025)), B50_uci = quantile(B50, probs=c(0.975)), 
                                                                     B50_low_lci = quantile(B50_low, probs=c(0.025)), B50_low_uci = quantile(B50_low, probs=c(0.975)), 
                                                                     B50_high_lci = quantile(B50_high, probs=c(0.025)), B50_high_uci = quantile(B50_high, probs=c(0.975)), 
                                                                     breadthBT_lci = quantile(breadthBT, probs=c(0.025)), breadthBT_uci = quantile(breadthBT, probs=c(0.975)), 
                                                                     x_minBT_lci = quantile(x_minBT, probs=c(0.025)), x_minBT_uci = quantile(x_minBT, probs=c(0.975)), 
                                                                     x_maxBT_lci = quantile(x_maxBT, probs=c(0.025)), x_maxBT_uci = quantile(x_maxBT, probs=c(0.975)), 
                                                                     max_RGR_lci = quantile(max_RGRBT, probs=c(0.025)), max_RGR_uci = quantile(max_RGRBT, probs=c(0.975)), 
                                                                     area_lci = quantile(area, probs=c(0.025)), area_uci = quantile(area, probs=c(0.975)))
mean_df_ci_bic$speciesn <- as.numeric(mean_df_ci_bic$species); mean_df_ci_bic <- mean_df_ci_bic[order(mean_df_ci_bic$speciesn),]
write.csv(mean_df_ci_bic, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_bic.csv")
mean_ci_table_bic <- mean_df_bic %>% right_join(mean_df_ci_bic, by="species")
mean_ci_table_bic <- as.data.frame(mean_ci_table_bic)
mean_ci_table_bic <- mean_ci_table_bic[,c("species", "speciesn.x", 
                                          "maximaBT", "maximaBT_lci", "maximaBT_uci", 
                                          "B50", "B50_lci", "B50_uci", 
                                          "B50_low", "B50_low_lci", "B50_low_uci", 
                                          "B50_high", "B50_high_lci", "B50_high_uci", 
                                          "breadthBT", "breadthBT_lci", "breadthBT_uci", 
                                          "x_minBT", "x_minBT_lci", "x_minBT_uci", 
                                          "x_maxBT", "x_maxBT_lci", "x_maxBT_uci", 
                                          "max_RGRBT", "max_RGR_lci", "max_RGR_uci", 
                                          "area", "area_lci", "area_uci")]
mean_ci_table_bic <- mean_ci_table_bic %>% mutate_at(vars(-c("species","speciesn.x")), funs(round(., 2)))
mean_ci_table_bic$ToptCI <- paste(mean_ci_table_bic$maximaBT, " [", mean_ci_table_bic$maximaBT_lci, ", ", mean_ci_table_bic$maximaBT_uci,"]", sep="")
mean_ci_table_bic$B50CI <- paste(mean_ci_table_bic$B50, " [", mean_ci_table_bic$B50_lci, ", ", mean_ci_table_bic$B50_uci,"]", sep="")
mean_ci_table_bic$B50_lowCI <- paste(mean_ci_table_bic$B50_low, " [", mean_ci_table_bic$B50_low_lci, ", ", mean_ci_table_bic$B50_low_uci,"]", sep="")
mean_ci_table_bic$B50_highCI <- paste(mean_ci_table_bic$B50_high, " [", mean_ci_table_bic$B50_high_lci, ", ", mean_ci_table_bic$B50_high_uci,"]", sep="")
mean_ci_table_bic$breadthCI <- paste(mean_ci_table_bic$breadthBT, " [", mean_ci_table_bic$breadthBT_lci, ", ", mean_ci_table_bic$breadthBT_uci,"]", sep="")
mean_ci_table_bic$x_minCI <- paste(mean_ci_table_bic$x_minBT, " [", mean_ci_table_bic$x_minBT_lci, ", ", mean_ci_table_bic$x_minBT_uci,"]", sep="")
mean_ci_table_bic$x_maxCI <- paste(mean_ci_table_bic$x_maxBT, " [", mean_ci_table_bic$x_maxBT_lci, ", ", mean_ci_table_bic$x_maxBT_uci,"]", sep="")
mean_ci_table_bic$max_RGRCI <- paste(mean_ci_table_bic$max_RGRBT, " [", mean_ci_table_bic$max_RGR_lci, ", ", mean_ci_table_bic$max_RGR_uci,"]", sep="")
mean_ci_table_bic$areaCI <- paste(mean_ci_table_bic$area, " [", mean_ci_table_bic$area_lci, ", ", mean_ci_table_bic$area_uci,"]", sep="")
mean_ci_table_bic <- mean_ci_table_bic[,c("species", "ToptCI", "B50CI", "B50_lowCI", "B50_highCI", "breadthCI", "x_minCI", "x_maxCI", "max_RGRCI", "areaCI")]
write.csv(mean_ci_table_bic, "Analysis-output/TPC/Mean-dfs/mean_ci_sxp_af_table_bic.csv")


tidy_perf_car <- performr::perform_df(model_fits_car, species_order = c(1:length(unique(car_dat$familyI)))) 
tidy_perf_car %<>% ungroup() %>% mutate(idx = 1:n()) %>% group_by(idx) %>% do({ mutate(.data = ., max_RGR = performance_mu(xs = .$maxima, .$shape1, .$shape2, .$stretch, .$x_min, .$x_max))})
#ssq_df_car <- bayes_p(stan_df = tidy_perf_car, raw_df = car_dat, raw_treatment = "dayTempc", raw_response = "RGRs", raw_group = "familyI")
tidy_perf_car %<>% ungroup() %>% bind_cols(optimum_breadth(par_df = tidy_perf_car,  prop_max = 0.5,  n_grid = 100))
tidy_perf_car$B50_low <- tidy_perf_car$opt_breadth_low + meansTot$dayTemp[which(meansTot$species=="car")] 
tidy_perf_car$B50_high <- tidy_perf_car$opt_breadth_high + meansTot$dayTemp[which(meansTot$species=="car")]  
tidy_perf_car$B50 <- tidy_perf_car$B50_high - tidy_perf_car$B50_low 
draws_car <- rstan::extract(model_fits_car); ndraws_car <- length(draws_car$lp__)
x_seq_car = seq(min(draws_car$x_min),max(draws_car$x_max),length.out = 100)
poly_draws_car <- sample(1:ndraws_car, 100)
creds_car <-  predict_interval(x = x_seq_car, spp = species, par_df = tidy_perf_car, x_draws = poly_draws_car, p = c(0.95, 0.5))
tidy_perf_car$x_minBT <- tidy_perf_car$x_min + meansTot$dayTemp[which(meansTot$species=="car")] 
tidy_perf_car$x_maxBT <- tidy_perf_car$x_max + meansTot$dayTemp[which(meansTot$species=="car")] 
tidy_perf_car$maximaBT <- tidy_perf_car$maxima + meansTot$dayTemp[which(meansTot$species=="car")] 
tidy_perf_car$breadthBT <- tidy_perf_car$x_maxBT-tidy_perf_car$x_minBT
tidy_perf_car$max_RGRBT <- tidy_perf_car$max_RGR * meansTot$RGR[which(meansTot$species=="car")] 
write.csv(tidy_perf_car, "Analysis-output/TPC/Models/tidy_perf_sxp_af_car.csv")
creds_car$x <- creds_car$x + meansTot$dayTemp[which(meansTot$species=="car")] ; creds_car$mu <- creds_car$mu * meansTot$RGR[which(meansTot$species=="car")] 
creds_car$upper <- creds_car$upper * meansTot$RGR[which(meansTot$species=="car")] ; creds_car$lower <- creds_car$lower * meansTot$RGR[which(meansTot$species=="car")] 
write.csv(creds_car, "Analysis-output/TPC/Creds/creds_sxp_af_car.csv")
mean_df_car <- tidy_perf_car %>% group_by(species) %>% summarise_if(is.numeric, .funs = c(mean))
mean_df_car$speciesn <- as.numeric(mean_df_car$species); mean_df_car <- mean_df_car[order(mean_df_car$speciesn),]
write.csv(mean_df_car, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_car.csv")
mean_df_ci_car <- tidy_perf_car %>%  group_by(species) %>% summarise(maximaBT_lci = quantile(maximaBT, probs=c(0.025)), maximaBT_uci = quantile(maximaBT, probs=c(0.975)), 
                                                                     B50_lci = quantile(B50, probs=c(0.025)), B50_uci = quantile(B50, probs=c(0.975)), 
                                                                     B50_low_lci = quantile(B50_low, probs=c(0.025)), B50_low_uci = quantile(B50_low, probs=c(0.975)), 
                                                                     B50_high_lci = quantile(B50_high, probs=c(0.025)), B50_high_uci = quantile(B50_high, probs=c(0.975)), 
                                                                     breadthBT_lci = quantile(breadthBT, probs=c(0.025)), breadthBT_uci = quantile(breadthBT, probs=c(0.975)), 
                                                                     x_minBT_lci = quantile(x_minBT, probs=c(0.025)), x_minBT_uci = quantile(x_minBT, probs=c(0.975)), 
                                                                     x_maxBT_lci = quantile(x_maxBT, probs=c(0.025)), x_maxBT_uci = quantile(x_maxBT, probs=c(0.975)), 
                                                                     max_RGR_lci = quantile(max_RGRBT, probs=c(0.025)), max_RGR_uci = quantile(max_RGRBT, probs=c(0.975)), 
                                                                     area_lci = quantile(area, probs=c(0.025)), area_uci = quantile(area, probs=c(0.975)))
mean_df_ci_car$speciesn <- as.numeric(mean_df_ci_car$species); mean_df_ci_car <- mean_df_ci_car[order(mean_df_ci_car$speciesn),]
write.csv(mean_df_ci_car, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_car.csv")
mean_ci_table_car <- mean_df_car %>% right_join(mean_df_ci_car, by="species")
mean_ci_table_car <- as.data.frame(mean_ci_table_car)
mean_ci_table_car <- mean_ci_table_car[,c("species", "speciesn.x", 
                                          "maximaBT", "maximaBT_lci", "maximaBT_uci", 
                                          "B50", "B50_lci", "B50_uci", 
                                          "B50_low", "B50_low_lci", "B50_low_uci", 
                                          "B50_high", "B50_high_lci", "B50_high_uci", 
                                          "breadthBT", "breadthBT_lci", "breadthBT_uci", 
                                          "x_minBT", "x_minBT_lci", "x_minBT_uci", 
                                          "x_maxBT", "x_maxBT_lci", "x_maxBT_uci", 
                                          "max_RGRBT", "max_RGR_lci", "max_RGR_uci", 
                                          "area", "area_lci", "area_uci")]
mean_ci_table_car <- mean_ci_table_car %>% mutate_at(vars(-c("species","speciesn.x")), funs(round(., 2)))
mean_ci_table_car$ToptCI <- paste(mean_ci_table_car$maximaBT, " [", mean_ci_table_car$maximaBT_lci, ", ", mean_ci_table_car$maximaBT_uci,"]", sep="")
mean_ci_table_car$B50CI <- paste(mean_ci_table_car$B50, " [", mean_ci_table_car$B50_lci, ", ", mean_ci_table_car$B50_uci,"]", sep="")
mean_ci_table_car$B50_lowCI <- paste(mean_ci_table_car$B50_low, " [", mean_ci_table_car$B50_low_lci, ", ", mean_ci_table_car$B50_low_uci,"]", sep="")
mean_ci_table_car$B50_highCI <- paste(mean_ci_table_car$B50_high, " [", mean_ci_table_car$B50_high_lci, ", ", mean_ci_table_car$B50_high_uci,"]", sep="")
mean_ci_table_car$breadthCI <- paste(mean_ci_table_car$breadthBT, " [", mean_ci_table_car$breadthBT_lci, ", ", mean_ci_table_car$breadthBT_uci,"]", sep="")
mean_ci_table_car$x_minCI <- paste(mean_ci_table_car$x_minBT, " [", mean_ci_table_car$x_minBT_lci, ", ", mean_ci_table_car$x_minBT_uci,"]", sep="")
mean_ci_table_car$x_maxCI <- paste(mean_ci_table_car$x_maxBT, " [", mean_ci_table_car$x_maxBT_lci, ", ", mean_ci_table_car$x_maxBT_uci,"]", sep="")
mean_ci_table_car$max_RGRCI <- paste(mean_ci_table_car$max_RGRBT, " [", mean_ci_table_car$max_RGR_lci, ", ", mean_ci_table_car$max_RGR_uci,"]", sep="")
mean_ci_table_car$areaCI <- paste(mean_ci_table_car$area, " [", mean_ci_table_car$area_lci, ", ", mean_ci_table_car$area_uci,"]", sep="")
mean_ci_table_car <- mean_ci_table_car[,c("species", "ToptCI", "B50CI", "B50_lowCI", "B50_highCI", "breadthCI", "x_minCI", "x_maxCI", "max_RGRCI", "areaCI")]
write.csv(mean_ci_table_car, "Analysis-output/TPC/Mean-dfs/mean_ci_sxp_af_table_car.csv")


tidy_perf_eas <- performr::perform_df(model_fits_eas, species_order = c(1:length(unique(eas_dat$familyI)))) 
tidy_perf_eas %<>% ungroup() %>% mutate(idx = 1:n()) %>% group_by(idx) %>% do({ mutate(.data = ., max_RGR = performance_mu(xs = .$maxima, .$shape1, .$shape2, .$stretch, .$x_min, .$x_max))})
#ssq_df_eas <- bayes_p(stan_df = tidy_perf_eas, raw_df = eas_dat, raw_treatment = "dayTempc", raw_response = "RGRs", raw_group = "familyI")
tidy_perf_eas %<>% ungroup() %>% bind_cols(optimum_breadth(par_df = tidy_perf_eas,  prop_max = 0.5,  n_grid = 100))
tidy_perf_eas$B50_low <- tidy_perf_eas$opt_breadth_low + meansTot$dayTemp[which(meansTot$species=="eas")] 
tidy_perf_eas$B50_high <- tidy_perf_eas$opt_breadth_high + meansTot$dayTemp[which(meansTot$species=="eas")]  
tidy_perf_eas$B50 <- tidy_perf_eas$B50_high - tidy_perf_eas$B50_low 
draws_eas <- rstan::extract(model_fits_eas); ndraws_eas <- length(draws_eas$lp__)
x_seq_eas = seq(min(draws_eas$x_min),max(draws_eas$x_max),length.out = 100)
poly_draws_eas <- sample(1:ndraws_eas, 100)
creds_eas <-  predict_interval(x = x_seq_eas, spp = species, par_df = tidy_perf_eas, x_draws = poly_draws_eas, p = c(0.95, 0.5))
tidy_perf_eas$x_minBT <- tidy_perf_eas$x_min + meansTot$dayTemp[which(meansTot$species=="eas")] 
tidy_perf_eas$x_maxBT <- tidy_perf_eas$x_max + meansTot$dayTemp[which(meansTot$species=="eas")] 
tidy_perf_eas$maximaBT <- tidy_perf_eas$maxima + meansTot$dayTemp[which(meansTot$species=="eas")] 
tidy_perf_eas$breadthBT <- tidy_perf_eas$x_maxBT-tidy_perf_eas$x_minBT
tidy_perf_eas$max_RGRBT <- tidy_perf_eas$max_RGR * meansTot$RGR[which(meansTot$species=="eas")] 
write.csv(tidy_perf_eas, "Analysis-output/TPC/Models/tidy_perf_sxp_af_eas.csv")
creds_eas$x <- creds_eas$x + meansTot$dayTemp[which(meansTot$species=="eas")] ; creds_eas$mu <- creds_eas$mu * meansTot$RGR[which(meansTot$species=="eas")] 
creds_eas$upper <- creds_eas$upper * meansTot$RGR[which(meansTot$species=="eas")] ; creds_eas$lower <- creds_eas$lower * meansTot$RGR[which(meansTot$species=="eas")] 
write.csv(creds_eas, "Analysis-output/TPC/Creds/creds_sxp_af_eas.csv")
mean_df_eas <- tidy_perf_eas %>% group_by(species) %>% summarise_if(is.numeric, .funs = c(mean))
mean_df_eas$speciesn <- as.numeric(mean_df_eas$species); mean_df_eas <- mean_df_eas[order(mean_df_eas$speciesn),]
write.csv(mean_df_eas, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_eas.csv")
mean_df_ci_eas <- tidy_perf_eas %>%  group_by(species) %>% summarise(maximaBT_lci = quantile(maximaBT, probs=c(0.025)), maximaBT_uci = quantile(maximaBT, probs=c(0.975)), 
                                                                     B50_lci = quantile(B50, probs=c(0.025)), B50_uci = quantile(B50, probs=c(0.975)), 
                                                                     B50_low_lci = quantile(B50_low, probs=c(0.025)), B50_low_uci = quantile(B50_low, probs=c(0.975)), 
                                                                     B50_high_lci = quantile(B50_high, probs=c(0.025)), B50_high_uci = quantile(B50_high, probs=c(0.975)), 
                                                                     breadthBT_lci = quantile(breadthBT, probs=c(0.025)), breadthBT_uci = quantile(breadthBT, probs=c(0.975)), 
                                                                     x_minBT_lci = quantile(x_minBT, probs=c(0.025)), x_minBT_uci = quantile(x_minBT, probs=c(0.975)), 
                                                                     x_maxBT_lci = quantile(x_maxBT, probs=c(0.025)), x_maxBT_uci = quantile(x_maxBT, probs=c(0.975)), 
                                                                     max_RGR_lci = quantile(max_RGRBT, probs=c(0.025)), max_RGR_uci = quantile(max_RGRBT, probs=c(0.975)), 
                                                                     area_lci = quantile(area, probs=c(0.025)), area_uci = quantile(area, probs=c(0.975)))
mean_df_ci_eas$speciesn <- as.numeric(mean_df_ci_eas$species); mean_df_ci_eas <- mean_df_ci_eas[order(mean_df_ci_eas$speciesn),]
write.csv(mean_df_ci_eas, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_eas.csv")
mean_ci_table_eas <- mean_df_eas %>% right_join(mean_df_ci_eas, by="species")
mean_ci_table_eas <- as.data.frame(mean_ci_table_eas)
mean_ci_table_eas <- mean_ci_table_eas[,c("species", "speciesn.x", 
                                          "maximaBT", "maximaBT_lci", "maximaBT_uci", 
                                          "B50", "B50_lci", "B50_uci", 
                                          "B50_low", "B50_low_lci", "B50_low_uci", 
                                          "B50_high", "B50_high_lci", "B50_high_uci", 
                                          "breadthBT", "breadthBT_lci", "breadthBT_uci", 
                                          "x_minBT", "x_minBT_lci", "x_minBT_uci", 
                                          "x_maxBT", "x_maxBT_lci", "x_maxBT_uci", 
                                          "max_RGRBT", "max_RGR_lci", "max_RGR_uci", 
                                          "area", "area_lci", "area_uci")]
mean_ci_table_eas <- mean_ci_table_eas %>% mutate_at(vars(-c("species","speciesn.x")), funs(round(., 2)))
mean_ci_table_eas$ToptCI <- paste(mean_ci_table_eas$maximaBT, " [", mean_ci_table_eas$maximaBT_lci, ", ", mean_ci_table_eas$maximaBT_uci,"]", sep="")
mean_ci_table_eas$B50CI <- paste(mean_ci_table_eas$B50, " [", mean_ci_table_eas$B50_lci, ", ", mean_ci_table_eas$B50_uci,"]", sep="")
mean_ci_table_eas$B50_lowCI <- paste(mean_ci_table_eas$B50_low, " [", mean_ci_table_eas$B50_low_lci, ", ", mean_ci_table_eas$B50_low_uci,"]", sep="")
mean_ci_table_eas$B50_highCI <- paste(mean_ci_table_eas$B50_high, " [", mean_ci_table_eas$B50_high_lci, ", ", mean_ci_table_eas$B50_high_uci,"]", sep="")
mean_ci_table_eas$breadthCI <- paste(mean_ci_table_eas$breadthBT, " [", mean_ci_table_eas$breadthBT_lci, ", ", mean_ci_table_eas$breadthBT_uci,"]", sep="")
mean_ci_table_eas$x_minCI <- paste(mean_ci_table_eas$x_minBT, " [", mean_ci_table_eas$x_minBT_lci, ", ", mean_ci_table_eas$x_minBT_uci,"]", sep="")
mean_ci_table_eas$x_maxCI <- paste(mean_ci_table_eas$x_maxBT, " [", mean_ci_table_eas$x_maxBT_lci, ", ", mean_ci_table_eas$x_maxBT_uci,"]", sep="")
mean_ci_table_eas$max_RGRCI <- paste(mean_ci_table_eas$max_RGRBT, " [", mean_ci_table_eas$max_RGR_lci, ", ", mean_ci_table_eas$max_RGR_uci,"]", sep="")
mean_ci_table_eas$areaCI <- paste(mean_ci_table_eas$area, " [", mean_ci_table_eas$area_lci, ", ", mean_ci_table_eas$area_uci,"]", sep="")
mean_ci_table_eas <- mean_ci_table_eas[,c("species", "ToptCI", "B50CI", "B50_lowCI", "B50_highCI", "breadthCI", "x_minCI", "x_maxCI", "max_RGRCI", "areaCI")]
write.csv(mean_ci_table_eas, "Analysis-output/TPC/Mean-dfs/mean_ci_sxp_af_table_eas.csv")



tidy_perf_fil <- performr::perform_df(model_fits_fil, species_order = c(1:length(unique(fil_dat$familyI)))) 
tidy_perf_fil %<>% ungroup() %>% mutate(idx = 1:n()) %>% group_by(idx) %>% do({ mutate(.data = ., max_RGR = performance_mu(xs = .$maxima, .$shape1, .$shape2, .$stretch, .$x_min, .$x_max))})
#ssq_df_fil <- bayes_p(stan_df = tidy_perf_fil, raw_df = fil_dat, raw_treatment = "dayTempc", raw_response = "RGRs", raw_group = "familyI")
tidy_perf_fil %<>% ungroup() %>% bind_cols(optimum_breadth(par_df = tidy_perf_fil,  prop_max = 0.5,  n_grid = 100))
tidy_perf_fil$B50_low <- tidy_perf_fil$opt_breadth_low + meansTot$dayTemp[which(meansTot$species=="fil")] 
tidy_perf_fil$B50_high <- tidy_perf_fil$opt_breadth_high + meansTot$dayTemp[which(meansTot$species=="fil")]  
tidy_perf_fil$B50 <- tidy_perf_fil$B50_high - tidy_perf_fil$B50_low 
draws_fil <- rstan::extract(model_fits_fil); ndraws_fil <- length(draws_fil$lp__)
x_seq_fil = seq(min(draws_fil$x_min),max(draws_fil$x_max),length.out = 100)
poly_draws_fil <- sample(1:ndraws_fil, 100)
creds_fil <-  predict_interval(x = x_seq_fil, spp = species, par_df = tidy_perf_fil, x_draws = poly_draws_fil, p = c(0.95, 0.5))
tidy_perf_fil$x_minBT <- tidy_perf_fil$x_min + meansTot$dayTemp[which(meansTot$species=="fil")] 
tidy_perf_fil$x_maxBT <- tidy_perf_fil$x_max + meansTot$dayTemp[which(meansTot$species=="fil")] 
tidy_perf_fil$maximaBT <- tidy_perf_fil$maxima + meansTot$dayTemp[which(meansTot$species=="fil")] 
tidy_perf_fil$breadthBT <- tidy_perf_fil$x_maxBT-tidy_perf_fil$x_minBT
tidy_perf_fil$max_RGRBT <- tidy_perf_fil$max_RGR * meansTot$RGR[which(meansTot$species=="fil")] 
write.csv(tidy_perf_fil, "Analysis-output/TPC/Models/tidy_perf_sxp_af_fil.csv")
creds_fil$x <- creds_fil$x + meansTot$dayTemp[which(meansTot$species=="fil")] ; creds_fil$mu <- creds_fil$mu * meansTot$RGR[which(meansTot$species=="fil")] 
creds_fil$upper <- creds_fil$upper * meansTot$RGR[which(meansTot$species=="fil")] ; creds_fil$lower <- creds_fil$lower * meansTot$RGR[which(meansTot$species=="fil")] 
write.csv(creds_fil, "Analysis-output/TPC/Creds/creds_sxp_af_fil.csv")
mean_df_fil <- tidy_perf_fil %>% group_by(species) %>% summarise_if(is.numeric, .funs = c(mean))
mean_df_fil$speciesn <- as.numeric(mean_df_fil$species); mean_df_fil <- mean_df_fil[order(mean_df_fil$speciesn),]
write.csv(mean_df_fil, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_fil.csv")
mean_df_ci_fil <- tidy_perf_fil %>%  group_by(species) %>% summarise(maximaBT_lci = quantile(maximaBT, probs=c(0.025)), maximaBT_uci = quantile(maximaBT, probs=c(0.975)), 
                                                                     B50_lci = quantile(B50, probs=c(0.025)), B50_uci = quantile(B50, probs=c(0.975)), 
                                                                     B50_low_lci = quantile(B50_low, probs=c(0.025)), B50_low_uci = quantile(B50_low, probs=c(0.975)), 
                                                                     B50_high_lci = quantile(B50_high, probs=c(0.025)), B50_high_uci = quantile(B50_high, probs=c(0.975)), 
                                                                     breadthBT_lci = quantile(breadthBT, probs=c(0.025)), breadthBT_uci = quantile(breadthBT, probs=c(0.975)), 
                                                                     x_minBT_lci = quantile(x_minBT, probs=c(0.025)), x_minBT_uci = quantile(x_minBT, probs=c(0.975)), 
                                                                     x_maxBT_lci = quantile(x_maxBT, probs=c(0.025)), x_maxBT_uci = quantile(x_maxBT, probs=c(0.975)), 
                                                                     max_RGR_lci = quantile(max_RGRBT, probs=c(0.025)), max_RGR_uci = quantile(max_RGRBT, probs=c(0.975)), 
                                                                     area_lci = quantile(area, probs=c(0.025)), area_uci = quantile(area, probs=c(0.975)))
mean_df_ci_fil$speciesn <- as.numeric(mean_df_ci_fil$species); mean_df_ci_fil <- mean_df_ci_fil[order(mean_df_ci_fil$speciesn),]
write.csv(mean_df_ci_fil, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_fil.csv")
mean_ci_table_fil <- mean_df_fil %>% right_join(mean_df_ci_fil, by="species")
mean_ci_table_fil <- as.data.frame(mean_ci_table_fil)
mean_ci_table_fil <- mean_ci_table_fil[,c("species", "speciesn.x", 
                                          "maximaBT", "maximaBT_lci", "maximaBT_uci", 
                                          "B50", "B50_lci", "B50_uci", 
                                          "B50_low", "B50_low_lci", "B50_low_uci", 
                                          "B50_high", "B50_high_lci", "B50_high_uci", 
                                          "breadthBT", "breadthBT_lci", "breadthBT_uci", 
                                          "x_minBT", "x_minBT_lci", "x_minBT_uci", 
                                          "x_maxBT", "x_maxBT_lci", "x_maxBT_uci", 
                                          "max_RGRBT", "max_RGR_lci", "max_RGR_uci", 
                                          "area", "area_lci", "area_uci")]
mean_ci_table_fil <- mean_ci_table_fil %>% mutate_at(vars(-c("species","speciesn.x")), funs(round(., 2)))
mean_ci_table_fil$ToptCI <- paste(mean_ci_table_fil$maximaBT, " [", mean_ci_table_fil$maximaBT_lci, ", ", mean_ci_table_fil$maximaBT_uci,"]", sep="")
mean_ci_table_fil$B50CI <- paste(mean_ci_table_fil$B50, " [", mean_ci_table_fil$B50_lci, ", ", mean_ci_table_fil$B50_uci,"]", sep="")
mean_ci_table_fil$B50_lowCI <- paste(mean_ci_table_fil$B50_low, " [", mean_ci_table_fil$B50_low_lci, ", ", mean_ci_table_fil$B50_low_uci,"]", sep="")
mean_ci_table_fil$B50_highCI <- paste(mean_ci_table_fil$B50_high, " [", mean_ci_table_fil$B50_high_lci, ", ", mean_ci_table_fil$B50_high_uci,"]", sep="")
mean_ci_table_fil$breadthCI <- paste(mean_ci_table_fil$breadthBT, " [", mean_ci_table_fil$breadthBT_lci, ", ", mean_ci_table_fil$breadthBT_uci,"]", sep="")
mean_ci_table_fil$x_minCI <- paste(mean_ci_table_fil$x_minBT, " [", mean_ci_table_fil$x_minBT_lci, ", ", mean_ci_table_fil$x_minBT_uci,"]", sep="")
mean_ci_table_fil$x_maxCI <- paste(mean_ci_table_fil$x_maxBT, " [", mean_ci_table_fil$x_maxBT_lci, ", ", mean_ci_table_fil$x_maxBT_uci,"]", sep="")
mean_ci_table_fil$max_RGRCI <- paste(mean_ci_table_fil$max_RGRBT, " [", mean_ci_table_fil$max_RGR_lci, ", ", mean_ci_table_fil$max_RGR_uci,"]", sep="")
mean_ci_table_fil$areaCI <- paste(mean_ci_table_fil$area, " [", mean_ci_table_fil$area_lci, ", ", mean_ci_table_fil$area_uci,"]", sep="")
mean_ci_table_fil <- mean_ci_table_fil[,c("species", "ToptCI", "B50CI", "B50_lowCI", "B50_highCI", "breadthCI", "x_minCI", "x_maxCI", "max_RGRCI", "areaCI")]
write.csv(mean_ci_table_fil, "Analysis-output/TPC/Mean-dfs/mean_ci_sxp_af_table_fil.csv")


tidy_perf_flo <- performr::perform_df(model_fits_flo, species_order = c(1:length(unique(flo_dat$familyI)))) 
tidy_perf_flo %<>% ungroup() %>% mutate(idx = 1:n()) %>% group_by(idx) %>% do({ mutate(.data = ., max_RGR = performance_mu(xs = .$maxima, .$shape1, .$shape2, .$stretch, .$x_min, .$x_max))})
#ssq_df_flo <- bayes_p(stan_df = tidy_perf_flo, raw_df = flo_dat, raw_treatment = "dayTempc", raw_response = "RGRs", raw_group = "familyI")
tidy_perf_flo %<>% ungroup() %>% bind_cols(optimum_breadth(par_df = tidy_perf_flo,  prop_max = 0.5,  n_grid = 100))
tidy_perf_flo$B50_low <- tidy_perf_flo$opt_breadth_low + meansTot$dayTemp[which(meansTot$species=="flo")] 
tidy_perf_flo$B50_high <- tidy_perf_flo$opt_breadth_high + meansTot$dayTemp[which(meansTot$species=="flo")]  
tidy_perf_flo$B50 <- tidy_perf_flo$B50_high - tidy_perf_flo$B50_low 
draws_flo <- rstan::extract(model_fits_flo); ndraws_flo <- length(draws_flo$lp__)
x_seq_flo = seq(min(draws_flo$x_min),max(draws_flo$x_max),length.out = 100)
poly_draws_flo <- sample(1:ndraws_flo, 100)
creds_flo <-  predict_interval(x = x_seq_flo, spp = species, par_df = tidy_perf_flo, x_draws = poly_draws_flo, p = c(0.95, 0.5))
tidy_perf_flo$x_minBT <- tidy_perf_flo$x_min + meansTot$dayTemp[which(meansTot$species=="flo")] 
tidy_perf_flo$x_maxBT <- tidy_perf_flo$x_max + meansTot$dayTemp[which(meansTot$species=="flo")] 
tidy_perf_flo$maximaBT <- tidy_perf_flo$maxima + meansTot$dayTemp[which(meansTot$species=="flo")] 
tidy_perf_flo$breadthBT <- tidy_perf_flo$x_maxBT-tidy_perf_flo$x_minBT
tidy_perf_flo$max_RGRBT <- tidy_perf_flo$max_RGR * meansTot$RGR[which(meansTot$species=="flo")] 
write.csv(tidy_perf_flo, "Analysis-output/TPC/Models/tidy_perf_sxp_af_flo.csv")
creds_flo$x <- creds_flo$x + meansTot$dayTemp[which(meansTot$species=="flo")] ; creds_flo$mu <- creds_flo$mu * meansTot$RGR[which(meansTot$species=="flo")] 
creds_flo$upper <- creds_flo$upper * meansTot$RGR[which(meansTot$species=="flo")] ; creds_flo$lower <- creds_flo$lower * meansTot$RGR[which(meansTot$species=="flo")] 
write.csv(creds_flo, "Analysis-output/TPC/Creds/creds_sxp_af_flo.csv")
mean_df_flo <- tidy_perf_flo %>% group_by(species) %>% summarise_if(is.numeric, .funs = c(mean))
mean_df_flo$speciesn <- as.numeric(mean_df_flo$species); mean_df_flo <- mean_df_flo[order(mean_df_flo$speciesn),]
write.csv(mean_df_flo, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_flo.csv")
mean_df_ci_flo <- tidy_perf_flo %>%  group_by(species) %>% summarise(maximaBT_lci = quantile(maximaBT, probs=c(0.025)), maximaBT_uci = quantile(maximaBT, probs=c(0.975)), 
                                                                     B50_lci = quantile(B50, probs=c(0.025)), B50_uci = quantile(B50, probs=c(0.975)), 
                                                                     B50_low_lci = quantile(B50_low, probs=c(0.025)), B50_low_uci = quantile(B50_low, probs=c(0.975)), 
                                                                     B50_high_lci = quantile(B50_high, probs=c(0.025)), B50_high_uci = quantile(B50_high, probs=c(0.975)), 
                                                                     breadthBT_lci = quantile(breadthBT, probs=c(0.025)), breadthBT_uci = quantile(breadthBT, probs=c(0.975)), 
                                                                     x_minBT_lci = quantile(x_minBT, probs=c(0.025)), x_minBT_uci = quantile(x_minBT, probs=c(0.975)), 
                                                                     x_maxBT_lci = quantile(x_maxBT, probs=c(0.025)), x_maxBT_uci = quantile(x_maxBT, probs=c(0.975)), 
                                                                     max_RGR_lci = quantile(max_RGRBT, probs=c(0.025)), max_RGR_uci = quantile(max_RGRBT, probs=c(0.975)), 
                                                                     area_lci = quantile(area, probs=c(0.025)), area_uci = quantile(area, probs=c(0.975)))
mean_df_ci_flo$speciesn <- as.numeric(mean_df_ci_flo$species); mean_df_ci_flo <- mean_df_ci_flo[order(mean_df_ci_flo$speciesn),]
write.csv(mean_df_ci_flo, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_flo.csv")
mean_ci_table_flo <- mean_df_flo %>% right_join(mean_df_ci_flo, by="species")
mean_ci_table_flo <- as.data.frame(mean_ci_table_flo)
mean_ci_table_flo <- mean_ci_table_flo[,c("species", "speciesn.x", 
                                          "maximaBT", "maximaBT_lci", "maximaBT_uci", 
                                          "B50", "B50_lci", "B50_uci", 
                                          "B50_low", "B50_low_lci", "B50_low_uci", 
                                          "B50_high", "B50_high_lci", "B50_high_uci", 
                                          "breadthBT", "breadthBT_lci", "breadthBT_uci", 
                                          "x_minBT", "x_minBT_lci", "x_minBT_uci", 
                                          "x_maxBT", "x_maxBT_lci", "x_maxBT_uci", 
                                          "max_RGRBT", "max_RGR_lci", "max_RGR_uci", 
                                          "area", "area_lci", "area_uci")]
mean_ci_table_flo <- mean_ci_table_flo %>% mutate_at(vars(-c("species","speciesn.x")), funs(round(., 2)))
mean_ci_table_flo$ToptCI <- paste(mean_ci_table_flo$maximaBT, " [", mean_ci_table_flo$maximaBT_lci, ", ", mean_ci_table_flo$maximaBT_uci,"]", sep="")
mean_ci_table_flo$B50CI <- paste(mean_ci_table_flo$B50, " [", mean_ci_table_flo$B50_lci, ", ", mean_ci_table_flo$B50_uci,"]", sep="")
mean_ci_table_flo$B50_lowCI <- paste(mean_ci_table_flo$B50_low, " [", mean_ci_table_flo$B50_low_lci, ", ", mean_ci_table_flo$B50_low_uci,"]", sep="")
mean_ci_table_flo$B50_highCI <- paste(mean_ci_table_flo$B50_high, " [", mean_ci_table_flo$B50_high_lci, ", ", mean_ci_table_flo$B50_high_uci,"]", sep="")
mean_ci_table_flo$breadthCI <- paste(mean_ci_table_flo$breadthBT, " [", mean_ci_table_flo$breadthBT_lci, ", ", mean_ci_table_flo$breadthBT_uci,"]", sep="")
mean_ci_table_flo$x_minCI <- paste(mean_ci_table_flo$x_minBT, " [", mean_ci_table_flo$x_minBT_lci, ", ", mean_ci_table_flo$x_minBT_uci,"]", sep="")
mean_ci_table_flo$x_maxCI <- paste(mean_ci_table_flo$x_maxBT, " [", mean_ci_table_flo$x_maxBT_lci, ", ", mean_ci_table_flo$x_maxBT_uci,"]", sep="")
mean_ci_table_flo$max_RGRCI <- paste(mean_ci_table_flo$max_RGRBT, " [", mean_ci_table_flo$max_RGR_lci, ", ", mean_ci_table_flo$max_RGR_uci,"]", sep="")
mean_ci_table_flo$areaCI <- paste(mean_ci_table_flo$area, " [", mean_ci_table_flo$area_lci, ", ", mean_ci_table_flo$area_uci,"]", sep="")
mean_ci_table_flo <- mean_ci_table_flo[,c("species", "ToptCI", "B50CI", "B50_lowCI", "B50_highCI", "breadthCI", "x_minCI", "x_maxCI", "max_RGRCI", "areaCI")]
write.csv(mean_ci_table_flo, "Analysis-output/TPC/Mean-dfs/mean_ci_sxp_af_table_flo.csv")


tidy_perf_gut <- performr::perform_df(model_fits_gut, species_order = c(1:length(unique(gut_dat$familyI)))) 
tidy_perf_gut %<>% ungroup() %>% mutate(idx = 1:n()) %>% group_by(idx) %>% do({ mutate(.data = ., max_RGR = performance_mu(xs = .$maxima, .$shape1, .$shape2, .$stretch, .$x_min, .$x_max))})
#ssq_df_gut <- bayes_p(stan_df = tidy_perf_gut, raw_df = gut_dat, raw_treatment = "dayTempc", raw_response = "RGRs", raw_group = "familyI")
tidy_perf_gut %<>% ungroup() %>% bind_cols(optimum_breadth(par_df = tidy_perf_gut,  prop_max = 0.5,  n_grid = 100))
tidy_perf_gut$B50_low <- tidy_perf_gut$opt_breadth_low + meansTot$dayTemp[which(meansTot$species=="gut")] 
tidy_perf_gut$B50_high <- tidy_perf_gut$opt_breadth_high + meansTot$dayTemp[which(meansTot$species=="gut")]  
tidy_perf_gut$B50 <- tidy_perf_gut$B50_high - tidy_perf_gut$B50_low 
draws_gut <- rstan::extract(model_fits_gut); ndraws_gut <- length(draws_gut$lp__)
x_seq_gut = seq(min(draws_gut$x_min),max(draws_gut$x_max),length.out = 100)
poly_draws_gut <- sample(1:ndraws_gut, 100)
creds_gut <-  predict_interval(x = x_seq_gut, spp = species, par_df = tidy_perf_gut, x_draws = poly_draws_gut, p = c(0.95, 0.5))
tidy_perf_gut$x_minBT <- tidy_perf_gut$x_min + meansTot$dayTemp[which(meansTot$species=="gut")] 
tidy_perf_gut$x_maxBT <- tidy_perf_gut$x_max + meansTot$dayTemp[which(meansTot$species=="gut")] 
tidy_perf_gut$maximaBT <- tidy_perf_gut$maxima + meansTot$dayTemp[which(meansTot$species=="gut")] 
tidy_perf_gut$breadthBT <- tidy_perf_gut$x_maxBT-tidy_perf_gut$x_minBT
tidy_perf_gut$max_RGRBT <- tidy_perf_gut$max_RGR * meansTot$RGR[which(meansTot$species=="gut")] 
write.csv(tidy_perf_gut, "Analysis-output/TPC/Models/tidy_perf_sxp_af_gut.csv")
creds_gut$x <- creds_gut$x + meansTot$dayTemp[which(meansTot$species=="gut")] ; creds_gut$mu <- creds_gut$mu * meansTot$RGR[which(meansTot$species=="gut")] 
creds_gut$upper <- creds_gut$upper * meansTot$RGR[which(meansTot$species=="gut")] ; creds_gut$lower <- creds_gut$lower * meansTot$RGR[which(meansTot$species=="gut")] 
write.csv(creds_gut, "Analysis-output/TPC/Creds/creds_sxp_af_gut.csv")
mean_df_gut <- tidy_perf_gut %>% group_by(species) %>% summarise_if(is.numeric, .funs = c(mean))
mean_df_gut$speciesn <- as.numeric(mean_df_gut$species); mean_df_gut <- mean_df_gut[order(mean_df_gut$speciesn),]
write.csv(mean_df_gut, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_gut.csv")
mean_df_ci_gut <- tidy_perf_gut %>%  group_by(species) %>% summarise(maximaBT_lci = quantile(maximaBT, probs=c(0.025)), maximaBT_uci = quantile(maximaBT, probs=c(0.975)), 
                                                                     B50_lci = quantile(B50, probs=c(0.025)), B50_uci = quantile(B50, probs=c(0.975)), 
                                                                     B50_low_lci = quantile(B50_low, probs=c(0.025)), B50_low_uci = quantile(B50_low, probs=c(0.975)), 
                                                                     B50_high_lci = quantile(B50_high, probs=c(0.025)), B50_high_uci = quantile(B50_high, probs=c(0.975)), 
                                                                     breadthBT_lci = quantile(breadthBT, probs=c(0.025)), breadthBT_uci = quantile(breadthBT, probs=c(0.975)), 
                                                                     x_minBT_lci = quantile(x_minBT, probs=c(0.025)), x_minBT_uci = quantile(x_minBT, probs=c(0.975)), 
                                                                     x_maxBT_lci = quantile(x_maxBT, probs=c(0.025)), x_maxBT_uci = quantile(x_maxBT, probs=c(0.975)), 
                                                                     max_RGR_lci = quantile(max_RGRBT, probs=c(0.025)), max_RGR_uci = quantile(max_RGRBT, probs=c(0.975)), 
                                                                     area_lci = quantile(area, probs=c(0.025)), area_uci = quantile(area, probs=c(0.975)))
mean_df_ci_gut$speciesn <- as.numeric(mean_df_ci_gut$species); mean_df_ci_gut <- mean_df_ci_gut[order(mean_df_ci_gut$speciesn),]
write.csv(mean_df_ci_gut, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_gut.csv")
mean_ci_table_gut <- mean_df_gut %>% right_join(mean_df_ci_gut, by="species")
mean_ci_table_gut <- as.data.frame(mean_ci_table_gut)
mean_ci_table_gut <- mean_ci_table_gut[,c("species", "speciesn.x", 
                                          "maximaBT", "maximaBT_lci", "maximaBT_uci", 
                                          "B50", "B50_lci", "B50_uci", 
                                          "B50_low", "B50_low_lci", "B50_low_uci", 
                                          "B50_high", "B50_high_lci", "B50_high_uci", 
                                          "breadthBT", "breadthBT_lci", "breadthBT_uci", 
                                          "x_minBT", "x_minBT_lci", "x_minBT_uci", 
                                          "x_maxBT", "x_maxBT_lci", "x_maxBT_uci", 
                                          "max_RGRBT", "max_RGR_lci", "max_RGR_uci", 
                                          "area", "area_lci", "area_uci")]
mean_ci_table_gut <- mean_ci_table_gut %>% mutate_at(vars(-c("species","speciesn.x")), funs(round(., 2)))
mean_ci_table_gut$ToptCI <- paste(mean_ci_table_gut$maximaBT, " [", mean_ci_table_gut$maximaBT_lci, ", ", mean_ci_table_gut$maximaBT_uci,"]", sep="")
mean_ci_table_gut$B50CI <- paste(mean_ci_table_gut$B50, " [", mean_ci_table_gut$B50_lci, ", ", mean_ci_table_gut$B50_uci,"]", sep="")
mean_ci_table_gut$B50_lowCI <- paste(mean_ci_table_gut$B50_low, " [", mean_ci_table_gut$B50_low_lci, ", ", mean_ci_table_gut$B50_low_uci,"]", sep="")
mean_ci_table_gut$B50_highCI <- paste(mean_ci_table_gut$B50_high, " [", mean_ci_table_gut$B50_high_lci, ", ", mean_ci_table_gut$B50_high_uci,"]", sep="")
mean_ci_table_gut$breadthCI <- paste(mean_ci_table_gut$breadthBT, " [", mean_ci_table_gut$breadthBT_lci, ", ", mean_ci_table_gut$breadthBT_uci,"]", sep="")
mean_ci_table_gut$x_minCI <- paste(mean_ci_table_gut$x_minBT, " [", mean_ci_table_gut$x_minBT_lci, ", ", mean_ci_table_gut$x_minBT_uci,"]", sep="")
mean_ci_table_gut$x_maxCI <- paste(mean_ci_table_gut$x_maxBT, " [", mean_ci_table_gut$x_maxBT_lci, ", ", mean_ci_table_gut$x_maxBT_uci,"]", sep="")
mean_ci_table_gut$max_RGRCI <- paste(mean_ci_table_gut$max_RGRBT, " [", mean_ci_table_gut$max_RGR_lci, ", ", mean_ci_table_gut$max_RGR_uci,"]", sep="")
mean_ci_table_gut$areaCI <- paste(mean_ci_table_gut$area, " [", mean_ci_table_gut$area_lci, ", ", mean_ci_table_gut$area_uci,"]", sep="")
mean_ci_table_gut <- mean_ci_table_gut[,c("species", "ToptCI", "B50CI", "B50_lowCI", "B50_highCI", "breadthCI", "x_minCI", "x_maxCI", "max_RGRCI", "areaCI")]
write.csv(mean_ci_table_gut, "Analysis-output/TPC/Mean-dfs/mean_ci_sxp_af_table_gut.csv")


tidy_perf_lac <- performr::perform_df(model_fits_lac, species_order = c(1:length(unique(lac_dat$familyI)))) 
tidy_perf_lac %<>% ungroup() %>% mutate(idx = 1:n()) %>% group_by(idx) %>% do({ mutate(.data = ., max_RGR = performance_mu(xs = .$maxima, .$shape1, .$shape2, .$stretch, .$x_min, .$x_max))})
#ssq_df_lac <- bayes_p(stan_df = tidy_perf_lac, raw_df = lac_dat, raw_treatment = "dayTempc", raw_response = "RGRs", raw_group = "familyI")
tidy_perf_lac %<>% ungroup() %>% bind_cols(optimum_breadth(par_df = tidy_perf_lac,  prop_max = 0.5,  n_grid = 100))
tidy_perf_lac$B50_low <- tidy_perf_lac$opt_breadth_low + meansTot$dayTemp[which(meansTot$species=="lac")] 
tidy_perf_lac$B50_high <- tidy_perf_lac$opt_breadth_high + meansTot$dayTemp[which(meansTot$species=="lac")]  
tidy_perf_lac$B50 <- tidy_perf_lac$B50_high - tidy_perf_lac$B50_low 
tidy_perf_lac %<>% ungroup() %>% bind_cols(optimum_breadth(par_df = tidy_perf_lac,  prop_max = 0.5,  n_grid = 100))
draws_lac <- rstan::extract(model_fits_lac); ndraws_lac <- length(draws_lac$lp__)
x_seq_lac = seq(min(draws_lac$x_min),max(draws_lac$x_max),length.out = 100)
poly_draws_lac <- sample(1:ndraws_lac, 100)
creds_lac <-  predict_interval(x = x_seq_lac, spp = species, par_df = tidy_perf_lac, x_draws = poly_draws_lac, p = c(0.95, 0.5))
tidy_perf_lac$x_minBT <- tidy_perf_lac$x_min + meansTot$dayTemp[which(meansTot$species=="lac")] 
tidy_perf_lac$x_maxBT <- tidy_perf_lac$x_max + meansTot$dayTemp[which(meansTot$species=="lac")] 
tidy_perf_lac$maximaBT <- tidy_perf_lac$maxima + meansTot$dayTemp[which(meansTot$species=="lac")] 
tidy_perf_lac$breadthBT <- tidy_perf_lac$x_maxBT-tidy_perf_lac$x_minBT
tidy_perf_lac$max_RGRBT <- tidy_perf_lac$max_RGR * meansTot$RGR[which(meansTot$species=="lac")] 
write.csv(tidy_perf_lac, "Analysis-output/TPC/Models/tidy_perf_sxp_af_lac.csv")
creds_lac$x <- creds_lac$x + meansTot$dayTemp[which(meansTot$species=="lac")] ; creds_lac$mu <- creds_lac$mu * meansTot$RGR[which(meansTot$species=="lac")] 
creds_lac$upper <- creds_lac$upper * meansTot$RGR[which(meansTot$species=="lac")] ; creds_lac$lower <- creds_lac$lower * meansTot$RGR[which(meansTot$species=="lac")] 
write.csv(creds_lac, "Analysis-output/TPC/Creds/creds_sxp_af_lac.csv")
mean_df_lac <- tidy_perf_lac %>% group_by(species) %>% summarise_if(is.numeric, .funs = c(mean))
mean_df_lac$speciesn <- as.numeric(mean_df_lac$species); mean_df_lac <- mean_df_lac[order(mean_df_lac$speciesn),]
write.csv(mean_df_lac, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_lac.csv")
mean_df_ci_lac <- tidy_perf_lac %>%  group_by(species) %>% summarise(maximaBT_lci = quantile(maximaBT, probs=c(0.025)), maximaBT_uci = quantile(maximaBT, probs=c(0.975)), 
                                                                     B50_lci = quantile(B50, probs=c(0.025)), B50_uci = quantile(B50, probs=c(0.975)), 
                                                                     B50_low_lci = quantile(B50_low, probs=c(0.025)), B50_low_uci = quantile(B50_low, probs=c(0.975)), 
                                                                     B50_high_lci = quantile(B50_high, probs=c(0.025)), B50_high_uci = quantile(B50_high, probs=c(0.975)), 
                                                                     breadthBT_lci = quantile(breadthBT, probs=c(0.025)), breadthBT_uci = quantile(breadthBT, probs=c(0.975)), 
                                                                     x_minBT_lci = quantile(x_minBT, probs=c(0.025)), x_minBT_uci = quantile(x_minBT, probs=c(0.975)), 
                                                                     x_maxBT_lci = quantile(x_maxBT, probs=c(0.025)), x_maxBT_uci = quantile(x_maxBT, probs=c(0.975)), 
                                                                     max_RGR_lci = quantile(max_RGRBT, probs=c(0.025)), max_RGR_uci = quantile(max_RGRBT, probs=c(0.975)), 
                                                                     area_lci = quantile(area, probs=c(0.025)), area_uci = quantile(area, probs=c(0.975)))
mean_df_ci_lac$speciesn <- as.numeric(mean_df_ci_lac$species); mean_df_ci_lac <- mean_df_ci_lac[order(mean_df_ci_lac$speciesn),]
write.csv(mean_df_ci_lac, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_lac.csv")
mean_ci_table_lac <- mean_df_lac %>% right_join(mean_df_ci_lac, by="species")
mean_ci_table_lac <- as.data.frame(mean_ci_table_lac)
mean_ci_table_lac <- mean_ci_table_lac[,c("species", "speciesn.x", 
                                          "maximaBT", "maximaBT_lci", "maximaBT_uci", 
                                          "B50", "B50_lci", "B50_uci", 
                                          "B50_low", "B50_low_lci", "B50_low_uci", 
                                          "B50_high", "B50_high_lci", "B50_high_uci", 
                                          "breadthBT", "breadthBT_lci", "breadthBT_uci", 
                                          "x_minBT", "x_minBT_lci", "x_minBT_uci", 
                                          "x_maxBT", "x_maxBT_lci", "x_maxBT_uci", 
                                          "max_RGRBT", "max_RGR_lci", "max_RGR_uci", 
                                          "area", "area_lci", "area_uci")]
mean_ci_table_lac <- mean_ci_table_lac %>% mutate_at(vars(-c("species","speciesn.x")), funs(round(., 2)))
mean_ci_table_lac$ToptCI <- paste(mean_ci_table_lac$maximaBT, " [", mean_ci_table_lac$maximaBT_lci, ", ", mean_ci_table_lac$maximaBT_uci,"]", sep="")
mean_ci_table_lac$B50CI <- paste(mean_ci_table_lac$B50, " [", mean_ci_table_lac$B50_lci, ", ", mean_ci_table_lac$B50_uci,"]", sep="")
mean_ci_table_lac$B50_lowCI <- paste(mean_ci_table_lac$B50_low, " [", mean_ci_table_lac$B50_low_lci, ", ", mean_ci_table_lac$B50_low_uci,"]", sep="")
mean_ci_table_lac$B50_highCI <- paste(mean_ci_table_lac$B50_high, " [", mean_ci_table_lac$B50_high_lci, ", ", mean_ci_table_lac$B50_high_uci,"]", sep="")
mean_ci_table_lac$breadthCI <- paste(mean_ci_table_lac$breadthBT, " [", mean_ci_table_lac$breadthBT_lci, ", ", mean_ci_table_lac$breadthBT_uci,"]", sep="")
mean_ci_table_lac$x_minCI <- paste(mean_ci_table_lac$x_minBT, " [", mean_ci_table_lac$x_minBT_lci, ", ", mean_ci_table_lac$x_minBT_uci,"]", sep="")
mean_ci_table_lac$x_maxCI <- paste(mean_ci_table_lac$x_maxBT, " [", mean_ci_table_lac$x_maxBT_lci, ", ", mean_ci_table_lac$x_maxBT_uci,"]", sep="")
mean_ci_table_lac$max_RGRCI <- paste(mean_ci_table_lac$max_RGRBT, " [", mean_ci_table_lac$max_RGR_lci, ", ", mean_ci_table_lac$max_RGR_uci,"]", sep="")
mean_ci_table_lac$areaCI <- paste(mean_ci_table_lac$area, " [", mean_ci_table_lac$area_lci, ", ", mean_ci_table_lac$area_uci,"]", sep="")
mean_ci_table_lac <- mean_ci_table_lac[,c("species", "ToptCI", "B50CI", "B50_lowCI", "B50_highCI", "breadthCI", "x_minCI", "x_maxCI", "max_RGRCI", "areaCI")]
write.csv(mean_ci_table_lac, "Analysis-output/TPC/Mean-dfs/mean_ci_sxp_af_table_lac.csv")


tidy_perf_nor <- performr::perform_df(model_fits_nor, species_order = c(1:length(unique(nor_dat$familyI)))) 
tidy_perf_nor %<>% ungroup() %>% mutate(idx = 1:n()) %>% group_by(idx) %>% do({ mutate(.data = ., max_RGR = performance_mu(xs = .$maxima, .$shape1, .$shape2, .$stretch, .$x_min, .$x_max))})
#ssq_df_nor <- bayes_p(stan_df = tidy_perf_nor, raw_df = nor_dat, raw_treatment = "dayTempc", raw_response = "RGRs", raw_group = "familyI")
tidy_perf_nor %<>% ungroup() %>% bind_cols(optimum_breadth(par_df = tidy_perf_nor,  prop_max = 0.5,  n_grid = 100))
tidy_perf_nor$B50_low <- tidy_perf_nor$opt_breadth_low + meansTot$dayTemp[which(meansTot$species=="nor")] 
tidy_perf_nor$B50_high <- tidy_perf_nor$opt_breadth_high + meansTot$dayTemp[which(meansTot$species=="nor")]  
tidy_perf_nor$B50 <- tidy_perf_nor$B50_high - tidy_perf_nor$B50_low 
draws_nor <- rstan::extract(model_fits_nor); ndraws_nor <- length(draws_nor$lp__)
x_seq_nor = seq(min(draws_nor$x_min),max(draws_nor$x_max),length.out = 100)
poly_draws_nor <- sample(1:ndraws_nor, 100)
creds_nor <-  predict_interval(x = x_seq_nor, spp = species, par_df = tidy_perf_nor, x_draws = poly_draws_nor, p = c(0.95, 0.5))
tidy_perf_nor$x_minBT <- tidy_perf_nor$x_min + meansTot$dayTemp[which(meansTot$species=="nor")] 
tidy_perf_nor$x_maxBT <- tidy_perf_nor$x_max + meansTot$dayTemp[which(meansTot$species=="nor")] 
tidy_perf_nor$maximaBT <- tidy_perf_nor$maxima + meansTot$dayTemp[which(meansTot$species=="nor")] 
tidy_perf_nor$breadthBT <- tidy_perf_nor$x_maxBT-tidy_perf_nor$x_minBT
tidy_perf_nor$max_RGRBT <- tidy_perf_nor$max_RGR * meansTot$RGR[which(meansTot$species=="nor")] 
write.csv(tidy_perf_nor, "Analysis-output/TPC/Models/tidy_perf_sxp_af_nor.csv")
creds_nor$x <- creds_nor$x + meansTot$dayTemp[which(meansTot$species=="nor")] ; creds_nor$mu <- creds_nor$mu * meansTot$RGR[which(meansTot$species=="nor")] 
creds_nor$upper <- creds_nor$upper * meansTot$RGR[which(meansTot$species=="nor")] ; creds_nor$lower <- creds_nor$lower * meansTot$RGR[which(meansTot$species=="nor")] 
write.csv(creds_nor, "Analysis-output/TPC/Creds/creds_sxp_af_nor.csv")
mean_df_nor <- tidy_perf_nor %>% group_by(species) %>% summarise_if(is.numeric, .funs = c(mean))
mean_df_nor$speciesn <- as.numeric(mean_df_nor$species); mean_df_nor <- mean_df_nor[order(mean_df_nor$speciesn),]
write.csv(mean_df_nor, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_nor.csv")
mean_df_ci_nor <- tidy_perf_nor %>%  group_by(species) %>% summarise(maximaBT_lci = quantile(maximaBT, probs=c(0.025)), maximaBT_uci = quantile(maximaBT, probs=c(0.975)), 
                                                                     B50_lci = quantile(B50, probs=c(0.025)), B50_uci = quantile(B50, probs=c(0.975)), 
                                                                     B50_low_lci = quantile(B50_low, probs=c(0.025)), B50_low_uci = quantile(B50_low, probs=c(0.975)), 
                                                                     B50_high_lci = quantile(B50_high, probs=c(0.025)), B50_high_uci = quantile(B50_high, probs=c(0.975)), 
                                                                     breadthBT_lci = quantile(breadthBT, probs=c(0.025)), breadthBT_uci = quantile(breadthBT, probs=c(0.975)), 
                                                                     x_minBT_lci = quantile(x_minBT, probs=c(0.025)), x_minBT_uci = quantile(x_minBT, probs=c(0.975)), 
                                                                     x_maxBT_lci = quantile(x_maxBT, probs=c(0.025)), x_maxBT_uci = quantile(x_maxBT, probs=c(0.975)), 
                                                                     max_RGR_lci = quantile(max_RGRBT, probs=c(0.025)), max_RGR_uci = quantile(max_RGRBT, probs=c(0.975)), 
                                                                     area_lci = quantile(area, probs=c(0.025)), area_uci = quantile(area, probs=c(0.975)))
mean_df_ci_nor$speciesn <- as.numeric(mean_df_ci_nor$species); mean_df_ci_nor <- mean_df_ci_nor[order(mean_df_ci_nor$speciesn),]
write.csv(mean_df_ci_nor, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_nor.csv")
mean_ci_table_nor <- mean_df_nor %>% right_join(mean_df_ci_nor, by="species")
mean_ci_table_nor <- as.data.frame(mean_ci_table_nor)
mean_ci_table_nor <- mean_ci_table_nor[,c("species", "speciesn.x", 
                                          "maximaBT", "maximaBT_lci", "maximaBT_uci", 
                                          "B50", "B50_lci", "B50_uci", 
                                          "B50_low", "B50_low_lci", "B50_low_uci", 
                                          "B50_high", "B50_high_lci", "B50_high_uci", 
                                          "breadthBT", "breadthBT_lci", "breadthBT_uci", 
                                          "x_minBT", "x_minBT_lci", "x_minBT_uci", 
                                          "x_maxBT", "x_maxBT_lci", "x_maxBT_uci", 
                                          "max_RGRBT", "max_RGR_lci", "max_RGR_uci", 
                                          "area", "area_lci", "area_uci")]
mean_ci_table_nor <- mean_ci_table_nor %>% mutate_at(vars(-c("species","speciesn.x")), funs(round(., 2)))
mean_ci_table_nor$ToptCI <- paste(mean_ci_table_nor$maximaBT, " [", mean_ci_table_nor$maximaBT_lci, ", ", mean_ci_table_nor$maximaBT_uci,"]", sep="")
mean_ci_table_nor$B50CI <- paste(mean_ci_table_nor$B50, " [", mean_ci_table_nor$B50_lci, ", ", mean_ci_table_nor$B50_uci,"]", sep="")
mean_ci_table_nor$B50_lowCI <- paste(mean_ci_table_nor$B50_low, " [", mean_ci_table_nor$B50_low_lci, ", ", mean_ci_table_nor$B50_low_uci,"]", sep="")
mean_ci_table_nor$B50_highCI <- paste(mean_ci_table_nor$B50_high, " [", mean_ci_table_nor$B50_high_lci, ", ", mean_ci_table_nor$B50_high_uci,"]", sep="")
mean_ci_table_nor$breadthCI <- paste(mean_ci_table_nor$breadthBT, " [", mean_ci_table_nor$breadthBT_lci, ", ", mean_ci_table_nor$breadthBT_uci,"]", sep="")
mean_ci_table_nor$x_minCI <- paste(mean_ci_table_nor$x_minBT, " [", mean_ci_table_nor$x_minBT_lci, ", ", mean_ci_table_nor$x_minBT_uci,"]", sep="")
mean_ci_table_nor$x_maxCI <- paste(mean_ci_table_nor$x_maxBT, " [", mean_ci_table_nor$x_maxBT_lci, ", ", mean_ci_table_nor$x_maxBT_uci,"]", sep="")
mean_ci_table_nor$max_RGRCI <- paste(mean_ci_table_nor$max_RGRBT, " [", mean_ci_table_nor$max_RGR_lci, ", ", mean_ci_table_nor$max_RGR_uci,"]", sep="")
mean_ci_table_nor$areaCI <- paste(mean_ci_table_nor$area, " [", mean_ci_table_nor$area_lci, ", ", mean_ci_table_nor$area_uci,"]", sep="")
mean_ci_table_nor <- mean_ci_table_nor[,c("species", "ToptCI", "B50CI", "B50_lowCI", "B50_highCI", "breadthCI", "x_minCI", "x_maxCI", "max_RGRCI", "areaCI")]
write.csv(mean_ci_table_nor, "Analysis-output/TPC/Mean-dfs/mean_ci_sxp_af_table_nor.csv")


tidy_perf_par <- performr::perform_df(model_fits_par, species_order = c(1:length(unique(par_dat$familyI)))) 
tidy_perf_par %<>% ungroup() %>% mutate(idx = 1:n()) %>% group_by(idx) %>% do({ mutate(.data = ., max_RGR = performance_mu(xs = .$maxima, .$shape1, .$shape2, .$stretch, .$x_min, .$x_max))})
#ssq_df_par <- bayes_p(stan_df = tidy_perf_par, raw_df = par_dat, raw_treatment = "dayTempc", raw_response = "RGRs", raw_group = "familyI")
tidy_perf_par %<>% ungroup() %>% bind_cols(optimum_breadth(par_df = tidy_perf_par,  prop_max = 0.5,  n_grid = 100))
tidy_perf_par$B50_low <- tidy_perf_par$opt_breadth_low + meansTot$dayTemp[which(meansTot$species=="par")] 
tidy_perf_par$B50_high <- tidy_perf_par$opt_breadth_high + meansTot$dayTemp[which(meansTot$species=="par")]  
tidy_perf_par$B50 <- tidy_perf_par$B50_high - tidy_perf_par$B50_low 
draws_par <- rstan::extract(model_fits_par); ndraws_par <- length(draws_par$lp__)
x_seq_par = seq(min(draws_par$x_min),max(draws_par$x_max),length.out = 100)
poly_draws_par <- sample(1:ndraws_par, 100)
creds_par <-  predict_interval(x = x_seq_par, spp = species, par_df = tidy_perf_par, x_draws = poly_draws_par, p = c(0.95, 0.5))
tidy_perf_par$x_minBT <- tidy_perf_par$x_min + meansTot$dayTemp[which(meansTot$species=="par")] 
tidy_perf_par$x_maxBT <- tidy_perf_par$x_max + meansTot$dayTemp[which(meansTot$species=="par")] 
tidy_perf_par$maximaBT <- tidy_perf_par$maxima + meansTot$dayTemp[which(meansTot$species=="par")] 
tidy_perf_par$breadthBT <- tidy_perf_par$x_maxBT-tidy_perf_par$x_minBT
tidy_perf_par$max_RGRBT <- tidy_perf_par$max_RGR * meansTot$RGR[which(meansTot$species=="par")] 
write.csv(tidy_perf_par, "Analysis-output/TPC/Models/tidy_perf_sxp_af_par.csv")
creds_par$x <- creds_par$x + meansTot$dayTemp[which(meansTot$species=="par")] ; creds_par$mu <- creds_par$mu * meansTot$RGR[which(meansTot$species=="par")] 
creds_par$upper <- creds_par$upper * meansTot$RGR[which(meansTot$species=="par")] ; creds_par$lower <- creds_par$lower * meansTot$RGR[which(meansTot$species=="par")] 
write.csv(creds_par, "Analysis-output/TPC/Creds/creds_sxp_af_par.csv")
mean_df_par <- tidy_perf_par %>% group_by(species) %>% summarise_if(is.numeric, .funs = c(mean))
mean_df_par$speciesn <- as.numeric(mean_df_par$species); mean_df_par <- mean_df_par[order(mean_df_par$speciesn),]
write.csv(mean_df_par, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_par.csv")
mean_df_ci_par <- tidy_perf_par %>%  group_by(species) %>% summarise(maximaBT_lci = quantile(maximaBT, probs=c(0.025)), maximaBT_uci = quantile(maximaBT, probs=c(0.975)), 
                                                                     B50_lci = quantile(B50, probs=c(0.025)), B50_uci = quantile(B50, probs=c(0.975)), 
                                                                     B50_low_lci = quantile(B50_low, probs=c(0.025)), B50_low_uci = quantile(B50_low, probs=c(0.975)), 
                                                                     B50_high_lci = quantile(B50_high, probs=c(0.025)), B50_high_uci = quantile(B50_high, probs=c(0.975)), 
                                                                     breadthBT_lci = quantile(breadthBT, probs=c(0.025)), breadthBT_uci = quantile(breadthBT, probs=c(0.975)), 
                                                                     x_minBT_lci = quantile(x_minBT, probs=c(0.025)), x_minBT_uci = quantile(x_minBT, probs=c(0.975)), 
                                                                     x_maxBT_lci = quantile(x_maxBT, probs=c(0.025)), x_maxBT_uci = quantile(x_maxBT, probs=c(0.975)), 
                                                                     max_RGR_lci = quantile(max_RGRBT, probs=c(0.025)), max_RGR_uci = quantile(max_RGRBT, probs=c(0.975)), 
                                                                     area_lci = quantile(area, probs=c(0.025)), area_uci = quantile(area, probs=c(0.975)))
mean_df_ci_par$speciesn <- as.numeric(mean_df_ci_par$species); mean_df_ci_par <- mean_df_ci_par[order(mean_df_ci_par$speciesn),]
write.csv(mean_df_ci_par, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_par.csv")
mean_ci_table_par <- mean_df_par %>% right_join(mean_df_ci_par, by="species")
mean_ci_table_par <- as.data.frame(mean_ci_table_par)
mean_ci_table_par <- mean_ci_table_par[,c("species", "speciesn.x", 
                                          "maximaBT", "maximaBT_lci", "maximaBT_uci", 
                                          "B50", "B50_lci", "B50_uci", 
                                          "B50_low", "B50_low_lci", "B50_low_uci", 
                                          "B50_high", "B50_high_lci", "B50_high_uci", 
                                          "breadthBT", "breadthBT_lci", "breadthBT_uci", 
                                          "x_minBT", "x_minBT_lci", "x_minBT_uci", 
                                          "x_maxBT", "x_maxBT_lci", "x_maxBT_uci", 
                                          "max_RGRBT", "max_RGR_lci", "max_RGR_uci", 
                                          "area", "area_lci", "area_uci")]
mean_ci_table_par <- mean_ci_table_par %>% mutate_at(vars(-c("species","speciesn.x")), funs(round(., 2)))
mean_ci_table_par$ToptCI <- paste(mean_ci_table_par$maximaBT, " [", mean_ci_table_par$maximaBT_lci, ", ", mean_ci_table_par$maximaBT_uci,"]", sep="")
mean_ci_table_par$B50CI <- paste(mean_ci_table_par$B50, " [", mean_ci_table_par$B50_lci, ", ", mean_ci_table_par$B50_uci,"]", sep="")
mean_ci_table_par$B50_lowCI <- paste(mean_ci_table_par$B50_low, " [", mean_ci_table_par$B50_low_lci, ", ", mean_ci_table_par$B50_low_uci,"]", sep="")
mean_ci_table_par$B50_highCI <- paste(mean_ci_table_par$B50_high, " [", mean_ci_table_par$B50_high_lci, ", ", mean_ci_table_par$B50_high_uci,"]", sep="")
mean_ci_table_par$breadthCI <- paste(mean_ci_table_par$breadthBT, " [", mean_ci_table_par$breadthBT_lci, ", ", mean_ci_table_par$breadthBT_uci,"]", sep="")
mean_ci_table_par$x_minCI <- paste(mean_ci_table_par$x_minBT, " [", mean_ci_table_par$x_minBT_lci, ", ", mean_ci_table_par$x_minBT_uci,"]", sep="")
mean_ci_table_par$x_maxCI <- paste(mean_ci_table_par$x_maxBT, " [", mean_ci_table_par$x_maxBT_lci, ", ", mean_ci_table_par$x_maxBT_uci,"]", sep="")
mean_ci_table_par$max_RGRCI <- paste(mean_ci_table_par$max_RGRBT, " [", mean_ci_table_par$max_RGR_lci, ", ", mean_ci_table_par$max_RGR_uci,"]", sep="")
mean_ci_table_par$areaCI <- paste(mean_ci_table_par$area, " [", mean_ci_table_par$area_lci, ", ", mean_ci_table_par$area_uci,"]", sep="")
mean_ci_table_par <- mean_ci_table_par[,c("species", "ToptCI", "B50CI", "B50_lowCI", "B50_highCI", "breadthCI", "x_minCI", "x_maxCI", "max_RGRCI", "areaCI")]
write.csv(mean_ci_table_par, "Analysis-output/TPC/Mean-dfs/mean_ci_sxp_af_table_par.csv")


tidy_perf_ver <- performr::perform_df(model_fits_ver, species_order = c(1:length(unique(ver_dat$familyI)))) 
tidy_perf_ver %<>% ungroup() %>% mutate(idx = 1:n()) %>% group_by(idx) %>% do({ mutate(.data = ., max_RGR = performance_mu(xs = .$maxima, .$shape1, .$shape2, .$stretch, .$x_min, .$x_max))})
#ssq_df_ver <- bayes_p(stan_df = tidy_perf_ver, raw_df = ver_dat, raw_treatment = "dayTempc", raw_response = "RGRs", raw_group = "familyI")
tidy_perf_ver %<>% ungroup() %>% bind_cols(optimum_breadth(par_df = tidy_perf_ver,  prop_max = 0.5,  n_grid = 100))
tidy_perf_ver$B50_low <- tidy_perf_ver$opt_breadth_low + meansTot$dayTemp[which(meansTot$species=="ver")] 
tidy_perf_ver$B50_high <- tidy_perf_ver$opt_breadth_high + meansTot$dayTemp[which(meansTot$species=="ver")]  
tidy_perf_ver$B50 <- tidy_perf_ver$B50_high - tidy_perf_ver$B50_low 
draws_ver <- rstan::extract(model_fits_ver); ndraws_ver <- length(draws_ver$lp__)
x_seq_ver = seq(min(draws_ver$x_min),max(draws_ver$x_max),length.out = 100)
poly_draws_ver <- sample(1:ndraws_ver, 100)
creds_ver <-  predict_interval(x = x_seq_ver, spp = species, par_df = tidy_perf_ver, x_draws = poly_draws_ver, p = c(0.95, 0.5))
tidy_perf_ver$x_minBT <- tidy_perf_ver$x_min + meansTot$dayTemp[which(meansTot$species=="ver")] 
tidy_perf_ver$x_maxBT <- tidy_perf_ver$x_max + meansTot$dayTemp[which(meansTot$species=="ver")] 
tidy_perf_ver$maximaBT <- tidy_perf_ver$maxima + meansTot$dayTemp[which(meansTot$species=="ver")] 
tidy_perf_ver$breadthBT <- tidy_perf_ver$x_maxBT-tidy_perf_ver$x_minBT
tidy_perf_ver$max_RGRBT <- tidy_perf_ver$max_RGR * meansTot$RGR[which(meansTot$species=="ver")] 
write.csv(tidy_perf_ver, "Analysis-output/TPC/Models/tidy_perf_sxp_af_ver.csv")
creds_ver$x <- creds_ver$x + meansTot$dayTemp[which(meansTot$species=="ver")] ; creds_ver$mu <- creds_ver$mu * meansTot$RGR[which(meansTot$species=="ver")] 
creds_ver$upper <- creds_ver$upper * meansTot$RGR[which(meansTot$species=="ver")] ; creds_ver$lower <- creds_ver$lower * meansTot$RGR[which(meansTot$species=="ver")] 
write.csv(creds_ver, "Analysis-output/TPC/Creds/creds_sxp_af_ver.csv")
mean_df_ver <- tidy_perf_ver %>% group_by(species) %>% summarise_if(is.numeric, .funs = c(mean))
mean_df_ver$speciesn <- as.numeric(mean_df_ver$species); mean_df_ver <- mean_df_ver[order(mean_df_ver$speciesn),]
write.csv(mean_df_ver, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ver.csv")
mean_df_ci_ver <- tidy_perf_ver %>%  group_by(species) %>% summarise(maximaBT_lci = quantile(maximaBT, probs=c(0.025)), maximaBT_uci = quantile(maximaBT, probs=c(0.975)), 
                                                                     B50_lci = quantile(B50, probs=c(0.025)), B50_uci = quantile(B50, probs=c(0.975)), 
                                                                     B50_low_lci = quantile(B50_low, probs=c(0.025)), B50_low_uci = quantile(B50_low, probs=c(0.975)), 
                                                                     B50_high_lci = quantile(B50_high, probs=c(0.025)), B50_high_uci = quantile(B50_high, probs=c(0.975)), 
                                                                     breadthBT_lci = quantile(breadthBT, probs=c(0.025)), breadthBT_uci = quantile(breadthBT, probs=c(0.975)), 
                                                                     x_minBT_lci = quantile(x_minBT, probs=c(0.025)), x_minBT_uci = quantile(x_minBT, probs=c(0.975)), 
                                                                     x_maxBT_lci = quantile(x_maxBT, probs=c(0.025)), x_maxBT_uci = quantile(x_maxBT, probs=c(0.975)), 
                                                                     max_RGR_lci = quantile(max_RGRBT, probs=c(0.025)), max_RGR_uci = quantile(max_RGRBT, probs=c(0.975)), 
                                                                     area_lci = quantile(area, probs=c(0.025)), area_uci = quantile(area, probs=c(0.975)))
mean_df_ci_ver$speciesn <- as.numeric(mean_df_ci_ver$species); mean_df_ci_ver <- mean_df_ci_ver[order(mean_df_ci_ver$speciesn),]
write.csv(mean_df_ci_ver, "Analysis-output/TPC/Mean-dfs/mean_df_sxp_af_ci_ver.csv")
mean_ci_table_ver <- mean_df_ver %>% right_join(mean_df_ci_ver, by="species")
mean_ci_table_ver <- as.data.frame(mean_ci_table_ver)
mean_ci_table_ver <- mean_ci_table_ver[,c("species", "speciesn.x", 
                                          "maximaBT", "maximaBT_lci", "maximaBT_uci", 
                                          "B50", "B50_lci", "B50_uci", 
                                          "B50_low", "B50_low_lci", "B50_low_uci", 
                                          "B50_high", "B50_high_lci", "B50_high_uci", 
                                          "breadthBT", "breadthBT_lci", "breadthBT_uci", 
                                          "x_minBT", "x_minBT_lci", "x_minBT_uci", 
                                          "x_maxBT", "x_maxBT_lci", "x_maxBT_uci", 
                                          "max_RGRBT", "max_RGR_lci", "max_RGR_uci", 
                                          "area", "area_lci", "area_uci")]
mean_ci_table_ver <- mean_ci_table_ver %>% mutate_at(vars(-c("species","speciesn.x")), funs(round(., 2)))
mean_ci_table_ver$ToptCI <- paste(mean_ci_table_ver$maximaBT, " [", mean_ci_table_ver$maximaBT_lci, ", ", mean_ci_table_ver$maximaBT_uci,"]", sep="")
mean_ci_table_ver$B50CI <- paste(mean_ci_table_ver$B50, " [", mean_ci_table_ver$B50_lci, ", ", mean_ci_table_ver$B50_uci,"]", sep="")
mean_ci_table_ver$B50_lowCI <- paste(mean_ci_table_ver$B50_low, " [", mean_ci_table_ver$B50_low_lci, ", ", mean_ci_table_ver$B50_low_uci,"]", sep="")
mean_ci_table_ver$B50_highCI <- paste(mean_ci_table_ver$B50_high, " [", mean_ci_table_ver$B50_high_lci, ", ", mean_ci_table_ver$B50_high_uci,"]", sep="")
mean_ci_table_ver$breadthCI <- paste(mean_ci_table_ver$breadthBT, " [", mean_ci_table_ver$breadthBT_lci, ", ", mean_ci_table_ver$breadthBT_uci,"]", sep="")
mean_ci_table_ver$x_minCI <- paste(mean_ci_table_ver$x_minBT, " [", mean_ci_table_ver$x_minBT_lci, ", ", mean_ci_table_ver$x_minBT_uci,"]", sep="")
mean_ci_table_ver$x_maxCI <- paste(mean_ci_table_ver$x_maxBT, " [", mean_ci_table_ver$x_maxBT_lci, ", ", mean_ci_table_ver$x_maxBT_uci,"]", sep="")
mean_ci_table_ver$max_RGRCI <- paste(mean_ci_table_ver$max_RGRBT, " [", mean_ci_table_ver$max_RGR_lci, ", ", mean_ci_table_ver$max_RGR_uci,"]", sep="")
mean_ci_table_ver$areaCI <- paste(mean_ci_table_ver$area, " [", mean_ci_table_ver$area_lci, ", ", mean_ci_table_ver$area_uci,"]", sep="")
mean_ci_table_ver <- mean_ci_table_ver[,c("species", "ToptCI", "B50CI", "B50_lowCI", "B50_highCI", "breadthCI", "x_minCI", "x_maxCI", "max_RGRCI", "areaCI")]
write.csv(mean_ci_table_ver, "Analysis-output/TPC/Mean-dfs/mean_ci_sxp_af_table_ver.csv")






# compute bayesian p value for each species and export dataframe
# bps <- data.frame(group=c(1:10),
#                   sp=c("bic", "car", "eas", "fil", "flo", "gut", "lac", "nor", "par", "ver"),
#                   ps=rep(NA,10))
# bps$ps[1] <- ssq_df_bic %>% summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# bps$ps[2] <- ssq_df_car %>% summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# bps$ps[3] <- ssq_df_eas %>% summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# bps$ps[4] <- ssq_df_fil %>% summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# bps$ps[5] <- ssq_df_flo %>% summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# bps$ps[6] <- ssq_df_gut %>% summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# bps$ps[7] <- ssq_df_lac %>% summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# bps$ps[8] <- ssq_df_nor %>% summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# bps$ps[9] <- ssq_df_par %>% summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# bps$ps[10] <- ssq_df_ver %>% summarise(b_pval=mean(ssq_obs > ssq_pseudo))
# 
# bps$ps <- as.numeric(bps$ps)
# write.csv(bps, file="Analysis-output/TPC/Bayesian-p-values_sxp-af.csv")





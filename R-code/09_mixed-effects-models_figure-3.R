#### PROJECT: Mimulus niche breadth partitioning
#### PURPOSE: Test hypotheses about how family-level breadth contributes to 
####          population-level breadth and create figure 3
#### AUTHOR/LAST UPDATE: EMC/2022-08-19


## load packages ##
library("lme4")
library("car")
library("lsmeans")
library("tidyverse")
library("gridExtra")
library("lsmeans")
library("lmerTest")
library("MuMIn")
library("lattice")

## read data ##
tpc <- read.csv("Analysis-output/family-mean-estimates-sxp-af.csv") 
var <- read.csv("Processed-data/breadth-and-var-by-species.csv")
# note: "breadth-var-by-species.csv" was created manually in excel

## calculate breadth ##
tpc$t.breadth2 <- tpc$t.max2 - tpc$t.min2

## join data ##
var.data <- select(var, species, pop.breadth)
tpc <- left_join(tpc, var.data, by="species")


#############################################
# population breadth ~ family level breadth #
#############################################

## run model and get summary output ##
pop.breadth.lmer <- lmer(pop.breadth ~ t.breadth2 + (1 | pair), data=tpc)
summary(pop.breadth.lmer)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: pop.breadth ~ t.breadth2 + (1 | pair)
# Data: tpc
# 
# REML criterion at convergence: 725
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -3.02492 -0.25744  0.00588  0.44081  1.84974 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# pair     (Intercept) 23.22    4.819   
# Residual              1.11    1.053   
# Number of obs: 235, groups:  pair, 5
# 
# Fixed effects:
#                 Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)     28.47168    2.27824   4.90365  12.497 6.62e-05 ***
#   t.breadth2    0.06223    0.03367 231.72130   1.849   0.0658 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr)
# t.breadth2 -0.323

## calculate p-value for two tailed test ##
m1.pval <- summary(pop.breadth.lmer)$coefficients[,"Pr(>|t|)"][2]
m1.pval
# 0.06579085 

## calculate marginal and conditional r2 ##
r.squaredGLMM(pop.breadth.lmer)
#           R2m       R2c
# [1,] 0.004975914 0.9546151



#############################################
## population breadth ~ variation in T-opt ## 
#############################################

## run model and get summary output ##
pop.breadth.opt.lmer <- lmer(pop.breadth ~ opt.var + (1 | pair), data=var)
summary(pop.breadth.opt.lmer)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: pop.breadth ~ opt.var + (1 | pair)
# Data: var
# 
# REML criterion at convergence: 46.4
# 
# Scaled residuals: 
# Min       1Q   Median       3Q      Max 
# -0.84821 -0.62645 -0.02471  0.51413  1.12570 
# 
#  #Random effects:
# Groups   Name        Variance Std.Dev.
# pair     (Intercept) 22.131   4.704   
# Residual              1.809   1.345   
# Number of obs: 10, groups:  pair, 5
#
# Fixed effects:
# Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)  25.5181     2.7668  7.5429   9.223 2.27e-05 ***
# opt.var       0.5084     0.2130  4.3721   2.386     0.07 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr)
# opt.var -0.631

## calculate p-value for two tailed test ##
m2.pval <- summary(pop.breadth.opt.lmer)$coefficients[,"Pr(>|t|)"][2]
m2.pval
# 0.07000248 

## calculate marginal and conditional r2 ##
r.squaredGLMM(pop.breadth.opt.lmer)
# R2m      R2c
# [1,] 0.08662951 0.9309972


#################################################
## population breadth ~ variation in T-breadth ##
#################################################

## run model and get summary output ##
pop.breadth.var.lmer <- lmer(pop.breadth ~ breadth.var + (1 | pair), data=var)
summary(pop.breadth.var.lmer)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: pop.breadth ~ breadth.var + (1 | pair)
# Data: var
#
# REML criterion at convergence: 49.3
# 
# Scaled residuals: 
# Min       1Q   Median       3Q      Max 
# -1.21460 -0.44997  0.06626  0.46405  1.08802 
#
# Random effects:
# Groups   Name        Variance Std.Dev.
# pair     (Intercept) 22.284   4.721   
# Residual              4.299   2.073   
# Number of obs: 10, groups:  pair, 5
# 
# Fixed effects:
# Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)  27.7796     4.0259  7.7109   6.900 0.000148 ***
# breadth.var   0.2747     0.4854  4.7384   0.566 0.597141    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr)
# breadth.var -0.836

## calculate p-value for two tailed test ##
m3.pval <- summary(pop.breadth.var.lmer)$coefficients[,"Pr(>|t|)"][2]
m3.pval
# 0.5971406  

## calculate marginal and conditional r2 ##
r.squaredGLMM(pop.breadth.var.lmer)
# R2m       R2c
# [1,] 0.01425338 0.8405845


############################################################
## Figure 3: Results that correspond to hypotheses/models ##
############################################################

se1 <- 0.3

h1 <- ggplot(tpc, aes(y=pop.breadth, x=t.breadth2)) + 
  geom_point(aes(x=t.breadth2, y=pop.breadth, fill=species), shape=21, size=4, alpha=0.8, color="gray25") + 
  ggtitle("A") + theme_bw() + 
  xlab(expression(family ~italic(T[breadth]))) + ylab("  ") +
  #geom_point(data=var, aes(y=pop.breadth, x=fam.breadth), shape=3, size=4, stroke=1, color="black") + 
  scale_fill_manual(name="Species", values=c("mediumpurple3","thistle1","goldenrod2","lightgoldenrod1","green4","palegreen1","royalblue3","lightskyblue1","violetred3","pink1"), 
                     breaks=c("par","car","ver","eas", "bic", "fil","flo", "nor","lac", "gut"), 
                     labels=c("M. parishii","M. cardinalis","M. verbanaceous","M. eastwoodiae", "M. bicolor", "M. filicaulis","M. floribundus", "M. norrisii","M. laciniatus", "M. guttatus")) + 
  theme(legend.position = "bottom", text=element_text(size=23),  legend.text=element_text(size=20, face="italic"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.direction = "vertical") + 
  guides(fill=guide_legend(override.aes=list(size=5))) + geom_abline(slope=1 , intercept=0, color="gray", size=2) + 
  geom_point(data=var, aes(x=fam.breadth, y=pop.breadth), fill="black", shape=23, size=3) + 
  geom_errorbar(data=var, aes(x=fam.breadth, y=pop.breadth, xmin=fam.breadth.cmin, xmax=fam.breadth.cmax, width=se1, alpha=0.75), show.legend=FALSE) + 
  geom_errorbar(data=var, aes(x=fam.breadth, y=pop.breadth, ymin=pop.breadth.cmin, ymax=pop.breadth.cmax, width=se1, alpha=0.75), show.legend=FALSE) + 
  geom_line(aes(y = predict(pop.breadth.lmer), color=pair), size=1, show.legend=FALSE) + 
  scale_color_manual(values=c("mediumpurple4", "darkgoldenrod","darkgreen","royalblue4", "violetred4"), 
                     breaks=c("car-par","ver-eas", "bic-fil", "flo-nor","gut-lac")) + guides(fill=guide_legend(ncol=2))


h2 <- ggplot(var, aes(y=pop.breadth, x=opt.var)) + 
  theme_bw() + 
  theme(legend.position = "none", text=element_text(size=23), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab(expression(family ~italic(T[opt])~ variance)) + 
  ylab(" ") + ggtitle("B") + 
  scale_fill_manual(name="Species", values=c("mediumpurple3","thistle1","goldenrod2","lightgoldenrod1","green4","palegreen1","royalblue3","lightskyblue1","violetred3","pink1"), 
                    breaks=c("par","car","ver","eas", "bic", "fil","flo", "nor","lac", "gut"), 
                    labels=c("M. parishii","M. cardinalis","M. verbanaceous","M. eastwoodiae", "M. bicolor", "M. filicaulis","M. floribundus", "M. norrisii","M. laciniatus", "M. guttatus")) + 
  geom_errorbar(data=var, aes(x=opt.var, y=pop.breadth, xmin=opt.var.cmin, xmax=opt.var.cmax, width=se1, alpha=0.75)) + 
  geom_errorbar(data=var, aes(x=opt.var, y=pop.breadth, ymin=pop.breadth.cmin, ymax=pop.breadth.cmax, width=se1, alpha=0.75)) + 
  geom_point(aes(x=opt.var, y=pop.breadth, fill=species), shape=21, size=5, stroke=1, color="gray25") + 
  geom_line(aes(y = predict(pop.breadth.opt.lmer), color=pair), size=1) + 
  scale_color_manual(values=c("mediumpurple4", "darkgoldenrod","darkgreen","royalblue4", "violetred4"), 
                     breaks=c("car-par","ver-eas", "bic-fil", "flo-nor","gut-lac"))
#  scale_color_manual(values=c("mediumpurple3", "goldenrod2","green4","royalblue3", "violetred3"), breaks=c("car-par","ver-eas", "bic-fil", "flo-nor","gut-lac"))

h3_lines = data.frame(x = c(min(var$breadth.var), max(var$breadth.var)),
                      y = c(b3 + m3 * min(var$breadth.var), b3 + m3 * max(var$breadth.var)))
 
h3 <- ggplot(var, aes(x=breadth.var, y=pop.breadth)) + theme_bw() + 
  theme(text=element_text(size=23), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab(expression(family ~italic(T[breadth])~ variance)) + 
  ylab(" ") + ggtitle("C") + 
    scale_fill_manual(name="Species", values=c("mediumpurple3","thistle1","goldenrod2","lightgoldenrod1","green4","palegreen1","royalblue3","lightskyblue1","violetred3","pink1"), 
                      breaks=c("par","car","ver","eas", "bic", "fil","flo", "nor","lac", "gut"), 
                      labels=c("M. parishii","M. cardinalis","M. verbanaceous","M. eastwoodiae", "M. bicolor", "M. filicaulis","M. floribundus", "M. norrisii","M. laciniatus", "M. guttatus")) + 
  geom_errorbar(data=var, aes(x=breadth.var, y=pop.breadth, xmin=breadth.var.cmin, xmax=breadth.var.cmax, width=se1, alpha=0.75)) + 
  geom_errorbar(data=var, aes(x=breadth.var, y=pop.breadth, ymin=pop.breadth.cmin, ymax=pop.breadth.cmax, width=se1, alpha=0.75)) + 
  geom_point(aes(x=breadth.var, y=pop.breadth, fill=species), shape=21, size=5, stroke=1, color="gray25") + 
  geom_line(aes(y = predict(pop.breadth.var.lmer), color=pair), size=1, linetype="dashed") + 
  scale_color_manual(values=c("mediumpurple4", "darkgoldenrod","darkgreen","royalblue4", "violetred4"), 
                     breaks=c("car-par","ver-eas", "bic-fil", "flo-nor","gut-lac")) + xlim(NA, 17)
  
## extract legend ##
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(h1)

text = expression(paste(population ~italic(T[breadth])))
lab <- ggplot() + annotate("text", x = 4, y = 25, size=10, angle=90, label=text) + theme_void()

text = ""
nolab <- ggplot() + annotate("text", x = 4, y = 25, size=23, label = text) + 
  theme_void()

## assemble and export multipanel figure ##
## note: this script is used to create figure 3 for the manuscript as of 2022-05-01
pdf("Figures/3_Figure-3_hypothesis-testing-results_vertical.pdf",  height=18.5, width=6.5)
grid.arrange(arrangeGrob(nolab, h1+theme(legend.position='hidden'), lab, h2+theme(legend.position='hidden'),
                        nolab, h3+theme(legend.position='hidden'), nolab, mylegend, ncol=2, nrow=4,
                         heights=c(5.5,5.5,5.5,2), widths=c(0.4,6), padding = 0.1))
dev.off()


#####
## testing w/o pair
#####



## population-level breadth ~ family-level breadth ##

pop.breadth.lm <- lm(pop.breadth ~ t.breadth2, data=tpc)
summary(pop.breadth.lm)
# Call:
# lm(formula = pop.breadth ~ t.breadth2, data = tpc)
#
# Residuals:
#   Min     1Q Median     3Q    Max 
# -5.176 -1.435  0.051  1.329  4.975 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 15.59287    0.57132   27.29   <2e-16 ***
#   t.breadth2   0.66339    0.02323   28.56   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 1.991 on 233 degrees of freedom
# Multiple R-squared:  0.7778,	Adjusted R-squared:  0.7769 
# F-statistic: 815.7 on 1 and 233 DF,  p-value: < 2.2e-16



## Population-level breadth ~ Topt variation ##

pop.breadth.opt.lm <- lm(pop.breadth ~ opt.var, data=var)
summary(pop.breadth.opt.lm)

# Call:
# lm(formula = pop.breadth ~ opt.var, data = var)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -5.429 -2.828 -1.899  3.540  6.871 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  22.3009     4.5900   4.859  0.00126 **
#   opt.var       0.9009     0.5298   1.701  0.12744   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.711 on 8 degrees of freedom
# Multiple R-squared:  0.2655,	Adjusted R-squared:  0.1737 
# F-statistic: 2.892 on 1 and 8 DF,  p-value: 0.1274



## Population-level breadth ~ Tbreadth variation ##

pop.breadth.var.lm <- lm(pop.breadth ~ breadth.var, data=var)
summary(pop.breadth.var.lm)

# Call:
# lm(formula = pop.breadth ~ breadth.var, data = var)
#
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.9414 -3.1154 -0.1153  2.4801  5.3846 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  18.9406     4.3471   4.357  0.00242 **
#   breadth.var   1.5498     0.5992   2.587  0.03228 * 
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 4.057 on 8 degrees of freedom
# Multiple R-squared:  0.4554,	Adjusted R-squared:  0.3874 
# F-statistic:  6.69 on 1 and 8 DF,  p-value: 0.03228


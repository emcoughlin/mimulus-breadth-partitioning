#### PROJECT: Mimulus niche breadth partitioning
#### PURPOSE: Illustrate hypotheses about how family-level breadth contributes to population-level breadth
#### AUTHOR/LAST UPDATED: RCW & EMC/2021-07-29


# load necessary packages
library(devtools)
library(performr)
library(ggplot2)
library(gridExtra)


##################################################   
# Families of widespread species are generalists #
##################################################

# build a widespread species tpc 
shp1 <- 2.1 # When shape1 is larger than shape2, the curve will skew right
shp2 <- 2.5 # When shape2 is larger than shape1, the curve will skew left
lower <- 4 # set lower limit of curve
upper <- 51 # set upper limit of curve
xs <- seq(lower, upper, length.out=100) # temperatures for the tpc
tpc <- performance_mu(xs = xs, # Temperatures for the tpc
                      shape1 = shp1, 
                      shape2 = shp2, 
                      stretch = 0.51, # Dictates the maximum expected value of the response trait
                      x_min = lower, # Location along the environmental axis left of the optimum where the response trait falls to 0
                      x_max = upper) # Location along the environmental axis right of the optimum where the response trait falls to 0
tpcdf_sp1 <- data.frame(x= xs, y=tpc)

# build 3 family level tpcs
lower <- c(6, 7, 8)
upper <- c(47, 48, 49)
stretch <- c(0.49, 0.49, 0.49)

xs <- seq(lower[1], upper[1], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1, shape2 = shp2, 
                      stretch = stretch[1], x_min = lower[1], x_max = upper[1])
tpcdf_fam1 <- data.frame(x= xs, y=tpc)

xs <- seq(lower[2], upper[2], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1, shape2 = shp2, 
                      stretch = stretch[2], x_min = lower[2], x_max = upper[2])
tpcdf_fam2 <- data.frame(x= xs, y=tpc)

xs <- seq(lower[3], upper[3], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1, shape2 = shp2, 
                      stretch = stretch[3], x_min = lower[3], x_max = upper[3])
tpcdf_fam3 <- data.frame(x= xs, y=tpc)

# build a restricted species tpc 
shp1 <- 2.1 # When shape1 is larger than shape2, the curve will skew right
shp2 <- 2.5 # When shape2 is larger than shape1, the curve will skew left
lower <- 14 # set lower limit of curve
upper <- 41 # set upper limit of curve
xs <- seq(lower, upper, length.out=100) # temperatures for the tpc
tpc <- performance_mu(xs = xs, # Temperatures for the tpc
                      shape1 = shp1, 
                      shape2 = shp2, 
                      stretch = 0.51, # Dictates the maximum expected value of the response trait
                      x_min = lower, # Location along the environmental axis left of the optimum where the response trait falls to 0
                      x_max = upper) # Location along the environmental axis right of the optimum where the response trait falls to 0
tpcdf_sp2 <- data.frame(x= xs, y=tpc)

# build 3 family level tpcs
lower <- c(16, 17, 18)
upper <- c(37, 38, 39)
stretch <- c(0.49, 0.49, 0.49)

xs <- seq(lower[1], upper[1], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1, shape2 = shp2, 
                      stretch = stretch[1], x_min = lower[1], x_max = upper[1])
tpcdf_fam4 <- data.frame(x= xs, y=tpc)

xs <- seq(lower[2], upper[2], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1, shape2 = shp2, 
                      stretch = stretch[2], x_min = lower[2], x_max = upper[2])
tpcdf_fam5 <- data.frame(x= xs, y=tpc)

xs <- seq(lower[3], upper[3], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1, shape2 = shp2, 
                      stretch = stretch[3], x_min = lower[3], x_max = upper[3])
tpcdf_fam6 <- data.frame(x= xs, y=tpc)

# now plot species and family curves together

tpc_a <- ggplot() +
  labs(x=expression(paste("Environment")), y="Performance", title="A") + 
  scale_x_continuous(limits=c(0,55), expand = c(0, 0), breaks = NULL) +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0), breaks = NULL) +
  # narrow curves
  geom_line(data=tpcdf_fam4, aes(x, y), color="deepskyblue2", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam5, aes(x, y), color="deepskyblue3", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam6, aes(x, y), color="deepskyblue4", alpha=0.5, size=1) +
  # broad curves
  geom_line(data=tpcdf_fam1, aes(x, y), color="red2", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam2, aes(x, y), color="red3", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam3, aes(x, y), color="red4", alpha=0.5, size=1) +
  # broad breadth
  annotate(geom="segment", x=13.85, y=.4, xend=41.9, yend=.4, color="red", size=1.5, arrow=arrow(ends="both", length=unit(.3,"cm"), angle=25), lineend="round", linejoin="mitre") + 
  # narrow breadth
  annotate(geom="segment", x=19.8, y=.36, xend=35.6, yend=.36, color="deepskyblue", size=1.5, arrow=arrow(ends="both", length=unit(.3,"cm"), angle=25), lineend="round", linejoin="mitre") + 
  # general settings
  guides(fill=FALSE, color=FALSE) + theme_classic() + 
  theme(panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA), text=element_text(size=30), 
        strip.background=element_rect(colour=NA,fill=NA), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position = "bottom") 


##########################################################   
# Widespread species have among-family variation in Topt #
##########################################################

# build 3 family level tpcs
shp1 <- c(5, 3.3, 2.8) 
shp2 <- c(9, 7, 9) 
lower <- c(1,3,4)
upper <- c(53,55,57)
stretch <- c(0.2625, 0.328, 0.317)

xs <- seq(lower[1], upper[1], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1[1], shape2 = shp2[1], stretch = stretch[1], x_min = lower[1], x_max = upper[1])
tpcdf_fam1 <- data.frame(x= xs, y=tpc)
xs <- seq(lower[2], upper[2], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1[2], shape2 = shp2[2], stretch = stretch[2], x_min = lower[2], x_max = upper[2])
tpcdf_fam2 <- data.frame(x= xs, y=tpc)
xs <- seq(lower[3], upper[3], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1[3], shape2 = shp2[3], stretch = stretch[3], x_min = lower[3], x_max = upper[3])
tpcdf_fam3 <- data.frame(x= xs, y=tpc)

# now plot species and family curves together
tpc_b <- ggplot() +
  labs(x=expression(paste("Environment")), y="", title="B") + 
  scale_x_continuous(limits=c(0,55), expand = c(0, 0), breaks = NULL) +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0), breaks = NULL) +
  # narrow curves
  geom_line(data=tpcdf_fam4, aes(x, y), color="deepskyblue2", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam5, aes(x, y), color="deepskyblue3", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam6, aes(x, y), color="deepskyblue4", alpha=0.5, size=1) +
  # broad curves
  geom_line(data=tpcdf_fam1, aes(x, y), color="red2", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam2, aes(x, y), color="red3", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam3, aes(x, y), color="red4", alpha=0.5, size=1) +
  # wide breadth
  annotate(geom="segment", x=15, y=.4, xend=40, yend=.4, color="red", size=1.5, arrow=arrow(ends="both", length=unit(.3,"cm"), angle=25), lineend="round", linejoin="mitre") + 
  # narrow breadth
  annotate(geom="segment", x=19.8, y=.36, xend=35.6, yend=.36, color="deepskyblue", size=1.5, arrow=arrow(ends="both", length=unit(.3,"cm"), angle=25), lineend="round", linejoin="mitre") + 
  # general settings
  guides(fill=FALSE, color=FALSE) + theme_classic() + 
  theme(panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA), text=element_text(size=30), 
        strip.background=element_rect(colour=NA,fill=NA), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())        


#########################################################
# Widespread species have among-family variation in b50 #
#########################################################

# build 3 family level tpcs
shp1 <- c(2.1, 3, 8) 
shp2 <- c(2.5, 5, 200) 
lower <- c(5.5,5,10)
upper <- c(49.5,50,45)
stretch <- c(0.49, 0.365, 0.14)

xs <- seq(lower[1], upper[1], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1[1], shape2 = shp2[1], 
                      stretch = stretch[1], x_min = lower[1], x_max = upper[1])
tpcdf_fam1 <- data.frame(x= xs, y=tpc)

xs <- seq(lower[2], upper[2], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1[2], shape2 = shp2[2], 
                      stretch = stretch[2], x_min = lower[2], x_max = upper[2])
tpcdf_fam2 <- data.frame(x= xs, y=tpc)

xs <- seq(lower[3], upper[3], length.out=100)
tpc <- performance_mu(xs = xs, shape1 = shp1[3], shape2 = shp2[3], 
                      stretch = stretch[3], x_min = lower[3], x_max = upper[3])
tpcdf_fam3 <- data.frame(x= xs, y=tpc)


# now plot species and family curves together
tpc_c <- ggplot() +
  labs(x=expression(paste("Environment")), y="", title="C") + 
  scale_x_continuous(limits=c(0,55), expand = c(0, 0), breaks = NULL) +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0), breaks = NULL) +
  # narrow curves
  geom_line(data=tpcdf_fam4, aes(x, y), color="deepskyblue2", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam5, aes(x, y), color="deepskyblue3", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam6, aes(x, y), color="deepskyblue4", alpha=0.5, size=1) +
  # broad curves
  geom_line(data=tpcdf_fam1, aes(x, y), color="red2", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam2, aes(x, y), color="red3", alpha=0.5, size=1) +
  geom_line(data=tpcdf_fam3, aes(x, y), color="red4", alpha=0.5, size=1) +
  # broad breadth
  annotate(geom="segment", x=13.85, y=.4, xend=41.9, yend=.4, color="red", size=1.5, arrow=arrow(ends="both", length=unit(.3,"cm"), angle=25), lineend="round", linejoin="mitre") + 
  # narrow breadth
  annotate(geom="segment", x=19.8, y=.36, xend=35.6, yend=.36, color="deepskyblue", size=1.5, arrow=arrow(ends="both", length=unit(.3,"cm"), angle=25), lineend="round", linejoin="mitre") + 
  # general settings
  guides(fill=FALSE, color=FALSE) + theme_classic() + 
  theme(panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA), text=element_text(size=30), 
        strip.background=element_rect(colour=NA,fill=NA), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())
        #panel.border = element_rect(colour = "black", fill=NA, size=1.5))


################# 
# Create legend #
#################   

leg <- ggplot() +
  annotate(geom="text", x=c(1,2), y=0.5, xend=35.6, yend=.36, color="black", label=c("Narrow population-level niche", "Broad population-level niche"), size=8) + 
  annotate(geom="segment", x=c(0.5, 1.51), y=0.5, xend=c(0.6, 1.61), yend=0.5, color=c("deepskyblue", "red"), size=1.5, lineend="round", linejoin="mitre") + 
  scale_x_continuous(limits=c(0,3), expand = c(0, 0), breaks = NULL) +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0), breaks = NULL) +
  theme(panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA), text=element_text(size=0), 
        strip.background=element_rect(colour=NA,fill=NA), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())


#################
# Export figure #
#################

# note that this is script to make Fig. 1 in manuscript as of 20210729

pdf("Figures/1_Figure-1_hypotheses.pdf", height=7, width=18)
lay <- rbind(c(1,2,3),
             c(4,4,4))
figure <- grid.arrange(tpc_a, tpc_b, tpc_c, leg, ncol=3, nrow=2,
                       layout_matrix = lay, heights=c(7,1))
dev.off()


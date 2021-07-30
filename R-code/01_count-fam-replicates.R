#### PROJECT: Mimulus niche breadth partitioning
#### PURPOSE: Counting families and replicates per species per treatment
#### AUTHOR/LAST UPDATE: EMC/2021-07-29


library("tidyverse")

# import #
temp <- read.csv("Raw-data/thermal_data.csv")

# count individuals per species #
indiv <- temp %>% select(species, family) %>% group_by(species) %>% count()

indiv
# bic 326
# car 691
# eas 1311
# fil 179
# flo 242
# gut 175
# lac 170
# nor 268
# par 1327
# ver 762

# count family replicates per species per temperature #
reps <- temp %>% select(species, family, dayTemp) %>% group_by(species, dayTemp, family) %>% count()

summary(reps)
#     species             dayTemp          family            n        
# Length:1843        Min.   :15.00   Min.   :  1.0   Min.   :1.000  
# Class :character   1st Qu.:25.00   1st Qu.: 67.0   1st Qu.:2.000  
# Mode  :character   Median :35.00   Median :127.0   Median :3.000  
# Mean   :32.57   Mean   :142.2   Mean   :2.958  
# 3rd Qu.:45.00   3rd Qu.:234.0   3rd Qu.:4.000  
# Max.   :50.00   Max.   :306.0   Max.   :4.000  

# 1-4 replicates per species per temperature #

# export #
write.csv(indiv, "Processed-data/individuals-per-species.csv")
write.csv(reps, "Processed-data/replicates-per-species-per-temp.csv")

#Script for building models w/env predictors

remove(list = ls())

library(tidyverse)
library(gridExtra)
library(plotrix)
library(lme4)
library(MuMIn)
library(glmmTMB)
library(DHARMa)

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)
msat <- read.csv("output/msatloci_assembled.csv", stringsAsFactors = FALSE)
mtdna_env <- read.csv("output/mtdna_env.csv", stringsAsFactors = FALSE)
msat_env <- read.csv("output/msat_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)

#merge dataframes
mtdna <- merge(mtdna, mtdna_env[, c('X', 'sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstrange', 
                                    'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                    'chloroA.BO_chlomax', 'chloroA.BO_chlomin')], all.x = TRUE)
msat <- cbind(msat[, -1], msat_env[, c('sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstrange', 
                                       'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                       'chloroA.BO_chlomax', 'chloroA.BO_chlomin')]) #merge not working for some reason, cbind bc in same order
mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                  'Centroid')], all.x = TRUE)
msat <- merge(msat, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                'Centroid')], all.x = TRUE)

#clean up dataframes
#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)

####################################################################################################################

######## Building models for mtdna Hd ########
#following best model structure from model.R script

######## Cleaning mtdna Hd dataset ########
#subset mtdna to remove He = NA columns
mtdna_small_He <- subset(mtdna_small, mtdna_small$He != "NA")

##### scaling bp & making sure numeric ####
mtdna_small_He$bp <- as.numeric(mtdna_small_He$bp)
mtdna_small_He$bp_scale <- scale(as.numeric(mtdna_small_He$bp))

#### calculating position in range ####
mtdna_small_He$Centroid <- as.numeric(mtdna_small_He$Centroid)
mtdna_small_He$Half_RangeSize <- as.numeric(mtdna_small_He$Half_RangeSize)

mtdna_small_He$range_position <- NA #create column to fill in

for (i in 1:nrow(mtdna_small_He)) { #calculate distance from range center as percentage (0-1 both sides of centroid, 1 = all the way at range edge)
  mtdna_small_He$range_position[i] <- abs((mtdna_small_He$lat[i] - mtdna_small_He$Centroid[i])/mtdna_small_He$Half_RangeSize[i])
}

#check those >1 and round to 1 (slight discrepancy btwn aquamaps and reality)
mtdna_check <- subset(mtdna_small_He, mtdna_small_He$range_position > 1) #often right at aquamaps limit, round to 1 and keep
mtdna_small_He$range_position[mtdna_small_He$range_position > 1] <- 1

#subset to only those with range_position
mtdna_small_He <- subset(mtdna_small_He, range_position != "NA")

#calculate success and failures
mtdna_small_He$success <- round(mtdna_small_He$He*mtdna_small_He$n)
mtdna_small_He$failure<- round((1 - mtdna_small_He$He)*mtdna_small_He$n)

######## Build null model (only includes nuisance variables) ########
#variables to include: bp, position in spp range, (1|MarkerName), (1|Family/Genus/spp), (1|Source), (1|Site)
#binomial model

mtdna_hd_binomial_null <- glmer(cbind(success, failure) ~ bp_scale + range_position + (1|Family/Genus/spp) + (1|Source) + 
                         (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail")

#checking fit with DHARMa
mtdna_hd_binomial_null_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_null, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_hd_binomial_null_sim)
plotResiduals(mtdna_hd_binomial_null_sim)

#against specific predictors
plotResiduals(mtdna_hd_binomial_null_sim, mtdna_small_He$bp)
plotResiduals(mtdna_hd_binomial_null_sim, mtdna_small_He$range_position)

######## Build sst models (include sst variables) ########
#variables to include: bp, position in spp range, sst predictor, (1|MarkerName), (1|Family/Genus/spp), (1|Source), (1|Site)
#binomial model

#subset to only those with sst data
mtdna_small_He <- subset(mtdna_small_He, sst.BO_sstmean != "NA") #should look at where these are...

#scale sst data (for convergence issues)
mtdna_small_He$sstmean_scale <- scale(as.numeric(mtdna_small_He$sst.BO_sstmean))

#scale?
mtdna_hd_binomial_sstmean <- glmer(cbind(success, failure) ~ bp_scale + range_position + sstmean_scale + (1|Family/Genus/spp) + (1|Source) + 
                                  (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#Script for building models w/env predictors

remove(list = ls())

library(tidyverse)
library(gridExtra)
library(plotrix)
library(lme4)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(effects)
library(sjPlot)

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)
msat <- read.csv("output/msatloci_assembled.csv", stringsAsFactors = FALSE)
mtdna_env <- read.csv("output/mtdna_env.csv", stringsAsFactors = FALSE)
msat_env <- read.csv("output/msat_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)

#merge dataframes
mtdna <- merge(mtdna, mtdna_env[, c('X', 'sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin', 
                                    'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                    'chloroA.BO_chlomax', 'chloroA.BO_chlomin')], all.x = TRUE)
msat <- cbind(msat[, -1], msat_env[, c('sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin', 
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

######## mtDNA Hd models ########

#### Clean up dataframe ####
#subset mtdna to remove He = NA columns
mtdna_small_He <- subset(mtdna_small, mtdna_small$He != "NA")

#scaling bp & making sure numeric
mtdna_small_He$bp <- as.numeric(mtdna_small_He$bp)
  mtdna_small_He$bp_scale <- scale(as.numeric(mtdna_small_He$bp))

#calculating position in range
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

#calculating successes and failures
mtdna_small_He$success <- round(mtdna_small_He$He*mtdna_small_He$n)
mtdna_small_He$failure<- round((1 - mtdna_small_He$He)*mtdna_small_He$n)

#### Build sst models (include sst variables) ####
#following best model structure from model.R script

## log transform sst data ##

##subset to only those with sst data
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$sst.BO_sstmean != "NA") #shouldn't remove any

#log transform sst data
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$sst.BO_sstmean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$sst.BO_sstrange != 0)
  mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$sst.BO_sstmax != 0)
  mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$sst.BO_sstmin != 0)

mtdna_small_He$logsstmean <- log10(mtdna_small_He$sst.BO_sstmean)
  mtdna_small_He$logsstrange <- log10(mtdna_small_He$sst.BO_sstrange)
  mtdna_small_He$logsstmax <- log10(mtdna_small_He$sst.BO_sstmax)
  mtdna_small_He$logsstmin <- log10(mtdna_small_He$sst.BO_sstmin)

#remove logsst = NA columns
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$logsstmean != "NaN")
  mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$logsstrange != "NaN")
  mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$logsstmax != "NaN")
  mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$logsstmin != "NaN")

## sst mean model ##
mtdna_hd_binomial_sstmean <- glmer(cbind(success, failure) ~ bp_scale + range_position + logsstmean + 
                                     (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                                   family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_sstmean_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_sstmean, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_sstmean_sim)
plotResiduals(mtdna_hd_binomial_sstmean_sim)
  plotResiduals(mtdna_hd_binomial_sstmean_sim, mtdna_small_He$sstmean_scale)

#marginal effects
plot_model(mtdna_hd_binomial_sstmean, type = "pred", pred.type = "re",
           terms = "logsstmean [all]")

## sst range model ##
mtdna_hd_binomial_sstrange <- glmer(cbind(success, failure) ~ bp_scale + range_position + logsstrange + 
                                      (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                                    family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                   control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_sstrange_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_sstrange, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_sstrange_sim)
plotResiduals(mtdna_hd_binomial_sstrange_sim)
  plotResiduals(mtdna_hd_binomial_sstrange_sim, mtdna_small_He$sstrange_scale)

#marginal effects
plot_model(mtdna_hd_binomial_sstrange, type = "pred", pred.type = "re",
           terms = "logsstrange [all]")

## sst max model ##
mtdna_hd_binomial_sstmax <- glmer(cbind(success, failure) ~ bp_scale + range_position + logsstmax + 
                                    (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                                  family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                    control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_sstmax_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_sstmax, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_sstmax_sim)
plotResiduals(mtdna_hd_binomial_sstmax_sim)
  plotResiduals(mtdna_hd_binomial_sstmax_sim, mtdna_small_He$sstmax_scale)

#marginal effects
plot_model(mtdna_hd_binomial_sstmean, type = "pred", pred.type = "re",
           terms = "logsstrange [all]")

## sst min model ##
mtdna_hd_binomial_sstmin <- glmer(cbind(success, failure) ~ bp_scale + range_position + logsstmin + 
                                    (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                                  family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                    control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_sstmin_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_sstmin, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_sstmin_sim)
plotResiduals(mtdna_hd_binomial_sstmin_sim)
  plotResiduals(mtdna_hd_binomial_sstmin_sim, mtdna_small_He$sstmin_scale)

#marginal effects
plot_model(mtdna_hd_binomial_sstmean, type = "pred", pred.type = "re",
           terms = "logsstmin [all]")

######## Build oxygen models (include mean dissolved oxygen variable) ########
#variables to include: bp, position in spp range, dissolved oxygen predictor, (1|MarkerName), (1|Family/Genus/spp), (1|Source), (1|Site)
#binomial model

#logtransform dissox data

#### log transform dissolved ox ####
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$BO_dissox != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

mtdna_small_He$logdissox <- log10(mtdna_small_He$BO_dissox)

#remove logdissox = NA columns
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$logdissox != "NaN")

##### diss oxy mean model ####
## diss oxy w/rp ##
mtdna_hd_binomial_dissox <- glmer(cbind(success, failure) ~ bp_scale + range_position + logdissox + (1|Family/Genus/spp) + (1|Source) + 
                                     (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                   control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_dissox_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_dissox, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_dissox_sim)
plotResiduals(mtdna_hd_binomial_dissox_sim)

#against dissox
plotResiduals(mtdna_hd_binomial_dissox_sim, mtdna_small_He$dissox_scale)

#look at partial residuals
dissox_eff <- effect("logdissox", residuals = TRUE, mtdna_hd_binomial_dissox)
plot(dissox_eff, smooth.residuals = TRUE)

## diss oxy w/out rp ##
mtdna_hd_norp_binomial_dissox <- glmer(cbind(success, failure) ~ bp_scale + logdissox + (1|Family/Genus/spp) + (1|Source) + 
                                    (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_norp_binomial_dissox_sim <- simulateResiduals(fittedModel = mtdna_hd_norp_binomial_dissox, n = 1000, plot = F)
plotQQunif(mtdna_hd_norp_binomial_dissox_sim)
plotResiduals(mtdna_hd_norp_binomial_dissox_sim)

#against dissox
plotResiduals(mtdna_hd_norp_binomial_dissox_sim, mtdna_small_He$dissox_scale)

#look at partial residuals
dissox_eff <- effect("logdissox", residuals = TRUE, mtdna_hd_norp_binomial_dissox)
plot(dissox_eff, smooth.residuals = TRUE)

######## Build chlorophyll A models (include chlorophyll A variables) ########
#variables to include: bp, position in spp range, chlorophyll A predictor, (1|MarkerName), (1|Family/Genus/spp), (1|Source), (1|Site)
#binomial model

#logtransform chlorophyll A data

#### log transform chlorophyll A ####
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$sst.BO_sstmean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$chloroA.BO_chlorange != 0)
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$chloroA.BO_chlomax != 0)
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$chloroA.BO_chlomin != 0)

mtdna_small_He$logchlomean <- log10(mtdna_small_He$sst.BO_sstmean)
mtdna_small_He$logchlorange <- log10(mtdna_small_He$chloroA.BO_chlorange)
mtdna_small_He$logchlomax <- log10(mtdna_small_He$chloroA.BO_chlomax)
mtdna_small_He$logchlomin <- log10(mtdna_small_He$chloroA.BO_chlomin)

#remove logchlo = NA columns
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$logchlomean != "Inf" | mtdna_small_He$logchlomean != "NaN")
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$logchlorange != "Inf")
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$logchlomax != "Inf")
mtdna_small_He <- subset(mtdna_small_He, mtdna_small_He$logchlomin != "Inf")
                         
##### chloroA mean model ####
## chloroA mean w/rp ##
mtdna_hd_binomial_chloromean <- glmer(cbind(success, failure) ~ bp_scale + range_position + logsstmean + I(logsstmean^2) + (1|Family/Genus/spp) + (1|Source) + 
                                    (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_chloromean_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_chloromean, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_chloromean_sim)
plotResiduals(mtdna_hd_binomial_chloromean_sim)

#against chloroA mean
plotResiduals(mtdna_hd_binomial_chloromean_sim, mtdna_small_He$logchlomean)

#look at partial residuals
chloromean_eff <- effect("logchlomean", residuals = TRUE, mtdna_hd_binomial_chloromean)
plot(chloromean_eff, smooth.residuals = TRUE)

## chloroA mean w/out rp ##
mtdna_hd_norp_binomial_chloromean <- glmer(cbind(success, failure) ~ bp_scale + logchlomean + I(logchlomean^2) + (1|Family/Genus/spp) + (1|Source) + 
                                        (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                      control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_norp_binomial_chloromean_sim <- simulateResiduals(fittedModel = mtdna_hd_norp_binomial_chloromean, n = 1000, plot = F)
plotQQunif(mtdna_hd_norp_binomial_chloromean_sim)
plotResiduals(mtdna_hd_norp_binomial_chloromean_sim)

#against chloroA mean
plotResiduals(mtdna_hd_norp_binomial_chloromean_sim, mtdna_small_He$logchlomean)

#look at partial residuals
chloromean_eff <- effect("logchlomean", residuals = TRUE, mtdna_hd_norp_binomial_chloromean)
plot(chloromean_eff, smooth.residuals = TRUE)

##### chloroA range model ####
## chloroA range w/rp ##
mtdna_hd_binomial_chlororange <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlorange + I(logchlorange^2) + (1|Family/Genus/spp) + (1|Source) + 
                                        (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                      control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_chlororange_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_chlororange, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_chlororange_sim)
plotResiduals(mtdna_hd_binomial_chlororange_sim)

#against chloroA range
plotResiduals(mtdna_hd_binomial_chlororange_sim, mtdna_small_He$logchlorange)

#look at partial residuals
chlororange_eff <- effect("logchlorange", residuals = TRUE, mtdna_hd_binomial_chlororange)
plot(chlororange_eff, smooth.residuals = TRUE)

## chloroA range w/out rp ##
mtdna_hd_norp_binomial_chlororange <- glmer(cbind(success, failure) ~ bp_scale + logchlorange + I(logchlorange^2) + (1|Family/Genus/spp) + (1|Source) + 
                                         (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                       control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_norp_binomial_chlororange_sim <- simulateResiduals(fittedModel = mtdna_hd_norp_binomial_chlororange, n = 1000, plot = F)
plotQQunif(mtdna_hd_norp_binomial_chlororange_sim)
plotResiduals(mtdna_hd_norp_binomial_chlororange_sim)

#against chloroA range
plotResiduals(mtdna_hd_norp_binomial_chlororange_sim, mtdna_small_He$logchlorange)

#look at partial residuals
chloromax_eff <- effect("logchlorange", residuals = TRUE, mtdna_hd_norp_binomial_chlororange)
plot(chloromax_eff, smooth.residuals = TRUE)

##### chloroA max model ####
## chloroA max w/rp ##
mtdna_hd_binomial_chloromax <- glmer(cbind(success, failure) ~ bp_scale + logchlomax + I(logchlomax^2) + range_position + (1|Family/Genus/spp) + (1|Source) + 
                                         (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                       control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_chloromax_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_chloromax, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_chloromax_sim)
plotResiduals(mtdna_hd_binomial_chloromax_sim)

#against chloroA max
plotResiduals(mtdna_hd_binomial_chloromax_sim, mtdna_small_He$logchlomax)

#look at partial residuals
chloromax_eff <- effect("logchlomax", residuals = TRUE, mtdna_hd_binomial_chloromax)
plot(chloromax_eff, smooth.residuals = TRUE)

## chloroA max w/outrp ##
mtdna_hd_norp_binomial_chloromax <- glmer(cbind(success, failure) ~ bp_scale + logchlomax + I(logchlomax^2) + (1|Family/Genus/spp) + (1|Source) + 
                                       (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                     control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_norp_binomial_chloromax_sim <- simulateResiduals(fittedModel = mtdna_hd_norp_binomial_chloromax, n = 1000, plot = F)
plotQQunif(mtdna_hd_norp_binomial_chloromax_sim)
plotResiduals(mtdna_hd_norp_binomial_chloromax_sim)

#against chloroA max
plotResiduals(mtdna_hd_norp_binomial_chloromax_sim, mtdna_small_He$logchlomax)

#look at partial residuals
chloromax_eff <- effect("logchlomax", residuals = TRUE, mtdna_hd_norp_binomial_chloromax)
plot(chloromax_eff, smooth.residuals = TRUE)

##### chloroA min model ####
## chloroA min w/rp ##
mtdna_hd_binomial_chloromin <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomin + I(logchlomin^2) + (1|Family/Genus/spp) + (1|Source) + 
                                       (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                     control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_chloromin_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_chloromin, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_chloromin_sim)
plotResiduals(mtdna_hd_binomial_chloromin_sim)

#against chloroA min
plotResiduals(mtdna_hd_binomial_chloromin_sim, mtdna_small_He$logchlomin)

#look at partial residuals
chloromin_eff <- effect("logchlomin", residuals = TRUE, mtdna_hd_binomial_chloromin)
plot(chloromin_eff, smooth.residuals = TRUE)

## chloroA min w/out rp ##
mtdna_hd_norp_binomial_chloromin <- glmer(cbind(success, failure) ~ bp_scale + logchlomin + I(logchlomin^2) + (1|Family/Genus/spp) + (1|Source) + 
                                       (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                                     control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_norp_binomial_chloromin_sim <- simulateResiduals(fittedModel = mtdna_hd_norp_binomial_chloromin, n = 1000, plot = F)
plotQQunif(mtdna_hd_norp_binomial_chloromin_sim)
plotResiduals(mtdna_hd_norp_binomial_chloromin_sim)

#against chloroA min
plotResiduals(mtdna_hd_norp_binomial_chloromin_sim, mtdna_small_He$logchlomin)

#look at partial residuals
chloromin_eff <- effect("logchlomin", residuals = TRUE, mtdna_hd_norp_binomial_chloromin)
plot(chloromin_eff, smooth.residuals = TRUE)

##############################################################################################################

######## Building models for mtdna pi ########
#following best model structure from model.R script

######## Cleaning mtdna Hd dataset ########

#subset mtdna to remove Pi = NA columns
mtdna_small_pi <- subset(mtdna_small, mtdna_small$Pi != "NA")

#### log transform pi ####
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$Pi != 0) #if any zeros will screw up log transformation (log10(0) is undefined, also probably shouldn't be there anyway)
mtdna_small_pi$logpi <- log10(mtdna_small_pi$Pi)

#remove logpi = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logpi != "Inf")

#### chalculating position in range ####
#fix character type
mtdna_small_pi$Centroid <- as.numeric(mtdna_small_pi$Centroid)
mtdna_small_pi$Half_RangeSize <- as.numeric(mtdna_small_pi$Half_RangeSize)

mtdna_small_pi$range_position <- NA #create column to fill in

for (i in 1:nrow(mtdna_small_pi)) { #gcalculate distance from range center as percentage (0-1 both sides of centroid, 1 = all the way at range edge)
  mtdna_small_pi$range_position[i] <- abs((mtdna_small_pi$lat[i] - mtdna_small_pi$Centroid[i])/mtdna_small_pi$Half_RangeSize[i])
}

#check those >1 and round to 1 (slight discrepancy btwn aquamaps and reality)
mtdna_check <- subset(mtdna_small_pi, mtdna_small_pi$range_position > 1) #often right at aquamaps limit, round to 1 and keep
mtdna_small_pi$range_position[mtdna_small_pi$range_position > 1] <- 1

#subset to only those with range_position
mtdna_small_pi <- subset(mtdna_small_pi, range_position != "NA")

######## build null model (only includes nuisance variables) ########
#variables to include: position in spp range, (1|MarkerName), (1|Family/Genus/spp), (1|Source), (1|Site)
#linear model

## null w/rp ##
mtdna_pi_linear_null <- lmer(logpi ~ range_position + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                               (1|Site), REML = FALSE, data = mtdna_small_pi, 
                             na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_null)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_linear_null_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_null, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_linear_null_sim)
plotResiduals(mtdna_pi_linear_null_sim)

#against specific predictors
plotResiduals(mtdna_pi_linear_null_sim, mtdna_small_pi$range_position)

## null w/out rp ##
mtdna_pi_norp_linear_null <- lmer(logpi ~ (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                               (1|Site), REML = FALSE, data = mtdna_small_pi, 
                             na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_norp_linear_null)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_norp_linear_null_sim <- simulateResiduals(fittedModel = mtdna_pi_norp_linear_null, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_norp_linear_null_sim)
plotResiduals(mtdna_pi_norp_linear_null_sim)

#against specific predictors
plotResiduals(mtdna_pi_norp_linear_null_sim, mtdna_small_pi$range_position)

######## Build sst models (include sst variables) ########
#variables to include: position in spp range, sst predictor, (1|MarkerName), (1|Family/Genus/spp), (1|Source), (1|Site)
#linear model

#subset to only those with sst data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$sst.BO_sstmean != "NA") #should look at where these are...

#### log transform sst data ####
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$sst.BO_sstmean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$sst.BO_sstrange != 0)
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$sst.BO_sstmax != 0)
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$sst.BO_sstmin != 0)

mtdna_small_pi$logsstmean <- log10(mtdna_small_pi$sst.BO_sstmean)
mtdna_small_pi$logsstrange <- log10(mtdna_small_pi$sst.BO_sstrange)
mtdna_small_pi$logsstmax <- log10(mtdna_small_pi$sst.BO_sstmax)
mtdna_small_pi$logsstmin <- log10(mtdna_small_pi$sst.BO_sstmin)

#remove logsst = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstmean != "NaN")
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstrange != "NaN")
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstmax != "NaN")
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstmin != "NaN")

#### sst mean models ####
## sst mean w/rp ##
mtdna_pi_linear_sstmean <- lmer(logpi ~ range_position + logsstmean + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                               (1|Site), REML = FALSE, data = mtdna_small_pi, 
                             na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_sstmean)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_linear_sstmean_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_sstmean, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_linear_sstmean_sim)
plotResiduals(mtdna_pi_linear_sstmean_sim)

#against sstmean
plotResiduals(mtdna_pi_linear_sstmean_sim, mtdna_small_pi$sstmean_scale)

#look at partial residuals
sstmean_eff <- effect("sstmean_scale", residuals = TRUE, mtdna_pi_linear_sstmean)
plot(sstmean_eff, smooth.residuals = TRUE)

## sst mean w/out rp ##
mtdna_pi_norp_linear_sstmean <- lmer(logpi ~ logsstmean + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                  (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

mtdna_pi_norp_linear_sstmean_CP <- lmer(logpi ~ logsstmean + Pelagic_Coastal + Pelagic_Coastal:logsstmean + 
                                          (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                       (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                     na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_norp_linear_sstmean)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_norp_linear_sstmean_sim <- simulateResiduals(fittedModel = mtdna_pi_norp_linear_sstmean, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_norp_linear_sstmean_sim)
plotResiduals(mtdna_pi_norp_linear_sstmean_sim)

#against sstmean
plotResiduals(mtdna_pi_norp_linear_sstmean_sim, mtdna_small_pi$sstmean_scale)

#look at partial residuals
sstmean_eff <- effect("logsstmean", residuals = TRUE, mtdna_pi_norp_linear_sstmean_CP)
plot(sstmean_eff, smooth.residuals = TRUE)

#marginal effects
CP_int_eff <- plot_model(mtdna_pi_norp_linear_sstmean_CP, type = "pred", terms = c("logsstmean [all]", "Pelagic_Coastal"))
sst_eff <- plot_model(mtdna_pi_norp_linear_sstmean_CP, type = "pred", terms = "logsstmean [all]")

#### sst range models ####
## sst range w/rp ##
mtdna_pi_linear_sstrange <- lmer(logpi ~ range_position + logsstrange + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                  (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_sstrange)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_linear_sstrange_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_sstrange, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_linear_sstrange_sim)
plotResiduals(mtdna_pi_linear_sstrange_sim)

#against sstrange
plotResiduals(mtdna_pi_linear_sstrange_sim, mtdna_small_pi$sstrange_scale)

#look at partial residuals
sstrange_eff <- effect("sstrange_scale", residuals = TRUE, mtdna_pi_linear_sstrange)
plot(sstrange_eff, smooth.residuals = TRUE)

## sst range w/out rp ##
mtdna_pi_norp_linear_sstrange <- lmer(logpi ~ logsstrange + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                   (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                 na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_norp_linear_sstrange)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_norp_linear_sstrange_sim <- simulateResiduals(fittedModel = mtdna_pi_norp_linear_sstrange, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_norp_linear_sstrange_sim)
plotResiduals(mtdna_pi_norp_linear_sstrange_sim)

#against sstrange
plotResiduals(mtdna_pi_norp_linear_sstrange_sim, mtdna_small_pi$sstrange_scale)

#look at partial residuals
sstrange_eff <- effect("sstrange_scale", residuals = TRUE, mtdna_pi_norp_linear_sstrange)
plot(sstrange_eff, smooth.residuals = TRUE)

#### sst max models ####
## sst max w/rp ##
mtdna_pi_linear_sstmax <- lmer(logpi ~ range_position + logsstmax + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                   (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                 na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_sstmax)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_linear_sstmax_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_sstmax, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_linear_sstmax_sim)
plotResiduals(mtdna_pi_linear_sstmax_sim)

#against sstmax
plotResiduals(mtdna_pi_linear_sstmaxsim, mtdna_small_pi$sstmax_scale)

#look at partial residuals
sstmax_eff <- effect("sstmax_scale", residuals = TRUE, mtdna_pi_linear_sstmax)
plot(sstmax_eff, smooth.residuals = TRUE)

## sst max w/out rp ##
mtdna_pi_norp_linear_sstmax <- lmer(logpi ~ logsstmax + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_norp_linear_sstmax)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_norp_linear_sstmax_sim <- simulateResiduals(fittedModel = mtdna_pi_norp_linear_sstmax, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_norp_linear_sstmax_sim)
plotResiduals(mtdna_pi_norp_linear_sstmax_sim)

#against sstmax
plotResiduals(mtdna_pi_norp_linear_sstmaxsim, mtdna_small_pi$sstmax_scale)

#look at partial residuals
sstmax_eff <- effect("sstmax_scale", residuals = TRUE, mtdna_pi_norp_linear_sstmax)
plot(sstmax_eff, smooth.residuals = TRUE)

#### sst min models ####
## sst min w/rp ##
mtdna_pi_linear_sstmin <- lmer(logpi ~ range_position + logsstmin + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_sstmin)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_linear_sstmin_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_sstmin, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_linear_sstmin_sim)
plotResiduals(mtdna_pi_linear_sstmin_sim)

#against sstmin
plotResiduals(mtdna_pi_linear_sstminsim, mtdna_small_pi$sstmin_scale)

#look at partial residuals
sstmin_eff <- effect("logsstmin", residuals = TRUE, mtdna_pi_linear_sstmin)
plot(sstmin_eff, smooth.residuals = TRUE)

## sst min w/out rp ##
mtdna_pi_norp_linear_sstmin <- lmer(logpi ~ logsstmin + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_norp_linear_sstmin)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_norp_linear_sstmin_sim <- simulateResiduals(fittedModel = mtdna_pi_norp_linear_sstmin, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_norp_linear_sstmin_sim)
plotResiduals(mtdna_pi_norp_linear_sstmin_sim)

#against sstmin
plotResiduals(mtdna_pi_norp_linear_sstminsim, mtdna_small_pi$sstmin_scale)

#look at partial residuals
sstmin_eff <- effect("logsstmin", residuals = TRUE, mtdna_pi_norp_linear_sstmin)
plot(sstmin_eff, smooth.residuals = TRUE)

######## Build oxygen models (include mean dissolved oxygen variable) ########
#variables to include: position in spp range, dissolved oxygen predictor, (1|MarkerName), (1|Family/Genus/spp), (1|Source), (1|Site)
#linear model

#log transform dissolved ox

#### log transform diss ox data ####
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$BO_dissox != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

mtdna_small_pi$logdissox <- log10(mtdna_small_pi$BO_dissox)

#remove logdissox = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logdissox != "NaN")

##### diss oxy mean model ####
## diss oxy mean w/rp ##
mtdna_pi_linear_dissox <- lmer(logpi ~ range_position + logdissox + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_dissox)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_linear_dissox_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_dissox, n = 1000, plot = F)
plotQQunif(mtdna_pi_linear_dissox_sim)
plotResiduals(mtdna_pi_linear_dissox_sim)

#against dissox
plotResiduals(mtdna_pi_linear_dissox_sim, mtdna_small_pi$dissox_scale)

#look at partial residuals
dissox_eff <- effect("logdissox", residuals = TRUE, mtdna_pi_linear_dissox)
plot(dissox_eff, smooth.residuals = TRUE)

## diss oxy mean w/out rp ##
mtdna_pi_norp_linear_dissox <- lmer(logpi ~ logdissox + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_norp_linear_dissox)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_norp_linear_dissox_sim <- simulateResiduals(fittedModel = mtdna_pi_norp_linear_dissox, n = 1000, plot = F)
plotQQunif(mtdna_pi_norp_linear_dissox_sim)
plotResiduals(mtdna_pi_norp_linear_dissox_sim)

#against dissox
plotResiduals(mtdna_pi_norp_linear_dissox_sim, mtdna_small_pi$dissox_scale)

#look at partial residuals
dissox_eff <- effect("logdissox", residuals = TRUE, mtdna_pi_norp_linear_dissox)
plot(dissox_eff, smooth.residuals = TRUE)

######## Build chlorophyll A models (include chlorophyll A variables) ########
#variables to include: position in spp range, chlorophyll A predictor, (1|MarkerName), (1|Family/Genus/spp), (1|Source), (1|Site)
#linear model

#logtransform chlorophyll A data

#### log transform chlorophyll A ####
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlorange != 0)
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomax != 0)
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomin != 0)

mtdna_small_pi$logchlomean <- log10(mtdna_small_pi$chloroA.BO_chlomean)
mtdna_small_pi$logchlorange <- log10(mtdna_small_pi$chloroA.BO_chlorange)
mtdna_small_pi$logchlomax <- log10(mtdna_small_pi$chloroA.BO_chlomax)
mtdna_small_pi$logchlomin <- log10(mtdna_small_pi$chloroA.BO_chlomin)

#remove logchlo = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlomean != "Inf")
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlorange != "Inf")
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlomax != "Inf")
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlomin != "Inf")

##### chloroA mean model ####
## chloroA mean w/rp ##
mtdna_pi_linear_chloromean <- lmer(logpi ~ range_position + logchlomean + I(logchlomean^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_chloromean)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_linear_chloromean_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_chloromean, n = 1000, plot = F)
plotQQunif(mtdna_pi_linear_chloromean_sim)
plotResiduals(mtdna_pi_linear_chloromean_sim)

#against chloroA mean
plotResiduals(mtdna_pi_linear_chloromean_sim, mtdna_small_pi$logchlomean)

#look at partial residuals
chloromean_eff <- effect("logchlomean", residuals = TRUE, mtdna_pi_linear_chloromean)
plot(chloromean_eff, smooth.residuals = TRUE)

## chloroA mean w/out rp ##
mtdna_pi_norp_linear_chloromean <- lmer(logpi ~ logchlomean + I(logchlomean^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                     (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                   na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_norp_linear_chloromean)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_norp_linear_chloromean_sim <- simulateResiduals(fittedModel = mtdna_pi_norp_linear_chloromean, n = 1000, plot = F)
plotQQunif(mtdna_pi_norp_linear_chloromean_sim)
plotResiduals(mtdna_pi_norp_linear_chloromean_sim)

#against chloroA mean
plotResiduals(mtdna_pi_norp_linear_chloromean_sim, mtdna_small_pi$logchlomean)

#look at partial residuals
chloromean_eff <- effect("logchlomean", residuals = TRUE, mtdna_pi_norp_linear_chloromean)
plot(chloromean_eff, smooth.residuals = TRUE)

##### chloroA range model ####
## chloroA range w/rp ##
mtdna_pi_linear_chlororange <- lmer(logpi ~ range_position + logchlorange + I(logchlorange^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                     (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                   na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_chlororange)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_linear_chlororange_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_chlororange, n = 1000, plot = F)
plotQQunif(mtdna_pi_linear_chlororange_sim)
plotResiduals(mtdna_pi_linear_chlororange_sim)

#against chloroA range
plotResiduals(mtdna_pi_linear_chlororange_sim, mtdna_small_pi$logchlorange)

#look at partial residuals
chlororange_eff <- effect("logchlorange", residuals = TRUE, mtdna_pi_linear_chlororange)
plot(chlororange_eff, smooth.residuals = TRUE)

## chloroA range w/out rp ##
mtdna_pi_norp_linear_chlororange <- lmer(logpi ~ logchlorange + I(logchlorange^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                      (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                    na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_norp_linear_chlororange)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_norp_linear_chlororange_sim <- simulateResiduals(fittedModel = mtdna_pi_norp_linear_chlororange, n = 1000, plot = F)
plotQQunif(mtdna_pi_norp_linear_chlororange_sim)
plotResiduals(mtdna_pi_norp_linear_chlororange_sim)

#against chloroA range
plotResiduals(mtdna_pi_norp_linear_chlororange_sim, mtdna_small_pi$logchlorange)

#look at partial residuals
chlororange_eff <- effect("logchlorange", residuals = TRUE, mtdna_pi_norp_linear_chlororange)
plot(chlororange_eff, smooth.residuals = TRUE)

##### chloroA max model ####
## chloroA max w/rp ##
mtdna_pi_linear_chloromax <- lmer(logpi ~ range_position + logchlomax + I(logchlomax^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                      (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                    na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_chloromax)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_linear_chloromax_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_chloromax, n = 1000, plot = F)
plotQQunif(mtdna_pi_linear_chloromax_sim)
plotResiduals(mtdna_pi_linear_chloromax_sim)

#against chloroA max
plotResiduals(mtdna_pi_linear_chloromax_sim, mtdna_small_pi$logchlomax)

#look at partial residuals
chloromax_eff <- effect("logchlomax", residuals = TRUE, mtdna_pi_linear_chloromax)
plot(chloromax_eff, smooth.residuals = TRUE)

## chloroA max w/out rp ##
mtdna_pi_norp_linear_chloromax <- lmer(logpi ~ logchlomax + I(logchlomax^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                    (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                  na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_norp_linear_chloromax)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_norp_linear_chloromax_sim <- simulateResiduals(fittedModel = mtdna_pi_norp_linear_chloromax, n = 1000, plot = F)
plotQQunif(mtdna_pi_norp_linear_chloromax_sim)
plotResiduals(mtdna_pi_norp_linear_chloromax_sim)

#against chloroA max
plotResiduals(mtdna_pi_norp_linear_chloromax_sim, mtdna_small_pi$logchlomax)

#look at partial residuals
chloromax_eff <- effect("logchlomax", residuals = TRUE, mtdna_pi_norp_linear_chloromax)
plot(chloromax_eff, smooth.residuals = TRUE)

##### chloroA min model ####
## chloroA min w/rp ##
mtdna_pi_linear_chloromin <- lmer(logpi ~ range_position + logchlomin + I(logchlomin^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                    (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                  na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_chloromin)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_linear_chloromin_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_chloromin, n = 1000, plot = F)
plotQQunif(mtdna_pi_linear_chloromin_sim)
plotResiduals(mtdna_pi_linear_chloromin_sim)

#against chloroA min
plotResiduals(mtdna_pi_linear_chloromin_sim, mtdna_small_pi$logchlomin)

#look at partial residuals
chloromin_eff <- effect("logchlomin", residuals = TRUE, mtdna_pi_linear_chloromin)
plot(chloromin_eff, smooth.residuals = TRUE)

## chloroA min w/out rp ##
mtdna_pi_norp_linear_chloromin <- lmer(logpi ~ logchlomin + I(logchlomin^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                    (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                  na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_norp_linear_chloromin)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#checking fit with DHARMa
mtdna_pi_norp_linear_chloromin_sim <- simulateResiduals(fittedModel = mtdna_pi_norp_linear_chloromin, n = 1000, plot = F)
plotQQunif(mtdna_pi_norp_linear_chloromin_sim)
plotResiduals(mtdna_pi_norp_linear_chloromin_sim)

#against chloroA min
plotResiduals(mtdna_pi_norp_linear_chloromin_sim, mtdna_small_pi$logchlomin)

#look at partial residuals
chloromin_eff <- effect("logchlomin", residuals = TRUE, mtdna_pi_norp_linear_chloromin)
plot(chloromin_eff, smooth.residuals = TRUE)

##############################################################################################################

######## Building models for msat He ########
#following best model structure from model.R script

######## Cleaning msat He dataset ########

#### removing missing n values ####
msat$n <- as.numeric(msat$n)
msat <- subset(msat, msat$n != "NA")

#### removing missing repeat values ####
msat$Repeat <- as.numeric(msat$Repeat)
msat <- subset(msat, msat$Repeat != "NA")

#### calculating position in range ####
#fix character type
msat$Centroid <- as.numeric(msat$Centroid)
msat$Half_RangeSize <- as.numeric(msat$Half_RangeSize)

msat$range_position <- NA #create column to fill in

for (i in 1:nrow(msat)) { #gcalculate distance from range center as percentage (0-1 both sides of centroid, 1 = all the way at range edge)
  msat$range_position[i] <- abs((msat$lat[i] - msat$Centroid[i])/msat$Half_RangeSize[i])
}

#check those >1 and round to 1 (slight discrepancy btwn aquamaps and reality)
msat_check <- subset(msat, msat$range_position > 1) #often right at aquamaps limit, round to 1 and keep
msat_check$latcomp <- msat_check$lat
msat$range_position[msat$range_position > 1] <- 1

#subset to only those with range_position
msat <- subset(msat, range_position != "NA")

#### Calculating successes & failures ####
msat$success <- round(msat$He*msat$n)
msat$failure <- round((1 - msat$He)*msat$n)

#### Adding random effect for every unit ####
#correcting for overdispersion/heteroskedasticity
msat$ID <- (1:24182)

######## Build null model (only includes nuisance variables) ########
#variables to include: PrimerNote, CrossSpp, position in spp range, (1|Family/Genus/spp), (1|Source), (1|Site), (1|ID)
#binomial model

## null model w/rp ##
msat_he_binomial_null <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + (1|Family/Genus/spp) + (1|Source) + 
                              (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                              control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_null_sim <- simulateResiduals(fittedModel = msat_he_binomial_null, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(msat_he_binomial_null_sim)
plotResiduals(msat_he_binomial_null_sim)

#against specific predictors
plotResiduals(msat_he_binomial_null_sim, msat$PrimerNote)
plotResiduals(msat_he_binomial_null_sim, msat$CrossSpp)
plotResiduals(msat_he_binomial_null_sim, msat$Repeat)
plotResiduals(msat_he_binomial_null_sim, msat$range_position)

## null model w/out rp ##
msat_he_binomial_null <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + (1|Family/Genus/spp) + (1|Source) + 
                                 (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                               control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_norp_binomial_null_sim <- simulateResiduals(fittedModel = msat_he_norp_binomial_null, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(msat_he_norp_binomial_null_sim)
plotResiduals(msat_he_norp_binomial_null_sim)

#against specific predictors
plotResiduals(msat_he_norp_binomial_null_sim, msat$PrimerNote)
plotResiduals(msat_he_norp_binomial_null_sim, msat$CrossSpp)
plotResiduals(msat_he_norp_binomial_null_sim, msat$Repeat)
plotResiduals(msat_he_norp_binomial_null_sim, msat$range_position)

######## Build sst models (include sst variables) ########
#variables to include: PrimerNote, CrossSpp, position in spp range, sst predictor, (1|Family/Genus/spp), (1|Source), (1|Site), (1|ID)
#binomial model

#subset to only those with sst data
msat <- subset(msat, sst.BO_sstmean != "NA") #should look at where these are...

#### log transform sst data ####
msat <- subset(msat, msat$sst.BO_sstmean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
msat <- subset(msat, msat$sst.BO_sstrange != 0)
msat <- subset(msat, msat$sst.BO_sstmax != 0)
msat <- subset(msat, msat$sst.BO_sstmin != 0)

msat$logsstmean <- log10(msat$sst.BO_sstmean)
msat$logsstrange <- log10(msat$sst.BO_sstrange)
msat$logsstmax <- log10(msat$sst.BO_sstmax)
msat$logsstmin <- log10(msat$sst.BO_sstmin)

#remove logsst = NA columns
msat <- subset(msat, msat$logsstmean != "NaN")
msat <- subset(msat, msat$logsstrange != "NaN")
msat <- subset(msat, msat$logsstmax != "NaN")
msat <- subset(msat, msat$logsstmin != "NaN")

msat$lat_scale <- scale(msat$lat)

##### sst mean model ####
## sst mean w/rp ##
msat_he_binomial_sstmean <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + logsstmean + (1|Family/Genus/spp) + (1|Source) + 
                                 (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                 control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_sstmean_sim <- simulateResiduals(fittedModel = msat_he_binomial_sstmean, n = 1000, plot = F)
plotQQunif(msat_he_binomial_sstmean_sim)
plotResiduals(msat_he_binomial_sstmean_sim)

#against sstmean
plotResiduals(msat_he_binomial_sstmean_sim, msat$sstmean_scale)

#look at partial residuals
sstmean_eff <- effect("logsstmean", residuals = TRUE, msat_he_binomial_sstmean)
plot(sstmean_eff, smooth.residuals = TRUE)

## sst mean w/out rp ##
msat_he_norp_binomial_sstmean <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + logsstmean + (1|Family/Genus/spp) + (1|Source) + 
                                    (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa"))

msat_he_norp_binomial_sstmean_CP <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + logsstmean + 
                                            Pelagic_Coastal + Pelagic_Coastal:logsstmean + (1|Family/Genus/spp) + (1|Source) + 
                                         (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                       control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_norp_binomial_sstmean_sim <- simulateResiduals(fittedModel = msat_he_norp_binomial_sstmean, n = 1000, plot = F)
plotQQunif(msat_he_norp_binomial_sstmean_sim)
plotResiduals(msat_he_norp_binomial_sstmean_sim)

#against sstmean
plotResiduals(msat_he_norp_binomial_sstmean_sim, msat$sstmean_scale)

#look at partial residuals
sstmean_eff <- effect("logsstmean", residuals = TRUE, msat_he_norp_binomial_sstmean)
plot(sstmean_eff, smooth.residuals = TRUE)

#marginal effects
CP_int_eff <- plot_model(msat_he_norp_binomial_sstmean_CP, type = "pred", terms = c("logsstmean [all]", "Pelagic_Coastal"))
sst_eff <- plot_model(msat_he_norp_binomial_sstmean_CP, type = "pred", terms = "logsstmean [all]")

##### sst range model ####
## sst range w/rp ##
msat_he_binomial_sstrange <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + logsstrange + (1|Family/Genus/spp) + (1|Source) + 
                                    (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                   control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_sstrange_sim <- simulateResiduals(fittedModel = msat_he_binomial_sstrange, n = 1000, plot = F)
plotQQunif(msat_he_binomial_sstrange_sim)
plotResiduals(msat_he_binomial_sstrange_sim)

#against sstrange
plotResiduals(msat_he_binomial_sstrange_sim, msat$sstrange_scale)

#look at partial residuals
sstrange_eff <- effect("sstrange_scale", residuals = TRUE, msat_he_binomial_sstrange)
plot(sstrange_eff, smooth.residuals = TRUE)

## sst range w/out rp ##
msat_he_norp_binomial_sstrange <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + logsstrange + (1|Family/Genus/spp) + (1|Source) + 
                                     (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                   control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_norp_binomial_sstrange_sim <- simulateResiduals(fittedModel = msat_he_norp_binomial_sstrange, n = 1000, plot = F)
plotQQunif(msat_he_norp_binomial_sstrange_sim)
plotResiduals(msat_he_norp_binomial_sstrange_sim)

#against sstrange
plotResiduals(msat_he_norp_binomial_sstrange_sim, msat$sstrange_scale)

#look at partial residuals
sstrange_eff <- effect("sstrange_scale", residuals = TRUE, msat_he_norp_binomial_sstrange)
plot(sstrange_eff, smooth.residuals = TRUE)

##### sst max model ####
## sst max w/rp ##
msat_he_binomial_sstmax <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + logsstmax + (1|Family/Genus/spp) + (1|Source) + 
                                     (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                   control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_sstmax_sim <- simulateResiduals(fittedModel = msat_he_binomial_sstmax, n = 1000, plot = F)
plotQQunif(msat_he_binomial_sstmax_sim)
plotResiduals(msat_he_binomial_sstmax_sim)

#against sstmax
plotResiduals(msat_he_binomial_sstmax_sim, msat$sstmax_scale)

#look at partial residuals
sstmax_eff <- effect("sstmax_scale", residuals = TRUE, msat_he_binomial_sstmax)
plot(sstmax_eff, smooth.residuals = TRUE)

## sst max w/out rp ##
msat_he_norp_binomial_sstmax <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + logsstmax + (1|Family/Genus/spp) + (1|Source) + 
                                   (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                 control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_norp_binomial_sstmax_sim <- simulateResiduals(fittedModel = msat_he_norp_binomial_sstmax, n = 1000, plot = F)
plotQQunif(msat_he_norp_binomial_sstmax_sim)
plotResiduals(msat_he_norp_binomial_sstmax_sim)

#against sstmax
plotResiduals(msat_he_norp_binomial_sstmax_sim, msat$sstmax_scale)

#look at partial residuals
sstmax_eff <- effect("sstmax_scale", residuals = TRUE, msat_he_norp_binomial_sstmax)
plot(sstmax_eff, smooth.residuals = TRUE)

##### sst min model ####
## sst min w/rp ##
msat_he_binomial_sstmin <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + logsstmin + (1|Family/Genus/spp) + (1|Source) + 
                                   (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                 control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_sstmin_sim <- simulateResiduals(fittedModel = msat_he_binomial_sstmin, n = 1000, plot = F)
plotQQunif(msat_he_binomial_sstmin_sim)
plotResiduals(msat_he_binomial_sstmin_sim)

#against sstmin
plotResiduals(msat_he_binomial_sstmin_sim, msat$sstmin_scale)

#look at partial residuals
sstmin_eff <- effect("sstmin_scale", residuals = TRUE, msat_he_binomial_sstmin)
plot(sstmin_eff, smooth.residuals = TRUE)

## sst min w/out rp ##
msat_he_norp_binomial_sstmin <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + logsstmin + (1|Family/Genus/spp) + (1|Source) + 
                                   (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                 control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_norp_binomial_sstmin_sim <- simulateResiduals(fittedModel = msat_he_norp_binomial_sstmin, n = 1000, plot = F)
plotQQunif(msat_he_norp_binomial_sstmin_sim)
plotResiduals(msat_he_norp_binomial_sstmin_sim)

#against sstmin
plotResiduals(msat_he_norp_binomial_sstmin_sim, msat$sstmin_scale)

#look at partial residuals
sstmin_eff <- effect("sstmin_scale", residuals = TRUE, msat_he_norp_binomial_sstmin)
plot(sstmin_eff, smooth.residuals = TRUE)

######## Build oxygen models (include mean dissolved oxygen variable) ########
#variables to include: PrimerNote, CrossSpp, position in spp range, dissolved oxygen predictor, (1|Family/Genus/spp), (1|Source), (1|Site), (1|ID)
#binomial model

#log transform data

#### log transform dissox data ####
msat <- subset(msat, msat$BO_dissox != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

msat$logdissox <- log10(msat$BO_dissox)

#remove logdissox = NA columns
msat <- subset(msat, msat$logdissox != "NaN")

##### diss oxy mean model ####
## diss oxy mean w/rp ##
msat_he_binomial_dissox <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + logdissox + (1|Family/Genus/spp) + (1|Source) + 
                                   (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                 control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_dissox_sim <- simulateResiduals(fittedModel = msat_he_binomial_dissox, n = 1000, plot = F)
plotQQunif(msat_he_binomial_dissox_sim)
plotResiduals(msat_he_binomial_dissox_sim)

#against dissox
plotResiduals(msat_he_binomial_dissox_sim, msat$dissox_scale)

#look at partial residuals
dissox_eff <- effect("logdissox", residuals = TRUE, msat_he_binomial_dissox)
plot(dissox_eff, smooth.residuals = TRUE)

## diss oxy mean w/out rp ##
msat_he_norp_binomial_dissox <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + logdissox + (1|Family/Genus/spp) + (1|Source) + 
                                   (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                 control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_norp_binomial_dissox_sim <- simulateResiduals(fittedModel = msat_he_norp_binomial_dissox, n = 1000, plot = F)
plotQQunif(msat_he_norp_binomial_dissox_sim)
plotResiduals(msat_he_norp_binomial_dissox_sim)

#against dissox
plotResiduals(msat_he_norp_binomial_dissox_sim, msat$dissox_scale)

#look at partial residuals
dissox_eff <- effect("logdissox", residuals = TRUE, msat_he_norp_binomial_dissox)
plot(dissox_eff, smooth.residuals = TRUE)

######## Build chlorophyll A models (include chlorophyll A variables) ########
#variables to include: PrimerNote, CrossSpp, position in spp range, chlorophyll A predictor, (1|Family/Genus/spp), (1|Source), (1|Site), (1|ID)
#binomial model

#logtransform chlorophyll A data

#### log transform chlorophyll A ####
msat <- subset(msat, msat$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
msat <- subset(msat, msat$chloroA.BO_chlorange != 0)
msat <- subset(msat, msat$chloroA.BO_chlomax != 0)
msat <- subset(msat, msat$chloroA.BO_chlomin != 0)

msat$logchlomean <- log10(msat$chloroA.BO_chlomean)
msat$logchlorange <- log10(msat$chloroA.BO_chlorange)
msat$logchlomax <- log10(msat$chloroA.BO_chlomax)
msat$logchlomin <- log10(msat$chloroA.BO_chlomin)

#remove logchlo = NA columns
msat <- subset(msat, msat$logchlomean != "Inf")
msat <- subset(msat, msat$logchlorange != "Inf")
msat <- subset(msat, msat$logchlomax != "Inf")
msat <- subset(msat, msat$logchlomin != "Inf")

##### chloroA mean model ####
## chloroA mean w/rp ##
msat_he_binomial_chloromean <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + logchlomean + I(logchlomean^2) + (1|Family/Genus/spp) + (1|Source) + 
                                   (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                 control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_chloromean_sim <- simulateResiduals(fittedModel = msat_he_binomial_chloromean, n = 1000, plot = F)
plotQQunif(msat_he_binomial_chloromean_sim)
plotResiduals(msat_he_binomial_chloromean_sim)

#against chloroA mean
plotResiduals(msat_he_binomial_chloromean_sim, msat$logchlomean)

#look at partial residuals
chloromean_eff <- effect("logchlomean", residuals = TRUE, msat_he_binomial_chloromean)
plot(chloromean_eff, smooth.residuals = TRUE)

## chloroA mean w/out rp ##
msat_he_norp_binomial_chloromean <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + logchlomean + I(logchlomean^2) + (1|Family/Genus/spp) + (1|Source) + 
                                       (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                     control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_norp_binomial_chloromean_sim <- simulateResiduals(fittedModel = msat_he_norp_binomial_chloromean, n = 1000, plot = F)
plotQQunif(msat_he_orp_binomial_chloromean_sim)
plotResiduals(msat_he_norp_binomial_chloromean_sim)

#against chloroA mean
plotResiduals(msat_he_norp_binomial_chloromean_sim, msat$logchlomean)

#look at partial residuals
chloromean_eff <- effect("logchlomean", residuals = TRUE, msat_he_norp_binomial_chloromean)
plot(chloromean_eff, smooth.residuals = TRUE)

##### chloroA range model ####
## chloroA range w/rp ##
msat_he_binomial_chlororange <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + logchlorange + I(logchlorange^2) + (1|Family/Genus/spp) + (1|Source) + 
                                       (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                     control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_chlororange_sim <- simulateResiduals(fittedModel = msat_he_binomial_chlororange, n = 1000, plot = F)
plotQQunif(msat_he_binomial_chlororange_sim)
plotResiduals(msat_he_binomial_chlororange_sim)

#against chloroA range
plotResiduals(msat_he_binomial_chlororange_sim, msat$logchlorange)

#look at partial residuals
chlororange_eff <- effect("logchlorange", residuals = TRUE, msat_he_binomial_chlororange)
plot(chlororange_eff, smooth.residuals = TRUE)

## chloroA range w/out rp ##
msat_he_norp_binomial_chlororange <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + logchlorange + I(logchlorange^2) + (1|Family/Genus/spp) + (1|Source) + 
                                        (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                      control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_norp_binomial_chlororange_sim <- simulateResiduals(fittedModel = msat_he_norp_binomial_chlororange, n = 1000, plot = F)
plotQQunif(msat_he_norp_binomial_chlororange_sim)
plotResiduals(msat_he_norp_binomial_chlororange_sim)

#against chloroA range
plotResiduals(msat_he_norp_binomial_chlororange_sim, msat$logchlorange)

#look at partial residuals
chlororange_eff <- effect("logchlorange", residuals = TRUE, msat_he_norp_binomial_chlororange)
plot(chlororange_eff, smooth.residuals = TRUE)

##### chloroA max model ####
## chloroA max w/rp ##
msat_he_binomial_chloromax <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + logchlomax + I(logchlomax^2) + (1|Family/Genus/spp) + (1|Source) + 
                                        (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                      control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_chloromax_sim <- simulateResiduals(fittedModel = msat_he_binomial_chloromax, n = 1000, plot = F)
plotQQunif(msat_he_binomial_chloromax_sim)
plotResiduals(msat_he_binomial_chloromax_sim)

#against chloroA max
plotResiduals(msat_he_binomial_chloromax_sim, msat$logchlomax)

#look at partial residuals
chloromax_eff <- effect("logchlomax", residuals = TRUE, msat_he_binomial_chloromax)
plot(chloromax_eff, smooth.residuals = TRUE)

## chloroA max w/out rp ##
msat_he_norp_binomial_chloromax <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + logchlomax + I(logchlomax^2) + (1|Family/Genus/spp) + (1|Source) + 
                                      (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                    control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_norp_binomial_chloromax_sim <- simulateResiduals(fittedModel = msat_he_norp_binomial_chloromax, n = 1000, plot = F)
plotQQunif(msat_he_norp_binomial_chloromax_sim)
plotResiduals(msat_he_norp_binomial_chloromax_sim)

#against chloroA max
plotResiduals(msat_he_norp_binomial_chloromax_sim, msat$logchlomax)

#look at partial residuals
chloromax_eff <- effect("logchlomax", residuals = TRUE, msat_he_norp_binomial_chloromax)
plot(chloromax_eff, smooth.residuals = TRUE)

##### chloroA min model ####
## chloroA min w/rp ##
msat_he_binomial_chloromin <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + logchlomin + I(logchlomin^2) + (1|Family/Genus/spp) + (1|Source) + 
                                      (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                    control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_chloromin_sim <- simulateResiduals(fittedModel = msat_he_binomial_chloromin, n = 1000, plot = F)
plotQQunif(msat_he_binomial_chloromin_sim)
plotResiduals(msat_he_binomial_chloromin_sim)

#against chloroA min
plotResiduals(msat_he_binomial_chloromin_sim, msat$logchlomin)

chloromin_eff <- effect("logchlomin", residuals = TRUE, msat_he_binomial_chloromin)
plot(chloromin_eff, smooth.residuals = TRUE)

## chloroA min w/out rp ##
msat_he_norp_binomial_chloromin <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + logchlomin + I(logchlomin^2) + (1|Family/Genus/spp) + (1|Source) + 
                                      (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                                    control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_norp_binomial_chloromin_sim <- simulateResiduals(fittedModel = msat_he_norp_binomial_chloromin, n = 1000, plot = F)
plotQQunif(msat_he_norp_binomial_chloromin_sim)
plotResiduals(msat_he_norp_binomial_chloromin_sim)

#against chloroA min
plotResiduals(msat_he_norp_binomial_chloromin_sim, msat$logchlomin)

chloromin_eff <- effect("logchlomin", residuals = TRUE, msat_he_norp_binomial_chloromin)
plot(chloromin_eff, smooth.residuals = TRUE)

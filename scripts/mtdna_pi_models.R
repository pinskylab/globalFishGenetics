#Script for building mtDNA pi models

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
mtdna_env <- read.csv("output/mtdna_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)

#merge dataframes
mtdna <- merge(mtdna, mtdna_env[, c('X', 'sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin', 
                                    'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                    'chloroA.BO_chlomax', 'chloroA.BO_chlomin')], all.x = TRUE)
mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 
                                  'Half_RangeSize', 'Centroid')], all.x = TRUE)

#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)

######################################################################################################################

######## Clean up dataframe ########
#do this at beginning with all variables, so that all models run on same dataset (for comparisons)

#### Calculate nuisance variables ####
#subset mtdna to remove Pi = NA columns
mtdna_small_pi <- subset(mtdna_small, mtdna_small$Pi != "NA")

#add position in range
#fix character type
mtdna_small_pi$Centroid <- as.numeric(mtdna_small_pi$Centroid)
mtdna_small_pi$Half_RangeSize <- as.numeric(mtdna_small_pi$Half_RangeSize)

mtdna_small_pi$range_position <- NA #create column to fill in

for (i in 1:nrow(mtdna_small_pi)) { #calculate distance from range center as percentage (0-1 both sides of centroid, 1 = all the way at range edge)
  mtdna_small_pi$range_position[i] <- abs((mtdna_small_pi$lat[i] - mtdna_small_pi$Centroid[i])/mtdna_small_pi$Half_RangeSize[i])
}

#check those >1 and round to 1 (slight discrepancy btwn aquamaps and reality)
mtdna_check <- subset(mtdna_small_pi, mtdna_small_pi$range_position > 1) #often right at aquamaps limit, round to 1 and keep
mtdna_small_pi$range_position[mtdna_small_pi$range_position > 1] <- 1

#### Log transform pi ####
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$Pi != 0) #if any zeros will screw up log transformation (log10(0) is undefined, also probably shouldn't be there anyway)
  mtdna_small_pi$logpi <- log10(mtdna_small_pi$Pi)

#remove logpi = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logpi != "Inf")

#### Calculate latitude, longitude variables ####
#calculate abslat
mtdna_small_pi$abslat <- abs(mtdna_small_pi$lat)

#scale geographic variables
mtdna_small_pi$lat_scale <- as.numeric(scale(mtdna_small_pi$lat))
mtdna_small_pi$abslat_scale <- as.numeric(scale(mtdna_small_pi$abslat))

#convert lon to radians
mtdna_small_pi$lon_360 <- mtdna_small_pi$lon + 180 #convert (-180,180) to (0,360)
  mtdna_small_pi$lon_rad <- (2*pi*mtdna_small_pi$lon_360)/360

#### Calculate environmental variables ####
## log transform sst data ##
#subset to only those with sst data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$sst.BO_sstmean != "NA") #shouldn't remove any

mtdna_small_pi$logsstmean <- log10(mtdna_small_pi$sst.BO_sstmean)
  mtdna_small_pi$logsstrange <- log10(mtdna_small_pi$sst.BO_sstrange)
  mtdna_small_pi$logsstmax <- log10(mtdna_small_pi$sst.BO_sstmax)
  mtdna_small_pi$logsstmin <- log10(mtdna_small_pi$sst.BO_sstmin)

#remove logsst = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstmean != "NaN")
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstrange != "NaN")
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstmax != "NaN")
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstmin != "NaN")

## log transform dissox data ##
#subset to only those with dissox data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$BO_dissox != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

mtdna_small_pi$logdissox <- log10(mtdna_small_pi$BO_dissox)

#remove logdissox = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logdissox != "NaN")

## log transform chlorophyll A ##
#subset to only those with chloroA data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlorange != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomax != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomin != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

mtdna_small_pi$logchlomean <- log10(mtdna_small_pi$sst.BO_sstmean)
  mtdna_small_pi$logchlorange <- log10(mtdna_small_pi$chloroA.BO_chlorange)
  mtdna_small_pi$logchlomax <- log10(mtdna_small_pi$chloroA.BO_chlomax)
  mtdna_small_pi$logchlomin <- log10(mtdna_small_pi$chloroA.BO_chlomin)

#remove logchlo = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlomean != "Inf" | 
                           mtdna_small_pi$logchlomean != "NaN")
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlorange != "Inf" | 
                           mtdna_small_pi$logchlorange != "NaN")
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlomax != "Inf" | 
                           mtdna_small_pi$logchlomax != "NaN")
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlomin != "Inf" | 
                           mtdna_small_pi$logchlomin != "NaN")

#############################################################################################################

######## Null model ########
#no bp or range_position, models with those perform worse and they aren't significant predictors
#spp not nested here bc get singular boundary fit otherwise

null_model_pi <- lmer(logpi ~  (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                      REML = FALSE, data = mtdna_small_pi, 
                      na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood

#pull p-values
coefs <- data.frame(coef(summary(null_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
null_model_pi_sim_output <- simulateResiduals(null_model_pi, plot = F)
plotQQunif(null_model_pi_sim_output)
plotResiduals(null_model_pi_sim_output)
testDispersion(null_model_pi_sim_output)

###################################################################################################################

######## Latitude & longitude models ########

#### lat model ####
lat_model_pi <- lmer(logpi ~ lat_scale + I(lat_scale^2) + (1|Family/Genus) +   
                       (1|Source) + (1|MarkerName), REML = FALSE, data = mtdna_small_pi, 
                     na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(lat_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
lat_model_pi_sim_output <- simulateResiduals(lat_model_pi, plot = F)
plotQQunif(lat_model_pi_sim_output) 
plotResiduals(lat_model_pi_sim_output)

#marginal effects
plot_model(lat_model_pi, type = "pred", pred.type = "re",
           terms = "lat_scale [all]")

#### abslat model ####
abslat_model_pi <- lmer(logpi ~ abslat + (1|Family/Genus) + 
                          (1|Source) + (1|MarkerName), REML = FALSE, data = mtdna_small_pi, 
                        na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(abslat_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
abslat_model_pi_sim_output <- simulateResiduals(abslat_model_pi, plot = F)
plotQQunif(abslat_model_pi_sim_output) 
plotResiduals(abslat_model_pi_sim_output)

#marginal effects
plot_model(abslat_model_pi, type = "pred", pred.type = "re",
           terms = "abslat [all]")

#### lon model ####
lon_model_pi <- lmer(logpi ~ sin(lon_rad) + cos(lon_rad) + (1|Family/Genus) + (1|Source) + 
                       (1|MarkerName), REML = FALSE, data = mtdna_small_pi, 
                     na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(lon_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
lon_model_pi_sim_output <- simulateResiduals(lon_model_pi, plot = F)
plotQQunif(lon_model_pi_sim_output) 
plotResiduals(lon_model_pi_sim_output)

#marginal effects
plot_model(lon_model_pi, type = "pred", pred.type = "re",
           terms = "lon_rad [all]")

#### lat & lon model ####
lat_lon_model_pi <- lmer(logpi ~ lat_scale + I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + 
                           (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                         REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", 
                         control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(lat_lon_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
lat_lon_model_pi_sim_output <- simulateResiduals(lat_lon_model_pi, plot = F)
plotQQunif(lat_lon_model_pi_sim_output) 
plotResiduals(lat_lon_model_pi_sim_output)
  plotResiduals(lat_lon_model_pi_sim_output, mtdna_small_pi$lat_scale)
  plotResiduals(lat_lon_model_pi_sim_output, mtdna_small_pi$lon_rad)

#marginal effects
plot_model(lat_lon_model_pi, type = "pred", pred.type = "re",
           terms = "lat_scale [all]")
plot_model(lat_lon_model_pi, type = "pred", pred.type = "re",
           terms = "lon_rad [all]")

#### abslat & lon model ####
abslat_lon_model_pi <- lmer(logpi ~ abslat_scale + sin(lon_rad) + cos(lon_rad) + 
                              (1|Family/Genus) + (1|Source) + (1|MarkerName), REML = FALSE, 
                            data = mtdna_small_pi, na.action = "na.fail", 
                            control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(abslat_lon_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
abslat_lon_model_pi_sim_output <- simulateResiduals(abslat_lon_model_pi, plot = F)
plotQQunif(abslat_lon_model_pi_sim_output) 
plotResiduals(abslat_lon_model_pi_sim_output)
  plotResiduals(abslat_lon_model_pi_sim_output, mtdna_small_pi$abslat_scale)
  plotResiduals(abslat_lon_model_pi_sim_output, mtdna_small_pi$lon_rad)

#marginal effects
plot_model(abslat_lon_model_pi, type = "pred", pred.type = "re",
           terms = "abslat_scale [all]")
plot_model(abslat_lon_model_pi, type = "pred", pred.type = "re",
           terms = "lon_rad [all]")

####################################################################################################################

######## Environmental models ########

#### sst mean model ####
mtdna_pi_linear_sstmean <- lmer(logpi ~ logsstmean + (1|Family/Genus) + (1|Source) + (1|MarkerName), REML = FALSE, data = mtdna_small_pi, 
                                na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_sstmean)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
mtdna_pi_linear_sstmean_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_sstmean, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_linear_sstmean_sim)
plotResiduals(mtdna_pi_linear_sstmean_sim)
  plotResiduals(mtdna_pi_linear_sstmean_sim, mtdna_small_pi$logsstmean)

#marginal effects
plot_model(mtdna_pi_linear_sstmean, type = "pred", pred.type = "re",
           terms = "logsstmean [all]")

#### sst range model ####
mtdna_pi_linear_sstrange <- lmer(logpi ~ logsstrange + (1|Family/Genus) + 
                                   (1|Source) + (1|MarkerName), REML = FALSE, data = mtdna_small_pi, 
                                 na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_sstrange)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
mtdna_pi_linear_sstrange_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_sstrange, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_linear_sstrange_sim)
plotResiduals(mtdna_pi_linear_sstrange_sim)
  plotResiduals(mtdna_pi_linear_sstrange_sim, mtdna_small_pi$logsstrange)

#marginal effects
plot_model(mtdna_pi_linear_sstrange, type = "pred", pred.type = "re",
           terms = "logsstrange [all]")

#### sst max model ####
mtdna_pi_linear_sstmax <- lmer(logpi ~ logsstmax + (1|Family/Genus) + (1|Source) + (1|MarkerName), REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_sstmax)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
mtdna_pi_linear_sstmax_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_sstmax, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_linear_sstmax_sim)
plotResiduals(mtdna_pi_linear_sstmax_sim)
  plotResiduals(mtdna_pi_linear_sstmaxsim, mtdna_small_pi$logsstmax)

#marginal effects
plot_model(mtdna_pi_linear_sstmax, type = "pred", pred.type = "re",
           terms = "logsstmax [all]")

#### sst min model ####
mtdna_pi_linear_sstmin <- lmer(logpi ~ logsstmin + (1|Family/Genus) + (1|Source) + (1|MarkerName), REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_sstmin)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
mtdna_pi_linear_sstmin_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_sstmin, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_pi_linear_sstmin_sim)
plotResiduals(mtdna_pi_linear_sstmin_sim)
  plotResiduals(mtdna_pi_linear_sstminsim, mtdna_small_pi$logsstmin)

#marginal effects
plot_model(mtdna_pi_linear_sstmin, type = "pred", pred.type = "re",
           terms = "logsstmin [all]")

#### diss oxy mean model ####
mtdna_pi_linear_dissox <- lmer(logpi ~ logdissox + (1|Family/Genus) + (1|Source) + (1|MarkerName), REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_dissox)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
mtdna_pi_linear_dissox_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_dissox, n = 1000, plot = F)
plotQQunif(mtdna_pi_linear_dissox_sim)
plotResiduals(mtdna_pi_linear_dissox_sim)
  plotResiduals(mtdna_pi_linear_dissox_sim, mtdna_small_pi$logdissox)

#marginal effects
plot_model(mtdna_pi_linear_dissox, type = "pred", pred.type = "re",
           terms = "logdissox [all]")

#### chloroA mean model ####
mtdna_pi_linear_chloromean <- lmer(logpi ~ logchlomean + I(logchlomean^2) + (1|Family/Genus) + 
                                     (1|Source) + (1|MarkerName), 
                                   REML = FALSE, data = mtdna_small_pi, 
                                   na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_chloromean)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
mtdna_pi_linear_chloromean_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_chloromean, n = 1000, plot = F)
plotQQunif(mtdna_pi_linear_chloromean_sim)
plotResiduals(mtdna_pi_linear_chloromean_sim)
  plotResiduals(mtdna_pi_linear_chloromean_sim, mtdna_small_pi$logchlomean)

#marginal effects
plot_model(mtdna_pi_linear_chloromean, type = "pred", pred.type = "re",
           terms = "logchlomean [all]")

#### chloroA range model ####
mtdna_pi_linear_chlororange <- lmer(logpi ~ logchlorange + I(logchlorange^2) + (1|Family/Genus) + 
                                      (1|Source) + (1|MarkerName), REML = FALSE, 
                                    data = mtdna_small_pi, na.action = "na.fail", 
                                    control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_chlororange)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
mtdna_pi_linear_chlororange_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_chlororange, n = 1000, plot = F)
plotQQunif(mtdna_pi_linear_chlororange_sim)
plotResiduals(mtdna_pi_linear_chlororange_sim)
  plotResiduals(mtdna_pi_linear_chlororange_sim, mtdna_small_pi$logchlorange)

#marginal effects
plot_model(mtdna_pi_linear_chlororange, type = "pred", pred.type = "re",
           terms = "logchlorange [all]")


#### chloroA max model ####
mtdna_pi_linear_chloromax <- lmer(logpi ~  logchlomax + I(logchlomax^2) + (1|Family/Genus) + 
                                    (1|Source) + (1|MarkerName), REML = FALSE, 
                                  data = mtdna_small_pi, na.action = "na.fail", 
                                  control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_chloromax)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
mtdna_pi_linear_chloromax_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_chloromax, n = 1000, plot = F)
plotQQunif(mtdna_pi_linear_chloromax_sim)
plotResiduals(mtdna_pi_linear_chloromax_sim)
  plotResiduals(mtdna_pi_linear_chloromax_sim, mtdna_small_pi$logchlomax)

#marginal effects
plot_model(mtdna_pi_linear_chloromax, type = "pred", pred.type = "re",
           terms = "logchlomax [all]")

#### chloroA min model ####
mtdna_pi_linear_chloromin <- lmer(logpi ~  logchlomin + I(logchlomin^2) + (1|Family/Genus) + (1|Source) + 
                                    (1|MarkerName), REML = FALSE, data = mtdna_small_pi, 
                                  na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(mtdna_pi_linear_chloromin)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
mtdna_pi_linear_chloromin_sim <- simulateResiduals(fittedModel = mtdna_pi_linear_chloromin, n = 1000, plot = F)
plotQQunif(mtdna_pi_linear_chloromin_sim)
plotResiduals(mtdna_pi_linear_chloromin_sim)
  plotResiduals(mtdna_pi_linear_chloromin_sim, mtdna_small_pi$logchlomin)

#marginal effects
plot_model(mtdna_pi_linear_chloromin, type = "pred", pred.type = "re",
           terms = "logchlomin [all]")

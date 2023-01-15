#Script for building msat He models

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
msat <- read.csv("output/msatloci_assembled.csv", stringsAsFactors = FALSE)
msat_env <- read.csv("output/msat_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)

#merge dataframes
msat <- cbind(msat[, -1], msat_env[, c('sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin', 
                                       'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                       'chloroA.BO_chlomax', 'chloroA.BO_chlomin')]) #merge not working for some reason, cbind bc in same order
msat <- merge(msat, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                'Centroid')], all.x = TRUE)

####################################################################################################################

######## Clean up dataframe ########
#do this at beginning with all variables, so that all models run on same dataset (for comparisons)

#### Calculate nuisance variables ####
#remove missing n values
msat$n <- as.numeric(msat$n)
  msat <- subset(msat, msat$n != "NA")

#remove missing He
msat <- subset(msat, msat$He != "NA")

#remove missing CrossSpp
msat <- subset(msat, msat$CrossSpp!= "NA")

#add range position
#fix character type
msat$Centroid <- as.numeric(msat$Centroid)
msat$Half_RangeSize <- as.numeric(msat$Half_RangeSize)

msat$range_position <- NA #create column to fill in

for (i in 1:nrow(msat)) { #calculate distance from range center as percentage (0-1 both sides of centroid, 1 = all the way at range edge)
  msat$range_position[i] <- abs((msat$lat[i] - msat$Centroid[i])/msat$Half_RangeSize[i])
}

#check those >1 and round to 1 (slight discrepancy btwn aquamaps and reality)
msat_check <- subset(msat, msat$range_position > 1) #often right at aquamaps limit, round to 1 and keep
msat$range_position[msat$range_position > 1] <- 1

#add random effect for every unit
msat$ID <- (1:28142)

#### Calculate success and failure ####
msat$success <- round(msat$He*msat$n) #number of heterozygotes
msat$failure<- round((1 - msat$He)*msat$n) #number of homozygotes

#### Calculate latitude, longitude variables ####
#calculate abslat
msat$abslat <- abs(msat$lat)

#scale geographic variables
msat$lat_scale <- as.numeric(scale(msat$lat))
msat$abslat_scale <- as.numeric(scale(msat$abslat))

#convert lon to radians
msat$lon_360 <- msat$lon + 180 #convert (-180,180) to (0,360)
msat$lon_rad <- (2*pi*msat$lon_360)/360

#### Calculate environmental variables ####
## log transform sst data ##
#subset to only those with sst data
msat <- subset(msat, msat$sst.BO_sstmean != "NA") #shouldn't remove any

msat$logsstmean <- log10(msat$sst.BO_sstmean)
  msat$logsstrange <- log10(msat$sst.BO_sstrange)
  msat$logsstmax <- log10(msat$sst.BO_sstmax)
  msat$logsstmin <- log10(msat$sst.BO_sstmin)

#remove logsst = NA columns
msat <- subset(msat, msat$logsstmean != "NaN")
  msat <- subset(msat, msat$logsstrange != "NaN")
  msat <- subset(msat, msat$logsstmax != "NaN")
  msat <- subset(msat, msat$logsstmin != "NaN")

## log transform dissox data ##
#subset to only those with dissox data
msat <- subset(msat, msat$BO_dissox != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

msat$logdissox <- log10(msat$BO_dissox)

#remove logdissox = NA columns
msat <- subset(msat, msat$logdissox != "NaN")

## log transform chlorophyll A ##
#subset to only those with chloroA data
msat <- subset(msat, msat$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat <- subset(msat, msat$chloroA.BO_chlorange != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat <- subset(msat, msat$chloroA.BO_chlomax != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat <- subset(msat, msat$chloroA.BO_chlomin != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

msat$logchlomean <- log10(msat$sst.BO_sstmean)
  msat$logchlorange <- log10(msat$chloroA.BO_chlorange)
  msat$logchlomax <- log10(msat$chloroA.BO_chlomax)
  msat$logchlomin <- log10(msat$chloroA.BO_chlomin)

#remove logchlo = NA columns
msat <- subset(msat, msat$logchlomean != "Inf" | 
                 msat$logchlomean != "NaN")
msat <- subset(msat, msat$logchlorange != "Inf" | 
                 msat$logchlorange != "NaN")
msat <- subset(msat, msat$logchlomax != "Inf" | 
                 msat$logchlomax != "NaN")
msat <- subset(msat, msat$logchlomin != "Inf" | 
                 msat$logchlomin != "NaN")

#############################################################################################################

######## Null model ########
#no repeat and no range_position

binomial_msat_null <- glmer(cbind(success, failure) ~  CrossSpp + (1|Family/Genus/spp) + 
                              (1|Source) + (1|ID), family = binomial, data = msat, 
                            na.action = "na.fail", control = glmerControl(optimizer = "bobyqa")) #for some reason these ones keep throwing missing data error?

#checking fit with DHARMa
binomial_msat_null_sim <- simulateResiduals(fittedModel = binomial_msat_null, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_null_sim)
plotResiduals(binomial_msat_null_sim)
testDispersion(binomial_msat_null_sim)

#against CrossSpp
plotResiduals(binomial_msat_null_sim, msat$CrossSpp)

#marginal effects
plot_model(binomial_msat_null, type = "pred", pred.type = "re",
           terms = "CrossSpp [all]")

#### lat model ####
binomial_msat_lat <- glmer(cbind(success, failure) ~ CrossSpp + lat_scale + 
                             I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|ID), 
                           family = binomial, data = msat, na.action = "na.fail", 
                           control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_msat_lat_sim <- simulateResiduals(fittedModel = binomial_msat_lat, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_lat_sim)
plotResiduals(binomial_msat_lat_sim)
  plotResiduals(binomial_msat_lat_sim, msat$lat_scale)

#marginal effects
plot_model(binomial_msat_lat, type = "pred", pred.type = "re",
           terms = "lat_scale [all]")

### abslat model ####
binomial_msat_abslat <- glmer(cbind(success, failure) ~ CrossSpp + abslat_scale + 
                                (1|Family/Genus/spp) + (1|Source) + 
                                (1|ID), family = binomial, data = msat, na.action = "na.fail", 
                              control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_msat_abslat_sim <- simulateResiduals(fittedModel = binomial_msat_abslat, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_abslat_sim)
plotResiduals(binomial_msat_abslat_sim)
  plotResiduals(binomial_msat_abslat_sim, msat$abslat_scale)

#marginal effects
plot_model(binomial_msat_abslat, type = "pred", pred.type = "re",
           terms = "abslat_scale [all]")

#### lon model ####
binomial_msat_lon <- glmer(cbind(success, failure) ~ CrossSpp + sin(lon_rad) + cos(lon_rad) + 
                             (1|Family/Genus/spp) + (1|Source) + (1|ID), 
                           family = binomial, data = msat, na.action = "na.fail", 
                           control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_msat_lon_sim <- simulateResiduals(fittedModel = binomial_msat_lon, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_lon_sim)
plotResiduals(binomial_msat_lon_sim)
  plotResiduals(binomial_msat_lon_sim, msat$lon_rad)

#marginal effects
plot_model(binomial_msat_lon, type = "pred", pred.type = "re",
           terms = "lon_rad [all]")

#### lat & lon model ####
binomial_msat_lat_lon <- glmer(cbind(success, failure) ~  CrossSpp + lat_scale + I(lat_scale^2) + 
                                 sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|ID), 
                               family = binomial, data = msat, na.action = "na.fail",
                               control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_msat_lat_lon_sim <- simulateResiduals(fittedModel = binomial_msat_lat_lon, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_lat_lon_sim)
plotResiduals(binomial_msat_lat_lon_sim)
  plotResiduals(binomial_msat_lat_lon_sim, msat$lon_rad)
  plotResiduals(binomial_msat_lat_lon_sim, msat$lat_scale)

#marginal effects
plot_model(binomial_msat_lat_lon, type = "pred", pred.type = "re",
           terms = "lon_rad [all]")
plot_model(binomial_msat_lat_lon, type = "pred", pred.type = "re",
           terms = "lat_scale [all]")

#### abslat & lon model ####
binomial_msat_abslat_lon <- glmer(cbind(success, failure) ~ CrossSpp + abslat_scale +  
                                    cos(lon_rad) + sin(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|ID), 
                                  family = binomial, data = msat, na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_msat_abslat_lon_sim <- simulateResiduals(fittedModel = binomial_msat_abslat_lon, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_abslat_lon_sim)
plotResiduals(binomial_msat_abslat_lon_sim)
  plotResiduals(binomial_msat_abslat_lon_sim, msat$lon_rad)
  plotResiduals(binomial_msat_abslat_lon_sim, msat$abslat_scale)

#marginal effects
plot_model(binomial_msat_abslat_lon, type = "pred", pred.type = "re",
           terms = "lon_rad [all]")
plot_model(binomial_msat_abslat_lon, type = "pred", pred.type = "re",
           terms = "abslat_scale [all]")

###################################################################################################################

######## Environmental models ########

##### sst mean model ####
msat_he_binomial_sstmean <- glmer(cbind(success, failure) ~  CrossSpp + logsstmean + (1|Family/Genus/spp) + 
                                    (1|Source) + (1|ID), family = binomial, data = msat, 
                                  na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_sstmean_sim <- simulateResiduals(fittedModel = msat_he_binomial_sstmean, n = 1000, plot = F)
plotQQunif(msat_he_binomial_sstmean_sim)
plotResiduals(msat_he_binomial_sstmean_sim)
  plotResiduals(msat_he_binomial_sstmean_sim, msat$logsstmean)

plot_model(msat_he_binomial_sstmean, type = "pred", pred.type = "re",
           terms = "logsstmean [all]")

#### sst range model ####
msat_he_binomial_sstrange <- glmer(cbind(success, failure) ~ CrossSpp + logsstrange + (1|Family/Genus/spp) + 
                                     (1|Source) + (1|ID), family = binomial, data = msat, 
                                   na.action = "na.fail", 
                                   control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_sstrange_sim <- simulateResiduals(fittedModel = msat_he_binomial_sstrange, n = 1000, plot = F)
plotQQunif(msat_he_binomial_sstrange_sim)
plotResiduals(msat_he_binomial_sstrange_sim)
  plotResiduals(msat_he_binomial_sstrange_sim, msat$logsstrange)

plot_model(msat_he_binomial_sstrange, type = "pred", pred.type = "re",
           terms = "logsstrange [all]")

#### sst max model ####
## sst max w/rp ##
msat_he_binomial_sstmax <- glmer(cbind(success, failure) ~ CrossSpp + logsstmax + (1|Family/Genus/spp) + 
                                   (1|Source) + (1|ID), family = binomial, data = msat, 
                                 na.action = "na.fail", 
                                 control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_sstmax_sim <- simulateResiduals(fittedModel = msat_he_binomial_sstmax, n = 1000, plot = F)
plotQQunif(msat_he_binomial_sstmax_sim)
plotResiduals(msat_he_binomial_sstmax_sim)
  plotResiduals(msat_he_binomial_sstmax_sim, msat$logsstmax)

plot_model(msat_he_binomial_sstmax, type = "pred", pred.type = "re",
           terms = "logsstmax[all]")

#### sst min model ####
msat_he_binomial_sstmin <- glmer(cbind(success, failure) ~ CrossSpp + logsstmin + (1|Family/Genus/spp) + 
                                   (1|Source) + (1|ID), family = binomial, data = msat, 
                                 na.action = "na.fail", 
                                 control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_sstmin_sim <- simulateResiduals(fittedModel = msat_he_binomial_sstmin, n = 1000, plot = F)
plotQQunif(msat_he_binomial_sstmin_sim)
plotResiduals(msat_he_binomial_sstmin_sim)
  plotResiduals(msat_he_binomial_sstmin_sim, msat$logsstmin)

plot_model(msat_he_binomial_sstmin, type = "pred", pred.type = "re",
           terms = "logsstmin [all]")

#### diss oxy mean model ####
msat_he_binomial_dissox <- glmer(cbind(success, failure) ~ CrossSpp + logdissox + (1|Family/Genus/spp) + 
                                   (1|Source) + (1|ID), family = binomial, data = msat, 
                                 na.action = "na.fail", 
                                 control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_dissox_sim <- simulateResiduals(fittedModel = msat_he_binomial_dissox, n = 1000, plot = F)
plotQQunif(msat_he_binomial_dissox_sim)
plotResiduals(msat_he_binomial_dissox_sim)
  plotResiduals(msat_he_binomial_dissox_sim, msat$logdissox)

plot_model(msat_he_binomial_dissox, type = "pred", pred.type = "re",
           terms = "logdissox [all]")

#### chloroA mean model ####
msat_he_binomial_chloromean <- glmer(cbind(success, failure) ~ CrossSpp + logchlomean + I(logchlomean^2) + 
                                       (1|Family/Genus/spp) + (1|Source) + (1|ID), family = binomial, 
                                     data = msat, na.action = "na.fail", 
                                     control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_chloromean_sim <- simulateResiduals(fittedModel = msat_he_binomial_chloromean, n = 1000, plot = F)
plotQQunif(msat_he_binomial_chloromean_sim)
plotResiduals(msat_he_binomial_chloromean_sim)
  plotResiduals(msat_he_binomial_chloromean_sim, msat$logchlomean)

plot_model(msat_he_binomial_chloromean, type = "pred", pred.type = "re",
           terms = "logchlomean [all]")

#### chloroA range model ####
msat_he_binomial_chlororange <- glmer(cbind(success, failure) ~ CrossSpp + logchlorange + I(logchlorange^2) + 
                                        (1|Family/Genus/spp) + (1|Source) + (1|ID), family = binomial, 
                                      data = msat, na.action = "na.fail", 
                                      control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_chlororange_sim <- simulateResiduals(fittedModel = msat_he_binomial_chlororange, n = 1000, plot = F)
plotQQunif(msat_he_binomial_chlororange_sim)
plotResiduals(msat_he_binomial_chlororange_sim)
  plotResiduals(msat_he_binomial_chlororange_sim, msat$logchlorange)

plot_model(msat_he_binomial_chlororange, type = "pred", pred.type = "re",
           terms = "logchlorange [all]")

#### chloroA max model ####
msat_he_binomial_chloromax <- glmer(cbind(success, failure) ~ CrossSpp + logchlomax + I(logchlomax^2) + 
                                      (1|Family/Genus/spp) + (1|Source) + (1|ID), family = binomial, 
                                    data = msat, na.action = "na.fail", 
                                    control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_chloromax_sim <- simulateResiduals(fittedModel = msat_he_binomial_chloromax, n = 1000, plot = F)
plotQQunif(msat_he_binomial_chloromax_sim)
plotResiduals(msat_he_binomial_chloromax_sim)
  plotResiduals(msat_he_binomial_chloromax_sim, msat$logchlomax)

plot_model(msat_he_binomial_chloromax, type = "pred", pred.type = "re",
           terms = "logchlomax [all]")

#### chloroA min model ####
msat_he_binomial_chloromin <- glmer(cbind(success, failure) ~ CrossSpp + logchlomin + I(logchlomin^2) + 
                                      (1|Family/Genus/spp) + (1|Source) + (1|ID), family = binomial, 
                                    data = msat, na.action = "na.fail", 
                                    control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
msat_he_binomial_chloromin_sim <- simulateResiduals(fittedModel = msat_he_binomial_chloromin, n = 1000, plot = F)
plotQQunif(msat_he_binomial_chloromin_sim)
plotResiduals(msat_he_binomial_chloromin_sim)
  plotResiduals(msat_he_binomial_chloromin_sim, msat$logchlomin)

plot_model(msat_he_binomial_chloromin, type = "pred", pred.type = "re",
           terms = "logchlomin [all]")

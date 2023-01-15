#Script for building mtDNA Hd models

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
#subset mtdna to remove Hd = NA columns
mtdna_small_hd <- subset(mtdna_small, mtdna_small$He != "NA")

#add position in range
mtdna_small_hd$Centroid <- as.numeric(mtdna_small_hd$Centroid) #change to numeric for calculations
mtdna_small_hd$Half_RangeSize <- as.numeric(mtdna_small_hd$Half_RangeSize)

mtdna_small_hd$range_position <- NA #create column to fill in

for (i in 1:nrow(mtdna_small_hd)) { #calculate distance from range center as percentage (0-1 both sides of centroid, 1 = all the way at range edge)
  mtdna_small_hd$range_position[i] <- abs((mtdna_small_hd$lat[i] - 
                                             mtdna_small_hd$Centroid[i])/
                                            mtdna_small_hd$Half_RangeSize[i])
}

#check those >1 and round to 1 (slight discrepancy btwn aquamaps and reality)
mtdna_check <- subset(mtdna_small_hd, mtdna_small_hd$range_position > 1) #often right at aquamaps limit, round to 1 and keep
mtdna_small_hd$range_position[mtdna_small_hd$range_position > 1] <- 1

#subset to only those with range_position
mtdna_small_hd <- subset(mtdna_small_hd, range_position != "NA")

#scale bp
mtdna_small_hd$bp_scale <- as.numeric(scale(as.numeric(mtdna_small_hd$bp)))

#### Calculate success and failure ####
mtdna_small_hd$success <- round(mtdna_small_hd$He*mtdna_small_hd$n) #essentially, number of heterozygotes
mtdna_small_hd$failure <- round((1 - mtdna_small_hd$He)*mtdna_small_hd$n) #number of homozygotes

#### Calculate latitude, longitude variables ####
#calculate abslat
mtdna_small_hd$abslat <- abs(mtdna_small_hd$lat)

#scale geographic variables
mtdna_small_hd$lat_scale <- as.numeric(scale(mtdna_small_hd$lat))
mtdna_small_hd$abslat_scale <- as.numeric(scale(mtdna_small_hd$abslat))

#convert lon to radians
mtdna_small_hd$lon_360 <- mtdna_small_hd$lon + 180 #convert (-180,180) to (0,360)
mtdna_small_hd$lon_rad <- (2*pi*mtdna_small_hd$lon_360)/360

#### Calculate environmental variables ####
## log transform sst data ##
#subset to only those with sst data
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$sst.BO_sstmean != "NA") #shouldn't remove any

mtdna_small_hd$logsstmean <- log10(mtdna_small_hd$sst.BO_sstmean)
  mtdna_small_hd$logsstrange <- log10(mtdna_small_hd$sst.BO_sstrange)
  mtdna_small_hd$logsstmax <- log10(mtdna_small_hd$sst.BO_sstmax)
  mtdna_small_hd$logsstmin <- log10(mtdna_small_hd$sst.BO_sstmin)

#remove logsst = NA columns
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logsstmean != "NaN")
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logsstrange != "NaN")
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logsstmax != "NaN")
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logsstmin != "NaN")

## log transform dissox data ##
#subset to only those with dissox data
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$BO_dissox != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

mtdna_small_hd$logdissox <- log10(mtdna_small_hd$BO_dissox)

#remove logdissox = NA columns
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logdissox != "NaN")

## log transform chlorophyll A ##
#subset to only those with chloroA data
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlorange != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomax != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomin != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

mtdna_small_hd$logchlomean <- log10(mtdna_small_hd$sst.BO_sstmean)
  mtdna_small_hd$logchlorange <- log10(mtdna_small_hd$chloroA.BO_chlorange)
  mtdna_small_hd$logchlomax <- log10(mtdna_small_hd$chloroA.BO_chlomax)
  mtdna_small_hd$logchlomin <- log10(mtdna_small_hd$chloroA.BO_chlomin)

#remove logchlo = NA columns
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logchlomean != "Inf" | 
                           mtdna_small_hd$logchlomean != "NaN")
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logchlorange != "Inf" | 
                             mtdna_small_hd$logchlorange != "NaN")
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logchlomax != "Inf" | 
                             mtdna_small_hd$logchlomax != "NaN")
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logchlomin != "Inf" | 
                             mtdna_small_hd$logchlomin != "NaN")

#############################################################################################################

######## Null model ########

binomial_null <- glmer(cbind(success, failure) ~ bp_scale + range_position + (1|Family/Genus/spp) + (1|Source) + 
                         (1|Site) + (1|MarkerName), family = binomial, 
                       data = mtdna_small_hd, na.action = "na.fail", 
                       control = glmerControl(optimizer = "bobyqa"))

#check fit with DHARMa
binomial_null_sim <- simulateResiduals(fittedModel = binomial_null, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_null_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_null_sim) #residuals against predicted value -- looking for uniformity

#against bp & range position
plotResiduals(binomial_null_sim, mtdna_small_hd$bp_scale)
plotResiduals(binomial_null_sim, mtdna_small_hd$range_position)

#test dispersion for binomial model
testDispersion(binomial_null_sim)

#marginal effects
plot_model(binomial_null, type = "pred", pred.type = "re",
           terms = "range_position [all]")
plot_model(binomial_null, type = "pred", pred.type = "re",
           terms = "bp_scale [all]")

#plot partial residuals
#bp_eff <- effect("bp_scale", residuals = TRUE, binomial_null)
#  plot(bp_eff, smooth.residuals = TRUE)
#x_bp <- as.data.frame(effects_bp_eff) #if want to turn these into ggplot

###################################################################################################################

######## Latitude & longitude models ########

#### lat model ####
binomial_lat <- glmer(cbind(success, failure) ~ bp_scale + range_position + lat_scale + 
                        I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                      family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                      control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_lat_sim <- simulateResiduals(fittedModel = binomial_lat, n = 1000, plot = F)
plotQQunif(binomial_lat_sim)
plotResiduals(binomial_lat_sim)
  plotResiduals(binomial_lat_sim, mtdna_small_hd$lat_scale)

#marginal effects
plot_model(binomial_lat, type = "pred", pred.type = "re",
           terms = "lat_scale [all]")

#### abslat model ####
binomial_abslat <- glmer(cbind(success, failure) ~ bp_scale + abslat_scale + 
                           (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                         family = binomial, data = mtdna_small_hd, 
                         na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_abslat_sim <- simulateResiduals(fittedModel = binomial_abslat, n = 1000, plot = F)
plotResiduals(binomial_abslat_sim)
  plotResiduals(binomial_abslat_sim, mtdna_small_hd$abslat_scale)

#marginal effects
plot_model(binomial_abslat, type = "pred", pred.type = "re",
           terms = "abslat_scale [all]")

#### lon model ####
binomial_lon <- glmer(cbind(success, failure) ~ bp_scale + range_position + sin(lon_rad) + cos(lon_rad) + 
                        (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                      family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                      control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_lon_sim <- simulateResiduals(fittedModel = binomial_lon, n = 1000, plot = F)
plotResiduals(binomial_lon_sim)
  plotResiduals(binomial_lon_sim, mtdna_small_hd$lon_scale)

#marginal effects
plot_model(binomial_lon, type = "pred", pred.type = "re",
           terms = "lon_rad [all]")

#### lat & lon model ####
binomial_lat_lon <- glmer(cbind(success, failure) ~ bp_scale + range_position + lat_scale + 
                            I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + 
                            (1|Source) + (1|Site) + (1|MarkerName), family = binomial, 
                          data = mtdna_small_hd, na.action = "na.fail", 
                          control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_lat_lon_sim <- simulateResiduals(fittedModel = binomial_lat_lon, n = 1000, plot = F)
plotQQunif(binomial_lat_lon_sim)
plotResiduals(binomial_lat_lon_sim)
  plotResiduals(binomial_lat_lon_sim, mtdna_small_hd$lat_scale)
  plotResiduals(binomial_lat_lon_sim, mtdna_small_hd$lon_rad)

#marginal effects
plot_model(binomial_lat_lon, type = "pred", pred.type = "re",
           terms = "lat_scale [all]")
plot_model(binomial_lat_lon, type = "pred", pred.type = "re",
           terms = "lon_rad [all]")

#### abslat & lon model ###
binomial_abslat_lon <- glmer(cbind(success, failure) ~ bp_scale + range_position + abslat_scale + 
                               sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + 
                               (1|Source) + (1|Site) + (1|MarkerName), family = binomial, 
                             data = mtdna_small_hd, na.action = "na.fail", 
                             control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_abslat_lon_sim <- simulateResiduals(fittedModel = binomial_abslat_lon, n = 1000, plot = F)
plotQQunif(binomial_abslat_lon_sim)
plotResiduals(binomial_abslat_lon_sim)
  plotResiduals(binomial_abslat_lon_sim, mtdna_small_hd$abslat_scale)
  plotResiduals(binomial_abslat_lon_sim, mtdna_small_hd$lon_rad)

#marginal effects
plot_model(binomial_abslat_lon, type = "pred", pred.type = "re",
           terms = "abslat_scale [all]")
plot_model(binomial_abslat_lon, type = "pred", pred.type = "re",
           terms = "lon_rad [all]")

###################################################################################################################

######## Environmental models ########

#### sst mean model ####
mtdna_hd_binomial_sstmean <- glmer(cbind(success, failure) ~ bp_scale + range_position + logsstmean + 
                                     (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                                   family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                                   control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_sstmean_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_sstmean, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_sstmean_sim)
plotResiduals(mtdna_hd_binomial_sstmean_sim)
  plotResiduals(mtdna_hd_binomial_sstmean_sim, mtdna_small_hd$sstmean_scale)

#marginal effects
plot_model(mtdna_hd_binomial_sstmean, type = "pred", pred.type = "re",
           terms = "logsstmean [all]")

#### sst range model ####
mtdna_hd_binomial_sstrange <- glmer(cbind(success, failure) ~ bp_scale + range_position + logsstrange + 
                                      (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                                    family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                                    control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_sstrange_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_sstrange, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_sstrange_sim)
plotResiduals(mtdna_hd_binomial_sstrange_sim)
  plotResiduals(mtdna_hd_binomial_sstrange_sim, mtdna_small_hd$sstrange_scale)

#marginal effects
plot_model(mtdna_hd_binomial_sstrange, type = "pred", pred.type = "re",
           terms = "logsstrange [all]")

#### sst max model ####
mtdna_hd_binomial_sstmax <- glmer(cbind(success, failure) ~ bp_scale + range_position + logsstmax + 
                                    (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                                  family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_sstmax_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_sstmax, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_sstmax_sim)
plotResiduals(mtdna_hd_binomial_sstmax_sim)
  plotResiduals(mtdna_hd_binomial_sstmax_sim, mtdna_small_hd$sstmax_scale)

#marginal effects
plot_model(mtdna_hd_binomial_sstmean, type = "pred", pred.type = "re",
           terms = "logsstrange [all]")

#### sst min model ####
mtdna_hd_binomial_sstmin <- glmer(cbind(success, failure) ~ bp_scale + range_position + logsstmin + 
                                    (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                                  family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_sstmin_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_sstmin, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_sstmin_sim)
plotResiduals(mtdna_hd_binomial_sstmin_sim)
  plotResiduals(mtdna_hd_binomial_sstmin_sim, mtdna_small_He$sstmin_scale)

#marginal effects
plot_model(mtdna_hd_binomial_sstmean, type = "pred", pred.type = "re",
           terms = "logsstmin [all]")

#### diss oxy mean model ####
mtdna_hd_binomial_dissox <- glmer(cbind(success, failure) ~ bp_scale + range_position + logdissox + 
                                    (1|Family/Genus/spp) + (1|Source) + (1|Site) + (1|MarkerName), 
                                  family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_dissox_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_dissox, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_dissox_sim)
plotResiduals(mtdna_hd_binomial_dissox_sim)
  plotResiduals(mtdna_hd_binomial_dissox_sim, mtdna_small_hd$dissox_scale)

#marginal effects
plot_model(mtdna_hd_binomial_dissox, type = "pred", pred.type = "re",
           terms = "logdissox [all]")

#### chloroA mean model ####
mtdna_hd_binomial_chloromean <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomean + 
                                        I(logchlomean^2) + (1|Family/Genus/spp) + (1|Source) + 
                                        (1|Site) + (1|MarkerName), family = binomial, 
                                      data = mtdna_small_hd, na.action = "na.fail", 
                                      control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_chloromean_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_chloromean, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_chloromean_sim)
plotResiduals(mtdna_hd_binomial_chloromean_sim)
  plotResiduals(mtdna_hd_binomial_chloromean_sim, mtdna_small_hd$logchlomean)

#marginal effects
plot_model(mtdna_hd_binomial_chloromean, type = "pred", pred.type = "re",
           terms = "logchlomean [all]")

##### chloroA range model ####
mtdna_hd_binomial_chlororange <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlorange + 
                                         I(logchlorange^2) + (1|Family/Genus/spp) + (1|Source) + 
                                         (1|Site) + (1|MarkerName), family = binomial, 
                                       data = mtdna_small_hd, na.action = "na.fail", 
                                       control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_chlororange_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_chlororange, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_chlororange_sim)
plotResiduals(mtdna_hd_binomial_chlororange_sim)
  plotResiduals(mtdna_hd_binomial_chlororange_sim, mtdna_small_hdlogchlorange)

#marginal effects
plot_model(mtdna_hd_binomial_chlororange, type = "pred", pred.type = "re",
           terms = "logchlorange [all]")

#### chloroA max model ####
mtdna_hd_binomial_chloromax <- glmer(cbind(success, failure) ~ bp_scale + logchlomax + I(logchlomax^2) + 
                                       range_position + (1|Family/Genus/spp) + (1|Source) + 
                                       (1|Site) + (1|MarkerName), family = binomial, 
                                     data = mtdna_small_hd, na.action = "na.fail", 
                                     control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_chloromax_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_chloromax, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_chloromax_sim)
plotResiduals(mtdna_hd_binomial_chloromax_sim)
  plotResiduals(mtdna_hd_binomial_chloromax_sim, mtdna_small_hd$logchlomax)

#marginal effects
plot_model(mtdna_hd_binomial_chloromax, type = "pred", pred.type = "re",
           terms = "logchlomax [all]")

#### chloroA min model ####
mtdna_hd_binomial_chloromin <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomin + 
                                       I(logchlomin^2) + (1|Family/Genus/spp) + (1|Source) + 
                                       (1|Site) + (1|MarkerName), family = binomial, 
                                     data = mtdna_small_hd, na.action = "na.fail", 
                                     control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#checking fit with DHARMa
mtdna_hd_binomial_chloromin_sim <- simulateResiduals(fittedModel = mtdna_hd_binomial_chloromin, n = 1000, plot = F)
plotQQunif(mtdna_hd_binomial_chloromin_sim)
plotResiduals(mtdna_hd_binomial_chloromin_sim)
  plotResiduals(mtdna_hd_binomial_chloromin_sim, mtdna_small_hdlogchlomin)

#marginal effects
plot_model(mtdna_hd_binomial_chloromin, type = "pred", pred.type = "re",
           terms = "logchlomin [all]")

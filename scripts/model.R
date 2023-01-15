#Script for building models

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
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)
  mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                    'Centroid')], all.x = TRUE)
  msat <- merge(msat, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                  'Centroid')], all.x = TRUE)

#clean up dataframes
#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)
  
######################################################################################################################

######## mtDNA Hd models ########
  
#### Clean up dataframe ####

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

#calculate abslat
mtdna_small_hd$abslat <- abs(mtdna_small_hd$lat)

#scale geographic variables
mtdna_small_hd$lat_scale <- as.numeric(scale(mtdna_small_hd$lat))
mtdna_small_hd$abslat_scale <- as.numeric(scale(mtdna_small_hd$abslat))

#scale bp
mtdna_small_hd$bp_scale <- as.numeric(scale(as.numeric(mtdna_small_hd$bp)))

#convert lon to radians
mtdna_small_hd$lon_360 <- mtdna_small_hd$lon + 180 #convert (-180,180) to (0,360)
mtdna_small_hd$lon_rad <- (2*pi*mtdna_small_hd$lon_360)/360

#calculate success and failure
mtdna_small_hd$success <- round(mtdna_small_hd$He*mtdna_small_hd$n) #essentially, number of heterozygotes
mtdna_small_hd$failure <- round((1 - mtdna_small_hd$He)*mtdna_small_hd$n) #number of homozygotes

#### null model ####
#null model (no abslat/lat/lon)
binomial_null <- glmer(cbind(success, failure) ~ bp_scale + range_position + (1|Family/Genus/spp) + (1|Source) + 
                         (1|Site) + (1|MarkerName), family = binomial, 
                       data = mtdna_small_He, na.action = "na.fail", 
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

##############################################################################################################

######## mtdna pi models ########

#### Clean up dataframe ####

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

#log transform and check again
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$Pi != 0) #if any zeros will screw up log transformation (log10(0) is undefined, also probably shouldn't be there anyway)
  mtdna_small_pi$logpi <- log10(mtdna_small_pi$Pi)

#remove logpi = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logpi != "Inf")

#calculate abslat
mtdna_small_pi$abslat <- abs(mtdna_small_pi$lat)

#scale geographic variables
mtdna_small_pi$lat_scale <- scale(mtdna_small_pi$lat)
  mtdna_small_pi$lat_scale <- as.numeric(mtdna_small_pi$lat_scale)
mtdna_small_pi$abslat_scale <- scale(mtdna_small_pi$abslat)
  mtdna_small_pi$abslat_scale <- as.numeric(mtdna_small_pi$abslat_scale)

#convert lon to radians
mtdna_small_pi$lon_360 <- mtdna_small_pi$lon + 180 #convert (-180,180) to (0,360)
  mtdna_small_pi$lon_rad <- (2*pi*mtdna_small_pi$lon_360)/360

#### null model ####
#null model (no abslat/lat/lon)
#no bp or range_position, models with those perform worse and they aren't significant predictors
#spp not nested here bc get singular boundary fit otherwise
null_model_pi <- lmer(logpi ~  (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                        (1|Site), REML = FALSE, data = mtdna_small_pi, 
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

#### lat model ###
lat_model_pi <- lmer(logpi ~ lat_scale + I(lat_scale^2) + (1|Family/Genus) + (1|Site) +  
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

#### abslat model ###
abslat_model_pi <- lmer(logpi ~ abslat_scale + (1|Family/Genus) + 
                          (1|Source) + (1|MarkerName) + (1|Site), REML = FALSE, data = mtdna_small_pi, 
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
           terms = "abslat_scale [all]")

### lon model ###
lon_model_pi <- lmer(logpi ~ sin(lon_rad) + cos(lon_rad) + (1|Family/Genus) + (1|Source) + 
                       (1|MarkerName) + (1|Site), REML = FALSE, data = mtdna_small_pi, 
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

#### lat & lon model ###
lat_lon_model_pi <- lmer(logpi ~ lat_scale + I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + 
                           (1|Family/Genus) + (1|Source) + (1|MarkerName) + (1|Site), 
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
abslat_lon_model_pi <- lmer(logpi ~ abslat_scale + I(abslat_scale^2) + sin(lon_rad) + cos(lon_rad) + 
                              (1|Family/Genus) + (1|Source) + (1|MarkerName) + (1|Site), REML = FALSE, 
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

##############################################################################################################

######## msat he models ########

#### Clean up dataframe ####

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

#calculate abslat
msat$abslat <- abs(msat$lat)

#scale geographic variables
msat$lat_scale <- as.numeric(scale(msat$lat))
msat$abslat_scale <- as.numeric(scale(msat$abslat))

#convert lon to radians
msat$lon_360 <- msat$lon + 180 #convert (-180,180) to (0,360)
msat$lon_rad <- (2*pi*msat$lon_360)/360

#calculate success and failure
msat$success <- round(msat$He*msat$n) #number of heterozygotes
msat$failure<- round((1 - msat$He)*msat$n) #number of homozygotes

#add random effect for every unit
msat$ID <- (1:28152)

#### null model ####
#null model (no abslat/lat/lon)
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
                                  I(lat_scale^2) + I(lat_scale^3) + (1|Family/Genus/spp) + (1|Source) + (1|ID), 
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

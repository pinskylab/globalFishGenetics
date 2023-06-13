################################################### Script to build msat He models ########################################################

#Nuclear microsatellite expected heterozygosity (He) data
#Binomial generalized linear mixed effect models for He
#Success = heterozygote, failure = homozygote
#Check model fits and for spatial autocorrelation in residuals with DHARMa

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(tidyverse) #v.2.0.0
library(lme4) #v.1.1-31
library(DHARMa) #v.0.4.6
library(sjPlot) #v.2.8.12
library(splines) #v.4.2.2

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

################################################################################################################################

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

#subset to only those with range_position
msat <- subset(msat, range_position != "NA")

#add random effect for every unit
msat$ID <- (1:26733)

#### Calculate success and failure ####
msat$success <- round(msat$He*msat$n) #number of heterozygotes
msat$failure<- round((1 - msat$He)*msat$n) #number of homozygotes

#### Calculate latitude, longitude variables ####
#calculate abslat
msat$abslat <- abs(msat$lat)

#scale geographic variables
msat$lat_scale <- as.numeric(scale(msat$lat))
msat$abslat_scale <- as.numeric(scale(msat$abslat))
msat$lon_scale <- as.numeric(scale(msat$lon))

#### Calculate environmental variables ####
## log transform chlorophyll A ##
#subset to only those with chloroA data
msat <- subset(msat, msat$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat <- subset(msat, msat$chloroA.BO_chlorange != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat <- subset(msat, msat$chloroA.BO_chlomax != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat <- subset(msat, msat$chloroA.BO_chlomin != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

msat$logchlomean <- log10(msat$chloroA.BO_chlomean)
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

#### Create coordinate dataframe for SAC tests ####
#grouping factor for residuals -- need to identify which ones have the same lat/lon
msat$coords <- paste(msat$lat, "_", msat$lon, sep = "")
  coords <- as.data.frame(unique(msat$coords))
  colnames(coords) <- c("coords")

#pull out unique coordinate for SAC test
coords_unique <- coords %>% separate(coords, sep = "_", c("lat_unique", "lon_unique"))
  x_unique <- coords_unique$lat_unique
  y_unique <- coords_unique$lon_unique

#############################################################################################################################

######## Null model ########

null_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + 
                         (1|Family/Genus) + (1|Source) + (1|ID),
                       data = msat, family = binomial,
                       na.action = "na.fail", nAGQ = 0, #nAGQ = 0 & #nAGQ = 1 qualitatively the same, so using 0 bc converges much faster
                       control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
null_model_he_sim <- simulateResiduals(fittedModel = null_model_he, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(null_model_he_sim)
plotResiduals(null_model_he_sim)
  plotResiduals(null_model_he_sim, msat$CrossSpp)
  plotResiduals(null_model_he_sim, msat$range_position)

#test for SAC
sim_recalc <- recalculateResiduals(null_model_he_sim, group = msat$coords) #need to group same lat/lons together
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

############################################################################################################################
  
######## Latitude & longitude models ########
  
#### lat model ####
lat_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + 
                        lat_scale + I(lat_scale^2) + (1|Family/Genus) + 
                        (1|Source) + (1|ID), 
                      family = binomial, data = msat, 
                      na.action = "na.fail", nAGQ = 0,
                      control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
lat_model_he_sim <- simulateResiduals(fittedModel = lat_model_he, plot = F)
plotQQunif(lat_model_he_sim)
plotResiduals(lat_model_he_sim)
  plotResiduals(lat_model_he_sim, msat$lat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(lat_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

### abslat model ####
abslat_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + abslat_scale + 
                           (1|Family/Genus) + (1|Source) + (1|ID), 
                         family = binomial, data = msat, 
                         na.action = "na.fail", nAGQ = 0,
                         control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
abslat_model_he_sim <- simulateResiduals(fittedModel = abslat_model_he, plot = F)
plotQQunif(abslat_model_he_sim)
plotResiduals(abslat_model_he_sim)
  plotResiduals(abslat_model_he_sim, msat$abslat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(abslat_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### lon model ####
lon_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + bs(lon_scale) + 
                        (1|Family/Genus) + (1|Source) + (1|ID), 
                      family = binomial, data = msat, 
                      na.action = "na.fail", nAGQ = 0,
                      control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
lon_model_he_sim <- simulateResiduals(fittedModel = lon_model_he, plot = F)
plotQQunif(lon_model_he_sim)
plotResiduals(lon_model_he_sim)
  plotResiduals(lon_model_he_sim, msat$lon_scale)

#test for SAC
sim_recalc <- recalculateResiduals(lon_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### lat & lon model ####
lat_lon_model_he <- glmer(cbind(success, failure) ~  CrossSpp + range_position + lat_scale + 
                            I(lat_scale^2) + bs(lon_scale) + (1|Family/Genus) + 
                            (1|Source) + (1|ID), 
                          family = binomial, data = msat, 
                          na.action = "na.fail", nAGQ = 0,
                          control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
lat_lon_model_he_sim <- simulateResiduals(fittedModel = lat_lon_model_he, plot = F)
plotQQunif(lat_lon_model_he_sim)
plotResiduals(lat_lon_model_he_sim)
  plotResiduals(lat_lon_model_he_sim, msat$lon_scale)
  plotResiduals(lat_lon_model_he_sim, msat$lat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(lat_lon_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### abslat & lon model ####
abslat_lon_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + 
                               abslat_scale +  bs(lon_scale) + (1|Family/Genus) + 
                               (1|Source) + (1|ID), 
                             family = binomial, data = msat, 
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
abslat_lon_model_he_sim <- simulateResiduals(fittedModel = abslat_lon_model_he, plot = F)
plotQQunif(abslat_lon_model_he_sim)
plotResiduals(abslat_lon_model_he_sim)
  plotResiduals(abslat_lon_model_he_sim, msat$lon_scale)
  plotResiduals(abslat_lon_model_he_sim, msat$abslat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(abslat_lon_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

###################################################################################################################

######## Environmental models ########

##### sst mean model ####
sstmean_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + sst.BO_sstmean + 
                            (1|Family/Genus) + (1|Source) + (1|ID), 
                          family = binomial, data = msat, 
                          na.action = "na.fail", nAGQ = 0,
                          control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
sstmean_model_he_sim <- simulateResiduals(fittedModel = sstmean_model_he, plot = F)
plotQQunif(sstmean_model_he_sim)
plotResiduals(sstmean_model_he_sim)
  plotResiduals(sstmean_model_he_sim, msat$sst.BO_sstmean)

#test for SAC
sim_recalc <- recalculateResiduals(sstmean_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
#### sst range model ####
sstrange_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + sst.BO_sstrange + 
                             (1|Family/Genus) + (1|Source) + (1|ID), 
                           family = binomial, data = msat, 
                           na.action = "na.fail", nAGQ = 0,
                           control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
sstrange_model_he_sim <- simulateResiduals(fittedModel = sstrange_model_he, plot = F)
plotQQunif(sstrange_model_he_sim)
plotResiduals(sstrange_model_he_sim)
  plotResiduals(sstrange_model_he_sim, msat$sst.BO_sstrange)

#test for SAC
sim_recalc <- recalculateResiduals(sstrange_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### sst max model ####
sstmax_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + sst.BO_sstmax + 
                           (1|Family/Genus) + (1|Source) + (1|ID), 
                         family = binomial, data = msat, 
                         na.action = "na.fail", nAGQ = 0,
                         control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
sstmax_model_he_sim <- simulateResiduals(fittedModel = sstmax_model_he, plot = F)
plotQQunif(sstmax_model_he_sim)
plotResiduals(sstmax_model_he_sim)
  plotResiduals(sstmax_model_he_sim, msat$sst.BO_sstmax)

#test for SAC
sim_recalc <- recalculateResiduals(sstmax_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### sst min model ####
sstmin_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + sst.BO_sstmin + 
                           (1|Family/Genus) + (1|Source) + (1|ID), 
                         family = binomial, data = msat, 
                         na.action = "na.fail", nAGQ = 0,
                         control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
sstmin_model_he_sim <- simulateResiduals(fittedModel = sstmin_model_he, plot = F)
plotQQunif(sstmin_model_he_sim)
plotResiduals(sstmin_model_he_sim)
  plotResiduals(sstmin_model_he_sim, msat$sst.BO_sstmin)

#test for SAC
sim_recalc <- recalculateResiduals(sstmin_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

### chloro mean model ####
chlomean_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + 
                             logchlomean + I(logchlomean^2) + 
                             (1|Family/Genus) + (1|Source) + (1|ID), 
                           family = binomial, data = msat, 
                           na.action = "na.fail", nAGQ = 0,
                           control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
chlomean_model_he_sim <- simulateResiduals(fittedModel = chlomean_model_he, plot = F)
plotQQunif(chlomean_model_he_sim)
plotResiduals(chlomean_model_he_sim)
  plotResiduals(chlomean_model_he_sim, msat$logchlomean)

#test for SAC
sim_recalc <- recalculateResiduals(chlomean_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### chloro range model ####
chlorange_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + 
                              logchlorange + I(logchlorange^2) + 
                              (1|Family/Genus) + (1|Source) + (1|ID), 
                            family = binomial, data = msat, 
                            na.action = "na.fail", nAGQ = 0,
                            control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
chlorange_model_he_sim <- simulateResiduals(fittedModel = chlorange_model_he, plot = F)
plotQQunif(chlorange_model_he_sim)
plotResiduals(chlorange_model_he_sim)
  plotResiduals(chlorange_model_he_sim, msat$logchlorange)
  
#test for SAC
sim_recalc <- recalculateResiduals(chlorange_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### chloro max model ####
chlomax_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + 
                            logchlomax + I(logchlomax^2) + (1|Family/Genus) + 
                            (1|Source) + (1|ID), 
                          family = binomial, data = msat, 
                          na.action = "na.fail", nAGQ = 0,
                          control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
chlomax_model_he_sim <- simulateResiduals(fittedModel = chlomax_model_he, plot = F)
plotQQunif(chlomax_model_he_sim)
plotResiduals(chlomax_model_he_sim)
  plotResiduals(chlomax_model_he_sim, msat$logchlomax)

#test for SAC
sim_recalc <- recalculateResiduals(chlomax_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### chloro min model ####
chlomin_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + 
                            logchlomin + I(logchlomin^2) + (1|Family/Genus) + 
                            (1|Source) + (1|ID), 
                          family = binomial, data = msat, 
                          na.action = "na.fail", nAGQ = 0,
                          control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
chlomin_model_he_sim <- simulateResiduals(fittedModel = chlomin_model_he, plot = F)
plotQQunif(chlomin_model_he_sim)
plotResiduals(chlomin_model_he_sim)
  plotResiduals(chlomin_model_he_sim, msat$logchlomin)

#test for SAC
sim_recalc <- recalculateResiduals(chlomin_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
################################################### Script to build msat He models ########################################################

#Nuclear microsatellite expected heterozygosity (He) data
#Beta generalized linear mixed effect models for He
#All independent variables scaled & centered (except chlorophyll, which is log-transformed)
#Check model fits and for spatial autocorrelation in residuals with DHARMa

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(tidyverse) #v.2.0.0
library(glmmTMB) #1.1.7
library(DHARMa) #v.0.4.6
library(splines) #v.4.2.2
library(performance) #0.10.4
library(MuMIn) #1.47.5

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
  msat$CrossSpp_scale <- as.numeric(scale(msat$CrossSpp))

#### Add range position ####
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

#scale range position
msat$range_pos_scale <- as.numeric(scale(msat$range_position))

#### Calculate latitude, longitude variables ####
#calculate abslat
msat$abslat <- abs(msat$lat)

#scale geographic variables
msat$lat_scale <- as.numeric(scale(msat$lat))
msat$abslat_scale <- as.numeric(scale(msat$abslat))
msat$lon_scale <- as.numeric(scale(msat$lon))

#### Calculate environmental variables ####
#scale SST variables
msat$sstmean_scale <- as.numeric(scale(msat$sst.BO_sstmean))

## log transform chlorophyll A ##
#subset to only those with chloroA data
msat <- subset(msat, msat$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat$logchlomean <- log10(msat$chloroA.BO_chlomean)

#remove logchlo = NA columns
msat <- subset(msat, msat$logchlomean != "Inf" | 
                 msat$logchlomean != "NaN")

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
  
beta_null_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                (1|Family/Genus) + (1|Source), 
                              data = msat, family = ordbeta, 
                              na.action = "na.fail")
  
#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(beta_null_model_he)

#checking fit with DHARMa
null_model_he_sim <- simulateResiduals(fittedModel = beta_null_model_he, plot = F) #creates "DHARMa" residuals from simulations
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
beta_lat_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                               lat_scale + I(lat_scale^2) +
                               (1|Family/Genus) + (1|Source) + 
                               (0 + lat_scale|Family), 
                              data = msat, family = ordbeta, 
                              na.action = "na.fail")  

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013) --> for random slopes, has to be this (extended version of Nakagawa pseudo-rsquared)
r.squaredGLMM(beta_lat_model_he)

#checking fit with DHARMa
lat_model_he_sim <- simulateResiduals(fittedModel = beta_lat_model_he, plot = F)
plotQQunif(lat_model_he_sim)
plotResiduals(lat_model_he_sim)
  plotResiduals(lat_model_he_sim, msat$lat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(lat_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### abslat model ####
beta_abslat_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + abslat_scale + 
                                  (1|Family/Genus) + (1|Source) + 
                                  (0 + abslat_scale|Family), 
                                data = msat, family = ordbeta, 
                                na.action = "na.fail") 

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(beta_abslat_model_he)

#checking fit with DHARMa
abslat_model_he_sim <- simulateResiduals(fittedModel = beta_abslat_model_he, plot = F)
plotQQunif(abslat_model_he_sim)
plotResiduals(abslat_model_he_sim)
  plotResiduals(abslat_model_he_sim, msat$abslat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(abslat_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### lon model ####
beta_lon_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + bs(lon_scale) +
                               (1|Family/Genus) + (1|Source) + 
                               (0 + lon_scale|Family), 
                             data = msat, family = ordbeta, 
                             na.action = "na.fail") 

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(beta_lon_model_he)

#checking fit with DHARMa
lon_model_he_sim <- simulateResiduals(fittedModel = beta_lon_model_he, plot = F)
plotQQunif(lon_model_he_sim)
plotResiduals(lon_model_he_sim)
  plotResiduals(lon_model_he_sim, msat$lon_scale)

#test for SAC
sim_recalc <- recalculateResiduals(lon_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
#### lat & lon model ####
beta_lat_lon_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                   lat_scale + I(lat_scale^2) + bs(lon_scale) +
                                   (1|Family/Genus) + (1|Source) + 
                                   (0 + lat_scale|Family) + (0 + lon_scale|Family), 
                                 data = msat, family = ordbeta, 
                                 na.action = "na.fail") 

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
 r.squaredGLMM(beta_lat_lon_model_he)

#checking fit with DHARMa
lat_lon_model_he_sim <- simulateResiduals(fittedModel = beta_lat_lon_model_he, plot = F)
plotQQunif(lat_lon_model_he_sim)
plotResiduals(lat_lon_model_he_sim)
  plotResiduals(lat_lon_model_he_sim, msat$lon_scale)
  plotResiduals(lat_lon_model_he_sim, msat$lat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(lat_lon_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
#### abslat & lon model ####
beta_abslat_lon_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                      abslat_scale + bs(lon_scale) +
                                      (1|Family/Genus) + (1|Source) + 
                                      (0 + abslat_scale|Family) + (0 + lon_scale|Family), 
                                    data = msat, family = ordbeta, 
                                    na.action = "na.fail") 

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(beta_abslat_lon_model_he)

#checking fit with DHARMa
abslat_lon_model_he_sim <- simulateResiduals(fittedModel = beta_abslat_lon_model_he, plot = F)
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
beta_sstmean_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + sstmean_scale +
                                   (1|Family/Genus) + (1|Source) + 
                                   (0 + sstmean_scale|Family),
                                 data = msat, family = ordbeta, 
                                 na.action = "na.fail")   
  
#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(beta_sstmean_model_he)

#checking fit with DHARMa
sstmean_model_he_sim <- simulateResiduals(fittedModel = beta_sstmean_model_he, plot = F)
plotQQunif(sstmean_model_he_sim)
plotResiduals(sstmean_model_he_sim)
  plotResiduals(sstmean_model_he_sim, msat$sst.BO_sstmean)

#test for SAC
sim_recalc <- recalculateResiduals(sstmean_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
### chloro mean model ####
beta_chlomean_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale +
                                    logchlomean + I(logchlomean^2) + 
                                    (1|Family/Genus) + (1|Source) + 
                                    (0 + logchlomean|Family),
                                   data = msat, family = ordbeta, 
                                   na.action = "na.fail")  

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(beta_chlomean_model_he)
  
#checking fit with DHARMa
chlomean_model_he_sim <- simulateResiduals(fittedModel = beta_chlomean_model_he, plot = F)
plotQQunif(chlomean_model_he_sim)
plotResiduals(chlomean_model_he_sim)
  plotResiduals(chlomean_model_he_sim, msat$logchlomean)

#test for SAC
sim_recalc <- recalculateResiduals(chlomean_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
### SST mean & chloro mean model ####
beta_sstmean_chlomean_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                            sstmean_scale + logchlomean + I(logchlomean^2) + 
                                            (1|Family/Genus) + (1|Source) + 
                                            (0 + sstmean_scale|Family) + (0 + logchlomean|Family),
                                          data = msat, family = ordbeta, 
                                          na.action = "na.fail")  
  
#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(beta_sstmean_chlomean_model_he)
  
#checking fit with DHARMa
sstmean_chlomean_model_he_sim <- simulateResiduals(fittedModel = beta_sstmean_chlomean_model_he, plot = F)
plotQQunif(sstmean_chlomean_model_he_sim)
plotResiduals(sstmean_chlomean_model_he_sim)
  plotResiduals(sstmean_chlomean_model_he_sim, msat$logchlomean)

#test for SAC
sim_recalc <- recalculateResiduals(sstmean_chlomean_model_he_sim, group = msat$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
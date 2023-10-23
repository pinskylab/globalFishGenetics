################################################### Script to build mtDNA Hd models ########################################################

#Mitochondrial haplotype diversity (Hd) data
#Beta generalized linear mixed effect models for Hd
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

#subset mtdna to remove Hd = NA columns
mtdna_small_hd <- subset(mtdna_small, mtdna_small$He != "NA")

#### Calculate nuisance variables ####
#scale bp
mtdna_small_hd$bp_scale <- as.numeric(scale(as.numeric(mtdna_small_hd$bp)))

#### Add range position ####
#fix character type
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

#scale range_position
mtdna_small_hd$range_pos_scale <- as.numeric(scale(mtdna_small_hd$range_position))

#### Calculate latitude, longitude variables ####
#calculate abslat
mtdna_small_hd$abslat <- abs(mtdna_small_hd$lat)

#scale geographic variables
mtdna_small_hd$lat_scale <- as.numeric(scale(mtdna_small_hd$lat))
mtdna_small_hd$abslat_scale <- as.numeric(scale(mtdna_small_hd$abslat))
mtdna_small_hd$lon_scale <- as.numeric(scale(mtdna_small_hd$lon))

#### Calculate environmental variables ####
#scale SST variables
mtdna_small_hd$sstmean_scale <- as.numeric(scale(mtdna_small_hd$sst.BO_sstmean))

## log transform chlorophyll A ##
#subset to only those with chloroA data
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd$logchlomean <- log10(mtdna_small_hd$chloroA.BO_chlomean)

#remove logchlo = NA columns
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logchlomean != "Inf" | 
                           mtdna_small_hd$logchlomean != "NaN")
  
#### Create coordinate dataframe for SAC tests ####
#grouping factor for residuals -- need to identify which ones have the same lat/lon
mtdna_small_hd$coords <- paste(mtdna_small_hd$lat, "_", mtdna_small_hd$lon, sep = "")
  coords <- as.data.frame(unique(mtdna_small_hd$coords))
  colnames(coords) <- c("coords")
  
#pull out unique coordinate for SAC test
coords_unique <- coords %>% separate(coords, sep = "_", c("lat_unique", "lon_unique"))
  x_unique <- coords_unique$lat_unique
  y_unique <- coords_unique$lon_unique
    
#########################################################################################################################
  
######## Null model ########
  
beta_null_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + 
                                (1|Family/Genus) + (1|Source) + (1|MarkerName),
                              data = mtdna_small_hd, family = ordbeta, #use ordbeta (ordinal beta family) bc have a lot of 1s
                              na.action = "na.fail")  

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r2_nakagawa(beta_null_model_hd)

#check fit with DHARMa
null_model_hd_sim <- simulateResiduals(fittedModel = beta_null_model_hd, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(null_model_hd_sim) #QQplot
plotResiduals(null_model_hd_sim) #residuals against predicted value -- looking for uniformity
  plotResiduals(null_model_hd_sim, mtdna_small_hd$bp_scale)
  plotResiduals(null_model_hd_sim, mtdna_small_hd$range_position)
  
#test for SAC
sim_recalc <- recalculateResiduals(null_model_hd_sim, group = mtdna_small_hd$coords) #need to group same lat/lons together
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
###################################################################################################################

######## Latitude & longitude models ########

#### lat model ####
beta_lat_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + 
                               lat_scale + I(lat_scale^2) + (1|Family/Genus) +
                               (1|Source) + (1|MarkerName),
                             data = mtdna_small_hd, family = ordbeta, 
                             na.action = "na.fail")  

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r2_nakagawa(beta_lat_model_hd)
  
#check fit with DHARMa
lat_model_hd_sim <- simulateResiduals(fittedModel = beta_lat_model_hd, plot = F)
plotQQunif(lat_model_hd_sim)
plotResiduals(lat_model_hd_sim)
  plotResiduals(lat_model_hd_sim, mtdna_small_hd$lat_scale)
  
#test for SAC
sim_recalc <- recalculateResiduals(lat_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### abslat model ####
beta_abslat_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + 
                                  abslat_scale + (1|Family/Genus) +
                                  (1|Source) + (1|MarkerName),
                                data = mtdna_small_hd, family = ordbeta, 
                                na.action = "na.fail")  

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r2_nakagawa(beta_abslat_model_hd)

#check fit with DHARMa
abslat_model_hd_sim <- simulateResiduals(fittedModel = beta_abslat_model_hd, plot = F)
plotQQunif(abslat_model_hd_sim)
plotResiduals(abslat_model_hd_sim)
  plotResiduals(abslat_model_hd_sim, mtdna_small_hd$abslat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(abslat_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
#### lon model ####
beta_lon_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + 
                               bs(lon_scale) + (1|Family/Genus) + 
                               (1|Source) + (1|MarkerName),
                             data = mtdna_small_hd, family = ordbeta, 
                             na.action = "na.fail")  

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r2_nakagawa(beta_lon_model_hd)
  
#checking fit with DHARMa
lon_model_hd_sim <- simulateResiduals(fittedModel = beta_lon_model_hd, plot = F)
plotQQunif(lon_model_hd_sim)
plotResiduals(lon_model_hd_sim)
  plotResiduals(lon_model_hd_sim, mtdna_small_hd$lon_scale)
  
#test for SAC
sim_recalc <- recalculateResiduals(lon_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### lat & lon model ####
beta_lat_lon_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                                   I(lat_scale^2) + bs(lon_scale) + (1|Family/Genus) + 
                                   (1|Source) + (1|MarkerName),
                                 data = mtdna_small_hd, family = ordbeta, 
                                 na.action = "na.fail")  

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r2_nakagawa(beta_lat_lon_model_hd)
  
#checking fit with DHARMa
lat_lon_model_hd_sim <- simulateResiduals(fittedModel = beta_lat_lon_model_hd, plot = F)
plotQQunif(lat_lon_model_hd_sim)
plotResiduals(lat_lon_model_hd_sim)
  plotResiduals(lat_lon_model_hd_sim, mtdna_small_hd$lon_scale)
  plotResiduals(lat_lon_model_hd_sim, mtdna_small_hd$lat_scale)
  
#test for SAC
sim_recalc <- recalculateResiduals(lat_lon_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
#### abslat & lon model ####
beta_abslat_lon_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale + 
                                      bs(lon_scale) + (1|Family/Genus) + 
                                      (1|Source) + (1|MarkerName),
                                    data = mtdna_small_hd, family = ordbeta, 
                                    na.action = "na.fail")  

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r2_nakagawa(beta_abslat_lon_model_hd)
  
#checking fit with DHARMa
abslat_lon_model_hd_sim <- simulateResiduals(fittedModel = absbeta_lat_lon_model_hd, plot = F)
plotQQunif(abslat_lon_model_hd_sim)
plotResiduals(abslat_lon_model_hd_sim)
  plotResiduals(abslat_lon_model_hd_sim, mtdna_small_hd$lon_scale)
  plotResiduals(abslat_lon_model_hd_sim, mtdna_small_hd$abslat_scale)
  
#test for SAC
sim_recalc <- recalculateResiduals(abslat_lon_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
###################################################################################################################

######## Environmental models ########

#### sst mean model ####
beta_sstmean_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                   (1|Family/Genus) + (1|Source) + (1|MarkerName),
                                 data = mtdna_small_hd, family = ordbeta, 
                                 na.action = "na.fail")  
  
#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r2_nakagawa(beta_sstmean_model_hd)
  
#checking fit with DHARMa
sstmean_model_hd_sim <- simulateResiduals(fittedModel = beta_sstmean_model_hd, plot = F)
plotQQunif(sstmean_model_hd_sim)
plotResiduals(sstmean_model_hd_sim)
  plotResiduals(sstmean_model_hd_sim, mtdna_small_hd$sst.BO_sstmean)
  
#test for SAC
sim_recalc <- recalculateResiduals(sstmean_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### chloro mean model ####
beta_chlomean_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean +
                                    I(logchlomean^2) + (1|Family/Genus) + 
                                    (1|Source) + (1|MarkerName),
                                  data = mtdna_small_hd, family = ordbeta, 
                                  na.action = "na.fail")  
  
#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r2_nakagawa(beta_chlomean_model_hd)

#checking fit with DHARMa
chlomean_model_hd_sim <- simulateResiduals(fittedModel = beta_chlomean_model_hd, plot = F)
plotQQunif(chlomean_model_hd_sim)
plotResiduals(chlomean_model_hd_sim)
  plotResiduals(chlomean_model_hd_sim, mtdna_small_hd$logchlomean)
  
#test for SAC
sim_recalc <- recalculateResiduals(chlomean_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
beta_sstmean_chlomean_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                            logchlomean + I(logchlomean^2) + (1|Family/Genus) + 
                                            (1|Source) + (1|MarkerName),
                                          data = mtdna_small_hd, family = ordbeta, 
                                          na.action = "na.fail")  

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r2_nakagawa(beta_sstmean_chlomean_model_hd)

#checking fit with DHARMa
sstmean_chlomean_model_hd_sim <- simulateResiduals(fittedModel = beta_sstmean_chlomean_model_hd, plot = F)
plotQQunif(sstmean_chlomean_model_hd_sim)
plotResiduals(sstmean_chlomean_model_hd_sim)
  plotResiduals(sstmean_chlomean_model_hd_sim, mtdna_small_hd$logchlomean)
  plotResiduals(sstmean_chlomean_model_hd_sim, mtdna_small_hd$sstmean_scale)

#test for SAC
sim_recalc <- recalculateResiduals(sstmean_chlomean_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

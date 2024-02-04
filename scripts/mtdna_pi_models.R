################################################### Script to build mtDNA pi models ########################################################

#Mitochondrial average pairwise diversity (pi) data
#Linear generalized linear mixed effect models for pi
#All independent variables scaled & centered (except chlorophyll, which is log-transformed)
#Pi log-transformed
#Check model fits and for spatial autocorrelation in residuals with DHARMa

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(tidyverse) #v.2.0.0
library(lme4) #v.1.1-31
library(DHARMa) #v.0.4.6
library(splines) #v.4.2.2
library(performance) #0.10.4
library(MuMIn) #1.47.5

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

#### Add range position ####
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

#subset to only those with range_position
mtdna_small_pi <- subset(mtdna_small_pi, range_position != "NA")

#scale range position
mtdna_small_pi$range_pos_scale <- as.numeric(scale(mtdna_small_pi$range_position))

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
mtdna_small_pi$lon_scale <- as.numeric(scale(mtdna_small_pi$lon))

#### Calculate environmental variables ####
#scale SST variables
mtdna_small_pi$sstmean_scale <- as.numeric(scale(mtdna_small_pi$sst.BO_sstmean))

## log transform chlorophyll A ##
#subset to only those with chloroA data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi$logchlomean <- log10(mtdna_small_pi$chloroA.BO_chlomean)

#remove logchlo = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlomean != "Inf" | 
                           mtdna_small_pi$logchlomean != "NaN")

#### Create coordinate dataframe for SAC tests ####
#grouping factor for residuals -- need to identify which ones have the same lat/lon
mtdna_small_pi$coords <- paste(mtdna_small_pi$lat, "_", mtdna_small_pi$lon, sep = "")
  coords <- as.data.frame(unique(mtdna_small_pi$coords))
    colnames(coords) <- c("coords")

#pull out unique coordinate for SAC test
coords_unique <- coords %>% separate(coords, sep = "_", c("lat_unique", "lon_unique"))
  x_unique <- coords_unique$lat_unique
  y_unique <- coords_unique$lon_unique

#############################################################################################################
  
######## Null model ########
  
null_model_pi <- lmer(logpi ~ range_pos_scale + (1|Family/Genus) + 
                        (1|Source) + (1|MarkerName), 
                      REML = FALSE, data = mtdna_small_pi,
                      na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(null_model_pi)
 
#pull p-values
coefs <- data.frame(coef(summary(null_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs
 
#checking fit with DHARMa
null_model_pi_sim_output <- simulateResiduals(null_model_pi, plot = F)
plotQQunif(null_model_pi_sim_output)
plotResiduals(null_model_pi_sim_output)
  
#test for SAC
sim_recalc <- recalculateResiduals(null_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords) #need to group same lat/lons together
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

###################################################################################################################

######## Latitude & longitude models ########

#### lat model ####
lat_model_pi <- lmer(logpi ~ range_pos_scale + lat_scale + I(lat_scale^2) +
                       (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                       (0 + lat_scale|Family), 
                     REML = FALSE, data = mtdna_small_pi, 
                     na.action = "na.fail",
                     control = lmerControl(optimizer = "bobyqa"))

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(lat_model_pi)
  
#pull p-values
coefs <- data.frame(coef(summary(lat_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
lat_model_pi_sim_output <- simulateResiduals(lat_model_pi, plot = F)
plotQQunif(lat_model_pi_sim_output) 
plotResiduals(lat_model_pi_sim_output)
  plotResiduals(lat_model_pi_output, mtdna_small_pi$lat_scale)
  
#test for SAC
sim_recalc <- recalculateResiduals(lat_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### abslat model ####
abslat_model_pi <- lmer(logpi ~ range_pos_scale + abslat_scale + 
                          (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                          (0 + abslat_scale|Family),
                        REML = FALSE, data = mtdna_small_pi, 
                        na.action = "na.fail", 
                        control = lmerControl(optimizer = "bobyqa"))

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(abslat_model_pi)
  
#pull p-values
coefs <- data.frame(coef(summary(abslat_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
abslat_model_pi_sim_output <- simulateResiduals(abslat_model_pi, plot = F)
plotQQunif(abslat_model_pi_sim_output) 
plotResiduals(abslat_model_pi_sim_output)
  plotResiduals(abslat_model_pi_output, mtdna_small_pi$abslat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(abslat_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### lon model ####
lon_model_pi <- lmer(logpi ~ range_pos_scale + bs(lon_scale) + 
                       (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                       (0 + lon_scale|Family), 
                     REML = FALSE, data = mtdna_small_pi, 
                     na.action = "na.fail",
                     control = lmerControl(optimizer = "bobyqa"))

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(lon_model_pi)
  
#pull p-values
coefs <- data.frame(coef(summary(lon_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
lon_model_pi_sim_output <- simulateResiduals(lon_model_pi, plot = F)
plotQQunif(lon_model_pi_sim_output) 
plotResiduals(lon_model_pi_sim_output)
  plotResiduals(lon_model_pi_output, mtdna_small_pi$lon_scale)

#test for SAC
sim_recalc <- recalculateResiduals(lon_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
#### lat & lon model ####
lat_lon_model_pi <- lmer(logpi ~ range_pos_scale + lat_scale + 
                           I(lat_scale^2) + bs(lon_scale) + 
                           (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                           (0 + lat_scale|Family) + (0 + lon_scale|Family), 
                         REML = FALSE, data = mtdna_small_pi, 
                         na.action = "na.fail",
                         control = lmerControl(optimizer = "bobyqa"))

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(lat_lon_model_pi)
  
#pull p-values
coefs <- data.frame(coef(summary(lat_lon_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
lat_lon_model_pi_sim_output <- simulateResiduals(lat_lon_model_pi, plot = F)
plotQQunif(lat_lon_model_pi_sim_output) 
plotResiduals(lat_lon_model_pi_sim_output)
  plotResiduals(lat_lon_model_pi_sim_output, mtdna_small_pi$lon_scale)
  plotResiduals(lat_lon_model_pi_sim_output, mtdna_small_pi$lat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(lat_lon_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
#### abslat & lon model #### 
abslat_lon_model_pi <- lmer(logpi ~ range_pos_scale + abslat_scale + bs(lon_scale) +
                              (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                              (0 + abslat_scale|Family) + (0 + lon_scale|Family), 
                            REML = FALSE, data = mtdna_small_pi, 
                            na.action = "na.fail",
                            control = lmerControl(optimizer = "bobyqa"))

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(abslat_lon_model_pi)
  
#pull p-values
coefs <- data.frame(coef(summary(abslat_lon_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
abslat_lon_model_pi_sim_output <- simulateResiduals(abslat_lon_model_pi, plot = F)
plotQQunif(abslat_lon_model_pi_sim_output) 
plotResiduals(abslat_lon_model_pi_sim_output)
  plotResiduals(abslat_lon_model_pi_sim_output, mtdna_small_pi$lon_scale)
  plotResiduals(abslat_lon_model_pi_sim_output, mtdna_small_pi$abslat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(abslat_lon_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

####################################################################################################################

######## Environmental models ########

#### sst mean model ####
sstmean_model_pi <- lmer(logpi ~ range_pos_scale + sstmean_scale + 
                           (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                           (0 + sstmean_scale|Family), 
                         REML = FALSE, data = mtdna_small_pi, 
                         na.action = "na.fail",
                         control = lmerControl(optimizer = "bobyqa"))
 
#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(sstmean_model_pi)
   
#pull p-values
coefs <- data.frame(coef(summary(sstmean_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
sstmean_model_pi_sim_output <- simulateResiduals(fittedModel = sstmean_model_pi,plot = F)
plotQQunif(sstmean_model_pi_sim_output)
plotResiduals(sstmean_model_pi_sim_output)
  plotResiduals(sstmean_model_pi_sim_output, mtdna_small_pi$sst.BO_sstmean)

#test for SAC
sim_recalc <- recalculateResiduals(sstmean_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### chloro mean model ####
chlomean_model_pi <- lmer(logpi ~ range_pos_scale + logchlomean + I(logchlomean^2) + 
                            (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                            (0 + logchlomean|Family), 
                          REML = FALSE, data = mtdna_small_pi, 
                          na.action = "na.fail",
                          control = lmerControl(optimizer = "bobyqa"))

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(chlomean_model_pi)

#pull p-values
coefs <- data.frame(coef(summary(chlomean_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
chlomean_model_pi_sim_output <- simulateResiduals(fittedModel = chlomean_model_pi, plot = F)
plotQQunif(chlomean_model_pi_sim_output)
plotResiduals(chlomean_model_pi_sim_output)
  plotResiduals(chlomean_model_pi_sim_output, mtdna_small_pi$logchlomean)

#test for SAC
sim_recalc <- recalculateResiduals(chlomean_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
#### SST mean & chloro mean model #### --> DOESN'T CONVERGE (singular fit)
sstmean_chlomean_model_pi <- lmer(logpi ~ range_pos_scale + sstmean_scale + 
                                    logchlomean + I(logchlomean^2) + 
                                    (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                                    (0 + sstmean_scale|Family) + (0 + logchlomean|Family), 
                                  REML = FALSE, data = mtdna_small_pi, 
                                  na.action = "na.fail",
                                  control = lmerControl(optimizer = "bobyqa"))

#calculate pseudo-rsquared (Nakagawa & Schielzeth 2013)
r.squaredGLMM(sstmean_chlomean_model_pi)

#pull p-values
coefs <- data.frame(coef(summary(sstmean_chlomean_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
sstmean_chlomean_model_pi_sim_output <- simulateResiduals(fittedModel = sstmean_chlomean_model_pi, plot = F)
plotQQunif(sstmean_chlomean_model_pi_sim_output)
plotResiduals(sstmean_chlomean_model_pi_sim_output)
  plotResiduals(sstmean_chlomean_model_pi_sim_output, mtdna_small_pi$logchlomean)

#test for SAC
sim_recalc <- recalculateResiduals(sstmean_chlomean_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
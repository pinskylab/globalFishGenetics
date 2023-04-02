#Script for building mtDNA pi models

remove(list = ls())

library(tidyverse)
library(lme4)
library(DHARMa)
library(sjPlot)
library(splines)
library(spdep)
library(spatialreg)

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)
mtdna_env <- read.csv("output/mtdna_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)
load("mtdna_pi_ME.Rdata") #for ME weights, so don't have to re-run every time

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

#subset to only those with range_position
mtdna_small_pi <- subset(mtdna_small_pi, range_position != "NA")

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
## log transform chlorophyll A ##
#subset to only those with chloroA data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlorange != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomax != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomin != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

mtdna_small_pi$logchlomean <- log10(mtdna_small_pi$chloroA.BO_chlomean)
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

#### SAC variables ####
#to account for spatial autocorrelation

## Spatial random effect ##
#bin in 1 degree lat/lon bands
mtdna_small_pi$latlonbins <- paste(round(mtdna_small_pi$lat, 0), "_", 
                                   round(mtdna_small_pi$lon, 0), sep = "")

## Moran's eigenvectors ##
#create matrix of geographic coordinates
lonlat_mat <- cbind(mtdna_small_pi$lon, mtdna_small_pi$lat)

#groups into nearest neighbor groups
pi_nb4 <- knearneigh(x = lonlat_mat, k = 4, longlat = TRUE) #creates matrix with the indices of points belonging to the sets of k nearest neighbors
  pi_nb4 <- knn2nb(pi_nb4) #converts knn object to neighbors list (e.g. these indices in matrix group together spatially as nearest neighbors)
  pi_nb4 <- make.sym.nb(pi_nb4) #checks for symmetry/transitivity

#calculate weights
pi_wt4 <- nb2listw(pi_nb4, style = "W") #creates spatial weights

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
  
null_model_pi <- lmer(logpi ~ range_position + (1|Family/Genus) + 
                        (1|Source) + (1|MarkerName), 
                      REML = FALSE, data = mtdna_small_pi,
                      na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
  
## with SAC ##
#create MI eigenvectors
set.seed(8484) #because bootstraps for ME

#null_ME <- ME(logpi ~ range_position,
#              data = mtdna_small_pi, listw = pi_wt4) #ME is a permutation bootstrap on Moran's I for regression residuals and picks which ones account for most autocorrelation

#models
null_model_pi_SA_RE <- lmer(logpi ~ range_position + (1|Family/Genus) + 
                              (1|Source) + (1|MarkerName) + (1|latlonbins), 
                            REML = FALSE, data = mtdna_small_pi,
                            na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

null_model_pi_SA_ME <- lmer(logpi ~ range_position + (1|Family/Genus) + 
                              (1|Source) + (1|MarkerName) + fitted(null_ME), 
                            REML = FALSE, data = mtdna_small_pi,
                            na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

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
lat_model_pi <- lmer(logpi ~ range_position + lat_scale + I(lat_scale^2) + 
                       (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                     REML = FALSE, data = mtdna_small_pi, 
                     na.action = "na.fail",
                     control = lmerControl(optimizer = "bobyqa"))
  
## with SAC ##
#create ME eigenvectors
set.seed(8484)
  
#lat_ME <- ME(logpi ~ range_position + lat_scale + I(lat_scale^2),
#             data = mtdna_small_pi, listw = pi_wt4)


#models
lat_model_pi_SA_RE <- lmer(logpi ~ range_position + lat_scale + I(lat_scale^2) + 
                             (1|Family/Genus) + (1|Source) + (1|MarkerName) + (1|latlonbins), 
                           REML = FALSE, data = mtdna_small_pi, 
                           na.action = "na.fail",
                           control = lmerControl(optimizer = "bobyqa"))

lat_model_pi_SA_ME <- lmer(logpi ~ range_position + lat_scale + I(lat_scale^2) + 
                             (1|Family/Genus) + (1|Source) + (1|MarkerName) + fitted(lat_ME), 
                           REML = FALSE, data = mtdna_small_pi, 
                           na.action = "na.fail",
                           control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(lat_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
lat_model_pi_sim_output <- simulateResiduals(lat_model_pi, plot = F)
plotQQunif(lat_model_pi_sim_output) 
plotResiduals(lat_model_pi_sim_output)
  plotResiduals(sstmin_model_pi_output, mtdna_small_pi$lat_scale)
  
#test for SAC
sim_recalc <- recalculateResiduals(lat_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### abslat model ####
abslat_model_pi <- lmer(logpi ~ range_position + abslat_scale + (1|Family/Genus) + 
                          (1|Source) + (1|MarkerName), 
                        REML = FALSE, data = mtdna_small_pi, 
                        na.action = "na.fail",
                        control = lmerControl(optimizer = "bobyqa"))

## with SAC ##
#create ME eigenvectors
set.seed(8484)
  
#abslat_ME <- ME(logpi ~ range_position + abslat_scale,
#                data = mtdna_small_pi, listw = pi_wt4)

#models
abslat_model_pi_SA_RE <- lmer(logpi ~ range_position + abslat_scale + (1|Family/Genus) +
                                (1|Source) + (1|MarkerName) + (1|latlonbins), 
                              REML = FALSE, data = mtdna_small_pi, 
                              na.action = "na.fail",
                              control = lmerControl(optimizer = "bobyqa"))

abslat_model_pi_SA_ME <- lmer(logpi ~ range_position + abslat_scale + (1|Family/Genus) +
                                (1|Source) + (1|MarkerName) + fitted(abslat_ME), 
                              REML = FALSE, data = mtdna_small_pi, 
                              na.action = "na.fail",
                              control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(abslat_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
abslat_model_pi_sim_output <- simulateResiduals(abslat_model_pi, plot = F)
plotQQunif(abslat_model_pi_sim_output) 
plotResiduals(abslat_model_pi_sim_output)
  plotResiduals(sstmin_model_pi_output, mtdna_small_pi$abslat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(abslat_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### lon model ####
lon_model_pi <- lmer(logpi ~ range_position + bs(lon_scale) + (1|Family/Genus) + 
                       (1|Source) + (1|MarkerName), 
                     REML = FALSE, data = mtdna_small_pi, 
                     na.action = "na.fail",
                     control = lmerControl(optimizer = "bobyqa"))
  
## with SAC ##
#create ME eigenvectors
set.seed(8484)

#lon_ME <- ME(logpi ~ range_position + bs(lon_scale),
#             data = mtdna_small_pi, listw = pi_wt4)

#models
lon_model_pi_SA_RE <- lmer(logpi ~ range_position + bs(lon_scale) + (1|Family/Genus) + 
                             (1|Source) + (1|MarkerName) + (1|latlonbins), 
                           REML = FALSE, data = mtdna_small_pi, 
                           na.action = "na.fail",
                           control = lmerControl(optimizer = "bobyqa"))

lon_model_pi_SA_ME <- lmer(logpi ~ range_position + bs(lon_scale) + (1|Family/Genus) + 
                             (1|Source) + (1|MarkerName) + fitted(lon_ME), 
                           REML = FALSE, data = mtdna_small_pi, 
                           na.action = "na.fail",
                           control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(lon_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
lon_model_pi_sim_output <- simulateResiduals(lon_model_pi, plot = F)
plotQQunif(lon_model_pi_sim_output) 
plotResiduals(lon_model_pi_sim_output)
  plotResiduals(sstmin_model_pi_output, mtdna_small_pi$lon_scale)

#test for SAC
sim_recalc <- recalculateResiduals(lon_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### lat & lon model ####
lat_lon_model_pi <- lmer(logpi ~ range_position + lat_scale + I(lat_scale^2) + 
                           bs(lon_scale) + (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                         REML = FALSE, data = mtdna_small_pi, 
                         na.action = "na.fail",
                         control = lmerControl(optimizer = "bobyqa"))
  
## with SAC ##
#create ME eigenvectors
set.seed(8484)

#lat_lon_ME <- ME(logpi ~ range_position + lat_scale + 
#                   I(lat_scale^2) + bs(lon_scale),
#                 data = mtdna_small_pi, listw = pi_wt4)

#models
lat_lon_model_pi_SA_RE <- lmer(logpi ~ range_position + lat_scale + I(lat_scale^2) + 
                                 bs(lon_scale) + (1|Family/Genus) + (1|Source) + 
                                 (1|MarkerName) + (1|latlonbins), 
                               REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail",
                               control = lmerControl(optimizer = "bobyqa"))

lat_lon_model_pi_SA_ME <- lmer(logpi ~ range_position + lat_scale + I(lat_scale^2) + 
                                 bs(lon_scale) + (1|Family/Genus) + (1|Source) + 
                                 (1|MarkerName) + fitted(lat_lon_ME), 
                               REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail",
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
  plotResiduals(lat_lon_model_pi_sim_output, mtdna_small_pi$lon_scale)
  
#test for SAC
sim_recalc <- recalculateResiduals(lat_lon_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### abslat & lon model ####
abslat_lon_model_pi <- lmer(logpi ~ range_position + abslat_scale + bs(lon_scale) + 
                              (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                            REML = FALSE, data = mtdna_small_pi, 
                            na.action = "na.fail",
                            control = lmerControl(optimizer = "bobyqa"))
  
## with SAC ##
#create ME eigenvectors
set.seed(8484)
  
#abslat_lon_ME <- ME(logpi ~ range_position + abslat_scale + bs(lon_scale),
#                    data = mtdna_small_pi, listw = pi_wt4)



## abslat lon model with SA ##
abslat_lon_model_pi_SA_RE <- lmer(logpi ~ range_position + abslat_scale + 
                                    bs(lon_scale) + (1|Family/Genus) + 
                                    (1|Source) + (1|MarkerName) + (1|latlonbins), 
                                  REML = FALSE, data = mtdna_small_pi, 
                                  na.action = "na.fail",
                                  control = lmerControl(optimizer = "bobyqa"))

abslat_lon_model_pi_SA_ME <- lmer(logpi ~ range_position + abslat_scale + 
                                    bs(lon_scale) + (1|Family/Genus) + 
                                    (1|Source) + (1|MarkerName) + fitted(abslat_lon_ME), 
                                  REML = FALSE, data = mtdna_small_pi, 
                                  na.action = "na.fail",
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
  plotResiduals(abslat_lon_model_pi_sim_output, mtdna_small_pi$lon_scale)
  
#test for SAC
sim_recalc <- recalculateResiduals(abslat_lon_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

####################################################################################################################

######## Environmental models ########

#### sst mean model ####
sstmean_model_pi <- lmer(logpi ~ range_position + sst.BO_sstmean + (1|Family/Genus) + 
                           (1|Source) + (1|MarkerName), 
                         REML = FALSE, data = mtdna_small_pi, 
                         na.action = "na.fail",
                         control = lmerControl(optimizer = "bobyqa"))
## with SAC ##
#create ME eigenvectors
set.seed(8484)
  
#sstmean_ME <- ME(logpi ~ range_position + sst.BO_sstmean,
#                 data = mtdna_small_pi, listw = pi_wt4)
  
#models
sstmean_model_pi_SA_RE <- lmer(logpi ~ range_position + sst.BO_sstmean + (1|Family/Genus) + 
                                 (1|Source) + (1|MarkerName) + (1|latlonbins), 
                               REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail",
                               control = lmerControl(optimizer = "bobyqa"))

sstmean_model_pi_SA_ME <- lmer(logpi ~ range_position + sst.BO_sstmean + (1|Family/Genus) + 
                                 (1|Source) + (1|MarkerName) + fitted(sstmean_ME), 
                               REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail",
                               control = lmerControl(optimizer = "bobyqa"))
  
#pull p-values
coefs <- data.frame(coef(summary(sstmean_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
sstmean_model_pi_sim_output <- simulateResiduals(fittedModel = sstmean_model_pi, n = 1000, plot = F)
plotQQunif(sstmean_model_pi_sim_output)
plotResiduals(sstmean_model_pi_sim_output)
  plotResiduals(sstmean_model_pi_sim_output, mtdna_small_pi$sst.BO_sstmean)

#test for SAC
sim_recalc <- recalculateResiduals(sstmean_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### sst range model ####
sstrange_model_pi <- lmer(logpi ~ range_position + sst.BO_sstrange + (1|Family/Genus) + 
                            (1|Source) + (1|MarkerName), 
                          REML = FALSE, data = mtdna_small_pi, 
                          na.action = "na.fail",
                          control = lmerControl(optimizer = "bobyqa"))

## with SAC ##
#create ME eigenvectors
set.seed(8484)
  
#sstrange_ME <- ME(logpi ~ range_position + sst.BO_sstrange,
#                   data = mtdna_small_pi, listw = pi_wt4)

#models
sstrange_model_pi_SA_RE <- lmer(logpi ~ range_position + sst.BO_sstrange + (1|Family/Genus) + 
                                  (1|Source) + (1|MarkerName) + (1|latlonbins), 
                                REML = FALSE, data = mtdna_small_pi, 
                                na.action = "na.fail",
                                control = lmerControl(optimizer = "bobyqa"))

sstrange_model_pi_SA_ME <- lmer(logpi ~ range_position + sst.BO_sstrange + (1|Family/Genus) + 
                                  (1|Source) + (1|MarkerName) + fitted(sstrange_ME), 
                                REML = FALSE, data = mtdna_small_pi, 
                                na.action = "na.fail",
                                control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(sstrange_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
sstrange_model_pi_sim_output <- simulateResiduals(fittedModel = sstrange_model_pi, n = 1000, plot = F)
plotQQunif(sstrange_model_pi_sim_output)
plotResiduals(sstrange_model_pi_sim_output)
  plotResiduals(sstrange_model_pi_sim_output, mtdna_small_pi$sst.BO_sstrange)
  
#test for SAC
sim_recalc <- recalculateResiduals(sstrange_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### sst max model ####
sstmax_model_pi <- lmer(logpi ~ range_position + sst.BO_sstmax + (1|Family/Genus) + 
                          (1|Source) + (1|MarkerName), 
                        REML = FALSE, data = mtdna_small_pi, 
                        na.action = "na.fail", 
                        control = lmerControl(optimizer = "bobyqa"))

## with SAC ##
#create ME eigenvectors
set.seed(8484)
  
#sstmax_ME <- ME(logpi ~ range_position + sst.BO_sstmax,
#                data = mtdna_small_pi, listw = pi_wt4)
 
#models 
sstmax_model_pi_SA_RE <- lmer(logpi ~ range_position + sst.BO_sstmax + (1|Family/Genus) + 
                                (1|Source) + (1|MarkerName) + (1|latlonbins), 
                              REML = FALSE, data = mtdna_small_pi, 
                              na.action = "na.fail", 
                              control = lmerControl(optimizer = "bobyqa"))

sstmax_model_pi_SA_ME <- lmer(logpi ~ range_position + sst.BO_sstmax + (1|Family/Genus) + 
                                (1|Source) + (1|MarkerName) + fitted(sstmax_ME), 
                              REML = FALSE, data = mtdna_small_pi, 
                              na.action = "na.fail", 
                              control = lmerControl(optimizer = "bobyqa"))
  
#pull p-values
coefs <- data.frame(coef(summary(sstmax_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs
  
#checking fit with DHARMa
sstmax_model_pi_sim_output <- simulateResiduals(fittedModel = sstmax_model_pi, n = 1000, plot = F)
plotQQunif(sstmax_model_pi_sim_output)
plotResiduals(sstmax_model_pi_sim_output)
  plotResiduals(sstmax_model_pi_sim_output, mtdna_small_pi$sst.BO_sstmax)

#test for SAC
sim_recalc <- recalculateResiduals(sstmax_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### sst min model ####
sstmin_model_pi <- lmer(logpi ~ range_position + sst.BO_sstmin + (1|Family/Genus) + 
                          (1|Source) + (1|MarkerName), 
                        REML = FALSE, data = mtdna_small_pi, 
                        na.action = "na.fail",
                        control = lmerControl(optimizer = "bobyqa"))

## with SAC ##
#create ME eigenvectors
set.seed(8484)
  
#sstmin_ME <- ME(logpi ~ range_position + sst.BO_sstmin,
#                data = mtdna_small_pi, listw = pi_wt4)

#models
sstmin_model_pi_SA_RE <- lmer(logpi ~ range_position + sst.BO_sstmin + (1|Family/Genus) + 
                                (1|Source) + (1|MarkerName) + (1|latlonbins), 
                              REML = FALSE, data = mtdna_small_pi, 
                              na.action = "na.fail",
                              control = lmerControl(optimizer = "bobyqa"))

sstmin_model_pi_SA_ME <- lmer(logpi ~ range_position + sst.BO_sstmin + (1|Family/Genus) + 
                                (1|Source) + (1|MarkerName) + fitted(sstmin_ME), 
                              REML = FALSE, data = mtdna_small_pi, 
                              na.action = "na.fail",
                              control = lmerControl(optimizer = "bobyqa"))
  
#pull p-values
coefs <- data.frame(coef(summary(sstmin_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs
  
#checking fit with DHARMa
sstmin_model_pi_sim_output <- simulateResiduals(fittedModel = sstmin_model_pi, n = 1000, plot = F)
plotQQunif(sstmin_model_pi_sim_output)
plotResiduals(sstmin_model_pi_sim_output)
  plotResiduals(sstmin_model_pi_sim_output, mtdna_small_pi$sst.BO_sstmin)
  
#test for SAC
sim_recalc <- recalculateResiduals(sstmin_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### chloro mean model ####
chlomean_model_pi <- lmer(logpi ~ range_position + logchlomean + I(logchlomean^2) + 
                            (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                          REML = FALSE, data = mtdna_small_pi, 
                          na.action = "na.fail",
                          control = lmerControl(optimizer = "bobyqa"))
  
## with SAC ##
#create ME eigenvectors
set.seed(8484)

#chlomean_ME <- ME(logpi ~ range_position + logchlomean + I(logchlomean^2),
#                  data = mtdna_small_pi, listw = pi_wt4)

#models
chlomean_model_pi_SA_RE <- lmer(logpi ~ range_position + logchlomean + I(logchlomean^2) + 
                                  (1|Family/Genus) + (1|Source) + (1|MarkerName) + (1|latlonbins), 
                                REML = FALSE, data = mtdna_small_pi, 
                                na.action = "na.fail",
                                control = lmerControl(optimizer = "bobyqa"))

chlomean_model_pi_SA_ME <- lmer(logpi ~ range_position + logchlomean + I(logchlomean^2) + 
                                  (1|Family/Genus) + (1|Source) + (1|MarkerName) + fitted(chlomean_ME), 
                             REML = FALSE, data = mtdna_small_pi, 
                             na.action = "na.fail",
                             control = lmerControl(optimizer = "bobyqa"))

#pull p-values
coefs <- data.frame(coef(summary(chlomean_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs

#checking fit with DHARMa
chlomean_model_pi_sim_output <- simulateResiduals(fittedModel = chlomean_model_pi, n = 1000, plot = F)
plotQQunif(chlomean_model_pi_sim_output)
plotResiduals(chlomean_model_pi_sim_output)
  plotResiduals(chlomean_model_pi_sim_output, mtdna_small_pi$logchlomean)

#test for SAC
sim_recalc <- recalculateResiduals(chlomean_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
#### chloro range model ####
chlorange_model_pi <- lmer(logpi ~ range_position + logchlorange + I(logchlorange^2) + 
                             (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                           REML = FALSE, data = mtdna_small_pi, 
                           na.action = "na.fail",
                           control = lmerControl(optimizer = "bobyqa"))
  
## with SAC ##
#create ME eigenvectors
set.seed(8484)

#chlorange_ME <- ME(logpi ~ range_position + logchlorange + I(logchlorange^2),
#                    data = mtdna_small_pi, listw = pi_wt4)

#models
chlorange_model_pi_SA_RE <- lmer(logpi ~ range_position + logchlorange + I(logchlorange^2) + 
                                   (1|Family/Genus) + (1|Source) + (1|MarkerName) + (1|latlonbins), 
                              REML = FALSE, data = mtdna_small_pi, 
                              na.action = "na.fail",
                              control = lmerControl(optimizer = "bobyqa"))

chlorange_model_pi_SA_ME <- lmer(logpi ~ range_position + logchlorange + I(logchlorange^2) + 
                                   (1|Family/Genus) + (1|Source) + (1|MarkerName) + fitted(chlorange_ME), 
                                 REML = FALSE, data = mtdna_small_pi, 
                                 na.action = "na.fail",
                                 control = lmerControl(optimizer = "bobyqa"))
  
#pull p-values
coefs <- data.frame(coef(summary(chlorange_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs
  
#checking fit with DHARMa
chlorange_model_pi_sim_output <- simulateResiduals(fittedModel = chlorange_model_pi, n = 1000, plot = F)
plotQQunif(chlorange_model_pi_sim_output)
plotResiduals(chlorange_model_pi_sim_output)
  plotResiduals(chlorange_model_pi_sim_output, mtdna_small_pi$logchlorange)
  
#test for SAC
sim_recalc <- recalculateResiduals(chlorange_model_pi_sim_output, 
                                    group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### chloro max model ####
chlomax_model_pi <- lmer(logpi ~ range_position + logchlomax + I(logchlomax^2) + 
                           (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                         REML = FALSE, data = mtdna_small_pi, 
                         na.action = "na.fail",
                         control = lmerControl(optimizer = "bobyqa"))
## with SAC ##
#create ME eigenvectors
set.seed(8484)

#chlomax_ME <- ME(logpi ~ range_position + logchlomax + I(logchlomax^2),
#                 data = mtdna_small_pi, listw = pi_wt4)

#models
chlomax_model_pi_SA_RE <- lmer(logpi ~ range_position + logchlomax + I(logchlomax^2) + 
                                 (1|Family/Genus) + (1|Source) + (1|MarkerName) + (1|latlonbins), 
                            REML = FALSE, data = mtdna_small_pi, 
                            na.action = "na.fail",
                            control = lmerControl(optimizer = "bobyqa"))

chlomax_model_pi_SA_ME <- lmer(logpi ~ range_position + logchlomax + I(logchlomax^2) + 
                                 (1|Family/Genus) + (1|Source) + (1|MarkerName) + fitted(chlomax_ME), 
                               REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail",
                               control = lmerControl(optimizer = "bobyqa"))
  
#pull p-values
coefs <- data.frame(coef(summary(chlomax_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs
  
#checking fit with DHARMa
chlomax_model_pi_sim_output <- simulateResiduals(fittedModel = chlomax_model_pi, n = 1000, plot = F)
plotQQunif(chlomax_model_pi_sim_output)
plotResiduals(chlomax_model_pi_sim_output)
  plotResiduals(chlomax_model_pi_sim_output, mtdna_small_pi$logchlomax)
  
#test for SAC
sim_recalc <- recalculateResiduals(chlomax_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
  
#### chloro min model ####
chlomin_model_pi <- lmer(logpi ~ range_position + logchlomin + I(logchlomin^2) + 
                           (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                          REML = FALSE, data = mtdna_small_pi, 
                          na.action = "na.fail",
                          control = lmerControl(optimizer = "bobyqa"))
## with SAC ##
#create ME eigenvectors
set.seed(8484)

#chlomin_ME <- ME(logpi ~ range_position + logchlomin + I(logchlomin^2),
#                 data = mtdna_small_pi, listw = pi_wt4)

#models
chlomin_model_pi_SA_RE <- lmer(logpi ~ range_position + logchlomin + I(logchlomin^2) + 
                                 (1|Family/Genus) + (1|Source) + (1|MarkerName) + (1|latlonbins), 
                               REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail",
                               control = lmerControl(optimizer = "bobyqa"))

chlomin_model_pi_SA_ME <- lmer(logpi ~ range_position + logchlomin + I(logchlomin^2) + 
                                 (1|Family/Genus) + (1|Source) + (1|MarkerName) + fitted(chlomin_ME), 
                               REML = FALSE, data = mtdna_small_pi, 
                                na.action = "na.fail",
                                control = lmerControl(optimizer = "bobyqa"))
  
#pull p-values
coefs <- data.frame(coef(summary(chlomin_model_pi)))
  coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
  coefs
  
#checking fit with DHARMa
chlomin_model_pi_sim_output <- simulateResiduals(fittedModel = chlomin_model_pi, n = 1000, plot = F)
plotQQunif(chlomin_model_pi_sim_output)
plotResiduals(chlomin_model_pi_sim_output)
  plotResiduals(chlomin_model_pi_sim_output, mtdna_small_pi$logchlomin)
  
#test for SAC
sim_recalc <- recalculateResiduals(chlomin_model_pi_sim_output, 
                                   group = mtdna_small_pi$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
######################################################################################################################
  
######## Save all ME ########
#save all ME as an Rdata file so don't have to recreate every time
#save(list = c("null_ME", "lat_ME", "abslat_ME", "lon_ME", "lat_lon_ME",
#              "abslat_lon_ME", "sstmean_ME", "sstrange_ME", "sstmax_ME", 
#              "sstmin_ME", "chlomean_ME", "chlorange_ME", "chlomax_ME", "chlomin_ME"), 
#       file = "mtdna_pi_ME.Rdata")
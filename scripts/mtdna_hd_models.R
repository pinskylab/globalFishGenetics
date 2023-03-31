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
library(spdep)
library(spatialreg)

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)
mtdna_env <- read.csv("output/mtdna_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)
load("mtdna_Hd_ME.Rdata") #for ME weights, so don't have to re-run every time

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
mtdna_small_hd$lon_scale <- as.numeric(scale(mtdna_small_hd$lon))

#convert lon to radians
mtdna_small_hd$lon_360 <- mtdna_small_hd$lon + 180 #convert (-180,180) to (0,360)
  mtdna_small_hd$lon_rad <- (2*pi*mtdna_small_hd$lon_360)/360

#### Calculate environmental variables ####
## log transform chlorophyll A ##
#subset to only those with chloroA data
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlorange != 0)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomax != 0)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomin != 0)

mtdna_small_hd$logchlomean <- log10(mtdna_small_hd$chloroA.BO_chlomean)
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

  #### Calculate weights ####
  #to account for spatial autocorrelation
  
  #create matrix of geographic coordinates
  lonlat_mat <- cbind(mtdna_small_hd$lon, mtdna_small_hd$lat)
  
  #groups into nearest neighbor groups
  hd_nb4 <- knearneigh(x = lonlat_mat, k = 4, longlat = TRUE) #creates matrix with the indices of points belonging to the sets of k nearest neighbors
  hd_nb4 <- knn2nb(hd_nb4) #converts knn object to neighbors list (e.g. these indices in matrix group together spatially as nearest neighbors)
  hd_nb4 <- make.sym.nb(hd_nb4) #checks for symmetry/transitivity
  #calculate weights
  hd_wt4 <- nb2listw(hd_nb4, style = "W") #creates spatial weights
  
#### Create coordinate dataframe for SAC tests ####
#need to create grouping factor for observations with identical lat/lon
mtdna_small_hd$latlonbins <- paste(round(mtdna_small_hd$lat, 0), "_", round(mtdna_small_hd$lon, 0), sep = "") #grouping factor for 1x1 lat/lon bins
mtdna_small_hd$coords <- paste(mtdna_small_hd$lat, "_", mtdna_small_hd$lon, sep = "") #grouping factor for residuals
  coords <- as.data.frame(unique(mtdna_small_hd$coords))
  colnames(coords) <- c("coords")

#pull out unique coordinate for SAC test
coords_unique <- coords %>% separate(coords, sep = "_", c("lat_unique", "lon_unique"))
  x_unique <- coords_unique$lat_unique
  y_unique <- coords_unique$lon_unique
  
  ac <- autocov_dist(mtdna_small_hd$He, lonlat_mat, nbs = 100000, style = "B", longlat = TRUE)
    mtdna_small_hd$ac <- ac
    mtdna_small_hd$ac[mtdna_small_hd$ac == "Inf"] <- 1
    
  
#########################################################################################################################

######## Null model ########

## null model ##
null_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                         (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                       data = mtdna_small_hd, family = binomial,
                       na.action = "na.fail", nAGQ = 0, #nAGQ = 0 & #nAGQ = 1 qualitatively the same, so using 0 bc conveges much faster
                       control = glmerControl(optimizer = "bobyqa"))

## null model with SA ##
null_model_hd_SA_AC <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                            (1|Family/Genus) + (1|Source) + (1|MarkerName) + (1|latlonbins), 
                          data = mtdna_small_hd, family = binomial,
                          na.action = "na.fail", nAGQ = 0,
                          control = glmerControl(optimizer = "bobyqa"))

null_model_hd_SA_ME <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                               (1|Family/Genus) + (1|Source) + (1|MarkerName) + fitted(null_ME), 
                             data = mtdna_small_hd, family = binomial,
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa"))

#check fit with DHARMa
null_model_hd_sim <- simulateResiduals(fittedModel = null_model_hd_SA_AC, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(null_model_hd_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(null_model_hd_sim) #residuals against predicted value -- looking for uniformity
  plotResiduals(null_model_hd_sim, mtdna_small_hd$bp_scale)
  plotResiduals(null_model_hd_sim, mtdna_small_hd$range_position)

 resids <- residuals(null_model_hd_sim)
  
#test for SAC
sim_recalc <- recalculateResiduals(null_model_hd_sim, group = mtdna_small_hd$coords) #need to group same lat/lons together
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
###################################################################################################################

######## Latitude & longitude models ########

#### lat model ####
lat_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                           lat_scale + I(lat_scale^2) + (1|Family/Genus) + 
                           (1|Source) + (1|MarkerName), 
                         data = mtdna_small_hd, family = binomial,
                         na.action = "na.fail", nAGQ = 0,
                         control = glmerControl(optimizer = "bobyqa"))
  
## lat model with SA ##
lat_model_hd_SA_RE <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                           lat_scale + I(lat_scale^2) + (1|Family/Genus) + 
                           (1|Source) + (1|MarkerName) + (1|latlonbins), 
                         data = mtdna_small_hd, family = binomial,
                         na.action = "na.fail", nAGQ = 0,
                         control = glmerControl(optimizer = "bobyqa"))

lat_model_hd_SA_ME <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                              lat_scale + I(lat_scale^2) + (1|Family/Genus) + 
                              (1|Source) + (1|MarkerName) + fitted(lat_ME), 
                            data = mtdna_small_hd, family = binomial,
                            na.action = "na.fail", nAGQ = 0,
                            control = glmerControl(optimizer = "bobyqa"))
  
#check fit with DHARMa
lat_model_hd_sim <- simulateResiduals(fittedModel = lat_model_hd, n = 1000, plot = F)
plotQQunif(lat_model_hd_sim)
plotResiduals(lat_model_hd_sim)
  plotResiduals(lat_model_hd_sim, mtdna_small_hd$lat_scale)
  
#test for SAC
sim_recalc <- recalculateResiduals(lat_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### abslat model ####
abslat_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                                abslat_scale + (1|Family/Genus) + (1|Source) + 
                                (1|MarkerName), 
                              data = mtdna_small_hd, family = binomial,
                              na.action = "na.fail", nAGQ = 0,
                              control = glmerControl(optimizer = "bobyqa"))
## abslat model with SA ##
abslat_model_hd_SA_RE <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                              abslat_scale + (1|Family/Genus) + (1|Source) + 
                              (1|MarkerName) + (1|latlonbins),
                            data = mtdna_small_hd, family = binomial,
                            na.action = "na.fail", nAGQ = 0,
                            control = glmerControl(optimizer = "bobyqa"))

abslat_model_hd_SA_ME <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                                 abslat_scale + (1|Family/Genus) + (1|Source) + 
                                 (1|MarkerName) + fitted(abslat_ME),
                               data = mtdna_small_hd, family = binomial,
                               na.action = "na.fail", nAGQ = 0,
                               control = glmerControl(optimizer = "bobyqa"))
  
#check fit with DHARMa
abslat_model_hd_sim <- simulateResiduals(fittedModel = abslat_model_hd, n = 1000, plot = F)
plotQQunif(abslat_model_hd_sim)
plotResiduals(abslat_model_hd_sim)
  plotResiduals(abslat_model_hd_sim, mtdna_small_hd$abslat_scale)

#test for SAC
sim_recalc <- recalculateResiduals(abslat_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
#### lon model ####
lon_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                          bs(lon_scale) + (1|Family/Genus) + (1|Source) + 
                          (1|MarkerName),
                        family = binomial, data = mtdna_small_hd, 
                        na.action = "na.fail", nAGQ = 0,
                        control = glmerControl(optimizer = "bobyqa"))

## lon model with SA ##
lon_model_hd_SA_RE <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                          bs(lon_scale) + (1|Family/Genus) + (1|Source) + 
                          (1|MarkerName) + (1|latlonbins),
                        family = binomial, data = mtdna_small_hd, 
                        na.action = "na.fail", nAGQ = 0,
                        control = glmerControl(optimizer = "bobyqa"))

lon_model_hd_SA_ME <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                              bs(lon_scale) + (1|Family/Genus) + (1|Source) + 
                              (1|MarkerName) + fitted(lon_ME),
                            family = binomial, data = mtdna_small_hd, 
                            na.action = "na.fail", nAGQ = 0,
                            control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
lon_model_hd_sim <- simulateResiduals(fittedModel = lon_model_hd_SA, n = 1000, plot = F)
plotQQunif(lon_model_hd_sim)
plotResiduals(lon_model_hd_sim)
  plotResiduals(lon_model_hd_sim, mtdna_small_hd$lon_scale)
  

  sim_recalc <- recalculateResiduals(lon_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### lat & lon model ####
lat_lon_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + bs(lon_scale) + 
                               lat_scale + I(lat_scale^2) +(1|Family/Genus) + (1|Source) + 
                               (1|MarkerName),
                             family = binomial, data = mtdna_small_hd, 
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa")) 

## lat lon model with SA ## 
#account for and remove spatial autocorrelation from residuals
set.seed(8484)

#create eigenvectors  
lat_lon_ME <- ME(cbind(success, failure) ~ bp_scale + range_position + 
                   lat_scale + I(lat_scale^2) + bs(lon_scale),
               data = mtdna_small_hd, family = binomial, listw = hd_wt4)
  
#model with SA
lat_lon_model_hd_SA_RE <- glmer(cbind(success, failure) ~ bp_scale + range_position + bs(lon_scale) + 
                               lat_scale + I(lat_scale^2) +(1|Family/Genus) + (1|Source) + 
                               (1|MarkerName) + (1|latlonbins),
                             family = binomial, data = mtdna_small_hd, 
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa"))

lat_lon_model_hd_SA_ME <- glmer(cbind(success, failure) ~ bp_scale + range_position + bs(lon_scale) + 
                                  lat_scale + I(lat_scale^2) +(1|Family/Genus) + (1|Source) + 
                                  (1|MarkerName) + fitted(lat_lon_ME),
                                family = binomial, data = mtdna_small_hd, 
                                na.action = "na.fail", nAGQ = 0,
                                control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
lat_lon_model_hd_sim <- simulateResiduals(fittedModel = lat_lon_model_hd_SA, n = 1000, plot = F)
plotQQunif(lat_lon_model_hd_sim)
  plotResiduals(lat_lon_model_hd_sim)
  plotResiduals(lat_lon_model_hd_sim, mtdna_small_hd$lon_scale)
  plotResiduals(lat_lon_model_hd_sim, mtdna_small_hd$lat_scale)
  
  sim_recalc <- recalculateResiduals(lat_lon_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### abslat & lon model ###
abslat_lon_model_hd  <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                                        bs(lon_scale) + abslat_scale + (1|Family/Genus) + 
                                        (1|Source) + (1|MarkerName), 
                                      family = binomial, data = mtdna_small_hd, 
                                      na.action = "na.fail", nAGQ = 0, 
                                      control = glmerControl(optimizer = "bobyqa"))

## abslat lon model with SA ##
#account for and remove spatial autocorrelation from residuals
set.seed(8484)

#create eigenvectors
abslat_lon_ME <- ME(cbind(success, failure) ~ bp_scale + range_position + 
                     abslat_scale + bs(lon_scale),
                   data = mtdna_small_hd, family = binomial, listw = hd_wt4)
  
#model with SA
abslat_lon_model_hd_SA_RE  <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                                      bs(lon_scale) + abslat_scale + (1|Family/Genus) + 
                                      (1|Source) + (1|MarkerName) + (1|latlonbins), 
                                    family = binomial, data = mtdna_small_hd, 
                                    na.action = "na.fail", nAGQ = 0, 
                                    control = glmerControl(optimizer = "bobyqa"))

abslat_lon_model_hd_SA_ME  <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                                   bs(lon_scale) + abslat_scale + (1|Family/Genus) + 
                                   (1|Source) + (1|MarkerName) + fitted(abslat_lon_ME), 
                                 family = binomial, data = mtdna_small_hd, 
                                 na.action = "na.fail", nAGQ = 0, 
                                 control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
abslat_lon_model_hd_sim <- simulateResiduals(fittedModel = abslat_lon_model_hd_SA, n = 1000, plot = F)
plotQQunif(abslat_lon_model_hd_sim)
plotResiduals(abslat_lon_model_hd_sim)
  plotResiduals(abslat_lon_model_hd_sim, mtdna_small_hd$lon_scale)
  plotResiduals(abslat_lon_model_hd_sim, mtdna_small_hd$abslat_scale)

###################################################################################################################

######## Environmental models ########

#### sst mean model ####
#account for and remove spatial autocorrelation from residuals
sstmean_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                              sst.BO_sstmean + (1|Family/Genus) + (1|Source) + 
                              (1|MarkerName), 
                            family = binomial, data = mtdna_small_hd, 
                            na.action = "na.fail", nAGQ = 0,
                            control = glmerControl(optimizer = "bobyqa"))
set.seed(8484)
  
sstmean_ME <- ME(cbind(success, failure) ~ bp_scale + range_position + sst.BO_sstmean,
                      data = mtdna_small_hd, family = binomial, listw = hd_wt4)

## sstmean model with SA ##
sstmean_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                               sst.BO_sstmean + (1|Family/Genus) + (1|Source) + 
                               (1|MarkerName), 
                             family = binomial, data = mtdna_small_hd, 
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa"))

## sstmean model with SA ##
sstmean_model_hd_SA <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                               sst.BO_sstmean + (1|Family/Genus) + (1|Source) + 
                               (1|MarkerName) + fitted(sstmean_ME), 
                             family = binomial, data = mtdna_small_hd, 
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
sstmean_model_hd_sim <- simulateResiduals(fittedModel = sstmean_model_hd_SA, n = 1000, plot = F)
plotQQunif(sstmean_model_hd_sim)
plotResiduals(sstmean_model_hd_sim)
  plotResiduals(sstmean_model_hd_sim, mtdna_small_hd$sst.BO_sstmean)
  
  sim_recalc <- recalculateResiduals(sstmean_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### sst range model ####
sstrange_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                                  sst.BO_sstrange + (1|Family/Genus) + (1|Source) + 
                                  (1|MarkerName), 
                                family = binomial, data = mtdna_small_hd, 
                                na.action = "na.fail", nAGQ = 0,
                                control = glmerControl(optimizer = "bobyqa"))
#account for and remove spatial autocorrelation from residuals
set.seed(8484)
  
sstrange_ME <- ME(cbind(success, failure) ~ bp_scale + range_position + sst.BO_sstrange,
                   data = mtdna_small_hd, family = binomial, listw = hd_wt4)
  
## sstrange model with SA ##
sstrange_model_hd_SA <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                                sst.BO_sstrange + (1|Family/Genus) + (1|Source) + 
                                (1|MarkerName) + fitted(sstrange_ME), 
                              family = binomial, data = mtdna_small_hd, 
                              na.action = "na.fail", nAGQ = 0,
                              control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
sstrange_model_hd_sim <- simulateResiduals(fittedModel = sstrange_model_hd_SA, n = 1000, plot = F)
plotQQunif(sstrange_model_hd_sim)
plotResiduals(sstrange_model_hd_sim)
  plotResiduals(sstrange_model_hd_sim, mtdna_small_hd$sst.BO_sstrange)

#### sst max model ####
sstmax_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                              sst.BO_sstmax + (1|Family/Genus) + (1|Source) + 
                                (1|MarkerName), 
                              family = binomial, data = mtdna_small_hd, 
                              na.action = "na.fail", nAGQ = 0,
                              control = glmerControl(optimizer = "bobyqa"))
#account for and remove spatial autocorrelation from residuals
set.seed(8484)
  
sstmax_ME <- ME(cbind(success, failure) ~ bp_scale + range_position + sst.BO_sstmax,
                    data = mtdna_small_hd, family = binomial, listw = hd_wt4)
  
## sstmax model with SA ##
sstmax_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                              sst.BO_sstmax + (1|Family/Genus) + (1|Source) + 
                              (1|MarkerName) + fitted(sstmax_ME), 
                            family = binomial, data = mtdna_small_hd, 
                            na.action = "na.fail", nAGQ = 0,
                            control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
sstmax_model_hd_sim <- simulateResiduals(fittedModel = sstmax_model_hd_SA, n = 1000, plot = F)
plotQQunif(sstmax_model_hd_sim)
plotResiduals(sstmax_model_hd_sim)
  plotResiduals(sstmax_model_hd_sim, mtdna_small_hd$sst.BO_sstmax)

#### sst min model ####
  sstmin_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                                sst.BO_sstmin + (1|Family/Genus) + (1|Source) + 
                                (1|MarkerName), 
                              family = binomial, data = mtdna_small_hd, 
                              na.action = "na.fail", nAGQ = 0,
                              control = glmerControl(optimizer = "bobyqa"))
#account for and remove spatial autocorrelation from residuals
set.seed(8484)
  
sstmin_ME <- ME(cbind(success, failure) ~ bp_scale + range_position + sst.BO_sstmin,
                  data = mtdna_small_hd, family = binomial, listw = hd_wt4)
  
## sstmin model with SA ##
sstmin_model_hd_SA <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                              sst.BO_sstmin + (1|Family/Genus) + (1|Source) + 
                              (1|MarkerName) + fitted(sstmin_ME), 
                            family = binomial, data = mtdna_small_hd, 
                            na.action = "na.fail", nAGQ = 0,
                            control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
sstmin_model_hd_sim <- simulateResiduals(fittedModel = sstmin_model_hd_SA, n = 1000, plot = F)
plotQQunif(sstmin_model_hd_sim)
plotResiduals(sstmin_model_hd_sim)
  plotResiduals(sstmin_model_hd_sim, mtdna_small_hd$sst.BO_sstmin)

#### chloro mean model ####
  chlomean_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomean + 
                                  I(logchlomean^2) + (1|Family/Genus) + (1|Source) + 
                                  (1|MarkerName), 
                                family = binomial, data = mtdna_small_hd, 
                                na.action = "na.fail", nAGQ = 0,
                                control = glmerControl(optimizer = "bobyqa"))
  
#account for and remove spatial autocorrelation from residuals
set.seed(8484)
  
chlomean_ME <- ME(cbind(success, failure) ~ bp_scale + range_position + 
                    logchlomean + I(logchlomean^2),
                  data = mtdna_small_hd, family = binomial, listw = hd_wt4)
  
## chlomean model with SA ##
chlomean_model_hd_SA <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomean + 
                                I(logchlomean^2) + (1|Family/Genus) + (1|Source) + 
                                (1|MarkerName) + fitted(chlomean_ME), 
                              family = binomial, data = mtdna_small_hd, 
                              na.action = "na.fail", nAGQ = 0,
                              control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
chlomean_model_hd_sim <- simulateResiduals(fittedModel = chlomean_model_hd_SA, n = 1000, plot = F)
plotQQunif(chlomean_model_hd_sim)
plotResiduals(chlomean_model_hd_sim)
  plotResiduals(chlomean_model_hd_sim, mtdna_small_hd$logchlomean)

##### chloro range model ####
  chlorange_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlorange + 
                                   I(logchlorange^2) + (1|Family/Genus) + (1|Source) + 
                                   (1|MarkerName), 
                                 family = binomial, data = mtdna_small_hd, 
                                 na.action = "na.fail", nAGQ = 0,
                                 control = glmerControl(optimizer = "bobyqa"))
  
#account for and remove spatial autocorrelation from residuals
set.seed(8484)
  
chlorange_ME <- ME(cbind(success, failure) ~ bp_scale + range_position + 
                      logchlorange + I(logchlorange^2),
                    data = mtdna_small_hd, family = binomial, listw = hd_wt4)
  
## chlorange model with SA ##
chlorange_model_hd_SA <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlorange + 
                                 I(logchlorange^2) + (1|Family/Genus) + (1|Source) + 
                                (1|MarkerName) + fitted(chlorange_ME), 
                               family = binomial, data = mtdna_small_hd, 
                               na.action = "na.fail", nAGQ = 0,
                               control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
chlorange_model_hd_sim <- simulateResiduals(fittedModel = chlorange_model_hd_SA, n = 1000, plot = F)
plotQQunif(chlorange_model_hd_sim)
  plotResiduals(chlorange_model_hd_sim)
  plotResiduals(chlorange_model_hd_sim, mtdna_small_hd$logchlorange)
  
  sim_recalc <- recalculateResiduals(chlorange_model_hd_sim, group = mtdna_small_hd$coords)
  testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)

#### chloro max model ####
   chlomax_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomax + 
                              I(logchlomax^2) + (1|Family/Genus) + (1|Source) + 
                              (1|MarkerName), 
                            family = binomial, data = mtdna_small_hd, 
                            na.action = "na.fail", nAGQ = 0,
                            control = glmerControl(optimizer = "bobyqa"))
  
#account for and remove spatial autocorrelation from residuals
set.seed(8484)
  
chlomax_ME <- ME(cbind(success, failure) ~ bp_scale + range_position + 
                       logchlomax + I(logchlomax^2),
                     data = mtdna_small_hd, family = binomial, listw = hd_wt4)
  
## chlomax model with SA ##
chlomax_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomax + 
                               I(logchlomax^2) + (1|Family/Genus) + (1|Source) + 
                               (1|MarkerName), 
                             family = binomial, data = mtdna_small_hd, 
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa"))

mtdna_test <- mtdna_small_hd[1:20,]
mtdna_small_hd$coords <- paste(round(mtdna_small_hd$lat, 0), "_", round(mtdna_small_hd$lon, 0), sep = "")
coords <- as.data.frame(unique(mtdna_small_hd$coords))
  colnames(coords) <- c("coords")
coords_unique <- coords %>% separate(coords, sep = "_", c("lat_unique", "lon_unique"))
x_unique <- coords_unique$lat_unique
y_unique <- coords_unique$lon_unique


#checking fit with DHARMa
chlomax_model_hd_sim <- simulateResiduals(fittedModel = chlomax_model_hd, n = 1000, plot = F)
plotQQunif(chlomax_model_hd_sim)
plotResiduals(chlomax_model_hd_sim)
  plotResiduals(chlomax_model_hd_sim, mtdna_small_hd$logchlomax)

mtdna_test <- mtdna_small_hd[1:20,]
  
sim_recalc <- recalculateResiduals(chlomin_model_hd_sim, group = mtdna_small_hd$coords)
testSpatialAutocorrelation(sim_recalc, x = x_unique, y = y_unique)
  
#### chloro min model ####
chlomin_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomin + 
                               I(logchlomin^2) + (1|Family/Genus) + (1|Source) + 
                               (1|MarkerName), 
                             family = binomial, data = mtdna_small_hd, 
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa"))

#account for and remove spatial autocorrelation from residuals
set.seed(8484)
  
chlomin_ME <- ME(cbind(success, failure) ~ bp_scale + range_position + 
                     logchlomin + I(logchlomin^2),
                   data = mtdna_small_hd, family = binomial, listw = hd_wt4)
  
## chlomin model with SA ##
chlomin_model_hd_SA <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomin + 
                               I(logchlomin^2) + (1|Family/Genus) + (1|Source) + 
                               (1|MarkerName) + fitted(chlomin_ME), 
                             family = binomial, data = mtdna_small_hd, 
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa"))
  
#checking fit with DHARMa
chlomin_model_hd_sim <- simulateResiduals(fittedModel = chlomin_model_hd_SA, n = 1000, plot = F)
plotQQunif(chlomin_model_hd_sim)
plotResiduals(chlomin_model_hd_sim)
  plotResiduals(chlomin_model_hd_sim, mtdna_small_hd$logchlomin)
  
######################################################################################################################
  
######## Save all ME ########
#save all ME as an Rdata file so don't have to recreate every time
save(list = c("null_ME", "lat_ME", "abslat_ME", "lon_ME", "lat_lon_ME", 
              "abslat_lon_ME", "sstmean_ME", "sstrange_ME", "sstmax_ME", 
              "sstmin_ME", "chlomean_ME", "chlorange_ME", "chlomax_ME", "chlomin_ME"), 
     file = "mtdna_Hd_ME.Rdata")
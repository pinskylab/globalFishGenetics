################################################## Bootstrap pi models #########################################################

#bootstrapping (10,000 iterations/model) to get 95% CIs for coefficients
#do relationships change (e.g. - to +?)
#how sensitive are relationships/models to data input?

#set up to run on Wahab (Old Dominion University HPC)

################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(lme4) #v.1.1-31
library(splines) #v.4.2.2

#functions
boot_coef <- function(data){ #to pull out fixed effects and sigma for bootstrapping
  s <- sigma(data)
  c(getME(data, "fixef"), sigma = s)
}

#read in data
mtdna <- read.csv("mtdna_assembled.csv", stringsAsFactors = FALSE)
mtdna_env <- read.csv("mtdna_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("spp_combined_info.csv", stringsAsFactors = FALSE)

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

#############################################################################################################################

######## Null model ########
  
null_model_pi <- lmer(logpi ~ range_position + (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                      REML = FALSE, data = mtdna_small_pi, 
                      na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
  
## bootstrap confidence interval for coefficients ##  
pi_null_boot <- bootMer(null_model_pi, FUN = boot_coef, nsim = 10000)
`boot_coef_null <- as.data.frame(pi_null_boot$t) 
  
#pull out CIs
lower_ci_null <- as.data.frame(apply(pi_null_boot$t, 2, quantile, 0.025))
upper_ci_null<- as.data.frame(apply(pi_null_boot$t, 2, quantile, 0.975))
  
cis_null <- cbind(lower_ci_null, upper_ci_null)
  cis_null$fixef <- rownames(cis_null)
  cis_null$model <- "null"
  colnames(cis_null) <- c("2.5ci", "97.5ci", "fixef", "model")
  
################################################################################################################
  
######### Latitude & longitude models ########

#### lat model ####
lat_model_pi <- lmer(logpi ~ range_position + lat_scale + I(lat_scale^2) + (1|Family/Genus) +   
                       (1|Source) + (1|MarkerName), 
                     REML = FALSE, data = mtdna_small_pi, 
                     na.action = "na.fail", 
                     control = lmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
pi_lat_boot <- bootMer(lat_model_pi, FUN = boot_coef, nsim = 10000)
  boot_coef_lat <- as.data.frame(pi_lat_boot$t) 

#pull out CIs
lower_ci_lat <- as.data.frame(apply(pi_lat_boot$t, 2, quantile, 0.025))
upper_ci_lat <- as.data.frame(apply(pi_lat_boot$t, 2, quantile, 0.975))

cis_lat <- cbind(lower_ci_lat, upper_ci_lat)
  cis_lat$fixef <- rownames(cis_lat)
  cis_lat$model <- "lat"
  colnames(cis_lat) <- c("2.5ci", "97.5ci", "fixef", "model")

#### abslat model ####
abslat_model_pi <- lmer(logpi ~ range_position + abslat_scale + (1|Family/Genus) + 
                          (1|Source) + (1|MarkerName), 
                        REML = FALSE, data = mtdna_small_pi, 
                        na.action = "na.fail", 
                        control = lmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
pi_abslat_boot <- bootMer(abslat_model_pi, FUN = boot_coef, nsim = 10000)
  boot_coef_abslat <- as.data.frame(pi_abslat_boot$t) 

#pull out CIs
lower_ci_abslat <- as.data.frame(apply(pi_abslat_boot$t, 2, quantile, 0.025))
upper_ci_abslat <- as.data.frame(apply(pi_abslat_boot$t, 2, quantile, 0.975))

cis_abslat <- cbind(lower_ci_abslat, upper_ci_abslat)
  cis_abslat$fixef <- rownames(cis_abslat)
  cis_abslat$model <- "abslat"
  colnames(cis_abslat) <- c("2.5ci", "97.5ci", "fixef", "model")

#### lon model ####
lon_model_pi <- lmer(logpi ~ range_position + bs(lon_scale) + 
                       (1|Family/Genus) + (1|Source) + (1|MarkerName),
                     REML = FALSE, data = mtdna_small_pi, 
                     na.action = "na.fail", 
                     control = lmerControl(optimizer = "bobyqa"))
  
## bootstrap confidence interval for coefficients ##  
pi_lon_boot <- bootMer(lon_model_pi, FUN = boot_coef, nsim = 10000)
  boot_coef_lon <- as.data.frame(pi_lon_boot$t) 

#pull out CIs
lower_ci_lon <- as.data.frame(apply(pi_lon_boot$t, 2, quantile, 0.025))
upper_ci_lon <- as.data.frame(apply(pi_lon_boot$t, 2, quantile, 0.975))

cis_lon <- cbind(lower_ci_lon, upper_ci_lon)
  cis_lon$fixef <- rownames(cis_lon)
  cis_lon$model <- "lon"
  colnames(cis_lon) <- c("2.5ci", "97.5ci", "fixef", "model")
  
###############################################################################################################################

######## Environmental models ########

#### sst mean model ####
mtdna_pi_linear_sstmean <- lmer(logpi ~ range_position + sst.BO_sstmean + 
                                  (1|Family/Genus) + (1|Source) + (1|MarkerName),
                                REML = FALSE, data = mtdna_small_pi, 
                                na.action = "na.fail", 
                                control = lmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
pi_sstmean_boot <- bootMer(mtdna_pi_linear_sstmean, FUN = boot_coef, nsim = 10000)
  boot_coef_sstmean <- as.data.frame(pi_sstmean_boot$t) 

#pull out CIs
lower_ci_sstmean <- as.data.frame(apply(pi_sstmean_boot$t, 2, quantile, 0.025))
upper_ci_sstmean <- as.data.frame(apply(pi_sstmean_boot$t, 2, quantile, 0.975))

cis_sstmean <- cbind(lower_ci_sstmean, upper_ci_sstmean)
  cis_sstmean$fixef <- rownames(cis_sstmean)
  cis_sstmean$model <- "sstmean"
  colnames(cis_sstmean) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### sst range model ####
mtdna_pi_linear_sstrange <- lmer(logpi ~ range_position + sst.BO_sstrange + 
                                   (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                                REML = FALSE, data = mtdna_small_pi, 
                                na.action = "na.fail", 
                                control = lmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
pi_sstrange_boot <- bootMer(mtdna_pi_linear_sstrange, FUN = boot_coef, nsim = 10000)
boot_coef_sstrange <- as.data.frame(pi_sstrange_boot$t) 

#pull out CIs
lower_ci_sstrange <- as.data.frame(apply(pi_sstrange_boot$t, 2, quantile, 0.025))
upper_ci_sstrange <- as.data.frame(apply(pi_sstrange_boot$t, 2, quantile, 0.975))

cis_sstrange <- cbind(lower_ci_sstrange, upper_ci_sstrange)
  cis_sstrange$fixef <- rownames(cis_sstrange)
  cis_sstrange$model <- "sstrange"
  colnames(cis_sstrange) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### sst max model ####
mtdna_pi_linear_sstmax <- lmer(logpi ~ range_position + sst.BO_sstmax + 
                                 (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                               REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", 
                               control = lmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
pi_sstmax_boot <- bootMer(mtdna_pi_linear_sstmax, FUN = boot_coef, nsim = 10000)
  boot_coef_sstmax <- as.data.frame(pi_sstmax_boot$t) 

#pull out CIs
lower_ci_sstmax <- as.data.frame(apply(pi_sstmax_boot$t, 2, quantile, 0.025))
upper_ci_sstmax <- as.data.frame(apply(pi_sstmax_boot$t, 2, quantile, 0.975))

cis_sstmax <- cbind(lower_ci_sstmax, upper_ci_sstmax)
  cis_sstmax$fixef <- rownames(cis_sstmax)
  cis_sstmax$model <- "sstmax"
  colnames(cis_sstmax) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### sst min model ####
mtdna_pi_linear_sstmin <- lmer(logpi ~ range_position + sst.BO_sstmin + 
                                 (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                               REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", 
                               control = lmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
pi_sstmin_boot <- bootMer(mtdna_pi_linear_sstmin, FUN = boot_coef, nsim = 10000)
  boot_coef_sstmin <- as.data.frame(pi_sstmin_boot$t) 

#pull out CIs
lower_ci_sstmin <- as.data.frame(apply(pi_sstmin_boot$t, 2, quantile, 0.025))
upper_ci_sstmin <- as.data.frame(apply(pi_sstmin_boot$t, 2, quantile, 0.975))

cis_sstmin <- cbind(lower_ci_sstmin, upper_ci_sstmin)
  cis_sstmin$fixef <- rownames(cis_sstmin)
  cis_sstmin$model <- "sstmin"
  colnames(cis_sstmin) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### chloroA mean model ####
mtdna_pi_linear_chloromean <- lmer(logpi ~ range_position + logchlomean + I(logchlomean^2) + 
                                     (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                                   REML = FALSE, data = mtdna_small_pi, 
                                   na.action = "na.fail", 
                                   control = lmerControl(optimizer = "bobyqa"))
  
## bootstrap confidence interval for coefficients ##  
pi_chloromean_boot <- bootMer(mtdna_pi_linear_chloromean, FUN = boot_coef, nsim = 10000)
  boot_coef_chloromean <- as.data.frame(pi_chloromean_boot$t) 

#pull out CIs
lower_ci_chloromean <- as.data.frame(apply(pi_chloromean_boot$t, 2, quantile, 0.025))
upper_ci_chloromean <- as.data.frame(apply(pi_chloromean_boot$t, 2, quantile, 0.975))

cis_chloromean <- cbind(lower_ci_chloromean, upper_ci_chloromean)
  cis_chloromean$fixef <- rownames(cis_chloromean)
  cis_chloromean$model <- "chloromean"
  colnames(cis_chloromean) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### chloroA range model ####
mtdna_pi_linear_chlororange <- lmer(logpi ~ range_position + logchlorange + I(logchlorange^2) + 
                                      (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                                    REML = FALSE, data = mtdna_small_pi, 
                                    na.action = "na.fail", 
                                    control = lmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
pi_chlororange_boot <- bootMer(mtdna_pi_linear_chlororange, FUN = boot_coef, nsim = 10000)
  boot_coef_chlororange <- as.data.frame(pi_chlororange_boot$t) 

#pull out CIs
lower_ci_chlororange <- as.data.frame(apply(pi_chlororange_boot$t, 2, quantile, 0.025))
upper_ci_chlororange <- as.data.frame(apply(pi_chlororange_boot$t, 2, quantile, 0.975))

cis_chlororange <- cbind(lower_ci_chlororange, upper_ci_chlororange)
  cis_chlororange$fixef <- rownames(cis_chlororange)
  cis_chlororange$model <- "chlororange"
  colnames(cis_chlororange) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### chloroA max model ####
mtdna_pi_linear_chloromax <- lmer(logpi ~ range_position + logchlomax + I(logchlomax^2) + 
                                    (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                                  REML = FALSE, data = mtdna_small_pi, 
                                  na.action = "na.fail", 
                                  control = lmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
pi_chloromax_boot <- bootMer(mtdna_pi_linear_chloromax, FUN = boot_coef, nsim = 10000)
  boot_coef_chloromax <- as.data.frame(pi_chloromax_boot$t) 

#pull out CIs
lower_ci_chloromax <- as.data.frame(apply(pi_chloromax_boot$t, 2, quantile, 0.025))
upper_ci_chloromax <- as.data.frame(apply(pi_chloromax_boot$t, 2, quantile, 0.975))

cis_chloromax <- cbind(lower_ci_chloromax, upper_ci_chloromax)
  cis_chloromax$fixef <- rownames(cis_chloromax)
  cis_chloromax$model <- "chloromax"
  colnames(cis_chloromax) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### chloroA min model ####
mtdna_pi_linear_chloromin <- lmer(logpi ~ range_position + logchlomin + I(logchlomin^2) + 
                                    (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                                  REML = FALSE, data = mtdna_small_pi, 
                                  na.action = "na.fail", 
                                  control = lmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
pi_chloromin_boot <- bootMer(mtdna_pi_linear_chloromin, FUN = boot_coef, nsim = 10000)
  boot_coef_chloromin <- as.data.frame(pi_chloromin_boot$t) 

#pull out CIs
lower_ci_chloromin <- as.data.frame(apply(pi_chloromin_boot$t, 2, quantile, 0.025))
upper_ci_chloromin <- as.data.frame(apply(pi_chloromin_boot$t, 2, quantile, 0.975))

cis_chloromin <- cbind(lower_ci_chloromin, upper_ci_chloromin)
  cis_chloromin$fixef <- rownames(cis_chloromin)
  cis_chloromin$model <- "chloromax"
  colnames(cis_chloromin) <- c("2.5ci", "97.5ci", "fixef", "model")
  
###############################################################################################################################

######## Merge into one dataframe and write out ########

#rbind cis together
cis_pi <- as.data.frame(rbind(cis_lat, cis_abslat, cis_lon, 
                              cis_sstmean, cis_sstrange, cis_sstmax, cis_sstmin, 
                              cis_chloromean, cis_chlororange, 
                              cis_chloromax, cis_chloromin))

write.csv(cis_pi, "cis_pi.csv")
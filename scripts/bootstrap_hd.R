################################################## Bootstrap Hd models #########################################################

#bootstrapping (1000 iterations/model) to get 95% CIs for coefficients
#do relationships change (e.g. - to +?)
#how sensitive are relationships/models to data input?

################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(glmmTMB) #1.1.7
library(lme4) #v.1.1-31
library(splines) #v.4.2.2

#read in data
mtdna <- read.csv("data/mtdna_assembled.csv", stringsAsFactors = FALSE)
mtdna_env <- read.csv("data/mtdna_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("data/spp_combined_info.csv", stringsAsFactors = FALSE)

#merge dataframes
mtdna <- merge(mtdna, mtdna_env[, c('X', 'sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin', 
                                    'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                    'chloroA.BO_chlomax', 'chloroA.BO_chlomin')], all.x = TRUE)
mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 
                                  'Half_RangeSize', 'Centroid')], all.x = TRUE)

#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)

#################################################################################################################################

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

#scale range position
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

################################################################################################################################
  
######## Null model ########
  
null_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + 
                           (1|Family/Genus) + (1|Source) + (1|MarkerName),
                         data = mtdna_small_hd, family = ordbeta,
                         na.action = "na.fail")  
  
##bootstrap confidence interval for coefficients ##  
hd_null_boot <- bootMer(null_model_hd, FUN = function(x) fixef(x)$cond, nsim = 1000)
  boot_coef_null <- as.data.frame(hd_null_boot$t) 

#pull out CIs
lower_ci_null <- as.data.frame(apply(hd_null_boot$t, 2, quantile, 0.025))
upper_ci_null <- as.data.frame(apply(hd_null_boot$t, 2, quantile, 0.975))

cis_null <- cbind(lower_ci_null, upper_ci_null)
  cis_null$fixef <- rownames(cis_null)
  cis_null$model <- "lat"
  colnames(cis_null) <- c("2.5ci", "97.5ci", "fixef", "model")

#############################################################################################################################

######### Latitude & longitude models ########

#### lat model ####
lat_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                          I(lat_scale^2) + (1|Family/Genus) + 
                          (1|Source) + (1|MarkerName),
                        data = mtdna_small_hd, family = ordbeta,
                        na.action = "na.fail")  

## bootstrap confidence interval for coefficients ##  
hd_lat_boot <- bootMer(lat_model_hd, FUN = function(x) fixef(x)$cond, nsim = 1000)
  boot_coef_lat <- as.data.frame(hd_lat_boot$t) 

#pull out CIs
lower_ci_lat <- as.data.frame(apply(hd_lat_boot$t, 2, quantile, 0.025))
upper_ci_lat <- as.data.frame(apply(hd_lat_boot$t, 2, quantile, 0.975))

cis_lat <- cbind(lower_ci_lat, upper_ci_lat)
  cis_lat$fixef <- rownames(cis_lat)
  cis_lat$model <- "lat"
  colnames(cis_lat) <- c("2.5ci", "97.5ci", "fixef", "model")

#### abslat model ####
abslat_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale +
                             (1|Family/Genus) + (1|Source) + (1|MarkerName),
                           data = mtdna_small_hd, family = ordbeta,
                           na.action = "na.fail")  

## bootstrap confidence interval for coefficients ##  
hd_abslat_boot <- bootMer(abslat_model_hd, FUN = function(x) fixef(x)$cond, nsim = 1000)
  boot_coef_abslat <- as.data.frame(hd_abslat_boot$t) 

#pull out CIs
lower_ci_abslat <- as.data.frame(apply(hd_abslat_boot$t, 2, quantile, 0.025))
upper_ci_abslat <- as.data.frame(apply(hd_abslat_boot$t, 2, quantile, 0.975))

cis_abslat <- cbind(lower_ci_abslat, upper_ci_abslat)
  cis_abslat$fixef <- rownames(cis_abslat)
  cis_abslat$model <- "abslat"
  colnames(cis_abslat) <- c("2.5ci", "97.5ci", "fixef", "model")

#### lon model ####
lon_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + bs(lon_scale) +
                          (1|Family/Genus) + (1|Source) + (1|MarkerName),
                        data = mtdna_small_hd, family = ordbeta,
                        na.action = "na.fail")  

## bootstrap confidence interval for coefficients ##  
hd_lon_boot <- bootMer(lon_model_hd, FUN = function(x) fixef(x)$cond, nsim = 1000)
boot_coef_lon <- as.data.frame(hd_lon_boot$t) 

#pull out CIs
lower_ci_lon <- as.data.frame(apply(hd_lon_boot$t, 2, quantile, 0.025))
upper_ci_lon <- as.data.frame(apply(hd_lon_boot$t, 2, quantile, 0.975))

cis_lon <- cbind(lower_ci_lon, upper_ci_lon)
  cis_lon$fixef <- rownames(cis_lon)
  cis_lon$model <- "lon"
  colnames(cis_lon) <- c("2.5ci", "97.5ci", "fixef", "model")
  
################################################################################################################################
  
######## Environmental models ########
  
#### sst mean model ####
sstmean_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                              (1|Family/Genus) + (1|Source) + (1|MarkerName),
                            data = mtdna_small_hd, family = ordbeta,
                            na.action = "na.fail")  
  
## bootstrap confidence interval for coefficients ##  
hd_sstmean_boot <- bootMer(sstmean_model_hd, FUN = function(x) fixef(x)$cond, nsim = 1000)
boot_coef_sstmean <- as.data.frame(hd_sstmean_boot$t) 
  
#pull out CIs
lower_ci_sstmean <- as.data.frame(apply(hd_sstmean_boot$t, 2, quantile, 0.025))
upper_ci_sstmean <- as.data.frame(apply(hd_sstmean_boot$t, 2, quantile, 0.975))
  
cis_sstmean <- cbind(lower_ci_sstmean, upper_ci_sstmean)
  cis_sstmean$fixef <- rownames(cis_sstmean)
  cis_sstmean$model <- "sstmean"
  colnames(cis_sstmean) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### chloroA mean model ####
chlomean_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean + 
                               I(logchlomean^2) + (1|Family/Genus) + 
                               (1|Source) + (1|MarkerName),
                             data = mtdna_small_hd, family = ordbeta,
                             na.action = "na.fail")  
  
## bootstrap confidence interval for coefficients ##  
hd_chloromean_boot <- bootMer(chlomean_model_hd, FUN = function(x) fixef(x)$cond, nsim = 1000)
boot_coef_chloromean <- as.data.frame(hd_chloromean_boot$t) 
  
#pull out CIs
lower_ci_chloromean <- as.data.frame(apply(hd_chloromean_boot$t, 2, quantile, 0.025))
upper_ci_chloromean <- as.data.frame(apply(hd_chloromean_boot$t, 2, quantile, 0.975))
  
cis_chloromean <- cbind(lower_ci_chloromean, upper_ci_chloromean)
  cis_chloromean$fixef <- rownames(cis_chloromean)
  cis_chloromean$model <- "chloromean"
  colnames(cis_chloromean) <- c("2.5ci", "97.5ci", "fixef", "model")
  
###############################################################################################################################
  
######## Merge into one dataframe and write out ########
  
#rbind cis together
cis_hd <- as.data.frame(rbind(cis_null, cis_lat, cis_abslat, cis_lon, 
                              cis_sstmean, cis_chloromean)) 
  
write.csv(cis_hd, "output/cis_hd.csv")
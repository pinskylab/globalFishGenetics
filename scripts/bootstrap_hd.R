################################################## Bootstrap Hd models #########################################################

#bootstrapping (10,000 iterations/model) to get 95% CIs for coefficients
#do relationships change (e.g. - to +?)
#how sensitive are relationships/models to data input?

#set up to run on Wahab (Old Dominion University HPC)

################################################################################################################################

#set-up
remove(list = ls())

#load libraries
library(lme4)
library(splines)

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

#################################################################################################################################

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

#### Calculate environmental variables ####
## log transform chlorophyll A ##
#subset to only those with chloroA data
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlorange != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomax != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomin != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

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

################################################################################################################################
  
######## Null model ########
  
binomial_null <- glmer(cbind(success, failure) ~ bp_scale + range_position + 
                         (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                       family = binomial, nAGQ = 0,
                       data = mtdna_small_hd, na.action = "na.fail", 
                       control = glmerControl(optimizer = "bobyqa"))
  
##bootstrap confidence interval for coefficients ##  
hd_null_boot <- bootMer(binomial_null, FUN = boot_coef, nsim = 10000)
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
binomial_lat <- glmer(cbind(success, failure) ~ bp_scale + range_position + lat_scale + 
                        I(lat_scale^2) + (1|Family/Genus) + 
                        (1|Source) + (1|MarkerName), 
                      family = binomial, nAGQ = 0,
                      data = mtdna_small_hd, na.action = "na.fail", 
                      control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
hd_lat_boot <- bootMer(binomial_lat, FUN = boot_coef, nsim = 10000)
  boot_coef_lat <- as.data.frame(hd_lat_boot$t) 

#pull out CIs
lower_ci_lat <- as.data.frame(apply(hd_lat_boot$t, 2, quantile, 0.025))
upper_ci_lat <- as.data.frame(apply(hd_lat_boot$t, 2, quantile, 0.975))

cis_lat <- cbind(lower_ci_lat, upper_ci_lat)
  cis_lat$fixef <- rownames(cis_lat)
  cis_lat$model <- "lat"
  colnames(cis_lat) <- c("2.5ci", "97.5ci", "fixef", "model")

#### abslat model ####
binomial_abslat <- glmer(cbind(success, failure) ~ bp_scale + range_position + abslat_scale + 
                         (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                      family = binomial, nAGQ = 0,
                      data = mtdna_small_hd, na.action = "na.fail", 
                      control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
hd_abslat_boot <- bootMer(binomial_abslat, FUN = boot_coef, nsim = 10000)
  boot_coef_abslat <- as.data.frame(hd_abslat_boot$t) 

#pull out CIs
lower_ci_abslat <- as.data.frame(apply(hd_abslat_boot$t, 2, quantile, 0.025))
upper_ci_abslat <- as.data.frame(apply(hd_abslat_boot$t, 2, quantile, 0.975))

cis_abslat <- cbind(lower_ci_abslat, upper_ci_abslat)
  cis_abslat$fixef <- rownames(cis_abslat)
  cis_abslat$model <- "abslat"
  colnames(cis_abslat) <- c("2.5ci", "97.5ci", "fixef", "model")

#### lon model ####
binomial_lon <- glmer(cbind(success, failure) ~ bp_scale + range_position + bs(lon_scale) + 
                        (1|Family/Genus) + (1|Source) + (1|MarkerName),
                      family = binomial, nAGQ = 0,
                      data = mtdna_small_hd, na.action = "na.fail", 
                      control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
hd_lon_boot <- bootMer(binomial_lon, FUN = boot_coef, nsim = 10000)
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
mtdna_hd_binomial_sstmean <- glmer(cbind(success, failure) ~ bp_scale + range_position + sst.BO_sstmean + 
                                     (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                                   family = binomial, nAGQ = 0,
                                   data = mtdna_small_hd, na.action = "na.fail", 
                                   control = glmerControl(optimizer = "bobyqa"))
  
## bootstrap confidence interval for coefficients ##  
hd_sstmean_boot <- bootMer(mtdna_hd_binomial_sstmean, FUN = boot_coef, nsim = 10000)
boot_coef_sstmean <- as.data.frame(hd_sstmean_boot$t) 
  
#pull out CIs
lower_ci_sstmean <- as.data.frame(apply(hd_sstmean_boot$t, 2, quantile, 0.025))
upper_ci_sstmean <- as.data.frame(apply(hd_sstmean_boot$t, 2, quantile, 0.975))
  
cis_sstmean <- cbind(lower_ci_sstmean, upper_ci_sstmean)
  cis_sstmean$fixef <- rownames(cis_sstmean)
  cis_sstmean$model <- "sstmean"
  colnames(cis_sstmean) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### sst range model ####
mtdna_hd_binomial_sstrange <- glmer(cbind(success, failure) ~ bp_scale + range_position + sst.BO_sstrange + 
                                      (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                                    family = binomial, nAGQ = 0,
                                    data = mtdna_small_hd, na.action = "na.fail", 
                                    control = glmerControl(optimizer = "bobyqa"))
  
## bootstrap confidence interval for coefficients ##  
hd_sstrange_boot <- bootMer(mtdna_hd_binomial_sstmrange, FUN = boot_coef, nsim = 10000)
boot_coef_sstrange <- as.data.frame(hd_sstrange_boot$t) 
  
#pull out CIs
lower_ci_sstrange <- as.data.frame(apply(hd_sstrange_boot$t, 2, quantile, 0.025))
upper_ci_sstrange <- as.data.frame(apply(hd_sstrange_boot$t, 2, quantile, 0.975))
  
cis_sstrange <- cbind(lower_ci_sstrange, upper_ci_sstrange)
  cis_sstrange$fixef <- rownames(cis_sstrange)
  cis_sstrange$model <- "sstrange"
  colnames(cis_sstrange) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### sst max model ####
mtdna_hd_binomial_sstmax <- glmer(cbind(success, failure) ~ bp_scale + range_position + sst.BO_sstmax + 
                                    (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                                  family = binomial, nAGQ = 0,
                                  data = mtdna_small_hd, na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa"))
  
## bootstrap confidence interval for coefficients ##  
hd_sstmax_boot <- bootMer(mtdna_hd_binomial_sstmax, FUN = boot_coef, nsim = 10000)
boot_coef_sstmax <- as.data.frame(hd_sstmax_boot$t) 
  
#pull out CIs
lower_ci_sstmax <- as.data.frame(apply(hd_sstmax_boot$t, 2, quantile, 0.025))
upper_ci_sstmax <- as.data.frame(apply(hd_sstmax_boot$t, 2, quantile, 0.975))
  
cis_sstmax <- cbind(lower_ci_sstmax, upper_ci_sstmax)
  cis_sstmax$fixef <- rownames(cis_sstmax)
  cis_sstmax$model <- "sstmax"
  colnames(cis_sstmax) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### sst min model ####
mtdna_hd_binomial_sstmin <- glmer(cbind(success, failure) ~ bp_scale + range_position + sst.BO_sstmin + 
                                    (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                                  family = binomial, nAGQ = 0,
                                  data = mtdna_small_hd, na.action = "na.fail", 
                                  control = glmerControl(optimizer = "bobyqa"))
  
## bootstrap confidence interval for coefficients ##  
hd_sstmin_boot <- bootMer(mtdna_hd_binomial_sstmin, FUN = boot_coef, nsim = 10000)
boot_coef_sstmin <- as.data.frame(hd_sstmin_boot$t) 
  
#pull out CIs
lower_ci_sstmin <- as.data.frame(apply(hd_sstmin_boot$t, 2, quantile, 0.025))
upper_ci_sstmin <- as.data.frame(apply(hd_sstmin_boot$t, 2, quantile, 0.975))
  
cis_sstmin <- cbind(lower_ci_sstmin, upper_ci_sstmin)
  cis_sstmin$fixef <- rownames(cis_sstmin)
  cis_sstmin$model <- "sstmin"
  colnames(cis_sstmin) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### chloroA mean model ####
mtdna_hd_binomial_chloromean <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomean + 
                                        I(logchlomean^2) + (1|Family/Genus) + 
                                        (1|Source) + (1|MarkerName), 
                                      family = binomial, nAGQ = 0,
                                      data = mtdna_small_hd, na.action = "na.fail", 
                                      control = glmerControl(optimizer = "bobyqa"))
  
## bootstrap confidence interval for coefficients ##  
hd_chloromean_boot <- bootMer(mtdna_hd_binomial_chloromean, FUN = boot_coef, nsim = 10000)
boot_coef_chloromean <- as.data.frame(hd_chloromean_boot$t) 
  
#pull out CIs
lower_ci_chloromean <- as.data.frame(apply(hd_chloromean_boot$t, 2, quantile, 0.025))
upper_ci_chloromean <- as.data.frame(apply(hd_chloromean_boot$t, 2, quantile, 0.975))
  
cis_chloromean <- cbind(lower_ci_chloromean, upper_ci_chloromean)
  cis_chloromean$fixef <- rownames(cis_chloromean)
  cis_chloromean$model <- "chloromean"
  colnames(cis_chloromean) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### chloroA range model ####
mtdna_hd_binomial_chlororange <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlorange + 
                                         I(logchlorange^2) + (1|Family/Genus) + 
                                         (1|Source) + (1|MarkerName), 
                                       family = binomial, nAGQ = 0,
                                       data = mtdna_small_hd, na.action = "na.fail", 
                                       control = glmerControl(optimizer = "bobyqa"))
  
## bootstrap confidence interval for coefficients ##  
hd_chlororange_boot <- bootMer(mtdna_hd_binomial_chlororange, FUN = boot_coef, nsim = 10000)
boot_coef_chlororange <- as.data.frame(hd_chlororange_boot$t) 
  
#pull out CIs
lower_ci_chlororange <- as.data.frame(apply(hd_chlororange_boot$t, 2, quantile, 0.025))
upper_ci_chlororange <- as.data.frame(apply(hd_chlororange_boot$t, 2, quantile, 0.975))
  
cis_chlororange <- cbind(lower_ci_chlororange, upper_ci_chlororange)
  cis_chlororange$fixef <- rownames(cis_chlororange)
  cis_chlororange$model <- "chlororange"
  colnames(cis_chlororange) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### chloroA max model ####
mtdna_hd_binomial_chloromax <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomax + 
                                       I(logchlomax^2) + (1|Family/Genus) + 
                                       (1|Source) + (1|MarkerName), 
                                     family = binomial, nAGQ = 0,
                                     data = mtdna_small_hd, na.action = "na.fail", 
                                     control = glmerControl(optimizer = "bobyqa"))
  
## bootstrap confidence interval for coefficients ##  
hd_chloromax_boot <- bootMer(mtdna_hd_binomial_chloromax, FUN = boot_coef, nsim = 10000)
boot_coef_chloromax <- as.data.frame(hd_chloromax_boot$t) 
  
#pull out CIs
lower_ci_chloromax <- as.data.frame(apply(hd_chloromax_boot$t, 2, quantile, 0.025))
upper_ci_chloromax <- as.data.frame(apply(hd_chloromax_boot$t, 2, quantile, 0.975))
  
cis_chloromax <- cbind(lower_ci_chloromax, upper_ci_chloromax)
  cis_chloromax$fixef <- rownames(cis_chloromax)
  cis_chloromax$model <- "chloromax"
  colnames(cis_chloromax) <- c("2.5ci", "97.5ci", "fixef", "model")
  
#### chloroA min model ####
mtdna_hd_binomial_chloromin <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomin + 
                                       I(logchlomin^2) + (1|Family/Genus) + 
                                       (1|Source) + (1|MarkerName), 
                                     family = binomial, nAGQ = 0,
                                     data = mtdna_small_hd, na.action = "na.fail", 
                                     control = glmerControl(optimizer = "bobyqa"))
  
## bootstrap confidence interval for coefficients ##  
hd_chloromin_boot <- bootMer(mtdna_hd_binomial_chloromin, FUN = boot_coef, nsim = 10000)
boot_coef_chloromin <- as.data.frame(hd_chloromin_boot$t) 
  
#pull out CIs
lower_ci_chloromin <- as.data.frame(apply(hd_chloromin_boot$t, 2, quantile, 0.025))
upper_ci_chloromin <- as.data.frame(apply(hd_chloromin_boot$t, 2, quantile, 0.975))
  
cis_chloromin <- cbind(lower_ci_chloromin, upper_ci_chloromin)
  cis_chloromin$fixef <- rownames(cis_chloromin)
  cis_chloromin$model <- "chloromax"
  colnames(cis_chloromin) <- c("2.5ci", "97.5ci", "fixef", "model")
  
###############################################################################################################################
  
######## Merge into one dataframe and write out ########
  
#rbind cis together
cis_hd <- as.data.frame(rbind(cis_null, cis_lat, cis_abslat, cis_lon, 
                              cis_sstmean, cis_sstrange, cis_sstmax, cis_sstmin, 
                              cis_chloromean, cis_chlororange, cis_chloromax, cis_chloromin))
  
write.csv(cis_hd, "cis_hd.csv")
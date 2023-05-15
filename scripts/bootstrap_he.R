################################################## Bootstrap He models #########################################################

#bootstrapping (10,000 iterations/model) to get 95% CIs for coefficients
#do relationships change (e.g. - to +?)
#how sensitive are relationships/models to data input?

#set up to run on Wahab (Old Dominion University HPC)
#Here, reading out each csv (not taking 95% CIs)

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
msat <- read.csv("msatloci_assembled.csv", stringsAsFactors = FALSE)
msat_env <- read.csv("msat_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("spp_combined_info.csv", stringsAsFactors = FALSE)

#merge dataframes
msat <- merge(msat, msat_env[, c('X', 'sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin', 
                                    'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                    'chloroA.BO_chlomax', 'chloroA.BO_chlomin')], all.x = TRUE)
msat <- merge(msat, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 
                                  'Half_RangeSize', 'Centroid')], all.x = TRUE)

#################################################################################################################################

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

##################################################################################################################################

######## Null model ########

binomial_msat_null <- glmer(cbind(success, failure) ~  range_position + CrossSpp + 
                              (1|Family/Genus) + (1|Source) + (1|ID), 
                            family = binomial, data = msat, nAGQ = 0,
                            na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

##bootstrap confidence interval for coefficients ##  
he_null_boot <- bootMer(binomial_msat_null, FUN = boot_coef, nsim = 10000)
  boot_coef_null <- as.data.frame(he_null_boot$t) 

write.csv(boot_coef_null, "/scratch/r3clark/div_bootstrap/boot_he_null.csv")

#####################################################################################################################################

######### Latitude & longitude models ########

#### lat model ####
binomial_msat_lat <- glmer(cbind(success, failure) ~ range_position + CrossSpp + lat_scale + 
                             I(lat_scale^2) + (1|Family/Genus) + (1|Source) + (1|ID), 
                           family = binomial, nAGQ = 0, 
                           data = msat, na.action = "na.fail", 
                           control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
he_lat_boot <- bootMer(binomial_msat_lat, FUN = boot_coef, nsim = 10000)
  boot_coef_lat <- as.data.frame(he_lat_boot$t) 

write.csv(boot_coef_lat, "/scratch/r3clark/div_bootstrap/boot_he_lat.csv")

#### abslat model ####
binomial_msat_abslat <- glmer(cbind(success, failure) ~ range_position + CrossSpp + abslat_scale + 
                            (1|Family/Genus) + (1|Source) + (1|ID), 
                           family = binomial, nAGQ = 0,
                           data = msat, na.action = "na.fail", 
                           control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
he_abslat_boot <- bootMer(binomial_msat_abslat, FUN = boot_coef, nsim = 10000)
  boot_coef_abslat <- as.data.frame(he_abslat_boot$t) 

write.csv(boot_coef_abslat, "/scratch/r3clark/div_bootstrap/boot_he_abslat.csv")
  
#### lon model ####
binomial_msat_lon <- glmer(cbind(success, failure) ~ range_position + CrossSpp + bs(lon_scale) + 
                             (1|Family/Genus) + (1|Source) + (1|ID),
                           family = binomial, nAGQ = 0,
                           data = msat, na.action = "na.fail", 
                           control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
he_lon_boot <- bootMer(binomial_msat_lon, FUN = boot_coef, nsim = 10000)
  boot_coef_lon <- as.data.frame(he_lon_boot$t) 

write.csv(boot_coef_lon, "/scratch/r3clark/div_bootstrap/boot_he_lon.csv")

#####################################################################################################################################

######## Environmental models ########

#### sst mean model ####
msat_he_binomial_sstmean <- glmer(cbind(success, failure) ~  range_position + CrossSpp + sst.BO_sstmean +
                                    (1|Family/Genus) + (1|Source) + (1|ID), 
                                  family = binomial, data = msat, 
                                  na.action = "na.fail", nAGQ = 0,
                                  control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
he_sstmean_boot <- bootMer(msat_he_binomial_sstmean, FUN = boot_coef, nsim = 10000)
boot_coef_sstmean <- as.data.frame(he_sstmean_boot$t) 

write.csv(boot_coef_sstmean, "/scratch/r3clark/div_bootstrap/boot_he_sstmean.csv")

#### sst range model ####
msat_he_binomial_sstrange <- glmer(cbind(success, failure) ~  range_position + CrossSpp + sst.BO_sstrange + 
                                     (1|Family/Genus) + (1|Source) + (1|ID), 
                                   family = binomial, data = msat, 
                                   na.action = "na.fail", nAGQ = 0,
                                   control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
he_sstrange_boot <- bootMer(msat_he_binomial_sstrange, FUN = boot_coef, nsim = 10000)
boot_coef_sstrange <- as.data.frame(he_sstrange_boot$t) 

write.csv(boot_coef_sstrange, "/scratch/r3clark/div_bootstrap/boot_he_sstrange.csv")

#### sst max model ####
msat_he_binomial_sstmax <- glmer(cbind(success, failure) ~  range_position + CrossSpp + sst.BO_sstmax + 
                                   (1|Family/Genus) + (1|Source) + (1|ID), 
                                 family = binomial, data = msat, 
                                 na.action = "na.fail", nAGQ = 0,
                                 control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
he_sstmax_boot <- bootMer(msat_he_binomial_sstmax, FUN = boot_coef, nsim = 10000)
boot_coef_sstmax <- as.data.frame(he_sstmax_boot$t) 

write.csv(boot_coef_sstmax, "/scratch/r3clark/div_bootstrap/boot_he_sstmax.csv")

#### sst min model ####
msat_he_binomial_sstmin <- glmer(cbind(success, failure) ~  range_position + CrossSpp + sst.BO_sstmin + 
                                   (1|Family/Genus) + (1|Source) + (1|ID), 
                                 family = binomial, data = msat, 
                                 na.action = "na.fail", nAGQ = 0,
                                 control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
he_sstmin_boot <- bootMer(msat_he_binomial_sstmin, FUN = boot_coef, nsim = 10000)
boot_coef_sstmin <- as.data.frame(he_sstmin_boot$t) 

write.csv(boot_coef_sstmin, "/scratch/r3clark/div_bootstrap/boot_he_sstmin.csv")

#### chloroA mean model ####
msat_he_binomial_chloromean <- glmer(cbind(success, failure) ~ range_position + CrossSpp + 
                                       logchlomean + I(logchlomean^2) + 
                                       (1|Family/Genus) + (1|Source) + (1|ID), 
                                     family = binomial, nAGQ = 0,
                                     data = msat, na.action = "na.fail", 
                                     control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
he_chloromean_boot <- bootMer(msat_he_binomial_chloromean, FUN = boot_coef, nsim = 10000)
boot_coef_chloromean <- as.data.frame(he_chloromean_boot$t) 

write.csv(boot_coef_chloromean, "/scratch/r3clark/div_bootstrap/boot_he_chloromean.csv")

#### chloroA range model ####
msat_he_binomial_chlororange <- glmer(cbind(success, failure) ~ range_position + CrossSpp + 
                                        logchlorange + I(logchlorange^2) + 
                                        (1|Family/Genus) + (1|Source) + (1|ID), 
                                      family = binomial, nAGQ = 0,
                                      data = msat, na.action = "na.fail", 
                                      control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
he_chlororange_boot <- bootMer(msat_he_binomial_chlororange, FUN = boot_coef, nsim = 10000)
boot_coef_chlororange <- as.data.frame(he_chlororange_boot$t) 

write.csv(boot_coef_chlororange, "/scratch/r3clark/div_bootstrap/boot_he_chlororange.csv")

#### chloroA max model ####
msat_he_binomial_chloromax <- glmer(cbind(success, failure) ~ range_position + CrossSpp +
                                      logchlomax + I(logchlomax^2) + 
                                      (1|Family/Genus) + (1|Source) + (1|ID), 
                                    family = binomial, nAGQ = 0,
                                    data = msat, na.action = "na.fail", 
                                    control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
he_chloromax_boot <- bootMer(msat_he_binomial_chloromax, FUN = boot_coef, nsim = 10000)
boot_coef_chloromax <- as.data.frame(he_chloromax_boot$t) 

write.csv(boot_coef_chloromax, "/scratch/r3clark/div_bootstrap/boot_he_chloromax.csv")

#### chloroA min model ####
msat_he_binomial_chloromin <- glmer(cbind(success, failure) ~ range_position + CrossSpp + 
                                      logchlomin + I(logchlomin^2) + 
                                      (1|Family/Genus) + (1|Source) + (1|ID), 
                                    family = binomial, nAGQ = 0,
                                    data = msat, na.action = "na.fail", 
                                    control = glmerControl(optimizer = "bobyqa"))

## bootstrap confidence interval for coefficients ##  
he_chloromin_boot <- bootMer(msat_he_binomial_chloromin, FUN = boot_coef, nsim = 10000)
boot_coef_chloromin <- as.data.frame(he_chloromin_boot$t) 

write.csv(boot_coef_chloromin, "/scratch/r3clark/div_bootstrap/boot_he_chloromin.csv")
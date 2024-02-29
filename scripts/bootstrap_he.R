################################################## Bootstrap He models #########################################################

#bootstrapping (1000 iterations/model) to get 95% CIs for coefficients
#do relationships change (e.g. - to +?)
#how sensitive are relationships/models to data input?

#here, reading out each csv (not taking 95% CIs)

################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(glmmTMB) #1.1.7
library(lme4) #v.1.1-31
library(splines) #v.4.2.2

#read in data
msat <- read.csv("data/msatloci_assembled.csv", stringsAsFactors = FALSE)
msat_env <- read.csv("data/msat_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("data/spp_combined_info.csv", stringsAsFactors = FALSE)

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

#scale range_position
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

##################################################################################################################################

######### Latitude & longitude models ########

#### lat model ####
lat_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + lat_scale + 
                          I(lat_scale^2) + (1|Family/Genus) + (1|Source),
                        data = msat, family = ordbeta, na.action = "na.fail")

##bootstrap confidence interval for coefficients ##  
he_lat_boot <- bootMer(lat_model_he, FUN = function(x) fixef(x)$cond, nsim = 1000)
boot_coef_lat <- as.data.frame(he_lat_boot$t) 

write.csv(boot_coef_lat, "output/boot_he_lat.csv")

#### abslat model ####
abslat_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + abslat_scale + 
                             (1|Family/Genus) + (1|Source),
                           data = msat, family = ordbeta, na.action = "na.fail")

##bootstrap confidence interval for coefficients ##  
he_abslat_boot <- bootMer(abslat_model_he, FUN = function(x) fixef(x)$cond, nsim = 1000)
boot_coef_abslat <- as.data.frame(he_abslat_boot$t) 

write.csv(boot_coef_abslat, "output/boot_he_abslat.csv")
  
#### lon model ####
lon_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + bs(lon_scale) + 
                          (1|Family/Genus) + (1|Source),
                        data = msat, family = ordbeta, na.action = "na.fail")

##bootstrap confidence interval for coefficients ##  
he_lon_boot <- bootMer(lon_model_he, FUN = function(x) fixef(x)$cond, nsim = 1000)
boot_coef_lon <- as.data.frame(he_lon_boot$t) 

write.csv(boot_coef_lon, "output/boot_he_lon.csv")

#####################################################################################################################################

######## Environmental models ########

#### sst mean model ####
sstmean_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + sstmean_scale + 
                              (1|Family/Genus) + (1|Source),
                            data = msat, family = ordbeta, na.action = "na.fail")

##bootstrap confidence interval for coefficients ##  
he_sstmean_boot <- bootMer(sstmean_model_he, FUN = function(x) fixef(x)$cond, nsim = 1000)
boot_coef_sstmean <- as.data.frame(he_sstmean_boot$t) 

write.csv(boot_coef_sstmean, "output/boot_he_sstmean.csv")

#### chloroA mean model ####
chlomean_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + logchlomean + 
                               I(logchlomean^2) + (1|Family/Genus) + (1|Source),
                             data = msat, family = ordbeta, na.action = "na.fail")

##bootstrap confidence interval for coefficients ##  
he_chlomean_boot <- bootMer(chlomean_model_he, FUN = function(x) fixef(x)$cond, nsim = 1000)
boot_coef_chlomean <- as.data.frame(he_chlomean_boot$t) 

write.csv(boot_coef_chlomean, "output/boot_he_chlomean.csv")
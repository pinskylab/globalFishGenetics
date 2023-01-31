######################## Script for Building Regression Figures for msat ##############################

#Uses final msat he binomial regression models
#Predicts he at certain values of predictor variable of interest
#Calculates mean he of raw data (binned every X units)
#Plots two together

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(tidyverse)
library(plotrix)
library(lme4)
library(data.table)
library(DescTools)
library(sjPlot)

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

####################################################################################################################

######## Clean up dataframe ########

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

#add random effect for every unit
msat$ID <- (1:28142)

#### Calculate success and failure ####
msat$success <- round(msat$He*msat$n) #number of heterozygotes
msat$failure<- round((1 - msat$He)*msat$n) #number of homozygotes

#### Calculate latitude, longitude variables ####
#calculate abslat
msat$abslat <- abs(msat$lat)

#scale geographic variables
msat$lat_scale <- as.numeric(scale(msat$lat))
msat$abslat_scale <- as.numeric(scale(msat$abslat))

#convert lon to radians
msat$lon_360 <- msat$lon + 180 #convert (-180,180) to (0,360)
  msat$lon_rad <- (2*pi*msat$lon_360)/360

#############################################################################################################

######### Abslat figures #######

abslat_model_he <- glm(cbind(success, failure) ~ CrossSpp + abslat_scale, 
                       family = binomial, data = msat, na.action = "na.fail")

#### Predict ####
#marginal effects
abslat_eff <- plot_model(abslat_model_he, type = "pred",
                         terms = "abslat_scale [all]")

#pull out marginal effects dataframe
abslat_eff_data <- as.data.frame(abslat_eff$data)

#unscale lat
#use same scaled:center & scaled:scale from original data
abslat_scale <- scale(msat$abslat) #bc had to convert to numeric to run model/calculate marginal effects
abslat_eff_data$abslat <- (abslat_eff_data$x * attr(abslat_scale, "scaled:scale")) + 
  attr(abslat_scale, "scaled:center")

#### Calculate means from raw data ####
## Round abslat DOWN to nearest multiple of 10 ##
#create column to fill with the rounded abslat
msat$abslat_round <- NA

#fill round column in
for(i in 1:nrow(msat)) {
  cat(paste(i, " ", sep = ''))
  {msat$abslat_round[i] <- DescTools::RoundTo(msat$abslat[i], 
                                              multiple = 10, FUN = floor)} #rounding DOWN (e.g. abslat 10-20 gets value of 10)
}

msat$abslat_round <- as.factor(msat$abslat_round)

## Calculate means and SE within each 10 degree band ##
msat <- data.table(msat) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
abslat_he_mean <- msat[, mean(He), by = abslat_round]
  colnames(abslat_he_mean) <- c("abslat", "he_mean")

#calculate SE in each 10 degree band
abslat_he_SE <- msat[, std.error(He), by = abslat_round]
  colnames(abslat_he_SE) <- c("abslat", "he_SE")

#count observations in each 10 degree band
abslat_he_count <- msat[, .N, by = abslat_round]
  colnames(abslat_he_count) <- c("abslat", "he_count")

#merge mean and SE dataframes together
abslat_he_binned_means <- list(abslat_he_mean, abslat_he_SE, abslat_he_count) %>% 
  reduce(full_join, by = "abslat")
abslat_he_binned_means$abslat <- as.numeric(as.character(abslat_he_binned_means$abslat))
  abslat_he_binned_means <- abslat_he_binned_means[order(abslat), ]
  abslat_he_binned_means$X <- abslat_he_binned_means$abslat + 5 #for plotting, plot in MIDDLE of 10 degree band

#calculate error bars (standard error)
abslat_he_binned_means$mean_lowerSE <- abslat_he_binned_means$he_mean - 
  abslat_he_binned_means$he_SE
abslat_he_binned_means$mean_upperSE <- abslat_he_binned_means$he_mean + 
  abslat_he_binned_means$he_SE

#### Plot abslat ####
#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_he_abslat_plot_both <- ggplot() + 
  geom_line(data = abslat_eff_data, 
            aes(x = abslat, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = abslat_eff_data, 
              aes(x = abslat, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = abslat_he_binned_means, 
             aes(x = X, y = he_mean, color = "10-degree binned means", size = he_count), shape = "square") + 
  geom_errorbar(data = abslat_he_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("absolute latitude") + ylab("msat He") + labs(color = "Legend") + 
  scale_y_continuous(limits = c(0, 1.0)) + 
  scale_color_manual(values = colors) + 
  scale_size_continuous(breaks = c(500, 1000, 2000, 5000), 
                        range = c(10, 20))
mtdna_he_abslat_plot_annotated_both <- mtdna_he_abslat_plot_both + theme_bw() + #coord_flip() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_he_abslat_plot_annotated_both

#############################################################################################################

######### Lat figures #######

lat_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + lat_scale + I(lat_scale^2), 
                    family = binomial, data = mtdna_small_hd, na.action = "na.fail")

#### Predict ####
#marginal effects
lat_eff <- plot_model(lat_model_hd, type = "pred",
                      terms = "lat_scale [all]")

#pull out marginal effects dataframe
lat_eff_data <- as.data.frame(lat_eff$data)

#unscale lat
#use same scaled:center & scaled:scale from original data
lat_scale <- scale(mtdna_small_hd$lat) #bc had to convert to numeric to run model/calculate marginal effects
lat_eff_data$lat <- (lat_eff_data$x * attr(lat_scale, "scaled:scale")) + 
  attr(lat_scale, "scaled:center")

#### Calculate means from raw data ####
## Round lat DOWN to nearest multiple of 10 ##
#create column to fill with the rounded lat
mtdna_small_hd$lat_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_hd)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_hd$lat_round[i] <- DescTools::RoundTo(mtdna_small_hd$lat[i], 
                                                     multiple = 10, FUN = floor)} #rounding DOWN (e.g. lat 10-20 gets value of 10)
}

mtdna_small_hd$lat_round <- as.factor(mtdna_small_hd$lat_round)

## Calculate means and SE within each 10 degree band ##
mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
lat_hd_mean <- mtdna_small_hd[, mean(He), by = lat_round]
colnames(lat_hd_mean) <- c("lat", "hd_mean")

#calculate SE in each 10 degree band
lat_hd_SE <- mtdna_small_hd[, std.error(He), by = lat_round]
colnames(lat_hd_SE) <- c("lat", "hd_SE")

#count observations in each 10 degree band
lat_hd_count <- mtdna_small_hd[, .N, by = lat_round]
colnames(lat_hd_count) <- c("lat", "hd_count")

#merge mean and SE dataframes together
lat_hd_binned_means <- list(lat_hd_mean, lat_hd_SE, lat_hd_count) %>% 
  reduce(full_join, by = "lat")
lat_hd_binned_means$lat <- as.numeric(as.character(lat_hd_binned_means$lat))
lat_hd_binned_means <- lat_hd_binned_means[order(lat), ]
lat_hd_binned_means$X <- lat_hd_binned_means$lat + 5 #for plotting, plot in MIDDLE of 10 degree band

#calculate error bars (standard error)
lat_hd_binned_means$mean_lowerSE <- lat_hd_binned_means$hd_mean - 
  lat_hd_binned_means$hd_SE
lat_hd_binned_means$mean_upperSE <- lat_hd_binned_means$hd_mean + 
  lat_hd_binned_means$hd_SE

#### Plot lat ####
#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_hd_lat_plot_both <- ggplot() + 
  geom_line(data = lat_eff_data, 
            aes(x = lat, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = lat_eff_data, 
              aes(x = lat, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = lat_hd_binned_means, 
             aes(x = X, y = hd_mean, color = "10-degree binned means", size = hd_count), shape = "square") + 
  geom_errorbar(data = lat_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("latitude") + ylab("mtdna Hd") + labs(color = "Legend") + 
  scale_y_continuous(limits = c(0, 1.0)) + 
  scale_color_manual(values = colors) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_hd_lat_plot_annotated_both <- mtdna_hd_lat_plot_both + theme_bw() + #coord_flip() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_lat_plot_annotated_both

#############################################################################################################

######### Lon figures #######

lon_model_he <- glm(cbind(success, failure) ~ CrossSpp + sin(lon_rad) + cos(lon_rad), 
                    family = binomial, data = msat, na.action = "na.fail")

#### Predict ####
#marginal effects
lon_eff <- plot_model(lon_model_he, type = "pred",
                      terms = "lon_rad [all]")

#pull out marginal effects dataframe
lon_eff_data <- as.data.frame(lon_eff$data)
  lon_eff_data$lon <- ((360*lon_eff_data$x)/(2*pi))-180

#### Calculate means from raw data ####
## Round lon DOWN to nearest multiple of 10 ##
#create column to fill with the rounded lon
msat$lon_round <- NA

#fill round column in
for(i in 1:nrow(msat)) {
  cat(paste(i, " ", sep = ''))
  {msat$lon_round[i] <- DescTools::RoundTo(msat$lon[i], 
                                                     multiple = 10, FUN = floor)} #rounding DOWN (e.g. lat 10-20 gets value of 10)
}

  msat$lon_round <- as.factor(msat$lon_round)

## Calculate means and SE within each 10 degree band ##
msat <- data.table(msat) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
lon_he_mean <- msat[, mean(He), by = lon_round]
  colnames(lon_he_mean) <- c("lon", "he_mean")

#calculate SE in each 10 degree band
lon_he_SE <- msat[, std.error(He), by = lon_round]
  colnames(lon_he_SE) <- c("lon", "he_SE")

#count observations in each 10 degree band
lon_he_count <- msat[, .N, by = lon_round]
  colnames(lon_he_count) <- c("lon", "he_count")

#merge mean and SE dataframes together
lon_he_binned_means <- list(lon_he_mean, lon_he_SE, lon_he_count) %>% 
  reduce(full_join, by = "lon")
lon_he_binned_means$lon <- as.numeric(as.character(lon_he_binned_means$lon))
  lon_he_binned_means <- lon_he_binned_means[order(lon), ]
  lon_he_binned_means$X <- lon_he_binned_means$lon + 5 #for plotting, plot in MIDDLE of 10 degree band

#calculate error bars (standard error)
lon_he_binned_means$mean_lowerSE <- lon_he_binned_means$he_mean - 
  lon_he_binned_means$he_SE
lon_he_binned_means$mean_upperSE <- lon_he_binned_means$he_mean + 
  lon_he_binned_means$he_SE

#### Plot lon ####
#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_he_lon_plot_both <- ggplot() + 
  geom_line(data = lon_eff_data, 
            aes(x = lon, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = lon_eff_data, 
              aes(x = lon, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = lon_he_binned_means, 
             aes(x = X, y = he_mean, color = "10-degree binned means", size = he_count), shape = "square") + 
  geom_errorbar(data = lon_he_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("longitude") + ylab("msat He") + labs(color = "Legend") + 
  scale_y_continuous(limits = c(0, 1.0)) + 
  scale_color_manual(values = colors) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_he_lon_plot_annotated_both <- mtdna_he_lon_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_he_lon_plot_annotated_both

#############################################################################################################

######## Calculate environmental variables ########

#### log transform sst data ####
#subset to only those with sst data
msat <- subset(msat, msat$sst.BO_sstmean != "NA") #shouldn't remove any

msat$logsstmean <- log10(msat$sst.BO_sstmean)
  msat$logsstrange <- log10(msat$sst.BO_sstrange)
  msat$logsstmax <- log10(msat$sst.BO_sstmax)
  msat$logsstmin <- log10(msat$sst.BO_sstmin)

#remove logsst = NA columns
msat <- subset(msat, msat$logsstmean != "NaN")
  msat <- subset(msat, msat$logsstrange != "NaN")
  msat <- subset(msat, msat$logsstmax != "NaN")
  msat <- subset(msat, msat$logsstmin != "NaN")

#### log transform dissox data ####
#subset to only those with dissox data
msat <- subset(msat, msat$BO_dissox != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

msat$logdissox <- log10(msat$BO_dissox)

#remove logdissox = NA columns
msat <- subset(msat, msat$logdissox != "NaN")

#### log transform chlorophyll A ####
#subset to only those with chloroA data
msat <- subset(msat, msat$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat <- subset(msat, msat$chloroA.BO_chlorange != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat <- subset(msat, msat$chloroA.BO_chlomax != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat <- subset(msat, msat$chloroA.BO_chlomin != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

msat$logchlomean <- log10(msat$sst.BO_sstmean)
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

##################################################################################################################

######### SST mean figures #######

SSTmean_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + logsstmean, 
                        family = binomial, data = mtdna_small_hd, na.action = "na.fail")

#### Predict ####
#marginal effects
SSTmean_eff <- plot_model(SSTmean_model_hd, type = "pred",
                          terms = "logsstmean [all]")

#pull out marginal effects dataframe
SSTmean_eff_data <- as.data.frame(SSTmean_eff$data)
SSTmean_eff_data$sstmean <- 10^(SSTmean_eff_data$x)

#### Calculate means from raw data ####
## Round lon DOWN to nearest multiple of 5 ##
#create column to fill with the rounded SSTmean
mtdna_small_hd$SSTmean_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_hd)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_hd$SSTmean_round[i] <- DescTools::RoundTo(mtdna_small_hd$sst.BO_sstmean[i], 
                                                         multiple = 5, FUN = floor)} #rounding DOWN
}

mtdna_small_hd$SSTmean_round <- as.factor(mtdna_small_hd$SSTmean_round)

## Calculate means and SE within each 5 degree band ##
mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate mean in each 5 degree band
SSTmean_hd_mean <- mtdna_small_hd[, mean(He), by = SSTmean_round]
colnames(SSTmean_hd_mean) <- c("SSTmean", "hd_mean")

#calculate SE in each 5 degree band
SSTmean_hd_SE <- mtdna_small_hd[, std.error(He), by = SSTmean_round]
colnames(SSTmean_hd_SE) <- c("SSTmean", "hd_SE")

#count observations in each 5 degree band
SSTmean_hd_count <- mtdna_small_hd[, .N, by = SSTmean_round]
colnames(SSTmean_hd_count) <- c("SSTmean", "hd_count")

#merge mean and SE dataframes together
SSTmean_hd_binned_means <- list(SSTmean_hd_mean, SSTmean_hd_SE, SSTmean_hd_count) %>% 
  reduce(full_join, by = "SSTmean")
SSTmean_hd_binned_means$SSTmean <- as.numeric(as.character(SSTmean_hd_binned_means$SSTmean))
SSTmean_hd_binned_means <- SSTmean_hd_binned_means[order(SSTmean), ]
SSTmean_hd_binned_means$X <- SSTmean_hd_binned_means$SSTmean + 2.5 #for plotting, plot in MIDDLE of 5 degree band

#calculate error bars (standard error)
SSTmean_hd_binned_means$mean_lowerSE <- SSTmean_hd_binned_means$hd_mean - 
  SSTmean_hd_binned_means$hd_SE
SSTmean_hd_binned_means$mean_upperSE <- SSTmean_hd_binned_means$hd_mean + 
  SSTmean_hd_binned_means$hd_SE

#### Plot SSTmean ####
#for legend
colors <- c("5-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_hd_SSTmean_plot_both <- ggplot() + 
  geom_line(data = SSTmean_eff_data, 
            aes(x = sstmean, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = SSTmean_eff_data, 
              aes(x = sstmean, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = SSTmean_hd_binned_means, 
             aes(x = X, y = hd_mean, color = "5-degree binned means", size = hd_count), shape = "square") + 
  geom_errorbar(data = SSTmean_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "5-degree binned means"), width = 1.75, size = 3) + 
  xlab("SST mean") + ylab("mtdna Hd") + labs(color = "Legend") + 
  scale_y_continuous(limits = c(0, 1.0)) + 
  scale_color_manual(values = colors) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_hd_SSTmean_plot_annotated_both <- mtdna_hd_SSTmean_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_SSTmean_plot_annotated_both

##################################################################################################################

######### SST min figures #######

SSTmin_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + logsstmin + 
                           (1|Family/Genus/spp) + (1|Source) + (1|MarkerName), 
                         family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                         control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#### Predict ####
#marginal effects
SSTmin_eff <- plot_model(SSTmin_model_hd, type = "pred",
                         terms = "logsstmin [all]")

#pull out marginal effects dataframe
SSTmin_eff_data <- as.data.frame(SSTmin_eff$data)
SSTmin_eff_data$sstmin <- 10^(SSTmin_eff_data$x)

#### Calculate means from raw data ####
## Round lon DOWN to nearest multiple of 5 ##
#create column to fill with the rounded SSTmin
mtdna_small_hd$SSTmin_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_hd)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_hd$SSTmin_round[i] <- DescTools::RoundTo(mtdna_small_hd$sst.BO_sstmin[i], 
                                                        multiple = 5, FUN = floor)} #rounding DOWN
}

mtdna_small_hd$SSTmin_round <- as.factor(mtdna_small_hd$SSTmin_round)

## Calculate means and SE within each 5 degree band ##
mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate mean in each 5 degree band
SSTmin_hd_mean <- mtdna_small_hd[, mean(He), by = SSTmin_round]
colnames(SSTmin_hd_mean) <- c("SSTmin", "hd_mean")

#calculate SE in each 5 degree band
SSTmin_hd_SE <- mtdna_small_hd[, std.error(He), by = SSTmin_round]
colnames(SSTmin_hd_SE) <- c("SSTmin", "hd_SE")

#count observations in each 5 degree band
SSTmin_hd_count <- mtdna_small_hd[, .N, by = SSTmin_round]
colnames(SSTmin_hd_count) <- c("SSTmin", "hd_count")

#merge mean and SE dataframes together
SSTmin_hd_binned_means <- list(SSTmin_hd_mean, SSTmin_hd_SE, SSTmin_hd_count) %>% 
  reduce(full_join, by = "SSTmin")
SSTmin_hd_binned_means$SSTmin <- as.numeric(as.character(SSTmin_hd_binned_means$SSTmin))
SSTmin_hd_binned_means <- SSTmin_hd_binned_means[order(SSTmin), ]
SSTmin_hd_binned_means$X <- SSTmin_hd_binned_means$SSTmin + 2.5 #for plotting, plot in MIDDLE of 5 degree band

#calculate error bars (standard error)
SSTmin_hd_binned_means$mean_lowerSE <- SSTmin_hd_binned_means$hd_mean - 
  SSTmin_hd_binned_means$hd_SE
SSTmin_hd_binned_means$mean_upperSE <- SSTmin_hd_binned_means$hd_mean + 
  SSTmin_hd_binned_means$hd_SE

#### Plot SSTmin ####
#for legend
colors <- c("5-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_hd_SSTmin_plot_both <- ggplot() + 
  geom_line(data = SSTmin_eff_data, 
            aes(x = sstmin, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = SSTmin_eff_data, 
              aes(x = sstmin, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = SSTmin_hd_binned_means, 
             aes(x = X, y = hd_mean, color = "5-degree binned means", size = hd_count), shape = "square") + 
  geom_errorbar(data = SSTmin_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "5-degree binned means"), width = 1.75, size = 3) + 
  xlab("SST min") + ylab("mtdna Hd") + labs(color = "Legend") + 
  scale_y_continuous(limits = c(0, 1.0)) + 
  scale_color_manual(values = colors) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_hd_SSTmin_plot_annotated_both <- mtdna_hd_SSTmin_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_SSTmin_plot_annotated_both

###########################################################################################################

######## Dissox figures ########

dissox_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + logdissox + 
                           (1|Family/Genus/spp) + (1|Source) + (1|MarkerName), 
                         family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                         control = glmerControl(optimizer = "bobyqa")) #had to add bobyqa to converge

#### Predict ####
#marginal effects
dissox_eff <- plot_model(dissox_model_hd, type = "pred", 
                         terms = "logdissox [all]")

#pull out marginal effects dataframe
dissox_eff_data <- as.data.frame(dissox_eff$data)
dissox_eff_data$dissox <- 10^(dissox_eff_data$x)

#### Calculate means from raw data ####
## Round dissox DOWN to nearest multiple of 5 ##
#create column to fill with the rounded dissox
mtdna_small_hd$dissox_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_hd)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_hd$dissox_round[i] <- DescTools::RoundTo(mtdna_small_hd$BO_dissox[i], 
                                                        multiple = 0.5, FUN = floor)} #rounding DOWN
}

mtdna_small_hd$dissox_round <- as.factor(mtdna_small_hd$dissox_round)

## Calculate means and SE within each 0.5 band ##
mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so that can use data.table functions

#calculate mean in each 0.5 band
dissox_hd_mean <- mtdna_small_hd[, mean(He), by = dissox_round]
colnames(dissox_hd_mean) <- c("dissox", "hd_mean")

#calculate SE in each 0.5 band  
dissox_hd_SE <- mtdna_small_hd[, std.error(He), by = dissox_round]
colnames(dissox_hd_SE) <- c("dissox", "hd_SE")

#count observations in each 0.5 band
dissox_hd_count <- mtdna_small_hd[, .N, by = dissox_round]
colnames(dissox_hd_count) <- c("dissox", "hd_count")

#merge mean and SE dataframes together
dissox_hd_binned_means <- list(dissox_hd_mean, dissox_hd_SE, dissox_hd_count) %>% 
  reduce(full_join, by = "dissox")
dissox_hd_binned_means$dissox <- as.numeric(as.character(dissox_hd_binned_means$dissox))
dissox_hd_binned_means <- dissox_hd_binned_means[order(dissox), ]
dissox_hd_binned_means$X <- dissox_hd_binned_means$dissox + 0.25 #for plotting, plot in MIDDLE of 0.5 band

#add error bars (standard error)
dissox_hd_binned_means$mean_lowerSE <- dissox_hd_binned_means$hd_mean - 
  dissox_hd_binned_means$hd_SE
dissox_hd_binned_means$mean_upperSE <- dissox_hd_binned_means$hd_mean + 
  dissox_hd_binned_means$hd_SE

#### Plot dissox ####
#for legend
colors <- c("Binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_hd_dissox_plot_both <- ggplot() + 
  geom_line(data = dissox_eff_data, 
            aes(x = dissox, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = dissox_eff_data, 
              aes(x = dissox, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = dissox_hd_binned_means, 
             aes(x = X, y = hd_mean, color = "Binned means", size = hd_count), shape = "square") + 
  geom_errorbar(data = dissox_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "Binned means"), width = 0.25, size = 3) + 
  xlab("Dissolved Oxygen Mean") + ylab("mtDNA Hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 1.0)) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_hd_dissox_plot_annotated_both <- mtdna_hd_dissox_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_dissox_plot_annotated_both

#############################################################################################################

######## ChloroA mean figures ########

chloroAmean_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlomean + 
                                I(logchlomean^2) + (1|Family/Genus/spp) + (1|Source) + 
                                (1|MarkerName), family = binomial, 
                              data = mtdna_small_hd, na.action = "na.fail", 
                              control = glmerControl(optimizer = "bobyqa"))
#### Predict ####
#marginal effects
chloroAmean_eff <- plot_model(chloroAmean_model_hd, type = "pred", 
                              terms = "logchlomean [all]")

#pull out marginal effects dataframe
chloroAmean_eff_data <- as.data.frame(chloroAmean_eff$data)
chloroAmean_eff_data$chlomean <- 10^(chloroAmean_eff_data$x)

#### Calculate means from raw data ####
## Round chlomean DOWN to nearest multiple of 10 ##
#create column to fill with the rounded chlomean
mtdna_small_hd$chloroAmean_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_hd)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_hd$chloroAmean_round[i] <- DescTools::RoundTo(mtdna_small_hd$chloroA.BO_chlomean[i], 
                                                             multiple = 10, FUN = floor)} #rounding DOWN
}

mtdna_small_hd$chloroAmean_round <- as.factor(mtdna_small_hd$chloroAmean_round)

## Calculate means and SE within each band ##
mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate mean in each band
chloroAmean_hd_mean <- mtdna_small_hd[, mean(He), by = chloroAmean_round]
colnames(chloroAmean_hd_mean) <- c("chloroAmean", "hd_mean")

#calculate SE in each band
chloroAmean_hd_SE <- mtdna_small_hd[, std.error(He), by = chloroAmean_round]
colnames(chloroAmean_hd_SE) <- c("chloroAmean", "hd_SE")

#count observations in each band
chloroAmean_hd_count <- mtdna_small_hd[, .N, by = chloroAmean_round]
colnames(chloroAmean_hd_count) <- c("chloroAmean", "hd_count")

#merge mean and SE dataframes together
chloroAmean_hd_binned_means <- list(chloroAmean_hd_mean, chloroAmean_hd_SE, chloroAmean_hd_count) %>% 
  reduce(full_join, by = "chloroAmean")
chloroAmean_hd_binned_means$chloroAmean <- as.numeric(as.character(chloroAmean_hd_binned_means$chloroAmean))
chloroAmean_hd_binned_means <- chloroAmean_hd_binned_means[order(chloroAmean), ]
chloroAmean_hd_binned_means$X <- chloroAmean_hd_binned_means$chloroAmean + 5 #for plotting, plot in MIDDLE of band

#calculate error bars (standard error)
chloroAmean_hd_binned_means$mean_lowerSE <- chloroAmean_hd_binned_means$hd_mean - 
  chloroAmean_hd_binned_means$hd_SE
chloroAmean_hd_binned_means$mean_upperSE <- chloroAmean_hd_binned_means$hd_mean + 
  chloroAmean_hd_binned_means$hd_SE

#### Plot chloroAmean ####
#for legend
colors <- c("Binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_hd_chloroAmean_plot_both <- ggplot() + 
  geom_line(data = chloroAmean_eff_data, 
            aes(x = chlomean, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = chloroAmean_eff_data, 
              aes(x = chlomean, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = chloroAmean_hd_binned_means, 
             aes(x = X, y = hd_mean, color = "Binned means", size = hd_count), shape = "square") + 
  geom_errorbar(data = chloroAmean_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "Binned means"), width = 1.75, size = 3) + 
  xlab("Chloro A mean") + ylab("mtDNA Hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 1.0)) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_hd_chloroAmean_plot_annotated_both <- mtdna_hd_chloroAmean_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_chloroAmean_plot_annotated_both

#############################################################################################################

######## ChloroA max figures ########

chloroAmax_model_hd <- glm(cbind(success, failure) ~ CrossSpp + logchlomax + 
                                 I(logchlomax^2), family = binomial, 
                               data = msat, na.action = "na.fail"
                           )
#### Predict ####
#marginal effects
chloroAmax_eff <- plot_model(chloroAmax_model_hd, type = "pred", 
                               terms = "logchlomax [all]")

#pull out marginal effects dataframe
chloroAmax_eff_data <- as.data.frame(chloroAmax_eff$data)
  chloroAmax_eff_data$chlomax <- 10^(chloroAmax_eff_data$x)

#### Calculate means from raw data ####
## Round chlomax DOWN to nearest multiple of 10 ##
#create column to fill with the rounded chlomean
msat$chloroAmax_round <- NA

#fill round column in
for(i in 1:nrow(msat)) {
  cat(paste(i, " ", sep = ''))
  {msat$chloroAmax_round[i] <- DescTools::RoundTo(msat$chloroA.BO_chlomax[i], 
                                                              multiple = 10, FUN = floor)} #rounding DOWN
}

msat$chloroAmax_round <- as.factor(msat$chloroAmax_round)

## Calculate means and SE within each band ##
msat <- data.table(msat) #make data.table so can use data.table functions

#calculate mean in each band
chloroAmax_he_mean <- msat[, mean(He), by = chloroAmax_round]
  colnames(chloroAmax_he_mean) <- c("chloroAmax", "he_mean")

#calculate SE in each band
chloroAmax_he_SE <- msat[, std.error(He), by = chloroAmax_round]
  colnames(chloroAmax_he_SE) <- c("chloroAmax", "he_SE")

#count observations in each band
chloroAmax_he_count <- msat[, .N, by = chloroAmax_round]
  colnames(chloroAmax_he_count) <- c("chloroAmax", "he_count")

#merge mean and SE dataframes together
chloroAmax_he_binned_means <- list(chloroAmax_he_mean, chloroAmax_he_SE, chloroAmax_he_count) %>% 
  reduce(full_join, by = "chloroAmax")
chloroAmax_he_binned_means$chloroAmax <- as.numeric(as.character(chloroAmax_he_binned_means$chloroAmax))
  chloroAmax_he_binned_means <- chloroAmax_he_binned_means[order(chloroAmax), ]
  chloroAmax_he_binned_means$X <- chloroAmax_he_binned_means$chloroAmax + 5 #for plotting, plot in MIDDLE of band

#calculate error bars (standard error)
chloroAmax_he_binned_means$mean_lowerSE <- chloroAmax_he_binned_means$he_mean - 
  chloroAmax_he_binned_means$he_SE
chloroAmax_he_binned_means$mean_upperSE <- chloroAmax_he_binned_means$he_mean + 
  chloroAmax_he_binned_means$he_SE

#### Plot chloroAmax ####
#for legend
colors <- c("Binned means" = "#3F6DAA", "Regression" = "black")

#plot
msat_he_chloroAmax_plot_both <- ggplot() + 
  geom_line(data = chloroAmax_eff_data, 
            aes(x = chlomax, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = chloroAmax_eff_data, 
              aes(x = chlomax, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = chloroAmax_he_binned_means, 
             aes(x = X, y = he_mean, color = "Binned means", size = he_count), shape = "square") + 
  geom_errorbar(data = chloroAmax_he_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "Binned means"), width = 1.75, size = 3) + 
  xlab("Chloro A max") + ylab("msat He") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 1.0)) + 
  scale_size_continuous(breaks = c(50, 1000, 2000, 5000), 
                        range = c(10, 20))
msat_he_chloroAmax_plot_annotated_both <- msat_he_chloroAmax_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
msat_he_chloroAmax_plot_annotated_both

chlomax_lat_plot <- ggplot() + 
  geom_point(data = msat, 
             aes(x = abslat, y = chloroA.BO_chlomax)) + 
  theme_bw()

cor(msat$chloroA.BO_chlomax, msat$abslat)

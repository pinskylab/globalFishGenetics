######################## Script for Building Regression Figures for Pi  ##############################

#Uses final mtdna pi linear regression models
#Predicts pi at certain values of predictor variable of interest
#Calculates mean pi of raw data (binned every X units)
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
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)
mtdna_env <- read.csv("output/mtdna_env.csv", stringsAsFactors = FALSE)

mtdna <- merge(mtdna, mtdna_env[, c('X', 'sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 
                                    'sst.BO_sstmin', 'BO_dissox', 'chloroA.BO_chlomean', 
                                    'chloroA.BO_chlorange', 'chloroA.BO_chlomax', 
                                    'chloroA.BO_chlomin')], all.x = TRUE)
mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 
                                  'Southernmost', 'Half_RangeSize', 'Centroid')], all.x = TRUE)

######################################################################################################################

######## Clean up dataframe ########

#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)

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

#convert lon to radians
mtdna_small_pi$lon_360 <- mtdna_small_pi$lon + 180 #convert (-180,180) to (0,360)
  mtdna_small_pi$lon_rad <- (2*pi*mtdna_small_pi$lon_360)/360

#############################################################################################################

######### Abslat figures #######

abslat_model_pi <- lm(logpi ~ abslat_scale, data = mtdna_small_pi, 
                      na.action = "na.fail")

#### Predict ####
#marginal effects
abslat_eff <- plot_model(abslat_model_pi, type = "pred",
                         terms = "abslat_scale [all]")

#pull out marginal effects dataframe
abslat_eff_data <- as.data.frame(abslat_eff$data)

#unscale lat
#use same scaled:center & scaled:scale from original data
abslat_scale <- scale(mtdna_small_pi$abslat) #bc had to convert to numeric to run model/calculate marginal effects
abslat_eff_data$abslat <- (abslat_eff_data$x * attr(abslat_scale, "scaled:scale")) + 
  attr(abslat_scale, "scaled:center")

#unlog pi
abslat_eff_data$unlog_pi <- 10^(abslat_eff_data$predicted)
  abslat_eff_data$unlog_conf.low <- 10^(abslat_eff_data$conf.low)
  abslat_eff_data$unlog_conf.high <- 10^(abslat_eff_data$conf.high)

#### Calculate means from raw data ####
## Round abslat DOWN to nearest multiple of 10 ##
#create column to fill with the rounded abslat
mtdna_small_pi$abslat_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$abslat_round[i] <- DescTools::RoundTo(mtdna_small_pi$abslat[i], 
                                                        multiple = 10, FUN = floor)} #rounding DOWN (e.g. abslat 10-20 gets value of 10)
}

mtdna_small_pi$abslat_round <- as.factor(mtdna_small_pi$abslat_round)

## Calculate means and SE within each 10 degree band ##
mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
abslat_pi_mean <- mtdna_small_pi[, mean(logpi), by = abslat_round]
  colnames(abslat_pi_mean) <- c("abslat", "pi_mean")
  abslat_pi_mean$pi_mean <- 10^(abslat_pi_mean$pi_mean) #unlog mean

#calculate SE in each 10 degree band
abslat_pi_SE <- mtdna_small_pi[, std.error(Pi), by = abslat_round]
  colnames(abslat_pi_SE) <- c("abslat", "pi_SE")

#count observations in each 10 degree band
abslat_pi_count <- mtdna_small_pi[, .N, by = abslat_round]
  colnames(abslat_pi_count) <- c("abslat", "pi_count")

#merge mean and SE dataframes together
abslat_pi_binned_means <- list(abslat_pi_mean, abslat_pi_SE, abslat_pi_count) %>% 
  reduce(full_join, by = "abslat")
abslat_pi_binned_means$abslat <- as.numeric(as.character(abslat_pi_binned_means$abslat))
  abslat_pi_binned_means <- abslat_pi_binned_means[order(abslat), ]
  abslat_pi_binned_means$X <- abslat_pi_binned_means$abslat + 5 #for plotting, plot in MIDDLE of 10 degree band
  
#calculate error bars (standard error)
abslat_pi_binned_means$mean_lowerSE <- abslat_pi_binned_means$pi_mean - 
  abslat_pi_binned_means$pi_SE
abslat_pi_binned_means$mean_upperSE <- abslat_pi_binned_means$pi_mean + 
  abslat_pi_binned_means$pi_SE

#### Plot abslat ####
#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_abslat_plot_both <- ggplot() + 
  geom_line(data = abslat_eff_data, 
            aes(x = abslat, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = abslat_eff_data, 
              aes(x = abslat, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = abslat_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "10-degree binned means", size = pi_count), shape = "square") + 
  geom_errorbar(data = abslat_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("absolute latitude") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_y_continuous(limits = c(0, 0.0075)) + 
  scale_color_manual(values = colors) +
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
             range = c(10, 20))
mtdna_pi_abslat_plot_annotated_both <- mtdna_pi_abslat_plot_both + theme_bw() + coord_flip() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_abslat_plot_annotated_both

#############################################################################################################

######### Lat figures #######

lat_model_pi <- lmer(logpi ~ lat_scale + I(lat_scale^2), data = mtdna_small_pi, 
                     na.action = "na.fail")

#### Predict ####
#marginal effects
lat_eff <- plot_model(lat_model_pi, type = "pred",
                         terms = "lat_scale [all]")

#pull out marginal effects dataframe
lat_eff_data <- as.data.frame(lat_eff$data)

#unscale lat
#use same scaled:center & scaled:scale from original data
lat_scale <- scale(mtdna_small_pi$lat) #bc had to convert to numeric to run model/calculate marginal effects
lat_eff_data$lat <- (lat_eff_data$x * attr(lat_scale, "scaled:scale")) + 
  attr(lat_scale, "scaled:center")

#unlog pi
lat_eff_data$unlog_pi <- 10^(lat_eff_data$predicted)
lat_eff_data$unlog_conf.low <- 10^(lat_eff_data$conf.low)
lat_eff_data$unlog_conf.high <- 10^(lat_eff_data$conf.high)

#### Calculate means from raw data ####
## Round lat DOWN to nearest multiple of 10 ##
#create column to fill with the rounded abslat
mtdna_small_pi$lat_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$lat_round[i] <- DescTools::RoundTo(mtdna_small_pi$lat[i], 
                                                     multiple = 10, FUN = floor)} #rounding DOWN (e.g. lat 10-12 gets value of 10)
}

mtdna_small_pi$lat_round <- as.factor(mtdna_small_pi$lat_round)

##Calculate means and SE within each 10 degree band ##

mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
lat_pi_mean <- mtdna_small_pi[, mean(logpi), by = list(lat_round)]
  colnames(lat_pi_mean) <- c("lat", "pi_mean")
  lat_pi_mean$pi_mean <- 10^(lat_pi_mean$pi_mean)

#calculate SE in each 10 degree band
lat_pi_SE <- mtdna_small_pi[, std.error(Pi), by = list(lat_round)] 
  colnames(lat_pi_SE) <- c("lat", "pi_SE")

#count observations in each 10 degree band
lat_pi_count <- mtdna_small_pi[, .N, by = lat_round]
  colnames(lat_pi_count) <- c("lat", "pi_count")

#merge mean and SE dataframes together
lat_pi_binned_means <- list(lat_pi_mean, lat_pi_SE, lat_pi_count) %>% 
  reduce(full_join, by = "lat")
lat_pi_binned_means$lat <- as.numeric(as.character(lat_pi_binned_means$lat))
  lat_pi_binned_means <- lat_pi_binned_means[order(lat), ]
  lat_pi_binned_means$X <- lat_pi_binned_means$lat + 5 #for plotting, plot in MIDDLE of 10 degree band

#calculate error bars (standard error)
lat_pi_binned_means$mean_lowerSE <- lat_pi_binned_means$pi_mean - 
  lat_pi_binned_means$pi_SE
lat_pi_binned_means$mean_upperSE <- lat_pi_binned_means$pi_mean + 
  lat_pi_binned_means$pi_SE

#### Plot lat ####
#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_lat_plot_both <- ggplot() + 
  geom_line(data = lat_eff_data, 
            aes(x = lat, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = lat_eff_data, 
              aes(x = lat, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = lat_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "10-degree binned means", size = pi_count), shape = "square") + 
  geom_errorbar(data = lat_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("latitude") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_y_continuous(limits = c(0, 0.01)) + 
  scale_color_manual(values = colors) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_pi_lat_plot_annotated_both <- mtdna_pi_lat_plot_both + theme_bw() + coord_flip() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_lat_plot_annotated_both

#####################################################################################################

######## Lon figures ########

lon_model_pi <- lmer(logpi ~ sin(lon_rad) + cos(lon_rad), data = mtdna_small_pi, 
                     na.action = "na.fail")

#### Predict ####
#marginal effects
lon_eff <- plot_model(lon_model_pi, type = "pred", 
                         terms = "lon_rad [all]")

#pull out marginal effects dataframe
lon_eff_data <- as.data.frame(lon_eff$data)
  lon_eff_data$lon <- ((360*lon_eff_data$x)/(2*pi))-180

#unlog pi
lon_eff_data$unlog_pi <- 10^(lon_eff_data$predicted)
  lon_eff_data$unlog_conf.low <- 10^(lon_eff_data$conf.low)
  lon_eff_data$unlog_conf.high <- 10^(lon_eff_data$conf.high)
  
#### Calculate means from raw data ####
## Round lon DOWN to nearest multiple of 10 ##
#create column to fill with the rounded lon
mtdna_small_pi$lon_round <- NA

#fil round column in
for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$lon_round[i] <- DescTools::RoundTo(mtdna_small_pi$lon[i], 
                                                     multiple = 10, FUN = floor)} #rounding DOWN (e.g. lon 10-12 gets value of 10)
}

mtdna_small_pi$lon_round <- as.factor(mtdna_small_pi$lon_round)

## Calculate means and SE within each 10 degree band ##
mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
lon_pi_mean <- mtdna_small_pi[, mean(logpi), by = lon_round]
  colnames(lon_pi_mean) <- c("lon", "pi_mean")
  lon_pi_mean$pi_mean <- 10^(lon_pi_mean$pi_mean) #unlog mean

#calculate SE in each 10 degree band
lon_pi_SE <- mtdna_small_pi[, std.error(Pi), by = lon_round]
  colnames(lon_pi_SE) <- c("lon", "pi_SE")

#count observations in each 10 degree band
lon_pi_count <- mtdna_small_pi[, .N, by = lon_round]
  colnames(lon_pi_count) <- c("lon", "pi_count")

#merge mean and SE dataframes together
lon_pi_binned_means <- list(lon_pi_mean, lon_pi_SE, lon_pi_count) %>% 
  reduce(full_join, by = "lon")
  lon_pi_binned_means$lon <- as.numeric(as.character(lon_pi_binned_means$lon))
  lon_pi_binned_means <- lon_pi_binned_means[order(lon), ]
  lon_pi_binned_means$X <- lon_pi_binned_means$lon + 5 #for plotting, plot in MIDDLE of 10 degree band
    lon_pi_binned_means$X[lon_pi_binned_means$lon == 180] <- 180 #these two are on the 180 line
  
#calculate error bars (standard error)
lon_pi_binned_means$mean_lowerSE <- lon_pi_binned_means$pi_mean - 
  lon_pi_binned_means$pi_SE
  lon_pi_binned_means$mean_lowerSE[lon_pi_binned_means$mean_lowerSE < 0] <- 0 #capping at 0 if goes negative
lon_pi_binned_means$mean_upperSE <- lon_pi_binned_means$pi_mean + 
  lon_pi_binned_means$pi_SE
  lon_pi_binned_means$mean_upperSE[lon_pi_binned_means$mean_upperSE > 0.015] <- 0.015 #capping at 0.015 for plotting

#### Plot lon ####
#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_lon_plot_both <- ggplot() + 
  geom_line(data = lon_eff_data, 
            aes(x = lon, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = lon_eff_data, 
              aes(x = lon, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = lon_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "10-degree binned means", size = pi_count), shape = "square") + 
  geom_errorbar(data = lon_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("longitude") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 0.015)) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_pi_lon_plot_annotated_both <- mtdna_pi_lon_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_lon_plot_annotated_both

#############################################################################################################

######## Calculate environmental variables ########

#### log transform sst data ####
#subset to only those with sst data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$sst.BO_sstmean != "NA") #shouldn't remove any

mtdna_small_pi$logsstmean <- log10(mtdna_small_pi$sst.BO_sstmean)
  mtdna_small_pi$logsstrange <- log10(mtdna_small_pi$sst.BO_sstrange)
  mtdna_small_pi$logsstmax <- log10(mtdna_small_pi$sst.BO_sstmax)
  mtdna_small_pi$logsstmin <- log10(mtdna_small_pi$sst.BO_sstmin)

#remove logsst = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstmean != "NaN")
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstrange != "NaN")
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstmax != "NaN")
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstmin != "NaN")

#### log transform dissox data ####
#subset to only those with dissox data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$BO_dissox != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

mtdna_small_pi$logdissox <- log10(mtdna_small_pi$BO_dissox)

#remove logdissox = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logdissox != "NaN")

#### log transform chlorophyll A ####
#subset to only those with chloroA data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlorange != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomax != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomin != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

mtdna_small_pi$logchlomean <- log10(mtdna_small_pi$sst.BO_sstmean)
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

#######################################################################################################

######## SST mean figures ########

SSTmean_model_pi <- lm(logpi ~ logsstmean, data = mtdna_small_pi, 
                        na.action = "na.fail")

#### Predict ####
#marginal effects
SSTmean_eff <- plot_model(SSTmean_model_pi, type = "pred", 
                      terms = "logsstmean [all]")

#pull out marginal effects dataframe
SSTmean_eff_data <- as.data.frame(SSTmean_eff$data)
  SSTmean_eff_data$sstmean <- 10^(SSTmean_eff_data$x)

#unlog pi
SSTmean_eff_data$unlog_pi <- 10^(SSTmean_eff_data$predicted)
  SSTmean_eff_data$unlog_conf.low <- 10^(SSTmean_eff_data$conf.low)
  SSTmean_eff_data$unlog_conf.high <- 10^(SSTmean_eff_data$conf.high)

#### Calculate means from raw data ####
## Round SSTmean DOWN to nearest multiple of 5 ##
#create column to fill with the rounded SSTmean
mtdna_small_pi$SSTmean_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$SSTmean_round[i] <- DescTools::RoundTo(mtdna_small_pi$sst.BO_sstmean[i], 
                                                        multiple = 5, FUN = floor)} #rounding DOWN (e.g. SSTmean 5-10 gets value of 5)
}

mtdna_small_pi$SSTmean_round <- as.factor(mtdna_small_pi$SSTmean_round)

## Calculate means and SE within each 5 degree band ##
mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so can use data.table functions

#calculate mean in each 5 degree band
SSTmean_pi_mean <- mtdna_small_pi[, mean(logpi), by = SSTmean_round]
  colnames(SSTmean_pi_mean) <- c("SSTmean", "pi_mean")
  SSTmean_pi_mean$pi_mean <- 10^(SSTmean_pi_mean$pi_mean) #unlog mean

#calculate SE in each 5 degree band
SSTmean_pi_SE <- mtdna_small_pi[, std.error(Pi), by = SSTmean_round]
  colnames(SSTmean_pi_SE) <- c("SSTmean", "pi_SE")

#count observations in each 5 degree band
SSTmean_pi_count <- mtdna_small_pi[, .N, by = SSTmean_round]
  colnames(SSTmean_pi_count) <- c("SSTmean", "pi_count")

#merge mean and SE dataframes together
SSTmean_pi_binned_means <- list(SSTmean_pi_mean, SSTmean_pi_SE, SSTmean_pi_count) %>% 
  reduce(full_join, by = "SSTmean")
  SSTmean_pi_binned_means$SSTmean <- as.numeric(as.character(SSTmean_pi_binned_means$SSTmean))
  SSTmean_pi_binned_means <- SSTmean_pi_binned_means[order(SSTmean), ]
  SSTmean_pi_binned_means$X <- SSTmean_pi_binned_means$SSTmean + 2.5 #for plotting, plot in MIDDLE of 5 degree band

#calculate error bars (standard error)
SSTmean_pi_binned_means$mean_lowerSE <- SSTmean_pi_binned_means$pi_mean - 
  SSTmean_pi_binned_means$pi_SE
SSTmean_pi_binned_means$mean_upperSE <- SSTmean_pi_binned_means$pi_mean + 
  SSTmean_pi_binned_means$pi_SE

#### Plot SSTmean ####
#for legend
colors <- c("5-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_SSTmean_plot_both <- ggplot() + 
  geom_line(data = SSTmean_eff_data, 
            aes(x = sstmean, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = SSTmean_eff_data, 
              aes(x = sstmean, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = SSTmean_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "5-degree binned means", size = pi_count), shape = "square") + 
  geom_errorbar(data = SSTmean_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "5-degree binned means"), width = 1.75, size = 3) + 
  xlab("SST mean") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 0.012)) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_pi_SSTmean_plot_annotated_both <- mtdna_pi_SSTmean_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_SSTmean_plot_annotated_both

#############################################################################################################

######## SST min models ########

SSTmin_model_pi <- lmer(logpi ~ logsstmin + (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                        REML = FALSE, data = mtdna_small_pi, 
                                   na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
SSTmin_eff <- plot_model(SSTmin_model_pi, type = "pred", 
                          terms = "logsstmin [all]")

#pull out marginal effects dataframe
SSTmin_eff_data <- as.data.frame(SSTmin_eff$data)
SSTmin_eff_data$sstmin <- 10^(SSTmin_eff_data$x)

#unlog pi
SSTmin_eff_data$unlog_pi <- 10^(SSTmin_eff_data$predicted)
SSTmin_eff_data$unlog_conf.low <- 10^(SSTmin_eff_data$conf.low)
SSTmin_eff_data$unlog_conf.high <- 10^(SSTmin_eff_data$conf.high)

#### Calculate means from raw data ####
## Round SSTmin DOWN to nearest multiple of 5 ##
#create column to fill with the rounded SSTmin
mtdna_small_pi$SSTmin_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$SSTmin_round[i] <- DescTools::RoundTo(mtdna_small_pi$sst.BO_sstmin[i], 
                                                         multiple = 5, FUN = floor)} #rounding DOWN (e.g. SSTmin 5-10 gets value of 5)
}

mtdna_small_pi$SSTmin_round <- as.factor(mtdna_small_pi$SSTmin_round)

## Calculate means and SE within each 5 degree band ##
mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so can use data.table functions

#calculate mean in each 5 degree band
SSTmin_pi_mean <- mtdna_small_pi[, mean(logpi), by = SSTmin_round]
colnames(SSTmin_pi_mean) <- c("SSTmin", "pi_mean")
SSTmin_pi_mean$pi_mean <- 10^(SSTmin_pi_mean$pi_mean) #unlog mean

#calculate SE in each 5 degree band
SSTmin_pi_SE <- mtdna_small_pi[, std.error(Pi), by = SSTmin_round]
colnames(SSTmin_pi_SE) <- c("SSTmin", "pi_SE")

#count observations in each 5 degree band
SSTmin_pi_count <- mtdna_small_pi[, .N, by = SSTmin_round]
colnames(SSTmin_pi_count) <- c("SSTmin", "pi_count")

#merge mean and SE dataframes together
SSTmin_pi_binned_means <- list(SSTmin_pi_mean, SSTmin_pi_SE, SSTmin_pi_count) %>% 
  reduce(full_join, by = "SSTmin")
SSTmin_pi_binned_means$SSTmin <- as.numeric(as.character(SSTmin_pi_binned_means$SSTmin))
  SSTmin_pi_binned_means <- SSTmin_pi_binned_means[order(SSTmin), ]
  SSTmin_pi_binned_means$X <- SSTmin_pi_binned_means$SSTmin + 2.5 #for plotting, plot in MIDDLE of 5 degree band

#calculate error bars (standard error)
SSTmin_pi_binned_means$mean_lowerSE <- SSTmin_pi_binned_means$pi_mean - 
  SSTmin_pi_binned_means$pi_SE
SSTmin_pi_binned_means$mean_upperSE <- SSTmin_pi_binned_means$pi_mean + 
  SSTmin_pi_binned_means$pi_SE

#### Plot SSTmin ####
#for legend
colors <- c("5-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_SSTmin_plot_both <- ggplot() + 
  geom_line(data = SSTmin_eff_data, 
            aes(x = sstmin, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = SSTmin_eff_data, 
              aes(x = sstmin, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = SSTmin_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "5-degree binned means", size = pi_count), shape = "square") + 
  geom_errorbar(data = SSTmin_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "5-degree binned means"), width = 1.75, size = 3) + 
  xlab("SST min") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 0.012)) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_pi_SSTmin_plot_annotated_both <- mtdna_pi_SSTmin_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_SSTmin_plot_annotated_both

###########################################################################################################

######## Dissox figures ########

dissox_model_pi <- lmer(logpi ~ logdissox + (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                        REML = FALSE, data = mtdna_small_pi, 
                        na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#### Predict ####
#marginal effects
dissox_eff <- plot_model(dissox_model_pi, type = "pred", 
                         terms = "logdissox [all]")

#pull out marginal effects dataframe
dissox_eff_data <- as.data.frame(dissox_eff$data)
dissox_eff_data$dissox <- 10^(dissox_eff_data$x)

#unlog pi
dissox_eff_data$unlog_pi <- 10^(dissox_eff_data$predicted)
dissox_eff_data$unlog_conf.low <- 10^(dissox_eff_data$conf.low)
dissox_eff_data$unlog_conf.high <- 10^(dissox_eff_data$conf.high)

#### Calculate means from raw data ####
## Round dissox DOWN to nearest multiple of 5 ##
#create column to fill with the rounded dissox
mtdna_small_pi$dissox_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$dissox_round[i] <- DescTools::RoundTo(mtdna_small_pi$BO_dissox[i], 
                                                        multiple = 0.5, FUN = floor)} #rounding DOWN
}

  mtdna_small_pi$dissox_round <- as.factor(mtdna_small_pi$dissox_round)

## Calculate means and SE within each 0.5 band ##
mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so that can use data.table functions

#calculate mean in each 0.5 band
dissox_pi_mean <- mtdna_small_pi[, mean(logpi), by = dissox_round]
  colnames(dissox_pi_mean) <- c("dissox", "logpi_mean")
  dissox_pi_mean$pi_mean <- 10^(dissox_pi_mean$logpi_mean) #unlog mean

#calculate SE in each 0.5 band  
dissox_pi_SE <- mtdna_small_pi[, std.error(Pi), by = dissox_round]
  colnames(dissox_pi_SE) <- c("dissox", "pi_SE")
  
#count observations in each 0.5 band
dissox_pi_count <- mtdna_small_pi[, .N, by = dissox_round]
  colnames(dissox_pi_count) <- c("dissox", "pi_count")

#merge mean and SE dataframes together
dissox_pi_binned_means <- list(dissox_pi_mean, dissox_pi_SE, dissox_pi_count) %>% 
  reduce(full_join, by = "dissox")
  dissox_pi_binned_means$dissox <- as.numeric(as.character(dissox_pi_binned_means$dissox))
  dissox_pi_binned_means <- dissox_pi_binned_means[order(dissox), ]
  dissox_pi_binned_means$X <- dissox_pi_binned_means$dissox + 0.25 #for plotting, plot in MIDDLE of 0.5 band
  
#add error bars (standard error)
dissox_pi_binned_means$mean_lowerSE <- dissox_pi_binned_means$pi_mean - 
  dissox_pi_binned_means$pi_SE
dissox_pi_binned_means$mean_upperSE <- dissox_pi_binned_means$pi_mean + 
  dissox_pi_binned_means$pi_SE

#### Plot dissox ####
#for legend
colors <- c("Binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_dissox_plot_both <- ggplot() + 
  geom_line(data = dissox_eff_data, 
            aes(x = dissox, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = dissox_eff_data, 
              aes(x = dissox, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = dissox_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "Binned means", size = pi_count), shape = "square") + 
  geom_errorbar(data = dissox_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "Binned means"), width = 0.5, size = 3) + 
  xlab("Dissolved Oxygen Mean") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 0.012)) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_pi_dissox_plot_annotated_both <- mtdna_pi_dissox_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_dissox_plot_annotated_both

#############################################################################################################

######## ChloroA mean figures ########

chloroAmean_model_pi <- lmer(logpi ~ logchlomean + I(logchlomean^2) + (1|Family/Genus) + (1|Source) + 
                              (1|MarkerName), REML = FALSE, data = mtdna_small_pi, 
                            na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#### Predict ####
#marginal effects
chloroAmean_eff <- plot_model(chloroAmean_model_pi, type = "pred", 
                         terms = "logchlomean [all]")

#pull out marginal effects dataframe
chloroAmean_eff_data <- as.data.frame(chloroAmean_eff$data)
  chloroAmean_eff_data$chlomean <- 10^(chloroAmean_eff_data$x)

#unlog pi
chloroAmean_eff_data$unlog_pi <- 10^(chloroAmean_eff_data$predicted)
  chloroAmean_eff_data$unlog_conf.low <- 10^(chloroAmean_eff_data$conf.low)
  chloroAmean_eff_data$unlog_conf.high <- 10^(chloroAmean_eff_data$conf.high)

#### Calculate means from raw data ####
## Round chlomean DOWN to nearest multiple of 10 ##
#create column to fill with the rounded chlomean
mtdna_small_pi$chloroAmean_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$chloroAmean_round[i] <- DescTools::RoundTo(mtdna_small_pi$chloroA.BO_chlomean[i], 
                                                        multiple = 1, FUN = floor)} #rounding DOWN
}

  mtdna_small_pi$chloroAmean_round <- as.factor(mtdna_small_pi$chloroAmean_round)

## Calculate means and SE within each band ##
mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so can use data.table functions

#calculate mean in each band
chloroAmean_pi_mean <- mtdna_small_pi[, mean(logpi), by = chloroAmean_round]
  colnames(chloroAmean_pi_mean) <- c("chloroAmean", "pi_mean")
  chloroAmean_pi_mean$pi_mean <- 10^(chloroAmean_pi_mean$pi_mean) #unlog mean

#calculate SE in each band
  chloroAmean_pi_SE <- mtdna_small_pi[, std.error(Pi), by = chloroAmean_round]
  colnames(chloroAmean_pi_SE) <- c("chloroAmean", "pi_SE")

#count observations in each band
chloroAmean_pi_count <- mtdna_small_pi[, .N, by = chloroAmean_round]
  colnames(chloroAmean_pi_count) <- c("chloroAmean", "pi_count")

#merge mean and SE dataframes together
chloroAmean_pi_binned_means <- list(chloroAmean_pi_mean, chloroAmean_pi_SE, chloroAmean_pi_count) %>% 
  reduce(full_join, by = "chloroAmean")
  chloroAmean_pi_binned_means$chloroAmean <- as.numeric(as.character(chloroAmean_pi_binned_means$chloroAmean))
  chloroAmean_pi_binned_means <- chloroAmean_pi_binned_means[order(chloroAmean), ]
  chloroAmean_pi_binned_means$X <- chloroAmean_pi_binned_means$chloroAmean + 0.5 #for plotting, plot in MIDDLE of band

#calculate error bars (standard error)
chloroAmean_pi_binned_means$mean_lowerSE <- chloroAmean_pi_binned_means$pi_mean - 
  chloroAmean_pi_binned_means$pi_SE
chloroAmean_pi_binned_means$mean_upperSE <- chloroAmean_pi_binned_means$pi_mean + 
  chloroAmean_pi_binned_means$pi_SE

#### Plot chloroAmean ####
#for legend
colors <- c("Binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_chloroAmean_plot_both <- ggplot() + 
  geom_line(data = chloroAmean_eff_data, 
            aes(x = chlomean, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = chloroAmean_eff_data, 
              aes(x = chlomean, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = chloroAmean_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "Binned means", size = pi_count), shape = "square") + 
  geom_errorbar(data = chloroAmean_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "Binned means"), width = 1.75, size = 3) + 
  xlab("Chloro A mean") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 0.012)) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_pi_chloroAmean_plot_annotated_both <- mtdna_pi_chloroAmean_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_chloroAmean_plot_annotated_both

#############################################################################################################

######## ChloroA range figures ########

chloroArange_model_pi <- lmer(logpi ~ logchlorange + I(logchlorange^2) + (1|Family/Genus) + (1|Source) + 
                               (1|MarkerName), REML = FALSE, data = mtdna_small_pi, 
                             na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#### Predict ####
#marginal effects
chloroArange_eff <- plot_model(chloroArange_model_pi, type = "pred", 
                              terms = "logchlorange [all]")

#pull out marginal effects dataframe
chloroArange_eff_data <- as.data.frame(chloroArange_eff$data)
chloroArange_eff_data$chlorange <- 10^(chloroArange_eff_data$x)

#unlog pi
chloroArange_eff_data$unlog_pi <- 10^(chloroArange_eff_data$predicted)
  chloroArange_eff_data$unlog_conf.low <- 10^(chloroArange_eff_data$conf.low)
  chloroArange_eff_data$unlog_conf.high <- 10^(chloroArange_eff_data$conf.high)

#### Calculate means from raw data ####
## Round chlomrange DOWN to nearest multiple of 10 ##
#create column to fill with the rounded chlomean
mtdna_small_pi$chloroArange_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$chloroArange_round[i] <- DescTools::RoundTo(mtdna_small_pi$chloroA.BO_chlorange[i], 
                                                             multiple = 1, FUN = floor)} #rounding DOWN
}

  mtdna_small_pi$chloroArange_round <- as.factor(mtdna_small_pi$chloroArange_round)

## Calculate means and SE within each band ##
mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so can use data.table functions

#calculate mean in each band
chloroArange_pi_mean <- mtdna_small_pi[, mean(logpi), by = chloroArange_round]
  colnames(chloroArange_pi_mean) <- c("chloroArange", "pi_mean")
  chloroArange_pi_mean$pi_mean <- 10^(chloroArange_pi_mean$pi_mean) #unlog mean

#calculate SE in each band
chloroArange_pi_SE <- mtdna_small_pi[, std.error(Pi), by = chloroArange_round]
  colnames(chloroArange_pi_SE) <- c("chloroArange", "pi_SE")

#count observations in each band
chloroArange_pi_count <- mtdna_small_pi[, .N, by = chloroArange_round]
  colnames(chloroArange_pi_count) <- c("chloroArange", "pi_count")

#merge mean and SE dataframes together
chloroArange_pi_binned_means <- list(chloroArange_pi_mean, chloroArange_pi_SE, chloroArange_pi_count) %>% 
  reduce(full_join, by = "chloroArange")
  chloroArange_pi_binned_means$chloroArange <- as.numeric(as.character(chloroArange_pi_binned_means$chloroArange))
  chloroArange_pi_binned_means <- chloroArange_pi_binned_means[order(chloroArange), ]
  chloroArange_pi_binned_means$X <- chloroArange_pi_binned_means$chloroArange + 0.5 #for plotting, plot in MIDDLE of band

#calculate error bars (standard error)
chloroArange_pi_binned_means$mean_lowerSE <- chloroArange_pi_binned_means$pi_mean - 
  chloroArange_pi_binned_means$pi_SE
chloroArange_pi_binned_means$mean_upperSE <- chloroArange_pi_binned_means$pi_mean + 
  chloroArange_pi_binned_means$pi_SE

#### Plot chloroArange ####
#for legend
colors <- c("Binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_chloroArange_plot_both <- ggplot() + 
  geom_line(data = chloroArange_eff_data, 
            aes(x = chlorange, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = chloroArange_eff_data, 
              aes(x = chlorange, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = chloroArange_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "Binned means", size = pi_count), shape = "square") + 
  geom_errorbar(data = chloroArange_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "Binned means"), width = 1.75, size = 3) + 
  xlab("Chloro A range") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 0.012)) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_pi_chloroArange_plot_annotated_both <- mtdna_pi_chloroArange_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_chloroArange_plot_annotated_both
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
msat <- read.csv("output/msatloci_assembled.csv", stringsAsFactors = FALSE) #should this be msatloci??
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)
msat_env <- read.csv("output/msat_env.csv", stringsAsFactors = FALSE) 

msat <- merge(msat, msat_env[, c('X', 'sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 
                                    'sst.BO_sstmin', 'BO_dissox', 'chloroA.BO_chlomean', 
                                    'chloroA.BO_chlorange', 'chloroA.BO_chlomax', 
                                    'chloroA.BO_chlomin')], all.x = TRUE)
msat <- merge(msat, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 
                                  'Southernmost', 'Half_RangeSize', 'Centroid')], all.x = TRUE)

###################################################################################################

######### Clean up dataframe ########

#remove missing n values
msat$n <- as.numeric(msat$n)
msat <- subset(msat, msat$n != "NA")

#remove missing He
msat <- subset(msat, msat$He != "NA")

#calculate abslat
msat$abslat <- abs(msat$lat)

###################################################################################################

######## Build lat, abslat & lon models ########

#### clean up data ####

#scale geographic variables
msat$lat_scale <- as.numeric(scale(msat$lat))
msat$abslat_scale <- as.numeric(scale(msat$abslat))

#convert lon to radians
msat$lon_360 <- msat$lon + 180 #convert (-180,180) to (0,360)
msat$lon_rad <- (2*pi*msat$lon_360)/360

#calculate success and failure
msat$success <- round(msat$He*msat$n) #number of heterozygotes
msat$failure<- round((1 - msat$He)*msat$n) #number of homozygotes

#### Models ####
#for prediction plots, do NOT want any random effects
#also not using different (faster) optimizer

lat_model_he <- glm(cbind(success, failure) ~ PrimerNote + CrossSpp + lat_scale + I(lat_scale^2) + I(lat_scale^3) + 
                      I(lat_scale^4) + sin(lon_rad) + cos(lon_rad), family = binomial, data = msat, 
                    na.action = "na.omit")

lat_CP_model_he <- glm(cbind(success, failure) ~ bp_scale + range_position + lat_scale + 
                         I(lat_scale^2) + Pelagic_Coastal + Pelagic_Coastal:lat_scale + 
                         Pelagic_Coastal:I(lat_scale^2), family = binomial, data = mtdna_small_hd, 
                       na.action = "na.fail")

abslat_model_he <- glm(cbind(success, failure) ~ PrimerNote + CrossSpp + abslat_scale + I(abslat_scale^2) + 
                         I(abslat_scale^3) + I(abslat_scale^4), family = binomial, data = msat, 
                       na.action = "na.omit")

abslat_CP_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + 
                            Pelagic_Coastal + abslat_scale + I(abslat_scale^2) + 
                            Pelagic_Coastal:abslat_scale + Pelagic_Coastal:I(abslat_scale^2), 
                          family = binomial, data = mtdna_small_hd, na.action = "na.fail")

lon_model_he <- glm(cbind(success, failure) ~ PrimerNote + CrossSpp + sin(lon_rad) + cos(lon_rad), 
                    family = binomial, data = msat, 
                    na.action = "na.omit")

lon_CP_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + sin(lon_rad) + 
                         cos(lon_rad) + Pelagic_Coastal + Pelagic_Coastal:sin(lon_rad) + 
                         Pelagic_Coastal:cos(lon_rad), family = binomial, 
                       data = mtdna_small_hd, na.action = "na.fail")

lat_lon_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + lat_scale + 
                          I(lat_scale^2) + sin(lon_rad) + 
                          cos(lon_rad), family = binomial, data = mtdna_small_hd, 
                        na.action = "na.fail")

lat_lon_CP_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + lat_scale + 
                             I(lat_scale^2) + sin(lon_rad) + 
                             cos(lon_rad) + Pelagic_Coastal + Pelagic_Coastal:lat_scale + 
                             Pelagic_Coastal:I(lat_scale^2) + Pelagic_Coastal:sin(lon_rad) + 
                             Pelagic_Coastal:cos(lon_rad), family = binomial, 
                           data = mtdna_small_hd, na.action = "na.fail")

#logtransform sst data
msat <- subset(msat, msat$sst.BO_sstmean != "NA") #subset to only rows with sstmean data

msat <- subset(msat, msat$sst.BO_sstmean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
msat$logsstmean <- log10(msat$sst.BO_sstmean)
msat <- subset(msat, msat$logsstmean != "NaN")

sstmean_model_he <- glm(cbind(success, failure) ~ PrimerNote + CrossSpp + logsstmean + I(logsstmean^2), 
                        family = binomial, data = msat, na.action = "na.omit")

sstmean_CP_model_he <- glm(cbind(success, failure) ~ bp_scale + range_position + logsstmean + Pelagic_Coastal + 
                             Pelagic_Coastal:logsstmean, family = binomial, data = mtdna_small_hd, na.action = "na.fail")

#############################################################################################################

######### Abs Lat model figures #######

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
CP_abslat_he_mean <- msat[, mean(He), by = list(abslat_round, Pelagic_Coastal)]
colnames(abslat_he_mean) <- c("abslat", "he_mean")
colnames(CP_abslat_he_mean) <- c("abslat", "Pelagic_Coastal", "hd_mean")

#calculate SE in each 10 degree band
abslat_he_SE <- msat[, std.error(He), by = abslat_round]
CP_abslat_hd_SE <- msat[, std.error(He), by = list(abslat_round, Pelagic_Coastal)]
colnames(abslat_he_SE) <- c("abslat", "he_SE")
colnames(CP_abslat_hd_SE) <- c("abslat", "Pelagic_Coastal", "hd_SE")

#count observations in each 10 degree band
abslat_he_count <- msat[, .N, by = abslat_round]
CP_abslat_hd_count <- msat[, .N, by = list(abslat_round, Pelagic_Coastal)]
colnames(abslat_he_count) <- c("abslat", "he_count")
colnames(CP_abslat_hd_count) <- c("abslat", "Pelagic_Coastal", "hd_count")

#merge mean and SE dataframes together
abslat_he_binned_means <- list(abslat_he_mean, abslat_he_SE, abslat_he_count) %>% 
  reduce(full_join, by = "abslat")
abslat_he_binned_means$abslat <- as.numeric(as.character(abslat_he_binned_means$abslat))
abslat_he_binned_means <- abslat_he_binned_means[order(abslat), ]
abslat_he_binned_means$X <- abslat_he_binned_means$abslat + 5 #for plotting, plot in MIDDLE of 10 degree band

CP_lat_hd_binned_means <- merge(CP_lat_hd_mean, CP_lat_hd_SE, by = c("lat", "Pelagic_Coastal"))
CP_lat_hd_binned_means <- merge(CP_lat_hd_binned_means, CP_lat_hd_count, by = c("lat", "Pelagic_Coastal"))
colnames(CP_lat_hd_binned_means) <- c("lat", "Pelagic_Coastal", "hd_mean", "hd_SE", "hd_count")
CP_lat_hd_binned_means$lat <- as.numeric(as.character(CP_lat_hd_binned_means$lat))
CP_lat_hd_binned_means <- CP_lat_hd_binned_means[order(lat), ]
CP_lat_hd_binned_means$X <- CP_lat_hd_binned_means$lat + 5 #for plotting, plot in MIDDLE of 10 degree band

#calculate error bars (standard error)
abslat_he_binned_means$mean_lowerSE <- abslat_he_binned_means$he_mean - 
  abslat_he_binned_means$he_SE
abslat_he_binned_means$mean_upperSE <- abslat_he_binned_means$he_mean + 
  abslat_he_binned_means$he_SE

CP_lat_hd_binned_means$mean_lowerSE <- CP_lat_hd_binned_means$hd_mean - 
  CP_lat_hd_binned_means$hd_SE
CP_lat_hd_binned_means$mean_upperSE <- CP_lat_hd_binned_means$hd_mean + 
  CP_lat_hd_binned_means$hd_SE

#add count bin for plotting (all < number, 2001 = >2000)
abslat_he_binned_means$count_bin <- c(2000, 4000, 5000, 5001, 5000, 5001, 4000, 2000)
abslat_he_binned_means$count_bin <- factor(abslat_he_binned_means$count_bin, levels = c(2000, 4000, 5000, 5001))

CP_lat_hd_binned_means$count_bin <- c(5, 5, 5, 5, 50, 50, 50, 50, 100, 50, 200, 5, 100, 50, 200, 
                                      50, 200, 50, 201, 50, 201, 50, 200, 50, 100, 50, 50, 50, 50)
CP_lat_hd_binned_means$count_bin <- factor(CP_lat_hd_binned_means$count_bin, levels = c(5, 50, 100, 200, 201))

#### Plot abslat ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
msat_he_abslat_plot_both <- ggplot() + 
  geom_line(data = abslat_eff_data, 
            aes(x = abslat, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = abslat_eff_data, 
              aes(x = abslat, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = abslat_he_binned_means, 
             aes(x = X, y = he_mean, color = "10-degree binned means", size = count_bin), shape = "square") + 
  geom_errorbar(data = abslat_he_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("absolute latitude") + ylab("nuclear He") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_size_manual(values = c(8, 12, 14, 16), 
                    labels = c("\u22642000", "\u22644000", "\u22645000", "\u22655000")) 
  #scale_y_continuous(limits = c(0.65, 0.775))
msat_he_abslat_plot_annotated_both <- msat_he_abslat_plot_both + theme_bw() + coord_flip() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 120, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "none", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
msat_he_abslat_plot_annotated_both

## With pelagic and coastal separate ##
#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_hd_lat_CP_plot_both <- ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey", size = 3) + 
  geom_smooth(data = CP_hd_lat_predict_df, 
              aes(x = X, y = predict_hd, color = Pelagic_Coastal), size = 6) + 
  geom_ribbon(data = CP_hd_lat_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = CP_lat_hd_binned_means, 
             aes(x = X, y = hd_mean, color = Pelagic_Coastal, size = count_bin), shape = "square") + 
  geom_errorbar(data = CP_lat_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 1.75, size = 3) + 
  xlab("latitude") + ylab("mtdna Hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP) + 
  scale_size_manual(values = c(6, 8, 12, 14, 16), 
                    labels = c("\u22645", "\u226450", "\u2264100", "\u2264200", "\u2265200")) + 
  scale_y_continuous(c(0.55, 0.8))
mtdna_hd_lat_CP_plot_annotated_both <- mtdna_hd_lat_CP_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 70, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 60, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 50),
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_lat_CP_plot_annotated_both

#############################################################################################################

######### Lat model figures #######

#### Predict ####

#marginal effects
lat_eff <- plot_model(lat_model_he, type = "pred", 
                      terms = "lat_scale [all]")

#pull out marginal effects dataframe
lat_eff_data <- as.data.frame(lat_eff$data)

#unscale lat
#use same scaled:center & scaled:scale from original data
lat_scale <- scale(msat$lat) #sometimes, original scale in dataframe doesn't store attributes
lat_eff_data$lat <- (lat_eff_data$x * attr(lat_scale, "scaled:scale")) + 
  attr(lat_scale, "scaled:center")

#### Calculate means from raw data ####

## Round lat DOWN to nearest multiple of 10 ##

#create column to fill with the rounded abslat
msat$lat_round <- NA

#fill round column in
for(i in 1:nrow(msat)) {
  cat(paste(i, " ", sep = ''))
  {msat$lat_round[i] <- DescTools::RoundTo(msat$lat[i], 
                                                   multiple = 10, FUN = floor)} #rounding DOWN (e.g. abslat 10-20 gets value of 10)
}

msat$lat_round <- as.factor(msat$lat_round)

## Calculate means and SE within each 10 degree band ##

msat <- data.table(msat) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
lat_he_mean <- msat[, mean(He), by = lat_round]
CP_lat_he_mean <- msat[, mean(He), by = list(lat_round, Pelagic_Coastal)]
  colnames(lat_he_mean) <- c("lat", "he_mean")
  colnames(CP_lat_he_mean) <- c("lat", "Pelagic_Coastal", "he_mean")

#calculate SE in each 10 degree band
lat_he_SE <- msat[, std.error(He), by = lat_round]
CP_lat_he_SE <- msat[, std.error(He), by = list(lat_round, Pelagic_Coastal)]
  colnames(lat_he_SE) <- c("lat", "he_SE")
  colnames(CP_lat_he_SE) <- c("lat", "Pelagic_Coastal", "he_SE")

#count observations in each 10 degree band
lat_he_count <- msat[, .N, by = lat_round]
CP_lat_he_count <- msat[, .N, by = list(lat_round, Pelagic_Coastal)]
  colnames(lat_he_count) <- c("lat", "he_count")
  colnames(CP_lat_he_count) <- c("lat", "Pelagic_Coastal", "he_count")

#merge mean and SE dataframes together
lat_he_binned_means <- list(lat_he_mean, lat_he_SE, lat_he_count) %>% 
  reduce(full_join, by = "lat")
  lat_he_binned_means$lat <- as.numeric(as.character(lat_he_binned_means$lat))
  lat_he_binned_means <- lat_he_binned_means[order(lat), ]
  lat_he_binned_means$X <- lat_he_binned_means$lat + 5 #for plotting, plot in MIDDLE of 10 degree band

CP_lat_hd_binned_means <- merge(CP_lat_hd_mean, CP_lat_hd_SE, by = c("lat", "Pelagic_Coastal"))
CP_lat_hd_binned_means <- merge(CP_lat_hd_binned_means, CP_lat_hd_count, by = c("lat", "Pelagic_Coastal"))
colnames(CP_lat_hd_binned_means) <- c("lat", "Pelagic_Coastal", "hd_mean", "hd_SE", "hd_count")
CP_lat_hd_binned_means$lat <- as.numeric(as.character(CP_lat_hd_binned_means$lat))
CP_lat_hd_binned_means <- CP_lat_hd_binned_means[order(lat), ]
CP_lat_hd_binned_means$X <- CP_lat_hd_binned_means$lat + 5 #for plotting, plot in MIDDLE of 10 degree band

#calculate error bars (standard error)
lat_he_binned_means$mean_lowerSE <- lat_he_binned_means$he_mean - 
  lat_he_binned_means$he_SE
lat_he_binned_means$mean_upperSE <- lat_he_binned_means$he_mean + 
  lat_he_binned_means$he_SE

CP_lat_hd_binned_means$mean_lowerSE <- CP_lat_hd_binned_means$hd_mean - 
  CP_lat_hd_binned_means$hd_SE
CP_lat_hd_binned_means$mean_upperSE <- CP_lat_hd_binned_means$hd_mean + 
  CP_lat_hd_binned_means$hd_SE

#add count bin for plotting (all < number, 5001 = >5000)
lat_he_binned_means$count_bin <- c(50, 500, 500, 500, 2000, 2000, 2000, 2000, 2000,
                                   2000, 5000, 5001, 5000, 5001, 5000, 500)
lat_he_binned_means$count_bin <- factor(lat_he_binned_means$count_bin, levels = c(50, 500, 2000, 5000, 5001))

CP_lat_hd_binned_means$count_bin <- c(5, 5, 5, 5, 50, 50, 50, 50, 100, 50, 200, 5, 100, 50, 200, 
                                      50, 200, 50, 201, 50, 201, 50, 200, 50, 100, 50, 50, 50, 50)
CP_lat_hd_binned_means$count_bin <- factor(CP_lat_hd_binned_means$count_bin, levels = c(5, 50, 100, 200, 201))

#### Plot lat ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
msat_he_lat_plot_both <- ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey", size = 3) + 
  geom_line(data = lat_eff_data, 
            aes(x = lat, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = lat_eff_data, 
              aes(x = lat, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = lat_he_binned_means, 
             aes(x = X, y = he_mean, color = "10-degree binned means", size = count_bin), shape = "square") + 
  geom_errorbar(data = lat_he_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("latitude") + ylab("nuclear hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_size_manual(values = c(6, 8, 12, 14, 16), 
                    labels = c("\u22645", "\u226450", "\u2264100", "\u2264200", "\u2265200"))
msat_he_lat_plot_annotated_both <- msat_he_lat_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 70, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 60, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 50), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
msat_he_lat_plot_annotated_both

## With pelagic and coastal separate ##
#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_hd_lat_CP_plot_both <- ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey", size = 3) + 
  geom_smooth(data = CP_hd_lat_predict_df, 
              aes(x = X, y = predict_hd, color = Pelagic_Coastal), size = 6) + 
  geom_ribbon(data = CP_hd_lat_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = CP_lat_hd_binned_means, 
             aes(x = X, y = hd_mean, color = Pelagic_Coastal, size = count_bin), shape = "square") + 
  geom_errorbar(data = CP_lat_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 1.75, size = 3) + 
  xlab("latitude") + ylab("mtdna Hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP) + 
  scale_size_manual(values = c(6, 8, 12, 14, 16), 
                    labels = c("\u22645", "\u226450", "\u2264100", "\u2264200", "\u2265200"))
mtdna_hd_lat_CP_plot_annotated_both <- mtdna_hd_lat_CP_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 70, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 60, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 50),
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_lat_CP_plot_annotated_both

#############################################################################################################

######### Lon model figures #######

#### Predict ####

## Build prediction dataframe ##
lon_eff <- plot_model(lon_model_he, type = "pred", terms = "lon_rad [all]")
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
                                                     multiple = 10, FUN = floor)} #rounding DOWN (e.g. abslat 10-20 gets value of 10)
}

msat$lon_round <- as.factor(msat$lon_round)

## Calculate means and SE within each 10 degree band ##

msat <- data.table(msat) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
lon_he_mean <- msat[, mean(He), by = lon_round]
CP_lon_he_mean <- msat[, mean(He), by = list(lon_round, Pelagic_Coastal)]
  colnames(lon_he_mean) <- c("lon", "he_mean")
  colnames(CP_lon_he_mean) <- c("lon", "Pelagic_Coastal", "he_mean")

#calculate SE in each 10 degree band
lon_he_SE <- msat[, std.error(He), by = lon_round]
CP_lon_he_SE <- msat[, std.error(He), by = list(lon_round, Pelagic_Coastal)]
  colnames(lon_he_SE) <- c("lon", "he_SE")
  colnames(CP_lon_he_SE) <- c("lon", "Pelagic_Coastal", "he_SE")

#count observations in each 10 degree band
lon_he_count <- msat[, .N, by = lon_round]
CP_lon_he_count <- msat[, .N, by = list(lon_round, Pelagic_Coastal)]
  colnames(lon_he_count) <- c("lon", "he_count")
  colnames(CP_lon_he_count) <- c("lon", "Pelagic_Coastal", "he_count")

#merge mean and SE dataframes together
lon_he_binned_means <- list(lon_he_mean, lon_he_SE, lon_he_count) %>% 
  reduce(full_join, by = "lon")
  lon_he_binned_means$lon <- as.numeric(as.character(lon_he_binned_means$lon))
  lon_he_binned_means <- lon_he_binned_means[order(lon), ]
  lon_he_binned_means$X <- lon_he_binned_means$lon + 5 #for plotting, plot in MIDDLE of 10 degree band
  lon_he_binned_means$X[lon_he_binned_means$lon == 180] <- 180 #these two are on the 180 line

CP_lon_he_binned_means <- merge(CP_lon_he_mean, CP_lon_he_SE, by = c("lon", "Pelagic_Coastal"))
CP_lon_he_binned_means <- merge(CP_lon_he_binned_means, CP_lon_he_count, by = c("lon", "Pelagic_Coastal"))
  colnames(CP_lon_he_binned_means) <- c("lon", "Pelagic_Coastal", "he_mean", "he_SE", "he_count")
  CP_lon_he_binned_means$lon <- as.numeric(as.character(CP_lon_he_binned_means$lon))
  CP_lon_he_binned_means <- CP_lon_he_binned_means[order(lon), ]
  CP_lon_he_binned_means$X <- CP_lon_he_binned_means$lon + 5 #for plotting, plot in MIDDLE of 10 degree band
  CP_lon_he_binned_means$X[CP_lon_he_binned_means$lon == 180] <- 180 #these two are on the 180 line

#calculate error bars (standard error)
lon_he_binned_means$mean_lowerSE <- lon_he_binned_means$he_mean - 
  lon_he_binned_means$he_SE
lon_he_binned_means$mean_upperSE <- lon_he_binned_means$he_mean + 
  lon_he_binned_means$he_SE

CP_lon_hd_binned_means$mean_lowerSE <- CP_lon_hd_binned_means$hd_mean - 
  CP_lon_hd_binned_means$hd_SE
CP_lon_hd_binned_means$mean_upperSE <- CP_lon_hd_binned_means$hd_mean + 
  CP_lon_hd_binned_means$hd_SE

#add count bin for plotting (all < number, 5001 = >5000)
lon_he_binned_means$count_bin <- c(1000, 1000, 1000, 1000, 1000, 5000, 1000, 1000, 1000,
                                   2000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 5000, 
                                   2000, 5000, 2000, 1000, 1000, 1000, 1000, 1000, 1000, 
                                   1000, 1000, 2000, 5000, 1000, 1000, 1000, 1000, 1000)
  lon_he_binned_means$count_bin <- factor(lon_he_binned_means$count_bin, levels = c(1000, 2000, 4000, 5000))

CP_lon_hd_binned_means$count_bin <- c(50, 5, 50, 100, 5, 50, 5, 50, 5, 50, 5, 50, 5, 5, 5, 
                                      50, 5, 50, 50, 100, 50, 50, 50, 50, 50, 50, 5, 50, 5, 50,
                                      50, 50, 50, 100, 5, 100, 5, 100, 50, 50, 5, 50, 5, 50, 5, 
                                      50, 50, 50, 50, 5, 50, 50, 5, 50, 5, 100, 50, 201, 50, 100, 
                                      50, 100, 50, 50, 5, 50, 50, 50, 50, 5)
CP_lon_hd_binned_means$count_bin <- factor(CP_lon_hd_binned_means$count_bin, levels = c(5, 50, 100, 201))

#### Plot lon ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
msat_he_lon_plot_both <- ggplot() + 
  geom_line(data = lon_eff_data, 
            aes(x = lon, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = lon_eff_data, 
              aes(x = lon, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = lon_he_binned_means, 
             aes(x = X, y = he_mean, color = "10-degree binned means", size = count_bin), shape = "square") + 
  geom_errorbar(data = lon_he_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("longitude") + ylab("nuclear He") + labs(color = "Legend") + 
  scale_color_manual(values = colors) +
  scale_size_manual(values = c(8, 12, 16), 
                    labels = c("\u22641000", "\u22642000", "\u22645000"))
msat_he_lon_plot_annotated_both <- msat_he_lon_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
msat_he_lon_plot_annotated_both

## With pelagic and coastal separate ##
#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_hd_lon_CP_plot_both <- ggplot() + 
  geom_smooth(data = lon_eff_data, 
              aes(x = lon, y = predicted, color = group), size = 6) + 
  geom_ribbon(data = lon_eff_data, 
              aes(x = lon, ymin = conf.low, ymax = conf.high, color = group), alpha = 0.1) + 
  geom_point(data = CP_lon_hd_binned_means, 
             aes(x = X, y = hd_mean, color = Pelagic_Coastal, size = count_bin), shape = "square") + 
  geom_errorbar(data = CP_lon_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 1.75, size = 3) + 
  xlab("longitude") + ylab("mtdna Hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP) + 
  scale_size_manual(values = c(6, 8, 12, 14, 16), 
                    labels = c("\u22645", "\u226450", "\u2264100", "\u2264200", "\u2265200"))
mtdna_hd_lon_CP_plot_annotated_both <- mtdna_hd_lon_CP_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 70, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 60, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 50),
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_lon_CP_plot_annotated_both

#############################################################################################################

######### SSTmean model figures #######

#### Predict ####

#### Predict ####

## Build prediction dataframe ##
SSTmean_eff <- plot_model(sstmean_model_he, type = "pred", terms = "logsstmean [all]")
SSTmean_eff_data <- as.data.frame(SSTmean_eff$data)
SSTmean_eff_data$sstmean <- 10^(SSTmean_eff_data$x)
SSTmean_eff_data <- SSTmean_eff_data[-1, ]

#### Calculate means from raw data ####

# Round SSTmean DOWN to nearest multiple of 5 ##

#create column to fill with the rounded SSTmin
msat$SSTmean_round <- NA

#fill round column in
for(i in 1:nrow(msat)) {
  cat(paste(i, " ", sep = ''))
  {msat$SSTmean_round[i] <- DescTools::RoundTo(msat$sst.BO_sstmean[i], 
                                                         multiple = 5, FUN = floor)} #rounding DOWN (e.g. SSTmin 5-10 gets value of 5)
}

msat$SSTmean_round <- as.factor(msat$SSTmean_round)

## Calculate means and SE within each 10 degree band ##

msat <- data.table(msat) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
SSTmean_he_mean <- msat[, mean(He), by = SSTmean_round]
colnames(SSTmean_he_mean) <- c("SSTmean", "he_mean")

#calculate SE in each 10 degree band
SSTmean_he_SE <- msat[, std.error(He), by = SSTmean_round]
colnames(SSTmean_he_SE) <- c("SSTmean", "he_SE")

#count observations in each 10 degree band
SSTmean_he_count <- msat[, .N, by = SSTmean_round]
colnames(SSTmean_he_count) <- c("SSTmean", "he_count")

#merge mean and SE dataframes together
SSTmean_he_binned_means <- list(SSTmean_he_mean, SSTmean_he_SE, SSTmean_he_count) %>% 
  reduce(full_join, by = "SSTmean")
SSTmean_he_binned_means$SSTmean <- as.numeric(as.character(SSTmean_he_binned_means$SSTmean))
SSTmean_he_binned_means <- SSTmean_he_binned_means[order(SSTmean), ]
SSTmean_he_binned_means$X <- SSTmean_he_binned_means$SSTmean + 2.5 #for plotting, plot in MIDDLE of 5 degree band

#calculate error bars (standard error)
SSTmean_he_binned_means$mean_lowerSE <- SSTmean_he_binned_means$he_mean - 
  SSTmean_he_binned_means$he_SE
SSTmean_he_binned_means$mean_upperSE <- SSTmean_he_binned_means$he_mean + 
  SSTmean_he_binned_means$he_SE

#add count bin for plotting (all < number, 5001 = >5000)
SSTmean_he_binned_means$count_bin <- c(1000, 5001, 5001, 5001, 5000, 5001, 1000)
SSTmean_he_binned_means$count_bin <- factor(SSTmean_he_binned_means$count_bin, levels = c(1000, 5000, 5001))


#### Plot SSTmean ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
msat_he_SSTmean_plot_both <- ggplot() + 
  geom_line(data = SSTmean_eff_data, 
            aes(x = sstmean, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = SSTmean_eff_data, 
              aes(x = sstmean, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = SSTmean_he_binned_means, 
             aes(x = X, y = he_mean, color = "10-degree binned means", size = count_bin), shape = "square") + 
  geom_errorbar(data = SSTmean_he_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("SST mean") + ylab("nuclear He") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_size_manual(values = c(8, 14, 16), 
                    labels = c("\u22641000", "\u22645000", "\u22655000")) + 
  scale_y_continuous(limits = c(0.65, 0.8))
msat_he_SSTmean_plot_annotated_both <- msat_he_SSTmean_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
msat_he_SSTmean_plot_annotated_both

## With pelagic and coastal separate ##
#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_hd_sstmean_CP_plot_both <- ggplot() + 
  geom_smooth(data = CP_hd_sstmean_predict_df, 
              aes(x = X, y = predict_hd, color = Pelagic_Coastal), size = 6, se = FALSE) + 
  geom_ribbon(data = CP_hd_sstmean_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = CP_sstmean_hd_binned_means, 
             aes(x = X, y = hd_mean, color = Pelagic_Coastal, size = count_bin), shape = "square") + 
  geom_errorbar(data = CP_sstmean_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 1.75, size = 3) + 
  xlab("SST mean") + ylab("mtdna Hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP) + 
  scale_size_manual(values = c(6, 8, 12, 14, 16), 
                    labels = c("\u22645", "\u226450", "\u2264100", "\u2264200", "\u2265200"))
mtdna_hd_sstmean_CP_plot_annotated_both <- mtdna_hd_sstmean_CP_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 70, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 60, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 50),
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_sstmean_CP_plot_annotated_both
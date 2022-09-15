#Predict w/GLMMs and plot with binned means

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
library(data.table)
library(DescTools)

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)
msat <- read.csv("output/msatloci_assembled.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)
mtdna_env <- read.csv("output/mtdna_env.csv", stringsAsFactors = FALSE)
msat_env <- read.csv("output/msat_env.csv", stringsAsFactors = FALSE)

mtdna <- merge(mtdna, mtdna_env[, c('X', 'sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin', 
                                    'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                    'chloroA.BO_chlomax', 'chloroA.BO_chlomin')], all.x = TRUE)
msat <- cbind(msat[, -1], msat_env[, c('sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin', 
                                       'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                       'chloroA.BO_chlomax', 'chloroA.BO_chlomin')]) #merge not working for some reason, cbind bc in same order
mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                  'Centroid')], all.x = TRUE)
msat <- merge(msat, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                'Centroid')], all.x = TRUE)
#clean up dataframes
#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)

#calculate abslat
mtdna_small$abslat <- abs(mtdna_small$lat)
msat$abslat <- abs(msat$lat)

########## pi
#subset mtdna to remove Pi = NA columns
mtdna_small_pi <- subset(mtdna_small, mtdna_small$Pi != "NA")

#log transform and check again
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$Pi != 0) #if any zeros will screw up log transformation (log10(0) is undefined, also probably shouldn't be there anyway)
mtdna_small_pi$logpi <- log10(mtdna_small_pi$Pi)

#remove logpi = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logpi != "Inf")

######## building lat, abslat & lon model ########

#### clean up data ####
#scale geographic variables
mtdna_small_pi$lat_scale <- scale(mtdna_small_pi$lat)
mtdna_small_pi$abslat_scale <- scale(mtdna_small_pi$abslat)

#convert lon to radians
mtdna_small_pi$lon_360 <- mtdna_small_pi$lon + 180 #convert (-180,180) to (0,360)
mtdna_small_pi$lon_rad <- (2*pi*mtdna_small_pi$lon_360)/360

lat_model_pi <- lmer(logpi ~ lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                       (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

abslat_model_pi <- lmer(logpi ~ abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                       (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

lon_model_pi <- lmer(logpi ~ sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                       (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

abslat_CP_model_pi <- lmer(logpi ~ abslat_scale + I(abslat_scale^2) + Pelagic_Coastal + Pelagic_Coastal:abslat_scale + 
                                  Pelagic_Coastal:I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                  (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

lon_CP_model_pi <- lmer(logpi ~ sin(lon_rad) + cos(lon_rad) + Pelagic_Coastal + Pelagic_Coastal:sin(lon_rad) + Pelagic_Coastal:cos(lon_rad) + 
                          (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + (1|Site), 
                        REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#predict for abslat
abslat <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
#need to scale this using the same scale mean and sd as in real dataset
abslat_scale <- (abslat - attr(mtdna_small_pi$abslat_scale, "scaled:center"))/attr(mtdna_small_pi$abslat_scale, 
                                                                             "scaled:scale")

pi_abslat_predict_df <- as.data.frame(abslat_scale)
  pi_abslat_predict_df$Family <- c(rep(0, times = 9))
  pi_abslat_predict_df$Genus <- c(rep(0, times = 9))
  pi_abslat_predict_df$spp <- c(rep(0, times = 9))
  pi_abslat_predict_df$Source <- c(rep(0, times = 9))
  pi_abslat_predict_df$MarkerName <- c(rep(0, times = 9))
  pi_abslat_predict_df$Site <- c(rep(0, times = 9))
pi_abslat_predict_df$X <- abslat #for plotting

Pel_pi_abslat_predict_df <- pi_abslat_predict_df
  Pel_pi_abslat_predict_df$Pelagic_Coastal <- c(rep("Pelagic", times = 9))

Coast_pi_abslat_predict_df <- pi_abslat_predict_df
  Coast_pi_abslat_predict_df$Pelagic_Coastal <- c(rep("Coastal", times = 9))

pi_predict_abslat <- predict(abslat_model_pi, pi_abslat_predict_df, re.form = NA)
  pi_abslat_predict_df$predict_pi <- 10^(pi_predict_abslat)
  
Pel_pi_predict_abslat <- predict(abslat_CP_model_pi, Pel_pi_abslat_predict_df, re.form = NA)
  Pel_pi_abslat_predict_df$predict_pi <- 10^(Pel_pi_predict_abslat)
  
Coast_pi_predict_abslat <- predict(abslat_model_pi, Coast_pi_abslat_predict_df, re.form = NA)
  Coast_pi_abslat_predict_df$predict_pi <- 10^(Coast_pi_predict_abslat)

#other variables -- just use average?
#not interested in effect of range on regression necessarily, just account for it?

#bootstrap for confidence intervals
#want to use bootMer but on predicted dataset
  
pi_abslat_boot <- bootMer(abslat_model_pi, 
                          FUN=function(x)predict(x, pi_abslat_predict_df, re.form = NA), nsim = 100)
pi_abslat_predict_df$lower_ci <- 10^(apply(pi_abslat_boot$t, 2, quantile, 0.025))
pi_abslat_predict_df$upper_ci <- 10^(apply(pi_abslat_boot$t, 2, quantile, 0.975))

Pel_pi_abslat_boot <- bootMer(abslat_CP_model_pi, 
                          FUN=function(x)predict(x, Pel_pi_abslat_predict_df, re.form = NA), nsim = 100)
Pel_pi_abslat_predict_df$lower_ci <- 10^(apply(Pel_pi_abslat_boot$t, 2, quantile, 0.025))
Pel_pi_abslat_predict_df$upper_ci <- 10^(apply(Pel_pi_abslat_boot$t, 2, quantile, 0.975))

Coast_pi_abslat_boot <- bootMer(abslat_CP_model_pi, 
                          FUN=function(x)predict(x, Coast_pi_abslat_predict_df, re.form = NA), nsim = 100)
Coast_pi_abslat_predict_df$lower_ci <- 10^(apply(Coast_pi_abslat_boot$t, 2, quantile, 0.025))
Coast_pi_abslat_predict_df$upper_ci <- 10^(apply(Coast_pi_abslat_boot$t, 2, quantile, 0.975))

CP_pi_abslat_predict_df <- rbind(Pel_pi_abslat_predict_df, Coast_pi_abslat_predict_df)

#add real means

#mtdna pi abslat 10 binning
mtdna_small_pi$abslat_round <- NA

for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$abslat_round[i] <- DescTools::RoundTo(mtdna_small_pi$abslat[i], multiple = 10, FUN = floor)}
}

mtdna_small_pi$abslat_round <- as.factor(mtdna_small_pi$abslat_round)

#subset by lat bands and then get mean and plot that

#abslat bands
mtdna_small_pi <- data.table(mtdna_small_pi)
abslat_logpi_mean <- mtdna_small_pi[, mean(logpi), by = list(abslat_round, Pelagic_Coastal)]
  colnames(abslat_logpi_mean) <- c("abslat", "Pelagic_Coastal", "logpi_mean")
  abslat_logpi_mean$pi_mean <- 10^(abslat_logpi_mean$logpi_mean)
abslat_pi_SE <- mtdna_small_pi[, std.error(Pi), by = list(abslat_round, Pelagic_Coastal)] #not quite sure if this is correct...
  colnames(abslat_pi_SE) <- c("abslat", "Pelagic_Coastal", "pi_SE")

abslat_pi_binned_means <- merge(abslat_logpi_mean, abslat_pi_SE, by = c("abslat", "Pelagic_Coastal"))
  abslat_pi_binned_means <- abslat_pi_binned_means[order(abslat), ]
  abslat_pi_binned_means$X <- c(5, 5, 15, 15, 25, 25, 35, 35, 45, 45, 55, 55, 65, 65, 75)

#calculate error bars (standard error)
abslat_pi_binned_means$mean_lowerSE <- abslat_pi_binned_means$pi_mean - abslat_pi_binned_means$pi_SE
abslat_pi_binned_means$mean_upperSE <- abslat_pi_binned_means$pi_mean + abslat_pi_binned_means$pi_SE

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_abslat_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = abslat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = pi_abslat_predict_df, aes(x = X, y = predict_pi, color = "Regression"), size = 2) + 
  geom_ribbon(data = pi_abslat_predict_df, aes(x = X, ymin = lower_ci, ymax = upper_ci, color = "Regression"), alpha = 0.1) + 
  geom_point(data = abslat_pi_binned_means, aes(x = X, y = pi_mean, color = "10-degree binned means"), size = 8, shape = "square") + 
  geom_errorbar(data = abslat_pi_binned_means, aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = "10-degree binned means"), 
                width = 0.75, size = 1.5) + 
  xlab("absolute latitude") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors)
mtdna_pi_abslat_plot_annotated_both <- mtdna_pi_abslat_plot_both + theme_bw() + 
  scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.03)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_abslat_plot_annotated_both


#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_pi_abslat_CP_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = abslat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = CP_pi_abslat_predict_df, aes(x = X, y = predict_pi, color = Pelagic_Coastal), size = 2) + 
  geom_ribbon(data = CP_pi_abslat_predict_df, aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = abslat_pi_binned_means, aes(x = X, y = pi_mean, color = Pelagic_Coastal), size = 8, shape = "square") + 
  geom_errorbar(data = abslat_pi_binned_means, aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 0.75, size = 1.5) + 
   xlab("absolute latitude") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP)
mtdna_pi_abslat_CP_plot_annotated_both <- mtdna_pi_abslat_CP_plot_both + theme_bw() + 
  scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.03)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_abslat_CP_plot_annotated_both

#predict for lat
lat <- c(-90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80)
#need to scale this using the same scale mean and sd as in real dataset
lat_scale <- (lat - attr(mtdna_small_pi$lat_scale, "scaled:center"))/attr(mtdna_small_pi$lat_scale, 
                                                                          "scaled:scale")

pi_lat_predict_df <- as.data.frame(lat_scale)
pi_lat_predict_df$Family <- c(rep(0, times = 9))
pi_lat_predict_df$Genus <- c(rep(0, times = 9))
pi_lat_predict_df$spp <- c(rep(0, times = 9))
pi_lat_predict_df$Source <- c(rep(0, times = 9))
pi_lat_predict_df$MarkerName <- c(rep(0, times = 9))
pi_lat_predict_df$Site <- c(rep(0, times = 9))
pi_lat_predict_df$X <- lat #for plotting

pi_predict_lat <- predict(lat_model_pi, pi_lat_predict_df, re.form = NA)
pi_lat_predict_df$predict_pi <- 10^(pi_predict_lat)

#other variables -- just use average?
#not interested in effect of range on regression necessarily, just account for it?

#bootstrap for confidence intervals
#want to use bootMer but on predicted dataset

pi_lat_boot <- bootMer(lat_model_pi, 
                       FUN=function(x)predict(x, pi_lat_predict_df, re.form = NA), nsim = 100)
pi_lat_predict_df$lower_ci <- 10^(apply(pi_lat_boot$t, 2, quantile, 0.025))
pi_lat_predict_df$upper_ci <- 10^(apply(pi_lat_boot$t, 2, quantile, 0.975))

#add real means
#mtdna pi lat 10 binning
mtdna_small_pi$lat_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)

for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$lat_round[i] <- DescTools::RoundTo(mtdna_small_pi$lat[i], multiple = 10, FUN = floor)}
}

mtdna_small_pi$lat_round <- as.factor(mtdna_small_pi$lat_round)

#subset by lat bands and then get mean and plot that

#lat bands
mtdna_small_pi <- data.table(mtdna_small_pi)
lat_logpi_mean <- mtdna_small_pi[, mean(logpi), by = list(lat_round)]
  colnames(lat_logpi_mean) <- c("lat", "logpi_mean")
  lat_logpi_mean$pi_mean <- 10^(lat_logpi_mean$logpi_mean)
lat_pi_SE <- mtdna_small_pi[, std.error(Pi), by = list(lat_round)] #not quite sure if this is correct...
  colnames(lat_pi_SE) <- c("lat", "pi_SE")

lat_pi_binned_means <- merge(lat_logpi_mean, lat_pi_SE, by = "lat")
  lat_pi_binned_means <- lat_pi_binned_means[order(lat), ]
  lat_pi_binned_means$X <- c(-75, -65, -55, -45, -35, -25, -15, -5, 
                             5, 15, 25, 35, 45, 55, 65, 75)

#Calculate error bars (standard error)
lat_pi_binned_means$mean_lowerSE <- lat_pi_binned_means$pi_mean - lat_pi_binned_means$pi_SE
lat_pi_binned_means$mean_upperSE <- lat_pi_binned_means$pi_mean + lat_pi_binned_means$pi_SE

#for legend
#colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_lat_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = lat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = pi_lat_predict_df, aes(x = X, y = predict_pi, color = "Regression"), size = 2) + 
  geom_ribbon(data = pi_lat_predict_df, aes(x = X, ymin = lower_ci, ymax = upper_ci), alpha = 0.1) + 
  geom_point(data = lat_pi_binned_means, aes(x = X, y = pi_mean, color = "10-degree binned means"), size = 8, shape = "square") + 
  geom_errorbar(data = lat_pi_binned_means, aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = "10-degree binned means"), 
                width = 0.75, size = 1.5) + 
  xlab("latitude") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors)
mtdna_pi_lat_plot_annotated_both <- mtdna_pi_lat_plot_both + theme_bw() + 
  #scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.02)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_lat_plot_annotated_both

#predict for lon
lon <- c(-180, -170, -160, -150, -140, -130, -120, -110, -100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 
         0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170)
  lon_360 <- (lon + 180)
  lon_rad <- (2*pi*lon)/360

pi_lon_predict_df <- as.data.frame(lon_rad)
pi_lon_predict_df$Family <- c(rep(0, times = 9))
pi_lon_predict_df$Genus <- c(rep(0, times = 9))
pi_lon_predict_df$spp <- c(rep(0, times = 9))
pi_lon_predict_df$Source <- c(rep(0, times = 9))
pi_lon_predict_df$MarkerName <- c(rep(0, times = 9))
pi_lon_predict_df$Site <- c(rep(0, times = 9))
pi_lon_predict_df$X <- lon #for plotting

Pel_pi_lon_predict_df <- pi_lon_predict_df
  Pel_pi_lon_predict_df$Pelagic_Coastal <- c(rep("Pelagic", times = 9))

Coast_pi_lon_predict_df <- pi_lon_predict_df
  Coast_pi_lon_predict_df$Pelagic_Coastal <- c(rep("Coastal", times = 9))

pi_predict_lon <- predict(lon_model_pi, pi_lon_predict_df, re.form = NA)
  pi_lon_predict_df$predict_pi <- 10^(pi_predict_lon)

Pel_pi_predict_lon <- predict(lon_CP_model_pi, Pel_pi_lon_predict_df, re.form = NA)
  Pel_pi_lon_predict_df$predict_pi <- 10^(Pel_pi_predict_lon)
  
Coast_pi_predict_lon <- predict(lon_CP_model_pi, Coast_pi_lon_predict_df, re.form = NA)
  Coast_pi_lon_predict_df$predict_pi <- 10^(Coast_pi_predict_lon)

#other variables -- just use average?
#not interested in effect of range on regression necessarily, just account for it?

#bootstrap for confidence intervals
#want to use bootMer but on predicted dataset

pi_lon_boot <- bootMer(lon_model_pi, 
                       FUN=function(x)predict(x, pi_lon_predict_df, re.form = NA), nsim = 100)
  pi_lon_predict_df$lower_ci <- 10^(apply(pi_lon_boot$t, 2, quantile, 0.025))
  pi_lon_predict_df$upper_ci <- 10^(apply(pi_lon_boot$t, 2, quantile, 0.975))

Pel_pi_lon_boot <- bootMer(lon_CP_model_pi, 
                       FUN=function(x)predict(x, Pel_pi_lon_predict_df, re.form = NA), nsim = 100)
  Pel_pi_lon_predict_df$lower_ci <- 10^(apply(Pel_pi_lon_boot$t, 2, quantile, 0.025))
  Pel_pi_lon_predict_df$upper_ci <- 10^(apply(Pel_pi_lon_boot$t, 2, quantile, 0.975))

Coast_pi_lon_boot <- bootMer(lon_CP_model_pi, 
                       FUN=function(x)predict(x, Coast_pi_lon_predict_df, re.form = NA), nsim = 100)
  Coast_pi_lon_predict_df$lower_ci <- 10^(apply(Coast_pi_lon_boot$t, 2, quantile, 0.025))
  Coast_pi_lon_predict_df$upper_ci <- 10^(apply(Coast_pi_lon_boot$t, 2, quantile, 0.975))
  
CP_pi_lon_predict_df <- rbind(Pel_pi_lon_predict_df, Coast_pi_lon_predict_df)

#add real means
#mtdna pi lon 10 binning
mtdna_small_pi$lon_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)

for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$lon_round[i] <- DescTools::RoundTo(mtdna_small_pi$lon[i], multiple = 10, FUN = floor)}
}

mtdna_small_pi$lon_round <- as.factor(mtdna_small_pi$lon_round)

#subset by lon bands and then get mean and plot that

#lon bands
mtdna_small_pi <- data.table(mtdna_small_pi)
lon_logpi_mean <- mtdna_small_pi[, mean(logpi), by = list(lon_round, Pelagic_Coastal)]
  colnames(lon_logpi_mean) <- c("lon", "Pelagic_Coastal", "logpi_mean")
  lon_logpi_mean$pi_mean <- 10^(lon_logpi_mean$logpi_mean)
lon_pi_SE <- mtdna_small_pi[, std.error(Pi), by = list(lon_round, Pelagic_Coastal)] #not quite sure if this is correct...
colnames(lon_pi_SE) <- c("lon", "Pelagic_Coastal", "pi_SE")

lon_pi_binned_means <- merge(lon_logpi_mean, lon_pi_SE, by = c("lon", "Pelagic_Coastal"))
  lon_pi_binned_means <- lon_pi_binned_means[order(lon), ]
  lon_pi_binned_means$X <- c(-175, -175, -165, -155, -155, -145, -145, -135, -135, -125, -125, -115, -115, -105, -105, 
                             -95, -95, -85, -85, -75, -75, -65, -65, -55, -55, -45, -45, -35, -35, -25, -25, -15, -15, -5, -5, 
                             5, 5, 15, 15, 25, 25, 35, 35, 45, 45, 55, 55, 75, 75, 85, 85, 95, 95, 105, 105, 115, 115, 
                             125, 125, 135, 135, 145, 145, 155, 155, 165, 165, 175, 175, 180) #180 is bc range goes to 180, so must be observation(s) at 180

#calculate error bars (standard error)
lon_pi_binned_means$mean_lowerSE <- lon_pi_binned_means$pi_mean - lon_pi_binned_means$pi_SE
lon_pi_binned_means$mean_upperSE <- lon_pi_binned_means$pi_mean + lon_pi_binned_means$pi_SE

#for legend
#colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_lon_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = lon, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = pi_lon_predict_df, aes(x = X, y = predict_pi, color = "Regression"), size = 2) + 
  geom_ribbon(data = pi_lon_predict_df, aes(x = X, ymin = lower_ci, ymax = upper_ci), alpha = 0.1) + 
  geom_point(data = lon_pi_binned_means, aes(x = X, y = pi_mean, color = "10-degree binned means"), size = 8, shape = "square") + 
  geom_errorbar(data = lon_pi_binned_means, aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = "10-degree binned means"), 
                width = 0.75, size = 1.5) + 
  xlab("longitude") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors)
mtdna_pi_lon_plot_annotated_both <- mtdna_pi_lon_plot_both + theme_bw() + 
  #scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.05)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_lon_plot_annotated_both

#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_pi_lon_CP_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = abslat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = CP_pi_lon_predict_df, aes(x = X, y = predict_pi, color = Pelagic_Coastal), size = 2) + 
  geom_ribbon(data = CP_pi_lon_predict_df, aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = lon_pi_binned_means, aes(x = X, y = pi_mean, color = Pelagic_Coastal), size = 8, shape = "square") + 
  geom_errorbar(data = lon_pi_binned_means, aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 0.75, size = 1.5) + 
  xlab("longitude") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP)
mtdna_pi_lon_CP_plot_annotated_both <- mtdna_pi_lon_CP_plot_both + theme_bw() + 
  #scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.03)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_lon_CP_plot_annotated_both

#CI --> Zoe did it this way
#used Ben Bolker's tutorial https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions

#predict for SST
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$sst.BO_sstmean != "NA") #should look at where these are...

#### log transform sst data ####
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$sst.BO_sstmin != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi$logsstmin <- log10(mtdna_small_pi$sst.BO_sstmin)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstmin != "NaN")

##### sst min model ####
SSTmin_model_pi <- lmer(logpi ~ logsstmin + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                               (1|Site), REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

SSTmin_CP_model_pi <- lmer(logpi ~ logsstmin + Pelagic_Coastal + Pelagic_Coastal:logsstmin + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                          (1|Site), REML = FALSE, data = mtdna_small_pi, 
                          na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?
  
#predict for SST min (range is 0.129 to 30.295 C)
SSTmin <- c(.01, 5, 10, 15, 20, 25, 30, 35) #starting at .01 bc log10 of 0 is Inf
#need to logtransform
logsstmin <- log10(SSTmin)

pi_SSTmin_predict_df <- as.data.frame(logsstmin)
pi_SSTmin_predict_df$Family <- c(rep(0, times = 8))
pi_SSTmin_predict_df$Genus <- c(rep(0, times = 8))
pi_SSTmin_predict_df$spp <- c(rep(0, times = 8))
pi_SSTmin_predict_df$Source <- c(rep(0, times = 8))
pi_SSTmin_predict_df$MarkerName <- c(rep(0, times = 8))
pi_SSTmin_predict_df$Site <- c(rep(0, times = 8))
pi_SSTmin_predict_df$X <- SSTmin #need this for plotting

Pel_pi_SSTmin_predict_df <- pi_SSTmin_predict_df
  Pel_pi_SSTmin_predict_df$Pelagic_Coastal <- c(rep("Pelagic", times = 8))
  
Coast_pi_SSTmin_predict_df <- pi_SSTmin_predict_df
  Coast_pi_SSTmin_predict_df$Pelagic_Coastal <- c(rep("Coastal", times = 8))  

pi_predict_SSTmin <- predict(SSTmin_model_pi, pi_SSTmin_predict_df, re.form = NA)
pi_SSTmin_predict_df$predict_pi <- 10^(pi_predict_SSTmin)

Pel_pi_predict_SSTmin <- predict(SSTmin_CP_model_pi, Pel_pi_SSTmin_predict_df, re.form = NA)
  Pel_pi_SSTmin_predict_df$predict_pi <- 10^(Pel_pi_predict_SSTmin)

Coast_pi_predict_SSTmin <- predict(SSTmin_CP_model_pi, Coast_pi_SSTmin_predict_df, re.form = NA)
  Coast_pi_SSTmin_predict_df$predict_pi <- 10^(Coast_pi_predict_SSTmin)

#bootstrap for confidence intervals
#want to use bootMer but on predicted dataset

pi_SSTmin_boot <- bootMer(SSTmin_model_pi, 
                       FUN=function(x)predict(x, pi_SSTmin_predict_df, re.form = NA), nsim = 100)
  pi_SSTmin_predict_df$lower_ci <- 10^(apply(pi_SSTmin_boot$t, 2, quantile, 0.025))
  pi_SSTmin_predict_df$upper_ci <- 10^(apply(pi_SSTmin_boot$t, 2, quantile, 0.975))

Pel_pi_SSTmin_boot <- bootMer(SSTmin_CP_model_pi, 
                          FUN=function(x)predict(x, Pel_pi_SSTmin_predict_df, re.form = NA), nsim = 100)
  Pel_pi_SSTmin_predict_df$lower_ci <- 10^(apply(Pel_pi_SSTmin_boot$t, 2, quantile, 0.025))
  Pel_pi_SSTmin_predict_df$upper_ci <- 10^(apply(Pel_pi_SSTmin_boot$t, 2, quantile, 0.975))

Coast_pi_SSTmin_boot <- bootMer(SSTmin_CP_model_pi, 
                          FUN=function(x)predict(x, Coast_pi_SSTmin_predict_df, re.form = NA), nsim = 100)
  Coast_pi_SSTmin_predict_df$lower_ci <- 10^(apply(Coast_pi_SSTmin_boot$t, 2, quantile, 0.025))
  Coast_pi_SSTmin_predict_df$upper_ci <- 10^(apply(Coast_pi_SSTmin_boot$t, 2, quantile, 0.975))

CP_pi_SSTmin_predict_df <- rbind(Pel_pi_SSTmin_predict_df, Coast_pi_SSTmin_predict_df)
  
#add real means
#mtdna pi SSTmin binning
mtdna_small_pi$SSTmin_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)

for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$SSTmin_round[i] <- DescTools::RoundTo(mtdna_small_pi$sst.BO_sstmin[i], multiple = 5, FUN = floor)}
}

mtdna_small_pi$SSTmin_round <- as.factor(mtdna_small_pi$SSTmin_round)

#SSTmin bands
mtdna_small_pi <- data.table(mtdna_small_pi)
SSTmin_logpi_mean <- mtdna_small_pi[, mean(logpi), by = list(SSTmin_round, Pelagic_Coastal)]
  colnames(SSTmin_logpi_mean) <- c("SSTmin", "Pelagic_Coastal", "logpi_mean")
  SSTmin_logpi_mean$pi_mean <- 10^(SSTmin_logpi_mean$logpi_mean)
SSTmin_pi_SE <- mtdna_small_pi[, std.error(Pi), by = list(SSTmin_round, Pelagic_Coastal)] #not quite sure if this is correct...
  colnames(SSTmin_pi_SE) <- c("SSTmin", "Pelagic_Coastal", "pi_SE")

SSTmin_pi_binned_means <- merge(SSTmin_logpi_mean, SSTmin_pi_SE, by = c("SSTmin", "Pelagic_Coastal"))
SSTmin_pi_binned_means <- SSTmin_pi_binned_means[order(SSTmin), ]
  SSTmin_pi_binned_means$X <- c(2.5, 2.5, 7.5, 7.5, 12.5, 12.5, 17.5, 17.5, 22.5, 22.5, 27.5, 27.5, 32.5) #because took mean every 5 units, so should plot in the middle of range, not the floor

#add error bars (standard error)
SSTmin_pi_binned_means$mean_lowerSE <- SSTmin_pi_binned_means$pi_mean - SSTmin_pi_binned_means$pi_SE
SSTmin_pi_binned_means$mean_upperSE <- SSTmin_pi_binned_means$pi_mean + SSTmin_pi_binned_means$pi_SE

#for legend
#colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_SSTmin_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = lat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = pi_SSTmin_predict_df, aes(x = X, y = predict_pi, color = "Regression"), size = 2) + 
  geom_ribbon(data = pi_SSTmin_predict_df, aes(x = X, ymin = lower_ci, ymax = upper_ci), alpha = 0.1) + 
  geom_point(data = SSTmin_pi_binned_means, aes(x = X, y = pi_mean, color = "10-degree binned means"), size = 8, shape = "square") + 
  geom_errorbar(data = SSTmin_pi_binned_means, aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = "10-degree binned means"), 
                width = 0.75, size = 1.5) + 
  xlab("SSTmin") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors)
mtdna_pi_SSTmin_plot_annotated_both <- mtdna_pi_SSTmin_plot_both + theme_bw() + 
  #scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.02)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_SSTmin_plot_annotated_both

#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_pi_SSTmin_CP_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = abslat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = CP_pi_SSTmin_predict_df, aes(x = X, y = predict_pi, color = Pelagic_Coastal), size = 2) + 
  geom_ribbon(data = CP_pi_SSTmin_predict_df, aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = SSTmin_pi_binned_means, aes(x = X, y = pi_mean, color = Pelagic_Coastal), size = 8, shape = "square") + 
  geom_errorbar(data = SSTmin_pi_binned_means, aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 0.75, size = 1.5) + 
  xlab("SSTmin") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP)
mtdna_pi_SSTmin_CP_plot_annotated_both <- mtdna_pi_SSTmin_CP_plot_both + theme_bw() + 
  #scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.03)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_SSTmin_CP_plot_annotated_both

#predict for dissox
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$BO_dissox != "NA") #should look at where these are...

#### log transform dissox data ####
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$BO_dissox != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
mtdna_small_pi$logdissox <- log10(mtdna_small_pi$BO_dissox)
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logdissox != "NaN")

##### dissox mean model ####
dissox_model_pi <- lmer(logpi ~ logdissox + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), REML = FALSE, data = mtdna_small_pi, 
                               na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

dissox_CP_model_pi <- lmer(logpi ~ logdissox + Pelagic_Coastal + Pelagic_Coastal:logdissox + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                          (1|Site), REML = FALSE, data = mtdna_small_pi, 
                        na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?


#predict for dissox (range is 4.213 to 8.426)
dissox <- c(4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5)
#need to logtransform
logdissox <- log10(dissox)

pi_dissox_predict_df <- as.data.frame(logdissox)
pi_dissox_predict_df$Family <- c(rep(0, times = 10))
pi_dissox_predict_df$Genus <- c(rep(0, times = 10))
pi_dissox_predict_df$spp <- c(rep(0, times = 10))
pi_dissox_predict_df$Source <- c(rep(0, times = 10))
pi_dissox_predict_df$MarkerName <- c(rep(0, times = 10))
pi_dissox_predict_df$Site <- c(rep(0, times = 10))
pi_dissox_predict_df$X <- dissox #need this for plotting

Pel_pi_dissox_predict_df <- pi_dissox_predict_df
  Pel_pi_dissox_predict_df$Pelagic_Coastal <- c(rep("Pelagic", times = 10))
  
Coast_pi_dissox_predict_df <- pi_dissox_predict_df
  Coast_pi_dissox_predict_df$Pelagic_Coastal <- c(rep("Coastal", times = 10))

pi_predict_dissox <- predict(dissox_model_pi, pi_dissox_predict_df, re.form = NA)
pi_dissox_predict_df$predict_pi <- 10^(pi_predict_dissox)

Pel_pi_predict_dissox <- predict(dissox_CP_model_pi, Pel_pi_dissox_predict_df, re.form = NA)
  Pel_pi_dissox_predict_df$predict_pi <- 10^(Pel_pi_predict_dissox)

Coast_pi_predict_dissox <- predict(dissox_CP_model_pi, Coast_pi_dissox_predict_df, re.form = NA)
  Coast_pi_dissox_predict_df$predict_pi <- 10^(Coast_pi_predict_dissox)

#bootstrap for confidence intervals
#want to use bootMer but on predicted dataset

pi_dissox_boot <- bootMer(dissox_model_pi, 
                          FUN=function(x)predict(x, pi_dissox_predict_df, re.form = NA), nsim = 100)
  pi_dissox_predict_df$lower_ci <- 10^(apply(pi_dissox_boot$t, 2, quantile, 0.025))
  pi_dissox_predict_df$upper_ci <- 10^(apply(pi_dissox_boot$t, 2, quantile, 0.975))
  
Pel_pi_dissox_boot <- bootMer(dissox_CP_model_pi, 
                          FUN=function(x)predict(x, Pel_pi_dissox_predict_df, re.form = NA), nsim = 100)
  Pel_pi_dissox_predict_df$lower_ci <- 10^(apply(Pel_pi_dissox_boot$t, 2, quantile, 0.025))
  Pel_pi_dissox_predict_df$upper_ci <- 10^(apply(Pel_pi_dissox_boot$t, 2, quantile, 0.975))
  
Coast_pi_dissox_boot <- bootMer(dissox_CP_model_pi, 
                           FUN=function(x)predict(x, Coast_pi_dissox_predict_df, re.form = NA), nsim = 100)
  Coast_pi_dissox_predict_df$lower_ci <- 10^(apply(Coast_pi_dissox_boot$t, 2, quantile, 0.025))
  Coast_pi_dissox_predict_df$upper_ci <- 10^(apply(Coast_pi_dissox_boot$t, 2, quantile, 0.975))

CP_pi_dissox_predict_df <- rbind(Pel_pi_dissox_predict_df, Coast_pi_dissox_predict_df)
  
#add real means
#mtdna pi dissox binning
mtdna_small_pi$dissox_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)

for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$dissox_round[i] <- DescTools::RoundTo(mtdna_small_pi$BO_dissox[i], multiple = 0.5, FUN = floor)}
}

mtdna_small_pi$dissox_round <- as.factor(mtdna_small_pi$dissox_round)

#dissox bands
mtdna_small_pi <- data.table(mtdna_small_pi)
dissox_logpi_mean <- mtdna_small_pi[, mean(logpi), by = list(dissox_round, Pelagic_Coastal)]
colnames(dissox_logpi_mean) <- c("dissox", "Pelagic_Coastal", "logpi_mean")
dissox_logpi_mean$pi_mean <- 10^(dissox_logpi_mean$logpi_mean)
dissox_pi_SE <- mtdna_small_pi[, std.error(Pi), by = list(dissox_round, Pelagic_Coastal)] #not quite sure if this is correct...
colnames(dissox_pi_SE) <- c("dissox", "Pelagic_Coastal", "pi_SE")

dissox_pi_binned_means <- merge(dissox_logpi_mean, dissox_pi_SE, by = c("dissox", "Pelagic_Coastal"))
dissox_pi_binned_means <- dissox_pi_binned_means[order(dissox), ]
dissox_pi_binned_means$X <- c(4.25, 4.25, 4.75, 4.75, 5.25, 5.25, 5.75, 5.75, 6.25, 6.25, 6.75, 6.75, 
                              7.25, 7.25, 7.75, 7.75, 8.25, 8.25) #because took mean every 5 units, so should plot in the middle of range, not the floor

#add error bars (standard error)
dissox_pi_binned_means$mean_lowerSE <- dissox_pi_binned_means$pi_mean - dissox_pi_binned_means$pi_SE
dissox_pi_binned_means$mean_upperSE <- dissox_pi_binned_means$pi_mean + dissox_pi_binned_means$pi_SE

#for legend
#colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_dissox_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = lat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = pi_dissox_predict_df, aes(x = X, y = predict_pi, color = "Regression"), size = 2) + 
  geom_ribbon(data = pi_dissox_predict_df, aes(x = X, ymin = lower_ci, ymax = upper_ci), alpha = 0.1) + 
  geom_point(data = dissox_pi_binned_means, aes(x = X, y = pi_mean, color = "10-degree binned means"), size = 8, shape = "square") + 
  geom_errorbar(data = dissox_pi_binned_means, aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = "10-degree binned means"), 
                width = 0.75, size = 1.5) + 
  xlab("Dissolved oxygen (mean)") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors)
mtdna_pi_dissox_plot_annotated_both <- mtdna_pi_dissox_plot_both + theme_bw() + 
  #scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.02)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_dissox_plot_annotated_both

#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_pi_dissox_CP_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = abslat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = CP_pi_dissox_predict_df, aes(x = X, y = predict_pi, color = Pelagic_Coastal), size = 2) + 
  geom_ribbon(data = CP_pi_dissox_predict_df, aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = dissox_pi_binned_means, aes(x = X, y = pi_mean, color = Pelagic_Coastal), size = 8, shape = "square") + 
  geom_errorbar(data = dissox_pi_binned_means, aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 0.05, size = 1.5) + 
  xlab("Dissolved oxygen (mean)") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP)
mtdna_pi_dissox_CP_plot_annotated_both <- mtdna_pi_dissox_CP_plot_both + theme_bw() + 
  #scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.03)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_dissox_CP_plot_annotated_both

#predict for chloroA
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomax != "NA") #should look at where these are...

#### log transform chloroA data ####
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomax != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
mtdna_small_pi$logchlomax <- log10(mtdna_small_pi$chloroA.BO_chlomax)
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlomax != "NaN")

##### chlomax mean model ####
chloroAmax_model_pi <- lmer(logpi ~ logchlomax + I(logchlomax^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                          (1|Site), REML = FALSE, data = mtdna_small_pi, 
                                        na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

chloroAmax_CP_model_pi <- lmer(logpi ~ logchlomax + I(logchlomax^2) + Pelagic_Coastal + Pelagic_Coastal:logchlomax + Pelagic_Coastal:I(logchlomax^2) 
                               + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                              (1|Site), REML = FALSE, data = mtdna_small_pi, 
                            na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?


#predict for chlomax (range is 0.037 to 61.783)
chlomax <- c(0.01, 10, 20, 30, 40, 50, 60) #starting at .01 bc log10 of 0 is Inf
#need to logtransform
logchlomax <- log10(chlomax)

pi_chlomax_predict_df <- as.data.frame(logchlomax)
pi_chlomax_predict_df$Family <- c(rep(0, times = 7))
pi_chlomax_predict_df$Genus <- c(rep(0, times = 7))
pi_chlomax_predict_df$spp <- c(rep(0, times = 7))
pi_chlomax_predict_df$Source <- c(rep(0, times = 7))
pi_chlomax_predict_df$MarkerName <- c(rep(0, times = 7))
pi_chlomax_predict_df$Site <- c(rep(0, times = 7))
pi_chlomax_predict_df$X <- chlomax #need this for plotting

Pel_pi_chlomax_predict_df <- pi_chlomax_predict_df
  Pel_pi_chlomax_predict_df$Pelagic_Coastal <- c(rep("Pelagic", times = 7))
  
Coast_pi_chlomax_predict_df <- pi_chlomax_predict_df
  Coast_pi_chlomax_predict_df$Pelagic_Coastal <- c(rep("Coastal", times = 7))

pi_predict_chlomax <- predict(chloroAmax_model_pi, pi_chlomax_predict_df, re.form = NA)
  pi_chlomax_predict_df$predict_pi <- 10^(pi_predict_chlomax)

Pel_pi_predict_chlomax <- predict(chloroAmax_CP_model_pi, Pel_pi_chlomax_predict_df, re.form = NA)
  Pel_pi_chlomax_predict_df$predict_pi <- 10^(Pel_pi_predict_chlomax)

Coast_pi_predict_chlomax <- predict(chloroAmax_CP_model_pi, Coast_pi_chlomax_predict_df, re.form = NA)
  Coast_pi_chlomax_predict_df$predict_pi <- 10^(Coast_pi_predict_chlomax)

#bootstrap for confidence intervals
#want to use bootMer but on predicted dataset

pi_chlomax_boot <- bootMer(chloroAmax_model_pi, 
                          FUN=function(x)predict(x, pi_chlomax_predict_df, re.form = NA), nsim = 100)
pi_chlomax_predict_df$lower_ci <- 10^(apply(pi_chlomax_boot$t, 2, quantile, 0.025))
pi_chlomax_predict_df$upper_ci <- 10^(apply(pi_chlomax_boot$t, 2, quantile, 0.975))

Pel_pi_chlomax_boot <- bootMer(chloroAmax_CP_model_pi, 
                           FUN=function(x)predict(x, Pel_pi_chlomax_predict_df, re.form = NA), nsim = 100)
  Pel_pi_chlomax_predict_df$lower_ci <- 10^(apply(Pel_pi_chlomax_boot$t, 2, quantile, 0.025))
  Pel_pi_chlomax_predict_df$upper_ci <- 10^(apply(Pel_pi_chlomax_boot$t, 2, quantile, 0.975))

Coast_pi_chlomax_boot <- bootMer(chloroAmax_CP_model_pi, 
                           FUN=function(x)predict(x, Coast_pi_chlomax_predict_df, re.form = NA), nsim = 100)
  Coast_pi_chlomax_predict_df$lower_ci <- 10^(apply(Coast_pi_chlomax_boot$t, 2, quantile, 0.025))
  Coast_pi_chlomax_predict_df$upper_ci <- 10^(apply(Coast_pi_chlomax_boot$t, 2, quantile, 0.975))

CP_pi_chlomax_predict_df <- rbind(Pel_pi_chlomax_predict_df, Coast_pi_chlomax_predict_df)
  
#add real means
#mtdna pi chlomax binning
mtdna_small_pi$chlomax_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)

for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$chlomax_round[i] <- DescTools::RoundTo(mtdna_small_pi$chloroA.BO_chlomax[i], multiple = 10, FUN = floor)}
}

mtdna_small_pi$chlomax_round <- as.factor(mtdna_small_pi$chlomax_round)

#chlomax bands
mtdna_small_pi <- data.table(mtdna_small_pi)
chlomax_logpi_mean <- mtdna_small_pi[, mean(logpi), by = list(chlomax_round, Pelagic_Coastal)]
colnames(chlomax_logpi_mean) <- c("chlomax", "Pelagic_Coastal", "logpi_mean")
chlomax_logpi_mean$pi_mean <- 10^(chlomax_logpi_mean$logpi_mean)
chlomax_pi_SE <- mtdna_small_pi[, std.error(Pi), by = list(chlomax_round, Pelagic_Coastal)] #not quite sure if this is correct...
colnames(chlomax_pi_SE) <- c("chlomax", "Pelagic_Coastal", "pi_SE")

chlomax_pi_binned_means <- merge(chlomax_logpi_mean, chlomax_pi_SE, by = c("chlomax", "Pelagic_Coastal"))
chlomax_pi_binned_means <- chlomax_pi_binned_means[order(chlomax), ]
chlomax_pi_binned_means$X <- c(5, 5, 15, 15, 25, 25, 35, 35, 45, 55, 60.5) #because took mean every 5 units, so should plot in the middle of range, not the floor --> except 60 as range of values only go to 61

#add error bars (standard error)
chlomax_pi_binned_means$mean_lowerSE <- chlomax_pi_binned_means$pi_mean - chlomax_pi_binned_means$pi_SE
chlomax_pi_binned_means$mean_upperSE <- chlomax_pi_binned_means$pi_mean + chlomax_pi_binned_means$pi_SE

#for legend
#colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_chlomax_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = lat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = pi_chlomax_predict_df, aes(x = X, y = predict_pi, color = "Regression"), size = 2) + 
  geom_ribbon(data = pi_chlomax_predict_df, aes(x = X, ymin = lower_ci, ymax = upper_ci), alpha = 0.1) + 
  geom_point(data = chlomax_pi_binned_means, aes(x = X, y = pi_mean, color = "10-degree binned means"), size = 8, shape = "square") + 
  geom_errorbar(data = chlomax_pi_binned_means, aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = "10-degree binned means"), 
                width = 0.75, size = 1.5) + 
  xlab("Chlorophyll A (max)") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors)
mtdna_pi_chlomax_plot_annotated_both <- mtdna_pi_chlomax_plot_both + theme_bw() + 
  #scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.02)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_chlomax_plot_annotated_both

#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_pi_chlomax_CP_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = abslat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = CP_pi_chlomax_predict_df, aes(x = X, y = predict_pi, color = Pelagic_Coastal), size = 2) + 
  geom_ribbon(data = CP_pi_chlomax_predict_df, aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = chlomax_pi_binned_means, aes(x = X, y = pi_mean, color = Pelagic_Coastal), size = 8, shape = "square") + 
  geom_errorbar(data = chlomax_pi_binned_means, aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 0.05, size = 1.5) + 
  xlab("Chlorophyll A (max)") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP)
mtdna_pi_chlomax_CP_plot_annotated_both <- mtdna_pi_chlomax_CP_plot_both + theme_bw() + 
  #scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.03)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_chlomax_CP_plot_annotated_both
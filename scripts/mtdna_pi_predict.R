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

###################################################################################################

######### Clean up dataframe ########

#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)

#calculate abslat
mtdna_small$abslat <- abs(mtdna_small$lat)

#subset mtdna to remove Pi = NA columns
mtdna_small_pi <- subset(mtdna_small, mtdna_small$Pi != "NA")

#log transform and check again
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$Pi != 0) #if any zeros will screw up log transformation (log10(0) is undefined, also probably shouldn't be there anyway)
  mtdna_small_pi$logpi <- log10(mtdna_small_pi$Pi)

#remove logpi = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logpi != "Inf")

###################################################################################################

######## Build lat, abslat & lon models ########

#### clean up data ####

#scale geographic variables
mtdna_small_pi$lat_scale <- scale(mtdna_small_pi$lat)
  mtdna_small_pi$lat_scale <- as.numeric(mtdna_small_pi$lat_scale)
mtdna_small_pi$abslat_scale <- scale(mtdna_small_pi$abslat)
  mtdna_small_pi$abslat_scale <- as.numeric(mtdna_small_pi$abslat_scale)

#convert lon to radians
mtdna_small_pi$lon_360 <- mtdna_small_pi$lon + 180 #convert (-180,180) to (0,360)
mtdna_small_pi$lon_rad <- (2*pi*mtdna_small_pi$lon_360)/360

#### Models ####

lat_model_pi <- lm(logpi ~ lat_scale + I(lat_scale^2), 
                     data = mtdna_small_pi, na.action = "na.fail")

abslat_model_pi <- lm(logpi ~ abslat_scale + I(abslat_scale^2), 
                        data = mtdna_small_pi, na.action = "na.fail")

#abslat_CP_model_pi <- lmer(logpi ~ abslat_scale + I(abslat_scale^2) + Pelagic_Coastal + 
 #                            Pelagic_Coastal:abslat_scale + Pelagic_Coastal:I(abslat_scale^2) + 
  #                           (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + (1|Site), 
   #                        REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", 
    #                       control = lmerControl(optimizer = "bobyqa"))

lon_model_pi <- lm(logpi ~ sin(lon_rad) + cos(lon_rad), data = mtdna_small_pi, 
                     na.action = "na.fail")

#lon_CP_model_pi <- lmer(logpi ~ sin(lon_rad) + cos(lon_rad) + Pelagic_Coastal + 
 #                         Pelagic_Coastal:sin(lon_rad) + Pelagic_Coastal:cos(lon_rad) + 
  #                        (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + (1|Site), 
   #                     REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", 
    #                    control = lmerControl(optimizer = "bobyqa"))

#############################################################################################################

######### Abslat figures #######

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

#add count bin for plotting (all < number, 201 = >200)
abslat_pi_binned_means$count_bin <- c(201, 201, 201, 201, 201, 100, 50, 50)
  abslat_pi_binned_means$count_bin <- factor(abslat_pi_binned_means$count_bin, levels = c(50, 100, 201))

#### Plot abslat ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_abslat_plot_both <- ggplot() + 
  geom_line(data = abslat_eff_data, 
            aes(x = abslat, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = abslat_eff_data, 
              aes(x = abslat, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = abslat_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "10-degree binned means", size = count_bin), shape = "square") + 
  geom_errorbar(data = abslat_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("absolute latitude") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 0.0075)) + 
  scale_size_manual(values = c(8, 12, 16), 
                    labels = c("\u226450", "\u2264100", "\u2265200"))
mtdna_pi_abslat_plot_annotated_both <- mtdna_pi_abslat_plot_both + coord_flip() + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "none", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_abslat_plot_annotated_both

## With pelagic and coastal separate ##
#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_pi_abslat_CP_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = abslat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = CP_pi_abslat_predict_df, 
            aes(x = X, y = predict_pi, color = Pelagic_Coastal), size = 2) + 
  geom_ribbon(data = CP_pi_abslat_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = abslat_pi_binned_means, 
             aes(x = X, y = pi_mean, color = Pelagic_Coastal), size = 8, shape = "square") + 
  geom_errorbar(data = abslat_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 0.75, size = 1.5) + 
   xlab("absolute latitude") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP)
mtdna_pi_abslat_CP_plot_annotated_both <- mtdna_pi_abslat_CP_plot_both + theme_bw() + 
  scale_y_continuous(limits = c(0, 0.015)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), 
        legend.title = element_blank())
mtdna_pi_abslat_CP_plot_annotated_both

#############################################################################################################

######### Lat figures #######

##### Predict ####

## Build prediction dataframe ##

#get lat values to predict at and scale
lat <- c(-90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 
         10, 20, 30, 40, 50, 60, 70, 80) #need to scale this using the same scale mean and sd as in real dataset
lat_scale <- (lat - attr(mtdna_small_pi$lat_scale, "scaled:center"))/attr(mtdna_small_pi$lat_scale, 
                                                                          "scaled:scale")
#build general prediction dataframe
#have to add a row for each variable in the model
pi_lat_predict_df <- as.data.frame(lat_scale)
  pi_lat_predict_df$Family <- c(rep(0, times = 9)) #random effects just set to zero as not included in predictions
  pi_lat_predict_df$Genus <- c(rep(0, times = 9))
  pi_lat_predict_df$spp <- c(rep(0, times = 9))
  pi_lat_predict_df$Source <- c(rep(0, times = 9))
  pi_lat_predict_df$MarkerName <- c(rep(0, times = 9))
  pi_lat_predict_df$Site <- c(rep(0, times = 9))
  pi_lat_predict_df$X <- lat #for plotting

## Predict for lat ##
  
pi_predict_lat <- predict(lat_model_pi, pi_lat_predict_df, re.form = NA) #re.form = NA means do NOT include random effects
  pi_lat_predict_df$predict_pi <- 10^(pi_predict_lat) #unlog pi

## Bootstrap for confidence intervals ##
#want to use bootMer on predicted dataset

pi_lat_boot <- bootMer(lat_model_pi, 
                       FUN=function(x)predict(x, pi_lat_predict_df, re.form = NA), 
                       nsim = 100)
  pi_lat_predict_df$lower_ci <- 10^(apply(pi_lat_boot$t, 2, quantile, 0.025))
  pi_lat_predict_df$upper_ci <- 10^(apply(pi_lat_boot$t, 2, quantile, 0.975))

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
lat_logpi_mean <- mtdna_small_pi[, mean(logpi), by = list(lat_round)]
  colnames(lat_logpi_mean) <- c("lat", "logpi_mean")
  lat_logpi_mean$pi_mean <- 10^(lat_logpi_mean$logpi_mean)

#calculate SE in each 10 degree band
lat_pi_SE <- mtdna_small_pi[, std.error(Pi), by = list(lat_round)] 
  colnames(lat_pi_SE) <- c("lat", "pi_SE")

#merge mean and SE dataframes together
lat_pi_binned_means <- merge(lat_logpi_mean, lat_pi_SE, by = "lat")
  lat_pi_binned_means <- lat_pi_binned_means[order(lat), ]
  lat_pi_binned_means$X <- c(-75, -65, -55, -45, -35, -25, -15, -5, 
                             5, 15, 25, 35, 45, 55, 65, 75) #for plotting, plot in MIDDLE of 10 degree band

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
  #geom_point(data = mtdna_small_pi, aes(x = lat, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = pi_lat_predict_df, 
            aes(x = X, y = predict_pi, color = "Regression"), size = 2) + 
  geom_ribbon(data = pi_lat_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci), alpha = 0.1) + 
  geom_point(data = lat_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "10-degree binned means"), size = 8, shape = "square") + 
  geom_errorbar(data = lat_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 0.75, size = 1.5) + 
  xlab("latitude") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors)
mtdna_pi_lat_plot_annotated_both <- mtdna_pi_lat_plot_both + theme_bw() + 
  scale_y_continuous(limits = c(0, 0.02)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), 
        legend.title = element_blank())
mtdna_pi_lat_plot_annotated_both

#####################################################################################################

######## Lon figures ########

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
  
#add count bin for plotting (all < number, 201 = >200)
lon_pi_binned_means$count_bin <- c(50, 50, 100, 50, 50, 100, 50, 50, 50, 100, 
                                   100, 100, 50, 50, 50, 50, 50, 100, 50, 
                                   100, 50, 50, 50, 50, 50, 50, 50, 50, 200, 201,
                                   100, 200, 50, 50, 100, 5)
lon_pi_binned_means$count_bin <- factor(lon_pi_binned_means$count_bin, levels = c(5, 50, 100, 200, 201))

#### Plot abslat ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_lon_plot_both <- ggplot() + 
  geom_line(data = lon_eff_data, 
            aes(x = lon, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = lon_eff_data, 
              aes(x = lon, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = lon_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "10-degree binned means", size = count_bin), shape = "square") + 
  geom_errorbar(data = lon_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("longitude") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 0.015)) + 
  scale_size_manual(values = c(6, 8, 10, 12, 16), 
                    labels = c("\u22645", "\u226450", "\u2264100", "\u2264200", "\u2265200"))
mtdna_pi_lon_plot_annotated_both <- mtdna_pi_lon_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_lon_plot_annotated_both

## With pelagic and coastal separate ##

#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_pi_lon_CP_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = lon, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = CP_pi_lon_predict_df, 
            aes(x = X, y = predict_pi, color = Pelagic_Coastal), size = 2) + 
  geom_ribbon(data = CP_pi_lon_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = lon_pi_binned_means, 
             aes(x = X, y = pi_mean, color = Pelagic_Coastal), size = 8, shape = "square") + 
  geom_errorbar(data = lon_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 0.75, size = 1.5) + 
  xlab("longitude") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP)
mtdna_pi_lon_CP_plot_annotated_both <- mtdna_pi_lon_CP_plot_both + theme_bw() + 
  scale_y_continuous(limits = c(0, 0.03)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), 
        legend.title = element_blank())
mtdna_pi_lon_CP_plot_annotated_both

#CI --> Zoe did it this way
#used Ben Bolker's tutorial https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions

#############################################################################################################

######## SST figures ########

#### Build SST models ####

## Clean up dataframe ##

#remove rows with SSTmean = NA
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$sst.BO_sstmean != "NA") #should look at where these are...

#log transform sst data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$sst.BO_sstmean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi$logsstmean <- log10(mtdna_small_pi$sst.BO_sstmean)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logsstmean != "NaN")

## SST mean models ##
  
SSTmean_model_pi <- lm(logpi ~ logsstmean, data = mtdna_small_pi, 
                        na.action = "na.fail") #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

SSTmin_CP_model_pi <- lmer(logpi ~ logsstmin + Pelagic_Coastal + Pelagic_Coastal:logsstmin + 
                             (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + (1|Site), 
                           REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", 
                           control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

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

#create column to fill with the rounded SSTmin
mtdna_small_pi$SSTmean_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$SSTmean_round[i] <- DescTools::RoundTo(mtdna_small_pi$sst.BO_sstmean[i], 
                                                        multiple = 5, FUN = floor)} #rounding DOWN (e.g. SSTmin 5-10 gets value of 5)
}

mtdna_small_pi$SSTmean_round <- as.factor(mtdna_small_pi$SSTmean_round)

## Calculate means and SE within each 10 degree band ##

mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
SSTmean_pi_mean <- mtdna_small_pi[, mean(logpi), by = SSTmean_round]
  colnames(SSTmean_pi_mean) <- c("SSTmean", "pi_mean")
  SSTmean_pi_mean$pi_mean <- 10^(SSTmean_pi_mean$pi_mean) #unlog mean

#calculate SE in each 10 degree band
SSTmean_pi_SE <- mtdna_small_pi[, std.error(Pi), by = SSTmean_round]
  colnames(SSTmean_pi_SE) <- c("SSTmean", "pi_SE")

#count observations in each 10 degree band
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

#add count bin for plotting (all < number, 201 = >200)
SSTmean_pi_binned_means$count_bin <- c(50, 100, 201, 201, 201, 201, 50)
SSTmean_pi_binned_means$count_bin <- factor(SSTmean_pi_binned_means$count_bin, levels = c(50, 100, 201))

#### Plot SSTmean ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_SSTmean_plot_both <- ggplot() + 
  geom_line(data = SSTmean_eff_data, 
            aes(x = sstmean, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = SSTmean_eff_data, 
              aes(x = sstmean, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = SSTmean_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "10-degree binned means", size = count_bin), shape = "square") + 
  geom_errorbar(data = SSTmean_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("SST mean") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 0.015)) + 
  scale_size_manual(values = c(8, 12, 16), 
                    labels = c("\u226450", "\u2264100", "\u2265200"))
mtdna_pi_SSTmean_plot_annotated_both <- mtdna_pi_SSTmean_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_SSTmean_plot_annotated_both

## With pelagic and coastal separately ##

#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_pi_SSTmin_CP_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = sst.BO_sstmin, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = CP_pi_SSTmin_predict_df, 
            aes(x = X, y = predict_pi, color = Pelagic_Coastal), size = 2) + 
  geom_ribbon(data = CP_pi_SSTmin_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = SSTmin_pi_binned_means, 
             aes(x = X, y = pi_mean, color = Pelagic_Coastal), size = 8, shape = "square") + 
  geom_errorbar(data = SSTmin_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 0.75, size = 1.5) + 
  xlab("SSTmin") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP)
mtdna_pi_SSTmin_CP_plot_annotated_both <- mtdna_pi_SSTmin_CP_plot_both + theme_bw() + 
  scale_y_continuous(limits = c(0, 0.03)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), 
        legend.title = element_blank())
mtdna_pi_SSTmin_CP_plot_annotated_both

###########################################################################################################

######## Dissox figures ########

#### Build dissox models ####

## Clean up dataframes ##

#remove rows with dissox = NA
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$BO_dissox != "NA") #should look at where these are...

#log transform dissox
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$BO_dissox != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi$logdissox <- log10(mtdna_small_pi$BO_dissox)
  mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logdissox != "NaN")

## dissox mean models ##
  
dissox_model_pi <- lmer(logpi ~ logdissox + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                          (1|Site), REML = FALSE, data = mtdna_small_pi, 
                        na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

dissox_CP_model_pi <- lmer(logpi ~ logdissox + Pelagic_Coastal + Pelagic_Coastal:logdissox + 
                             (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + (1|Site), 
                           REML = FALSE, data = mtdna_small_pi, 
                           na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#### Predict ####

## Build prediction data frame ##

#get dissox mean values to predict at and log transform
#predict for dissox (range is 4.213 to 8.426)
dissox <- c(4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5)
  logdissox <- log10(dissox)

#build general prediction dataframe
#have to add a row for each variable in the model
pi_dissox_predict_df <- as.data.frame(logdissox)
  pi_dissox_predict_df$Family <- c(rep(0, times = 10)) #random effects just set to zero as not included in predictions
  pi_dissox_predict_df$Genus <- c(rep(0, times = 10))
  pi_dissox_predict_df$spp <- c(rep(0, times = 10))
  pi_dissox_predict_df$Source <- c(rep(0, times = 10))
  pi_dissox_predict_df$MarkerName <- c(rep(0, times = 10))
  pi_dissox_predict_df$Site <- c(rep(0, times = 10))
  pi_dissox_predict_df$X <- dissox #need this for plotting

#build prediction dataframe for pelagics
Pel_pi_dissox_predict_df <- pi_dissox_predict_df
  Pel_pi_dissox_predict_df$Pelagic_Coastal <- c(rep("Pelagic", times = 10))

#build prediction dataframe for coastals
Coast_pi_dissox_predict_df <- pi_dissox_predict_df
  Coast_pi_dissox_predict_df$Pelagic_Coastal <- c(rep("Coastal", times = 10))

## Predict for dissox ##
  
pi_predict_dissox <- predict(dissox_model_pi, pi_dissox_predict_df, re.form = NA) #re.form = NA means do NOT include random effects
  pi_dissox_predict_df$predict_pi <- 10^(pi_predict_dissox) #unlog pi

Pel_pi_predict_dissox <- predict(dissox_CP_model_pi, Pel_pi_dissox_predict_df, re.form = NA)
  Pel_pi_dissox_predict_df$predict_pi <- 10^(Pel_pi_predict_dissox) #unlog pi

Coast_pi_predict_dissox <- predict(dissox_CP_model_pi, Coast_pi_dissox_predict_df, re.form = NA)
  Coast_pi_dissox_predict_df$predict_pi <- 10^(Coast_pi_predict_dissox) #unlog pi

## Bootstrap for confidence intervals ##
#want to use bootMer on predicted dataset

pi_dissox_boot <- bootMer(dissox_model_pi, 
                          FUN=function(x)predict(x, pi_dissox_predict_df, re.form = NA), 
                          nsim = 100)
  pi_dissox_predict_df$lower_ci <- 10^(apply(pi_dissox_boot$t, 2, quantile, 0.025))
  pi_dissox_predict_df$upper_ci <- 10^(apply(pi_dissox_boot$t, 2, quantile, 0.975))
  
Pel_pi_dissox_boot <- bootMer(dissox_CP_model_pi, 
                          FUN=function(x)predict(x, Pel_pi_dissox_predict_df, re.form = NA), 
                          nsim = 100)
  Pel_pi_dissox_predict_df$lower_ci <- 10^(apply(Pel_pi_dissox_boot$t, 2, quantile, 0.025))
  Pel_pi_dissox_predict_df$upper_ci <- 10^(apply(Pel_pi_dissox_boot$t, 2, quantile, 0.975))
  
Coast_pi_dissox_boot <- bootMer(dissox_CP_model_pi, 
                           FUN=function(x)predict(x, Coast_pi_dissox_predict_df, re.form = NA), 
                           nsim = 100)
  Coast_pi_dissox_predict_df$lower_ci <- 10^(apply(Coast_pi_dissox_boot$t, 2, quantile, 0.025))
  Coast_pi_dissox_predict_df$upper_ci <- 10^(apply(Coast_pi_dissox_boot$t, 2, quantile, 0.975))

#for CP, bind pelagic and coastal dataframes back together
CP_pi_dissox_predict_df <- rbind(Pel_pi_dissox_predict_df, Coast_pi_dissox_predict_df)
  
#### Calculate means from raw data ####

## Round dissox mean DOWN to nearest multiple of 0.5 ##

#create column to fill with the rounded dissox
mtdna_small_pi$dissox_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$dissox_round[i] <- DescTools::RoundTo(mtdna_small_pi$BO_dissox[i], 
                                                        multiple = 0.5, FUN = floor)} #rounding DOWN (e.g. dissox 1-1.5 gets value of 1)
}

mtdna_small_pi$dissox_round <- as.factor(mtdna_small_pi$dissox_round)

## Calculate means and SE within each 0.5 band ##
#written so calculates Pelagic & Coastal means separately

mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so that can use data.table functions

#calculate mean in each 0.5 band
dissox_logpi_mean <- mtdna_small_pi[, mean(logpi), by = list(dissox_round, Pelagic_Coastal)]
  colnames(dissox_logpi_mean) <- c("dissox", "Pelagic_Coastal", "logpi_mean")
  dissox_logpi_mean$pi_mean <- 10^(dissox_logpi_mean$logpi_mean) #unlog mean

#calculate SE in each 0.5 band  
dissox_pi_SE <- mtdna_small_pi[, std.error(Pi), by = list(dissox_round, Pelagic_Coastal)]
  colnames(dissox_pi_SE) <- c("dissox", "Pelagic_Coastal", "pi_SE")

#merge mean and SE dataframes together
dissox_pi_binned_means <- merge(dissox_logpi_mean, dissox_pi_SE, by = c("dissox", "Pelagic_Coastal"))
dissox_pi_binned_means <- dissox_pi_binned_means[order(dissox), ]
dissox_pi_binned_means$X <- c(4.25, 4.25, 4.75, 4.75, 5.25, 5.25, 5.75, 5.75, 6.25, 6.25, 6.75, 6.75, 
                              7.25, 7.25, 7.75, 7.75, 8.25, 8.25) #for plotting, plot in MIDDLE of 0.5 band

#add error bars (standard error)
dissox_pi_binned_means$mean_lowerSE <- dissox_pi_binned_means$pi_mean - 
  dissox_pi_binned_means$pi_SE
dissox_pi_binned_means$mean_upperSE <- dissox_pi_binned_means$pi_mean + 
  dissox_pi_binned_means$pi_SE

#### Plot dissox ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_dissox_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = BO_dissox, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = pi_dissox_predict_df, 
            aes(x = X, y = predict_pi, color = "Regression"), size = 2) + 
  geom_ribbon(data = pi_dissox_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci), alpha = 0.1) + 
  geom_point(data = dissox_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "10-degree binned means"), size = 8, shape = "square") + 
  geom_errorbar(data = dissox_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = "10-degree binned means"), 
                width = 0.75, size = 1.5) + 
  xlab("Dissolved oxygen (mean)") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors)
mtdna_pi_dissox_plot_annotated_both <- mtdna_pi_dissox_plot_both + theme_bw() + 
  scale_y_continuous(limits = c(0, 0.02)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), 
        legend.title = element_blank())
mtdna_pi_dissox_plot_annotated_both

## With pelagic and coastal separately ##

#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_pi_dissox_CP_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = BO_dissox, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = CP_pi_dissox_predict_df, 
            aes(x = X, y = predict_pi, color = Pelagic_Coastal), size = 2) + 
  geom_ribbon(data = CP_pi_dissox_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = dissox_pi_binned_means, 
             aes(x = X, y = pi_mean, color = Pelagic_Coastal), size = 8, shape = "square") + 
  geom_errorbar(data = dissox_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 0.05, size = 1.5) + 
  xlab("Dissolved oxygen (mean)") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP)
mtdna_pi_dissox_CP_plot_annotated_both <- mtdna_pi_dissox_CP_plot_both + theme_bw() + 
  scale_y_continuous(limits = c(0, 0.03)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), 
        legend.title = element_blank())
mtdna_pi_dissox_CP_plot_annotated_both

#############################################################################################################

######## ChloroA figures ########

#### Build ChloroA models ####

## Clean up data frame ##

#remove rows with chloroa = NA
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomax != "NA") #should look at where these are...

#log transform chloroA data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomax != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
mtdna_small_pi$logchlomax <- log10(mtdna_small_pi$chloroA.BO_chlomax)
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlomax != "NaN")

##chlomax models ##

chloroAmax_model_pi <- lmer(logpi ~ logchlomax + I(logchlomax^2) + (1|Family/Genus/spp) + (1|Source) + 
                              (1|MarkerName) + (1|Site), REML = FALSE, data = mtdna_small_pi, 
                            na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

chloroAmax_CP_model_pi <- lmer(logpi ~ logchlomax + I(logchlomax^2) + Pelagic_Coastal + 
                                 Pelagic_Coastal:logchlomax + Pelagic_Coastal:I(logchlomax^2) + 
                                 (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + (1|Site), 
                               REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", 
                               control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#### Predict ####

## Build prediction data frame ##

#get chlomax values to predict at and logtransform
#predict for chlomax (range is 0.037 to 61.783)
chlomax <- c(0.01, 10, 20, 30, 40, 50, 60) #starting at .01 bc log10 of 0 is Inf
  logchlomax <- log10(chlomax)

#build general prediction dataframe
#have to add a row for each variable in the model
pi_chlomax_predict_df <- as.data.frame(logchlomax)
  pi_chlomax_predict_df$Family <- c(rep(0, times = 7)) #random effects just set to zero as not included in predictions
  pi_chlomax_predict_df$Genus <- c(rep(0, times = 7))
  pi_chlomax_predict_df$spp <- c(rep(0, times = 7))
  pi_chlomax_predict_df$Source <- c(rep(0, times = 7))
  pi_chlomax_predict_df$MarkerName <- c(rep(0, times = 7))
  pi_chlomax_predict_df$Site <- c(rep(0, times = 7))
  pi_chlomax_predict_df$X <- chlomax #need this for plotting

#build prediction dataframe for pelagics
Pel_pi_chlomax_predict_df <- pi_chlomax_predict_df
  Pel_pi_chlomax_predict_df$Pelagic_Coastal <- c(rep("Pelagic", times = 7))
  
#build prediction dataframe for coastals
Coast_pi_chlomax_predict_df <- pi_chlomax_predict_df
  Coast_pi_chlomax_predict_df$Pelagic_Coastal <- c(rep("Coastal", times = 7))

## Predict for chlomax ##
  
pi_predict_chlomax <- predict(chloroAmax_model_pi, pi_chlomax_predict_df, re.form = NA) #re.form = NA means do NOT include random effects
  pi_chlomax_predict_df$predict_pi <- 10^(pi_predict_chlomax) #unlog pi

Pel_pi_predict_chlomax <- predict(chloroAmax_CP_model_pi, Pel_pi_chlomax_predict_df, re.form = NA)
  Pel_pi_chlomax_predict_df$predict_pi <- 10^(Pel_pi_predict_chlomax) #unlog pi

Coast_pi_predict_chlomax <- predict(chloroAmax_CP_model_pi, Coast_pi_chlomax_predict_df, re.form = NA)
  Coast_pi_chlomax_predict_df$predict_pi <- 10^(Coast_pi_predict_chlomax) #unlog pi

## Bootstrap for confidence intervals ##
#want to use bootMer on predicted dataset

pi_chlomax_boot <- bootMer(chloroAmax_model_pi, 
                          FUN=function(x)predict(x, pi_chlomax_predict_df, re.form = NA), 
                          nsim = 100)
  pi_chlomax_predict_df$lower_ci <- 10^(apply(pi_chlomax_boot$t, 2, quantile, 0.025))
  pi_chlomax_predict_df$upper_ci <- 10^(apply(pi_chlomax_boot$t, 2, quantile, 0.975))

Pel_pi_chlomax_boot <- bootMer(chloroAmax_CP_model_pi, 
                           FUN=function(x)predict(x, Pel_pi_chlomax_predict_df, re.form = NA), 
                           nsim = 100)
  Pel_pi_chlomax_predict_df$lower_ci <- 10^(apply(Pel_pi_chlomax_boot$t, 2, quantile, 0.025))
  Pel_pi_chlomax_predict_df$upper_ci <- 10^(apply(Pel_pi_chlomax_boot$t, 2, quantile, 0.975))

Coast_pi_chlomax_boot <- bootMer(chloroAmax_CP_model_pi, 
                           FUN=function(x)predict(x, Coast_pi_chlomax_predict_df, re.form = NA), 
                           nsim = 100)
  Coast_pi_chlomax_predict_df$lower_ci <- 10^(apply(Coast_pi_chlomax_boot$t, 2, quantile, 0.025))
  Coast_pi_chlomax_predict_df$upper_ci <- 10^(apply(Coast_pi_chlomax_boot$t, 2, quantile, 0.975))

#for CP, bind pelagic and coastal dataframes back together
CP_pi_chlomax_predict_df <- rbind(Pel_pi_chlomax_predict_df, Coast_pi_chlomax_predict_df)
  
#### Calculate means from raw data ####

## Round chlomax DOWN to nearest multiple of 10 ##

#create column to fill with the rounded chlomax
mtdna_small_pi$chlomax_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$chlomax_round[i] <- DescTools::RoundTo(mtdna_small_pi$chloroA.BO_chlomax[i], 
                                                         multiple = 10, FUN = floor)} #rounding DOWN (e.g. chlomax 10-20 gets value of 10)
}

mtdna_small_pi$chlomax_round <- as.factor(mtdna_small_pi$chlomax_round)

## Calculate means and SE within each 10 unit band ##
#written so calculates Pelagic & Coastal means separately

mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so that can use data.table functions

#calculate mean in each 10 unit band
chlomax_logpi_mean <- mtdna_small_pi[, mean(logpi), by = list(chlomax_round, Pelagic_Coastal)]
  colnames(chlomax_logpi_mean) <- c("chlomax", "Pelagic_Coastal", "logpi_mean")
  chlomax_logpi_mean$pi_mean <- 10^(chlomax_logpi_mean$logpi_mean) #unlog mean

#calculate SE in each 10 unit man
chlomax_pi_SE <- mtdna_small_pi[, std.error(Pi), by = list(chlomax_round, Pelagic_Coastal)]
  colnames(chlomax_pi_SE) <- c("chlomax", "Pelagic_Coastal", "pi_SE")

#merge mean and SE dataframes together
chlomax_pi_binned_means <- merge(chlomax_logpi_mean, chlomax_pi_SE, by = c("chlomax", "Pelagic_Coastal"))
chlomax_pi_binned_means <- chlomax_pi_binned_means[order(chlomax), ]
  chlomax_pi_binned_means$X <- c(5, 5, 15, 15, 25, 25, 35, 35, 45, 55, 60.5) #for plotting, plot in MIDDLE of 10 unit band

#add error bars (standard error)
chlomax_pi_binned_means$mean_lowerSE <- chlomax_pi_binned_means$pi_mean - 
  chlomax_pi_binned_means$pi_SE
chlomax_pi_binned_means$mean_upperSE <- chlomax_pi_binned_means$pi_mean + 
  chlomax_pi_binned_means$pi_SE

#### Plot chlomax ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_chlomax_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = chloroA.BO_chlomax, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = pi_chlomax_predict_df, 
            aes(x = X, y = predict_pi, color = "Regression"), size = 2) + 
  geom_ribbon(data = pi_chlomax_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci), alpha = 0.1) + 
  geom_point(data = chlomax_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "10-degree binned means"), size = 8, shape = "square") + 
  geom_errorbar(data = chlomax_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = "10-degree binned means"), 
                width = 0.75, size = 1.5) + 
  xlab("Chlorophyll A (max)") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors)
mtdna_pi_chlomax_plot_annotated_both <- mtdna_pi_chlomax_plot_both + theme_bw() + 
  scale_y_continuous(limits = c(0, 0.02)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), 
        legend.title = element_blank())
mtdna_pi_chlomax_plot_annotated_both

##With pelagic and coastal separately ##

#for legend
colors_CP <- c(Coastal = "#3F6DAA", Pelagic = "black")

#plot CP
mtdna_pi_chlomax_CP_plot_both <- ggplot() + 
  #geom_point(data = mtdna_small_pi, aes(x = chloroA.BO_chlomax, y = Pi, color = "Raw data"), 
  #           size = 4, alpha = 0.25) + 
  geom_line(data = CP_pi_chlomax_predict_df, 
            aes(x = X, y = predict_pi, color = Pelagic_Coastal), size = 2) + 
  geom_ribbon(data = CP_pi_chlomax_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci, color = Pelagic_Coastal), alpha = 0.1) + 
  geom_point(data = chlomax_pi_binned_means, 
             aes(x = X, y = pi_mean, color = Pelagic_Coastal), size = 8, shape = "square") + 
  geom_errorbar(data = chlomax_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, color = Pelagic_Coastal), 
                width = 0.05, size = 1.5) + 
  xlab("Chlorophyll A (max)") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors_CP)
mtdna_pi_chlomax_CP_plot_annotated_both <- mtdna_pi_chlomax_CP_plot_both + theme_bw() + 
  scale_y_continuous(limits = c(0, 0.03)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
mtdna_pi_chlomax_CP_plot_annotated_both
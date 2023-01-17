######################## Script for Building Regression Figures for mtDNA Hd  ##############################

#Uses final mtdna hd binomial regression models
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

###################################################################################################

######## Build lat, abslat & lon models ########
  
#### clean up data ####
  
#scale geographic variables
mtdna_small_hd$lat_scale <- as.numeric(scale(mtdna_small_hd$lat))
mtdna_small_hd$abslat_scale <- as.numeric(scale(mtdna_small_hd$abslat))

#scale bp
mtdna_small_hd$bp_scale <- as.numeric(scale(as.numeric(mtdna_small_hd$bp)))

#convert lon to radians
mtdna_small_hd$lon_360 <- mtdna_small_hd$lon + 180 #convert (-180,180) to (0,360)
  mtdna_small_hd$lon_rad <- (2*pi*mtdna_small_hd$lon_360)/360

#calculate success and failure
mtdna_small_hd$success <- round(mtdna_small_hd$He*mtdna_small_hd$n) #essentially, number of heterozygotes
  mtdna_small_hd$failure <- round((1 - mtdna_small_hd$He)*mtdna_small_hd$n) #number of homozygotes

#### Models ####
#for prediction plots, do NOT want any random effects
#also not using different (faster) optimizer

lat_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + lat_scale + 
                        I(lat_scale^2), family = binomial, data = mtdna_small_hd, 
                      na.action = "na.fail")

lat_CP_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + lat_scale + 
                        I(lat_scale^2) + Pelagic_Coastal + Pelagic_Coastal:lat_scale + 
                         Pelagic_Coastal:I(lat_scale^2), family = binomial, data = mtdna_small_hd, 
                      na.action = "na.fail")
  
abslat_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + abslat_scale + 
                           I(abslat_scale^2), family = binomial, data = mtdna_small_hd, 
                         na.action = "na.fail")

abslat_CP_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + 
                              Pelagic_Coastal + abslat_scale + I(abslat_scale^2) + 
                              Pelagic_Coastal:abslat_scale + Pelagic_Coastal:I(abslat_scale^2), 
                            family = binomial, data = mtdna_small_hd, na.action = "na.fail")

lon_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + sin(lon_rad) + cos(lon_rad), family = binomial, data = mtdna_small_hd, 
                      na.action = "na.fail")

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
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$sst.BO_sstmean != "NA") #subset to only rows with sstmean data

mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$sst.BO_sstmean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd$logsstmean <- log10(mtdna_small_hd$sst.BO_sstmean)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logsstmean != "NaN")

sstmean_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + logsstmean, 
                        family = binomial, data = mtdna_small_hd, na.action = "na.fail")

sstmean_CP_model_hd <- glm(cbind(success, failure) ~ bp_scale + range_position + logsstmean + Pelagic_Coastal + 
                                        Pelagic_Coastal:logsstmean, family = binomial, data = mtdna_small_hd, na.action = "na.fail")

#############################################################################################################

######### Abs Lat model figures #######

#### Predict ####

#marginal effects
abslat_eff <- plot_model(abslat_model_hd, type = "pred", #pred.type = "re" incorporates random effects
                      terms = "abslat_scale [all]")

#pull out marginal effects dataframe
abslat_eff_data <- as.data.frame(abslat_eff$data)

#unscale lat
#use same scaled:center & scaled:scale from original data
abslat_scale <- scale(mtdna_small_hd$abslat) #bc had to convert to numeric to run model/calculate marginal effects
abslat_eff_data$abslat <- (abslat_eff_data$x * attr(abslat_scale, "scaled:scale")) + 
  attr(abslat_scale, "scaled:center")

#### Calculate means from raw data ####

## Round abslat DOWN to nearest multiple of 10 ##

#create column to fill with the rounded abslat
mtdna_small_hd$abslat_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_hd)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_hd$abslat_round[i] <- DescTools::RoundTo(mtdna_small_hd$abslat[i], 
                                                        multiple = 10, FUN = floor)} #rounding DOWN (e.g. abslat 10-20 gets value of 10)
}

mtdna_small_hd$abslat_round <- as.factor(mtdna_small_hd$abslat_round)

## Calculate means and SE within each 10 degree band ##

mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
abslat_hd_mean <- mtdna_small_hd[, mean(He), by = abslat_round]
CP_abslat_hd_mean <- mtdna_small_hd[, mean(He), by = list(abslat_round, Pelagic_Coastal)]
  colnames(abslat_hd_mean) <- c("abslat", "hd_mean")
  colnames(CP_abslat_hd_mean) <- c("abslat", "Pelagic_Coastal", "hd_mean")

#calculate SE in each 10 degree band
abslat_hd_SE <- mtdna_small_hd[, std.error(He), by = abslat_round]
CP_abslat_hd_SE <- mtdna_small_hd[, std.error(He), by = list(abslat_round, Pelagic_Coastal)]
  colnames(abslat_hd_SE) <- c("abslat", "hd_SE")
  colnames(CP_abslat_hd_SE) <- c("abslat", "Pelagic_Coastal", "hd_SE")
  
#count observations in each 10 degree band
abslat_hd_count <- mtdna_small_hd[, .N, by = abslat_round]
CP_abslat_hd_count <- mtdna_small_hd[, .N, by = list(abslat_round, Pelagic_Coastal)]
  colnames(abslat_hd_count) <- c("abslat", "hd_count")
  colnames(CP_abslat_hd_count) <- c("abslat", "Pelagic_Coastal", "hd_count")

#merge mean and SE dataframes together
abslat_hd_binned_means <- list(abslat_hd_mean, abslat_hd_SE, abslat_hd_count) %>% 
  reduce(full_join, by = "abslat")
  abslat_hd_binned_means$abslat <- as.numeric(as.character(abslat_hd_binned_means$abslat))
  abslat_hd_binned_means <- abslat_hd_binned_means[order(abslat), ]
  abslat_hd_binned_means$X <- abslat_hd_binned_means$abslat + 5 #for plotting, plot in MIDDLE of 10 degree band

CP_lat_hd_binned_means <- merge(CP_lat_hd_mean, CP_lat_hd_SE, by = c("lat", "Pelagic_Coastal"))
  CP_lat_hd_binned_means <- merge(CP_lat_hd_binned_means, CP_lat_hd_count, by = c("lat", "Pelagic_Coastal"))
  colnames(CP_lat_hd_binned_means) <- c("lat", "Pelagic_Coastal", "hd_mean", "hd_SE", "hd_count")
  CP_lat_hd_binned_means$lat <- as.numeric(as.character(CP_lat_hd_binned_means$lat))
  CP_lat_hd_binned_means <- CP_lat_hd_binned_means[order(lat), ]
  CP_lat_hd_binned_means$X <- CP_lat_hd_binned_means$lat + 5 #for plotting, plot in MIDDLE of 10 degree band

#calculate error bars (standard error)
abslat_hd_binned_means$mean_lowerSE <- abslat_hd_binned_means$hd_mean - 
  abslat_hd_binned_means$hd_SE
abslat_hd_binned_means$mean_upperSE <- abslat_hd_binned_means$hd_mean + 
  abslat_hd_binned_means$hd_SE

CP_lat_hd_binned_means$mean_lowerSE <- CP_lat_hd_binned_means$hd_mean - 
  CP_lat_hd_binned_means$hd_SE
CP_lat_hd_binned_means$mean_upperSE <- CP_lat_hd_binned_means$hd_mean + 
  CP_lat_hd_binned_means$hd_SE

#add count bin for plotting (all < number, 201 = >200)
abslat_hd_binned_means$count_bin <- c(201, 201, 201, 201, 200, 100, 50, 50)
  abslat_hd_binned_means$count_bin <- factor(abslat_hd_binned_means$count_bin, levels = c(50, 100, 200, 201))

CP_lat_hd_binned_means$count_bin <- c(5, 5, 5, 5, 50, 50, 50, 50, 100, 50, 200, 5, 100, 50, 200, 
                                      50, 200, 50, 201, 50, 201, 50, 200, 50, 100, 50, 50, 50, 50)
  CP_lat_hd_binned_means$count_bin <- factor(CP_lat_hd_binned_means$count_bin, levels = c(5, 50, 100, 200, 201))

#### Plot abslat ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_hd_abslat_plot_both <- ggplot() + 
  geom_line(data = abslat_eff_data, 
            aes(x = abslat, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = abslat_eff_data, 
             aes(x = abslat, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = abslat_hd_binned_means, 
             aes(x = X, y = hd_mean, color = "10-degree binned means", size = count_bin), shape = "square") + 
  geom_errorbar(data = abslat_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("absolute latitude") + ylab("mtdna Hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_size_manual(values = c(8, 12, 14, 16), 
                    labels = c("\u226450", "\u2264100", "\u2264200", "\u2265200"))
mtdna_hd_abslat_plot_annotated_both <- mtdna_hd_abslat_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 70, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 60, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 50), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_abslat_plot_annotated_both

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

######### Lat model figures #######

#### Predict ####

#marginal effects
lat_eff <- plot_model(lat_model_hd, type = "pred", 
                         terms = "lat_scale [all]")

#pull out marginal effects dataframe
lat_eff_data <- as.data.frame(lat_eff$data)

#unscale lat
#use same scaled:center & scaled:scale from original data
lat_scale <- scale(mtdna_small_hd$lat) #sometimes, original scale in dataframe doesn't store attributes
lat_eff_data$lat <- (lat_eff_data$x * attr(lat_scale, "scaled:scale")) + 
  attr(lat_scale, "scaled:center")

#### Calculate means from raw data ####

## Round lat DOWN to nearest multiple of 10 ##

#create column to fill with the rounded abslat
mtdna_small_hd$lat_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_hd)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_hd$lat_round[i] <- DescTools::RoundTo(mtdna_small_hd$lat[i], 
                                                     multiple = 10, FUN = floor)} #rounding DOWN (e.g. abslat 10-20 gets value of 10)
}

mtdna_small_hd$lat_round <- as.factor(mtdna_small_hd$lat_round)

## Calculate means and SE within each 10 degree band ##

mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
lat_hd_mean <- mtdna_small_hd[, mean(He), by = lat_round]
CP_lat_hd_mean <- mtdna_small_hd[, mean(He), by = list(lat_round, Pelagic_Coastal)]
colnames(lat_hd_mean) <- c("lat", "hd_mean")
colnames(CP_lat_hd_mean) <- c("lat", "Pelagic_Coastal", "hd_mean")

#calculate SE in each 10 degree band
lat_hd_SE <- mtdna_small_hd[, std.error(He), by = lat_round]
CP_lat_hd_SE <- mtdna_small_hd[, std.error(He), by = list(lat_round, Pelagic_Coastal)]
colnames(lat_hd_SE) <- c("lat", "hd_SE")
colnames(CP_lat_hd_SE) <- c("lat", "Pelagic_Coastal", "hd_SE")

#count observations in each 10 degree band
lat_hd_count <- mtdna_small_hd[, .N, by = lat_round]
CP_lat_hd_count <- mtdna_small_hd[, .N, by = list(lat_round, Pelagic_Coastal)]
colnames(lat_hd_count) <- c("lat", "hd_count")
colnames(CP_lat_hd_count) <- c("lat", "Pelagic_Coastal", "hd_count")

#merge mean and SE dataframes together
lat_hd_binned_means <- list(lat_hd_mean, lat_hd_SE, lat_hd_count) %>% 
  reduce(full_join, by = "lat")
lat_hd_binned_means$lat <- as.numeric(as.character(lat_hd_binned_means$lat))
lat_hd_binned_means <- lat_hd_binned_means[order(lat), ]
lat_hd_binned_means$X <- lat_hd_binned_means$lat + 5 #for plotting, plot in MIDDLE of 10 degree band

CP_lat_hd_binned_means <- merge(CP_lat_hd_mean, CP_lat_hd_SE, by = c("lat", "Pelagic_Coastal"))
CP_lat_hd_binned_means <- merge(CP_lat_hd_binned_means, CP_lat_hd_count, by = c("lat", "Pelagic_Coastal"))
colnames(CP_lat_hd_binned_means) <- c("lat", "Pelagic_Coastal", "hd_mean", "hd_SE", "hd_count")
CP_lat_hd_binned_means$lat <- as.numeric(as.character(CP_lat_hd_binned_means$lat))
CP_lat_hd_binned_means <- CP_lat_hd_binned_means[order(lat), ]
CP_lat_hd_binned_means$X <- CP_lat_hd_binned_means$lat + 5 #for plotting, plot in MIDDLE of 10 degree band

#calculate error bars (standard error)
lat_hd_binned_means$mean_lowerSE <- lat_hd_binned_means$hd_mean - 
  lat_hd_binned_means$hd_SE
lat_hd_binned_means$mean_upperSE <- lat_hd_binned_means$hd_mean + 
  lat_hd_binned_means$hd_SE

CP_lat_hd_binned_means$mean_lowerSE <- CP_lat_hd_binned_means$hd_mean - 
  CP_lat_hd_binned_means$hd_SE
CP_lat_hd_binned_means$mean_upperSE <- CP_lat_hd_binned_means$hd_mean + 
  CP_lat_hd_binned_means$hd_SE

#add count bin for plotting (all < number, 201 = >200)
lat_hd_binned_means$count_bin <- c(5, 5, 5, 50, 100, 100, 200, 200, 200, 200, 201, 201, 200, 
                                   100, 50, 50)
lat_hd_binned_means$count_bin <- factor(lat_hd_binned_means$count_bin, levels = c(5, 50, 100, 200, 201))

CP_lat_hd_binned_means$count_bin <- c(5, 5, 5, 5, 50, 50, 50, 50, 100, 50, 200, 5, 100, 50, 200, 
                                      50, 200, 50, 201, 50, 201, 50, 200, 50, 100, 50, 50, 50, 50)
CP_lat_hd_binned_means$count_bin <- factor(CP_lat_hd_binned_means$count_bin, levels = c(5, 50, 100, 200, 201))

#### Plot lat ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")
mtdna_small_hd$X <- mtdna_small_hd$lat

#plot
mtdna_hd_lat_plot_both <- ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey", size = 3) + 
  geom_line(data = lat_eff_data, 
              aes(x = lat, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = lat_eff_data, 
              aes(x = lat, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = lat_hd_binned_means, 
             aes(x = X, y = hd_mean, color = "10-degree binned means", size = count_bin), shape = "square") + 
  geom_errorbar(data = lat_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("latitude") + ylab("mtdna Hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_size_manual(values = c(6, 8, 12, 14, 16), 
                    labels = c("\u22645", "\u226450", "\u2264100", "\u2264200", "\u2265200"))
mtdna_hd_lat_plot_annotated_both <- mtdna_hd_lat_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 70, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 60, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 50), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_lat_plot_annotated_both

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
lon_eff <- plot_model(lon_model_hd, type = "pred", terms = "lon_rad [all]")
  lon_eff_data <- as.data.frame(lon_eff$data)
  lon_eff_data$lon <- ((360*lon_eff_data$x)/(2*pi))-180
  
#### Calculate means from raw data ####

## Round lon DOWN to nearest multiple of 10 ##

#create column to fill with the rounded lon
mtdna_small_hd$lon_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_hd)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_hd$lon_round[i] <- DescTools::RoundTo(mtdna_small_hd$lon[i], 
                                                     multiple = 10, FUN = floor)} #rounding DOWN (e.g. abslat 10-20 gets value of 10)
}

mtdna_small_hd$lon_round <- as.factor(mtdna_small_hd$lon_round)

## Calculate means and SE within each 10 degree band ##

mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
lon_hd_mean <- mtdna_small_hd[, mean(He), by = lon_round]
  CP_lon_hd_mean <- mtdna_small_hd[, mean(He), by = list(lon_round, Pelagic_Coastal)]
  colnames(lon_hd_mean) <- c("lon", "hd_mean")
  colnames(CP_lon_hd_mean) <- c("lon", "Pelagic_Coastal", "hd_mean")

#calculate SE in each 10 degree band
lon_hd_SE <- mtdna_small_hd[, std.error(He), by = lon_round]
  CP_lon_hd_SE <- mtdna_small_hd[, std.error(He), by = list(lon_round, Pelagic_Coastal)]
  colnames(lon_hd_SE) <- c("lon", "hd_SE")
  colnames(CP_lon_hd_SE) <- c("lon", "Pelagic_Coastal", "hd_SE")

#count observations in each 10 degree band
lon_hd_count <- mtdna_small_hd[, .N, by = lon_round]
  CP_lon_hd_count <- mtdna_small_hd[, .N, by = list(lon_round, Pelagic_Coastal)]
  colnames(lon_hd_count) <- c("lon", "hd_count")
  colnames(CP_lon_hd_count) <- c("lon", "Pelagic_Coastal", "hd_count")

#merge mean and SE dataframes together
lon_hd_binned_means <- list(lon_hd_mean, lon_hd_SE, lon_hd_count) %>% 
  reduce(full_join, by = "lon")
  lon_hd_binned_means$lon <- as.numeric(as.character(lon_hd_binned_means$lon))
  lon_hd_binned_means <- lon_hd_binned_means[order(lon), ]
  lon_hd_binned_means$X <- lon_hd_binned_means$lon + 5 #for plotting, plot in MIDDLE of 10 degree band
    lon_hd_binned_means$X[lon_hd_binned_means$lon == 180] <- 180 #these two are on the 180 line

CP_lon_hd_binned_means <- merge(CP_lon_hd_mean, CP_lon_hd_SE, by = c("lon", "Pelagic_Coastal"))
CP_lon_hd_binned_means <- merge(CP_lon_hd_binned_means, CP_lon_hd_count, by = c("lon", "Pelagic_Coastal"))
  colnames(CP_lon_hd_binned_means) <- c("lon", "Pelagic_Coastal", "hd_mean", "hd_SE", "hd_count")
  CP_lon_hd_binned_means$lon <- as.numeric(as.character(CP_lon_hd_binned_means$lon))
  CP_lon_hd_binned_means <- CP_lon_hd_binned_means[order(lon), ]
  CP_lon_hd_binned_means$X <- CP_lon_hd_binned_means$lon + 5 #for plotting, plot in MIDDLE of 10 degree band
    CP_lon_hd_binned_means$X[CP_lon_hd_binned_means$lon == 180] <- 180 #these two are on the 180 line

#calculate error bars (standard error)
lon_hd_binned_means$mean_lowerSE <- lon_hd_binned_means$hd_mean - 
  lon_hd_binned_means$hd_SE
lon_hd_binned_means$mean_upperSE <- lon_hd_binned_means$hd_mean + 
  lon_hd_binned_means$hd_SE

CP_lon_hd_binned_means$mean_lowerSE <- CP_lon_hd_binned_means$hd_mean - 
  CP_lon_hd_binned_means$hd_SE
CP_lon_hd_binned_means$mean_upperSE <- CP_lon_hd_binned_means$hd_mean + 
  CP_lon_hd_binned_means$hd_SE

#add count bin for plotting (all < number, 201 = >200)
lon_hd_binned_means$count_bin <- c(50, 50, 100, 50, 50, 50, 50, 50, 50, 100, 100, 100, 
                                   50, 50, 50, 50, 50, 100, 100, 200, 50, 50, 50, 50, 50, 50, 
                                   50, 50, 100, 201, 100, 100, 50, 50, 50, 5)
  lon_hd_binned_means$count_bin <- factor(lon_hd_binned_means$count_bin, levels = c(5, 50, 100, 200, 201))

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
mtdna_hd_lon_plot_both <- ggplot() + 
  geom_line(data = lon_eff_data, 
              aes(x = lon, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = lon_eff_data, 
              aes(x = lon, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = lon_hd_binned_means, 
             aes(x = X, y = hd_mean, color = "10-degree binned means", size = count_bin), shape = "square") + 
  geom_errorbar(data = lon_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("longitude") + ylab("mtdna Hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors) +
  scale_size_manual(values = c(6, 8, 12, 14, 16), 
                    labels = c("\u22645", "\u226450", "\u2264100", "\u2264200", "\u2265200"))
mtdna_hd_lon_plot_annotated_both <- mtdna_hd_lon_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 70, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 60, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 50), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_lon_plot_annotated_both

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

## Build prediction dataframe ##

#get sst values to predict at and scale (range in dataset: 0.182, 32.856)
#predict for SST mean (range is 0.182 to 32.856 C)
SSTmean <- c(.01, 5, 10, 15, 20, 25, 30, 35) #starting at .01 bc log10 of 0 is Inf
  logsstmean <- log10(SSTmean)

#build general prediction dataframe
#have to add a row for each variable in the model
hd_sstmean_predict_df <- as.data.frame(logsstmean)
hd_sstmean_predict_df$bp_scale <- 1.36619e-16 #mean from full dataset
hd_sstmean_predict_df$range_position <- 0.437127 #mean from full dataset
hd_sstmean_predict_df$X <- SSTmean #for plotting

#build prediction dataframe for pelagics
Pel_hd_sstmean_predict_df <- hd_sstmean_predict_df
  Pel_hd_sstmean_predict_df$Pelagic_Coastal <- "Pelagic"

#build prediction dataframe for coastals
Coast_hd_sstmean_predict_df <- hd_sstmean_predict_df
  Coast_hd_sstmean_predict_df$Pelagic_Coastal <- "Coastal"

## Predict for SSTmean ##
#here, use type = "response" to get probability of success (probability of being heterozygote, ~hd)

hd_sstmean_predict_df$predict_hd <- predict(sstmean_model_hd, hd_sstmean_predict_df, type = "response")
Pel_hd_sstmean_predict_df$predict_hd <- predict(sstmean_CP_model_hd, Pel_hd_sstmean_predict_df, type = "response")
Coast_hd_sstmean_predict_df$predict_hd <- predict(sstmean_CP_model_hd, Coast_hd_sstmean_predict_df, type = "response")

## Add confidence intervals ##
#bootMer on predicted dataset (if mixed effects)
#Boot on predicted dataset (if fixed effects only)

hd_sstmean_boot <- Boot(sstmean_model_hd, f = function(x)predict(x, hd_sstmean_predict_df, 
                                                         type = "response"), R = 100)
  hd_sstmean_predict_df$lower_ci <- apply(hd_sstmean_boot$t, 2, quantile, 0.025)
  hd_sstmean_predict_df$upper_ci <- apply(hd_sstmean_boot$t, 2, quantile, 0.975)

Pel_hd_sstmean_boot <- Boot(sstmean_CP_model_hd, f = function(x)predict(x, Pel_hd_sstmean_predict_df, 
                                                                type = "response"), R = 100)
  Pel_hd_sstmean_predict_df$lower_ci <- apply(Pel_hd_sstmean_boot$t, 2, quantile, 0.025)
  Pel_hd_sstmean_predict_df$upper_ci <- apply(Pel_hd_sstmean_boot$t, 2, quantile, 0.975)

Coast_hd_sstmean_boot <- Boot(sstmean_CP_model_hd, f = function(x)predict(x, Coast_hd_sstmean_predict_df, 
                                                                  type = "response"), R = 100)
  Coast_hd_sstmean_predict_df$lower_ci <- apply(Coast_hd_sstmean_boot$t, 2, quantile, 0.025)
  Coast_hd_sstmean_predict_df$upper_ci <- apply(Coast_hd_sstmean_boot$t, 2, quantile, 0.975)


#for CP, bind pelagic and coastal dataframes back together
CP_hd_sstmean_predict_df <- rbind(Pel_hd_sstmean_predict_df, Coast_hd_sstmean_predict_df)

#### Calculate means from raw data ####

## Round SSTmean DOWN to nearest multiple of 5 ##

#create column to fill with the rounded lon
mtdna_small_hd$sstmean_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_hd)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_hd$sstmean_round[i] <- DescTools::RoundTo(mtdna_small_hd$sst.BO_sstmean[i], 
                                                     multiple = 5, FUN = floor)} #rounding DOWN (e.g. SSTmean 5-10 gets value of 5)
}

mtdna_small_hd$sstmean_round <- as.factor(mtdna_small_hd$sstmean_round)

## Calculate means and SE within each 5 degree band ##

mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate mean in each 5 degree band
sstmean_hd_mean <- mtdna_small_hd[, mean(He), by = sstmean_round]
CP_sstmean_hd_mean <- mtdna_small_hd[, mean(He), by = list(sstmean_round, Pelagic_Coastal)]
  colnames(sstmean_hd_mean) <- c("sstmean", "hd_mean")
  colnames(CP_sstmean_hd_mean) <- c("sstmean", "Pelagic_Coastal", "hd_mean")

#calculate SE in each 5 degree band
sstmean_hd_SE <- mtdna_small_hd[, std.error(He), by = sstmean_round]
CP_sstmean_hd_SE <- mtdna_small_hd[, std.error(He), by = list(sstmean_round, Pelagic_Coastal)]
  colnames(sstmean_hd_SE) <- c("sstmean", "hd_SE")
  colnames(CP_sstmean_hd_SE) <- c("sstmean", "Pelagic_Coastal", "hd_SE")

#count observations in each 5 degree band
sstmean_hd_count <- mtdna_small_hd[, .N, by = sstmean_round]
CP_sstmean_hd_count <- mtdna_small_hd[, .N, by = list(sstmean_round, Pelagic_Coastal)]
  colnames(sstmean_hd_count) <- c("sstmean", "hd_count")
  colnames(CP_sstmean_hd_count) <- c("sstmean", "Pelagic_Coastal", "hd_count")

#merge mean and SE dataframes together
sstmean_hd_binned_means <- list(sstmean_hd_mean, sstmean_hd_SE, sstmean_hd_count) %>% 
  reduce(full_join, by = "sstmean")
  sstmean_hd_binned_means$sstmean <- as.numeric(as.character(sstmean_hd_binned_means$sstmean))
  sstmean_hd_binned_means <- sstmean_hd_binned_means[order(sstmean), ]
  sstmean_hd_binned_means$X <- sstmean_hd_binned_means$sstmean + 2.5 #for plotting, plot in MIDDLE of 5 degree band

CP_sstmean_hd_binned_means <- merge(CP_sstmean_hd_mean, CP_sstmean_hd_SE, by = c("sstmean", "Pelagic_Coastal"))
CP_sstmean_hd_binned_means <- merge(CP_sstmean_hd_binned_means, CP_sstmean_hd_count, by = c("sstmean", "Pelagic_Coastal"))
  colnames(CP_sstmean_hd_binned_means) <- c("sstmean", "Pelagic_Coastal", "hd_mean", "hd_SE", "hd_count")
  CP_sstmean_hd_binned_means$sstmean <- as.numeric(as.character(CP_sstmean_hd_binned_means$sstmean))
  CP_sstmean_hd_binned_means <- CP_sstmean_hd_binned_means[order(sstmean), ]
  CP_sstmean_hd_binned_means$X <- CP_sstmean_hd_binned_means$sstmean + 2.5 #for plotting, plot in MIDDLE of 5 degree band

#calculate error bars (standard error)
sstmean_hd_binned_means$mean_lowerSE <- sstmean_hd_binned_means$hd_mean - 
  sstmean_hd_binned_means$hd_SE
sstmean_hd_binned_means$mean_upperSE <- sstmean_hd_binned_means$hd_mean + 
  sstmean_hd_binned_means$hd_SE

CP_sstmean_hd_binned_means$mean_lowerSE <- CP_sstmean_hd_binned_means$hd_mean - 
  CP_sstmean_hd_binned_means$hd_SE
CP_sstmean_hd_binned_means$mean_upperSE <- CP_sstmean_hd_binned_means$hd_mean + 
  CP_sstmean_hd_binned_means$hd_SE

#add count bin for plotting (all < number, 201 = >200)
sstmean_hd_binned_means$count_bin <- c(50, 100, 200, 201, 201, 201, 50)
sstmean_hd_binned_means$count_bin <- factor(sstmean_hd_binned_means$count_bin, levels = c(50, 100, 200, 201))

CP_sstmean_hd_binned_means$count_bin <- c(5, 50, 50, 50, 200, 50, 201, 50, 
                                          201, 50, 201, 100, 50, 5)
CP_sstmean_hd_binned_means$count_bin <- factor(CP_sstmean_hd_binned_means$count_bin, levels = c(5, 50, 100, 200, 201))

#### Plot SSTmean ####

## With pelagic and coastal together ##

#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_hd_sstmean_plot_both <- ggplot() + 
  geom_smooth(data = hd_sstmean_predict_df, 
              aes(x = X, y = predict_hd, color = "Regression"), size = 6, se = FALSE) + 
  geom_ribbon(data = hd_sstmean_predict_df, 
              aes(x = X, ymin = lower_ci, ymax = upper_ci, color = "Regression"), alpha = 0.1) + 
  geom_point(data = sstmean_hd_binned_means, 
             aes(x = X, y = hd_mean, color = "10-degree binned means", size = count_bin), shape = "square") + 
  geom_errorbar(data = sstmean_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("SST mean") + ylab("mtdna Hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors) +
  scale_size_manual(values = c(8, 12, 14, 16), 
                    labels = c("\u226450", "\u2264100", "\u2264200", "\u2265200"))
mtdna_hd_sstmean_plot_annotated_both <- mtdna_hd_sstmean_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 70, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 60, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 50), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_sstmean_plot_annotated_both

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
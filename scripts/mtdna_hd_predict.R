######################## Script for Building Regression Figures for mtDNA Hd  ##############################

#Uses final mtdna hd binomial regression models
#Predicts hd at certain values of predictor variable of interest
#Calculates mean hd of raw data (binned every X units)
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
library(splines)

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)
mtdna_env <- read.csv("output/mtdna_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)

#merge dataframes
mtdna <- merge(mtdna, mtdna_env[, c('X', 'sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin', 
                                    'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                    'chloroA.BO_chlomax', 'chloroA.BO_chlomin')], all.x = TRUE)
mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 
                                  'Half_RangeSize', 'Centroid')], all.x = TRUE)

#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)

######################################################################################################################

######## Clean up dataframe ########

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

#convert lon to radians
mtdna_small_hd$lon_360 <- mtdna_small_hd$lon + 180 #convert (-180,180) to (0,360)
mtdna_small_hd$lon_rad <- (2*pi*mtdna_small_hd$lon_360)/360

#######################################################################################################

######### Abslat figures #######

abslat_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + abslat_scale + 
                           (1|Family/Genus/spp) + (1|Source) + (1|MarkerName), 
                         family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                         control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
abslat_eff <- plot_model(abslat_model_hd, type = "pred",
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
abslat_hd_median <- mtdna_small_hd[, median(He), by = abslat_round]
  colnames(abslat_hd_median) <- c("abslat", "hd_median")
  abslat_hd_median$abslat <- as.numeric(as.character(abslat_hd_median$abslat))

#### Plot abslat ####

#just geom_point
mtdna_hd_abslat_plot_point <- ggplot() +
  geom_jitter(data = mtdna_small_hd, aes(x = as.numeric(as.character(abslat_round)) + 5, y = He), 
              fill = "#3F6DAA", color = "#3F6DAA", size = 8)+
  geom_line(data = abslat_eff_data,
            aes(x = abslat, y = predicted), color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = abslat_eff_data,
              aes(x = abslat, ymin = conf.low, ymax = conf.high), color ="black", alpha = 0.1)+
  scale_y_continuous(limits = c(0, 1))+
  scale_x_continuous(breaks = c(5, 15, 25, 35, 45, 55, 65, 75))+
  xlab("absolute latitude") + ylab("mtdna Hd") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 120),
        axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 120, color = "black"),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "none")
mtdna_hd_abslat_plot_point
  
#violin plot, with or without boxplot within
mtdna_small_hd$abslat_round_char <- as.factor(mtdna_small_hd$abslat_round)

mtdna_hd_abslat_plot_violin <- ggplot() +
  geom_violin(data= mtdna_small_hd, aes(x = abslat_round_char, y = He), 
              fill = "#3F6DAA", color = "#3F6DAA") + #factor, plotted as 0-8 numerically
  #geom_boxplot(data= mtdna_small_hd, aes(x = abslat_round_char, y = He), width = 0.2) +
  geom_point(data = abslat_hd_median, aes(x = (abslat + 10)/10, y = hd_median), 
             color = "darkblue", size = 24) + #to put on same scale as factors, divide by 10, adding 10 because 0 actually = factor of 1 (matches 0-10 round group) 
  geom_line(data = abslat_eff_data,
            aes(x = (abslat + 5)/10, y = predicted), col = "black", alpha = 0.3, linewidth = 10) + #dividing by 10 to put on same factor scale, have to add 5 because violin plots are midpoint of 10 degree bands
  geom_ribbon(data = abslat_eff_data,
              aes(x = (abslat + 5)/10, ymin = conf.low, ymax = conf.high), col = "black", alpha = 0.1)+
  scale_x_discrete(labels = c(5, 15, 25, 35, 45, 55, 65, 75)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("absolute latitude") + ylab("mtdna Hd") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 120),
        axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 120, color = "black"),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "none")
mtdna_hd_abslat_plot_violin

#############################################################################################################

######### Lat figures #######

lat_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + lat_scale + I(lat_scale^2) +
                           (1|Family/Genus/spp) + (1|Source) + (1|MarkerName), 
                         family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                         control = glmerControl(optimizer = "bobyqa"))

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
lat_hd_median <- mtdna_small_hd[, median(He), by = lat_round]
  colnames(lat_hd_median) <- c("lat", "hd_median")
  lat_hd_median$lat <- as.numeric(as.character(lat_hd_median$lat))

#### Plot lat ####

#just geom_point
mtdna_hd_lat_plot_point <- ggplot() +
  geom_jitter(data = mtdna_small_hd, aes(x = as.numeric(as.character(lat_round)) + 5, y = He), 
              fill = "#3F6DAA", color = "#3F6DAA", size = 8)+
  geom_line(data = lat_eff_data,
            aes(x = lat, y = predicted), color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = lat_eff_data,
              aes(x = lat, ymin = conf.low, ymax = conf.high), color ="black", alpha = 0.1)+
  scale_y_continuous(limits = c(0, 1))+
  scale_x_continuous(breaks = c(-75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75))+
  xlab("latitude") + ylab("mtdna Hd") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 120),
        axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 120, color = "black"),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "none")
mtdna_hd_lat_plot_point

#violin plot, with or without boxplot within
mtdna_small_hd$lat_round_char <- as.factor(mtdna_small_hd$lat_round)

mtdna_hd_lat_plot_violin <- ggplot() +
  geom_violin(data= mtdna_small_hd, aes(x = lat_round_char, y = He), 
              fill = "#3F6DAA", color = "#3F6DAA") + #factor, plotted as 0-8 numerically
  #geom_boxplot(data= mtdna_small_hd, aes(x = lat_round_char, y = He), width = 0.2) +
  geom_point(data = lat_hd_median, aes(x = (lat + 90)/10, y = hd_median), 
             color = "darkblue", size = 24) + #to put on same scale as factors, divide by 10, adding 90 because -80 actually = factor of 1 (matches -80 - -70 round group) 
  geom_line(data = lat_eff_data,
            aes(x = (lat + 85)/10, y = predicted), col = "black", alpha = 0.3, linewidth = 10) + #dividing by 10 to put on same factor scale, have to add 85 because violin plots are midpoint of 10 degree bands
  geom_ribbon(data = lat_eff_data,
              aes(x = (lat + 85)/10, ymin = conf.low, ymax = conf.high), col = "black", alpha = 0.1)+
  scale_x_discrete(labels = c(-75, -65, -55, -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("latitude") + ylab("mtdna Hd") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 120),
        axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 120, color = "black"),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "none")
mtdna_hd_lat_plot_violin

#############################################################################################################

######### Lon figures #######

lon_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + sin(lon_rad) + cos(lon_rad) +
                        (1|Family/Genus/spp) + (1|Source) + (1|MarkerName), 
                      family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                      control = glmerControl(optimizer = "bobyqa"))

lon_model_hd_spline <- glmer(cbind(success, failure) ~ bp_scale + range_position + ns(lon) + 
                          (1|Family/Genus/spp) + (1|Source) + (1|MarkerName),
                        family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                        control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
lon_eff <- plot_model(lon_model_test, type = "pred",
                      terms = "lon [all]")

#pull out marginal effects dataframe
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
                                                     multiple = 10, FUN = floor)} #rounding DOWN (e.g. lat 10-20 gets value of 10)
}

  mtdna_small_hd$lon_round <- as.factor(mtdna_small_hd$lon_round)

## Calculate means and SE within each 10 degree band ##
mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate mean in each 10 degree band
lon_hd_median <- mtdna_small_hd[, median(He), by = lon_round]
  colnames(lon_hd_median) <- c("lon", "hd_median")
  lon_hd_median$lon <- as.numeric(as.character(lon_hd_median$lon))
  lon_hd_median <- lon_hd_median[order(lon)]

#### Plot lon ####

#just geom_point
mtdna_hd_lon_plot_point <- ggplot() +
  geom_jitter(data = mtdna_small_hd, aes(x = as.numeric(as.character(lon_round)) + 5, y = He), 
              fill = "#3F6DAA", color = "#3F6DAA", size = 8)+
  geom_line(data = lon_eff_data,
            aes(x = x, y = predicted), color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = lon_eff_data,
              aes(x = x, ymin = conf.low, ymax = conf.high), color ="black", alpha = 0.1)+
  scale_y_continuous(limits = c(0, 1))+
  scale_x_continuous(breaks = c(-175, -165, -155, -145, -135, -125, -115, -105, -95, -85, -75, -65, -55, 
                                -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 
                                125, 135, 145, 155, 165, 175))+
  xlab("longitude") + ylab("mtdna Hd") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 120),
        axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 120, color = "black"),
        axis.line = element_line(linewidth = 4, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
mtdna_hd_lon_plot_point

#violin plot, with or without boxplot within
mtdna_small_hd$lon_round_char <- as.factor(mtdna_small_hd$lon_round)

mtdna_hd_lon_plot_violin <- ggplot() +
  #geom_violin(data= mtdna_small_hd, aes(x = lon_round, y = He), 
   #           fill = "#3F6DAA", color = "#3F6DAA") + #factor, plotted as 0-8 numerically
  #geom_boxplot(data= mtdna_small_hd, aes(x = lon_round_char, y = He), width = 0.2) +
  geom_point(data = lon_hd_median, aes(x = lon, y = hd_median), 
             color = "darkblue", size = 24) + #to put on same scale as factors, divide by 10, adding 190 because -180 actually = factor of 1 (matches -180 - -170 round group) 
  geom_line(data = lon_eff_data,
            aes(x = (x + 185)/10, y = predicted), col = "black", alpha = 0.3, linewidth = 10) + #dividing by 10 to put on same factor scale, have to add 85 because violin plots are midpoint of 10 degree bands
  geom_ribbon(data = lon_eff_data,
              aes(x = (x + 185)/10, ymin = conf.low, ymax = conf.high), col = "black", alpha = 0.1)+
  scale_x_discrete(labels = c(-175, -165, -155, -145, -135, -125, -115, -105, -95, -85, -75, -65, -55, 
                              -45, -35, -25, -15, -5, 5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 
                              125, 135, 145, 155, 165, 175)) +
  scale_y_continuous(limits = c(0.4, 1)) +
  xlab("longitude") + ylab("mtdna Hd") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 120),
        axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 100, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "none")
mtdna_hd_lon_plot_violin

#############################################################################################################

######## Calculate environmental variables ########

#### log transform sst data ####
#subset to only those with sst data
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$sst.BO_sstmean != "NA") #shouldn't remove any

mtdna_small_hd$logsstmean <- log10(mtdna_small_hd$sst.BO_sstmean)
  mtdna_small_hd$logsstrange <- log10(mtdna_small_hd$sst.BO_sstrange)
  mtdna_small_hd$logsstmax <- log10(mtdna_small_hd$sst.BO_sstmax)
  mtdna_small_hd$logsstmin <- log10(mtdna_small_hd$sst.BO_sstmin)

#remove logsst = NA columns
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logsstmean != "NaN")
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logsstrange != "NaN")
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logsstmax != "NaN")
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logsstmin != "NaN")

#### log transform dissox data ####
#subset to only those with dissox data
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$BO_dissox != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

mtdna_small_hd$logdissox <- log10(mtdna_small_hd$BO_dissox)

#remove logdissox = NA columns
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logdissox != "NaN")

#### log transform chlorophyll A ####
#subset to only those with chloroA data
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlorange != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomax != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomin != 0) #if any zeros will screw up log transformation (log10(0) is undefined)

mtdna_small_hd$logchlomean <- log10(mtdna_small_hd$sst.BO_sstmean)
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

##################################################################################################################

######### SST mean figures #######

SSTmean_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + sst.BO_sstmean +  
                        (1|Family/Genus/spp) + (1|Source) + (1|MarkerName), 
                        family = binomial, data = mtdna_small_hd, na.action = "na.fail", 
                        control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
SSTmean_eff <- plot_model(SSTmean_model_hd, type = "pred",
                      terms = "sst.BO_sstmean [all]")

#pull out marginal effects dataframe
SSTmean_eff_data <- as.data.frame(SSTmean_eff$data)

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

#calculate median in each 5 degree band
SSTmean_hd_median <- mtdna_small_hd[, median(He), by = SSTmean_round]
  colnames(SSTmean_hd_median) <- c("SSTmean", "hd_median")
  SSTmean_hd_median$SSTmean <- as.numeric(as.character(SSTmean_hd_median$SSTmean))

#### Plot SSTmean ####

#just geom_point
mtdna_hd_sstmean_plot_point <- ggplot() +
  geom_jitter(data = mtdna_small_hd, aes(x = as.numeric(as.character(SSTmean_round)) + 2.5, y = He), 
              fill = "#3F6DAA", color = "#3F6DAA", size = 8)+
  geom_line(data = SSTmean_eff_data,
            aes(x = x, y = predicted), color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = SSTmean_eff_data,
              aes(x = x, ymin = conf.low, ymax = conf.high), color ="black", alpha = 0.1)+
  scale_y_continuous(limits = c(0, 1))+
  scale_x_continuous(breaks = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5))+
  xlab("SST mean") + ylab("mtdna Hd") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 120),
        axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 120, color = "black"),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "none")
mtdna_hd_sstmean_plot_point
  
#violin plot, with or without boxplot within
mtdna_small_hd$sstmean_round_char <- as.factor(mtdna_small_hd$SSTmean_round)

mtdna_hd_sstmean_plot_violin <- ggplot() +
  geom_violin(data= mtdna_small_hd, aes(x = sstmean_round_char, y = He), 
              fill = "#3F6DAA", color = "#3F6DAA") + #factor, plotted as 0-8 numerically
  #geom_boxplot(data= mtdna_small_hd, aes(x = sstmean_round_char, y = He), width = 0.2) +
  geom_point(data = SSTmean_hd_median, aes(x = (SSTmean + 5)/5, y = hd_median), 
             color = "darkblue", size = 24) + #to put on same scale as factors, divide by 5, adding 5 because 0 actually = factor of 1 (matches 0 - 5 round group) 
  geom_line(data = SSTmean_eff_data,
            aes(x = (x + 2.5)/5, y = predicted), col = "black", alpha = 0.3, linewidth = 10) + #dividing by 10 to put on same factor scale, have to add 85 because violin plots are midpoint of 10 degree bands
  geom_ribbon(data = SSTmean_eff_data,
              aes(x = (x + 2.5)/5, ymin = conf.low, ymax = conf.high), col = "black", alpha = 0.1)+
  scale_x_discrete(labels = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("SST mean") + ylab("mtdna Hd") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 120),
        axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 120, color = "black"),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "none")
mtdna_hd_sstmean_plot_violin

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

chloroAmean_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + chloroA.BO_chlomean + 
                                I(chloroA.BO_chlomean^2) + (1|Family/Genus/spp) + (1|Source) + 
                                (1|MarkerName), family = binomial, 
                                      data = mtdna_small_hd, na.action = "na.fail", 
                                      control = glmerControl(optimizer = "bobyqa"))
#### Predict ####
#marginal effects
chloroAmean_eff <- plot_model(chloroAmean_model_hd, type = "pred", 
                              terms = "chloroA.BO_chlomean [all]")

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
chloroAmean_hd_median <- mtdna_small_hd[, median(He), by = chloroAmean_round]
  colnames(chloroAmean_hd_median) <- c("chloroAmean", "hd_median")
  chloroAmean_hd_median$chloroAmean <- as.numeric(as.character(chloroAmean_hd_median$chloroAmean))

#### Plot chloroAmean ####http://127.0.0.1:31005/graphics/plot_zoom_png?width=1274&height=890

#just geom_point
mtdna_hd_chloromean_plot_point <- ggplot() +
  geom_point(data = mtdna_small_hd, aes(x = chloroA.BO_chlomean, y = He), 
              fill = "#3F6DAA", color = "#3F6DAA", size = 8)+
  geom_line(data = chloroAmean_eff_data,
            aes(x = x, y = predicted), color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = chloroAmean_eff_data,
              aes(x = x, ymin = conf.low, ymax = conf.high), color ="black", alpha = 0.1)+
  scale_y_continuous(limits = c(0, 1))+
  scale_x_continuous(breaks = c(5, 15, 25, 35)) +
  xlab("chlorophyll A mean") + ylab("mtdna Hd") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 120),
        axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 120, color = "black"),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "none")
mtdna_hd_chloromean_plot_point
  
#violin plot, with or without boxplot within
mtdna_small_hd$chloroAmean_round_char <- as.factor(mtdna_small_hd$chloroAmean_round)

mtdna_hd_chloromean_plot_violin <- ggplot() +
  geom_violin(data= mtdna_small_hd, aes(x = chloroAmean_round_char, y = He), 
              fill = "#3F6DAA", color = "#3F6DAA") + #factor, plotted as 0-8 numerically
  #geom_boxplot(data= mtdna_small_hd, aes(x = chloroAmean_round_char, y = He), width = 0.2) +
  geom_point(data = chloroAmean_hd_median, aes(x = (chloroAmean + 10)/10, y = hd_median), 
             color = "darkblue", size = 24) + #to put on same scale as factors, divide by 10, adding 90 because -80 actually = factor of 1 (matches -80 - -70 round group) 
  geom_line(data = chloroAmean_eff_data,
            aes(x = (chlomean + 5)/10, y = predicted), col = "black", alpha = 0.3, linewidth = 10) + #dividing by 10 to put on same factor scale, have to add 85 because violin plots are midpoint of 10 degree bands
  geom_ribbon(data = chloroAmean_eff_data,
              aes(x = (chlomean + 5)/10, ymin = conf.low, ymax = conf.high), col = "black", alpha = 0.1)+
  scale_x_discrete(labels = c(5, 15, 25, 35)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("chlorophyll A mean") + ylab("mtdna Hd") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 120),
        axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 120, color = "black"),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "none")
mtdna_hd_chloromean_plot_violin

#############################################################################################################

######## ChloroA range figures ########

chloroArange_model_hd <- glmer(cbind(success, failure) ~ bp_scale + range_position + logchlorange + 
                                         I(logchlorange^2) + (1|Family/Genus/spp) + (1|Source) + 
                                         (1|MarkerName), family = binomial, 
                                       data = mtdna_small_hd, na.action = "na.fail", 
                                       control = glmerControl(optimizer = "bobyqa"))
#### Predict ####
#marginal effects
chloroArange_eff <- plot_model(chloroArange_model_hd, type = "pred", 
                               terms = "logchlorange [all]")

#pull out marginal effects dataframe
chloroArange_eff_data <- as.data.frame(chloroArange_eff$data)
  chloroArange_eff_data$chlorange <- 10^(chloroArange_eff_data$x)

#### Calculate means from raw data ####
## Round chlomrange DOWN to nearest multiple of 10 ##
#create column to fill with the rounded chlomean
mtdna_small_hd$chloroArange_round <- NA

#fill round column in
for(i in 1:nrow(mtdna_small_hd)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_hd$chloroArange_round[i] <- DescTools::RoundTo(mtdna_small_hd$chloroA.BO_chlorange[i], 
                                                              multiple = 10, FUN = floor)} #rounding DOWN
}

  mtdna_small_hd$chloroArange_round <- as.factor(mtdna_small_hd$chloroArange_round)

## Calculate means and SE within each band ##
mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate mean in each band
chloroArange_hd_mean <- mtdna_small_hd[, mean(He), by = chloroArange_round]
  colnames(chloroArange_hd_mean) <- c("chloroArange", "hd_mean")

#calculate SE in each band
chloroArange_hd_SE <- mtdna_small_hd[, std.error(He), by = chloroArange_round]
colnames(chloroArange_hd_SE) <- c("chloroArange", "hd_SE")

#count observations in each band
chloroArange_hd_count <- mtdna_small_hd[, .N, by = chloroArange_round]
  colnames(chloroArange_hd_count) <- c("chloroArange", "hd_count")

#merge mean and SE dataframes together
chloroArange_hd_binned_means <- list(chloroArange_hd_mean, chloroArange_hd_SE, chloroArange_hd_count) %>% 
    reduce(full_join, by = "chloroArange")
chloroArange_hd_binned_means$chloroArange <- as.numeric(as.character(chloroArange_hd_binned_means$chloroArange))
  chloroArange_hd_binned_means <- chloroArange_hd_binned_means[order(chloroArange), ]
  chloroArange_hd_binned_means$X <- chloroArange_hd_binned_means$chloroArange + 5 #for plotting, plot in MIDDLE of band

#calculate error bars (standard error)
chloroArange_hd_binned_means$mean_lowerSE <- chloroArange_hd_binned_means$hd_mean - 
  chloroArange_hd_binned_means$hd_SE
chloroArange_hd_binned_means$mean_upperSE <- chloroArange_hd_binned_means$hd_mean + 
  chloroArange_hd_binned_means$hd_SE

#### Plot chloroArange ####
#for legend
colors <- c("Binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_hd_chloroArange_plot_both <- ggplot() + 
  geom_line(data = chloroArange_eff_data, 
            aes(x = chlorange, y = predicted, color = "Regression"), size = 6) + 
  geom_ribbon(data = chloroArange_eff_data, 
              aes(x = chlorange, ymin = conf.low, ymax = conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = chloroArange_hd_binned_means, 
             aes(x = X, y = hd_mean, color = "Binned means", size = hd_count), shape = "square") + 
  geom_errorbar(data = chloroArange_hd_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "Binned means"), width = 1.75, size = 3) + 
  xlab("Chloro A range") + ylab("mtDNA Hd") + labs(color = "Legend") + 
  scale_color_manual(values = colors) + 
  scale_y_continuous(limits = c(0, 1.0)) + 
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_hd_chloroArange_plot_annotated_both <- mtdna_hd_chloroArange_plot_both + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_hd_chloroArange_plot_annotated_both
#################################### Script for Building Regression Figures for msat He #############################################

#Uses msat He binomial regression models
#Predicts He at certain values of predictor variable of interest
#Calculates mean or median He of raw data (binned every X units)
#Plots two together
 
##########################################################################################################################################
 
######## Set-up ########
 
remove(list = ls())
 
#load libraries
library(tidyverse)
library(lme4)
library(data.table)
library(sjPlot)
library(splines)
 
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
 
##############################################################################################################################
 
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

################################################################################################################################

######## Range position figure ########

null_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position +  
                         (1|Family/Genus) + (1|Source) + (1|MarkerName), 
                       family = binomial, data = msat, 
                       na.action = "na.fail", nAGQ = 0,
                       control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
rangepos_eff <- plot_model(null_model_he, type = "pred", 
                           terms = "range_position [all]")

#pull out marginal effects dataframe
rangepos_eff_data <- as.data.frame(rangepos_eff$data)

#### Plot range position ####
msat_he_rangepos_plot <- ggplot() +
  geom_line(data = rangepos_eff_data,
            aes(x = x, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = rangepos_eff_data,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              color ="black", alpha = 0.1) + #alpha makes this take way too long, so removing
  geom_rug(data = msat, mapping = aes(x = range_position), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 0.05, y = 0.985, label = "C", size =90) + 
  ylim(0.5, 1) + xlim(0, 1) +
  xlab("Range Position") + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 130, vjust = -1.8),
        axis.title.y = element_text(size = 130, vjust = 5),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 130, color = "black", margin = margin(t = 30)),
        axis.text.y = element_text(size = 130, color = "black", margin = margin(r = 30)),
        axis.line = element_line(linewidth = 4, color = "black"),
        plot.margin = unit(c(1,1.8,3,6), "cm"),
        legend.position = "none")
msat_he_rangepos_plot

###############################################################################################################################

######## CrossSpp figure ########

### Predict ####
#marginal effects
crosspp_eff <- plot_model(null_model_he, type = "pred", 
                          terms = "CrossSpp [all]")

#pull out marginal effects dataframe
crosspp_eff_data <- as.data.frame(crosspp_eff$data)

#### Plot CrossSpp ####
msat_he_crossspp_plot <- ggplot() +
  geom_line(data = crosspp_eff_data,
            aes(x = x, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = crosspp_eff_data,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              color ="black", alpha = 0.1) + 
  geom_rug(data = msat, mapping = aes(x = CrossSpp), 
           color = "#282828", linewidth = 2, inherit.aes = FALSE) + 
  ylim(0.5, 1) + xlim(0, 1) +
  xlab("Cross-Species Primer") + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 130, vjust = -1.8),
        axis.title.y = element_text(size = 130, vjust = 5),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 130, color = "black", margin = margin(t = 30)),
        axis.text.y = element_text(size = 130, color = "black", margin = margin(r = 30)),
        axis.line = element_line(linewidth = 4, color = "black"),
        plot.margin = unit(c(1,1.8,3,6), "cm"),
        legend.position = "none")
msat_he_crossspp_plot

###########################################################################################################################################

######### Abslat figure #######

abslat_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + abslat_scale +
                           (1|Family/Genus) + (1|Source) + (1|ID),
                         family = binomial, data = msat, 
                         na.action = "na.fail", nAGQ = 0,
                         control = glmerControl(optimizer = "bobyqa"))

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

#### Calculate medians from raw data ####
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

## Calculate medians within each 10 degree band ##
msat <- data.table(msat) #make data.table so can use data.table functions

#calculate median in each 10 degree band
abslat_he_median <- msat[, median(He), by = abslat_round]
  colnames(abslat_he_median) <- c("abslat", "he_median")
  abslat_he_median$abslat <- as.numeric(as.character(abslat_he_median$abslat))

#### Plot abslat ####

#violin plot
msat_he_abslat_plot_violin <- ggplot() +
  geom_violin(data= msat, aes(x = abslat_round, y = He),
              fill = "#768e92", color = "#768e92") + #factor, plotted as 0-8 numerically
  geom_point(data = abslat_he_median, aes(x = (abslat + 10)/10, y = he_median),
             color = "#3d4e50", size = 24) + #to put on same scale as factors, divide by 10, adding 10 because 0 actually = factor of 1 (matches 0-10 round group)
  geom_line(data = abslat_eff_data,
            aes(x = (abslat + 5)/10, y = predicted), 
            col = "black", alpha = 0.3, linewidth = 10) + #dividing by 10 to put on same factor scale, have to add 5 because violin plots are midpoint of 10 degree bands
  geom_ribbon(data = abslat_eff_data,
              aes(x = (abslat + 5)/10, ymin = conf.low, ymax = conf.high), 
              col = "black", alpha = 0.1) +
  annotate("text", x = 8, y = 0.97, label = "C", size = 90) +
  scale_x_discrete(labels = c(5, 15, 25, 35, 45, 55, 65, 75)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Absolute Latitude") + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 130, vjust = -1.6),
        axis.title.y = element_text(size = 130, vjust = 5),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 130, color = "black", margin = margin(t = 30)),
        axis.text.y = element_text(size = 130, color = "black", margin = margin(r = 30)),
        axis.line = element_line(linewidth = 4, color = "black"),
        plot.margin = unit(c(1,1.8,3,6), "cm"),
        legend.position = "none")
msat_he_abslat_plot_violin

################################################################################################################################

######### Lat figure #######

lat_model_he <- glmer(cbind(success, failure) ~ CrossSpp + range_position + 
                        lat_scale + I(lat_scale^2) +
                        (1|Family/Genus) + (1|Source) + (1|ID),
                      family = binomial, data = msat, 
                      na.action = "na.fail", nAGQ = 0,
                      control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
lat_eff <- plot_model(lat_model_he, type = "pred",
                      terms = "lat_scale [all]")

#pull out marginal effects dataframe
lat_eff_data <- as.data.frame(lat_eff$data)

#unscale lat
#use same scaled:center & scaled:scale from original data
lat_scale <- scale(msat$lat) #bc had to convert to numeric to run model/calculate marginal effects
lat_eff_data$lat <- (lat_eff_data$x * attr(lat_scale, "scaled:scale")) +
  attr(lat_scale, "scaled:center")

#### Calculate means from raw data ####
## Round lat DOWN to nearest multiple of 10 ##
#create column to fill with the rounded lat
msat$lat_round <- NA

#fill round column in
for(i in 1:nrow(msat)) {
  cat(paste(i, " ", sep = ''))
  {msat$lat_round[i] <- DescTools::RoundTo(msat$lat[i],
                                           multiple = 10, FUN = floor)} #rounding DOWN (e.g. lat 10-20 gets value of 10)
}
  msat$lat_round <- as.factor(msat$lat_round)

## Calculate medians and MADs within each 10 degree band ##
msat <- data.table(msat) #make data.table so can use data.table functions

#calculate median in each 10 degree band
lat_he_median <- msat[, median(He), by = lat_round]
  colnames(lat_he_median) <- c("lat", "he_median")

#calculate MAD in each 10 degree band
lat_he_MAD <- msat[, mad(He), by = lat_round] #median absolute deviation (less affected by outliers, measure of dispersion around median = spread of observations in a dataset)
  colnames(lat_he_MAD) <- c("lat", "he_MAD")

#count observations in each 10 degree band
lat_he_count <- msat[, .N, by = lat_round]
  colnames(lat_he_count) <- c("lat", "he_count")

#merge median and SE dataframes together
lat_he_binned_medians <- list(lat_he_median, lat_he_MAD, lat_he_count) %>%
  reduce(full_join, by = "lat")
lat_he_binned_medians$lat <- as.numeric(as.character(lat_he_binned_medians$lat))
lat_he_binned_medians <- lat_he_binned_medians[order(lat), ]
lat_he_binned_medians$X <- lat_he_binned_medians$lat + 5 #for plotting, plot in MIDDLE of 10 degree band

#calculate error bars (MAD, variance around median)
lat_he_binned_medians$median_lowerMAD <- lat_he_binned_medians$he_median -
  lat_he_binned_medians$he_MAD
lat_he_binned_medians$median_upperMAD <- lat_he_binned_medians$he_median +
  lat_he_binned_medians$he_MAD
  lat_he_binned_medians$median_upperMAD[lat_he_binned_medians$median_upperMAD > 1] <- 1 #bound at 1

#### Plot lat ####

msat_he_lat_plot <- ggplot() +
  geom_point(data = lat_he_binned_medians,
             aes(x = X, y = he_median), 
             color = "darkblue", shape = "circle", size = 24) +
  geom_errorbar(data = lat_he_binned_medians,
                aes(x = X, ymin = median_lowerMAD, ymax = median_upperMAD),
                color = "darkblue", width = 0, linewidth = 3) +
  geom_line(data = lat_eff_data,
            aes(x = lat, y = predicted), color ="black", 
            alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = lat_eff_data,
              aes(x = lat, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1)+
    annotate("text", x = -70, y = 0.97, label = "C", size = 90) +
    scale_x_continuous(breaks = seq(-80, 80, 20)) +
    scale_y_continuous(limits = c(0, 1)) +
    xlab("Latitude") + ylab(bquote(H[e]~"(nucDNA)")) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
          axis.title.x = element_text(size = 160, vjust = -1.6),
          axis.title.y = element_text(size = 160, vjust = 5),
          axis.ticks = element_line(color = "black", linewidth = 2),
          axis.text.x = element_text(size = 160, color = "black", margin = margin(t = 30)),
          axis.text.y = element_text(size = 160, color = "black", margin = margin(r = 30)),
          axis.line = element_line(linewidth = 4, color = "black"),
          plot.margin = unit(c(1,1.8,3,6), "cm"),
          legend.position = "none")
msat_he_lat_plot

##################################################################################################################################

######### Lon figure #######

lon_model_he_spline <- glmer(cbind(success, failure) ~ range_position + CrossSpp + bs(lon_scale) +
                               (1|Family/Genus) + (1|Source) + (1|ID),
                             family = binomial, data = msat, 
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
lon_eff <- plot_model(lon_model_he_spline, type = "eff",
                      terms = "lon_scale [all]")

#pull out marginal effects dataframe
lon_eff_data <- as.data.frame(lon_eff$data)

#unscale lat
#use same scaled:center & scaled:scale from original data
lon_scale <- scale(msat$lon) #bc had to convert to numeric to run model/calculate marginal effects
lon_eff_data$lon <- (lon_eff_data$x * attr(lon_scale, "scaled:scale")) +
  attr(lon_scale, "scaled:center")

#### Calculate medians from raw data ####
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

## Calculate medians and MADs within each 10 degree band ##
msat <- data.table(msat) #make data.table so can use data.table functions

#calculate median in each 10 degree band
lon_he_median <- msat[, median(He), by = lon_round]
  colnames(lon_he_median) <- c("lon", "he_median")

#calculate MAD in each 10 degree band
lon_he_MAD <- msat[, mad(He), by = lon_round] #median absolute deviation (less affected by outliers, measure of dispersion around median = spread of observations in a dataset)
  colnames(lon_he_MAD) <- c("lon", "he_MAD")

#count observations in each 10 degree band
lon_he_count <- msat[, .N, by = lon_round]
  colnames(lon_he_count) <- c("lon", "he_count")

#merge median and SE dataframes together
lon_he_binned_medians <- list(lon_he_median, lon_he_MAD, lon_he_count) %>%
  reduce(full_join, by = "lon")
lon_he_binned_medians$lon <- as.numeric(as.character(lon_he_binned_medians$lon))
lon_he_binned_medians <- lon_he_binned_medians[order(lon), ]
lon_he_binned_medians$X <- lon_he_binned_medians$lon + 5 #for plotting, plot in MIDDLE of 10 degree band

#calculate error bars (MAD, variance around median)
lon_he_binned_medians$median_lowerMAD <- lon_he_binned_medians$he_median -
  lon_he_binned_medians$he_MAD
lon_he_binned_medians$median_upperMAD <- lon_he_binned_medians$he_median +
  lon_he_binned_medians$he_MAD
  lon_he_binned_medians$median_upperMAD[lon_he_binned_medians$median_upperMAD > 1] <- 1 #bound at 1

#### Plot lon ####

msat_he_lon_plot <- ggplot() +
  geom_point(data = lon_he_binned_medians,
             aes(x = X, y = he_median), 
             color = "darkblue", shape = "circle", size = 24) +
  geom_errorbar(data = lon_he_binned_medians,
                aes(x = X, ymin = median_lowerMAD, ymax = median_upperMAD),
                color = "darkblue", width = 0, linewidth = 3) +
  geom_line(data = lon_eff_data,
            aes(x = lon, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = lon_eff_data,
              aes(x = lon, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  annotate("text", x = -175, y = 0.97, label = "C", size = 90) +
  scale_x_continuous(breaks = seq(-180, 180, 20)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Longitude") + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 160, vjust = -1.6),
        axis.title.y = element_text(size = 160, vjust = 5),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 160, color = "black", margin = margin(t = 30), 
                                   angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 160, color = "black", margin = margin(r = 30)),
        axis.line = element_line(linewidth = 4, color = "black"),
        plot.margin = unit(c(1,1.8,3,6), "cm"),
        legend.position = "none")
msat_he_lon_plot

################################################################################################################################
 
######## Calculate environmental variables ########
 
#### log transform chlorophyll ####
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

###################################################################################################################################

######### SST mean figure #######

SSTmean_model_he <- glmer(cbind(success, failure) ~ range_position + CrossSpp + sst.BO_sstmean +
                            (1|Family/Genus) + (1|Source) + (1|ID),
                          family = binomial, data = msat, 
                          na.action = "na.fail", nAGQ = 0,
                          control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
SSTmean_eff <- plot_model(SSTmean_model_he, type = "pred",
                          terms = "sst.BO_sstmean [all]")
  
#pull out marginal effects dataframe
SSTmean_eff_data <- as.data.frame(SSTmean_eff$data)

#### Plot sstmean ####

msat_he_sstmean_plot <- ggplot() +
  geom_line(data = SSTmean_eff_data, aes(x = x, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = SSTmean_eff_data,
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  geom_rug(data = msat, mapping = aes(x = sst.BO_sstmean), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 1, y = 0.842, label = "C", size = 100) + 
  ylim(0.65, 0.85) + xlim(0, 30) +
  xlab("Mean SST (°C)") + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 160, vjust = -2),
        axis.title.y = element_text(size = 160, vjust = 5),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 160, color = "black", margin = margin(t = 30)),
        axis.text.y = element_text(size = 160, color = "black", margin = margin(r = 30)),
        axis.line = element_line(linewidth = 4, color = "black"),
        plot.margin = unit(c(1,1.8,3,6), "cm"),
        legend.position = "none")
msat_he_sstmean_plot

#################################################################################################################################

######### SST range figure #######

SSTrange_model_he <- glmer(cbind(success, failure) ~ range_position + CrossSpp + sst.BO_sstrange +
                            (1|Family/Genus) + (1|Source) + (1|ID),
                          family = binomial, data = msat, 
                          na.action = "na.fail", nAGQ = 0,
                          control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
SSTrange_eff <- plot_model(SSTrange_model_he, type = "pred",
                          terms = "sst.BO_sstrange [all]")

#pull out marginal effects dataframe
SSTrange_eff_data <- as.data.frame(SSTrange_eff$data)

#### Plot sstrange ####

msat_he_sstrange_plot <- ggplot() +
  geom_line(data = SSTrange_eff_data, aes(x = x, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = SSTrange_eff_data,
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  geom_rug(data = msat, mapping = aes(x = sst.BO_sstrange), 
           color = "#282828",inherit.aes = FALSE) + 
  annotate("text", x = 1, y = 0.985, label = "C", size = 100) + 
  ylim(c(0.5, 1.0)) + xlim(0, 30) +
  xlab("SST Range (°C)") + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 160, vjust = -2),
        axis.title.y = element_text(size = 160, vjust = 5),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 160, color = "black", margin = margin(t = 30)),
        axis.text.y = element_text(size = 160, color = "black", margin = margin(r = 30)),
        axis.line = element_line(linewidth = 4, color = "black"),
        plot.margin = unit(c(1,1.8,3,6), "cm"),
        legend.position = "none")
msat_he_sstrange_plot

###############################################################################################################################

######## SST max figure #######

SSTmax_model_he <- glmer(cbind(success, failure) ~ range_position + CrossSpp + sst.BO_sstmax +
                           (1|Family/Genus) + (1|Source) + (1|ID),
                         family = binomial, data = msat, 
                         na.action = "na.fail", nAGQ = 0,
                        control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
SSTmax_eff <- plot_model(SSTmax_model_he, type = "pred",
                           terms = "sst.BO_sstmax [all]")

#pull out marginal effects dataframe
SSTmax_eff_data <- as.data.frame(SSTmax_eff$data)

#### Plot sstmax ####

msat_he_sstmax_plot <- ggplot() +
  geom_line(data = SSTmax_eff_data, aes(x = x, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = SSTmax_eff_data,
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  geom_rug(data = msat, mapping = aes(x = sst.BO_sstmax), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 1, y = 0.985, label = "F", size = 100) + 
  ylim(c(0.5, 1.0)) + xlim(0, 30) +
  xlab("Maximum SST (°C)") + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 160, vjust = -2),
        axis.title.y = element_text(size = 160, vjust = 5),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 160, color = "black", margin = margin(t = 30)),
        axis.text.y = element_text(size = 160, color = "black", margin = margin(r = 30)),
        axis.line = element_line(linewidth = 4, color = "black"),
        plot.margin = unit(c(1,1.8,3,6), "cm"),
        legend.position = "none")
msat_he_sstmax_plot

##################################################################################################################################

######## SST min figure #######

SSTmin_model_he <- glmer(cbind(success, failure) ~ range_position + CrossSpp + sst.BO_sstmin +
                           (1|Family/Genus) + (1|Source) + (1|ID),
                         family = binomial, data = msat, 
                         na.action = "na.fail", nAGQ = 0,
                         control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
SSTmin_eff <- plot_model(SSTmin_model_he, type = "pred",
                           terms = "sst.BO_sstmin [all]")

#pull out marginal effects dataframe
SSTmin_eff_data <- as.data.frame(SSTmin_eff$data)

#### Plot sstmin ####

msat_he_sstmin_plot <- ggplot() +
  geom_line(data = SSTmin_eff_data, aes(x = x, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = SSTmin_eff_data,
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  geom_rug(data = msat, mapping = aes(x = sst.BO_sstmin), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 1, y = 0.985, label = "I", size = 100) + 
  ylim(c(0.5, 1.0)) + xlim(0, 30) +
  xlab("Minimum SST (°C)") + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 160, vjust = -2),
        axis.title.y = element_text(size = 160, vjust = 5),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 160, color = "black", margin = margin(t = 30)),
        axis.text.y = element_text(size = 160, color = "black", margin = margin(r = 30)),
        axis.line = element_line(linewidth = 4, color = "black"),
        plot.margin = unit(c(1,1.8,3,6), "cm"),
        legend.position = "none")
msat_he_sstmin_plot

############################################################################################################################################

######## ChloroA mean figures ########

chloroAmean_model_he <- glmer(cbind(success, failure) ~ range_position + CrossSpp + 
                                logchlomean + I(logchlomean^2) +
                                (1|Family/Genus) + (1|Source) + (1|ID),
                              family = binomial, data = msat, 
                              na.action = "na.fail", nAGQ = 0,
                              control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
chloroAmean_eff <- plot_model(chloroAmean_model_he, type = "pred", 
                             terms = "logchlomean [all]")

#pull out marginal effects dataframe
chloroAmean_eff_data <- as.data.frame(chloroAmean_eff$data)
  chloroAmean_eff_data$chlomean <- 10^(chloroAmean_eff_data$x)

#### Plot chloroAmean ####

msat_he_chloromean_plot <- ggplot() +
  geom_line(data = chloroAmean_eff_data,
            aes(x = chlomean, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = chloroAmean_eff_data,
              aes(x = chlomean, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  geom_rug(data = msat, mapping = aes(x = chloroA.BO_chlomean), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 0.12, y = 0.842, label = "F", size = 100) + 
    ylim(0.65, 0.85) +
    scale_x_continuous(trans = "log10", limits = c(0.1, 10)) +
    xlab(bquote("Mean Chlorophyll"~(mg/m^3))) + ylab(bquote(H[e]~"(nucDNA)")) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
          axis.title.x = element_text(size = 160, vjust = -1.6),
          axis.title.y = element_text(size = 160, vjust = 5),
          axis.ticks = element_line(color = "black", linewidth = 2),
          axis.text.x = element_text(size = 160, color = "black", margin = margin(t = 30)),
          axis.text.y = element_text(size = 160, color = "black", margin = margin(r = 30)),
          axis.line = element_line(linewidth = 4, color = "black"),
          plot.margin = unit(c(1,1.8,3,6), "cm"),
          legend.position = "none")
msat_he_chloromean_plot

##########################################################################################################################################

######## ChloroA range figures ########

chloroArange_model_he <- glmer(cbind(success, failure) ~ range_position + CrossSpp + logchlorange + 
                                 I(logchlorange^2) +(1|Family/Genus) + (1|Source) + (1|ID),
                              family = binomial, data = msat,
                              na.action = "na.fail", nAGQ = 0,
                              control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
chloroArange_eff <- plot_model(chloroArange_model_he, type = "pred", 
                              terms = "logchlorange [all]")

#pull out marginal effects dataframe
chloroArange_eff_data <- as.data.frame(chloroArange_eff$data)
  chloroArange_eff_data$chlorange <- 10^(chloroArange_eff_data$x)

#### Plot chloroArange ####

msat_he_chlororange_plot <- ggplot() +
  geom_line(data = chloroArange_eff_data,
            aes(x = chlorange, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = chloroArange_eff_data,
              aes(x = chlorange, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  geom_rug(data = msat, mapping = aes(x = chloroA.BO_chlorange), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 0.12, y = 0.985, label = "C", size = 100) + 
  ylim(0.5, 1.0) +
  scale_x_continuous(trans = "log10", limits = c(0.1, 10)) +
  xlab(bquote("Chlorophyll Range"~(mg/m^3))) + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 160, vjust = -1.6),
        axis.title.y = element_text(size = 160, vjust = 5),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 160, color = "black", margin = margin(t = 30)),
        axis.text.y = element_text(size = 160, color = "black", margin = margin(r = 30)),
        axis.line = element_line(linewidth = 4, color = "black"),
        plot.margin = unit(c(1,1.8,3,6), "cm"),
        legend.position = "none")
msat_he_chlororange_plot

##############################################################################################################################################

######## ChloroA max figure ########

chloroAmax_model_he <- glmer(cbind(success, failure) ~ range_position + CrossSpp + logchlomax + 
                               I(logchlomax^2) + (1|Family/Genus) + (1|Source) + (1|ID),
                             family = binomial, data = msat, 
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
chloroAmax_eff <- plot_model(chloroAmax_model_he, type = "pred", 
                               terms = "logchlomax [all]")

#pull out marginal effects dataframe
chloroAmax_eff_data <- as.data.frame(chloroAmax_eff$data)
  chloroAmax_eff_data$chlomax <- 10^(chloroAmax_eff_data$x)

#### Plot chloroAmax ####

msat_he_chloromax_plot <- ggplot() +
  geom_line(data = chloroAmax_eff_data,
            aes(x = chlomax, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = chloroAmax_eff_data,
              aes(x = chlomax, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  geom_rug(data = msat, mapping = aes(x = chloroA.BO_chlomax), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 0.12, y = 0.985, label = "F", size = 100) + 
  ylim(0.5, 1.0) +
  scale_x_continuous(trans = "log10", limits = c(0.1, 10)) +
  xlab(bquote("Maximum Chlorophyll"~(mg/m^3))) + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 160, vjust = -1.6),
        axis.title.y = element_text(size = 160, vjust = 5),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 160, color = "black", margin = margin(t = 30)),
        axis.text.y = element_text(size = 160, color = "black", margin = margin(r = 30)),
        axis.line = element_line(linewidth = 4, color = "black"),
        plot.margin = unit(c(1,1.8,3,6), "cm"),
        legend.position = "none")
msat_he_chloromax_plot

############################################################################################################################################

######## ChloroA min figures ########

chloroAmin_model_he <- glmer(cbind(success, failure) ~ range_position + CrossSpp + logchlomin + 
                               I(logchlomin^2) + (1|Family/Genus) + (1|Source) + (1|ID),
                             family = binomial, data = msat, 
                             na.action = "na.fail", nAGQ = 0,
                             control = glmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
chloroAmin_eff <- plot_model(chloroAmin_model_he, type = "pred", 
                               terms = "logchlomin [all]")

#pull out marginal effects dataframe
chloroAmin_eff_data <- as.data.frame(chloroAmin_eff$data)
  chloroAmin_eff_data$chlomin <- 10^(chloroAmin_eff_data$x)

#### Plot chloroAmin ####

msat_he_chloromin_plot <- ggplot() +
    geom_line(data = chloroAmin_eff_data,
              aes(x = chlomin, y = predicted), 
              color ="black", alpha = 0.3, linewidth = 10) +
    geom_ribbon(data = chloroAmin_eff_data,
                aes(x = chlomin, ymin = conf.low, ymax = conf.high), 
                color ="black", alpha = 0.1) +
    geom_rug(data = msat, mapping = aes(x = chloroA.BO_chlomin), 
             color = "#282828", inherit.aes = FALSE) + 
    annotate("text", x = 0.12, y = 0.985, label = "I", size = 100) + 
    ylim(0.5, 1.0) +
    scale_x_continuous(trans = "log10", limits = c(0.1, 10)) +
    xlab(bquote("Minimum Chlorophyll"~(mg/m^3))) + ylab(bquote(H[e]~"(nucDNA)")) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
          axis.title.x = element_text(size = 160, vjust = -1.6),
          axis.title.y = element_text(size = 160, vjust = 5),
          axis.ticks = element_line(color = "black", linewidth = 2),
          axis.text.x = element_text(size = 160, color = "black", margin = margin(t = 30)),
          axis.text.y = element_text(size = 160, color = "black", margin = margin(r = 30)),
          axis.line = element_line(linewidth = 4, color = "black"),
          plot.margin = unit(c(1,1.8,3,6), "cm"),
          legend.position = "none")
msat_he_chloromin_plot
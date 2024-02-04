#################################### Script for Building Regression Figures for msat He #############################################

#Uses msat He beta regression models
#Predicts He at certain values of predictor variable of interest
#Calculates mean or median He of raw data (binned every X units)
#Plots two together
 
##########################################################################################################################################
 
######## Set-up ########
 
remove(list = ls())
 
#load libraries
library(tidyverse) #v.2.0.0
library(glmmTMB) #1.1.7
library(DHARMa) #v.0.4.6
library(sjPlot) #v.2.8.12
library(splines) #v.4.2.2
library(data.table) #1.14.8
 
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

################################################################################################################################

######## Range position figure ########

null_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                           (1|Family/Genus) + (1|Source), 
                         data = msat, family = ordbeta, 
                         na.action = "na.fail")

#### Predict ####
#marginal effects
rangepos_eff <- plot_model(null_model_he, type = "pred", 
                           terms = "range_pos_scale [all]",
                           allow.new.levels = TRUE)

#pull out marginal effects dataframe
rangepos_eff_data <- as.data.frame(rangepos_eff$data)

#unscale range_position
#use same scaled:center & scaled:scale from original data
range_pos_scale <- scale(msat$range_position) #bc had to convert to numeric to run model/calculate marginal effects
rangepos_eff_data$range_position <- (rangepos_eff_data$x * attr(range_pos_scale, "scaled:scale")) +
  attr(range_pos_scale, "scaled:center")

#### Plot range position ####
msat_he_rangepos_plot <- ggplot() +
  geom_line(data = rangepos_eff_data,
            aes(x = range_position, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = rangepos_eff_data,
              aes(x = range_position, ymin = conf.low, ymax = conf.high),
              color ="black", alpha = 0.1) + #alpha makes this take way too long, so removing
  geom_rug(data = msat, mapping = aes(x = range_position), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 0.05, y = 0.985, label = "(c)", size =90) + 
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
                          terms = "CrossSpp_scale [all]",
                          allow.new.levels = TRUE)

#pull out marginal effects dataframe
crosspp_eff_data <- as.data.frame(crosspp_eff$data)

#unscale CrossSpp
#use same scaled:center & scaled:scale from original data
CrossSpp_scale <- scale(msat$CrossSpp) #bc had to convert to numeric to run model/calculate marginal effects
crosspp_eff_data$CrossSpp <- (crosspp_eff_data$x * attr(CrossSpp_scale, "scaled:scale")) +
  attr(CrossSpp_scale, "scaled:center")

#### Plot CrossSpp ####
msat_he_crossspp_plot <- ggplot() +
  geom_line(data = crosspp_eff_data,
            aes(x = CrossSpp, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = crosspp_eff_data,
              aes(x = CrossSpp, ymin = conf.low, ymax = conf.high),
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

abslat_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + abslat_scale +
                             (1|Family/Genus) + (1|Source) + 
                             (0 + abslat_scale|Family), 
                           data = msat, family = ordbeta, 
                           na.action = "na.fail")
#### Predict ####
#marginal effects
abslat_eff <- plot_model(abslat_model_he, type = "pred",
                         terms = "abslat_scale [all]",
                         allow.new.levels = TRUE)

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
  annotate("text", x = 7.75, y = 0.97, label = "(c)", size = 90) +
  scale_x_discrete(labels = c(5, 15, 25, 35, 45, 55, 65, 75)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Absolute Latitude") + ylab(bquote(H[e]~"(nucDNA)")) + 
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
msat_he_abslat_plot_violin

################################################################################################################################

######### Lat figure #######

lat_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                          lat_scale + I(lat_scale^2) +
                          (1|Family/Genus) + (1|Source) + 
                          (0 + lat_scale|Family), 
                        data = msat, family = ordbeta, 
                        na.action = "na.fail")
#### Predict ####
#marginal effects
lat_eff <- plot_model(lat_model_he, type = "pred",
                      terms = "lat_scale [all]",
                      allow.new.levels = TRUE)

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
    annotate("text", x = -65, y = 0.97, label = "(c)", size = 90) +
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

lon_model_he_spline <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + bs(lon_scale) +
                                 (1|Family/Genus) + (1|Source) + 
                                 (0 + lon_scale|Family), 
                               data = msat, family = ordbeta, 
                               na.action = "na.fail")

#### Predict ####
#marginal effects
lon_eff <- plot_model(lon_model_he_spline, type = "pred",
                      terms = "lon_scale [all]",
                      allow.new.levels = TRUE)

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
  annotate("rect", xmin = 95, xmax = 165, ymin = 0, ymax = 1, #adding box highlighting coral triangle
           fill = "darkolivegreen", alpha = 0.4) + 
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
  annotate("text", x = -175, y = 0.97, label = "(c)", size = 90) +
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
 
#scale SST variables
msat$sstmean_scale <- as.numeric(scale(msat$sst.BO_sstmean))

#### log transform chlorophyll ####
#subset to only those with chloroA data
msat <- subset(msat, msat$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat$logchlomean <- log10(msat$chloroA.BO_chlomean)

#remove logchlo = NA columns
msat <- subset(msat, msat$logchlomean != "Inf" |
                 msat$logchlomean != "NaN")

###################################################################################################################################

######### SST mean figure #######

SSTmean_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + sstmean_scale + 
                              (1|Family/Genus) + (1|Source) + 
                              (0 + sstmean_scale|Family),
                            data = msat, family = ordbeta, 
                            na.action = "na.fail") 

#### Predict ####
#marginal effects
SSTmean_eff <- plot_model(SSTmean_model_he, type = "pred",
                          terms = "sstmean_scale [all]",
                          allow.new.levels = TRUE)
  
#pull out marginal effects dataframe
SSTmean_eff_data <- as.data.frame(SSTmean_eff$data)

#unscale SSTmean
#use same scaled:center & scaled:scale from original data
sstmean_scale <- scale(msat$sst.BO_sstmean) #bc had to convert to numeric to run model/calculate marginal effects
SSTmean_eff_data$SSTmean <- (SSTmean_eff_data$x * attr(sstmean_scale, "scaled:scale")) +
  attr(sstmean_scale, "scaled:center")

#### Plot sstmean ####

msat_he_sstmean_plot <- ggplot() +
  geom_line(data = SSTmean_eff_data, aes(x = SSTmean, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = SSTmean_eff_data,
              aes(x = SSTmean, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  geom_rug(data = msat, mapping = aes(x = sst.BO_sstmean), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 2, y = 0.985, label = "(c)", size = 100) + 
  ylim(0.5, 1.0) + xlim(0, 30) +
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

############################################################################################################################################

######## ChloroA mean figures ########

chloroAmean_model_he <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                  logchlomean + I(logchlomean^2) + 
                                  (1|Family/Genus) + (1|Source) + 
                                  (0 + logchlomean|Family),
                                data = msat, family = ordbeta, 
                                na.action = "na.fail")

#### Predict ####
#marginal effects
chloroAmean_eff <- plot_model(chloroAmean_model_he, type = "pred", 
                              terms = "logchlomean [all]",
                              allow.new.levels = TRUE)

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
  annotate("text", x = 0.135, y = 0.985, label = "(f)", size = 100) + 
    ylim(0.5, 1.0) +
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

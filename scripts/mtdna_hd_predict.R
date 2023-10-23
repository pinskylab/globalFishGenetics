######################################## Script for Building Regression Figures for mtDNA Hd  ########################################

#Uses mtDNA Hd beta regression models
#Predicts Hd at certain values of predictor variable of interest
#Calculates mean Hd of raw data (binned every X units)
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

#### Calculate nuisance variables ####
#scale bp
mtdna_small_hd$bp_scale <- as.numeric(scale(as.numeric(mtdna_small_hd$bp)))

#### Add range position ####
#fix character type
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

#scale range_position
mtdna_small_hd$range_pos_scale <- as.numeric(scale(mtdna_small_hd$range_position))

#### Calculate latitude, longitude variables ####
#calculate abslat
mtdna_small_hd$abslat <- abs(mtdna_small_hd$lat)

#scale geographic variables
mtdna_small_hd$lat_scale <- as.numeric(scale(mtdna_small_hd$lat))
mtdna_small_hd$abslat_scale <- as.numeric(scale(mtdna_small_hd$abslat))
mtdna_small_hd$lon_scale <- as.numeric(scale(mtdna_small_hd$lon))

#############################################################################################################

######## Range position figure ########

null_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + 
                           (1|Family/Genus) + (1|Source) + (1|MarkerName),
                         data = mtdna_small_hd, family = ordbeta,
                         na.action = "na.fail") 

#### Predict ####
#marginal effects
rangepos_eff <- plot_model(null_model_hd, type = "pred", 
                           terms = "range_pos_scale [all]",
                           allow.new.levels = TRUE)

#pull out marginal effects dataframe
rangepos_eff_data <- as.data.frame(rangepos_eff$data)

#unscale raange_position
#use same scaled:center & scaled:scale from original data
range_pos_scale <- scale(as.numeric(mtdna_small_hd$range_position)) #bc had to convert to numeric to run model/calculate marginal effects
rangepos_eff_data$range_position <- (rangepos_eff_data$x * attr(range_pos_scale, "scaled:scale")) + 
  attr(range_pos_scale, "scaled:center")

#### Plot range position ####

mtdna_hd_rangepos_plot <- ggplot() +
  geom_line(data = rangepos_eff_data,
            aes(x = range_position, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = rangepos_eff_data,
              aes(x = range_position, ymin = conf.low, ymax = conf.high),
              color ="black", alpha = 0.1) +
  geom_rug(data = mtdna_small_hd, mapping = aes(x = range_position), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 0.05, y = 0.985, label = "B", size =90) + 
  ylim(0.5, 1) + xlim(0, 1) +
  xlab("Range Position") + ylab(bquote(H[d]~"(mtDNA)")) + 
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
mtdna_hd_rangepos_plot

###########################################################################################################

######## bp figure ########

#### Predict ####
#marginal effects
bp_eff <- plot_model(null_model_hd, type = "pred", 
                     terms = "bp_scale [all]",
                     allow.new.levels = TRUE)

#pull out marginal effects dataframe
bp_eff_data <- as.data.frame(bp_eff$data)

#unscale bp
#use same scaled:center & scaled:scale from original data
bp_scale <- scale(as.numeric(mtdna_small_hd$bp)) #bc had to convert to numeric to run model/calculate marginal effects
bp_eff_data$bp <- (bp_eff_data$x * attr(bp_scale, "scaled:scale")) + 
  attr(bp_scale, "scaled:center")

#### Plot bp position ####

mtdna_hd_bp_plot <- ggplot() +
  geom_line(data = bp_eff_data,
            aes(x = bp, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = bp_eff_data,
              aes(x = bp, ymin = conf.low, ymax = conf.high),
              color ="black", alpha = 0.1) +
  geom_rug(data = mtdna_small_hd, mapping = aes(x = as.numeric(bp)), 
           color = "#282828", inherit.aes = FALSE) + 
  xlab("Base Pairs") + ylab(bquote(H[d]~"(mtDNA)")) + 
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
mtdna_hd_bp_plot

##########################################################################################################################

######### Abslat figure #######

abslat_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + 
                             abslat_scale + (1|Family/Genus) + 
                             (1|Source) + (1|MarkerName),
                           data = mtdna_small_hd, family = ordbeta, 
                           na.action = "na.fail") 

#### Predict ####
#marginal effects
abslat_eff <- plot_model(abslat_model_hd, type = "pred",
                      terms = "abslat_scale [all]",
                      allow.new.levels = TRUE)

#pull out marginal effects dataframe
abslat_eff_data <- as.data.frame(abslat_eff$data)

#unscale abslat
#use same scaled:center & scaled:scale from original data
abslat_scale <- scale(mtdna_small_hd$abslat) #bc had to convert to numeric to run model/calculate marginal effects
abslat_eff_data$abslat <- (abslat_eff_data$x * attr(abslat_scale, "scaled:scale")) + 
  attr(abslat_scale, "scaled:center")

#### Calculate medians from raw data ####
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

## Calculate medians within each 10 degree band ##
mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate median in each 10 degree band
abslat_hd_median <- mtdna_small_hd[, median(He), by = abslat_round]
  colnames(abslat_hd_median) <- c("abslat", "hd_median")
  abslat_hd_median$abslat <- as.numeric(as.character(abslat_hd_median$abslat))

#### Plot abslat ####
  
#violin plot
mtdna_hd_abslat_plot_violin <- ggplot() +
  geom_violin(data= mtdna_small_hd, aes(x = abslat_round, y = He), 
              fill = "#768e92", color = "#768e92") + #factor, plotted as 0-8 numerically
  geom_point(data = abslat_hd_median, aes(x = (abslat + 10)/10, y = hd_median), 
             color = "#3d4e50", size = 24) + #to put on same scale as factors, divide by 10, adding 10 because 0 actually = factor of 1 (matches 0-10 round group) 
  geom_line(data = abslat_eff_data,
            aes(x = (abslat + 5)/10, y = predicted), 
            col = "black", alpha = 0.3, linewidth = 10) + #dividing by 10 to put on same factor scale, have to add 5 because violin plots are midpoint of 10 degree bands
  geom_ribbon(data = abslat_eff_data,
              aes(x = (abslat + 5)/10, ymin = conf.low, ymax = conf.high), 
              col = "black", alpha = 0.1) +
  annotate("text", x = 8, y = 0.97, label = "B", size = 90) +
  scale_x_discrete(labels = c(5, 15, 25, 35, 45, 55, 65, 75)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Absolute Latitude") + ylab(bquote(H[d]~"(mtDNA)")) + 
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
mtdna_hd_abslat_plot_violin

#############################################################################################################

######### Lat figure #######

lat_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + 
                          lat_scale + I(lat_scale^2) + (1|Family/Genus) + 
                          (1|Source) + (1|MarkerName),
                        data = mtdna_small_hd, family = ordbeta, 
                        na.action = "na.fail")  

#### Predict ####
#marginal effects
lat_eff <- plot_model(lat_model_hd, type = "pred",
                      terms = "lat_scale [all]",
                      allow.new.levels = TRUE)

#pull out marginal effects dataframe
lat_eff_data <- as.data.frame(lat_eff$data)

#unscale lat
#use same scaled:center & scaled:scale from original data
lat_scale <- scale(mtdna_small_hd$lat) #bc had to convert to numeric to run model/calculate marginal effects
lat_eff_data$lat <- (lat_eff_data$x * attr(lat_scale, "scaled:scale")) + 
  attr(lat_scale, "scaled:center")

#### Calculate medians from raw data ####
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

## Calculate medians and MADs within each 10 degree band ##
mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions
  
#calculate median in each 10 degree band
lat_hd_median <- mtdna_small_hd[, median(He), by = lat_round]
  colnames(lat_hd_median) <- c("lat", "hd_median")
  
#calculate MAD in each 10 degree band
lat_hd_MAD <- mtdna_small_hd[, mad(He), by = lat_round] #median absolute deviation (less affected by outliers, measure of dispersion around median = spread of observations in a dataset)
  colnames(lat_hd_MAD) <- c("lat", "hd_MAD")
  
#count observations in each 10 degree band
lat_hd_count <- mtdna_small_hd[, .N, by = lat_round]
  colnames(lat_hd_count) <- c("lat", "hd_count")
  
#merge mean and SE dataframes together
lat_hd_binned_medians <- list(lat_hd_median, lat_hd_MAD, lat_hd_count) %>% 
  reduce(full_join, by = "lat")
  lat_hd_binned_medians$lat <- as.numeric(as.character(lat_hd_binned_medians$lat))
  lat_hd_binned_medians <- lat_hd_binned_medians[order(lat), ]
  lat_hd_binned_medians$X <- lat_hd_binned_medians$lat + 5 #for plotting, plot in MIDDLE of 10 degree band
  
#calculate error bars (MAD, variance around median)
lat_hd_binned_medians$median_lowerMAD <- lat_hd_binned_medians$hd_median - 
    lat_hd_binned_medians$hd_MAD
lat_hd_binned_medians$median_upperMAD <- lat_hd_binned_medians$hd_median + 
    lat_hd_binned_medians$hd_MAD
  lat_hd_binned_medians$median_upperMAD[lat_hd_binned_medians$median_upperMAD > 1] <- 1 #bound at 1
  
#### Plot lat ####

mtdna_hd_lat_plot <- ggplot() +
  geom_point(data = lat_hd_binned_medians, 
             aes(x = X, y = hd_median), color = "darkblue", shape = "circle", size = 24) + 
  geom_errorbar(data = lat_hd_binned_medians, 
                aes(x = X, ymin = median_lowerMAD, ymax = median_upperMAD), 
                color = "darkblue", width = 0, linewidth = 3) + 
  geom_line(data = lat_eff_data,
            aes(x = lat, y = predicted), color ="black", 
            alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = lat_eff_data,
              aes(x = lat, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  annotate("text", x = -70, y = 0.97, label = "B", size = 90) +
  scale_x_continuous(breaks = seq(-80, 80, 20)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Latitude") + ylab(bquote(H[d]~"(mtDNA)")) + 
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
mtdna_hd_lat_plot

#############################################################################################################

######### Lon figure #######

lon_model_hd_spline <- glmmTMB(He ~ bp_scale + range_pos_scale + 
                                 bs(lon_scale) + (1|Family/Genus) + 
                                 (1|Source) + (1|MarkerName),
                               data = mtdna_small_hd, family = ordbeta, 
                               na.action = "na.fail")  

#### Predict ####
#marginal effects
lon_eff <- plot_model(lon_model_hd_spline, type = "pred",
                      terms = "lon_scale [all]",
                      allow.new.levels = TRUE)

#pull out marginal effects dataframe
lon_eff_data <- as.data.frame(lon_eff$data)

#unscale lon
#use same scaled:center & scaled:scale from original data
lon_scale <- scale(mtdna_small_hd$lon) #bc had to convert to numeric to run model/calculate marginal effects
lon_eff_data$lon <- (lon_eff_data$x * attr(lon_scale, "scaled:scale")) + 
  attr(lon_scale, "scaled:center")

#### Calculate medians from raw data ####
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

## Calculate medians and MADs within each 10 degree band ##
mtdna_small_hd <- data.table(mtdna_small_hd) #make data.table so can use data.table functions

#calculate median in each 10 degree band
lon_hd_median <- mtdna_small_hd[, median(He), by = lon_round]
  colnames(lon_hd_median) <- c("lon", "hd_median")
  
#calculate MAD in each 10 degree band
lon_hd_MAD <- mtdna_small_hd[, mad(He), by = lon_round] #median absolute deviation (less affected by outliers, measure of dispersion around median = spread of observations in a dataset)
  colnames(lon_hd_MAD) <- c("lon", "hd_MAD")
  
#count observations in each 10 degree band
lon_hd_count <- mtdna_small_hd[, .N, by = lon_round]
  colnames(lon_hd_count) <- c("lon", "hd_count")
  
#merge median and MAD dataframes together
lon_hd_binned_medians <- list(lon_hd_median, lon_hd_MAD, lon_hd_count) %>% 
  reduce(full_join, by = "lon")
  lon_hd_binned_medians$lon <- as.numeric(as.character(lon_hd_binned_medians$lon))
  lon_hd_binned_medians <- lon_hd_binned_medians[order(lon), ]
lon_hd_binned_medians$X <- lon_hd_binned_medians$lon + 5 #for plotting, plot in MIDDLE of 10 degree band
lon_hd_binned_medians$lon[lon_hd_binned_medians$lon == 185] <- 180 #bc nothing above 180
  
#calculate error bars (MAD, variance around median)
lon_hd_binned_medians$median_lowerMAD <- lon_hd_binned_medians$hd_median - 
  lon_hd_binned_medians$hd_MAD
lon_hd_binned_medians$median_upperMAD <- lon_hd_binned_medians$hd_median + 
  lon_hd_binned_medians$hd_MAD
lon_hd_binned_medians$median_upperMAD[lon_hd_binned_medians$median_upperMAD > 1] <- 1 #bound at 1

#### Plot lon ####

mtdna_hd_lon_plot <- ggplot() +
  annotate("rect", xmin = 95, xmax = 165, ymin = 0, ymax = 1, #adding box highlighting coral triangle
           fill = "darkolivegreen", alpha = 0.4) + 
  geom_point(data = lon_hd_binned_medians, 
             aes(x = X, y = hd_median), 
             color = "darkblue", shape = "circle", size = 24) + 
  geom_errorbar(data = lon_hd_binned_medians, 
                aes(x = X, ymin = median_lowerMAD, ymax = median_upperMAD), 
                color = "darkblue", width = 0, linewidth = 3) + 
  geom_line(data = lon_eff_data,
            aes(x = lon, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = lon_eff_data,
              aes(x = lon, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  annotate("text", x = -175, y = 0.97, label = "B", size = 100) +
  scale_x_continuous(breaks = seq(-180, 180, 20)) +
  scale_y_continuous(limits = c(0, 1)) + 
  xlab("Longitude") + ylab(bquote(H[d]~"(mtDNA)")) + 
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
mtdna_hd_lon_plot

#############################################################################################################

######## Calculate environmental variables ########

#scale SST variables
mtdna_small_hd$sstmean_scale <- as.numeric(scale(mtdna_small_hd$sst.BO_sstmean))

#### log transform chlorophyll A ####
#subset to only those with chloroA data
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_hd$logchlomean <- log10(mtdna_small_hd$chloroA.BO_chlomean)

#remove logchlo = NA columns
mtdna_small_hd <- subset(mtdna_small_hd, mtdna_small_hd$logchlomean != "Inf" | 
                           mtdna_small_hd$logchlomean != "NaN")

##################################################################################################################

######### SST mean figure #######

SSTmean_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                              (1|Family/Genus) + (1|Source) + (1|MarkerName),
                            data = mtdna_small_hd, family = ordbeta, 
                            na.action = "na.fail")  
  
#### Predict ####
#marginal effects
SSTmean_eff <- plot_model(SSTmean_model_hd, type = "pred",
                      terms = "sstmean_scale [all]",
                      allow.new.levels = TRUE)

#pull out marginal effects dataframe
SSTmean_eff_data <- as.data.frame(SSTmean_eff$data)

#unscale SSTmean
#use same scaled:center & scaled:scale from original data
sstmean_scale <- scale(mtdna_small_hd$sst.BO_sstmean) #bc had to convert to numeric to run model/calculate marginal effects
SSTmean_eff_data$SSTmean <- (SSTmean_eff_data$x * attr(sstmean_scale, "scaled:scale")) +
  attr(sstmean_scale, "scaled:center")

#### Plot SSTmean ####

mtdna_hd_sstmean_plot <- ggplot() +
  geom_line(data = SSTmean_eff_data,
            aes(x = SSTmean, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = SSTmean_eff_data,
              aes(x = SSTmean, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  geom_rug(data = mtdna_small_hd, mapping = aes(x = sst.BO_sstmean), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 1, y = 0.985, label = "B", size = 100) + 
  ylim(c(0.5, 1.0)) + xlim(0, 30) +
  xlab("Mean SST (°C)") + ylab(bquote(H[d]~"(mtDNA)")) + 
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
mtdna_hd_sstmean_plot

###########################################################################################################

######## Chloro mean figure ########

chloromean_model_hd <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean +
                                 I(logchlomean^2) + (1|Family/Genus) + 
                                 (1|Source) + (1|MarkerName),
                               data = mtdna_small_hd, family = ordbeta, 
                               na.action = "na.fail")  

#### Predict ####
#marginal effects
chloromean_eff <- plot_model(chloromean_model_hd, type = "pred", 
                             terms = "logchlomean [all]",
                             allow.new.levels = TRUE)

#pull out marginal effects dataframe
chloromean_eff_data <- as.data.frame(chloromean_eff$data)
  chloromean_eff_data$chlomean <- 10^(chloromean_eff_data$x)

#### Plot chloromean ####

mtdna_hd_chloromean_plot <- ggplot() +
  geom_line(data = chloromean_eff_data,
            aes(x = chlomean, y = predicted), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = chloromean_eff_data,
              aes(x = chlomean, ymin = conf.low, ymax = conf.high), 
              color ="black", alpha = 0.1) +
  geom_rug(data = mtdna_small_hd, mapping = aes(x = chloroA.BO_chlomean), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 0.12, y = 0.985, label = "E", size = 100) + 
    ylim(0.5, 1.0) +
    scale_x_continuous(trans = "log10", limits = c(0.1, 10)) +
    xlab(bquote("Mean Chlorophyll"~(mg/m^3))) + ylab(bquote(H[d]~"(mtDNA)")) + 
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
mtdna_hd_chloromean_plot
################################### Script for Building Regression Figures for mtDNA pi  ############################################

#Uses mtDNA pi linear regression models
#Predicts pi at certain values of predictor variable of interest
#Calculates mean pi of raw data (binned every X units)
#Plots two together

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(tidyverse) #v.2.0.0
library(lme4) #v.1.1-31
library(DHARMa) #v.0.4.6
library(sjPlot) #v.2.8.12
library(splines) #v.4.2.2
library(data.table) #1.14.8

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

#### Add range position ####
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

#subset to only those with range_position
mtdna_small_pi <- subset(mtdna_small_pi, range_position != "NA")

#scale range_position
mtdna_small_pi$range_pos_scale <- as.numeric(scale(mtdna_small_pi$range_position))

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
mtdna_small_pi$lon_scale <- as.numeric(scale(mtdna_small_pi$lon))  

#############################################################################################################

######## Range position figure ########

null_model_pi <- lmer(logpi ~ range_pos_scale + (1|Family/Genus) + 
                        (1|Source) + (1|MarkerName),
                      REML = FALSE, data = mtdna_small_pi, 
                      na.action = "na.fail", 
                      control = lmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
rangepos_eff <- plot_model(null_model_pi, type = "pred", 
                           terms = "range_pos_scale [all]")

#pull out marginal effects dataframe
rangepos_eff_data <- as.data.frame(rangepos_eff$data)

#unscale raange_position
#use same scaled:center & scaled:scale from original data
range_pos_scale <- scale(as.numeric(mtdna_small_pi$range_position)) #bc had to convert to numeric to run model/calculate marginal effects
rangepos_eff_data$range_position <- (rangepos_eff_data$x * attr(range_pos_scale, "scaled:scale")) + 
  attr(range_pos_scale, "scaled:center")

#unlog pi
rangepos_eff_data$unlog_pi <- 10^(rangepos_eff_data$predicted)
rangepos_eff_data$unlog_conf.low <- 10^(rangepos_eff_data$conf.low)
rangepos_eff_data$unlog_conf.high <- 10^(rangepos_eff_data$conf.high)

#### Plot range position ####

mtdna_pi_rangepos_plot <- ggplot() +
  geom_line(data = rangepos_eff_data,
            aes(x = range_position, y = unlog_pi), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = rangepos_eff_data,
              aes(x = range_position, ymin = unlog_conf.low, ymax = unlog_conf.high),
              color ="black", alpha = 0.1) +
  geom_rug(data = mtdna_small_pi, mapping = aes(x = range_position), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 0.05, y = 0.0078, label = "(a)", size =90) + 
  ylim(0, 0.008) + xlim(0, 1) +
  xlab("Range Position") + ylab("π (mtDNA)") + 
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
mtdna_pi_rangepos_plot

#####################################################################################################################

######### Abslat figure #######

abslat_model_pi <- lmer(logpi ~ range_pos_scale + abslat_scale + (1|Family/Genus) +   
                          (1|Source) + (1|MarkerName) + 
                          (0 + abslat_scale|Family), 
                        REML = FALSE, data = mtdna_small_pi, 
                        na.action = "na.fail", 
                        control = lmerControl(optimizer = "bobyqa"))

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

#### Calculate medians from raw data ####
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

## Calculate medians within each 10 degree band ##
mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so can use data.table functions

#calculate median in each 10 degree band
abslat_pi_median <- mtdna_small_pi[, median(Pi), by = abslat_round]
  colnames(abslat_pi_median) <- c("abslat", "pi_median")
  abslat_pi_median$abslat <- as.numeric(as.character(abslat_pi_median$abslat))

#### Plot abslat ####

#violin plot
mtdna_pi_abslat_plot_violin <- ggplot() +
  geom_violin(data= mtdna_small_pi, aes(x = abslat_round, y = Pi), 
              fill = "#768e92", color = "#768e92") + #factor, plotted as 0-8 numerically
  geom_point(data = abslat_pi_median, aes(x = (abslat + 10)/10, y = pi_median), 
             color = "#3d4e50", size = 24) + #to put on same scale as factors, divide by 10, adding 10 because 0 actually = factor of 1 (matches 0-10 round group) 
  geom_line(data = abslat_eff_data,
            aes(x = (abslat + 5)/10, y = unlog_pi), col = "black", alpha = 0.3, linewidth = 10) + #dividing by 10 to put on same factor scale, have to add 5 because violin plots are midpoint of 10 degree bands
  geom_ribbon(data = abslat_eff_data,
              aes(x = (abslat + 5)/10, ymin = unlog_conf.low, ymax = unlog_conf.high), 
              col = "black", alpha = 0.1) +
  annotate("text", x = 7.75, y = 0.0145, label = "(a)", size = 90) +
  scale_y_continuous(limits = c(0, 0.015)) + 
  scale_x_discrete(labels = c(5, 15, 25, 35, 45, 55, 65, 75)) +
  xlab("Absolute Latitude") + ylab("π (mtDNA)") + 
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
mtdna_pi_abslat_plot_violin

#############################################################################################################

######### Lat figure #######

lat_model_pi <- lmer(logpi ~ range_pos_scale + lat_scale + I(lat_scale^2) + 
                       (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                       (0 + lat_scale|Family), 
                     REML = FALSE, data = mtdna_small_pi, 
                     na.action = "na.fail", 
                     control = lmerControl(optimizer = "bobyqa"))

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

#### Calculate medians from raw data ####
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

## Calculate medians and MADs within each 10 degree band ##
mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so can use data.table functions

#calculate median in each 10 degree band
lat_pi_median <- mtdna_small_pi[, median(Pi), by = lat_round]
  colnames(lat_pi_median) <- c("lat", "pi_median")

#calculate MAD in each 10 degree band
lat_pi_MAD <- mtdna_small_pi[, mad(Pi), by = lat_round] #median absolute deviation (less affected by outliers, measure of dispersion around median = spread of observations in a dataset)
  colnames(lat_pi_MAD) <- c("lat", "pi_MAD")

#count observations in each 10 degree band
lat_pi_count <- mtdna_small_pi[, .N, by = lat_round]
  colnames(lat_pi_count) <- c("lat", "pi_count")

#merge mean and SE dataframes together
lat_pi_binned_medians <- list(lat_pi_median, lat_pi_MAD, lat_pi_count) %>% 
  reduce(full_join, by = "lat")
lat_pi_binned_medians$lat <- as.numeric(as.character(lat_pi_binned_medians$lat))
lat_pi_binned_medians <- lat_pi_binned_medians[order(lat), ]
lat_pi_binned_medians$X <- lat_pi_binned_medians$lat + 5 #for plotting, plot in MIDDLE of 10 degree band

#calculate error bars (MAD, variance around median)
lat_pi_binned_medians$median_lowerMAD <- lat_pi_binned_medians$pi_median - 
  lat_pi_binned_medians$pi_MAD
  lat_pi_binned_medians$median_lowerMAD[lat_pi_binned_medians$median_lowerMAD < 0] <- 0.000001 #bound at ~0 for plotting
lat_pi_binned_medians$median_upperMAD <- lat_pi_binned_medians$pi_median + 
  lat_pi_binned_medians$pi_MAD

#### Plot lat ####

mtdna_pi_lat_plot <- ggplot() +
  geom_point(data = lat_pi_binned_medians, 
             aes(x = X, y = pi_median), color = "darkblue", shape = "circle", size = 24) + 
  geom_errorbar(data = lat_pi_binned_medians, 
                aes(x = X, ymin = median_lowerMAD, ymax = median_upperMAD), 
                color = "darkblue", width = 0, linewidth = 3) + 
  geom_line(data = lat_eff_data,
            aes(x = lat, y = unlog_pi), color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = lat_eff_data,
              aes(x = lat, ymin = unlog_conf.low, ymax = unlog_conf.high), 
              color ="black", alpha = 0.1) +
  annotate("text", x = -65, y = 0.0145, label = "(a)", size = 100) +
  scale_x_continuous(breaks = seq(-80, 80, 20)) +
  scale_y_continuous(limits = c(0, 0.015)) + 
  xlab("Latitude") + ylab("π (mtDNA)") +
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
mtdna_pi_lat_plot

#####################################################################################################

######## Lon figure ########

lon_model_pi_spline <- lmer(logpi ~ range_pos_scale + bs(lon_scale) + 
                               (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                              (0 + lon_scale|Family),
                            REML = FALSE, data = mtdna_small_pi, 
                            na.action = "na.fail", 
                            control = lmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
lon_eff <- plot_model(lon_model_pi_spline, type = "eff",
                         terms = "lon_scale [all]")

#pull out marginal effects dataframe
lon_eff_data <- as.data.frame(lon_eff$data)

#unscale lon
#use same scaled:center & scaled:scale from original data
lon_scale <- scale(mtdna_small_pi$lon) #bc had to convert to numeric to run model/calculate marginal effects
lon_eff_data$lon <- (lon_eff_data$x * attr(lon_scale, "scaled:scale")) + 
  attr(lon_scale, "scaled:center")

#unlog pi
lon_eff_data$unlog_pi <- 10^(lon_eff_data$predicted)
  lon_eff_data$unlog_conf.low <- 10^(lon_eff_data$conf.low)
  lon_eff_data$unlog_conf.high <- 10^(lon_eff_data$conf.high)
  
#### Calculate medians from raw data ####
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

## Calculate medians and MADs within each 10 degree band ##
mtdna_small_pi <- data.table(mtdna_small_pi) #make data.table so can use data.table functions

#calculate median in each 10 degree band
lon_pi_median <- mtdna_small_pi[, median(Pi), by = lon_round]
colnames(lon_pi_median) <- c("lon", "pi_median")

#calculate MAD in each 10 degree band
lon_pi_MAD <- mtdna_small_pi[, mad(Pi), by = lon_round] #median absolute deviation (less affected by outliers, measure of dispersion around median = spread of observations in a dataset)
colnames(lon_pi_MAD) <- c("lon", "pi_MAD")

#count observations in each 10 degree band
lon_pi_count <- mtdna_small_pi[, .N, by = lon_round]
colnames(lon_pi_count) <- c("lon", "pi_count")

#merge mean and SE dataframes together
lon_pi_binned_medians <- list(lon_pi_median, lon_pi_MAD, lon_pi_count) %>% 
  reduce(full_join, by = "lon")
lon_pi_binned_medians$lon <- as.numeric(as.character(lon_pi_binned_medians$lon))
lon_pi_binned_medians <- lon_pi_binned_medians[order(lon), ]
lon_pi_binned_medians$X <- lon_pi_binned_medians$lon + 5 #for plotting, plot in MIDDLE of 10 degree band
lon_pi_binned_medians$lon[lon_pi_binned_medians$lon == 185] <- 180 #bc nothing above 180

#calculate error bars (MAD, variance around median)
lon_pi_binned_medians$median_lowerMAD <- lon_pi_binned_medians$pi_median - 
  lon_pi_binned_medians$pi_MAD
  lon_pi_binned_medians$median_lowerMAD[lon_pi_binned_medians$median_lowerMAD < 0] <- 0.000001
lon_pi_binned_medians$median_upperMAD <- lon_pi_binned_medians$pi_median + 
  lon_pi_binned_medians$pi_MAD
  lon_pi_binned_medians$median_upperMAD[lon_pi_binned_medians$median_upperMAD > 0.015] <- 0.01499999 #for plotting

#### Plot lon ####

mtdna_pi_lon_plot <- ggplot() +
  annotate("rect", xmin = 95, xmax = 165, ymin = 0, ymax = 0.015, #adding box highlighting coral triangle
           fill = "darkolivegreen", alpha = 0.4) + 
  geom_point(data = lon_pi_binned_medians, 
             aes(x = X, y = pi_median), color = "darkblue", shape = "circle", size = 24) + 
  geom_errorbar(data = lon_pi_binned_medians, 
                aes(x = X, ymin = median_lowerMAD, ymax = median_upperMAD), 
                color = "darkblue", width = 0, linewidth = 3) + 
  geom_line(data = lon_eff_data,
            aes(x = lon, y = unlog_pi), color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = lon_eff_data,
              aes(x = lon, ymin = unlog_conf.low, ymax = unlog_conf.high), 
              color ="black", alpha = 0.1)  + 
  annotate("text", x = -175, y = 0.0145, label = "(a)", size = 100) +
  scale_x_continuous(breaks = seq(-180, 180, 20)) +
  scale_y_continuous(limits = c(0, 0.015)) + 
  xlab("Longitude") + ylab("π (mtDNA)") +
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
mtdna_pi_lon_plot

#############################################################################################################

######## Calculate environmental variables ########

#scale SST variables
mtdna_small_pi$sstmean_scale <- as.numeric(scale(mtdna_small_pi$sst.BO_sstmean))

#### log transform chlorophyll A ####
#subset to only those with chloroA data
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  mtdna_small_pi$logchlomean <- log10(mtdna_small_pi$chloroA.BO_chlomean)

#remove logchlo = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logchlomean != "Inf" | 
                           mtdna_small_pi$logchlomean != "NaN")

#######################################################################################################

######## SST mean figure ########

SSTmean_model_pi <- lmer(logpi ~ range_pos_scale + sstmean_scale + 
                           (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                           (0 + sstmean_scale|Family),
                         REML = FALSE, data = mtdna_small_pi, 
                         na.action = "na.fail", 
                         control = lmerControl(optimizer = "bobyqa"))

#### Predict ####
#marginal effects
SSTmean_eff <- plot_model(SSTmean_model_pi, type = "pred", 
                      terms = "sstmean_scale [all]")

#pull out marginal effects dataframe
SSTmean_eff_data <- as.data.frame(SSTmean_eff$data)

#unscale SSTmean
#use same scaled:center & scaled:scale from original data
sstmean_scale <- scale(mtdna_small_pi$sst.BO_sstmean) #bc had to convert to numeric to run model/calculate marginal effects
SSTmean_eff_data$SSTmean <- (SSTmean_eff_data$x * attr(sstmean_scale, "scaled:scale")) +
  attr(sstmean_scale, "scaled:center")

#unlog pi
SSTmean_eff_data$unlog_pi <- 10^(SSTmean_eff_data$predicted)
  SSTmean_eff_data$unlog_conf.low <- 10^(SSTmean_eff_data$conf.low)
  SSTmean_eff_data$unlog_conf.high <- 10^(SSTmean_eff_data$conf.high)

#### Plot SSTmean ####

mtdna_pi_sstmean_plot <- ggplot() +
  geom_line(data = SSTmean_eff_data,
            aes(x = SSTmean, y = unlog_pi), color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = SSTmean_eff_data,
              aes(x = SSTmean, ymin = unlog_conf.low, ymax = unlog_conf.high), 
              color ="black", alpha = 0.1) +
  geom_rug(data = mtdna_small_pi, mapping = aes(x = sst.BO_sstmean), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 2, y = 0.0078, label = "(a)", size = 100) + 
  ylim(0, 0.008) + xlim(0, 30) +
  xlab("Mean SST (°C)") + ylab("π (mtDNA)") + 
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
mtdna_pi_sstmean_plot

#############################################################################################################

######## ChloroA mean figure ########

chloroAmean_model_pi <- lmer(logpi ~ range_pos_scale + logchlomean + I(logchlomean^2) + 
                               (1|Family/Genus) + (1|Source) + (1|MarkerName) + 
                               (0 + logchlomean|Family),
                            REML = FALSE, data = mtdna_small_pi, 
                            na.action = "na.fail", 
                            control = lmerControl(optimizer = "bobyqa"))

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

#### plot chlomean ####

mtdna_pi_chloromean_plot <- ggplot() +
  geom_line(data = chloroAmean_eff_data,
            aes(x = chlomean, y = unlog_pi), 
            color ="black", alpha = 0.3, linewidth = 10) +
  geom_ribbon(data = chloroAmean_eff_data,
              aes(x = chlomean, ymin = unlog_conf.low, ymax = unlog_conf.high), 
              color ="black", alpha = 0.1)+
  geom_rug(data = mtdna_small_pi, mapping = aes(x = chloroA.BO_chlomean), 
           color = "#282828", inherit.aes = FALSE) + 
  annotate("text", x = 0.135, y = 0.0078, label = "(d)", size = 100) + 
  ylim(0, 0.008) +
  scale_x_continuous(trans = "log10", limits = c(0.1, 10)) +
  xlab(bquote("Mean Chlorophyll"~(mg/m^3))) + ylab("π (mtDNA)") + 
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
mtdna_pi_chloromean_plot

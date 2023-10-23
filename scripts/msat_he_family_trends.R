#################################### Script to look at family trends for msat he models ########################################################

#Uses msat he beta regression models
#Predicts he at certain values of predictor variable of interest
#Broken down by family

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

######## Count observations by family ########

msat_he_family <- c(unique(msat$Family)) #78 families

#how many observations for each family? really want to look at trends in the ones with enough data
msat <- as.data.table(msat)

he_family <- msat[, .N, by = .(Family)]
  he_family <- he_family[order(-N), ] #sort in descending order (families with most observations at top)

#10 families with most observations
#Clupeidae, Sebastidae, Sciaenidae, Gadidae, Serranidae, Pomacentridae, Sparidae, Scombridae, Lutjanidae, Pleuronectidae

#WORKING WITH THESE FAMILIES:
#bc fall into top 10 across all datasets and/or represent interesting life histories AND have a decent amount of data (>~30 obs)
#Scombridae, Lutjanidae, Serranidae, Pomacentridae, Sebastidae, Engraulidae, Gadidae, Syngnathidae, Rajidae, Carcharhinidae

###########################################################################################################

######### Abslat figure #######

#abslat models
abslat_model_he_Scombridae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                        abslat_scale + (1|Source),
                                      data = msat[msat$Family == "Scombridae", ], 
                                      na.action = "na.fail", family = ordbeta)

abslat_model_he_Lutjanidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                        abslat_scale + (1|Source),
                                      data = msat[msat$Family == "Lutjanidae", ], 
                                      na.action = "na.fail", family = ordbeta)

abslat_model_he_Serranidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                        abslat_scale + (1|Source),
                                      data = msat[msat$Family == "Serranidae", ], 
                                      na.action = "na.fail", family = ordbeta)

abslat_model_he_Pomacentridae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                           abslat_scale + (1|Source),
                                         data = msat[msat$Family == "Pomacentridae", ], 
                                         na.action = "na.fail", family = ordbeta)

abslat_model_he_Sebastidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                        abslat_scale + (1|Source),
                                      data = msat[msat$Family == "Sebastidae", ], 
                                      na.action = "na.fail", family = ordbeta)

abslat_model_he_Engraulidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                         abslat_scale + (1|Source),
                                       data = msat[msat$Family == "Engraulidae", ], 
                                       na.action = "na.fail", family =  ordbeta)

abslat_model_he_Gadidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                     abslat_scale + (1|Source),
                                   data = msat[msat$Family == "Gadidae", ], 
                                   na.action = "na.fail", family = ordbeta)

abslat_model_he_Syngnathidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                          abslat_scale + (1|Source),
                                        data = msat[msat$Family == "Syngnathidae", ], 
                                        na.action = "na.fail", family = ordbeta)

abslat_model_he_Rajidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                     abslat_scale + (1|Source),
                                   data = msat[msat$Family == "Rajidae", ], 
                                   na.action = "na.fail", family = ordbeta)

abslat_model_he_Carcharhinidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                            abslat_scale + (1|Source),
                                          data = msat[msat$Family == "Carcharhinidae", ], 
                                          na.action = "na.fail", family = ordbeta)

#### Predict ####
#marginal effects
abslat_eff_Scombridae <- plot_model(abslat_model_he_Scombridae, type = "pred", 
                                    terms = "abslat_scale [all]")
abslat_eff_Lutjanidae <- plot_model(abslat_model_he_Lutjanidae, type = "pred", 
                                    terms = "abslat_scale [all]")
abslat_eff_Serranidae <- plot_model(abslat_model_he_Serranidae, type = "pred", 
                                    terms = "abslat_scale [all]")
abslat_eff_Pomacentridae <- plot_model(abslat_model_he_Pomacentridae, type = "pred", 
                                       terms = "abslat_scale [all]")
abslat_eff_Sebastidae <- plot_model(abslat_model_he_Sebastidae, type = "pred", 
                                    terms = "abslat_scale [all]")
abslat_eff_Engraulidae <- plot_model(abslat_model_he_Engraulidae, type = "pred", 
                                     terms = "abslat_scale [all]")
abslat_eff_Gadidae <- plot_model(abslat_model_he_Gadidae, type = "pred", 
                                 terms = "abslat_scale [all]")
abslat_eff_Syngnathidae <- plot_model(abslat_model_he_Syngnathidae, type = "pred", 
                                      terms = "abslat_scale [all]")
abslat_eff_Rajidae <- plot_model(abslat_model_he_Rajidae, type = "pred", 
                                 terms = "abslat_scale [all]")
abslat_eff_Carcharhinidae <- plot_model(abslat_model_he_Carcharhinidae, type = "pred", 
                                        terms = "abslat_scale [all]")

#pull out marginal effects dataframes
abslat_eff_data_Scombridae <- as.data.frame(abslat_eff_Scombridae$data)
  abslat_eff_data_Scombridae$Family <- "Scombridae"
abslat_eff_data_Lutjanidae <- as.data.frame(abslat_eff_Lutjanidae$data)
  abslat_eff_data_Lutjanidae$Family <- "Lutjanidae"
abslat_eff_data_Serranidae <- as.data.frame(abslat_eff_Serranidae$data)
  abslat_eff_data_Serranidae$Family <- "Serranidae"
abslat_eff_data_Pomacentridae <- as.data.frame(abslat_eff_Pomacentridae$data)
  abslat_eff_data_Pomacentridae$Family <- "Pomacentridae"
abslat_eff_data_Sebastidae <- as.data.frame(abslat_eff_Sebastidae$data)
  abslat_eff_data_Sebastidae$Family <- "Sebastidae"
abslat_eff_data_Engraulidae <- as.data.frame(abslat_eff_Engraulidae$data)
  abslat_eff_data_Engraulidae$Family <- "Engraulidae"
abslat_eff_data_Gadidae <- as.data.frame(abslat_eff_Gadidae$data)
  abslat_eff_data_Gadidae$Family <- "Gadidae"
abslat_eff_data_Syngnathidae <- as.data.frame(abslat_eff_Syngnathidae$data)
  abslat_eff_data_Syngnathidae$Family <- "Syngnathidae"
abslat_eff_data_Rajidae <- as.data.frame(abslat_eff_Rajidae$data)
  abslat_eff_data_Rajidae$Family <- "Rajidae"
abslat_eff_data_Carcharhinidae <- as.data.frame(abslat_eff_Carcharhinidae$data)
  abslat_eff_data_Carcharhinidae$Family <- "Carcharhinidae"

#merge together
abslat_eff_data_allfam <- rbind(abslat_eff_data_Scombridae, abslat_eff_data_Lutjanidae,
                                abslat_eff_data_Serranidae, abslat_eff_data_Pomacentridae,
                                abslat_eff_data_Sebastidae, abslat_eff_data_Engraulidae,
                                abslat_eff_data_Gadidae, abslat_eff_data_Syngnathidae,
                                abslat_eff_data_Rajidae, abslat_eff_data_Carcharhinidae)

#use same scaled:center & scaled:scale from original data
abslat_scale <- scale(msat$abslat) #bc had to convert to numeric to run model/calculate marginal effects
  abslat_eff_data_allfam$abslat <- (abslat_eff_data_allfam$x * attr(abslat_scale, "scaled:scale")) +
    attr(abslat_scale, "scaled:center")

#### Plot abslat ####

msat_he_abslat_plot <- ggplot() +
  geom_line(data = abslat_eff_data_allfam,
            aes(x = abslat, y = predicted, color = Family), linewidth = 14) +
  annotate("text", x = 6, y = 0.97, label = "C", size = 100) +
  scale_color_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                                "#AA4499", "#44AA99", "#999933", "#882255", "#888888")) + 
  ylim(0.5, 1) + xlim(0, 90) +
  xlab("Absolute Latitude") + ylab(bquote(H[e]~"(nucDNA)")) + 
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
msat_he_abslat_plot

###########################################################################################################

######### Lat figure #######

#lat models
lat_model_he_Scombridae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                     lat_scale + I(lat_scale^2) + (1|Source),
                                   data = msat[msat$Family == "Scombridae", ], 
                                   na.action = "na.fail", family = ordbeta)

lat_model_he_Lutjanidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                     lat_scale + I(lat_scale^2) + (1|Source),
                                   data = msat[msat$Family == "Lutjanidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lat_model_he_Serranidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                     lat_scale + I(lat_scale^2) + (1|Source),
                                   data = msat[msat$Family == "Serranidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lat_model_he_Pomacentridae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                        lat_scale + I(lat_scale^2) + (1|Source),
                                      data = msat[msat$Family == "Pomacentridae", ], 
                                      na.action = "na.fail", family = ordbeta)

lat_model_he_Sebastidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                     lat_scale + I(lat_scale^2) + (1|Source),
                                   data = msat[msat$Family == "Sebastidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lat_model_he_Engraulidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                      lat_scale + I(lat_scale^2) + (1|Source),
                                    data = msat[msat$Family == "Engraulidae", ], 
                                    na.action = "na.fail", family =  ordbeta)

lat_model_he_Gadidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                  lat_scale + I(lat_scale^2) + (1|Source),
                                data = msat[msat$Family == "Gadidae", ], 
                                na.action = "na.fail", family = ordbeta)

lat_model_he_Syngnathidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                       lat_scale + I(lat_scale^2) + (1|Source),
                                     data = msat[msat$Family == "Syngnathidae", ], 
                                     na.action = "na.fail", family = ordbeta)

lat_model_he_Rajidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                  lat_scale + I(lat_scale^2) + (1|Source),
                                data = msat[msat$Family == "Rajidae", ], 
                                na.action = "na.fail", family = ordbeta)

lat_model_he_Carcharhinidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                         lat_scale + I(lat_scale^2) + (1|Source),
                                       data = msat[msat$Family == "Carcharhinidae", ], 
                                       na.action = "na.fail", family = ordbeta)

#### Predict ####
#marginal effects
lat_eff_Scombridae <- plot_model(lat_model_he_Scombridae, type = "pred", 
                                 terms = "lat_scale [all]")
lat_eff_Lutjanidae <- plot_model(lat_model_he_Lutjanidae, type = "pred", 
                                 terms = "lat_scale [all]")
lat_eff_Serranidae <- plot_model(lat_model_he_Serranidae, type = "pred", 
                                 terms = "lat_scale [all]")
lat_eff_Pomacentridae <- plot_model(lat_model_he_Pomacentridae, type = "pred", 
                                    terms = "lat_scale [all]")
lat_eff_Sebastidae <- plot_model(lat_model_he_Sebastidae, type = "pred", 
                                 terms = "lat_scale [all]")
lat_eff_Engraulidae <- plot_model(lat_model_he_Engraulidae, type = "pred", 
                                  terms = "lat_scale [all]")
lat_eff_Gadidae <- plot_model(lat_model_he_Gadidae, type = "pred", 
                              terms = "lat_scale [all]")
lat_eff_Syngnathidae <- plot_model(lat_model_he_Syngnathidae, type = "pred", 
                                   terms = "lat_scale [all]")
lat_eff_Rajidae <- plot_model(lat_model_he_Rajidae, type = "pred", 
                              terms = "lat_scale [all]")
lat_eff_Carcharhinidae <- plot_model(lat_model_he_Carcharhinidae, type = "pred", 
                                     terms = "lat_scale [all]")

#pull out marginal effects dataframes
lat_eff_data_Scombridae <- as.data.frame(lat_eff_Scombridae$data)
  lat_eff_data_Scombridae$Family <- "Scombridae"
lat_eff_data_Lutjanidae <- as.data.frame(lat_eff_Lutjanidae$data)
  lat_eff_data_Lutjanidae$Family <- "Lutjanidae"
lat_eff_data_Serranidae <- as.data.frame(lat_eff_Serranidae$data)
  lat_eff_data_Serranidae$Family <- "Serranidae"
lat_eff_data_Pomacentridae <- as.data.frame(lat_eff_Pomacentridae$data)
  lat_eff_data_Pomacentridae$Family <- "Pomacentridae"
lat_eff_data_Sebastidae <- as.data.frame(lat_eff_Sebastidae$data)
  lat_eff_data_Sebastidae$Family <- "Sebastidae"
lat_eff_data_Engraulidae <- as.data.frame(lat_eff_Engraulidae$data)
  lat_eff_data_Engraulidae$Family <- "Engraulidae"
lat_eff_data_Gadidae <- as.data.frame(lat_eff_Gadidae$data)
  lat_eff_data_Gadidae$Family <- "Gadidae"
lat_eff_data_Syngnathidae <- as.data.frame(lat_eff_Syngnathidae$data)
  lat_eff_data_Syngnathidae$Family <- "Syngnathidae"
lat_eff_data_Rajidae <- as.data.frame(lat_eff_Rajidae$data)
  lat_eff_data_Rajidae$Family <- "Rajidae"
lat_eff_data_Carcharhinidae <- as.data.frame(lat_eff_Carcharhinidae$data)
  lat_eff_data_Carcharhinidae$Family <- "Carcharhinidae"

#merge together
lat_eff_data_allfam <- rbind(lat_eff_data_Scombridae, lat_eff_data_Lutjanidae,
                             lat_eff_data_Serranidae, lat_eff_data_Pomacentridae,
                             lat_eff_data_Sebastidae, lat_eff_data_Engraulidae,
                             lat_eff_data_Gadidae, lat_eff_data_Syngnathidae,
                             lat_eff_data_Rajidae, lat_eff_data_Carcharhinidae)

#use same scaled:center & scaled:scale from original data
lat_scale <- scale(msat$lat) #bc had to convert to numeric to run model/calculate marginal effects
  lat_eff_data_allfam$lat <- (lat_eff_data_allfam$x * attr(lat_scale, "scaled:scale")) +
    attr(lat_scale, "scaled:center")

#### Plot lat ####

msat_he_lat_plot <- ggplot() +
  geom_line(data = lat_eff_data_allfam,
            aes(x = lat, y = predicted, color = Family), linewidth = 14) +
annotate("text", x = -70, y = 0.97, label = "C", size = 100) +
  scale_color_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                                "#AA4499", "#44AA99", "#999933", "#882255", "#888888")) + 
  ylim(0.5, 1) + xlim(-90, 90) +
  xlab("Latitude") + ylab(bquote(H[e]~"(nucDNA)")) + 
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
msat_he_lat_plot

###########################################################################################################

######### Lon figure #######

#lon models
lon_model_he_Scombridae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                     bs(lon_scale) + (1|Source),
                                   data = msat[msat$Family == "Scombridae", ], 
                                   na.action = "na.fail", family = ordbeta)

lon_model_he_Lutjanidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                     bs(lon_scale) + (1|Source),
                                   data = msat[msat$Family == "Lutjanidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lon_model_he_Serranidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                     bs(lon_scale) + (1|Source),
                                   data = msat[msat$Family == "Serranidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lon_model_he_Pomacentridae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                        bs(lon_scale) + (1|Source),
                                      data = msat[msat$Family == "Pomacentridae", ], 
                                      na.action = "na.fail", family = ordbeta)

lon_model_he_Sebastidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                     bs(lon_scale) + (1|Source),
                                   data = msat[msat$Family == "Sebastidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lon_model_he_Engraulidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                      bs(lon_scale) + (1|Source),
                                    data = msat[msat$Family == "Engraulidae", ], 
                                    na.action = "na.fail", family =  ordbeta)

lon_model_he_Gadidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                  bs(lon_scale) + (1|Source),
                                data = msat[msat$Family == "Gadidae", ], 
                                na.action = "na.fail", family = ordbeta)

lon_model_he_Syngnathidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                       bs(lon_scale) + (1|Source),
                                     data = msat[msat$Family == "Syngnathidae", ], 
                                     na.action = "na.fail", family = ordbeta)

lon_model_he_Rajidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                  bs(lon_scale) + (1|Source),
                                data = msat[msat$Family == "Rajidae", ], 
                                na.action = "na.fail", family = ordbeta)

lon_model_he_Carcharhinidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                         bs(lon_scale) + (1|Source),
                                       data = msat[msat$Family == "Carcharhinidae", ], 
                                       na.action = "na.fail", family = ordbeta)

#### Predict ####
#marginal effects
lon_eff_Scombridae <- plot_model(lon_model_he_Scombridae, type = "pred", 
                                 terms = "lon_scale [all]")
lon_eff_Lutjanidae <- plot_model(lon_model_he_Lutjanidae, type = "pred", 
                                 terms = "lon_scale [all]")
lon_eff_Serranidae <- plot_model(lon_model_he_Serranidae, type = "pred", 
                                 terms = "lon_scale [all]")
lon_eff_Pomacentridae <- plot_model(lon_model_he_Pomacentridae, type = "pred", 
                                    terms = "lon_scale [all]")
lon_eff_Sebastidae <- plot_model(lon_model_he_Sebastidae, type = "pred", 
                                 terms = "lon_scale [all]")
lon_eff_Engraulidae <- plot_model(lon_model_he_Engraulidae, type = "pred", 
                                  terms = "lon_scale [all]")
lon_eff_Gadidae <- plot_model(lon_model_he_Gadidae, type = "pred", 
                              terms = "lon_scale [all]")
lon_eff_Syngnathidae <- plot_model(lon_model_he_Syngnathidae, type = "pred", 
                                   terms = "lon_scale [all]")
lon_eff_Rajidae <- plot_model(lon_model_he_Rajidae, type = "pred", 
                              terms = "lon_scale [all]")
lon_eff_Carcharhinidae <- plot_model(lon_model_he_Carcharhinidae, type = "pred", 
                                     terms = "lon_scale [all]")

#pull out marginal effects dataframes
lon_eff_data_Scombridae <- as.data.frame(lon_eff_Scombridae$data)
  lon_eff_data_Scombridae$Family <- "Scombridae"
lon_eff_data_Lutjanidae <- as.data.frame(lon_eff_Lutjanidae$data)
  lon_eff_data_Lutjanidae$Family <- "Lutjanidae"
lon_eff_data_Serranidae <- as.data.frame(lon_eff_Serranidae$data)
  lon_eff_data_Serranidae$Family <- "Serranidae"
lon_eff_data_Pomacentridae <- as.data.frame(lon_eff_Pomacentridae$data)
  lon_eff_data_Pomacentridae$Family <- "Pomacentridae"
lon_eff_data_Sebastidae <- as.data.frame(lon_eff_Sebastidae$data)
  lon_eff_data_Sebastidae$Family <- "Sebastidae"
lon_eff_data_Engraulidae <- as.data.frame(lon_eff_Engraulidae$data)
  lon_eff_data_Engraulidae$Family <- "Engraulidae"
lon_eff_data_Gadidae <- as.data.frame(lon_eff_Gadidae$data)
  lon_eff_data_Gadidae$Family <- "Gadidae"
lon_eff_data_Syngnathidae <- as.data.frame(lon_eff_Syngnathidae$data)
  lon_eff_data_Syngnathidae$Family <- "Syngnathidae"
lon_eff_data_Rajidae <- as.data.frame(lon_eff_Rajidae$data)
  lon_eff_data_Rajidae$Family <- "Rajidae"
lon_eff_data_Carcharhinidae <- as.data.frame(lon_eff_Carcharhinidae$data)
  lon_eff_data_Carcharhinidae$Family <- "Carcharhinidae"
  
#merge together
lon_eff_data_allfam <- rbind(lon_eff_data_Scombridae, lon_eff_data_Lutjanidae,
                             lon_eff_data_Serranidae, lon_eff_data_Pomacentridae,
                             lon_eff_data_Sebastidae, lon_eff_data_Engraulidae,
                             lon_eff_data_Gadidae, lon_eff_data_Syngnathidae,
                             lon_eff_data_Rajidae, lon_eff_data_Carcharhinidae)

#use same scaled:center & scaled:scale from original data
lon_scale <- scale(msat$lon) #bc had to convert to numeric to run model/calculate marginal effects
  lon_eff_data_allfam$lon <- (lon_eff_data_allfam$x * attr(lon_scale, "scaled:scale")) +
    attr(lon_scale, "scaled:center")

#### Plot lon ####

msat_he_lon_plot <- ggplot() +
  annotate("rect", xmin = 95, xmax = 165, ymin = 0.5, ymax = 1, #adding box highlighting coral triangle
           fill = "darkolivegreen", alpha = 0.4) + 
  geom_line(data = lon_eff_data_allfam,
            aes(x = lon, y = predicted, color = Family), linewidth = 14) +
  annotate("text", x = -175, y = 0.97, label = "C", size = 100) +
  scale_color_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                                "#AA4499", "#44AA99", "#999933", "#882255", "#888888")) + 
  ylim(0.5, 1) + xlim(-180, 180) +
  xlab("Longitude") + ylab(bquote(H[e]~"(nucDNA)")) + 
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
msat_he_lon_plot

#############################################################################################################

######## Calculate environmental variables ########

#scale SST variables
msat$sstmean_scale <- as.numeric(scale(msat$sst.BO_sstmean))

#### log transform chlorophyll A ####
#subset to only those with chloroA data
msat <- subset(msat, msat$chloroA.BO_chlomean != 0) #if any zeros will screw up log transformation (log10(0) is undefined)
  msat$logchlomean <- log10(msat$chloroA.BO_chlomean)

#remove logchlo = NA columns
msat <- subset(msat, msat$logchlomean != "Inf" | 
                 msat$logchlomean != "NaN")

#######################################################################################################

######### SSTmean figure #######

#sstmean models
sstmean_model_he_Scombridae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                         sstmean_scale + (1|Source),
                                       data = msat[msat$Family == "Scombridae", ], 
                                       na.action = "na.fail", family = ordbeta)

sstmean_model_he_Lutjanidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                         sstmean_scale + (1|Source),
                                       data = msat[msat$Family == "Lutjanidae", ], 
                                       na.action = "na.fail", family = ordbeta)

sstmean_model_he_Serranidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                         sstmean_scale + (1|Source),
                                       data = msat[msat$Family == "Serranidae", ], 
                                       na.action = "na.fail", family = ordbeta)

sstmean_model_he_Pomacentridae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                            sstmean_scale + (1|Source),
                                          data = msat[msat$Family == "Pomacentridae", ], 
                                          na.action = "na.fail", family = ordbeta)

sstmean_model_he_Sebastidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                         sstmean_scale + (1|Source),
                                       data = msat[msat$Family == "Sebastidae", ], 
                                       na.action = "na.fail", family = ordbeta)

sstmean_model_he_Engraulidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                          sstmean_scale + (1|Source),
                                        data = msat[msat$Family == "Engraulidae", ], 
                                        na.action = "na.fail", family =  ordbeta)

sstmean_model_he_Gadidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                      sstmean_scale + (1|Source),
                                    data = msat[msat$Family == "Gadidae", ], 
                                    na.action = "na.fail", family = ordbeta)

sstmean_model_he_Syngnathidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                           sstmean_scale + (1|Source),
                                         data = msat[msat$Family == "Syngnathidae", ], 
                                         na.action = "na.fail", family = ordbeta)
    
sstmean_model_he_Rajidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                      sstmean_scale + (1|Source),
                                    data = msat[msat$Family == "Rajidae", ], 
                                    na.action = "na.fail", family = ordbeta)

sstmean_model_he_Carcharhinidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                             sstmean_scale + (1|Source),
                                           data = msat[msat$Family == "Carcharhinidae", ], 
                                           na.action = "na.fail", family = ordbeta)

#### Predict ####
#marginal effects
sstmean_eff_Scombridae <- plot_model(sstmean_model_he_Scombridae, type = "pred", 
                                     terms = "sstmean_scale [all]")
sstmean_eff_Lutjanidae <- plot_model(sstmean_model_he_Lutjanidae, type = "pred", 
                                     terms = "sstmean_scale [all]")
sstmean_eff_Serranidae <- plot_model(sstmean_model_he_Serranidae, type = "pred", 
                                     terms = "sstmean_scale [all]")
sstmean_eff_Pomacentridae <- plot_model(sstmean_model_he_Pomacentridae, type = "pred", 
                                        terms = "sstmean_scale [all]")
sstmean_eff_Sebastidae <- plot_model(sstmean_model_he_Sebastidae, type = "pred", 
                                     terms = "sstmean_scale [all]")
sstmean_eff_Engraulidae <- plot_model(sstmean_model_he_Engraulidae, type = "pred", 
                                      terms = "sstmean_scale [all]")
sstmean_eff_Gadidae <- plot_model(sstmean_model_he_Gadidae, type = "pred", 
                                  terms = "sstmean_scale [all]")
sstmean_eff_Syngnathidae <- plot_model(sstmean_model_he_Syngnathidae, type = "pred", 
                                       terms = "sstmean_scale [all]")
sstmean_eff_Rajidae <- plot_model(sstmean_model_he_Rajidae, type = "pred", 
                                  terms = "sstmean_scale [all]")
sstmean_eff_Carcharhinidae <- plot_model(sstmean_model_he_Carcharhinidae, type = "pred", 
                                         terms = "sstmean_scale [all]")

#pull out marginal effects dataframes
sstmean_eff_data_Scombridae <- as.data.frame(sstmean_eff_Scombridae$data)
  sstmean_eff_data_Scombridae$Family <- "Scombridae"
sstmean_eff_data_Lutjanidae <- as.data.frame(sstmean_eff_Lutjanidae$data)
  sstmean_eff_data_Lutjanidae$Family <- "Lutjanidae"
sstmean_eff_data_Serranidae <- as.data.frame(sstmean_eff_Serranidae$data)
  sstmean_eff_data_Serranidae$Family <- "Serranidae"
sstmean_eff_data_Pomacentridae <- as.data.frame(sstmean_eff_Pomacentridae$data)
  sstmean_eff_data_Pomacentridae$Family <- "Pomacentridae"
sstmean_eff_data_Sebastidae <- as.data.frame(sstmean_eff_Sebastidae$data)
  sstmean_eff_data_Sebastidae$Family <- "Sebastidae"
sstmean_eff_data_Engraulidae <- as.data.frame(sstmean_eff_Engraulidae$data)
  sstmean_eff_data_Engraulidae$Family <- "Engraulidae"
sstmean_eff_data_Gadidae <- as.data.frame(sstmean_eff_Gadidae$data)
  sstmean_eff_data_Gadidae$Family <- "Gadidae"
sstmean_eff_data_Syngnathidae <- as.data.frame(sstmean_eff_Syngnathidae$data)
  sstmean_eff_data_Syngnathidae$Family <- "Syngnathidae"
sstmean_eff_data_Rajidae <- as.data.frame(sstmean_eff_Rajidae$data)
  sstmean_eff_data_Rajidae$Family <- "Rajidae"
sstmean_eff_data_Carcharhinidae <- as.data.frame(sstmean_eff_Carcharhinidae$data)
  sstmean_eff_data_Carcharhinidae$Family <- "Carcharhinidae"

#fix Sebastidae dataframe (didn't make conf intervals??)
sstmean_eff_data_Rajidae$std.error <- "NA"
  sstmean_eff_data_Rajidae$conf.low <- "NA"
  sstmean_eff_data_Rajidae$conf.high <- "NA"

sstmean_eff_data_Rajidae <- sstmean_eff_data_Rajidae[, c("x", "predicted", 
                                                         "std.error", "conf.low",
                                                         "conf.high", "group", 
                                                         "group_col", "Family")]
  
#merge together
sstmean_eff_data_allfam <- rbind(sstmean_eff_data_Scombridae, sstmean_eff_data_Lutjanidae,
                                 sstmean_eff_data_Serranidae, sstmean_eff_data_Pomacentridae,
                                 sstmean_eff_data_Sebastidae, sstmean_eff_data_Engraulidae,
                                 sstmean_eff_data_Gadidae, sstmean_eff_data_Syngnathidae,
                                 sstmean_eff_data_Rajidae, sstmean_eff_data_Carcharhinidae)

#use same scaled:center & scaled:scale from original data
sstmean_scale <- scale(msat$sst.BO_sstmean) #bc had to convert to numeric to run model/calculate marginal effects
  sstmean_eff_data_allfam$sstmean <- (sstmean_eff_data_allfam$x * attr(sstmean_scale, "scaled:scale")) +
    attr(sstmean_scale, "scaled:center")

#### Plot SSTmean ####

msat_he_sstmean_plot <- ggplot() +
  geom_line(data = sstmean_eff_data_allfam,
            aes(x = sstmean, y = predicted, color = Family), linewidth = 14) +
  annotate("text", x = 1, y = 0.985, label = "C", size = 100) + 
  scale_color_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                                "#AA4499", "#44AA99", "#999933", "#882255", "#888888")) + 
  ylim(0.5, 1) + xlim(0, 30) +
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

###########################################################################################################

######### chloro mean figure #######

#chlomean models
chlomean_model_he_Scombridae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                          logchlomean + I(logchlomean^2) + (1|Source),
                                        data = msat[msat$Family == "Scombridae", ], 
                                        na.action = "na.fail", family = ordbeta)

chlomean_model_he_Lutjanidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                          logchlomean + I(logchlomean^2) + (1|Source),
                                        data = msat[msat$Family == "Lutjanidae", ], 
                                        na.action = "na.fail", family = ordbeta)

chlomean_model_he_Serranidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                          logchlomean + I(logchlomean^2) + (1|Source),
                                        data = msat[msat$Family == "Serranidae", ], 
                                        na.action = "na.fail", family = ordbeta)

chlomean_model_he_Pomacentridae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                             logchlomean + I(logchlomean^2) + (1|Source),
                                           data = msat[msat$Family == "Pomacentridae", ], 
                                           na.action = "na.fail", family = ordbeta)

chlomean_model_he_Sebastidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                          logchlomean + I(logchlomean^2) + (1|Source),
                                        data = msat[msat$Family == "Sebastidae", ], 
                                        na.action = "na.fail", family = ordbeta)

chlomean_model_he_Engraulidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                           logchlomean + I(logchlomean^2) + (1|Source),
                                         data = msat[msat$Family == "Engraulidae", ], 
                                         na.action = "na.fail", family =  ordbeta)

chlomean_model_he_Gadidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                       logchlomean + I(logchlomean^2) + (1|Source),
                                     data = msat[msat$Family == "Gadidae", ], 
                                     na.action = "na.fail", family = ordbeta)

chlomean_model_he_Syngnathidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                            logchlomean + I(logchlomean^2) + (1|Source),
                                          data = msat[msat$Family == "Syngnathidae", ], 
                                          na.action = "na.fail", family = ordbeta)

chlomean_model_he_Rajidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                       logchlomean + I(logchlomean^2) + (1|Source),
                                     data = msat[msat$Family == "Rajidae", ], 
                                     na.action = "na.fail", family = ordbeta)

chlomean_model_he_Carcharhinidae <- glmmTMB(He ~ CrossSpp_scale + range_pos_scale + 
                                              logchlomean + I(logchlomean^2) + (1|Source),
                                            data = msat[msat$Family == "Carcharhinidae", ], 
                                            na.action = "na.fail", family = ordbeta)

#### Predict ####
#marginal effects
chlomean_eff_Scombridae <- plot_model(chlomean_model_he_Scombridae, type = "pred", 
                                      terms = "logchlomean [all]")
chlomean_eff_Lutjanidae <- plot_model(chlomean_model_he_Lutjanidae, type = "pred", 
                                      terms = "logchlomean [all]")
chlomean_eff_Serranidae <- plot_model(chlomean_model_he_Serranidae, type = "pred", 
                                      terms = "logchlomean [all]")
chlomean_eff_Pomacentridae <- plot_model(chlomean_model_he_Pomacentridae, type = "pred", 
                                         terms = "logchlomean [all]")
chlomean_eff_Sebastidae <- plot_model(chlomean_model_he_Sebastidae, type = "pred", 
                                      terms = "logchlomean [all]")
chlomean_eff_Engraulidae <- plot_model(chlomean_model_he_Engraulidae, type = "pred", 
                                       terms = "logchlomean [all]")
chlomean_eff_Gadidae <- plot_model(chlomean_model_he_Gadidae, type = "pred", 
                                   terms = "logchlomean [all]")
chlomean_eff_Syngnathidae <- plot_model(chlomean_model_he_Syngnathidae, type = "pred", 
                                        terms = "logchlomean [all]")
chlomean_eff_Rajidae <- plot_model(chlomean_model_he_Rajidae, type = "pred", 
                                   terms = "logchlomean [all]")
chlomean_eff_Carcharhinidae <- plot_model(chlomean_model_he_Carcharhinidae, type = "pred", 
                                          terms = "logchlomean [all]")

#pull out marginal effects dataframes
chlomean_eff_data_Scombridae <- as.data.frame(chlomean_eff_Scombridae$data)
  chlomean_eff_data_Scombridae$Family <- "Scombridae"
chlomean_eff_data_Lutjanidae <- as.data.frame(chlomean_eff_Lutjanidae$data)
  chlomean_eff_data_Lutjanidae$Family <- "Lutjanidae"
chlomean_eff_data_Serranidae <- as.data.frame(chlomean_eff_Serranidae$data)
  chlomean_eff_data_Serranidae$Family <- "Serranidae"
chlomean_eff_data_Pomacentridae <- as.data.frame(chlomean_eff_Pomacentridae$data)
  chlomean_eff_data_Pomacentridae$Family <- "Pomacentridae"
chlomean_eff_data_Sebastidae <- as.data.frame(chlomean_eff_Sebastidae$data)
  chlomean_eff_data_Sebastidae$Family <- "Sebastidae"
chlomean_eff_data_Engraulidae <- as.data.frame(chlomean_eff_Engraulidae$data)
  chlomean_eff_data_Engraulidae$Family <- "Engraulidae"
chlomean_eff_data_Gadidae <- as.data.frame(chlomean_eff_Gadidae$data)
  chlomean_eff_data_Gadidae$Family <- "Gadidae"
chlomean_eff_data_Syngnathidae <- as.data.frame(chlomean_eff_Syngnathidae$data)
  chlomean_eff_data_Syngnathidae$Family <- "Syngnathidae"
chlomean_eff_data_Rajidae <- as.data.frame(chlomean_eff_Rajidae$data)
  chlomean_eff_data_Rajidae$Family <- "Rajidae"
chlomean_eff_data_Carcharhinidae <- as.data.frame(chlomean_eff_Carcharhinidae$data)
  chlomean_eff_data_Carcharhinidae$Family <- "Carcharhinidae"

#merge together
chlomean_eff_data_allfam <- rbind(chlomean_eff_data_Scombridae, chlomean_eff_data_Lutjanidae,
                                  chlomean_eff_data_Serranidae, chlomean_eff_data_Pomacentridae,
                                  chlomean_eff_data_Sebastidae, chlomean_eff_data_Engraulidae,
                                  chlomean_eff_data_Gadidae, chlomean_eff_data_Syngnathidae,
                                  chlomean_eff_data_Rajidae, chlomean_eff_data_Carcharhinidae)

#unlog chlorophyll mean
chlomean_eff_data_allfam$chlomean <- 10^(chlomean_eff_data_allfam$x)

#### Plot chlomean ####

msat_he_chlomean_plot <- ggplot() +
  geom_line(data = chlomean_eff_data_allfam,
            aes(x = chlomean, y = predicted, color = Family), linewidth = 14) +
  annotate("text", x = 0.12, y = 0.985, label = "F", size = 100) + 
  scale_color_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                                "#AA4499", "#44AA99", "#999933", "#882255", "#888888")) + 
  ylim(0.4, 1.0) +
  scale_x_continuous(trans = "log10", limits = c(0.1, 10)) +
  xlab(bquote("Mean Chlorophyll"~(mg/m^3))) + ylab(bquote(H[e]~"(nucDNA)")) + 
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
msat_he_chlomean_plot

#################################### Script to look at family trends for mtDNA hd models ########################################################

#Uses mtDNA hd beta regression models
#Predicts hd at certain values of predictor variable of interest
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

###########################################################################################################

######## Count observations by family ########

mtdna_hd_family <- c(unique(mtdna_small_hd$Family)) #74 families

#how many observations for each family? really want to look at trends in the ones with enough data
mtdna_small_hd <- as.data.table(mtdna_small_hd)

hd_family <- mtdna_small_hd[, .N, by = .(Family)]
  hd_family <- hd_family[order(-N), ] #sort in descending order (families with most observations at top)

#10 families with most observations
#Scombridae, Lutjanidae, Serranidae, Sciaenidae, Sparidae, Acanthuridae, Engraulidae, Rajidae, Gadidae, Pomacentridae

#WORKING WITH THESE FAMILIES:
#bc fall into top 10 across all datasets and/or represent interesting life histories AND have a decent amount of data (>~30 obs)
#Scombridae, Lutjanidae, Serranidae, Pomacentridae, Sebastidae, Engraulidae, Gadidae, Syngnathidae, Rajidae, Carcharhinidae

###########################################################################################################

######### Abslat figure #######

#abslat models
abslat_model_hd_Scombridae <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale + 
                                        (1|Source) + (1|MarkerName),
                                      data = mtdna_small_hd[mtdna_small_hd$Family == "Scombridae", ], 
                                      na.action = "na.fail", family = ordbeta)

abslat_model_hd_Lutjanidae <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale + 
                                        (1|Source) + (1|MarkerName),
                                      data = mtdna_small_hd[mtdna_small_hd$Family == "Lutjanidae", ], 
                                      na.action = "na.fail", family = ordbeta)

abslat_model_hd_Serranidae <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale + 
                                        (1|Source) + (1|MarkerName),
                                      data = mtdna_small_hd[mtdna_small_hd$Family == "Serranidae", ], 
                                      na.action = "na.fail", family = ordbeta)

abslat_model_hd_Pomacentridae <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale + 
                                           (1|Source) + (1|MarkerName),
                                         data = mtdna_small_hd[mtdna_small_hd$Family == "Pomacentridae", ], 
                                         na.action = "na.fail", family = ordbeta)

abslat_model_hd_Sebastidae <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale + 
                                        (1|Source) + (1|MarkerName),
                                      data = mtdna_small_hd[mtdna_small_hd$Family == "Sebastidae", ], 
                                      na.action = "na.fail", family = ordbeta)

abslat_model_hd_Engraulidae <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale + 
                                         (1|Source) + (1|MarkerName),
                                       data = mtdna_small_hd[mtdna_small_hd$Family == "Engraulidae", ], 
                                       na.action = "na.fail", family =  ordbeta)

abslat_model_hd_Gadidae <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale + 
                                     (1|Source) + (1|MarkerName),
                                   data = mtdna_small_hd[mtdna_small_hd$Family == "Gadidae", ], 
                                   na.action = "na.fail", family = ordbeta)

abslat_model_hd_Syngnathidae <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale + 
                                          (1|Source) + (1|MarkerName),
                                        data = mtdna_small_hd[mtdna_small_hd$Family == "Syngnathidae", ], 
                                        na.action = "na.fail", family = ordbeta)

abslat_model_hd_Rajidae <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale + 
                                     (1|Source) + (1|MarkerName),
                                   data = mtdna_small_hd[mtdna_small_hd$Family == "Rajidae", ], 
                                   na.action = "na.fail", family = ordbeta)

abslat_model_hd_Carcharhinidae <- glmmTMB(He ~ bp_scale + range_pos_scale + abslat_scale + 
                                            (1|Source) + (1|MarkerName),
                                          data = mtdna_small_hd[mtdna_small_hd$Family == "Carcharhinidae", ], 
                                          na.action = "na.fail", family = ordbeta)

#### Predict ####
#marginal effects
abslat_eff_Scombridae <- plot_model(abslat_model_hd_Scombridae, type = "pred", 
                                    terms = "abslat_scale [all]")
abslat_eff_Lutjanidae <- plot_model(abslat_model_hd_Lutjanidae, type = "pred", 
                                    terms = "abslat_scale [all]")
abslat_eff_Serranidae <- plot_model(abslat_model_hd_Serranidae, type = "pred", 
                                    terms = "abslat_scale [all]")
abslat_eff_Pomacentridae <- plot_model(abslat_model_hd_Pomacentridae, type = "pred", 
                                       terms = "abslat_scale [all]")
abslat_eff_Sebastidae <- plot_model(abslat_model_hd_Sebastidae, type = "pred", 
                                    terms = "abslat_scale [all]")
abslat_eff_Engraulidae <- plot_model(abslat_model_hd_Engraulidae, type = "pred", 
                                     terms = "abslat_scale [all]")
abslat_eff_Gadidae <- plot_model(abslat_model_hd_Gadidae, type = "pred", 
                                 terms = "abslat_scale [all]")
abslat_eff_Syngnathidae <- plot_model(abslat_model_hd_Syngnathidae, type = "pred", 
                                      terms = "abslat_scale [all]")
abslat_eff_Rajidae <- plot_model(abslat_model_hd_Rajidae, type = "pred", 
                                 terms = "abslat_scale [all]")
abslat_eff_Carcharhinidae <- plot_model(abslat_model_hd_Carcharhinidae, type = "pred", 
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
abslat_scale <- scale(mtdna_small_hd$abslat) #bc had to convert to numeric to run model/calculate marginal effects
abslat_eff_data_allfam$abslat <- (abslat_eff_data_allfam$x * attr(abslat_scale, "scaled:scale")) +
  attr(abslat_scale, "scaled:center")

#### Plot abslat ####

mtdna_hd_abslat_plot <- ggplot() +
  geom_line(data = abslat_eff_data_allfam,
            aes(x = abslat, y = predicted, color = Family), linewidth = 14) +
  annotate("text", x = 6, y = 0.97, label = "(b)", size = 100) +
  scale_color_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                                "#AA4499", "#44AA99", "#999933", "#882255", "#888888")) + 
  ylim(0.5, 1) + xlim(0, 90) +
  xlab("Absolute Latitude") + ylab(bquote(H[d]~"(mtDNA)")) + 
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
mtdna_hd_abslat_plot

###########################################################################################################

######### Lat figure #######

#lat models
lat_model_hd_Scombridae <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                                     I(lat_scale^2) + (1|Source) + (1|MarkerName),
                                   data = mtdna_small_hd[mtdna_small_hd$Family == "Scombridae", ], 
                                   na.action = "na.fail", family = ordbeta)

lat_model_hd_Lutjanidae <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                                     I(lat_scale^2) + (1|Source) + (1|MarkerName),
                                   data = mtdna_small_hd[mtdna_small_hd$Family == "Lutjanidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lat_model_hd_Serranidae <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                                     I(lat_scale^2) + (1|Source) + (1|MarkerName),
                                   data = mtdna_small_hd[mtdna_small_hd$Family == "Serranidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lat_model_hd_Pomacentridae <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                                        I(lat_scale^2) + (1|Source) + (1|MarkerName),
                                      data = mtdna_small_hd[mtdna_small_hd$Family == "Pomacentridae", ], 
                                      na.action = "na.fail", family = ordbeta)

lat_model_hd_Sebastidae <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                                     I(lat_scale^2) + (1|Source) + (1|MarkerName),
                                   data = mtdna_small_hd[mtdna_small_hd$Family == "Sebastidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lat_model_hd_Engraulidae <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                                      I(lat_scale^2) + (1|Source) + (1|MarkerName),
                                    data = mtdna_small_hd[mtdna_small_hd$Family == "Engraulidae", ], 
                                    na.action = "na.fail", family =  ordbeta)

lat_model_hd_Gadidae <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                                  I(lat_scale^2) + (1|Source) + (1|MarkerName),
                                data = mtdna_small_hd[mtdna_small_hd$Family == "Gadidae", ], 
                                na.action = "na.fail", family = ordbeta)

lat_model_hd_Syngnathidae <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                                       I(lat_scale^2) + (1|Source) + (1|MarkerName),
                                     data = mtdna_small_hd[mtdna_small_hd$Family == "Syngnathidae", ], 
                                     na.action = "na.fail", family = ordbeta)

lat_model_hd_Rajidae <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                                  I(lat_scale^2) + (1|Source) + (1|MarkerName),
                                data = mtdna_small_hd[mtdna_small_hd$Family == "Rajidae", ], 
                                na.action = "na.fail", family = ordbeta)

lat_model_hd_Carcharhinidae <- glmmTMB(He ~ bp_scale + range_pos_scale + lat_scale + 
                                         I(lat_scale^2) + (1|Source) + (1|MarkerName),
                                       data = mtdna_small_hd[mtdna_small_hd$Family == "Carcharhinidae", ], 
                                       na.action = "na.fail", family = ordbeta)

#### Predict ####
#marginal effects
lat_eff_Scombridae <- plot_model(lat_model_hd_Scombridae, type = "pred", 
                                 terms = "lat_scale [all]")
lat_eff_Lutjanidae <- plot_model(lat_model_hd_Lutjanidae, type = "pred", 
                                 terms = "lat_scale [all]")
lat_eff_Serranidae <- plot_model(lat_model_hd_Serranidae, type = "pred", 
                                 terms = "lat_scale [all]")
lat_eff_Pomacentridae <- plot_model(lat_model_hd_Pomacentridae, type = "pred", 
                                    terms = "lat_scale [all]")
lat_eff_Sebastidae <- plot_model(lat_model_hd_Sebastidae, type = "pred", 
                                 terms = "lat_scale [all]")
lat_eff_Engraulidae <- plot_model(lat_model_hd_Engraulidae, type = "pred", 
                                  terms = "lat_scale [all]")
lat_eff_Gadidae <- plot_model(lat_model_hd_Gadidae, type = "pred", 
                              terms = "lat_scale [all]")
lat_eff_Syngnathidae <- plot_model(lat_model_hd_Syngnathidae, type = "pred", 
                                   terms = "lat_scale [all]")
lat_eff_Rajidae <- plot_model(lat_model_hd_Rajidae, type = "pred", 
                              terms = "lat_scale [all]")
lat_eff_Carcharhinidae <- plot_model(lat_model_hd_Carcharhinidae, type = "pred", 
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
lat_scale <- scale(mtdna_small_hd$lat) #bc had to convert to numeric to run model/calculate marginal effects
  lat_eff_data_allfam$lat <- (lat_eff_data_allfam$x * attr(lat_scale, "scaled:scale")) +
    attr(lat_scale, "scaled:center")

#### Plot lat ####

mtdna_hd_lat_plot <- ggplot() +
  geom_line(data = lat_eff_data_allfam,
            aes(x = lat, y = predicted, color = Family), linewidth = 14) +
annotate("text", x = -70, y = 0.97, label = "(b)", size = 100) +
  scale_color_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                                "#AA4499", "#44AA99", "#999933", "#882255", "#888888")) + 
  ylim(0.5, 1) + xlim(-90, 90) +
  xlab("Latitude") + ylab(bquote(H[d]~"(mtDNA)")) + 
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
mtdna_hd_lat_plot

###########################################################################################################

######### Lon figure #######

#lon models
lon_model_hd_Scombridae <- glmmTMB(He ~ bp_scale + range_pos_scale + bs(lon_scale) + 
                                     (1|Source) + (1|MarkerName),
                                   data = mtdna_small_hd[mtdna_small_hd$Family == "Scombridae", ], 
                                   na.action = "na.fail", family = ordbeta)

lon_model_hd_Lutjanidae <- glmmTMB(He ~ bp_scale + range_pos_scale + bs(lon_scale) + 
                                     (1|Source) + (1|MarkerName),
                                   data = mtdna_small_hd[mtdna_small_hd$Family == "Lutjanidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lon_model_hd_Serranidae <- glmmTMB(He ~ bp_scale + range_pos_scale + bs(lon_scale) + 
                                     (1|Source) + (1|MarkerName),
                                   data = mtdna_small_hd[mtdna_small_hd$Family == "Serranidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lon_model_hd_Pomacentridae <- glmmTMB(He ~ bp_scale + range_pos_scale + bs(lon_scale) + 
                                        (1|Source) + (1|MarkerName),
                                      data = mtdna_small_hd[mtdna_small_hd$Family == "Pomacentridae", ], 
                                      na.action = "na.fail", family = ordbeta)

lon_model_hd_Sebastidae <- glmmTMB(He ~ bp_scale + range_pos_scale + bs(lon_scale) + 
                                     (1|Source) + (1|MarkerName),
                                   data = mtdna_small_hd[mtdna_small_hd$Family == "Sebastidae", ], 
                                   na.action = "na.fail", family = ordbeta)

lon_model_hd_Engraulidae <- glmmTMB(He ~ bp_scale + range_pos_scale + bs(lon_scale) + 
                                      (1|Source) + (1|MarkerName),
                                    data = mtdna_small_hd[mtdna_small_hd$Family == "Engraulidae", ], 
                                    na.action = "na.fail", family =  ordbeta)

lon_model_hd_Gadidae <- glmmTMB(He ~ bp_scale + range_pos_scale + bs(lon_scale) + 
                                  (1|Source) + (1|MarkerName),
                                data = mtdna_small_hd[mtdna_small_hd$Family == "Gadidae", ], 
                                na.action = "na.fail", family = ordbeta)

lon_model_hd_Syngnathidae <- glmmTMB(He ~ bp_scale + range_pos_scale + bs(lon_scale) + 
                                       (1|Source) + (1|MarkerName),
                                     data = mtdna_small_hd[mtdna_small_hd$Family == "Syngnathidae", ], 
                                     na.action = "na.fail", family = ordbeta)

lon_model_hd_Rajidae <- glmmTMB(He ~ bp_scale + range_pos_scale + bs(lon_scale) + 
                                  (1|Source) + (1|MarkerName),
                                data = mtdna_small_hd[mtdna_small_hd$Family == "Rajidae", ], 
                                na.action = "na.fail", family = ordbeta)

lon_model_hd_Carcharhinidae <- glmmTMB(He ~ bp_scale + range_pos_scale + bs(lon_scale) + 
                                         (1|Source) + (1|MarkerName),
                                       data = mtdna_small_hd[mtdna_small_hd$Family == "Carcharhinidae", ], 
                                       na.action = "na.fail", family = ordbeta)

#### Predict ####
#marginal effects
lon_eff_Scombridae <- plot_model(lon_model_hd_Scombridae, type = "pred", 
                                 terms = "lon_scale [all]")
lon_eff_Lutjanidae <- plot_model(lon_model_hd_Lutjanidae, type = "pred", 
                                 terms = "lon_scale [all]")
lon_eff_Serranidae <- plot_model(lon_model_hd_Serranidae, type = "pred", 
                                 terms = "lon_scale [all]")
lon_eff_Pomacentridae <- plot_model(lon_model_hd_Pomacentridae, type = "pred", 
                                    terms = "lon_scale [all]")
lon_eff_Sebastidae <- plot_model(lon_model_hd_Sebastidae, type = "pred", 
                                 terms = "lon_scale [all]")
lon_eff_Engraulidae <- plot_model(lon_model_hd_Engraulidae, type = "pred", 
                                  terms = "lon_scale [all]")
lon_eff_Gadidae <- plot_model(lon_model_hd_Gadidae, type = "pred", 
                              terms = "lon_scale [all]")
lon_eff_Syngnathidae <- plot_model(lon_model_hd_Syngnathidae, type = "pred", 
                                   terms = "lon_scale [all]")
lon_eff_Rajidae <- plot_model(lon_model_hd_Rajidae, type = "pred", 
                              terms = "lon_scale [all]")
lon_eff_Carcharhinidae <- plot_model(lon_model_hd_Carcharhinidae, type = "pred", 
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

#fix Syngnathidae dataframe (didn't make conf intervals??)
lon_eff_data_Syngnathidae$std.error <- "NA"
  lon_eff_data_Syngnathidae$conf.low <- "NA"
  lon_eff_data_Syngnathidae$conf.high <- "NA"

lon_eff_data_Syngnathidae <- lon_eff_data_Syngnathidae[, c("x", "predicted", 
                                                           "std.error", "conf.low",
                                                           "conf.high", "group", 
                                                           "group_col", "Family")]
  
#merge together
lon_eff_data_allfam <- rbind(lon_eff_data_Scombridae, lon_eff_data_Lutjanidae,
                             lon_eff_data_Serranidae, lon_eff_data_Pomacentridae,
                             lon_eff_data_Sebastidae, lon_eff_data_Engraulidae,
                             lon_eff_data_Gadidae, lon_eff_data_Syngnathidae,
                             lon_eff_data_Rajidae, lon_eff_data_Carcharhinidae)

#use same scaled:center & scaled:scale from original data
lon_scale <- scale(mtdna_small_hd$lon) #bc had to convert to numeric to run model/calculate marginal effects
lon_eff_data_allfam$lon <- (lon_eff_data_allfam$x * attr(lon_scale, "scaled:scale")) +
  attr(lon_scale, "scaled:center")

#### Plot lon ####

mtdna_hd_lon_plot <- ggplot() +
  annotate("rect", xmin = 95, xmax = 165, ymin = 0.5, ymax = 1, #adding box highlighting coral triangle
           fill = "darkolivegreen", alpha = 0.4) + 
  geom_line(data = lon_eff_data_allfam,
            aes(x = lon, y = predicted, color = Family), linewidth = 14) +
  annotate("text", x = -175, y = 0.97, label = "(b)", size = 100) +
  scale_color_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                                "#AA4499", "#44AA99", "#999933", "#882255", "#888888")) + 
  ylim(0.5, 1) + xlim(-180, 180) +
  xlab("Longitude") + ylab(bquote(H[d]~"(mtDNA)")) + 
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

#######################################################################################################

######### SSTmean figure #######

#sstmean models
sstmean_model_hd_Scombridae <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                         (1|Source) + (1|MarkerName),
                                       data = mtdna_small_hd[mtdna_small_hd$Family == "Scombridae", ], 
                                       na.action = "na.fail", family = ordbeta)

sstmean_model_hd_Lutjanidae <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                         (1|Source) + (1|MarkerName),
                                       data = mtdna_small_hd[mtdna_small_hd$Family == "Lutjanidae", ], 
                                       na.action = "na.fail", family = ordbeta)

sstmean_model_hd_Serranidae <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                         (1|Source) + (1|MarkerName),
                                       data = mtdna_small_hd[mtdna_small_hd$Family == "Serranidae", ], 
                                       na.action = "na.fail", family = ordbeta)

sstmean_model_hd_Pomacentridae <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                            (1|Source) + (1|MarkerName),
                                          data = mtdna_small_hd[mtdna_small_hd$Family == "Pomacentridae", ], 
                                          na.action = "na.fail", family = ordbeta)

sstmean_model_hd_Sebastidae <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                        (1|Source) + (1|MarkerName),
                                       data = mtdna_small_hd[mtdna_small_hd$Family == "Sebastidae", ], 
                                       na.action = "na.fail", family = ordbeta)

sstmean_model_hd_Engraulidae <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                          (1|Source) + (1|MarkerName),
                                        data = mtdna_small_hd[mtdna_small_hd$Family == "Engraulidae", ], 
                                        na.action = "na.fail", family =  ordbeta)

sstmean_model_hd_Gadidae <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                      (1|Source) + (1|MarkerName),
                                    data = mtdna_small_hd[mtdna_small_hd$Family == "Gadidae", ], 
                                    na.action = "na.fail", family = ordbeta)

sstmean_model_hd_Syngnathidae <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                           (1|Source) + (1|MarkerName),
                                         data = mtdna_small_hd[mtdna_small_hd$Family == "Syngnathidae", ], 
                                         na.action = "na.fail", family = ordbeta)

sstmean_model_hd_Rajidae <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                      (1|Source) + (1|MarkerName),
                                    data = mtdna_small_hd[mtdna_small_hd$Family == "Rajidae", ], 
                                    na.action = "na.fail", family = ordbeta)

sstmean_model_hd_Carcharhinidae <- glmmTMB(He ~ bp_scale + range_pos_scale + sstmean_scale + 
                                             (1|Source) + (1|MarkerName),
                                           data = mtdna_small_hd[mtdna_small_hd$Family == "Carcharhinidae", ], 
                                           na.action = "na.fail", family = ordbeta)

#### Predict ####
#marginal effects
sstmean_eff_Scombridae <- plot_model(sstmean_model_hd_Scombridae, type = "pred", 
                                     terms = "sstmean_scale [all]")
sstmean_eff_Lutjanidae <- plot_model(sstmean_model_hd_Lutjanidae, type = "pred", 
                                     terms = "sstmean_scale [all]")
sstmean_eff_Serranidae <- plot_model(sstmean_model_hd_Serranidae, type = "pred", 
                                     terms = "sstmean_scale [all]")
sstmean_eff_Pomacentridae <- plot_model(sstmean_model_hd_Pomacentridae, type = "pred", 
                                        terms = "sstmean_scale [all]")
sstmean_eff_Sebastidae <- plot_model(sstmean_model_hd_Sebastidae, type = "pred", 
                                     terms = "sstmean_scale [all]")
sstmean_eff_Engraulidae <- plot_model(sstmean_model_hd_Engraulidae, type = "pred", 
                                      terms = "sstmean_scale [all]")
sstmean_eff_Gadidae <- plot_model(sstmean_model_hd_Gadidae, type = "pred", 
                                  terms = "sstmean_scale [all]")
sstmean_eff_Syngnathidae <- plot_model(sstmean_model_hd_Syngnathidae, type = "pred", 
                                       terms = "sstmean_scale [all]")
sstmean_eff_Rajidae <- plot_model(sstmean_model_hd_Rajidae, type = "pred", 
                                  terms = "sstmean_scale [all]")
sstmean_eff_Carcharhinidae <- plot_model(sstmean_model_hd_Carcharhinidae, type = "pred", 
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

#merge together
sstmean_eff_data_allfam <- rbind(sstmean_eff_data_Scombridae, sstmean_eff_data_Lutjanidae,
                                 sstmean_eff_data_Serranidae, sstmean_eff_data_Pomacentridae,
                                 sstmean_eff_data_Sebastidae, sstmean_eff_data_Engraulidae,
                                 sstmean_eff_data_Gadidae, sstmean_eff_data_Syngnathidae,
                                 sstmean_eff_data_Rajidae, sstmean_eff_data_Carcharhinidae)

#use same scaled:center & scaled:scale from original data
sstmean_scale <- scale(mtdna_small_hd$sst.BO_sstmean) #bc had to convert to numeric to run model/calculate marginal effects
  sstmean_eff_data_allfam$sstmean <- (sstmean_eff_data_allfam$x * attr(sstmean_scale, "scaled:scale")) +
    attr(sstmean_scale, "scaled:center")

#### Plot SSTmean ####

mtdna_hd_sstmean_plot <- ggplot() +
  geom_line(data = sstmean_eff_data_allfam,
            aes(x = sstmean, y = predicted, color = Family), linewidth = 14) +
  annotate("text", x = 2, y = 0.985, label = "(b)", size = 100) + 
  scale_color_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                                "#AA4499", "#44AA99", "#999933", "#882255", "#888888")) + 
  ylim(0.5, 1) + xlim(0, 30) +
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

######### chloro mean figure #######

#chlomean models
chlomean_model_hd_Scombridae <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean + 
                                          I(logchlomean^2) + (1|Source) + (1|MarkerName),
                                        data = mtdna_small_hd[mtdna_small_hd$Family == "Scombridae", ], 
                                        na.action = "na.fail", family = ordbeta)

chlomean_model_hd_Lutjanidae <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean + 
                                          I(logchlomean^2) + (1|Source) + (1|MarkerName),
                                        data = mtdna_small_hd[mtdna_small_hd$Family == "Lutjanidae", ], 
                                        na.action = "na.fail", family = ordbeta)

chlomean_model_hd_Serranidae <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean + 
                                          I(logchlomean^2) + (1|Source) + (1|MarkerName),
                                        data = mtdna_small_hd[mtdna_small_hd$Family == "Serranidae", ], 
                                        na.action = "na.fail", family = ordbeta)

chlomean_model_hd_Pomacentridae <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean + 
                                             I(logchlomean^2) + (1|Source) + (1|MarkerName),
                                           data = mtdna_small_hd[mtdna_small_hd$Family == "Pomacentridae", ], 
                                           na.action = "na.fail", family = ordbeta)

chlomean_model_hd_Sebastidae <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean + 
                                          I(logchlomean^2) + (1|Source) + (1|MarkerName),
                                        data = mtdna_small_hd[mtdna_small_hd$Family == "Sebastidae", ], 
                                        na.action = "na.fail", family = ordbeta)

chlomean_model_hd_Engraulidae <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean + 
                                           I(logchlomean^2) + (1|Source) + (1|MarkerName),
                                         data = mtdna_small_hd[mtdna_small_hd$Family == "Engraulidae", ], 
                                         na.action = "na.fail", family =  ordbeta)

chlomean_model_hd_Gadidae <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean + 
                                       I(logchlomean^2) + (1|Source) + (1|MarkerName),
                                     data = mtdna_small_hd[mtdna_small_hd$Family == "Gadidae", ], 
                                     na.action = "na.fail", family = ordbeta)

chlomean_model_hd_Syngnathidae <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean + 
                                            I(logchlomean^2) + (1|Source) + (1|MarkerName),
                                          data = mtdna_small_hd[mtdna_small_hd$Family == "Syngnathidae", ], 
                                          na.action = "na.fail", family = ordbeta)

chlomean_model_hd_Rajidae <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean + 
                                       I(logchlomean^2) + (1|Source) + (1|MarkerName),
                                     data = mtdna_small_hd[mtdna_small_hd$Family == "Rajidae", ], 
                                     na.action = "na.fail", family = ordbeta)

chlomean_model_hd_Carcharhinidae <- glmmTMB(He ~ bp_scale + range_pos_scale + logchlomean + 
                                              I(logchlomean^2) + (1|Source) + (1|MarkerName),
                                            data = mtdna_small_hd[mtdna_small_hd$Family == "Carcharhinidae", ], 
                                            na.action = "na.fail", family = ordbeta)

#### Predict ####
#marginal effects
chlomean_eff_Scombridae <- plot_model(chlomean_model_hd_Scombridae, type = "pred", 
                                      terms = "logchlomean [all]")
chlomean_eff_Lutjanidae <- plot_model(chlomean_model_hd_Lutjanidae, type = "pred", 
                                      terms = "logchlomean [all]")
chlomean_eff_Serranidae <- plot_model(chlomean_model_hd_Serranidae, type = "pred", 
                                      terms = "logchlomean [all]")
chlomean_eff_Pomacentridae <- plot_model(chlomean_model_hd_Pomacentridae, type = "pred", 
                                         terms = "logchlomean [all]")
chlomean_eff_Sebastidae <- plot_model(chlomean_model_hd_Sebastidae, type = "pred", 
                                      terms = "logchlomean [all]")
chlomean_eff_Engraulidae <- plot_model(chlomean_model_hd_Engraulidae, type = "pred", 
                                       terms = "logchlomean [all]")
chlomean_eff_Gadidae <- plot_model(chlomean_model_hd_Gadidae, type = "pred", 
                                   terms = "logchlomean [all]")
chlomean_eff_Syngnathidae <- plot_model(chlomean_model_hd_Syngnathidae, type = "pred", 
                                        terms = "logchlomean [all]")
chlomean_eff_Rajidae <- plot_model(chlomean_model_hd_Rajidae, type = "pred", 
                                   terms = "logchlomean [all]")
chlomean_eff_Carcharhinidae <- plot_model(chlomean_model_hd_Carcharhinidae, type = "pred", 
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

#fix Sebastidae dataframe (didn't make conf intervals??)
chlomean_eff_data_Sebastidae$std.error <- "NA"
  chlomean_eff_data_Sebastidae$conf.low <- "NA"
  chlomean_eff_data_Sebastidae$conf.high <- "NA"

chlomean_eff_data_Sebastidae <- chlomean_eff_data_Sebastidae[, c("x", "predicted", 
                                                                 "std.error", "conf.low",
                                                                 "conf.high", "group", 
                                                                 "group_col", "Family")]

#fix Carcharhinidae dataframe (didn't make conf intervals??)
chlomean_eff_data_Carcharhinidae$std.error <- "NA"
  chlomean_eff_data_Carcharhinidae$conf.low <- "NA"
  chlomean_eff_data_Carcharhinidae$conf.high <- "NA"

chlomean_eff_data_Carcharhinidae <- chlomean_eff_data_Carcharhinidae[, c("x", "predicted", 
                                                                         "std.error", "conf.low",
                                                                         "conf.high", "group", 
                                                                         "group_col", "Family")]
  
#merge together
chlomean_eff_data_allfam <- rbind(chlomean_eff_data_Scombridae, chlomean_eff_data_Lutjanidae,
                                  chlomean_eff_data_Serranidae, chlomean_eff_data_Pomacentridae,
                                  chlomean_eff_data_Sebastidae, chlomean_eff_data_Engraulidae,
                                  chlomean_eff_data_Gadidae, chlomean_eff_data_Syngnathidae,
                                  chlomean_eff_data_Rajidae, chlomean_eff_data_Carcharhinidae)

#unlog chlorophyll mean
chlomean_eff_data_allfam$chlomean <- 10^(chlomean_eff_data_allfam$x)

#### Plot chlomean ####

mtdna_hd_chlomean_plot <- ggplot() +
  geom_line(data = chlomean_eff_data_allfam,
            aes(x = chlomean, y = predicted, color = Family), linewidth = 14) +
  annotate("text", x = 0.135, y = 0.985, label = "(e)", size = 100) + 
  scale_color_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                                "#AA4499", "#44AA99", "#999933", "#882255", "#888888")) + 
  ylim(0.4, 1.0) +
  scale_x_continuous(trans = "log10", limits = c(0.1, 10)) +
  xlab(bquote("Mean Chlorophyll"~(mg/m^3))) + ylab(bquote(H[d]~"(mtDNA)")) + 
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
mtdna_hd_chlomean_plot

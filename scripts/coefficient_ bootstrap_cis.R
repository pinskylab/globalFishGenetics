###################################### Script for Creating Coefficient 95% CIs  #######################################################

#read in data from the bootstrapping scripts
#calculate 95% CIs (where need be)
#create figures summarizing (supplemental figures for manuscript)

##########################################################################################################################################

######## Set-up ########
remove(list = ls())

#load libraries
library(here) #v.1.0.1
library(tidyverse) #v.2.0.0

#read in mtDNA data
pi_cis <- read.csv(here("output", "cis_pi.csv"))
  pi_cis <- pi_cis[,-1] #drop first --> same as fourth column
  pi_cis$metric <- c("pi")
hd_cis_env <- read.csv(here("output", "cis_hd_env.csv"))
  hd_cis_env <- hd_cis_env[,-1]
  hd_cis_env$metric <- c("hd")
hd_cis_geo <- read.csv(here("output", "cis_hd_geo.csv"))
  hd_cis_geo <- hd_cis_geo[,-1]
  hd_cis_geo$metric <- c("hd")

#read in msat data
he_null_cis <- read.csv(here("output", "boot_he_null.csv"))
he_abslat_cis <- read.csv(here("output", "boot_he_abslat.csv"))
he_lat_cis <- read.csv(here("output", "boot_he_lat.csv"))
he_lon_cis <- read.csv(here("output", "boot_he_lon.csv"))
he_sstmean_cis <- read.csv(here("output", "boot_he_sstmean.csv"))
he_chlomean_cis <- read.csv(here("output", "boot_he_chloromean.csv"))

##########################################################################################################################################

######## Clean up data ########

## create he dataframe (bc haven't calculated CIs here yet) ##
he_cis <- as.data.frame(c(rep("he", 11)))
  he_cis$X2.5ci <- c(quantile(he_null_cis$range_position, c(0.025)), #calculate quantiles for each model df
                     quantile(he_null_cis$CrossSpp, c(0.025)),
                     quantile(he_lat_cis$lat_scale, c(0.025)), 
                     quantile(he_lat_cis$I.lat_scale.2., c(0.025)), 
                     quantile(he_abslat_cis$abslat_scale, c(0.025)), 
                     quantile(he_lon_cis$bs.lon_scale.1, c(0.025)), 
                     quantile(he_lon_cis$bs.lon_scale.2, c(0.025)), 
                     quantile(he_lon_cis$bs.lon_scale.3, c(0.025)), 
                     quantile(he_sstmean_cis$sstmean_scale, c(0.025)), 
                     quantile(he_chlomean_cis$logchlomean, c(0.025)), 
                     quantile(he_chlomean_cis$I.logchlomean.2, c(0.025)))
  he_cis$X97.5ci <- c(quantile(he_null_cis$range_position, c(0.975)), 
                     quantile(he_null_cis$CrossSpp, c(0.975)),
                     quantile(he_lat_cis$lat_scale, c(0.975)), 
                     quantile(he_lat_cis$I.lat_scale.2., c(0.975)), 
                     quantile(he_abslat_cis$abslat_scale, c(0.975)), 
                     quantile(he_lon_cis$bs.lon_scale.1, c(0.975)), 
                     quantile(he_lon_cis$bs.lon_scale.2, c(0.975)), 
                     quantile(he_lon_cis$bs.lon_scale.3, c(0.975)), 
                     quantile(he_sstmean_cis$sstmean_scale, c(0.975)), 
                     quantile(he_chlomean_cis$logchlomean, c(0.975)), 
                     quantile(he_chlomean_cis$I.logchlomean.2, c(0.975)))
  he_cis$fixef <- c("range_position", "CrossSpp", "lat_scale", "I(lat_scale^2)", "abslat_scale",
                    "bs(lon_scale)1", "bs(lon_scale)2", "bs(lon_scale)3", "sst.BO_sstmean", 
                    "logchlomean", "I(logchlomean^2)")
  he_cis$model <- c("null", "null", "lat", "lat", "abslat", "lon", "lon", "lon",
                    "sstmean", "chloromean", "chloromean")
  colnames(he_cis) <- c("metric", "X2.5ci", "X97.5ci", "fixef", "model")

#merge all dataframes together  
all_cis <- rbind(pi_cis, hd_cis_geo, hd_cis_env, he_cis)

##################################################################################################################################

######## lat lon comparisons ########

#remove extra information
cis_latlon <- subset(all_cis, all_cis$model == "lat" | all_cis$model == "abslat"| 
                       all_cis$model == "lon")
  cis_latlon <- subset(cis_latlon, cis_latlon$fixef != "(Intercept)" & 
                         cis_latlon$fixef != "sigma" & cis_latlon$fixef != "range_pos_scale" & 
                         cis_latlon$fixef != "bp_scale" & cis_latlon$fixef != "range_position")

#add means (coefficients from "real" model)
#pi -> hd -> he: lat, lat^2, abslat, bslon1, bslon2, bslon3
cis_latlon$mean <- c(-0.038, -0.014, -0.045, 0.140, 0.282, 0.081, 
                     -0.144, -0.024, -0.103, -0.122, 1.160, 0.081, 
                     -0.053, -0.011, -0.030, 0.093, 0.003, 0.082)
 
## plot ## 
cis_latlon$fixef <- factor(cis_latlon$fixef, level = c("bs(lon_scale)3", "bs(lon_scale)2", "bs(lon_scale)1", 
                                                       "I(lat_scale^2)", "lat_scale", "abslat_scale"))

latlon_cis_plot <- ggplot(data = cis_latlon) + 
  geom_hline(aes(yintercept = 0), 
             color = "black", linewidth = 2, linetype = "dashed") +
  geom_point(aes(x = fixef, y = mean, color = metric), 
             position = position_dodge(width = 0.5), shape = "circle", size = 6) + 
  geom_errorbar(aes(x = fixef, ymin = X2.5ci, ymax = X97.5ci, color = metric), 
                position = position_dodge(width = 0.5), width = 0, linewidth = 2) + 
  scale_color_manual(values = c("#332288", "#88CCEE", "#117733")) +
  xlab("") + ylab("Coefficient") +
  scale_x_discrete(labels = c("b-spline lon 3", "b-spline lon 2", "b-spline lon 1", 
                              bquote(latitude^2), "latitude", "abs latitude")) +
  annotate("text", x = 6.25, y = -0.35, label = "(a)", size = 20) +
  coord_flip() +
  theme_minimal() + 
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 34),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.text = element_text(size = 34, color = "black"),
        legend.position = c(0.9, 0.9), 
        legend.text = element_text(size = 34),
        legend.title = element_blank(), 
        plot.margin = unit(c(1,1,1,1), "cm"),)
latlon_cis_plot

##################################################################################################################################

######## SST & chlorophyll comparisons ########

#remove extra information
cis_env <- subset(all_cis, all_cis$model == "sstmean" | all_cis$model == "chloromean")
  cis_env <- subset(cis_env, cis_env$fixef != "(Intercept)" & 
                      cis_env$fixef != "sigma" & cis_env$fixef != "bp_scale" & 
                      cis_env$fixef != "range_position" & cis_env$fixef != "range_pos_scale")
  
#add means (coefficients from "real" model)
#pi -> hd -> he: sstmean, chlomean, chlomean^2
cis_env$mean <- c(0.042, 0.022, -0.044, 
                  0.112, -0.022, -0.145, 
                  0.020, 0.027, -0.039)

## plot ##
cis_env$fixef <- factor(cis_env$fixef, level = c("I(logchlomean^2)", "logchlomean", "sstmean_scale"))

env_cis_plot <- ggplot(data = cis_env) + 
  geom_hline(aes(yintercept = 0), 
             color = "black", linewidth = 2, linetype = "dashed") +
  geom_point(aes(x = fixef, y = mean, color = metric), 
             position = position_dodge(width = 0.5), shape = "circle", size = 6) + 
  geom_errorbar(aes(x = fixef, ymin = X2.5ci, ymax = X97.5ci, color = metric), 
                position = position_dodge(width = 0.5), width = 0, linewidth = 2) + 
  scale_color_manual(values = c("#332288", "#88CCEE", "#117733")) +
  xlab("") + ylab("Coefficient") +
  scale_x_discrete(labels = c(bquote("chloro"~mean^2), "chloro mean", "SST mean")) +  
  annotate("text", x = 4.3, y = -0.225, label = "(b)", size = 20) +
  coord_flip() +
  theme_minimal() + 
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title = element_text(size = 34),
        axis.ticks = element_line(color = "black", linewidth = 1),
        axis.text = element_text(size = 34, color = "black"),
        legend.position = c(0.9, 0.9), 
        legend.text = element_text(size = 34),
        legend.title = element_blank(), 
        plot.margin = unit(c(1,1,1,1), "cm"),)
env_cis_plot 


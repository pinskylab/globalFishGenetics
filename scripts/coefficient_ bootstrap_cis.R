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
he_sstrange_cis <- read.csv(here("output", "boot_he_sstrange.csv"))
he_sstmax_cis <- read.csv(here("output", "boot_he_sstmax.csv"))
he_sstmin_cis <- read.csv(here("output", "boot_he_sstmin.csv"))
he_chlomean_cis <- read.csv(here("output", "boot_he_chloromean.csv"))
he_chlorange_cis <- read.csv(here("output", "boot_he_chlororange.csv"))
he_chlomax_cis <- read.csv(here("output", "boot_he_chloromax.csv"))
he_chlomin_cis <- read.csv(here("output", "boot_he_chloromin.csv"))

##########################################################################################################################################

######## Clean up data ########

## create he dataframe (bc haven't calculated CIs here yet) ##
he_cis <- as.data.frame(c(rep("he", 20)))
  he_cis$X2.5ci <- c(quantile(he_null_cis$range_position, c(0.025)), #calculate quantiles for each model df
                     quantile(he_null_cis$CrossSpp, c(0.025)),
                     quantile(he_lat_cis$lat_scale, c(0.025)), 
                     quantile(he_lat_cis$I.lat_scale.2., c(0.025)), 
                     quantile(he_abslat_cis$abslat_scale, c(0.025)), 
                     quantile(he_lon_cis$bs.lon_scale.1, c(0.025)), 
                     quantile(he_lon_cis$bs.lon_scale.2, c(0.025)), 
                     quantile(he_lon_cis$bs.lon_scale.3, c(0.025)), 
                     quantile(he_sstmean_cis$sst.BO_sstmean, c(0.025)), 
                     quantile(he_sstrange_cis$sst.BO_sstrange, c(0.025)), 
                     quantile(he_sstmax_cis$sst.BO_sstmax, c(0.025)), 
                     quantile(he_sstmin_cis$sst.BO_sstmin, c(0.025)), 
                     quantile(he_chlomean_cis$logchlomean, c(0.025)), 
                     quantile(he_chlomean_cis$I.logchlomean.2, c(0.025)),
                     quantile(he_chlorange_cis$logchlorange, c(0.025)), 
                     quantile(he_chlorange_cis$I.logchlorange.2, c(0.025)),
                     quantile(he_chlomax_cis$logchlomax, c(0.025)), 
                     quantile(he_chlomax_cis$I.logchlomax.2, c(0.025)),
                     quantile(he_chlomin_cis$logchlomin, c(0.025)), 
                     quantile(he_chlomin_cis$I.logchlomin.2, c(0.025)))
  he_cis$X97.5ci <- c(quantile(he_null_cis$range_position, c(0.975)), 
                     quantile(he_null_cis$CrossSpp, c(0.975)),
                     quantile(he_lat_cis$lat_scale, c(0.975)), 
                     quantile(he_lat_cis$I.lat_scale.2., c(0.975)), 
                     quantile(he_abslat_cis$abslat_scale, c(0.975)), 
                     quantile(he_lon_cis$bs.lon_scale.1, c(0.975)), 
                     quantile(he_lon_cis$bs.lon_scale.2, c(0.975)), 
                     quantile(he_lon_cis$bs.lon_scale.3, c(0.975)), 
                     quantile(he_sstmean_cis$sst.BO_sstmean, c(0.975)), 
                     quantile(he_sstrange_cis$sst.BO_sstrange, c(0.975)), 
                     quantile(he_sstmax_cis$sst.BO_sstmax, c(0.975)), 
                     quantile(he_sstmin_cis$sst.BO_sstmin, c(0.975)), 
                     quantile(he_chlomean_cis$logchlomean, c(0.975)), 
                     quantile(he_chlomean_cis$I.logchlomean.2, c(0.975)),
                     quantile(he_chlorange_cis$logchlorange, c(0.975)), 
                     quantile(he_chlorange_cis$I.logchlorange.2, c(0.975)),
                     quantile(he_chlomax_cis$logchlomax, c(0.975)), 
                     quantile(he_chlomax_cis$I.logchlomax.2, c(0.975)),
                     quantile(he_chlomin_cis$logchlomin, c(0.975)), 
                     quantile(he_chlomin_cis$I.logchlomin.2, c(0.975)))
  he_cis$fixef <- c("range_position", "CrossSpp", "lat_scale", "I(lat_scale^2)", "abslat_scale",
                    "bs(lon_scale)1", "bs(lon_scale)2", "bs(lon_scale)3", "sst.BO_sstmean", 
                    "sst.BO_sstrange", "sst.BO_sstmax", "sst.BO_sstmin",
                    "logchlomean", "I(logchlomean^2)", "logchlorange", "I(logchlorange^2)", 
                    "logchlomax", "I(logchlomax^2)", "logchlomin", "I(logchlomin^2)")
  he_cis$model <- c("null", "null", "lat", "lat", "lon", "lon", "lon", "abslat",
                    "sstmean", "sstrange", "sstmax", "sstmin", "chloromean", 
                    "chloromean", "chlororange", "chlororange", 
                    "chloromax", "chloromax", "chloromin", "chloromin")
  colnames(he_cis) <- c("metric", "X2.5ci", "X97.5ci", "fixef", "model")

#merge all dataframes together  
all_cis <- rbind(pi_cis, hd_cis_geo, hd_cis_env, he_cis)

##################################################################################################################################

######## lat lon comparisons ########

#remove extra information
cis_latlon <- subset(all_cis, all_cis$model == "lat" | all_cis$model == "abslat"| 
                       all_cis$model == "lon")
  cis_latlon <- subset(cis_latlon, cis_latlon$fixef != "(Intercept)" & 
                         cis_latlon$fixef != "sigma" & cis_latlon$fixef != "range_position" & 
                         cis_latlon$fixef != "bp_scale")

#add means (coefficients from "real" model)
#pi -> hd -> he: lat, lat^2, abslat, bslon1, bslon2, bslon3
cis_latlon$mean <- c(-0.038, -0.008, -0.046, 0.145, 0.310, 0.086, 
                     -0.270, 0.034, -0.158, 0.215, 1.500, 0.360, 
                     -0.024, -0.006, -0.020, 0.192, 0.128, 0.127)
 
## plot ## 
cis_latlon$fixef <- factor(cis_latlon$fixef, level = c("bs(lon_scale)3", "bs(lon_scale)2", "bs(lon_scale)1", 
                                                       "(lat_scale^2)", "lat_scale", "abslat_scale"))

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
  annotate("text", x = 6.25, y = -0.45, label = "A", size = 20) +
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

######## SST comparisons ########

#remove extra information
cis_sst <- subset(all_cis, all_cis$model == "sstmean" | all_cis$model == "sstrange"| 
                       all_cis$model == "sstmax" | all_cis$model == "sstmin")
  cis_sst <- subset(cis_sst, cis_sst$fixef != "(Intercept)" & 
                      cis_sst$fixef != "sigma" & cis_sst$fixef != "bp_scale" & 
                      cis_sst$fixef != "range_position")
  
#add means (coefficients from "real" model)
#pi -> hd -> he: mean, range, max, min
cis_sst$mean <- c(0.007, -0.003, 0.006, 0.005, 
                  0.031, -0.001, 0.030, 0.020,
                  0.0007, -0.005, -0.003, 0.003)

## plot ##
cis_sst$fixef <- factor(cis_sst$fixef, level = c("sst.BO_sstmin", "sst.BO_sstmax", 
                                                 "sst.BO_sstrange", "sst.BO_sstmean"))

sst_cis_plot <- ggplot(data = cis_sst) + 
  geom_hline(aes(yintercept = 0), 
             color = "black", linewidth = 2, linetype = "dashed") +
  geom_point(aes(x = fixef, y = mean, color = metric), 
             position = position_dodge(width = 0.5), shape = "circle", size = 6) + 
  geom_errorbar(aes(x = fixef, ymin = X2.5ci, ymax = X97.5ci, color = metric), 
                position = position_dodge(width = 0.5), width = 0, linewidth = 2) + 
  scale_color_manual(values = c("#332288", "#88CCEE", "#117733")) +
  xlab("") + ylab("Coefficient") +
  scale_x_discrete(labels = c("SST min", "SST max", "SST range", "SST mean")) +  
  annotate("text", x = 4.35, y = -0.0135, label = "B", size = 20) +
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
sst_cis_plot 

#################################################################################################################################

######## Chlorophyll  comparisons ########

#remove extra information
cis_chlo <- subset(all_cis, all_cis$model == "chloromean" | all_cis$model == "chlororange"| 
                     all_cis$model == "chloromax" | all_cis$model == "chloromin")
  cis_chlo <- subset(cis_chlo, cis_chlo$fixef != "(Intercept)" & 
                          cis_chlo$fixef != "sigma" & cis_chlo$fixef != "bp_scale" & 
                       cis_chlo$fixef != "range_position")

#fix wrong model designation
cis_chlo$model[cis_chlo$fixef == "logchlomin"] <- "chloromin"
  cis_chlo$model[cis_chlo$fixef == "I(logchlomin^2)"] <- "chloromin"

#add means (coefficients from "real" model)
#pi -> hd -> he: mean, mean^2, range, range^2, max, max^2, min, min^2
cis_chlo$mean <- c(0.024, -0.054, 0.019, -0.036, 0.045, -0.050, -0.001, -0.046,
                   0.099, -0.355, 0.041, -0.203, 0.214, -0.300, -0.078, -0.350, 
                   0.018, -0.065, 0.013, -0.031, 0.044, -0.051, -0.024, -0.055)
  
  
## plot ##
cis_chlo$fixef <- factor(cis_chlo$fixef, level = c("I(logchlomin^2)", "logchlomin", "I(logchlomax^2)",
                                                    "logchlomax", "I(logchlorange^2)", "logchlorange", 
                                                    "I(logchlomean^2)", "logchlomean"))

chlo_cis_plot <- ggplot(data = cis_chlo) + 
  geom_hline(aes(yintercept = 0), color = "black", linewidth = 2, linetype = "dashed") +
  geom_point(aes(x = fixef, y = mean, color = metric), 
             position = position_dodge(width = 0.5), shape = "circle", size = 6) + 
  geom_errorbar(aes(x = fixef, ymin = X2.5ci, ymax = X97.5ci, color = metric), 
                position = position_dodge(width = 0.5), width = 0, linewidth = 2) + 
  scale_color_manual(values = c("#332288", "#88CCEE", "#117733")) +
  xlab("") + ylab("Coefficient") +
  scale_x_discrete(labels = c(bquote("chloro"~min^2), "chloro min", 
                              bquote("chloro"~max^2), "chloro max",
                              bquote("chloro"~range^2), "chloro range", 
                              bquote("chloro"~mean^2), "chloro mean")) +
  annotate("text", x = 8.1, y = -0.425, label = "C", size = 20) +
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
chlo_cis_plot 

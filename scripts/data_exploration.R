#Script for cleaning/investigating datasets

remove(list = ls())

library(tidyverse)
library(gridExtra)
library(plotrix)
library(here)
library(splines)

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)
msat <- read.csv("output/msatloci_assembled.csv", stringsAsFactors = FALSE)
mtdna_env <- read.csv("output/mtdna_env.csv", stringsAsFactors = FALSE)
msat_env <- read.csv("output/msat_env.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)

#merge dataframes
mtdna <- merge(mtdna, mtdna_env[, c('X', 'sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin', 
                                    'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                    'chloroA.BO_chlomax', 'chloroA.BO_chlomin')], all.x = TRUE)
msat <- cbind(msat[, -1], msat_env[, c('sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin', 
                                       'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange', 
                                       'chloroA.BO_chlomax', 'chloroA.BO_chlomin')]) #merge not working for some reason, cbind bc in same order
mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                  'Centroid')], all.x = TRUE)
msat <- merge(msat, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                'Centroid')], all.x = TRUE)

#mtdna <- read.csv(here("output", "mtdna_assembled.csv"), stringsAsFactors = FALSE)

mtdna$abslat <- abs(mtdna$lat)
mtdna$abslon <- abs(mtdna$lon)

mtdna_nonapi <- subset(mtdna, mtdna$Pi != "NA")
mtdna_nonapi_small <- subset(mtdna_nonapi, mtdna_nonapi$Pi < 0.2)

mtdna_nonahe <- subset(mtdna, mtdna$He != "NA")
mtdna_nonahe$success <- round(mtdna_nonahe$He*mtdna_nonahe$n)
mtdna_nonahe$failure <- round((1 - mtdna_nonahe$He)*mtdna_nonahe$n)

#pi is linear regression model for now
mtdna_nonapi_small_lat_plot <- ggplot(data = mtdna_nonapi_small, aes(x = lat, y = Pi)) + geom_point(color = "darkgrey", size = 3, alpha = 0.5) + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "black", size = 2) + xlab("latitude") + ylab("mtdna pi")
mtdna_nonapi_small_lat_plot_annotated <- mtdna_nonapi_small_lat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdna_nonapi_small_lat_plot_annotated

mtdna_nonapi_small_plot <- ggplot(data = mtdna_nonapi_small, aes(x = abslat, y = Pi)) + geom_point(color = "darkgrey", size = 3, alpha = 0.5) + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "black", size = 2) + xlab("absolute latitude") + ylab("mtdna pi")
mtdna_nonapi_small_plot_annotated <- mtdna_nonapi_small_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdna_nonapi_small_plot_annotated
lm(mtdna_nonapi_small$Pi ~ mtdna_nonapi_small$abslat)

#binomial regression (changes relationship from linear)
mtdna_nonahe_lat_plot <- ggplot(data = mtdna_nonahe, aes(x = lat, y = He, succ = success, fail = failure)) + 
  geom_point(color = "darkgrey", size = 3, alpha = 0.5) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), color = "black", size = 2) + 
  xlab("latitude") + ylab("mtdna hd")
mtdna_nonahe_lat_plot_annotated <- mtdna_nonahe_lat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdna_nonahe_lat_plot_annotated

mtdna_nonahe_plot <- ggplot(data = mtdna_nonahe, aes(x = abslat, y = He, succ = success, fail = failure)) + 
  geom_point(color = "darkgrey", size = 3, alpha = 0.5) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), color = "black", size = 2) + 
  xlab("absolute latitude") + ylab("mtdna hd")
mtdna_nonahe_plot_annotated <- mtdna_nonahe_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdna_nonahe_plot_annotated
lm(mtdna_nonahe$He ~ mtdna_nonahe$abslat)

#msat stuff
#msat <- read.csv("output/msatloci_assembled.csv", stringsAsFactors = FALSE)

msat$abslat <- abs(msat$lat)
msat$abslon <- abs(msat$lon)
msat_nonan <- subset(msat, msat$n != "NA" & msat$n != "25-32")
  msat_nonan$n <- as.numeric(msat_nonan$n)
  class(msat_nonan$n) #check that is numeric
msat_nonan$success <- round(msat_nonan$He*msat_nonan$n)
msat_nonan$failure <- round((1 - msat_nonan$He)*msat_nonan$n)

msat_lat_plot <- ggplot(data = msat_nonan, aes(x = lat, y = He, succ = success, fail = failure)) + 
  geom_point(color = "darkgrey", size = 3, alpha = 0.5) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), color = "black", size = 2) + 
  xlab("latitude") + ylab("msat he")
msat_lat_plot_annotated <- msat_lat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
msat_lat_plot_annotated

msat_plot <- ggplot(data = msat_nonan, aes(x = abslat, y = He, succ = success, fail = failure)) + 
  geom_point(color = "darkgrey", size = 3, alpha = 0.5) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), color = "black", size = 2) + 
  xlab("absolute latitude") + ylab("msat he")
msat_plot_annotated <- msat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
msat_plot_annotated
lm(msat$He ~ msat$abslat)

div_lat_all_plot <- grid.arrange(mtdna_nonapi_small_lat_plot_annotated, mtdna_nonahe_lat_plot_annotated, 
                                    msat_onelocus_lat_plot_annotated, ncol = 3)
div_lat_all_plot

div_abslat_all_plot <- grid.arrange(mtdna_nonapi_small_plot_annotated, mtdna_nonahe_plot_annotated, 
                                 msat_onelocus_plot_annotated, ncol = 3)
div_abslat_all_plot

########################################################################################################

#lon plots

#linear for now
mtdna_nonapi_small_lon_plot <- ggplot(data = mtdna_nonapi_small, aes(x = lon, y = Pi)) + geom_point(color = "darkgrey", size = 3, alpha = 0.5) + 
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "black", size = 2) + xlab("longitude") + ylab("mtdna pi")
mtdna_nonapi_small_lon_plot_annotated <- mtdna_nonapi_small_lon_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdna_nonapi_small_lon_plot_annotated
lm(mtdna_nonapi_small$Pi ~ mtdna_nonapi_small$lon)

#binomial mtdna He v. lon
mtdna_nonahe_lon_plot <- ggplot(data = mtdna_nonahe, aes(x = lon, y = He, succ = success, fail = failure)) + 
  geom_point(color = "darkgrey", size = 3, alpha = 0.5) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), color = "black", size = 2) + 
  xlab("longitude") + ylab("mtdna hd")
mtdna_nonahe_lon_plot_annotated <- mtdna_nonahe_lon_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdna_nonahe_lon_plot_annotated
lm(mtdna_nonahe$He ~ mtdna_nonahe$lon)

#msat stuff

#binomial msat v. lon
msat_lon_plot <- ggplot(data = msat_nonan, aes(x = lon, y = He, succ = success, fail = failure)) + 
  geom_point(color = "darkgrey", size = 3, alpha = 0.5) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), color = "black", size = 2) + 
  xlab("longitude") + ylab("msat he")
msat_lon_plot_annotated <- msat_lon_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
msat_lon_plot_annotated
lm(msat$He ~ msat$lon)


div_lon_all_plot <- grid.arrange(mtdna_nonapi_small_lon_plot_annotated, mtdna_nonahe_lon_plot_annotated, 
                                 msat_onelocus_lon_plot_annotated, ncol = 3)
div_lon_all_plot

#bin observations into 10 deg latitudinal bands (and longitudinal bands) and plot mean (just for loop/lapply function OR round down and plot)
#Manel et al. --> used rgeos package to do equal area plotting and help with creating heatmap
  #need pick spatial scale to create equal area grid cells, then use lat/lon coordinates to put observations into grid cells and calculate means within grid cells
    #then re-plot with heat map

########################################################################################################

#fit splines

#msat
msat_abslat_spline_plot <- ggplot(data = msat, aes(x = abslat, y = He, color = Pelagic_Coastal)) + 
  geom_point(size = 3, alpha = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2) + 
  xlab("absolute latitude") + ylab("msat he") + scale_color_manual(values = c("#325168", "#7da0c0"))
msat_abslat_spline_plot_annotated <- msat_abslat_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
msat_abslat_spline_plot_annotated

msat_lat_spline_plot <- ggplot(data = msat, aes(x = lat, y = He, color = Pelagic_Coastal)) + 
  geom_point(size = 3, alpha = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2) + 
  xlab("latitude") + ylab("msat he") + scale_color_manual(values = c("#325168", "#7da0c0"))
msat_lat_spline_plot_annotated <- msat_lat_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
msat_lat_spline_plot_annotated

msat_lon_spline_plot <- ggplot(data = msat, aes(x = lon, y = He, color = Pelagic_Coastal)) + 
  geom_point(size = 3, alpha = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2) + 
  xlab("longitude") + ylab("msat he") + scale_color_manual(values = c("#325168", "#7da0c0"))
msat_lon_spline_plot_annotated <- msat_lon_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
msat_lon_spline_plot_annotated

#mtdna hd
mtdnahd_abslat_spline_plot <- ggplot(data = mtdna_nonahe, aes(x = abslat, y = He, color = Pelagic_Coastal)) + 
  geom_point(size = 3, alpha = 0.15) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2) + 
  xlab("absolute latitude") + ylab("mtdna hd") + scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnahd_abslat_spline_plot_annotated <- mtdnahd_abslat_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnahd_abslat_spline_plot_annotated

mtdnahd_lat_spline_plot <- ggplot(data = mtdna_nonahe, aes(x = lat, y = He, color = Pelagic_Coastal)) + 
  geom_point(size = 3, alpha = 0.15) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2) + 
  xlab("latitude") + ylab("mtdna hd") + ylim(0, 1.0) + scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnahd_lat_spline_plot_annotated <- mtdnahd_lat_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnahd_lat_spline_plot_annotated

mtdnahd_lon_spline_plot <- ggplot(data = mtdna_nonahe, aes(x = lon, y = He, color = Pelagic_Coastal)) + 
  geom_point(size = 3, alpha = 0.15) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2) + 
  xlab("longitude") + ylab("mtdna hd") + scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnahd_lon_spline_plot_annotated <- mtdnahd_lon_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnahd_lon_spline_plot_annotated

#mtdna pi
mtdnapi_abslat_spline_plot <- ggplot(data = mtdna_nonapi_small, aes(x = abslat, y = Pi, color = Pelagic_Coastal)) + 
  geom_point(size = 3, alpha = 0.15) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2) + 
  xlab("absolute latitude") + ylab("mtdna pi") + ylim(0, 0.05) + scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnapi_abslat_spline_plot_annotated <- mtdnapi_abslat_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnapi_abslat_spline_plot_annotated

mtdnapi_lat_spline_plot <- ggplot(data = mtdna_nonapi_small, aes(x = lat, y = Pi, color = Pelagic_Coastal)) + 
  geom_point(size = 3, alpha = 0.15) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2) + 
  xlab("latitude") + ylab("mtdna pi") + ylim(0, 0.05) + scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnapi_lat_spline_plot_annotated <- mtdnapi_lat_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnapi_lat_spline_plot_annotated

mtdnapi_lon_spline_plot <- ggplot(data = mtdna_nonapi_small, aes(x = lon, y = Pi, color = Pelagic_Coastal)) + 
  geom_point(size = 3, alpha = 0.15) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2) + 
  xlab("longitude") + ylab("mtdna pi") + ylim(0, 0.05) + scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnapi_lon_spline_plot_annotated <- mtdnapi_lon_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnapi_lon_spline_plot_annotated

#env sst mean
msat_sstmean_spline_plot <- ggplot(data = msat, aes(x = sst.BO_sstmean, y = He)) + 
  geom_point(color = "#7da0c0", size = 3, alpha = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2, color = "black") + 
  xlab("SST mean") + ylab("msat he") #+ scale_color_manual(values = c("#325168", "#7da0c0"))
msat_sstmean_spline_plot_annotated <- msat_sstmean_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
msat_sstmean_spline_plot_annotated

mtdnahe_sstmean_spline_plot <- ggplot(data = mtdna_nonahe, aes(x = sst.BO_sstmean, y = He)) + 
  geom_point(color = "#7da0c0", size = 3, alpha = 0.15) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2, color = "black") + 
  xlab("SST mean") + ylab("mtdna hd") #+ scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnahe_sstmean_spline_plot_annotated <- mtdnahe_sstmean_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnahe_sstmean_spline_plot_annotated

mtdnapi_sstmean_spline_plot <- ggplot(data = mtdna_nonapi_small, aes(x = sst.BO_sstmean, y = Pi)) + 
  geom_point(color = "#7da0c0", size = 3, alpha = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2, color = "black") + 
  xlab("SST mean") + ylab("mtdna pi") + ylim(0, 0.05) #+ scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnapi_sstmean_spline_plot_annotated <- mtdnapi_sstmean_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnapi_sstmean_spline_plot_annotated

#env DO
msat_DO_spline_plot <- ggplot(data = msat, aes(x = BO_dissox, y = He)) + 
  geom_point(color = "#7da0c0", size = 3, alpha = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2, color = "black") + 
  xlab("DO mean") + ylab("msat he")# + scale_color_manual(values = c("#325168", "#7da0c0"))
msat_DO_spline_plot_annotated <- msat_DO_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
msat_DO_spline_plot_annotated

mtdnahe_DO_spline_plot <- ggplot(data = mtdna_nonahe, aes(x = BO_dissox, y = He)) + 
  geom_point(color = "#7da0c0", size = 3, alpha = 0.15) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2, color = "black") + 
  xlab("DO mean") + ylab("mtdna hd") #+ scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnahe_DO_spline_plot_annotated <- mtdnahe_DO_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnahe_DO_spline_plot_annotated

mtdnapi_DO_spline_plot <- ggplot(data = mtdna_nonapi_small, aes(x = BO_dissox, y = Pi)) + 
  geom_point(color = "#7da0c0", size = 3, alpha = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2, color = "black") + 
  xlab("DO mean") + ylab("mtdna pi") + ylim(0, 0.05) #+ scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnapi_DO_spline_plot_annotated <- mtdnapi_DO_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnapi_DO_spline_plot_annotated

#env chloroA
msat_CA_spline_plot <- ggplot(data = msat, aes(x = chloroA.BO_chlomean, y = He)) + 
  geom_point(color = "#7da0c0", size = 3, alpha = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2, color = "black") + 
  xlab("chlorophyll A mean") + ylab("msat he")# + scale_color_manual(values = c("#325168", "#7da0c0"))
msat_CA_spline_plot_annotated <- msat_CA_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
msat_CA_spline_plot_annotated

mtdnahe_CA_spline_plot <- ggplot(data = mtdna_nonahe, aes(x = chloroA.BO_chlomean, y = He)) + 
  geom_point(color = "#7da0c0", size = 3, alpha = 0.15) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2, color = "black") + 
  xlab("chlorophyll A mean") + ylab("mtdna hd") #+ scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnahe_CA_spline_plot_annotated <- mtdnahe_CA_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnahe_CA_spline_plot_annotated

mtdnapi_CA_spline_plot <- ggplot(data = mtdna_nonapi_small, aes(x = chloroA.BO_chlomean, y = Pi)) + 
  geom_point(color = "#7da0c0", size = 3, alpha = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x), size = 2, color = "black") + 
  xlab("chlorophyll A mean") + ylab("mtdna pi") + ylim(0, 0.05) #+ scale_color_manual(values = c("#325168", "#7da0c0"))
mtdnapi_CA_spline_plot_annotated <- mtdnapi_CA_spline_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
mtdnapi_CA_spline_plot_annotated

##########################################################################################
#msat lat 10 degree binning

msat$lat_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)
msat$abslat_round <- NA

library(DescTools)

for(i in 1:nrow(msat)) {
  cat(paste(i, " ", sep = ''))
  {msat$lat_round[i] <- DescTools::RoundTo(msat$lat[i], multiple = 10, FUN = floor)}
}

for(i in 1:nrow(msat)) {
  cat(paste(i, " ", sep = ''))
  {msat$abslat_round[i] <- DescTools::RoundTo(msat$abslat[i], multiple = 10, FUN = floor)}
}
msat$abslat_round <- as.factor(msat$abslat_round)
msat$lat_round <- as.factor(msat$lat_round)

#subset by lat bands and then get mean and plot that

#abslat bands
msat_abslat0 <- subset(msat, abslat_round == 0)
msat_abslat0_mean <- mean(msat_abslat0$He)
msat_abslat0_sd <- std.error(msat_abslat0$He)
msat_abslat10 <- subset(msat, abslat_round == 10)
msat_abslat10_mean <- mean(msat_abslat10$He)
msat_abslat10_sd <- std.error(msat_abslat10$He)
msat_abslat20 <- subset(msat, abslat_round == 20)
msat_abslat20_mean <- mean(msat_abslat20$He)
msat_abslat20_sd <- std.error(msat_abslat20$He)
msat_abslat30 <- subset(msat, abslat_round == 30)
msat_abslat30_mean <- mean(msat_abslat30$He)
msat_abslat30_sd <- std.error(msat_abslat30$He)
msat_abslat40 <- subset(msat, abslat_round == 40)
msat_abslat40_mean <- mean(msat_abslat40$He)
msat_abslat40_sd <- std.error(msat_abslat40$He)
msat_abslat50 <- subset(msat, abslat_round == 50)
msat_abslat50_mean <- mean(msat_abslat50$He)
msat_abslat50_sd <- std.error(msat_abslat50$He)
msat_abslat60 <- subset(msat, abslat_round == 60)
msat_abslat60_mean <- mean(msat_abslat60$He)
msat_abslat60_sd <- std.error(msat_abslat60$He)
msat_abslat70 <- subset(msat, abslat_round == 70)
msat_abslat70_mean <- mean(msat_abslat70$He)
msat_abslat70_sd <- std.error(msat_abslat70$He)
msat_abslat80 <- subset(msat, abslat_round == 80)
msat_abslat80_mean <- mean(msat_abslat80$He)
msat_abslat80_sd <- std.error(msat_abslat80$He)

abslat_He_mean <- c(msat_abslat0_mean, msat_abslat10_mean, msat_abslat20_mean, msat_abslat30_mean, 
                 msat_abslat40_mean, msat_abslat50_mean, msat_abslat60_mean, msat_abslat70_mean, 
                 msat_abslat80_mean)
abslat_He_SD <- c(msat_abslat0_sd, msat_abslat10_sd, msat_abslat20_sd, msat_abslat30_sd, 
                    msat_abslat40_sd, msat_abslat50_sd, msat_abslat60_sd, msat_abslat70_sd, 
                    msat_abslat80_sd)
abslat_bin <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
abslat_He_bin_df <- data.frame(abslat_bin, abslat_He_mean, abslat_He_SD)
  colnames(abslat_He_bin_df) <- c("abslat", "He", "SE") #SD is really SE

msat_abslat_plot <- ggplot(data = abslat_He_bin_df, aes(x = abslat, y = He)) + 
  geom_line(size = 1.5, color = "grey") + 
  geom_point(size = 4) + 
  xlab("absolute latitude") + ylab("msat He") + 
  geom_errorbar(data = abslat_He_bin_df, aes(x = abslat, ymin = He-SE, ymax = He+SE), width = 4, size = 2)
  msat_abslat_plot_annotated <- msat_abslat_plot + theme_bw() + 
    scale_x_continuous(limits = c(-10, 80), breaks = c(0, 20, 40, 60, 80)) + 
    scale_y_continuous(limits = c(0, 1.0)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
msat_abslat_plot_annotated

msat_abslat_boxplot <- ggplot(data = msat, aes(y = He, x = abslat_round)) + 
  geom_jitter(shape = 16, position=position_jitter(0.2), alpha = 0.25, color ="#7da0c0") + 
  geom_boxplot(aes(y = He, x = abslat_round), lwd = 1.5, color = "#09263a")
msat_abslat_boxplot_annotated <- msat_abslat_boxplot + theme_bw() + xlab("absolute latitude") + ylab("msat He") +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 22, color = "black"), axis.title = element_text(size = 25), 
        legend.position = "none")
msat_abslat_boxplot_annotated

#lat bands
#msat_latneg90 <- subset(msat, lat_round == -90)
#msat_latneg90_mean <- mean(msat_latneg90$He)
#msat_latneg90_sd <- std.error(msat_latneg90$He)
#msat_latneg80 <- subset(msat, lat_round == -80)
#msat_latneg80_mean <- mean(msat_latneg80$He)
#msat_latneg80_sd <- std.error(msat_latneg80$He)
#msat_latneg70 <- subset(msat, lat_round == -70)
#msat_latneg70_mean <- mean(msat_latneg70$He)
#msat_latneg70_sd <- std.error(msat_latneg70$He)
#msat_latneg60 <- subset(msat, lat_round == -60)
#msat_latneg60_mean <- mean(msat_latneg60$He)
#msat_latneg60_sd <- std.error(msat_latneg60$He)
#msat_latneg50 <- subset(msat, lat_round == -50)
#msat_latneg50_mean <- mean(msat_latneg50$He)
#msat_latneg50_sd <- std.error(msat_latneg50$He)
#msat_latneg40 <- subset(msat, lat_round == -40)
#msat_latneg40_mean <- mean(msat_latneg40$He)
#msat_latneg40_sd <- std.error(msat_latneg40$He)
#msat_latneg30 <- subset(msat, lat_round == -30)
#msat_latneg30_mean <- mean(msat_latneg30$He)
#msat_latneg30_sd <- std.error(msat_latneg30$He)
#msat_latneg20 <- subset(msat, lat_round == -20)
#msat_latneg20_mean <- mean(msat_latneg20$He)
#msat_latneg20_sd <- std.error(msat_latneg20$He)
#msat_latneg10 <- subset(msat, lat_round == -10)
#msat_latneg10_mean <- mean(msat_latneg10$He)
#msat_latneg10_sd <- std.error(msat_latneg10$He)
#msat_lat0 <- subset(msat, lat_round == 0)
#msat_lat0_mean <- mean(msat_lat0$He)
#msat_lat0_sd <- std.error(msat_lat0$He)
#msat_lat10 <- subset(msat, lat_round == 10)
#msat_lat10_mean <- mean(msat_lat10$He)
#msat_lat10_sd <- std.error(msat_lat10$He)
#msat_lat20 <- subset(msat, lat_round == 20)
#msat_lat20_mean <- mean(msat_lat20$He)
#msat_lat20_sd <- std.error(msat_lat20$He)
#msat_lat30 <- subset(msat, lat_round == 30)
#msat_lat30_mean <- mean(msat_lat30$He)
#msat_lat30_sd <- std.error(msat_lat30$He)
#msat_lat40 <- subset(msat, lat_round == 40)
#msat_lat40_mean <- mean(msat_lat40$He)
#msat_lat40_sd <- std.error(msat_lat40$He)
#msat_lat50 <- subset(msat, lat_round == 50)
#msat_lat50_mean <- mean(msat_lat50$He)
#msat_lat50_sd <- std.error(msat_lat50$He)
#msat_lat60 <- subset(msat, lat_round == 60)
#msat_lat60_mean <- mean(msat_lat60$He)
#msat_lat60_sd <- std.error(msat_lat60$He)
#msat_lat70 <- subset(msat, lat_round == 70)
#msat_lat70_mean <- mean(msat_lat70$He)
#msat_lat70_sd <- std.error(msat_lat70$He)
#msat_lat80 <- subset(msat, lat_round == 80)
#msat_lat80_mean <- mean(msat_lat80$He)
#msat_lat80_sd <- std.error(msat_lat80$He)

#lat_He_mean <- c(msat_latneg90_mean, msat_latneg80_mean, msat_latneg70_mean, 
#                 msat_latneg60_mean, msat_latneg50_mean, msat_latneg40_mean, 
#                 msat_latneg30_mean, msat_latneg20_mean, msat_latneg10_mean, 
#                 msat_lat0_mean, msat_lat10_mean, msat_lat20_mean, msat_lat30_mean, 
#                 msat_lat40_mean, msat_lat50_mean, msat_lat60_mean, msat_lat70_mean, 
#                 msat_lat80_mean)
#lat_SD_mean <- c(msat_latneg90_sd, msat_latneg80_sd, msat_latneg70_sd, 
#                 msat_latneg60_sd, msat_latneg50_sd, msat_latneg40_sd, 
#                 msat_latneg30_sd, msat_latneg20_sd, msat_latneg10_sd, 
#                 msat_lat0_sd, msat_lat10_sd, msat_lat20_sd, msat_lat30_sd, 
#                 msat_lat40_sd, msat_lat50_sd, msat_lat60_sd, msat_lat70_sd, 
#                 msat_lat80_sd)
#lat_bin <- c(-90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80)
#lat_He_bin_df <- data.frame(lat_bin, lat_He_mean, lat_SD_mean)
#  colnames(lat_He_bin_df) <- c("lat", "He", "SD")

#msat_lat_bar_plot <- ggplot() + 
#  geom_bar(data = lat_He_bin_df, aes(x = lat, y = He), stat = "identity", color = "#00B8E5", fill = "#00B8E5") + xlab("latitude") + ylab("msat He") + 
#  #geom_smooth(data = msat_nonan, aes(x = lat, y = He, succ = success, fail = failure), method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x, color = "black", size = 2)
#  geom_errorbar(data = lat_He_bin_df, aes(x = lat, ymin = He-SD, ymax = He+SD), width = 4)
#msat_lat_bar_plot_annotated <- msat_lat_bar_plot + theme_bw() + coord_flip() + 
#  scale_x_continuous(limits = c(-100, 80), breaks = c(-80, -60, -40, -20, 0, 20, 40, 60, 80)) +
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
#msat_lat_bar_plot_annotated
#msat_lat_histogram <- ggplot(data = msat_nonan, aes(msat_nonan$lat)) + geom_histogram(binwidth = 10) + 
#  xlab("latitude") + ylab("count")
#msat_lat_histogram_annotated <- msat_lat_histogram + theme_bw() + 
#  scale_x_continuous(limits = c(-100, 80), breaks = c(-80, -60, -40, -20, 0, 20, 40, 60, 80)) + coord_flip() +
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
#msat_lat_histogram_annotated
#
#msat_lat_all_plot <- grid.arrange(msat_lat_bar_plot_annotated, msat_lat_histogram_annotated, ncol = 1)
#msat_lat_all_plot

msat_lat_boxplot <- ggplot(data = msat, aes(y = He, x = lat_round)) + 
  geom_jitter(shape = 16, position=position_jitter(0.2), alpha = 0.25, color ="#7da0c0") + 
  geom_boxplot(aes(y = He, x = lat_round), lwd = 1.5, color = "#09263a")
msat_lat_boxplot_annotated <- msat_lat_boxplot + theme_bw() + xlab("latitude") + ylab("msat He") +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 22, color = "black"), axis.title = element_text(size = 25), 
        legend.position = "none")
msat_lat_boxplot_annotated


#####################################################################################################

#mtdna he lat 10 binning

mtdna_nonahe$lat_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)
mtdna_nonahe$abslat_round <- NA

library(DescTools)

for(i in 1:nrow(mtdna_nonahe)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_nonahe$lat_round[i] <- DescTools::RoundTo(mtdna_nonahe$lat[i], multiple = 10, FUN = floor)}
}

for(i in 1:nrow(mtdna_nonahe)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_nonahe$abslat_round[i] <- DescTools::RoundTo(mtdna_nonahe$abslat[i], multiple = 10, FUN = floor)}
}

mtdna_nonahe$abslat_round <- as.factor(mtdna_nonahe$abslat_round)
mtdna_nonahe$lat_round <- as.factor(mtdna_nonahe$lat_round)

#subset by lat bands and then get mean and plot that

#abslat bands
#mtdna_abslat0 <- subset(mtdna_nonahe, abslat_round == 0)
#mtdna_abslat0_mean <- mean(mtdna_abslat0$He)
#mtdna_abslat0_se <- std.error(mtdna_abslat0$He)
#mtdna_abslat10 <- subset(mtdna_nonahe, abslat_round == 10)
#mtdna_abslat10_mean <- mean(mtdna_abslat10$He)
#mtdna_abslat10_se <- std.error(mtdna_abslat10$He)
#mtdna_abslat20 <- subset(mtdna_nonahe, abslat_round == 20)
#mtdna_abslat20_mean <- mean(mtdna_abslat20$He)
#mtdna_abslat20_se <- std.error(mtdna_abslat20$He)
#mtdna_abslat30 <- subset(mtdna_nonahe, abslat_round == 30)
#mtdna_abslat30_mean <- mean(mtdna_abslat30$He)
#mtdna_abslat30_se <- std.error(mtdna_abslat30$He)
#mtdna_abslat40 <- subset(mtdna_nonahe, abslat_round == 40)
#mtdna_abslat40_mean <- mean(mtdna_abslat40$He)
#mtdna_abslat40_se <- std.error(mtdna_abslat40$He)
#mtdna_abslat50 <- subset(mtdna_nonahe, abslat_round == 50)
#mtdna_abslat50_mean <- mean(mtdna_abslat50$He)
#mtdna_abslat50_se <- std.error(mtdna_abslat50$He)
#mtdna_abslat60 <- subset(mtdna_nonahe, abslat_round == 60)
#mtdna_abslat60_mean <- mean(mtdna_abslat60$He)
#mtdna_abslat60_se <- std.error(mtdna_abslat60$He)
#mtdna_abslat70 <- subset(mtdna_nonahe, abslat_round == 70)
#mtdna_abslat70_mean <- mean(mtdna_abslat70$He)
#mtdna_abslat70_se <- std.error(mtdna_abslat70$He)
#mtdna_abslat80 <- subset(mtdna_nonahe, abslat_round == 80)
#mtdna_abslat80_mean <- mean(mtdna_abslat80$He)
#mtdna_abslat80_se <- std.error(mtdna_abslat80$He)
#
#abslat_He_mean_mtdna <- c(mtdna_abslat0_mean, mtdna_abslat10_mean, mtdna_abslat20_mean, mtdna_abslat30_mean, 
#                    mtdna_abslat40_mean, mtdna_abslat50_mean, mtdna_abslat60_mean, mtdna_abslat70_mean, 
#                    mtdna_abslat80_mean)
#abslat_He_SE_mtdna <- c(mtdna_abslat0_se, mtdna_abslat10_se, mtdna_abslat20_se, mtdna_abslat30_se, 
#                          mtdna_abslat40_se, mtdna_abslat50_se, mtdna_abslat60_se, mtdna_abslat70_se, 
#                          mtdna_abslat80_se)
#abslat_bin <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
#abslat_He_bin_df_mtdna <- data.frame(abslat_bin, abslat_He_mean_mtdna, abslat_He_SE_mtdna)
#  colnames(abslat_He_bin_df_mtdna) <- c("abslat", "He", "SE")
#
#mtdna_abslat_bar_plot <- ggplot(data = abslat_He_bin_df_mtdna, aes(x = abslat, y = He, ymin = He-SE, ymax = He+SE)) + 
#  geom_bar(stat = "identity", color = "#00B8E5", fill = "#00B8E5") + xlab("absolute latitude") + ylab("mtdna Hd") + 
#  #geom_smooth(data = mtdna_nonahe, aes(x = abslat, y = He, succ = success, fail = failure), method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x, color = "black", size = 2)
#  geom_errorbar(width = 4)
#mtdna_abslat_bar_plot_annotated <- mtdna_abslat_bar_plot + theme_bw() + coord_flip() +
#  scale_x_continuous(limits = c(-10, 80), breaks = c(0, 20, 40, 60, 80)) + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
#mtdna_abslat_bar_plot_annotated
#mtdna_abslat_histogram <- ggplot(data = mtdna_nonahe, aes(mtdna_nonahe$abslat)) + geom_histogram(binwidth = 10) + 
#  xlab("absolute latitude") + ylab("count")
#mtdna_abslat_histogram_annotated <- mtdna_abslat_histogram + theme_bw() + 
#  scale_x_continuous(limits = c(-10, 80), breaks = c(0, 20, 40, 60, 80)) + coord_flip() + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
#mtdna_abslat_histogram_annotated
#
#mtdna_abslat_he_all_plot <- grid.arrange(mtdna_abslat_bar_plot_annotated, mtdna_abslat_histogram_annotated, ncol = 1)
#mtdna_abslat_he_all_plot

mtdnahe_abslat_boxplot <- ggplot(data = mtdna_nonahe, aes(y = He, x = abslat_round)) + 
  geom_boxplot(aes(y = He, x = abslat_round), lwd = 1.5, color = "#09263a") + 
  geom_jitter(shape = 16, position=position_jitter(0.2), alpha = 0.45, color ="#7da0c0") #+ 
 # geom_boxplot(aes(y = He, x = abslat_round), lwd = 1.5, color = "#09263a")
mtdnahe_abslat_boxplot_annotated <- mtdnahe_abslat_boxplot + theme_bw() + xlab("absolute latitude") + ylab("mtdna Hd") +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 22, color = "black"), axis.title = element_text(size = 25), 
        legend.position = "none")
mtdnahe_abslat_boxplot_annotated

#lat bands
#mtdna_latneg90 <- subset(mtdna_nonahe, lat_round == -90)
#mtdna_latneg90_mean <- mean(mtdna_latneg90$He)
#mtdna_latneg90_se <- std.error(mtdna_latneg90$He)
#mtdna_latneg80 <- subset(mtdna_nonahe, lat_round == -80)
#mtdna_latneg80_mean <- mean(mtdna_latneg80$He)
#mtdna_latneg80_se <- std.error(mtdna_latneg80$He)
#mtdna_latneg70 <- subset(mtdna_nonahe, lat_round == -70)
#mtdna_latneg70_mean <- mean(mtdna_latneg70$He)
#mtdna_latneg70_se <- std.error(mtdna_latneg70$He)
#mtdna_latneg60 <- subset(mtdna_nonahe, lat_round == -60)
#mtdna_latneg60_mean <- mean(mtdna_latneg60$He)
#mtdna_latneg60_se <- std.error(mtdna_latneg60$He)
#mtdna_latneg50 <- subset(mtdna_nonahe, lat_round == -50)
#mtdna_latneg50_mean <- mean(mtdna_latneg50$He)
#mtdna_latneg50_se <- std.error(mtdna_latneg50$He)
#mtdna_latneg40 <- subset(mtdna_nonahe, lat_round == -40)
#mtdna_latneg40_mean <- mean(mtdna_latneg40$He)
#mtdna_latneg40_se <- std.error(mtdna_latneg40$He)
#mtdna_latneg30 <- subset(mtdna_nonahe, lat_round == -30)
#mtdna_latneg30_mean <- mean(mtdna_latneg30$He)
#mtdna_latneg30_se <- std.error(mtdna_latneg30$He)
#mtdna_latneg20 <- subset(mtdna_nonahe, lat_round == -20)
#mtdna_latneg20_mean <- mean(mtdna_latneg20$He)
#mtdna_latneg20_se <- std.error(mtdna_latneg20$He)
#mtdna_latneg10 <- subset(mtdna_nonahe, lat_round == -10)
#mtdna_latneg10_mean <- mean(mtdna_latneg10$He)
#mtdna_latneg10_se <- std.error(mtdna_latneg10$He)
#mtdna_lat0 <- subset(mtdna_nonahe, lat_round == 0)
#mtdna_lat0_mean <- mean(mtdna_lat0$He)
#mtdna_lat0_se <- std.error(mtdna_lat0$He)
#mtdna_lat10 <- subset(mtdna_nonahe, lat_round == 10)
#mtdna_lat10_mean <- mean(mtdna_lat10$He)
#mtdna_lat10_se <- std.error(mtdna_lat10$He)
#mtdna_lat20 <- subset(mtdna_nonahe, lat_round == 20)
#mtdna_lat20_mean <- mean(mtdna_lat20$He)
#mtdna_lat20_se <- std.error(mtdna_lat20$He)
#mtdna_lat30 <- subset(mtdna_nonahe, lat_round == 30)
#mtdna_lat30_mean <- mean(mtdna_lat30$He)
#mtdna_lat30_se <- std.error(mtdna_lat30$He)
#mtdna_lat40 <- subset(mtdna_nonahe, lat_round == 40)
#mtdna_lat40_mean <- mean(mtdna_lat40$He)
#mtdna_lat40_se <- std.error(mtdna_lat40$He)
#mtdna_lat50 <- subset(mtdna_nonahe, lat_round == 50)
#mtdna_lat50_mean <- mean(mtdna_lat50$He)
#mtdna_lat50_se <- std.error(mtdna_lat50$He)
#mtdna_lat60 <- subset(mtdna_nonahe, lat_round == 60)
#mtdna_lat60_mean <- mean(mtdna_lat60$He)
#mtdna_lat60_se <- std.error(mtdna_lat60$He)
#mtdna_lat70 <- subset(mtdna_nonahe, lat_round == 70)
#mtdna_lat70_mean <- mean(mtdna_lat70$He)
#mtdna_lat70_se <- std.error(mtdna_lat70$He)
#mtdna_lat80 <- subset(mtdna_nonahe, lat_round == 80)
#mtdna_lat80_mean <- mean(mtdna_lat80$He)
#mtdna_lat80_se <- std.error(mtdna_lat80$He)
#
#lat_He_mean <- c(mtdna_latneg90_mean, mtdna_latneg80_mean, mtdna_latneg70_mean, 
#                 mtdna_latneg60_mean, mtdna_latneg50_mean, mtdna_latneg40_mean, 
#                 mtdna_latneg30_mean, mtdna_latneg20_mean, mtdna_latneg10_mean, 
#                 mtdna_lat0_mean, mtdna_lat10_mean, mtdna_lat20_mean, mtdna_lat30_mean, 
#                 mtdna_lat40_mean, mtdna_lat50_mean, mtdna_lat60_mean, mtdna_lat70_mean, 
#                 mtdna_lat80_mean)
#lat_He_SE <- c(mtdna_latneg90_se, mtdna_latneg80_se, mtdna_latneg70_se, 
#                 mtdna_latneg60_se, mtdna_latneg50_se, mtdna_latneg40_se, 
#                 mtdna_latneg30_se, mtdna_latneg20_se, mtdna_latneg10_se, 
#                 mtdna_lat0_se, mtdna_lat10_se, mtdna_lat20_se, mtdna_lat30_se, 
#                 mtdna_lat40_se, mtdna_lat50_se, mtdna_lat60_se, mtdna_lat70_se, 
#                 mtdna_lat80_se)
#lat_bin <- c(-90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80)
#lat_He_bin_df_mtdna <- data.frame(lat_bin, lat_He_mean, lat_He_SE)
#  colnames(lat_He_bin_df_mtdna) <- c("lat", "He", "SE")
#
#mtdna_lat_bar_plot <- ggplot(data = lat_He_bin_df_mtdna, aes(x = lat, y = He, ymin = He-SE, ymax = He+SE)) + 
#  geom_bar(stat = "identity", color = "#00B8E5", fill = "#00B8E5") + xlab("latitude") + ylab("mtdna Hd") + 
#  #geom_smooth(data = mtdna_nonahe, aes(x = lat, y = He, succ = success, fail = failure), method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x, color = "black", size = 2)
#  geom_errorbar(width = 4)
#mtdna_lat_bar_plot_annotated <- mtdna_lat_bar_plot + theme_bw() + coord_flip() +
#  scale_x_continuous(limits = c(-100, 80), breaks = c(-80, -60, -40, -20, 0, 20, 40, 60, 80)) + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
#mtdna_lat_bar_plot_annotated
#mtdna_lat_histogram <- ggplot(data = mtdna_nonahe, aes(mtdna_nonahe$lat)) + geom_histogram(binwidth = 10) + 
#  xlab("latitude") + ylab("count")
#mtdna_lat_histogram_annotated <- mtdna_lat_histogram + theme_bw() + 
#  scale_x_continuous(limits = c(-100, 80), breaks = c(-80, -60, -40, -20, 0, 20, 40, 60, 80)) + coord_flip() + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
#mtdna_lat_histogram_annotated
#
#mtdna_lat_he_all_plot <- grid.arrange(mtdna_lat_bar_plot_annotated, mtdna_lat_histogram_annotated, ncol = 1)
#mtdna_lat_he_all_plot

mtdnahe_lat_boxplot <- ggplot(data = mtdna_nonahe, aes(y = He, x = lat_round)) + 
  geom_boxplot(lwd = 1.5, color = "#09263a") + 
  geom_jitter(shape = 16, position=position_jitter(0.2), alpha = 0.45, color ="#7da0c0") #+ 
# geom_boxplot(aes(y = He, x = abslat_round), lwd = 1.5, color = "#09263a")
mtdnahe_lat_boxplot_annotated <- mtdnahe_lat_boxplot + theme_bw() + xlab("latitude") + ylab("mtdna Hd") +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 22, color = "black"), axis.title = element_text(size = 25), 
        legend.position = "none")
mtdnahe_lat_boxplot_annotated

#####################################################################################################

#mtdna pi lat 10 binning

mtdna_nonapi_small$lat_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)
mtdna_nonapi_small$abslat_round <- NA

library(DescTools)

for(i in 1:nrow(mtdna_nonapi_small)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_nonapi_small$lat_round[i] <- DescTools::RoundTo(mtdna_nonapi_small$lat[i], multiple = 10, FUN = floor)}
}

for(i in 1:nrow(mtdna_nonapi_small)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_nonapi_small$abslat_round[i] <- DescTools::RoundTo(mtdna_nonapi_small$abslat[i], multiple = 10, FUN = floor)}
}

mtdna_nonapi_small$abslat_round <- as.factor(mtdna_nonapi_small$abslat_round)
mtdna_nonapi_small$lat_round <- as.factor(mtdna_nonapi_small$lat_round)

#subset by lat bands and then get mean and plot that

#abslat bands
mtdna_pi_abslat0 <- subset(mtdna_nonapi_small, abslat_round == 0)
mtdna_pi_abslat0_mean <- mean(mtdna_pi_abslat0$Pi)
mtdna_pi_abslat0_se <- std.error(mtdna_pi_abslat0$Pi)
mtdna_pi_abslat10 <- subset(mtdna_nonapi_small, abslat_round == 10)
mtdna_pi_abslat10_mean <- mean(mtdna_pi_abslat10$Pi)
mtdna_pi_abslat10_se <- std.error(mtdna_pi_abslat10$Pi)
mtdna_pi_abslat20 <- subset(mtdna_nonapi_small, abslat_round == 20)
mtdna_pi_abslat20_mean <- mean(mtdna_pi_abslat20$Pi)
mtdna_pi_abslat20_se <- std.error(mtdna_pi_abslat20$Pi)
mtdna_pi_abslat30 <- subset(mtdna_nonapi_small, abslat_round == 30)
mtdna_pi_abslat30_mean <- mean(mtdna_pi_abslat30$Pi)
mtdna_pi_abslat30_se <- std.error(mtdna_pi_abslat30$Pi)
mtdna_pi_abslat40 <- subset(mtdna_nonapi_small, abslat_round == 40)
mtdna_pi_abslat40_mean <- mean(mtdna_pi_abslat40$Pi)
mtdna_pi_abslat40_se <- std.error(mtdna_pi_abslat40$Pi)
mtdna_pi_abslat50 <- subset(mtdna_nonapi_small, abslat_round == 50)
mtdna_pi_abslat50_mean <- mean(mtdna_pi_abslat50$Pi)
mtdna_pi_abslat50_se <- std.error(mtdna_pi_abslat50$Pi)
mtdna_pi_abslat60 <- subset(mtdna_nonapi_small, abslat_round == 60)
mtdna_pi_abslat60_mean <- mean(mtdna_pi_abslat60$Pi)
mtdna_pi_abslat60_se <- std.error(mtdna_pi_abslat60$Pi)
mtdna_pi_abslat70 <- subset(mtdna_nonapi_small, abslat_round == 70)
mtdna_pi_abslat70_mean <- mean(mtdna_pi_abslat70$Pi)
mtdna_pi_abslat70_se <- std.error(mtdna_pi_abslat70$Pi)
mtdna_pi_abslat80 <- subset(mtdna_nonapi_small, abslat_round == 80)
mtdna_pi_abslat80_mean <- mean(mtdna_pi_abslat80$Pi)
mtdna_pi_abslat80_se <- std.error(mtdna_pi_abslat80$Pi)
#
abslat_Pi_mean_mtdna <- c(mtdna_pi_abslat0_mean, mtdna_pi_abslat10_mean, mtdna_pi_abslat20_mean, mtdna_pi_abslat30_mean, 
                          mtdna_pi_abslat40_mean, mtdna_pi_abslat50_mean, mtdna_pi_abslat60_mean, mtdna_pi_abslat70_mean, 
                          mtdna_pi_abslat80_mean)
abslat_Pi_se_mtdna <- c(mtdna_pi_abslat0_se, mtdna_pi_abslat10_se, mtdna_pi_abslat20_se, mtdna_pi_abslat30_se, 
                          mtdna_pi_abslat40_se, mtdna_pi_abslat50_se, mtdna_pi_abslat60_se, mtdna_pi_abslat70_se, 
                          mtdna_pi_abslat80_se)
abslat_bin <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
abslat_Pi_bin_df_mtdna <- data.frame(abslat_bin, abslat_Pi_mean_mtdna, abslat_Pi_se_mtdna)
  colnames(abslat_Pi_bin_df_mtdna) <- c("abslat", "Pi", "SE")
#
#mtdna_pi_abslat_bar_plot <- ggplot(data = abslat_Pi_bin_df_mtdna, aes(x = abslat, y = Pi, ymin = Pi-SE, ymax = Pi+SE)) + 
#  geom_bar(stat = "identity", color = "#00B8E5", fill = "#00B8E5") + xlab("absolute latitude") + ylab("mtdna Pi") + 
#  #geom_smooth(data = mtdna_nonapi_small, aes(x = abslat, y = Pi), method = "lm", formula = y ~ x, color = "black", size = 2)
#  geom_errorbar(width = 4)
#mtdna_pi_abslat_bar_plot_annotated <- mtdna_pi_abslat_bar_plot + theme_bw() + coord_flip() + 
#  scale_x_continuous(limits = c(-10, 80), breaks = c(0, 20, 40, 60, 80)) + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
#mtdna_pi_abslat_bar_plot_annotated
#mtdna_pi_abslat_histogram <- ggplot(data = mtdna_nonapi_small, aes(mtdna_nonapi_small$abslat)) + geom_histogram(binwidth = 10) + 
#  xlab("absolute latitude") + ylab("count")
#mtdna_pi_abslat_histogram_annotated <- mtdna_pi_abslat_histogram + theme_bw() + 
#  scale_x_continuous(limits = c(-10, 80), breaks = c(0, 20, 40, 60, 80)) + coord_flip() + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
#mtdna_pi_abslat_histogram_annotated
#
#mtdna_abslat_pi_all_plot <- grid.arrange(mtdna_pi_abslat_bar_plot_annotated, mtdna_pi_abslat_histogram_annotated, ncol = 1)
#mtdna_abslat_pi_all_plot

mtdnapi_abslat_boxplot <- ggplot(data = mtdna_nonapi_small, aes(y = Pi, x = abslat_round)) + 
  geom_boxplot(lwd = 1.5, color = "#09263a") + 
  geom_jitter(shape = 16, position=position_jitter(0.2), alpha = 0.45, color ="#7da0c0") #+ 
# geom_boxplot(aes(y = He, x = abslat_round), lwd = 1.5, color = "#09263a")
mtdnapi_abslat_boxplot_annotated <- mtdnapi_abslat_boxplot + theme_bw() + xlab("absolute latitude") + ylab("mtdna pi") + ylim(0, 0.05) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 22, color = "black"), axis.title = element_text(size = 25), 
        legend.position = "none")
mtdnapi_abslat_boxplot_annotated

#lat bands
#mtdna_pi_latneg90 <- subset(mtdna_nonapi_small, lat_round == -90)
#mtdna_pi_latneg90_mean <- mean(mtdna_pi_latneg90$Pi)
#mtdna_pi_latneg90_se <- std.error(mtdna_pi_latneg90$Pi)
#mtdna_pi_latneg80 <- subset(mtdna_nonapi_small, lat_round == -80)
#mtdna_pi_latneg80_mean <- mean(mtdna_pi_latneg80$Pi)
#mtdna_pi_latneg80_se <- std.error(mtdna_pi_latneg80$Pi)
#mtdna_pi_latneg70 <- subset(mtdna_nonapi_small, lat_round == -70)
#mtdna_pi_latneg70_mean <- mean(mtdna_pi_latneg70$Pi)
#mtdna_pi_latneg70_se <- std.error(mtdna_pi_latneg70$Pi)
#mtdna_pi_latneg60 <- subset(mtdna_nonapi_small, lat_round == -60)
#mtdna_pi_latneg60_mean <- mean(mtdna_pi_latneg60$Pi)
#mtdna_pi_latneg60_se <- std.error(mtdna_pi_latneg60$Pi)
#mtdna_pi_latneg50 <- subset(mtdna_nonapi_small, lat_round == -50)
#mtdna_pi_latneg50_mean <- mean(mtdna_pi_latneg50$Pi)
#mtdna_pi_latneg50_se <- std.error(mtdna_pi_latneg50$Pi)
#mtdna_pi_latneg40 <- subset(mtdna_nonapi_small, lat_round == -40)
#mtdna_pi_latneg40_mean <- mean(mtdna_pi_latneg40$Pi)
#mtdna_pi_latneg40_se <- std.error(mtdna_pi_latneg40$Pi)
#mtdna_pi_latneg30 <- subset(mtdna_nonapi_small, lat_round == -30)
#mtdna_pi_latneg30_mean <- mean(mtdna_pi_latneg30$Pi)
#mtdna_pi_latneg30_se <- std.error(mtdna_pi_latneg30$Pi)
#mtdna_pi_latneg20 <- subset(mtdna_nonapi_small, lat_round == -20)
#mtdna_pi_latneg20_mean <- mean(mtdna_pi_latneg20$Pi)
#mtdna_pi_latneg20_se <- std.error(mtdna_pi_latneg20$Pi)
#mtdna_pi_latneg10 <- subset(mtdna_nonapi_small, lat_round == -10)
#mtdna_pi_latneg10_mean <- mean(mtdna_pi_latneg10$Pi)
#mtdna_pi_latneg10_se <- std.error(mtdna_pi_latneg10$Pi)
#mtdna_pi_lat0 <- subset(mtdna_nonapi_small, lat_round == 0)
#mtdna_pi_lat0_mean <- mean(mtdna_pi_lat0$Pi)
#mtdna_pi_lat0_se <- std.error(mtdna_pi_lat0$Pi)
#mtdna_pi_lat10 <- subset(mtdna_nonapi_small, lat_round == 10)
#mtdna_pi_lat10_mean <- mean(mtdna_pi_lat10$Pi)
#mtdna_pi_lat10_se <- std.error(mtdna_pi_lat10$Pi)
#mtdna_pi_lat20 <- subset(mtdna_nonapi_small, lat_round == 20)
#mtdna_pi_lat20_mean <- mean(mtdna_pi_lat20$Pi)
#mtdna_pi_lat20_se <- std.error(mtdna_pi_lat20$Pi)
#mtdna_pi_lat30 <- subset(mtdna_nonapi_small, lat_round == 30)
#mtdna_pi_lat30_mean <- mean(mtdna_pi_lat30$Pi)
#mtdna_pi_lat30_se <- std.error(mtdna_pi_lat30$Pi)
#mtdna_pi_lat40 <- subset(mtdna_nonapi_small, lat_round == 40)
#mtdna_pi_lat40_mean <- mean(mtdna_pi_lat40$Pi)
#mtdna_pi_lat40_se <- std.error(mtdna_pi_lat40$Pi)
#mtdna_pi_lat50 <- subset(mtdna_nonapi_small, lat_round == 50)
#mtdna_pi_lat50_mean <- mean(mtdna_pi_lat50$Pi)
#mtdna_pi_lat50_se <- std.error(mtdna_pi_lat50$Pi)
#mtdna_pi_lat60 <- subset(mtdna_nonapi_small, lat_round == 60)
#mtdna_pi_lat60_mean <- mean(mtdna_pi_lat60$Pi)
#mtdna_pi_lat60_se <- std.error(mtdna_pi_lat60$Pi)
#mtdna_pi_lat70 <- subset(mtdna_nonapi_small, lat_round == 70)
#mtdna_pi_lat70_mean <- mean(mtdna_pi_lat70$Pi)
#mtdna_pi_lat70_se <- std.error(mtdna_pi_lat70$Pi)
#mtdna_pi_lat80 <- subset(mtdna_nonapi_small, lat_round == 80)
#mtdna_pi_lat80_mean <- mean(mtdna_pi_lat80$Pi)
#mtdna_pi_lat80_se <- std.error(mtdna_pi_lat80$Pi)
#
#lat_Pi_mean <- c(mtdna_pi_latneg90_mean, mtdna_pi_latneg80_mean, mtdna_pi_latneg70_mean, 
#                 mtdna_pi_latneg60_mean, mtdna_pi_latneg50_mean, mtdna_pi_latneg40_mean, 
#                 mtdna_pi_latneg30_mean, mtdna_pi_latneg20_mean, mtdna_pi_latneg10_mean, 
#                 mtdna_pi_lat0_mean, mtdna_pi_lat10_mean, mtdna_pi_lat20_mean, mtdna_pi_lat30_mean, 
#                 mtdna_pi_lat40_mean, mtdna_pi_lat50_mean, mtdna_pi_lat60_mean, mtdna_pi_lat70_mean, 
#                 mtdna_pi_lat80_mean)
#lat_Pi_SE <- c(mtdna_pi_latneg90_se, mtdna_pi_latneg80_se, mtdna_pi_latneg70_se, 
#                 mtdna_pi_latneg60_se, mtdna_pi_latneg50_se, mtdna_pi_latneg40_se, 
#                 mtdna_pi_latneg30_se, mtdna_pi_latneg20_se, mtdna_pi_latneg10_se, 
#                 mtdna_pi_lat0_se, mtdna_pi_lat10_se, mtdna_pi_lat20_se, mtdna_pi_lat30_se, 
#                 mtdna_pi_lat40_se, mtdna_pi_lat50_se, mtdna_pi_lat60_se, mtdna_pi_lat70_se, 
#                 mtdna_pi_lat80_se)
#lat_bin <- c(-90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80)
#lat_Pi_bin_df_mtdna <- data.frame(lat_bin, lat_Pi_mean, lat_Pi_SE)
#  colnames(lat_Pi_bin_df_mtdna) <- c("lat", "Pi", "SE")
#
#mtdna_pi_lat_bar_plot <- ggplot(data = lat_Pi_bin_df_mtdna, aes(x = lat, y = Pi, ymin = Pi-SE, ymax = Pi+SE)) +  
#  geom_bar(stat = "identity", color = "#00B8E5", fill = "#00B8E5") + xlab("latitude") + ylab("mtdna Pi") + 
#  #geom_smooth(data = mtdna_nonapi_small, aes(x = lat, y = Pi), method = "lm", formula = y ~ x, color = "black", size = 2)
#  geom_errorbar(width = 4)
#mtdna_pi_lat_bar_plot_annotated <- mtdna_pi_lat_bar_plot + theme_bw() + coord_flip() + 
#  scale_x_continuous(limits = c(-100, 80), breaks = c(-80, -60, -40, -20, 0, 20, 40, 60, 80)) + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
#mtdna_pi_lat_bar_plot_annotated
#mtdna_pi_lat_histogram <- ggplot(data = mtdna_nonapi_small, aes(mtdna_nonapi_small$lat)) + geom_histogram(binwidth = 10) + 
#  xlab("latitude") + ylab("count")
#mtdna_pi_lat_histogram_annotated <- mtdna_pi_lat_histogram + theme_bw() + 
#  scale_x_continuous(limits = c(-100, 80), breaks = c(-80, -60, -40, -20, 0, 20, 40, 60, 80)) + coord_flip() +
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"))
#mtdna_pi_lat_histogram_annotated
#
#mtdna_lat_pi_all_plot <- grid.arrange(mtdna_pi_lat_bar_plot_annotated, mtdna_pi_lat_histogram_annotated, ncol = 1)
#mtdna_lat_pi_all_plot
#
#abslat_all_plot <- grid.arrange(mtdna_pi_abslat_bar_plot_annotated, mtdna_abslat_bar_plot_annotated, msat_abslat_bar_plot_annotated, 
#                                      ncol = 3)
#abslat_all_plot
#
#abslat_all_whist_plot <- grid.arrange(mtdna_pi_abslat_bar_plot_annotated, mtdna_abslat_bar_plot_annotated, msat_abslat_bar_plot_annotated, 
#                                mtdna_pi_abslat_histogram_annotated, mtdna_abslat_histogram_annotated, msat_abslat_histogram_annotated, 
#                                ncol = 3)
#abslat_all_whist_plot
#
#lat_all_plot <- grid.arrange(mtdna_pi_lat_bar_plot_annotated, mtdna_lat_bar_plot_annotated, msat_lat_bar_plot_annotated, 
#                                      ncol = 3)
#lat_all_plot
#
#lat_all_whist_plot <- grid.arrange(mtdna_pi_lat_bar_plot_annotated, mtdna_lat_bar_plot_annotated, msat_lat_bar_plot_annotated, 
#                                      mtdna_pi_lat_histogram_annotated, mtdna_lat_histogram_annotated, msat_lat_histogram_annotated, 
#                                      ncol = 3)
#lat_all_whist_plot

mtdnapi_lat_boxplot <- ggplot(data = mtdna_nonapi_small, aes(y = Pi, x = lat_round)) + 
  geom_boxplot(lwd = 1.5, color = "#09263a") + 
  geom_jitter(shape = 16, position=position_jitter(0.2), alpha = 0.45, color ="#7da0c0") #+ 
# geom_boxplot(aes(y = He, x = abslat_round), lwd = 1.5, color = "#09263a")
mtdnapi_lat_boxplot_annotated <- mtdnapi_lat_boxplot + theme_bw() + xlab("latitude") + ylab("mtdna pi") + ylim(0, 0.05) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 22, color = "black"), axis.title = element_text(size = 25), 
        legend.position = "none")
mtdnapi_lat_boxplot_annotated

########################################################################################################

#msat lon 10 degree binning

msat$lon_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)


library(DescTools)

for(i in 1:nrow(msat)) {
  cat(paste(i, " ", sep = ''))
  {msat$lon_round[i] <- DescTools::RoundTo(msat$lon[i], multiple = 10, FUN = floor)}
}

msat$lon_round <- as.factor(msat$lon_round)

#subset by lon bands and then get mean and plot that

#lon bands
#msat_lonneg180 <- subset(msat, lon_round == -180)
#msat_lonneg180_mean <- mean(msat_lonneg180$He)
#msat_lonneg180_se <- std.error(msat_lonneg180$He)
#msat_lonneg170 <- subset(msat, lon_round == -170)
#msat_lonneg170_mean <- mean(msat_lonneg170$He)
#msat_lonneg170_se <- std.error(msat_lonneg170$He)
#msat_lonneg160 <- subset(msat, lon_round == -160)
#msat_lonneg160_mean <- mean(msat_lonneg160$He)
#msat_lonneg160_se <- std.error(msat_lonneg160$He)
#msat_lonneg150 <- subset(msat, lon_round == -150)
#msat_lonneg150_mean <- mean(msat_lonneg150$He)
#msat_lonneg150_se <- std.error(msat_lonneg150$He)
#msat_lonneg140 <- subset(msat, lon_round == -140)
#msat_lonneg140_mean <- mean(msat_lonneg140$He)
#msat_lonneg140_se <- std.error(msat_lonneg140$He)
#msat_lonneg130 <- subset(msat, lon_round == -130)
#msat_lonneg130_mean <- mean(msat_lonneg130$He)
#msat_lonneg130_se <- std.error(msat_lonneg130$He)
#msat_lonneg120 <- subset(msat, lon_round == -120)
#msat_lonneg120_mean <- mean(msat_lonneg120$He)
#msat_lonneg120_se <- std.error(msat_lonneg120$He)
#msat_lonneg110 <- subset(msat, lon_round == -110)
#msat_lonneg110_mean <- mean(msat_lonneg110$He)
#msat_lonneg110_se <- std.error(msat_lonneg110$He)
#msat_lonneg100 <- subset(msat, lon_round == -100)
#msat_lonneg100_mean <- mean(msat_lonneg100$He)
#msat_lonneg100_se <- std.error(msat_lonneg100$He)
#msat_lonneg90 <- subset(msat, lon_round == -90)
#msat_lonneg90_mean <- mean(msat_lonneg90$He)
#msat_lonneg90_se <- std.error(msat_lonneg90$He)
#msat_lonneg80 <- subset(msat, lon_round == -80)
#msat_lonneg80_mean <- mean(msat_lonneg80$He)
#msat_lonneg80_se <- std.error(msat_lonneg80$He)
#msat_lonneg70 <- subset(msat, lon_round == -70)
#msat_lonneg70_mean <- mean(msat_lonneg70$He)
#msat_lonneg70_se <- std.error(msat_lonneg70$He)
#msat_lonneg60 <- subset(msat, lon_round == -60)
#msat_lonneg60_mean <- mean(msat_lonneg60$He)
#msat_lonneg60_se <- std.error(msat_lonneg60$He)
#msat_lonneg50 <- subset(msat, lon_round == -50)
#msat_lonneg50_mean <- mean(msat_lonneg50$He)
#msat_lonneg50_se <- std.error(msat_lonneg50$He)
#msat_lonneg40 <- subset(msat, lon_round == -40)
#msat_lonneg40_mean <- mean(msat_lonneg40$He)
#msat_lonneg40_se <- std.error(msat_lonneg40$He)
#msat_lonneg30 <- subset(msat, lon_round == -30)
#msat_lonneg30_mean <- mean(msat_lonneg30$He)
#msat_lonneg30_se <- std.error(msat_lonneg30$He)
#msat_lonneg20 <- subset(msat, lon_round == -20)
#msat_lonneg20_mean <- mean(msat_lonneg20$He)
#msat_lonneg20_se <- std.error(msat_lonneg20$He)
#msat_lonneg10 <- subset(msat, lon_round == -10)
#msat_lonneg10_mean <- mean(msat_lonneg10$He)
#msat_lonneg10_se <- std.error(msat_lonneg10$He)
#msat_lon0 <- subset(msat, lon_round == 0)
#msat_lon0_mean <- mean(msat_lon0$He)
#msat_lon0_se <- std.error(msat_lon0$He)
#msat_lon10 <- subset(msat, lon_round == 10)
#msat_lon10_mean <- mean(msat_lon10$He)
#msat_lon10_se <- std.error(msat_lon10$He)
#msat_lon20 <- subset(msat, lon_round == 20)
#msat_lon20_mean <- mean(msat_lon20$He)
#msat_lon20_se <- std.error(msat_lon20$He)
#msat_lon30 <- subset(msat, lon_round == 30)
#msat_lon30_mean <- mean(msat_lon30$He)
#msat_lon30_se <- std.error(msat_lon30$He)
#msat_lon40 <- subset(msat, lon_round == 40)
#msat_lon40_mean <- mean(msat_lon40$He)
#msat_lon40_se <- std.error(msat_lon40$He)
#msat_lon50 <- subset(msat, lon_round == 50)
#msat_lon50_mean <- mean(msat_lon50$He)
#msat_lon50_se <- std.error(msat_lon50$He)
#msat_lon60 <- subset(msat, lon_round == 60)
#msat_lon60_mean <- mean(msat_lon60$He)
#msat_lon60_se <- std.error(msat_lon60$He)
#msat_lon70 <- subset(msat, lon_round == 70)
#msat_lon70_mean <- mean(msat_lon70$He)
#msat_lon70_se <- std.error(msat_lon70$He)
#msat_lon80 <- subset(msat, lon_round == 80)
#msat_lon80_mean <- mean(msat_lon80$He)
#msat_lon80_se <- std.error(msat_lon80$He)
#msat_lon90 <- subset(msat, lon_round == 90)
#msat_lon90_mean <- mean(msat_lon90$He)
#msat_lon90_se <- std.error(msat_lon90$He)
#msat_lon100 <- subset(msat, lon_round == 100)
#msat_lon100_mean <- mean(msat_lon100$He)
#msat_lon100_se <- std.error(msat_lon100$He)
#msat_lon110 <- subset(msat, lon_round == 110)
#msat_lon110_mean <- mean(msat_lon110$He)
#msat_lon110_se <- std.error(msat_lon110$He)
#msat_lon120 <- subset(msat, lon_round == 120)
#msat_lon120_mean <- mean(msat_lon120$He)
#msat_lon120_se <- std.error(msat_lon120$He)
#msat_lon130 <- subset(msat, lon_round == 130)
#msat_lon130_mean <- mean(msat_lon130$He)
#msat_lon130_se <- std.error(msat_lon130$He)
#msat_lon140 <- subset(msat, lon_round == 140)
#msat_lon140_mean <- mean(msat_lon140$He)
#msat_lon140_se <- std.error(msat_lon140$He)
#msat_lon150 <- subset(msat, lon_round == 150)
#msat_lon150_mean <- mean(msat_lon150$He)
#msat_lon150_se <- std.error(msat_lon150$He)
#msat_lon160 <- subset(msat, lon_round == 160)
#msat_lon160_mean <- mean(msat_lon160$He)
#msat_lon160_se <- std.error(msat_lon160$He)
#msat_lon170 <- subset(msat, lon_round == 170)
#msat_lon170_mean <- mean(msat_lon170$He)
#msat_lon170_se <- std.error(msat_lon170$He)
#
#lon_He_mean <- c(msat_lonneg180_mean, msat_lonneg170_mean, msat_lonneg160_mean, 
#                 msat_lonneg150_mean, msat_lonneg140_mean, msat_lonneg130_mean, 
#                 msat_lonneg120_mean, msat_lonneg110_mean, msat_lonneg100_mean, 
#                 msat_lonneg90_mean, msat_lonneg80_mean, msat_lonneg70_mean, 
#                 msat_lonneg60_mean, msat_lonneg50_mean, msat_lonneg40_mean, 
#                 msat_lonneg30_mean, msat_lonneg20_mean, msat_lonneg10_mean, 
#                 msat_lon0_mean, msat_lon10_mean, msat_lon20_mean, msat_lon30_mean, 
#                 msat_lon40_mean, msat_lon50_mean, msat_lon60_mean, msat_lon70_mean, 
#                 msat_lon80_mean, msat_lon90_mean, msat_lon100_mean, msat_lon110_mean,
#                 msat_lon120_mean, msat_lon130_mean, msat_lon140_mean, msat_lon150_mean, 
#                 msat_lon160_mean, msat_lon170_mean)
#lon_He_SE <- c(msat_lonneg180_se, msat_lonneg170_se, msat_lonneg160_se, 
#                 msat_lonneg150_se, msat_lonneg140_se, msat_lonneg130_se, 
#                 msat_lonneg120_se, msat_lonneg110_se, msat_lonneg100_se, 
#                 msat_lonneg90_se, msat_lonneg80_se, msat_lonneg70_se, 
#                 msat_lonneg60_se, msat_lonneg50_se, msat_lonneg40_se, 
#                 msat_lonneg30_se, msat_lonneg20_se, msat_lonneg10_se, 
#                 msat_lon0_se, msat_lon10_se, msat_lon20_se, msat_lon30_se, 
#                 msat_lon40_se, msat_lon50_se, msat_lon60_se, msat_lon70_se, 
#                 msat_lon80_se, msat_lon90_se, msat_lon100_se, msat_lon110_se,
#                 msat_lon120_se, msat_lon130_se, msat_lon140_se, msat_lon150_se, 
#                 msat_lon160_se, msat_lon170_se)
#lon_bin <- c(-180, -170, -160, -150, -140, -130, -120, -110, -100, -
#               90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 
#             10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170)
#lon_He_bin_df <- data.frame(lon_bin, lon_He_mean, lon_He_SE)
#colnames(lon_He_bin_df) <- c("lon", "He", "SE")
#
#msat_lon_bar_plot <- ggplot(data = lon_He_bin_df, aes(x = lon, y = He, ymin = He-SE, ymax = He+SE)) + 
#  geom_bar(stat = "identity", color = "#00B8E5", fill = "#00B8E5") + xlab("longitude") + ylab("msat He") + 
#  #geom_smooth(data = msat_nonan, aes(x = lon, y = He, succ = success, fail = failure), method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x, color = "black", size = 2)
#  geom_errorbar(width = 4)
#msat_lon_bar_plot_annotated <- msat_lon_bar_plot + theme_bw() + 
#  scale_x_continuous(limits = c(-190, 180), breaks = c(-180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180)) +
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"))
#msat_lon_bar_plot_annotated
#msat_lon_histogram <- ggplot(data = msat_nonan, aes(msat_nonan$lon)) + geom_histogram(binwidth = 10) + 
#  xlab("longitude") + ylab("count")
#msat_lon_histogram_annotated <- msat_lon_histogram + theme_bw() + 
#  scale_x_continuous(limits = c(-190, 180), breaks = c(-180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180)) +
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"))
#msat_lon_histogram_annotated
#
#msat_lon_all_plot <- grid.arrange(msat_lon_bar_plot_annotated, msat_lon_histogram_annotated, ncol = 1)
#msat_lon_all_plot

msat_lon_boxplot <- ggplot(data = msat, aes(y = He, x = lon_round)) + 
 # geom_boxplot(lwd = 1.5, color = "#09263a") + 
  geom_jitter(shape = 16, position=position_jitter(0.2), alpha = 0.45, color ="#7da0c0") + 
  geom_boxplot(aes(y = He, x = lon_round), lwd = 1.5, color = "#09263a")
msat_lon_boxplot_annotated <- msat_lon_boxplot + theme_bw() + xlab("longitude") + ylab("msat He") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 22, color = "black"), axis.title = element_text(size = 25), 
        legend.position = "none")
msat_lon_boxplot_annotated

########################################################################################################

#mtdna he lon 10 degree binning

mtdna_nonahe$lon_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)

library(DescTools)

for(i in 1:nrow(mtdna_nonahe)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_nonahe$lon_round[i] <- DescTools::RoundTo(mtdna_nonahe$lon[i], multiple = 10, FUN = floor)}
}

mtdna_nonahe$lon_round <- as.factor(mtdna_nonahe$lon_round)

#subset by lon bands and then get mean and plot that

#lon bands
#mtdna_lonneg180 <- subset(mtdna_nonahe, lon_round == -180)
#mtdna_lonneg180_mean <- mean(mtdna_lonneg180$He)
#mtdna_lonneg180_se <- std.error(mtdna_lonneg180$He)
#mtdna_lonneg170 <- subset(mtdna_nonahe, lon_round == -170)
#mtdna_lonneg170_mean <- mean(mtdna_lonneg170$He)
#mtdna_lonneg170_se <- std.error(mtdna_lonneg170$He)
#mtdna_lonneg160 <- subset(mtdna_nonahe, lon_round == -160)
#mtdna_lonneg160_mean <- mean(mtdna_lonneg160$He)
#mtdna_lonneg160_se <- std.error(mtdna_lonneg160$He)
#mtdna_lonneg150 <- subset(mtdna_nonahe, lon_round == -150)
#mtdna_lonneg150_mean <- mean(mtdna_lonneg150$He)
#mtdna_lonneg150_se <- std.error(mtdna_lonneg150$He)
#mtdna_lonneg140 <- subset(mtdna_nonahe, lon_round == -140)
#mtdna_lonneg140_mean <- mean(mtdna_lonneg140$He)
#mtdna_lonneg140_se <- std.error(mtdna_lonneg140$He)
#mtdna_lonneg130 <- subset(mtdna_nonahe, lon_round == -130)
#mtdna_lonneg130_mean <- mean(mtdna_lonneg130$He)
#mtdna_lonneg130_se <- std.error(mtdna_lonneg130$He)
#mtdna_lonneg120 <- subset(mtdna_nonahe, lon_round == -120)
#mtdna_lonneg120_mean <- mean(mtdna_lonneg120$He)
#mtdna_lonneg120_se <- std.error(mtdna_lonneg120$He)
#mtdna_lonneg110 <- subset(mtdna_nonahe, lon_round == -110)
#mtdna_lonneg110_mean <- mean(mtdna_lonneg110$He)
#mtdna_lonneg110_se <- std.error(mtdna_lonneg110$He)
#mtdna_lonneg100 <- subset(mtdna_nonahe, lon_round == -100)
#mtdna_lonneg100_mean <- mean(mtdna_lonneg100$He)
#mtdna_lonneg100_se <- std.error(mtdna_lonneg100$He)
#mtdna_lonneg90 <- subset(mtdna_nonahe, lon_round == -90)
#mtdna_lonneg90_mean <- mean(mtdna_lonneg90$He)
#mtdna_lonneg90_se <- std.error(mtdna_lonneg90$He)
#mtdna_lonneg80 <- subset(mtdna_nonahe, lon_round == -80)
#mtdna_lonneg80_mean <- mean(mtdna_lonneg80$He)
#mtdna_lonneg80_se <- std.error(mtdna_lonneg80$He)
#mtdna_lonneg70 <- subset(mtdna_nonahe, lon_round == -70)
#mtdna_lonneg70_mean <- mean(mtdna_lonneg70$He)
#mtdna_lonneg70_se <- std.error(mtdna_lonneg70$He)
#mtdna_lonneg60 <- subset(mtdna_nonahe, lon_round == -60)
#mtdna_lonneg60_mean <- mean(mtdna_lonneg60$He)
#mtdna_lonneg60_se <- std.error(mtdna_lonneg60$He)
#mtdna_lonneg50 <- subset(mtdna_nonahe, lon_round == -50)
#mtdna_lonneg50_mean <- mean(mtdna_lonneg50$He)
#mtdna_lonneg50_se <- std.error(mtdna_lonneg50$He)
#mtdna_lonneg40 <- subset(mtdna_nonahe, lon_round == -40)
#mtdna_lonneg40_mean <- mean(mtdna_lonneg40$He)
#mtdna_lonneg40_se <- std.error(mtdna_lonneg40$He)
#mtdna_lonneg30 <- subset(mtdna_nonahe, lon_round == -30)
#mtdna_lonneg30_mean <- mean(mtdna_lonneg30$He)
#mtdna_lonneg30_se <- std.error(mtdna_lonneg30$He)
#mtdna_lonneg20 <- subset(mtdna_nonahe, lon_round == -20)
#mtdna_lonneg20_mean <- mean(mtdna_lonneg20$He)
#mtdna_lonneg20_se <- std.error(mtdna_lonneg20$He)
#mtdna_lonneg10 <- subset(mtdna_nonahe, lon_round == -10)
#mtdna_lonneg10_mean <- mean(mtdna_lonneg10$He)
#mtdna_lonneg10_se <- std.error(mtdna_lonneg10$He)
#mtdna_lon0 <- subset(mtdna_nonahe, lon_round == 0)
#mtdna_lon0_mean <- mean(mtdna_lon0$He)
#mtdna_lon0_se <- std.error(mtdna_lon0$He)
#mtdna_lon10 <- subset(mtdna_nonahe, lon_round == 10)
#mtdna_lon10_mean <- mean(mtdna_lon10$He)
#mtdna_lon10_se <- std.error(mtdna_lon10$He)
#mtdna_lon20 <- subset(mtdna_nonahe, lon_round == 20)
#mtdna_lon20_mean <- mean(mtdna_lon20$He)
#mtdna_lon20_se <- std.error(mtdna_lon20$He)
#mtdna_lon30 <- subset(mtdna_nonahe, lon_round == 30)
#mtdna_lon30_mean <- mean(mtdna_lon30$He)
#mtdna_lon30_se <- std.error(mtdna_lon30$He)
#mtdna_lon40 <- subset(mtdna_nonahe, lon_round == 40)
#mtdna_lon40_mean <- mean(mtdna_lon40$He)
#mtdna_lon40_se <- std.error(mtdna_lon40$He)
#mtdna_lon50 <- subset(mtdna_nonahe, lon_round == 50)
#mtdna_lon50_mean <- mean(mtdna_lon50$He)
#mtdna_lon50_se <- std.error(mtdna_lon50$He)
#mtdna_lon60 <- subset(mtdna_nonahe, lon_round == 60)
#mtdna_lon60_mean <- mean(mtdna_lon60$He)
#mtdna_lon60_se <- std.error(mtdna_lon60$He)
#mtdna_lon70 <- subset(mtdna_nonahe, lon_round == 70)
#mtdna_lon70_mean <- mean(mtdna_lon70$He)
#mtdna_lon70_se <- std.error(mtdna_lon70$He)
#mtdna_lon80 <- subset(mtdna_nonahe, lon_round == 80)
#mtdna_lon80_mean <- mean(mtdna_lon80$He)
#mtdna_lon80_se <- std.error(mtdna_lon80$He)
#mtdna_lon90 <- subset(mtdna_nonahe, lon_round == 90)
#mtdna_lon90_mean <- mean(mtdna_lon90$He)
#mtdna_lon90_se <- std.error(mtdna_lon90$He)
#mtdna_lon100 <- subset(mtdna_nonahe, lon_round == 100)
#mtdna_lon100_mean <- mean(mtdna_lon100$He)
#mtdna_lon100_se <- std.error(mtdna_lon100$He)
#mtdna_lon110 <- subset(mtdna_nonahe, lon_round == 110)
#mtdna_lon110_mean <- mean(mtdna_lon110$He)
#mtdna_lon110_se <- std.error(mtdna_lon110$He)
#mtdna_lon120 <- subset(mtdna_nonahe, lon_round == 120)
#mtdna_lon120_mean <- mean(mtdna_lon120$He)
#mtdna_lon120_se <- std.error(mtdna_lon120$He)
#mtdna_lon130 <- subset(mtdna_nonahe, lon_round == 130)
#mtdna_lon130_mean <- mean(mtdna_lon130$He)
#mtdna_lon130_se <- std.error(mtdna_lon130$He)
#mtdna_lon140 <- subset(mtdna_nonahe, lon_round == 140)
#mtdna_lon140_mean <- mean(mtdna_lon140$He)
#mtdna_lon140_se <- std.error(mtdna_lon140$He)
#mtdna_lon150 <- subset(mtdna_nonahe, lon_round == 150)
#mtdna_lon150_mean <- mean(mtdna_lon150$He)
#mtdna_lon150_se <- std.error(mtdna_lon150$He)
#mtdna_lon160 <- subset(mtdna_nonahe, lon_round == 160)
#mtdna_lon160_mean <- mean(mtdna_lon160$He)
#mtdna_lon160_se <- std.error(mtdna_lon160$He)
#mtdna_lon170 <- subset(mtdna_nonahe, lon_round == 170)
#mtdna_lon170_mean <- mean(mtdna_lon170$He)
#mtdna_lon170_se <- std.error(mtdna_lon170$He)
#
#lon_He_mean_mtdna <- c(mtdna_lonneg180_mean, mtdna_lonneg170_mean, mtdna_lonneg160_mean, 
#                 mtdna_lonneg150_mean, mtdna_lonneg140_mean, mtdna_lonneg130_mean, 
#                 mtdna_lonneg120_mean, mtdna_lonneg110_mean, mtdna_lonneg100_mean, 
#                 mtdna_lonneg90_mean, mtdna_lonneg80_mean, mtdna_lonneg70_mean, 
#                 mtdna_lonneg60_mean, mtdna_lonneg50_mean, mtdna_lonneg40_mean, 
#                 mtdna_lonneg30_mean, mtdna_lonneg20_mean, mtdna_lonneg10_mean, 
#                 mtdna_lon0_mean, mtdna_lon10_mean, mtdna_lon20_mean, mtdna_lon30_mean, 
#                 mtdna_lon40_mean, mtdna_lon50_mean, mtdna_lon60_mean, mtdna_lon70_mean, 
#                 mtdna_lon80_mean, mtdna_lon90_mean, mtdna_lon100_mean, mtdna_lon110_mean,
#                 mtdna_lon120_mean, mtdna_lon130_mean, mtdna_lon140_mean, mtdna_lon150_mean, 
#                 mtdna_lon160_mean, mtdna_lon170_mean)
#lon_He_SE_mtdna <- c(mtdna_lonneg180_se, mtdna_lonneg170_se, mtdna_lonneg160_se, 
#                       mtdna_lonneg150_se, mtdna_lonneg140_se, mtdna_lonneg130_se, 
#                       mtdna_lonneg120_se, mtdna_lonneg110_se, mtdna_lonneg100_se, 
#                       mtdna_lonneg90_se, mtdna_lonneg80_se, mtdna_lonneg70_se, 
#                       mtdna_lonneg60_se, mtdna_lonneg50_se, mtdna_lonneg40_se, 
#                       mtdna_lonneg30_se, mtdna_lonneg20_se, mtdna_lonneg10_se, 
#                       mtdna_lon0_se, mtdna_lon10_se, mtdna_lon20_se, mtdna_lon30_se, 
#                       mtdna_lon40_se, mtdna_lon50_se, mtdna_lon60_se, mtdna_lon70_se, 
#                       mtdna_lon80_se, mtdna_lon90_se, mtdna_lon100_se, mtdna_lon110_se,
#                       mtdna_lon120_se, mtdna_lon130_se, mtdna_lon140_se, mtdna_lon150_se, 
#                       mtdna_lon160_se, mtdna_lon170_se)
#lon_bin <- c(-180, -170, -160, -150, -140, -130, -120, -110, -100, -
#               90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 
#             10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170)
#lon_He_bin_df_mtdna <- data.frame(lon_bin, lon_He_mean_mtdna, lon_He_SE_mtdna)
#colnames(lon_He_bin_df_mtdna) <- c("lon", "He", "SE")
#
#mtdna_lon_bar_plot <- ggplot(data = lon_He_bin_df_mtdna, aes(x = lon, y = He, ymin = He-SE, ymax = He+SE)) + 
#  geom_bar(stat = "identity", color = "#00B8E5", fill = "#00B8E5") + xlab("longitude") + ylab("mtdna Hd") + 
#  #geom_smooth(data = mtdna_nonahe, aes(x = lon, y = He, succ = success, fail = failure), method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x, color = "black", size = 2)
#  geom_errorbar(width = 4)
#mtdna_lon_bar_plot_annotated <- mtdna_lon_bar_plot + theme_bw() + 
#  scale_x_continuous(limits = c(-190, 180), breaks = c(-180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180)) +
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"))
#mtdna_lon_bar_plot_annotated
#mtdna_lon_histogram <- ggplot(data = mtdna_nonahe, aes(mtdna_nonahe$lon)) + geom_histogram(binwidth = 10) + 
#  xlab("longitude") + ylab("count")
#mtdna_lon_histogram_annotated <- mtdna_lon_histogram + theme_bw() + 
#  scale_x_continuous(limits = c(-190, 180), breaks = c(-180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180)) +
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"))
#mtdna_lon_histogram_annotated
#
#mtdna_lon_all_plot <- grid.arrange(mtdna_lon_bar_plot_annotated, mtdna_lon_histogram_annotated, ncol = 1)
#mtdna_lon_all_plot

mtdnahd_lon_boxplot <- ggplot(data = mtdna_nonahe, aes(y = He, x = lon_round)) + 
  geom_boxplot(lwd = 1.5, color = "#09263a") + 
  geom_jitter(shape = 16, position=position_jitter(0.2), alpha = 0.45, color ="#7da0c0") #+ 
 # geom_boxplot(aes(y = He, x = lon_round), lwd = 1.5, color = "#09263a")
mtdnahd_lon_boxplot_annotated <- mtdnahd_lon_boxplot + theme_bw() + xlab("longitude") + ylab("mtdna Hd") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 22, color = "black"), axis.title = element_text(size = 25), 
        legend.position = "none")
mtdnahd_lon_boxplot_annotated

##################################################################################################################

#mtdna pi lon 10 degree binning

mtdna_nonapi_small$lon_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)

library(DescTools)

for(i in 1:nrow(mtdna_nonapi_small)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_nonapi_small$lon_round[i] <- DescTools::RoundTo(mtdna_nonapi_small$lon[i], multiple = 10, FUN = floor)}
}

mtdna_nonapi_small$lon_round <- as.factor(mtdna_nonapi_small$lon_round)


#subset by lon bands and then get mean and plot that

#lon bands
#mtdna_pi_lonneg180 <- subset(mtdna_nonapi_small, lon_round == -180)
#mtdna_pi_lonneg180_mean <- mean(mtdna_pi_lonneg180$Pi)
#mtdna_pi_lonneg180_se <- std.error(mtdna_pi_lonneg180$Pi)
#mtdna_pi_lonneg170 <- subset(mtdna_nonapi_small, lon_round == -170)
#mtdna_pi_lonneg170_mean <- mean(mtdna_pi_lonneg170$Pi)
#mtdna_pi_lonneg170_se <- std.error(mtdna_pi_lonneg170$Pi)
#mtdna_pi_lonneg160 <- subset(mtdna_nonapi_small, lon_round == -160)
#mtdna_pi_lonneg160_mean <- mean(mtdna_pi_lonneg160$Pi)
#mtdna_pi_lonneg160_se <- std.error(mtdna_pi_lonneg160$Pi)
#mtdna_pi_lonneg150 <- subset(mtdna_nonapi_small, lon_round == -150)
#mtdna_pi_lonneg150_mean <- mean(mtdna_pi_lonneg150$Pi)
#mtdna_pi_lonneg150_se <- std.error(mtdna_pi_lonneg150$Pi)
#mtdna_pi_lonneg140 <- subset(mtdna_nonapi_small, lon_round == -140)
#mtdna_pi_lonneg140_mean <- mean(mtdna_pi_lonneg140$Pi)
#mtdna_pi_lonneg140_se <- std.error(mtdna_pi_lonneg140$Pi)
#mtdna_pi_lonneg130 <- subset(mtdna_nonapi_small, lon_round == -130)
#mtdna_pi_lonneg130_mean <- mean(mtdna_pi_lonneg130$Pi)
#mtdna_pi_lonneg130_se <- std.error(mtdna_pi_lonneg130$Pi)
#mtdna_pi_lonneg120 <- subset(mtdna_nonapi_small, lon_round == -120)
#mtdna_pi_lonneg120_mean <- mean(mtdna_pi_lonneg120$Pi)
#mtdna_pi_lonneg120_se <- std.error(mtdna_pi_lonneg120$Pi)
#mtdna_pi_lonneg110 <- subset(mtdna_nonapi_small, lon_round == -110)
#mtdna_pi_lonneg110_mean <- mean(mtdna_pi_lonneg110$Pi)
#mtdna_pi_lonneg110_se <- std.error(mtdna_pi_lonneg110$Pi)
#mtdna_pi_lonneg100 <- subset(mtdna_nonapi_small, lon_round == -100)
#mtdna_pi_lonneg100_mean <- mean(mtdna_pi_lonneg100$Pi)
#mtdna_pi_lonneg100_se <- std.error(mtdna_pi_lonneg100$Pi)
#mtdna_pi_lonneg90 <- subset(mtdna_nonapi_small, lon_round == -90)
#mtdna_pi_lonneg90_mean <- mean(mtdna_pi_lonneg90$Pi)
#mtdna_pi_lonneg90_se <- std.error(mtdna_pi_lonneg90$Pi)
#mtdna_pi_lonneg80 <- subset(mtdna_nonapi_small, lon_round == -80)
#mtdna_pi_lonneg80_mean <- mean(mtdna_pi_lonneg80$Pi)
#mtdna_pi_lonneg80_se <- std.error(mtdna_pi_lonneg80$Pi)
#mtdna_pi_lonneg70 <- subset(mtdna_nonapi_small, lon_round == -70)
#mtdna_pi_lonneg70_mean <- mean(mtdna_pi_lonneg70$Pi)
#mtdna_pi_lonneg70_se <- std.error(mtdna_pi_lonneg70$Pi)
#mtdna_pi_lonneg60 <- subset(mtdna_nonapi_small, lon_round == -60)
#mtdna_pi_lonneg60_mean <- mean(mtdna_pi_lonneg60$Pi)
#mtdna_pi_lonneg60_se <- std.error(mtdna_pi_lonneg60$Pi)
#mtdna_pi_lonneg50 <- subset(mtdna_nonapi_small, lon_round == -50)
#mtdna_pi_lonneg50_mean <- mean(mtdna_pi_lonneg50$Pi)
#mtdna_pi_lonneg50_se <- std.error(mtdna_pi_lonneg50$Pi)
#mtdna_pi_lonneg40 <- subset(mtdna_nonapi_small, lon_round == -40)
#mtdna_pi_lonneg40_mean <- mean(mtdna_pi_lonneg40$Pi)
#mtdna_pi_lonneg40_se <- std.error(mtdna_pi_lonneg40$Pi)
#mtdna_pi_lonneg30 <- subset(mtdna_nonapi_small, lon_round == -30)
#mtdna_pi_lonneg30_mean <- mean(mtdna_pi_lonneg30$Pi)
#mtdna_pi_lonneg30_se <- std.error(mtdna_pi_lonneg30$Pi)
#mtdna_pi_lonneg20 <- subset(mtdna_nonapi_small, lon_round == -20)
#mtdna_pi_lonneg20_mean <- mean(mtdna_pi_lonneg20$Pi)
#mtdna_pi_lonneg20_se <- std.error(mtdna_pi_lonneg20$Pi)
#mtdna_pi_lonneg10 <- subset(mtdna_nonapi_small, lon_round == -10)
#mtdna_pi_lonneg10_mean <- mean(mtdna_pi_lonneg10$Pi)
#mtdna_pi_lonneg10_se <- std.error(mtdna_pi_lonneg10$Pi)
#mtdna_pi_lon0 <- subset(mtdna_nonapi_small, lon_round == 0)
#mtdna_pi_lon0_mean <- mean(mtdna_pi_lon0$Pi)
#mtdna_pi_lon0_se <- std.error(mtdna_pi_lon0$Pi)
#mtdna_pi_lon10 <- subset(mtdna_nonapi_small, lon_round == 10)
#mtdna_pi_lon10_mean <- mean(mtdna_pi_lon10$Pi)
#mtdna_pi_lon10_se <- std.error(mtdna_pi_lon10$Pi)
#mtdna_pi_lon20 <- subset(mtdna_nonapi_small, lon_round == 20)
#mtdna_pi_lon20_mean <- mean(mtdna_pi_lon20$Pi)
#mtdna_pi_lon20_se <- std.error(mtdna_pi_lon20$Pi)
#mtdna_pi_lon30 <- subset(mtdna_nonapi_small, lon_round == 30)
#mtdna_pi_lon30_mean <- mean(mtdna_pi_lon30$Pi)
#mtdna_pi_lon30_se <- std.error(mtdna_pi_lon30$Pi)
#mtdna_pi_lon40 <- subset(mtdna_nonapi_small, lon_round == 40)
#mtdna_pi_lon40_mean <- mean(mtdna_pi_lon40$Pi)
#mtdna_pi_lon40_se <- std.error(mtdna_pi_lon40$Pi)
#mtdna_pi_lon50 <- subset(mtdna_nonapi_small, lon_round == 50)
#mtdna_pi_lon50_mean <- mean(mtdna_pi_lon50$Pi)
#mtdna_pi_lon50_se <- std.error(mtdna_pi_lon50$Pi)
#mtdna_pi_lon60 <- subset(mtdna_nonapi_small, lon_round == 60)
#mtdna_pi_lon60_mean <- mean(mtdna_pi_lon60$Pi)
#mtdna_pi_lon60_se <- std.error(mtdna_pi_lon60$Pi)
#mtdna_pi_lon70 <- subset(mtdna_nonapi_small, lon_round == 70)
#mtdna_pi_lon70_mean <- mean(mtdna_pi_lon70$Pi)
#mtdna_pi_lon70_se <- std.error(mtdna_pi_lon70$Pi)
#mtdna_pi_lon80 <- subset(mtdna_nonapi_small, lon_round == 80)
#mtdna_pi_lon80_mean <- mean(mtdna_pi_lon80$Pi)
#mtdna_pi_lon80_se <- std.error(mtdna_pi_lon80$Pi)
#mtdna_pi_lon90 <- subset(mtdna_nonapi_small, lon_round == 90)
#mtdna_pi_lon90_mean <- mean(mtdna_pi_lon90$Pi)
#mtdna_pi_lon90_se <- std.error(mtdna_pi_lon90$Pi)
#mtdna_pi_lon100 <- subset(mtdna_nonapi_small, lon_round == 100)
#mtdna_pi_lon100_mean <- mean(mtdna_pi_lon100$Pi)
#mtdna_pi_lon100_se <- std.error(mtdna_pi_lon100$Pi)
#mtdna_pi_lon110 <- subset(mtdna_nonapi_small, lon_round == 110)
#mtdna_pi_lon110_mean <- mean(mtdna_pi_lon110$Pi)
#mtdna_pi_lon110_se <- std.error(mtdna_pi_lon110$Pi)
#mtdna_pi_lon120 <- subset(mtdna_nonapi_small, lon_round == 120)
#mtdna_pi_lon120_mean <- mean(mtdna_pi_lon120$Pi)
#mtdna_pi_lon120_se <- std.error(mtdna_pi_lon120$Pi)
#mtdna_pi_lon130 <- subset(mtdna_nonapi_small, lon_round == 130)
#mtdna_pi_lon130_mean <- mean(mtdna_pi_lon130$Pi)
#mtdna_pi_lon130_se <- std.error(mtdna_pi_lon130$Pi)
#mtdna_pi_lon140 <- subset(mtdna_nonapi_small, lon_round == 140)
#mtdna_pi_lon140_mean <- mean(mtdna_pi_lon140$Pi)
#mtdna_pi_lon140_se <- std.error(mtdna_pi_lon140$Pi)
#mtdna_pi_lon150 <- subset(mtdna_nonapi_small, lon_round == 150)
#mtdna_pi_lon150_mean <- mean(mtdna_pi_lon150$Pi)
#mtdna_pi_lon150_se <- std.error(mtdna_pi_lon150$Pi)
#mtdna_pi_lon160 <- subset(mtdna_nonapi_small, lon_round == 160)
#mtdna_pi_lon160_mean <- mean(mtdna_pi_lon160$Pi)
#mtdna_pi_lon160_se <- std.error(mtdna_pi_lon160$Pi)
#mtdna_pi_lon170 <- subset(mtdna_nonapi_small, lon_round == 170)
#mtdna_pi_lon170_mean <- mean(mtdna_pi_lon170$Pi)
#mtdna_pi_lon170_se <- std.error(mtdna_pi_lon170$Pi)
#
#lon_Pi_mean_mtdna <- c(mtdna_pi_lonneg180_mean, mtdna_pi_lonneg170_mean, mtdna_pi_lonneg160_mean, 
#                       mtdna_pi_lonneg150_mean, mtdna_pi_lonneg140_mean, mtdna_pi_lonneg130_mean, 
#                       mtdna_pi_lonneg120_mean, mtdna_pi_lonneg110_mean, mtdna_pi_lonneg100_mean, 
#                       mtdna_pi_lonneg90_mean, mtdna_pi_lonneg80_mean, mtdna_pi_lonneg70_mean, 
#                       mtdna_pi_lonneg60_mean, mtdna_pi_lonneg50_mean, mtdna_pi_lonneg40_mean, 
#                       mtdna_pi_lonneg30_mean, mtdna_pi_lonneg20_mean, mtdna_pi_lonneg10_mean, 
#                       mtdna_pi_lon0_mean, mtdna_pi_lon10_mean, mtdna_pi_lon20_mean, mtdna_pi_lon30_mean, 
#                       mtdna_pi_lon40_mean, mtdna_pi_lon50_mean, mtdna_pi_lon60_mean, mtdna_pi_lon70_mean, 
#                       mtdna_pi_lon80_mean, mtdna_pi_lon90_mean, mtdna_pi_lon100_mean, mtdna_pi_lon110_mean,
#                       mtdna_pi_lon120_mean, mtdna_pi_lon130_mean, mtdna_pi_lon140_mean, mtdna_pi_lon150_mean, 
#                       mtdna_pi_lon160_mean, mtdna_pi_lon170_mean)
#lon_Pi_SE_mtdna <- c(mtdna_pi_lonneg180_se, mtdna_pi_lonneg170_se, mtdna_pi_lonneg160_se, 
#                       mtdna_pi_lonneg150_se, mtdna_pi_lonneg140_se, mtdna_pi_lonneg130_se, 
#                       mtdna_pi_lonneg120_se, mtdna_pi_lonneg110_se, mtdna_pi_lonneg100_se, 
#                       mtdna_pi_lonneg90_se, mtdna_pi_lonneg80_se, mtdna_pi_lonneg70_se, 
#                       mtdna_pi_lonneg60_se, mtdna_pi_lonneg50_se, mtdna_pi_lonneg40_se, 
#                       mtdna_pi_lonneg30_se, mtdna_pi_lonneg20_se, mtdna_pi_lonneg10_se, 
#                       mtdna_pi_lon0_se, mtdna_pi_lon10_se, mtdna_pi_lon20_se, mtdna_pi_lon30_se, 
#                       mtdna_pi_lon40_se, mtdna_pi_lon50_se, mtdna_pi_lon60_se, mtdna_pi_lon70_se, 
#                       mtdna_pi_lon80_se, mtdna_pi_lon90_se, mtdna_pi_lon100_se, mtdna_pi_lon110_se,
#                       mtdna_pi_lon120_se, mtdna_pi_lon130_se, mtdna_pi_lon140_se, mtdna_pi_lon150_se, 
#                       mtdna_pi_lon160_se, mtdna_pi_lon170_se)
#lon_bin <- c(-180, -170, -160, -150, -140, -130, -120, -110, -100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 
#             10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170)
#lon_Pi_bin_df_mtdna <- data.frame(lon_bin, lon_Pi_mean_mtdna, lon_Pi_SE_mtdna)
#colnames(lon_Pi_bin_df_mtdna) <- c("lon", "Pi", "SE")
#
#mtdna_pi_lon_bar_plot <- ggplot(data = lon_Pi_bin_df_mtdna, aes(x = lon, y = Pi, ymin = Pi-SE, ymax = Pi+SE)) + 
#  geom_bar(stat = "identity", color = "#00B8E5", fill = "#00B8E5") + xlab("longitude") + ylab("mtdna Pi") + 
#  #geom_smooth(data = mtdna_nonapi_small, aes(x = lon, y = Pi), method = "lm", formula = y ~ x, color = "black", size = 2)
#  geom_errorbar(width = 4)
#mtdna_pi_lon_bar_plot_annotated <- mtdna_pi_lon_bar_plot + theme_bw() + 
#  scale_x_continuous(limits = c(-190, 180), breaks = c(-180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180)) +
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"))
#mtdna_pi_lon_bar_plot_annotated
#mtdna_pi_lon_histogram <- ggplot(data = mtdna_nonapi_small, aes(mtdna_nonapi_small$lon)) + geom_histogram(binwidth = 10) + 
#  xlab("longitude") + ylab("count")
#mtdna_pi_lon_histogram_annotated <- mtdna_pi_lon_histogram + theme_bw() + 
#  scale_x_continuous(limits = c(-190, 180), breaks = c(-180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180)) +
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"))
#mtdna_pi_lon_histogram_annotated
#
#mtdna_pi_lon_all_plot <- grid.arrange(mtdna_pi_lon_bar_plot_annotated, mtdna_pi_lon_histogram_annotated, ncol = 1)
#mtdna_pi_lon_all_plot
#
#lon_all_whist_plot <- grid.arrange(mtdna_pi_lon_bar_plot_annotated, mtdna_lon_bar_plot_annotated, msat_lon_bar_plot_annotated, 
#                                mtdna_pi_lon_histogram_annotated, mtdna_lon_histogram_annotated, msat_lon_histogram_annotated, ncol = 3)
#lon_all_whist_plot
#
#lon_all_plot <- grid.arrange(mtdna_pi_lon_bar_plot_annotated, mtdna_lon_bar_plot_annotated, msat_lon_bar_plot_annotated, 
#                             ncol = 3)
#lon_all_plot

mtdnapi_lon_boxplot <- ggplot(data = mtdna_nonapi_small, aes(y = Pi, x = lon_round)) + 
  geom_boxplot(lwd = 1.5, color = "#09263a") + 
  geom_jitter(shape = 16, position=position_jitter(0.2), alpha = 0.45, color ="#7da0c0") #+ 
# geom_boxplot(aes(y = He, x = lon_round), lwd = 1.5, color = "#09263a")
mtdnapi_lon_boxplot_annotated <- mtdnapi_lon_boxplot + theme_bw() + xlab("longitude") + ylab("mtdna pi") + ylim(0, 0.05) +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 22, color = "black"), axis.title = element_text(size = 25), 
        legend.position = "none")
mtdnapi_lon_boxplot_annotated
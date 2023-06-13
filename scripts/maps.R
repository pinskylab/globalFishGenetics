################################################### Script to plot sampling locations ########################################################

#Plot either each sampled population as a point OR combine as hexbins (every 3 degrees lat/lon)
#Maps of global chlorophyll & temp as well

#Plot size: 2000 x 1500

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(maps) #v.3.4.1
library(mapdata) #v.2.3.1
library(tidyverse) #v.2.0.0
library(ggthemes) #v.4.2.4
library(scales) #v.1.2.1

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv")
msat <- read.csv("output/msatloci_assembled.csv")
msat_env <- read.csv("output/msat_env.csv")

#merge dataframes
msat <- cbind(msat[, -1], msat_env[, c('sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin',
                                       'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange',
                                       'chloroA.BO_chlomax', 'chloroA.BO_chlomin')]) #merge not working for some reason, cbind bc in same order

#pull world map data
world_map <- map_data("world") %>% 
  filter(! long > 180)

################################################################################################################################################

######## Calculate tropical coverage ########

mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)

msat_tropics <- msat %>%
  filter(lat >= -23.4 & lat <= 23.4) %>%
  distinct(lat, lon, .keep_all = TRUE)

mtdna_tropics <- mtdna_small %>%
  filter(lat >= -23.4 & lat <= 23.4) %>%
  distinct(lat, lon, .keep_all = TRUE)

################################################################################################################################################

######## Create msat maps ########

#points only map
msat_he_point_plot <- world_map %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, color = "lightgray", fill = "lightgray") +
  expand_limits(x = world_map$long, y = world_map$lat) +  
  geom_point(data = msat, aes(x = lon, y = lat, color = He), 
             alpha = 0.5, size = 8, inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_hd, mapping = aes(x = lon),
           color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_hd, mapping = aes(y = lat), 
           color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  scale_color_gradient(low = "#1f2c32", high = "#e3ebed") +
  annotate("text", x = -170, y = 85, label = "C", size = 30) +
  ylim(c(-90, 90)) + labs(color = "He") + 
  xlab("Longitude") + ylab("Latitude") + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 50, color = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 50, color = "black", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "right", 
        legend.text = element_text(size = 50), 
        legend.title = element_text(size = 50), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(color = guide_colourbar(barwidth = unit(3, "cm"), barheight = unit(10, "cm")))
msat_he_point_plot

#hex bin mean map
msat_he_hexbin_plot <- world_map %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, color = "lightgray", fill = "lightgray") +
  expand_limits(x = world_map$long, y = world_map$lat) +  
  stat_summary_hex(data = msat, mapping= aes(x = lon, y = lat, z = He), 
                   fun = function(x) mean(x), binwidth = c(3, 3), #mean He (z) every 3 degrees lat and long bins
                   inherit.aes = FALSE) + 
  geom_rug(data = msat, mapping = aes(x = lon), 
           color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  geom_rug(data = msat, mapping = aes(y = lat), 
           color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  scale_fill_gradient(low = "#1f2c32", high = "#e3ebed") +
  annotate("text", x = -170, y = 85, label = "C", size = 30) +
  ylim(c(-90, 90)) + labs(fill = "He") + 
  xlab("Longitude") + ylab("Latitude") + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 50, color = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 50, color = "black", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "right", 
        legend.text = element_text(size = 50), 
        legend.title = element_text(size = 50), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(3, "cm"), barheight = unit(10, "cm")))
msat_he_hexbin_plot

###################################################################################################################################

######## Create mtdna maps ########

#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)
  mtdna_small_hd <- subset(mtdna_small, mtdna_small$He != "NA")
  mtdna_small_pi <- subset(mtdna_small, mtdna_small$Pi != "NA")

#### mtdna hd maps ####  

#points only map
mtdna_hd_point_plot <- world_map %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, color = "lightgray", fill = "lightgray") +
  expand_limits(x = world_map$long, y = world_map$lat) +  
  geom_point(data = mtdna_small_hd, aes(x = lon, y = lat, color = He), 
             alpha = 0.5, size = 8, inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_hd, mapping = aes(x = lon), 
           color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_hd, mapping = aes(y = lat), 
           color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  scale_color_gradient(low = "#1f2c32", high = "#e3ebed") +
  annotate("text", x = -170, y = 85, label = "B", size = 30) +
  ylim(c(-90, 90)) + labs(color = "Hd") + 
  xlab("Longitude") + ylab("Latitude") + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 50, color = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 50, color = "black", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "right", 
        legend.text = element_text(size = 50), 
        legend.title = element_text(size = 50), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(color = guide_colourbar(barwidth = unit(3, "cm"), barheight = unit(10, "cm")))
mtdna_hd_point_plot
  
#hex bin map
mtdna_hd_hexbin_plot <- world_map %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, color = "lightgray", fill = "lightgray") +
  expand_limits(x = world_map$long, y = world_map$lat) +  
  stat_summary_hex(data = mtdna_small_hd, mapping= aes(x = lon, y = lat, z = He), 
                   fun = function(x) mean(x), binwidth = c(3, 3), #mean Hd (z) every 3 degrees lat and long bins
                   inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_hd, mapping = aes(x = lon), color = "#282828", 
           alpha = 0.5, inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_hd, mapping = aes(y = lat), color = "#282828", 
           alpha = 0.5, inherit.aes = FALSE) + 
  scale_fill_gradient(low = "#1f2c32", high = "#e3ebed") + 
  annotate("text", x = -170, y = 85, label = "B", size = 30) +
  ylim(c(-90, 90)) + labs(fill = "Hd") + 
  xlab("Longitude") + ylab("Latitude") + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 50, color = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 50, color = "black", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "right", 
        legend.text = element_text(size = 50), 
        legend.title = element_text(size = 50), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(3, "cm"), barheight = unit(10, "cm")))
mtdna_hd_hexbin_plot

#### mtdna pi maps ####

#points only map
mtdna_pi_point_plot <- world_map %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, color = "lightgray", fill = "lightgray") +
  expand_limits(x = world_map$long, y = world_map$lat) +  
  geom_point(data = mtdna_small_pi, aes(x = lon, y = lat, color = Pi), 
             alpha = 0.5, size = 8, inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_hd, mapping = aes(x = lon), 
           color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_hd, mapping = aes(y = lat), 
           color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  scale_color_gradient(trans = "log10", low = "#1f2c32", high = "#e3ebed", 
                       labels = scales::label_log()) + 
  annotate("text", x = -170, y = 85, label = "A", size = 30) +
  ylim(c(-90, 90)) + labs(color = "π") + 
  xlab("Longitude") + ylab("Latitude") + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 50, color = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 50, color = "black", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "right", 
        legend.text = element_text(size = 50), 
        legend.title = element_text(size = 50), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(color = guide_colourbar(barwidth = unit(3, "cm"), barheight = unit(10, "cm")))
mtdna_pi_point_plot

#hex bin map
mtdna_pi_hexbin_plot <- world_map %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, color = "lightgray", fill = "lightgray") +
  expand_limits(x = world_map$long, y = world_map$lat) +  
  stat_summary_hex(data = mtdna_small_pi, mapping= aes(x = lon, y = lat, z = Pi), 
                   fun = function(x) mean(x), binwidth = c(3, 3), #mean Pi (z) every 3 degrees lat and long bins
                   inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_pi, mapping = aes(x = lon), 
           color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_pi, mapping = aes(y = lat), 
           color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  scale_fill_gradient(trans = "log10", low = "#1f2c32", high = "#e3ebed", 
                      labels = scales::label_log()) +
  annotate("text", x = -170, y = 85, label = "A", size = 30) +
  ylim(c(-90, 90)) + labs(fill = "π") + 
  xlab("Longitude") + ylab("Latitude") + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 50, color = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 50, color = "black", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "right", 
        legend.text = element_text(size = 50), 
        legend.title = element_text(size = 50), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(3, "cm"), barheight = unit(10, "cm")))
mtdna_pi_hexbin_plot

######################################################################################################################################

######## Create environmental data maps ########

#chlorophyll mean map
chlomean_hexbin_plot <- world_map %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, color = "lightgray", fill = "lightgray") +
  expand_limits(x = world_map$long, y = world_map$lat) +  
  stat_summary_hex(data = msat, mapping= aes(x = lon, y = lat, z = chloroA.BO_chlomean), 
                   fun = function(x) mean(x), binwidth = c(3, 3),
                   inherit.aes = FALSE) + 
  scale_fill_gradient(trans = "log10", low = "#c6cbc8", high = "#212924", 
                      labels = scales::label_log()) +
  annotate("text", x = -175, y = 85, label = "A", size = 30) +
  ylim(c(-90, 90)) + labs(fill = "Chlorophyll") + 
  xlab("Longitude") + ylab("Latitude") + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 50, color = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 50, color = "black", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "right", 
        legend.text = element_text(size = 30), 
        legend.title = element_text(size = 30), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(3, "cm"), barheight = unit(10, "cm")))
chlomean_hexbin_plot

#SST mean map
sstmean_hexbin_plot <- world_map %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, color = "lightgray", fill = "lightgray") +
  expand_limits(x = world_map$long, y = world_map$lat) +  
  stat_summary_hex(data = msat, mapping= aes(x = lon, y = lat, z = sst.BO_sstmean), 
                   fun = function(x) mean(x), binwidth = c(3, 3),
                   inherit.aes = FALSE) + 
  scale_fill_gradient(low = "#ebcb99", high = "#905800") +
  annotate("text", x = -175, y = 85, label = "B", size = 30) +
  ylim(c(-90, 90)) + labs(fill = "SST") + 
  xlab("Longitude") + ylab("Latitude") + 
  theme(panel.background = element_blank(), 
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 50, color = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 50, color = "black", margin = margin(r = 10)),
        axis.line = element_line(linewidth = 4, color = "black"),
        legend.position = "right", 
        legend.text = element_text(size = 30), 
        legend.title = element_text(size = 30), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(3, "cm"), barheight = unit(10, "cm")))
sstmean_hexbin_plot

#################################### Plot sampling locations ###########################################

#set-up

remove(list = ls())

#load libraries
library(maps)
library(mapdata)
library(tidyverse)

msat <- read.csv("output/latlon_msat_2020-12-16.csv")
mtdna <- read.csv("output/latlon_mtdna_2020-10-20.csv")
merged <- rbind(msat, mtdna)
unique_merged <- unique(merged)

#read in data
msat <- read.csv("output/msat_assembled_2.csv")
mtdna <- read.csv("output/mtdna_assembled_2.csv")

#pull world map data
geogr_data <- map_data('world')

################################################################################################################################################

######## Create msat map ########

msat_plot <- ggplot(geogr_data, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "lightgray", colour = "white") + 
  geom_point(data = msat, aes(x = lon, y = lat), size = 4, inherit.aes = FALSE)
#msat_plot
msat_plot_annotated <- msat_plot + xlab("Longitude (°)") + ylab("Latitude (°)") + theme_bw() + 
  ggtitle("msat locations") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), axis.title = element_text(size = 24, face = "bold"), 
        plot.title = element_text(size = 28, face = "bold"))
msat_plot_annotated

######## Create mtdna map ########

mtdna_plot <- ggplot(geogr_data, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "lightgray", colour = "white") + 
  geom_point(data = mtdna, aes(x = lon, y = lat), size = 4, inherit.aes = FALSE)
#mtdna_plot
mtdna_plot_annotated <- mtdna_plot + xlab("Longitude (°)") + ylab("Latitude (°)") + theme_bw() + 
  ggtitle("mtdna locations") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), axis.title = element_text(size = 24, face = "bold"), 
        plot.title = element_text(size = 28, face = "bold"))
mtdna_plot_annotated

######## Create mtdna map ########

mtdna_plot <- ggplot(geogr_data, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "lightgray", colour = "white") + 
  geom_point(data = mtdna, aes(x = lon, y = lat), size = 4, inherit.aes = FALSE)
#mtdna_plot
mtdna_plot_annotated <- mtdna_plot + xlab("Longitude (°)") + ylab("Latitude (°)") + theme_bw() + 
  ggtitle("mtdna locations") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), axis.title = element_text(size = 24, face = "bold"), 
        plot.title = element_text(size = 28, face = "bold"))
mtdna_plot_annotated

######## Create merged map ########

merged_plot <- ggplot(geogr_data, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "lightgray", colour = "white") + 
  geom_point(data = unique_merged, aes(x = lon, y = lat), size = 2, inherit.aes = FALSE)
#merged_plot
merged_plot_annotated <- merged_plot + xlab("Longitude (°)") + ylab("Latitude (°)") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), axis.title = element_text(size = 24, face = "bold"), 
        plot.title = element_text(size = 28, face = "bold"))
merged_plot_annotated

##########################################################################################################

######## Create msat histograms ########

#lat histogram
msat_lat_histogram <- ggplot(data = msat, aes(msat$lat)) + geom_histogram(binwidth = 10, color = "#00B8E5", fill = "#00B8E5") + 
  xlab("latitude") + ylab("count")
msat_lat_histogram_annotated <- msat_lat_histogram + theme_bw() + 
  scale_x_continuous(limits = c(-80, 80), breaks = c(-80, -60, -40, -20, 0, 20, 40, 60, 80)) + coord_flip() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 40, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 40, color = "black"))
msat_lat_histogram_annotated

#lon histogram
msat_lon_histogram <- ggplot(data = msat, aes(msat$lon)) + geom_histogram(binwidth = 10, color = "#00B8E5", fill = "#00B8E5") + 
  xlab("longitude") + ylab("count")
msat_lon_histogram_annotated <- msat_lon_histogram + theme_bw() + 
  scale_x_continuous(limits = c(-190, 180), breaks = c(-180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180)) +
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 40, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 40, color = "black"))
msat_lon_histogram_annotated

######## Create mtdna histograms ########

#lat histogram
mtdna_lat_histogram <- ggplot(data = mtdna, aes(mtdna$lat)) + geom_histogram(binwidth = 10, color = "#00B8E5", fill = "#00B8E5") + 
  xlab("latitude") + ylab("count")
mtdna_lat_histogram_annotated <- mtdna_lat_histogram + theme_bw() + 
  scale_x_continuous(limits = c(-80, 80), breaks = c(-80, -60, -40, -20, 0, 20, 40, 60, 80)) + coord_flip() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 40, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 40, color = "black"))
mtdna_lat_histogram_annotated

#lon histogram
mtdna_lon_histogram <- ggplot(data = mtdna, aes(mtdna$lon)) + geom_histogram(binwidth = 10, color = "#00B8E5", fill = "#00B8E5") + 
  xlab("longitude") + ylab("count")
mtdna_lon_histogram_annotated <- mtdna_lon_histogram + theme_bw() + 
  scale_x_continuous(limits = c(-190, 180), breaks = c(-180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180)) +
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 40, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 40, color = "black"))
mtdna_lon_histogram_annotated
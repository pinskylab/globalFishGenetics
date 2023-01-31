#################################### Plot sampling locations ###########################################

#set-up

remove(list = ls())

#load libraries
library(maps)
library(mapdata)
library(tidyverse)

#read in data
msat <- read.csv("output/latlon_msat_2023-01-31.csv")
mtdna <- read.csv("output/latlon_mtdna_2023-01-31.csv")
merged <- rbind(msat, mtdna)
  unique_merged <- unique(merged)

#pull world map data
world_map <- map_data("world") %>% 
  filter(! long > 180)

################################################################################################################################################

######## Create msat map ########

msat_plot <- world_map %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, color = "lightgray", fill = "lightgray") +
  expand_limits(x = world_map$long, y = world_map$lat) +  
  #ggplot() + 
  geom_point(data = msat, aes(x = lon, y = lat), color = "#626c76", alpha = 0.5, size = 8, 
             inherit.aes = FALSE) + 
  geom_rug(data = msat, mapping = aes(x = lon), color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  geom_rug(data = msat, mapping = aes(y = lat), color = "#282828", alpha = 0.5, inherit.aes = FALSE)
msat_plot_annotated <- msat_plot + ylim(c(-90, 90)) + theme_map()
msat_plot_annotated

######## Create mtdna map ########

mtdna_plot <- ggplot(world_map, aes(x = long, y = lat)) + 
  geom_polygon(fill = "lightgray", colour = "white") + 
  geom_point(data = mtdna, aes(x = lon, y = lat), size = 4, inherit.aes = FALSE) + geom_rug(data = mtdna)
#mtdna_plot
mtdna_plot_annotated <- mtdna_plot + xlab("Longitude (°)") + ylab("Latitude (°)") + theme_bw() + 
  ggtitle("mtdna locations") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), axis.title = element_text(size = 24, face = "bold"), 
        plot.title = element_text(size = 28, face = "bold"))
mtdna_plot_annotated

######## Create mtdna map ########

mtdna_plot <- world_map %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, color = "lightgray", fill = "lightgray") +
  expand_limits(x = world_map$long, y = world_map$lat) +  
  #ggplot() + 
  geom_point(data = mtdna, aes(x = lon, y = lat), color = "#626c76", alpha = 0.5, size = 8, 
             inherit.aes = FALSE) + 
  geom_rug(data = mtdna, mapping = aes(x = lon), color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  geom_rug(data = mtdna, mapping = aes(y = lat), color = "#282828", alpha = 0.5, inherit.aes = FALSE)
mtdna_plot_annotated <- mtdna_plot + ylim(c(-90, 90)) + theme_map()
mtdna_plot_annotated

######## Create merged map ########

merged_plot <- world_map %>% 
  ggplot(aes(map_id = region)) +
  geom_map(map = world_map, color = "lightgray", fill = "lightgray") +
  expand_limits(x = world_map$long, y = world_map$lat) +  
  #ggplot() + 
  geom_point(data = merged, aes(x = lon, y = lat), color = "#626c76", alpha = 0.5, size = 8, 
             inherit.aes = FALSE) + 
  geom_rug(data = merged, mapping = aes(x = lon), color = "#282828", alpha = 0.5, inherit.aes = FALSE) + 
  geom_rug(data = merged, mapping = aes(y = lat), color = "#282828", alpha = 0.5, inherit.aes = FALSE)
merged_plot_annotated <- merged_plot + ylim(c(-90, 90)) + theme_map()
merged_plot_annotated


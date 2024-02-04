################################################### Script to plot sampling locations ########################################################

#Plot either each sampled population as a point OR combine as hexbins (500x500 km equal-area grid cells)
#Mollweide projection
#Maps of global chlorophyll & temp as well

#Plot size: 6000 x 4000

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(sf) #v.1.0.13
library(tmap) #3.3.4
library(tidyverse) #v.2.0.0
library(data.table) #v.1.14.8

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv")
msat <- read.csv("output/msatloci_assembled.csv")
msat_env <- read.csv("output/msat_env.csv")

#merge dataframes
msat <- cbind(msat[, -1], msat_env[, c('sst.BO_sstmean', 'sst.BO_sstrange', 'sst.BO_sstmax', 'sst.BO_sstmin',
                                       'BO_dissox', 'chloroA.BO_chlomean', 'chloroA.BO_chlorange',
                                       'chloroA.BO_chlomax', 'chloroA.BO_chlomin')]) #merge not working for some reason, cbind bc in same order

#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)
  mtdna_small_hd <- subset(mtdna_small, mtdna_small$He != "NA")
  mtdna_small_pi <- subset(mtdna_small, mtdna_small$Pi != "NA")

#pull world map data
data("World") # load the background world outline
  World2 <- st_transform(World, "ESRI:54009") #project world outline to Mollweide

#ready Mollweide projection
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #the current coordinate reference systems (lat-lon)

################################################################################################################################################

######## Calculate tropical coverage ########
#want unique locations -- lat lon combinations/species

msat_tropics <- msat %>%
  filter(lat >= -23.4 & lat <= 23.4) %>%
  distinct(spp, lat, lon, .keep_all = TRUE)

mtdna_tropics <- mtdna_small %>%
  filter(lat >= -23.4 & lat <= 23.4) %>%
  distinct(spp, lat, lon, .keep_all = TRUE)

################################################################################################################################################

######## Calculate polar coverage ########
#want unique locations -- lat lon combinations/species

msat_polar <- msat %>%
  filter(lat <= -60 | lat >= 60) %>%
  distinct(spp, lat, lon, .keep_all = TRUE)

mtdna_polar <- mtdna_small %>%
  filter(lat <= -60 | lat >= 60) %>%
  distinct(spp, lat, lon, .keep_all = TRUE)

################################################################################################################################################

######## Create msat maps ########

#project lat & lon to Mollweide
msat_st <- st_as_sf(msat, coords = c("lon", "lat"), remove = FALSE, crs = projcrs) #make spatial object
msat_st_Moll <- st_transform(msat_st, "ESRI:54009") #project to Mollweide

#extract Mollweide coordinates in a data.frame
msat_coords <- data.frame(st_coordinates(msat_st_Moll)) #get X & Y coords from the projected points in Mollweide
msat_coords$He <- msat_st_Moll$He #copy over genetic diversity estimates

#points only map
msat_he_point_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  geom_point(data = msat_coords, aes(x = X, y = Y, color = He), 
             alpha = 0.5, size = 16, inherit.aes = FALSE)  + 
  geom_rug(data = msat_coords, mapping = aes(x = X), color = "#282828", 
           length = unit(0.03, "npc"), inherit.aes = FALSE) + 
  geom_rug(data = msat_coords, mapping = aes(y = Y), color = "#282828", 
           length = unit(0.015, "npc"), inherit.aes = FALSE) + #projection messes these lengths up - looks like y is 2x as big as x
  scale_color_gradient(low = "#1f2c32", high = "#e3ebed") + 
  annotate("text", x = -16000000, y = 8500000, label = "(c)", size = 100) +
  labs(color = "He") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(color = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
msat_he_point_plot

#hex bin mean map
msat_he_mean_hexbin_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  stat_summary_hex(data = msat_coords, aes(x = X, y = Y, z = He), 
                   fun = function(x) mean(x), binwidth = c(5e+05, 5e+05), #bin width and height in meters
                   inherit.aes = FALSE) + 
  geom_rug(data = msat_coords, mapping = aes(x = X), color = "#282828", 
           length = unit(0.03, "npc"), inherit.aes = FALSE) + 
  geom_rug(data = msat_coords, mapping = aes(y = Y), color = "#282828", 
           length = unit(0.015, "npc"), inherit.aes = FALSE) + #projection messes these lengths up - looks like y is 2x as big as x
  scale_fill_gradient(low = "#1f2c32", high = "#e3ebed") + 
  annotate("text", x = -16000000, y = 8500000, label = "(c)", size = 100) +
  labs(fill = "He") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
msat_he_mean_hexbin_plot

#hex bin SD map
msat_he_SD_hexbin_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  stat_summary_hex(data = msat_coords, aes(x = X, y = Y, z = He), 
                   fun = function(x) sd(x), binwidth = c(5e+05, 5e+05), #bin width and height in meters
                   inherit.aes = FALSE) + 
  scale_fill_gradient(low = "#1f2c32", high = "#e3ebed") + 
  annotate("text", x = -16000000, y = 8500000, label = "(c)", size = 100) +
  labs(fill = "SD") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
msat_he_SD_hexbin_plot

#hex bin count map
msat_he_count_hexbin_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  stat_summary_hex(data = msat_coords, aes(x = X, y = Y, z = He), 
                   fun = function(x) length(x), binwidth = c(5e+05, 5e+05), #bin width and height in meters
                   inherit.aes = FALSE) + 
  scale_fill_gradient(low = "#1f2c32", high = "#e3ebed") + 
  annotate("text", x = -16000000, y = 8500000, label = "(c)", size = 100) +
  labs(fill = "Count") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
msat_he_count_hexbin_plot
  
###################################################################################################################################

######## Create mtdna maps ########

#### mtdna hd maps ####  

#project lat & lon to Mollweide
mtdna_small_hd_st <- st_as_sf(mtdna_small_hd, coords = c("lon", "lat"), remove = FALSE, crs = projcrs) #make spatial object
mtdna_small_hd_st_Moll <- st_transform(mtdna_small_hd_st, "ESRI:54009") #project to Mollweide

#extract Mollweide coordinates in a data.frame
mtdna_small_hd_coords <- data.frame(st_coordinates(mtdna_small_hd_st_Moll)) #get X & Y coords from the projected points in Mollweide
mtdna_small_hd_coords$He <- mtdna_small_hd_st_Moll$He #copy over genetic diversity estimates

#points only map
mtdna_hd_point_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  geom_point(data = mtdna_small_hd_coords, aes(x = X, y = Y, color = He), 
             alpha = 0.5, size = 16, inherit.aes = FALSE)  + 
  geom_rug(data = mtdna_small_hd_coords, mapping = aes(x = X), color = "#282828", 
           length = unit(0.03, "npc"), inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_hd_coords, mapping = aes(y = Y), color = "#282828", 
           length = unit(0.015, "npc"), inherit.aes = FALSE) + #projection messes these lengths up - looks like y is 2x as big as x
  scale_color_gradient(low = "#1f2c32", high = "#e3ebed") + 
  annotate("text", x = -16000000, y = 8500000, label = "(b)", size = 100) +
  labs(color = "Hd") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(color = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
mtdna_hd_point_plot
  
#hex bin mean map
mtdna_hd_mean_hexbin_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  stat_summary_hex(data = mtdna_small_hd_coords, aes(x = X, y = Y, z = He), 
                   fun = function(x) mean(x), binwidth = c(5e+05, 5e+05), #bin width and height in meters
                   inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_hd_coords, mapping = aes(x = X), color = "#282828", 
           length = unit(0.03, "npc"), inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_hd_coords, mapping = aes(y = Y), color = "#282828", 
           length = unit(0.015, "npc"), inherit.aes = FALSE) + #projection messes these lengths up - looks like y is 2x as big as x
  scale_fill_gradient(low = "#1f2c32", high = "#e3ebed") + 
  annotate("text", x = -16000000, y = 8500000, label = "(b)", size = 100) +
  labs(fill = "Hd") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
mtdna_hd_mean_hexbin_plot

#hex bin SD map
mtdna_hd_SD_hexbin_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  stat_summary_hex(data = mtdna_small_hd_coords, aes(x = X, y = Y, z = He), 
                   fun = function(x) sd(x), binwidth = c(5e+05, 5e+05), #bin width and height in meters
                   inherit.aes = FALSE) + 
  scale_fill_gradient(low = "#1f2c32", high = "#e3ebed") + 
  annotate("text", x = -16000000, y = 8500000, label = "(b)", size = 100) +
  labs(fill = "SD") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
mtdna_hd_SD_hexbin_plot

#hex bin count map
mtdna_hd_count_hexbin_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  stat_summary_hex(data = mtdna_small_hd_coords, aes(x = X, y = Y, z = He), 
                   fun = function(x) length(x), binwidth = c(5e+05, 5e+05), #bin width and height in meters
                   inherit.aes = FALSE) + 
  scale_fill_gradient(low = "#1f2c32", high = "#e3ebed") + 
  annotate("text", x = -16000000, y = 8500000, label = "(b)", size = 100) +
  labs(fill = "Count") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
mtdna_hd_count_hexbin_plot

#### mtdna pi maps ####

#project lat & lon to Mollweide
mtdna_small_pi_st <- st_as_sf(mtdna_small_pi, coords = c("lon", "lat"), remove = FALSE, crs = projcrs) #make spatial object
mtdna_small_pi_st_Moll <- st_transform(mtdna_small_pi_st, "ESRI:54009") #project to Mollweide

#extract Mollweide coordinates in a data.frame
mtdna_small_pi_coords <- data.frame(st_coordinates(mtdna_small_pi_st_Moll)) #get X & Y coords from the projected points in Mollweide
mtdna_small_pi_coords$Pi <- mtdna_small_pi_st_Moll$Pi #copy over genetic diversity estimates

#points only map
mtdna_pi_point_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  geom_point(data = mtdna_small_pi_coords, aes(x = X, y = Y, color = Pi), 
             alpha = 0.5, size = 16, inherit.aes = FALSE)  + 
  geom_rug(data = mtdna_small_pi_coords, mapping = aes(x = X), color = "#282828", 
           length = unit(0.03, "npc"), inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_pi_coords, mapping = aes(y = Y), color = "#282828", 
           length = unit(0.015, "npc"), inherit.aes = FALSE) + #projection messes these lengths up - looks like y is 2x as big as x
  scale_color_gradient(trans = "log10", low = "#1f2c32", high = "#e3ebed", 
                       labels = scales::label_log()) + 
  annotate("text", x = -16000000, y = 8500000, label = "(a)", size = 100) +
  labs(color = "π") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 110), 
        legend.title = element_text(size = 110), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(color = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
mtdna_pi_point_plot

#hex bin mean map
mtdna_pi_mean_hexbin_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  stat_summary_hex(data = mtdna_small_pi_coords, aes(x = X, y = Y, z = Pi), 
                   fun = function(x) mean(x), binwidth = c(5e+05, 5e+05), #bin width and height in meters
                   inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_pi_coords, mapping = aes(x = X), color = "#282828", 
           length = unit(0.03, "npc"), inherit.aes = FALSE) + 
  geom_rug(data = mtdna_small_pi_coords, mapping = aes(y = Y), color = "#282828", 
           length = unit(0.015, "npc"), inherit.aes = FALSE) + #projection messes these lengths up - looks like y is 2x as big as x
  scale_fill_gradient(trans = "log10", low = "#1f2c32", high = "#e3ebed", 
                      labels = scales::label_log()) + 
  annotate("text", x = -16000000, y = 8500000, label = "(a)", size = 100) +
  labs(fill = "π") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
mtdna_pi_mean_hexbin_plot

#hex bin SD map
mtdna_pi_SD_hexbin_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  stat_summary_hex(data = mtdna_small_pi_coords, aes(x = X, y = Y, z = Pi), 
                   fun = function(x) sd(x), binwidth = c(5e+05, 5e+05), #bin width and height in meters
                   inherit.aes = FALSE) + 
  scale_fill_gradient(low = "#1f2c32", high = "#e3ebed") + 
  annotate("text", x = -16000000, y = 8500000, label = "(a)", size = 100) +
  labs(fill = "SD") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
mtdna_pi_SD_hexbin_plot

#hex bin count map
mtdna_pi_count_hexbin_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  stat_summary_hex(data = mtdna_small_pi_coords, aes(x = X, y = Y, z = Pi), 
                   fun = function(x) length(x), binwidth = c(5e+05, 5e+05), #bin width and height in meters
                   inherit.aes = FALSE) + 
  scale_fill_gradient(low = "#1f2c32", high = "#e3ebed") + 
  annotate("text", x = -16000000, y = 8500000, label = "(a)", size = 100) +
  labs(fill = "Count") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
mtdna_pi_count_hexbin_plot

######################################################################################################################################

######## Create environmental data maps ########

#extract Mollweide coordinates, using He database
msat_coords$chloromean <- msat_st_Moll$chloroA.BO_chlomean #copy over chlorophyll estimates
msat_coords$sstmean <- msat_st_Moll$sst.BO_sstmean #copy over SST estimates

#chlorophyll mean map
chlomean_mean_hexbin_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  stat_summary_hex(data = msat_coords, aes(x = X, y = Y, z = chloromean), 
                   fun = function(x) mean(x), binwidth = c(5e+05, 5e+05), #bin width and height in meters
                   inherit.aes = FALSE) + 
  geom_rug(data = msat_coords, mapping = aes(x = X), color = "#282828", 
           length = unit(0.03, "npc"), inherit.aes = FALSE) + 
  geom_rug(data = msat_coords, mapping = aes(y = Y), color = "#282828", 
           length = unit(0.015, "npc"), inherit.aes = FALSE) + #projection messes these lengths up - looks like y is 2x as big as x
  scale_fill_gradient(trans = "log10", low = "#c6cbc8", high = "#212924", 
                      labels = scales::label_log()) +
  annotate("text", x = -16000000, y = 8500000, label = "(a)", size = 100) +
  labs(fill = "Chlorophyll") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
chlomean_mean_hexbin_plot

#SST mean map
sstmean_mean_hexbin_plot <- ggplot(World2) +
  geom_sf(color = "lightgrey", fill = "lightgrey") + 
  stat_summary_hex(data = msat_coords, aes(x = X, y = Y, z = sstmean), 
                   fun = function(x) mean(x), binwidth = c(5e+05, 5e+05), #bin width and height in meters
                   inherit.aes = FALSE) + 
  geom_rug(data = msat_coords, mapping = aes(x = X), color = "#282828", 
           length = unit(0.03, "npc"), inherit.aes = FALSE) + 
  geom_rug(data = msat_coords, mapping = aes(y = Y), color = "#282828", 
           length = unit(0.015, "npc"), inherit.aes = FALSE) + #projection messes these lengths up - looks like y is 2x as big as x
  scale_fill_gradient(low = "#ebcb99", high = "#905800") +
  annotate("text", x = -16000000, y = 8500000, label = "(b)", size = 100) +
  labs(fill = "SST") + 
  theme_bw() + 
  theme(plot.margin = unit(c(8, 2, 10, 2), "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "black"), 
        axis.line = element_line(linewidth = 4, color = "black"), 
        axis.title = element_blank(), 
        legend.position = "right", 
        legend.text = element_text(size = 120), 
        legend.title = element_text(size = 120), 
        legend.justification = "center", 
        legend.title.align = 0.1) + 
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"), barheight = unit(18, "cm")))
sstmean_mean_hexbin_plot

#### SST mean & chlorophyll mean correlation plot ####

#corr coefficient (using msat bc largest database)
corr <- cor(x = msat$sst.BO_sstmean, y = msat$chloroA.BO_chlomean, 
            use = "complete.obs", method = "spearman") #r= -0.316

#plot
SST_chlo_plot <- ggplot() + 
  geom_point(data = msat, aes(x = sst.BO_sstmean, y = chloroA.BO_chlomean), 
             alpha = 0.5, size = 8, color = "#0E3C45") + 
  annotate("text", x = 28, y = 50, label = "r = -0.316", size = 12) +
  xlab("SST mean") + ylab("chlorophyll A mean") + 
  theme_minimal() + 
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 34),
        axis.title.y = element_text(size = 34),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"))
SST_chlo_plot
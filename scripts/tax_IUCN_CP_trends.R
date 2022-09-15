#Script for cleaning/investigating datasets

remove(list = ls())

library(tidyverse)
library(gridExtra)
library(plotrix)
library(here)
library(data.table)

#read in data
mtdna <- read.csv(here("output", "mtdna_assembled.csv"), stringsAsFactors = FALSE)
msat <- read.csv(here("output", "msatloci_assembled.csv"), stringsAsFactors = FALSE)
cp_info <- read.csv(here("output", "spp_combined_info.csv"), stringsAsFactors = FALSE)
  mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'IUCN_status', 'Genus', 'Family', 'Order', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                  'Centroid', 'IUCN_status')], all.x = TRUE)
  msat <- merge(msat, cp_info[, c('spp', 'Pelagic_Coastal', 'IUCN_status', 'Genus', 'Family', 'Order', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                'Centroid', 'IUCN_status')], all.x = TRUE)

#calculate abslat
mtdna$abslat <- abs(mtdna$lat)
msat$abslat <- abs(msat$lat)

#clean up dataframes
#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna <-subset(mtdna, as.numeric(mtdna$bp) < 2000)

#create mtdna pi dataset
mtdna_pi <- subset(mtdna, mtdna$Pi != "NA")
mtdna_pi <- subset(mtdna_pi, mtdna_pi$Pi != 0)

#create mtdna he dataset
mtdna_he <- subset(mtdna, mtdna$He != "NA")
mtdna_he$success <- round(mtdna_he$He*mtdna_he$n)
mtdna_he$failure<- round((1 - mtdna_he$He)*mtdna_he$n)

#create msat he dataset
msat$n <- as.numeric(msat$n)
  msat <- subset(msat, msat$n != "NA")
msat$success <- round(msat$He*msat$n)
msat$failure<- round((1 - msat$He)*msat$n)

################################################################################################

######## Coastal v. pelagic trends ########

mtdna_pi_abslat_plot_CP <- ggplot(data = mtdna_pi, aes(x = abslat, y = Pi, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("absolute latitude") + ylab("mtdna pi") + scale_color_manual(values = c("#999999", "#3333CC"))
mtdna_pi_abslat_plot_annotated_CP <- mtdna_pi_abslat_plot_CP + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
       axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
       legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_pi_abslat_plot_annotated_CP

################################################################################################

######## Order trends ########

#### mtdna pi trends ####

mtdna_pi_order <- c(unique(mtdna_pi$Order)) #27 orders

#how many observations for each order? really want to look at trends in the ones with enough data
mtdna_pi <- as.data.table(mtdna_pi)
pi_order <- mtdna_pi[, .N, by = .(Order)] #keep anything with over 50 observations 
#Perciformes, Scorpaeniformes, Beryciformes, Gadiformes, Carcharhiniformes, Pleuronectiformes, Clupeiformes

mtdna_pi_subset <- subset(mtdna_pi, mtdna_he$Order == "Perciformes" | mtdna_pi$Order == "Scorpaeniformes" | 
                            mtdna_pi$Order == "Beryciformes" | mtdna_pi$Order == "Gadiformes" | 
                            mtdna_pi$Order == "Carcharhiniformes" | mtdna_pi$Order == "Pleuronectiformes" | 
                            mtdna_pi$Order == "Clupeiformes" | mtdna_pi$Order == "Rajiformes" | 
                            mtdna_pi$Order == "Syngnathiformes")

mtdna_pi_abslat_plot_Order <- ggplot(data = mtdna_pi_subset, aes(x = abslat, y = Pi, color = Order)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = Order)) + 
  xlab("absolute latitude") + ylab("mtdna pi")
mtdna_pi_abslat_plot_annotated_Order <- mtdna_pi_abslat_plot_Order + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_pi_abslat_plot_annotated_Order

mtdna_pi_lat_plot_Order <- ggplot(data = mtdna_pi_subset, aes(x = lat, y = Pi, color = Order)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = Order)) + 
  xlab("latitude") + ylab("mtdna pi")
mtdna_pi_lat_plot_annotated_Order <- mtdna_pi_lat_plot_Order + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_pi_lat_plot_annotated_Order

mtdna_pi_lon_plot_Order <- ggplot(data = mtdna_pi_subset, aes(x = lon, y = Pi, color = Order)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = Order)) + 
  xlab("longitude") + ylab("mtdna pi")
mtdna_pi_lon_plot_annotated_Order <- mtdna_pi_lon_plot_Order + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_pi_lon_plot_annotated_Order

#### mtdna he trends ####

mtdna_he_order <- c(unique(mtdna_he$Order)) #28 orders

#how many observations for each order? really want to look at trends in the ones with enough data
mtdna_he <- as.data.table(mtdna_he)
mtdna_he_order <- mtdna_he[, .N, by = .(Order)] #keep anything with over 50 observations 
#Perciformes, Rajiformes, Scorpaeniformes, Beryciformes, Gadiformes, Carcharhiniformes, Pleuronectiformes, Clupeiformes, Syngnathiformes

mtdna_he_subset <- subset(mtdna_he, mtdna_he$Order == "Perciformes" | mtdna_he$Order == "Scorpaeniformes" | 
                           mtdna_he$Order == "Beryciformes" | mtdna_he$Order == "Gadiformes" | 
                           mtdna_he$Order == "Carcharhiniformes" | mtdna_he$Order == "Pleuronectiformes" | 
                           mtdna_he$Order == "Clupeiformes" | mtdna_he$Order == "Rajiformes" | 
                           mtdna_he$Order == "Syngnathiformes")

mtdna_he_abslat_plot_Order <- ggplot(data = mtdna_he_subset, aes(x = abslat, y = He, succ = success, fail = failure, color = Order)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2,  aes(linetype = Order)) + 
  xlab("absolute latitude") + ylab("mtdna he")
mtdna_he_abslat_plot_annotated_Order <- mtdna_he_abslat_plot_Order + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_he_abslat_plot_annotated_Order

mtdna_he_lat_plot_Order <- ggplot(data = mtdna_he_subset, aes(x = lat, y = He, succ = success, fail = failure, color = Order)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Order)) + 
  xlab("latitude") + ylab("mtdna he")
mtdna_he_lat_plot_annotated_Order <- mtdna_he_lat_plot_Order + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_he_lat_plot_annotated_Order

mtdna_he_lon_plot_Order <- ggplot(data = mtdna_he_subset, aes(x = lon, y = He, succ = success, fail = failure, color = Order)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Order)) + 
  xlab("longitude") + ylab("mtdna he")
mtdna_he_lon_plot_annotated_Order <- mtdna_he_lon_plot_Order + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_he_lon_plot_annotated_Order

#### msat he trends ####

msat_he_order <- c(unique(msat$Order)) #27 orders

#how many observations for each order? really want to look at trends in the ones with enough data
msat <- as.data.table(msat)
msat_he_order <- msat[, .N, by = .(Order)] #keep same orders as with mtdna dataset
#Perciformes, Rajiformes, Scorpaeniformes, Beryciformes, Gadiformes, Carcharhiniformes, Pleuronectiformes, Clupeiformes, Syngnathiformes

msat_he_subset <- subset(msat, msat$Order == "Perciformes" | msat$Order == "Scorpaeniformes" | 
                            msat$Order == "Beryciformes" | msat$Order == "Gadiformes" | 
                            msat$Order == "Carcharhiniformes" | msat$Order == "Pleuronectiformes" | 
                            msat$Order == "Clupeiformes" | msat$Order == "Rajiformes" | 
                            msat$Order == "Syngnathiformes")

msat_he_abslat_plot_Order <- ggplot(data = msat_he_subset, aes(x = abslat, y = He, succ = success, fail = failure, color = Order)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Order)) + 
  xlab("absolute latitude") + ylab("msat he")
msat_he_abslat_plot_annotated_Order <- msat_he_abslat_plot_Order + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_he_abslat_plot_annotated_Order

msat_he_lat_plot_Order <- ggplot(data = msat_he_subset, aes(x = lat, y = He, succ = success, fail = failure, color = Order)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Order)) + 
  xlab("latitude") + ylab("msat he")
msat_he_lat_plot_annotated_Order <- msat_he_lat_plot_Order + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_he_lat_plot_annotated_Order

msat_he_lon_plot_Order <- ggplot(data = msat_he_subset, aes(x = lon, y = He, succ = success, fail = failure, color = Order)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Order)) + 
  xlab("longitude") + ylab("msat he")
msat_he_lon_plot_annotated_Order <- msat_he_lon_plot_Order + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_he_lon_plot_annotated_Order

################################################################################################

######## IUCN trends ########

#### mtdna pi trends ####

#how many observations for each IUCN status
pi_IUCN <- mtdna_pi[, .N, by = .(IUCN_status)]

#plots
mtdna_pi_abslat_plot_IUCN <- ggplot(data = mtdna_pi, aes(x = abslat, y = Pi, color = IUCN_status)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = IUCN_status)) + 
  xlab("absolute latitude") + ylab("mtdna pi")
mtdna_pi_abslat_plot_annotated_IUCN <- mtdna_pi_abslat_plot_IUCN + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_pi_abslat_plot_annotated_IUCN

mtdna_pi_lat_plot_IUCN <- ggplot(data = mtdna_pi, aes(x = lat, y = Pi, color = IUCN_status)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = IUCN_status)) + 
  xlab("latitude") + ylab("mtdna pi")
mtdna_pi_lat_plot_annotated_IUCN <- mtdna_pi_lat_plot_IUCN + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_pi_lat_plot_annotated_IUCN

mtdna_pi_lon_plot_IUCN <- ggplot(data = mtdna_pi, aes(x = lon, y = Pi, color = IUCN_status)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = IUCN_status)) + 
  xlab("longitude") + ylab("mtdna pi")
mtdna_pi_lon_plot_annotated_IUCN <- mtdna_pi_lon_plot_IUCN + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_pi_lon_plot_annotated_IUCN

#### mtdna he trends ####

#how many observations for each IUCN status
he_IUCN <- mtdna_he[, .N, by = .(IUCN_status)]

mtdna_he_abslat_plot_IUCN <- ggplot(data = mtdna_he, aes(x = abslat, y = He, succ = success, fail = failure, color = IUCN_status)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = IUCN_status)) + 
  xlab("absolute latitude") + ylab("mtdna he")
mtdna_he_abslat_plot_annotated_IUCN <- mtdna_he_abslat_plot_IUCN + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_he_abslat_plot_annotated_IUCN

mtdna_he_lat_plot_IUCN <- ggplot(data = mtdna_he, aes(x = lat, y = He, succ = success, fail = failure, color = IUCN_Status)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = IUCN_status)) + 
  xlab("latitude") + ylab("mtdna he")
mtdna_he_lat_plot_annotated_IUCN <- mtdna_he_lat_plot_IUCN + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_he_lat_plot_annotated_IUCN

mtdna_he_lon_plot_IUCN <- ggplot(data = mtdna_he, aes(x = lon, y = He, succ = success, fail = failure, color = IUCN_status)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = IUCN_status)) + 
  xlab("longitude") + ylab("mtdna he")
mtdna_he_lon_plot_annotated_IUCN <- mtdna_he_lon_plot_IUCN + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_he_lon_plot_annotated_IUCN

#### msat he trends ####

#how many observations for each IUCN status
msat_he_IUCN <- msat[, .N, by = .(IUCN_status)]

msat_he_abslat_plot_IUCN <- ggplot(data = msat_he_subset, aes(x = abslat, y = He, succ = success, fail = failure, color = IUCN_status)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = IUCN_status)) + 
  xlab("absolute latitude") + ylab("msat he")
msat_he_abslat_plot_annotated_IUCN <- msat_he_abslat_plot_IUCN + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_he_abslat_plot_annotated_IUCN

msat_he_lat_plot_IUCN <- ggplot(data = msat_he_subset, aes(x = lat, y = He, succ = success, fail = failure, color = IUCN_status)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = IUCN_status)) + 
  xlab("latitude") + ylab("msat he")
msat_he_lat_plot_annotated_IUCN <- msat_he_lat_plot_IUCN + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_he_lat_plot_annotated_IUCN

msat_he_lon_plot_IUCN <- ggplot(data = msat_he_subset, aes(x = lon, y = He, succ = success, fail = failure, color = IUCN_status)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = IUCN_status)) + 
  xlab("longitude") + ylab("msat he")
msat_he_lon_plot_annotated_IUCN <- msat_he_lon_plot_IUCN + theme_bw() + ylim(0, 1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_he_lon_plot_annotated_IUCN

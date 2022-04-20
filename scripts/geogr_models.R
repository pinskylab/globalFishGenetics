remove(list = ls())

library(tidyverse)
library(gridExtra)

#read in data
mtdna <- read.csv("Output/mtdna_assembled.csv", stringsAsFactors = FALSE)
msat <- read.csv("Output/msatloci_assembled.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("Output/spp_combined_info.csv", stringsAsFactors = FALSE)


mtdna$abslat <- abs(mtdna$lat)
mtdna$abslon <- abs(mtdna$lon)

msat$abslat <- abs(msat$lat)
msat$abslon <- abs(msat$lon)

msat_cp <- merge(msat, cp_info[, c('spp', 'Pelagic_Coastal')], all.x = TRUE)
#msat_cp$Pelagic.Coastal.f <- factor(msat_cp$Pelagic.Coastal)
#is.factor(msat_cp$Pelagic.Coastal.f)

mtdna_cp <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal')], all.x = TRUE)
#mtdna_cp$Pelagic.Coastal.f <- factor(mtdna_cp$Pelagic.Coastal)
#is.factor(mtdna_cp$Pelagic.Coastal.f)

######################################################################################################

######## mtdna Pi plots ########

mtdna_nonapi_CP <- subset(mtdna_cp, mtdna_cp$Pi != "NA")
mtdna_nonapi_small_CP <- subset(mtdna_nonapi_CP, mtdna_nonapi_CP$Pi < 1)

#plots
#mtdna_nonapi_small_abslat_plot_CP <- ggplot(data = mtdna_nonapi_small_CP, aes(x = abslat, y = Pi, color = Pelagic_Coastal)) + 
#  geom_point(alpha = 0.3, size = 3) + 
#  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
#  xlab("absolute latitude") + ylab("mtdna pi") + scale_color_manual(values = c("#999999", "#3333CC"))
#mtdna_nonapi_small_abslat_plot_annotated_CP <- mtdna_nonapi_small_abslat_plot_CP + theme_bw() + ylim(0, 0.2) + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
#        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
#mtdna_nonapi_small_abslat_plot_annotated_CP

mtdna_nonapi_small_abslat_plot <- ggplot(data = mtdna_nonapi_small_CP, aes(x = abslat, y = Pi)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2) + 
  xlab("absolute latitude") + ylab("mtdna pi")
mtdna_nonapi_small_abslat_plot_annotated <- mtdna_nonapi_small_abslat_plot + theme_bw() + ylim(0, 0.1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonapi_small_abslat_plot_annotated

#mtdna_nonapi_small_lat_plot_CP <- ggplot(data = mtdna_nonapi_small_CP, aes(x = lat, y = Pi, color = Pelagic_Coastal)) + 
#  geom_point(alpha = 0.3, size = 3) + 
#  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
#  xlab("latitude") + ylab("mtdna pi") + scale_color_manual(values = c("#999999", "#3333CC"))
#mtdna_nonapi_small_lat_plot_annotated_CP <- mtdna_nonapi_small_lat_plot_CP + theme_bw() + ylim(0, 0.2) + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
#        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
#mtdna_nonapi_small_lat_plot_annotated_CP

mtdna_nonapi_small_lat_plot <- ggplot(data = mtdna_nonapi_small_CP, aes(x = lat, y = Pi)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2) + 
  xlab("latitude") + ylab("mtdna pi")
mtdna_nonapi_small_lat_plot_annotated <- mtdna_nonapi_small_lat_plot + theme_bw() + ylim(0, 0.1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonapi_small_lat_plot_annotated

#mtdna_nonapi_small_lon_plot_CP <- ggplot(data = mtdna_nonapi_small_CP, aes(x = lon, y = Pi, color = Pelagic_Coastal)) + 
#  geom_point(alpha = 0.3, size = 3) + 
#  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
#  xlab("longitude") + ylab("mtdna pi") + scale_color_manual(values = c("#999999", "#3333CC"))
#mtdna_nonapi_small_lon_plot_annotated_CP <- mtdna_nonapi_small_lon_plot_CP + theme_bw() + ylim(0, 0.2) + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"),
#        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
#mtdna_nonapi_small_lon_plot_annotated_CP

mtdna_nonapi_small_lon_plot <- ggplot(data = mtdna_nonapi_small_CP, aes(x = lon, y = Pi)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2) + 
  xlab("longitude") + ylab("mtdna pi")
mtdna_nonapi_small_lon_plot_annotated <- mtdna_nonapi_small_lon_plot + theme_bw() + ylim(0, 0.1) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"),
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonapi_small_lon_plot_annotated

######################################################################################################

######## mtdna He plots ########

mtdna_nonahe_CP <- subset(mtdna_cp, mtdna_cp$He != "NA")
mtdna_nonahe_CP$success <- round(mtdna_nonahe_CP$He*mtdna_nonahe_CP$n)
mtdna_nonahe_CP$failure <- round((1 - mtdna_nonahe_CP$He)*mtdna_nonahe_CP$n)


#plots
#mtdna_nonahe_abslat_plot_CP <- ggplot(data = mtdna_nonahe_CP, aes(x = abslat, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
#  geom_point(alpha = 0.3, size = 3) + 
#  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
#  xlab("absolute latitude") + ylab("mtdna hd") + scale_color_manual(values = c("#999999", "#3333CC"))
#mtdna_nonahe_abslat_plot_annotated_CP <- mtdna_nonahe_abslat_plot_CP + theme_bw() + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
#        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
#mtdna_nonahe_abslat_plot_annotated_CP

mtdna_nonahe_abslat_plot <- ggplot(data = mtdna_nonahe_CP, aes(x = abslat, y = He, succ = success, fail = failure)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2) + 
  xlab("absolute latitude") + ylab("mtdna hd")
mtdna_nonahe_abslat_plot_annotated <- mtdna_nonahe_abslat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonahe_abslat_plot_annotated

#mtdna_nonahe_lat_plot_CP <- ggplot(data = mtdna_nonahe_CP, aes(x = lat, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
#  geom_point(alpha = 0.3, size = 3) + 
#  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
#  xlab("latitude") + ylab("mtdna hd") + scale_color_manual(values = c("#999999", "#3333CC"))
#mtdna_nonahe_lat_plot_annotated_CP <- mtdna_nonahe_lat_plot_CP + theme_bw() + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
#        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
#mtdna_nonahe_lat_plot_annotated_CP

mtdna_nonahe_lat_plot <- ggplot(data = mtdna_nonahe_CP, aes(x = lat, y = He, succ = success, fail = failure)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2) + 
  xlab("latitude") + ylab("mtdna hd")
mtdna_nonahe_lat_plot_annotated <- mtdna_nonahe_lat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonahe_lat_plot_annotated

#mtdna_nonahe_lon_plot_CP <- ggplot(data = mtdna_nonahe_CP, aes(x = lon, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
#  geom_point(alpha = 0.3, size = 3) + 
#  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
#  xlab("longitude") + ylab("mtdna hd") + scale_color_manual(values = c("#999999", "#3333CC"))
#mtdna_nonahe_lon_plot_annotated_CP <- mtdna_nonahe_lon_plot_CP + theme_bw() + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
#        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
#mtdna_nonahe_lon_plot_annotated_CP

mtdna_nonahe_lon_plot <- ggplot(data = mtdna_nonahe_CP, aes(x = lon, y = He, succ = success, fail = failure)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2) + 
  xlab("longitude") + ylab("mtdna hd")
mtdna_nonahe_lon_plot_annotated <- mtdna_nonahe_lon_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonahe_lon_plot_annotated

######################################################################################################

######## msat He plots ########

msat_nonan <- subset(msat_cp, msat_cp$n != "NA" & msat_cp$n != "25-32")
msat_nonan$n <- as.numeric(msat_nonan$n)
  class(msat_nonan$n) #check that is numeric
msat_nonan$success <- round(msat_nonan$He*msat_nonan$n)
msat_nonan$failure <- round((1 - msat_nonan$He)*msat_nonan$n)

#plots
#msat_nonan_abslat_plot_CP <- ggplot(data = msat_nonan, aes(x = abslat, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
#  geom_point(alpha = 0.3, size = 1) + 
#  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
#  xlab("absolute latitude") + ylab("msat he") + scale_color_manual(values = c("#999999", "#3333CC"))
#msat_nonan_abslat_plot_annotated_CP <- msat_nonan_abslat_plot_CP + theme_bw() + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
#        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
#msat_nonan_abslat_plot_annotated_CP

msat_nonan_abslat_plot <- ggplot(data = msat_nonan, aes(x = abslat, y = He, succ = success, fail = failure)) + 
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2) + 
  xlab("absolute latitude") + ylab("msat he")
msat_nonan_abslat_plot_annotated <- msat_nonan_abslat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_nonan_abslat_plot_annotated

#msat_nonan_lat_plot_CP <- ggplot(data = msat_nonan, aes(x = lat, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
#  geom_point(alpha = 0.3, size = 1) + 
#  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
#  xlab("latitude") + ylab("msat he") + scale_color_manual(values = c("#999999", "#3333CC"))
#msat_nonan_lat_plot_annotated_CP <- msat_nonan_lat_plot_CP + theme_bw() + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
#        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
#msat_nonan_lat_plot_annotated_CP

msat_nonan_lat_plot <- ggplot(data = msat_nonan, aes(x = lat, y = He, succ = success, fail = failure)) + 
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2) + 
  xlab("latitude") + ylab("msat he")
msat_nonan_lat_plot_annotated <- msat_nonan_lat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_nonan_lat_plot_annotated

#msat_nonan_lon_plot_CP <- ggplot(data = msat_nonan, aes(x = lon, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
#  geom_point(alpha = 0.3, size = 1) + 
#  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
#  xlab("longitude") + ylab("msat he") + scale_color_manual(values = c("#999999", "#3333CC"))
#msat_nonan_lon_plot_annotated_CP <- msat_nonan_lon_plot_CP + theme_bw() + 
#  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
#        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
#        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
#msat_nonan_lon_plot_annotated_CP

msat_nonan_lon_plot <- ggplot(data = msat_nonan, aes(x = lon, y = He, succ = success, fail = failure)) + 
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2) + 
  xlab("longitude") + ylab("msat he")
msat_nonan_lon_plot_annotated <- msat_nonan_lon_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_nonan_lon_plot_annotated

###############################################################################################################

######## Pulling plots together in grid form ########

abslat_all_quad_plot_CP <- grid.arrange(mtdna_nonapi_small_abslat_plot_annotated_CP, mtdna_nonahe_abslat_plot_annotated_CP, 
                                     msat_nonan_abslat_plot_annotated_CP, ncol = 3)
abslat_all_quad_plot_CP

lat_all_quad_plot <- grid.arrange(mtdna_nonapi_small_lat_plot_annotated_CP, mtdna_nonahe_lat_plot_annotated_CP, 
                                  msat_nonan_lat_plot_annotated_CP, ncol = 3)
lat_all_quad_plot

lon_all_quad_plot <- grid.arrange(mtdna_nonapi_small_lon_plot_annotated_CP, mtdna_nonahe_lon_plot_annotated_CP, 
                                     msat_nonan_lon_plot_annotated_CP, ncol = 3)
lon_all_quad_plot

mtdna_pi_all_quad_plot <- grid.arrange(mtdna_nonapi_small_lat_plot_annotated,
                                       mtdna_nonapi_small_lon_plot_annotated, ncol = 2)
mtdna_pi_all_quad_plot

mtdna_hd_all_quad_plot <- grid.arrange(mtdna_nonahe_lat_plot_annotated, 
                                       mtdna_nonahe_lon_plot_annotated, ncol = 2)
mtdna_hd_all_quad_plot

msat_all_quad_plot <- grid.arrange(msat_nonan_lat_plot_annotated, 
                                       msat_nonan_lon_plot_annotated, ncol = 2)
msat_all_quad_plot

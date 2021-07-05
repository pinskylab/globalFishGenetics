remove(list = ls())

library(tidyverse)
library(gridExtra)

#read in data
mtdna <- read.csv("output/mtdna_assembled_2.csv", stringsAsFactors = FALSE)
msat <- read.csv("output/msatloci.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)


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

######## mtdna Pi models & plots ########

mtdna_nonapi_CP <- subset(mtdna_cp, mtdna_cp$Pi != "NA")
mtdna_nonapi_small_CP <- subset(mtdna_nonapi_CP, mtdna_nonapi_CP$Pi < 1)

#regression model
mtdna_nonapi_small_abslat_glm <- glm(mtdna_nonapi_small$Pi ~ mtdna_nonapi_small$abslat)
summary(mtdna_nonapi_small_abslat_glm)

mtdna_nonapi_small_abslat_glm_quad <- glm(mtdna_nonapi_small$Pi ~ mtdna_nonapi_small$abslat + I((mtdna_nonapi_small$abslat)^2))
summary(mtdna_nonapi_small_abslat_glm_quad)

mtdna_nonapi_small_lat_glm <- glm(mtdna_nonapi_small$Pi ~ mtdna_nonapi_small$lat)
summary(mtdna_nonapi_small_lat_glm)

mtdna_nonapi_small_lat_glm_quad <- glm(mtdna_nonapi_small$Pi ~ mtdna_nonapi_small$lat + I((mtdna_nonapi_small$lat)^2))
summary(mtdna_nonapi_small_lat_glm_quad)

mtdna_nonapi_small_abslon_glm <- glm(mtdna_nonapi_small$Pi ~ mtdna_nonapi_small$abslon)
summary(mtdna_nonapi_small_abslon_glm)

mtdna_nonapi_small_abslon_glm_quad <- glm(mtdna_nonapi_small$Pi ~ mtdna_nonapi_small$abslon + I((mtdna_nonapi_small$abslon)^2))
summary(mtdna_nonapi_small_abslon_glm_quad)

mtdna_nonapi_small_lon_glm <- glm(mtdna_nonapi_small$Pi ~ mtdna_nonapi_small$lon)
summary(mtdna_nonapi_small_lon_glm)

mtdna_nonapi_small_lon_glm_quad <- glm(mtdna_nonapi_small$Pi ~ mtdna_nonapi_small$lon + I((mtdna_nonapi_small$lon)^2))
summary(mtdna_nonapi_small_lon_glm_quad)

#plots
mtdna_nonapi_small_abslat_plot_CP <- ggplot(data = mtdna_nonapi_small_CP, aes(x = abslat, y = Pi, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("absolute latitude") + ylab("mtdna pi") + scale_color_manual(values = c("#999999", "#3333CC"))
mtdna_nonapi_small_abslat_plot_annotated_CP <- mtdna_nonapi_small_abslat_plot_CP + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonapi_small_abslat_plot_annotated_CP

mtdna_nonapi_small_abslat_plot <- ggplot(data = mtdna_nonapi_small_CP, aes(x = abslat, y = Pi)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2) + 
  xlab("absolute latitude") + ylab("mtdna pi")
mtdna_nonapi_small_abslat_plot_annotated <- mtdna_nonapi_small_abslat_plot + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonapi_small_abslat_plot_annotated

mtdna_nonapi_small_lat_plot_CP <- ggplot(data = mtdna_nonapi_small_CP, aes(x = lat, y = Pi, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("latitude") + ylab("mtdna pi") + scale_color_manual(values = c("#999999", "#3333CC"))
mtdna_nonapi_small_lat_plot_annotated_CP <- mtdna_nonapi_small_lat_plot_CP + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonapi_small_lat_plot_annotated_CP

mtdna_nonapi_small_lat_plot <- ggplot(data = mtdna_nonapi_small_CP, aes(x = lat, y = Pi)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2) + 
  xlab("latitude") + ylab("mtdna pi")
mtdna_nonapi_small_lat_plot_annotated <- mtdna_nonapi_small_lat_plot + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonapi_small_lat_plot_annotated

mtdna_nonapi_small_abslon_plot_CP <- ggplot(data = mtdna_nonapi_small_CP, aes(x = abslon, y = Pi, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("absolute longitude") + ylab("mtdna pi") + scale_color_manual(values = c("#999999", "#3333CC"))
mtdna_nonapi_small_abslon_plot_annotated <- mtdna_nonapi_small_abslon_plot_CP + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonapi_small_abslon_plot_annotated_CP

mtdna_nonapi_small_lon_plot_CP <- ggplot(data = mtdna_nonapi_small_CP, aes(x = lon, y = Pi, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("longitude") + ylab("mtdna pi") + scale_color_manual(values = c("#999999", "#3333CC"))
mtdna_nonapi_small_lon_plot_annotated_CP <- mtdna_nonapi_small_lon_plot_CP + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"),
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonapi_small_lon_plot_annotated_CP

mtdna_nonapi_small_lon_plot <- ggplot(data = mtdna_nonapi_small_CP, aes(x = lon, y = Pi)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", formula = y ~ x + I(x^2), size = 2) + 
  xlab("longitude") + ylab("mtdna pi")
mtdna_nonapi_small_lon_plot_annotated <- mtdna_nonapi_small_lon_plot + theme_bw() + ylim(0, 0.2) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"),
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonapi_small_lon_plot_annotated

######################################################################################################

######## mtdna He models & plots ########

mtdna_nonahe_CP <- subset(mtdna_cp, mtdna_cp$He != "NA")
mtdna_nonahe_CP$success <- round(mtdna_nonahe_CP$He*mtdna_nonahe_CP$n)
mtdna_nonahe_CP$failure <- round((1 - mtdna_nonahe_CP$He)*mtdna_nonahe_CP$n)

#binomial regression model
mtdna_nonahe_abslat_glm <- glm(cbind(mtdna_nonahe$success, mtdna_nonahe$failure) ~ mtdna_nonahe$abslat, family = "binomial")
summary(mtdna_nonahe_abslat_glm)

mtdna_nonahe_abslat_glm_quad <- glm(cbind(mtdna_nonahe$success, mtdna_nonahe$failure) ~ mtdna_nonahe$abslat + I((mtdna_nonahe$abslat)^2), family = "binomial")
summary(mtdna_nonahe_abslat_glm_quad)

mtdna_nonahe_lat_glm <- glm(cbind(mtdna_nonahe$success, mtdna_nonahe$failure) ~ mtdna_nonahe$lat, family = "binomial")
summary(mtdna_nonahe_lat_glm)

mtdna_nonahe_lat_glm_quad <- glm(cbind(mtdna_nonahe$success, mtdna_nonahe$failure) ~ mtdna_nonahe$lat + I((mtdna_nonahe$lat)^2), family = "binomial")
summary(mtdna_nonahe_lat_glm_quad)

mtdna_nonahe_abslon_glm <- glm(cbind(mtdna_nonahe$success, mtdna_nonahe$failure) ~ mtdna_nonahe$abslon, family = "binomial")
summary(mtdna_nonahe_abslon_glm)

mtdna_nonahe_abslon_glm_quad <- glm(cbind(mtdna_nonahe$success, mtdna_nonahe$failure) ~ mtdna_nonahe$abslon + I((mtdna_nonahe$abslon)^2), family = "binomial")
summary(mtdna_nonahe_abslon_glm_quad)

mtdna_nonahe_lon_glm <- glm(cbind(mtdna_nonahe$success, mtdna_nonahe$failure) ~ mtdna_nonahe$lon, family = "binomial")
summary(mtdna_nonahe_lon_glm)

mtdna_nonahe_lon_glm_quad <- glm(cbind(mtdna_nonahe$success, mtdna_nonahe$failure) ~ mtdna_nonahe$lon + I((mtdna_nonahe$lon)^2), family = "binomial")
summary(mtdna_nonahe_lon_glm_quad)

#plots
mtdna_nonahe_abslat_plot_CP <- ggplot(data = mtdna_nonahe_CP, aes(x = abslat, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("absolute latitude") + ylab("mtdna hd") + scale_color_manual(values = c("#999999", "#3333CC"))
mtdna_nonahe_abslat_plot_annotated_CP <- mtdna_nonahe_abslat_plot_CP + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonahe_abslat_plot_annotated_CP

mtdna_nonahe_abslat_plot <- ggplot(data = mtdna_nonahe_CP, aes(x = abslat, y = He, succ = success, fail = failure)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2) + 
  xlab("absolute latitude") + ylab("mtdna hd")
mtdna_nonahe_abslat_plot_annotated <- mtdna_nonahe_abslat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonahe_abslat_plot_annotated

mtdna_nonahe_lat_plot_CP <- ggplot(data = mtdna_nonahe_CP, aes(x = lat, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("latitude") + ylab("mtdna hd") + scale_color_manual(values = c("#999999", "#3333CC"))
mtdna_nonahe_lat_plot_annotated_CP <- mtdna_nonahe_lat_plot_CP + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonahe_lat_plot_annotated_CP

mtdna_nonahe_lat_plot <- ggplot(data = mtdna_nonahe_CP, aes(x = lat, y = He, succ = success, fail = failure)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2) + 
  xlab("latitude") + ylab("mtdna hd")
mtdna_nonahe_lat_plot_annotated <- mtdna_nonahe_lat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonahe_lat_plot_annotated

mtdna_nonahe_abslon_plot_CP <- ggplot(data = mtdna_nonahe_CP, aes(x = abslon, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 3) +  
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("absolute longitude") + ylab("mtdna hd") + scale_color_manual(values = c("#999999", "#3333CC"))
mtdna_nonahe_abslon_plot_annotated_CP <- mtdna_nonahe_abslon_plot_CP + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonahe_abslon_plot_annotated_CP

mtdna_nonahe_lon_plot_CP <- ggplot(data = mtdna_nonahe_CP, aes(x = lon, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 3) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("longitude") + ylab("mtdna hd") + scale_color_manual(values = c("#999999", "#3333CC"))
mtdna_nonahe_lon_plot_annotated_CP <- mtdna_nonahe_lon_plot_CP + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
mtdna_nonahe_lon_plot_annotated_CP

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

######## msat He models & plots ########

msat_nonan <- subset(msat_cp, msat_cp$n != "NA" & msat_cp$n != "25-32")
msat_nonan$n <- as.numeric(msat_nonan$n)
  class(msat_nonan$n) #check that is numeric
msat_nonan$success <- round(msat_nonan$He*msat_nonan$n)
msat_nonan$failure <- round((1 - msat_nonan$He)*msat_nonan$n)

#binomial regression model
msat_nonan_abslat_glm <- glm(cbind(msat_nonan$success, msat_nonan$failure) ~ msat_nonan$abslat, family = "binomial")
summary(msat_nonan_abslat_glm)

msat_nonan_abslat_glm_quad <- glm(cbind(msat_nonan$success, msat_nonan$failure) ~ msat_nonan$abslat + I((msat_nonan$abslat)^2), family = "binomial")
summary(msat_nonan_abslat_glm_quad)

msat_nonan_abslat_glm_PC <- glm(cbind(msat_nonan$success, msat_nonan$failure) ~ msat_nonan$abslat + msat_nonan$Pelagic.Coastal.f, family = "binomial")
summary(msat_nonan_abslat_glm)

msat_nonan_abslat_glm_quad_PC <- glm(cbind(msat_nonan$success, msat_nonan$failure) ~ msat_nonan$abslat + I((msat_nonan$abslat)^2) + msat_nonan$Pelagic.Coastal.f, family = "binomial")
summary(msat_nonan_abslat_glm_quad)

msat_nonan_lat_glm <- glm(cbind(msat_nonan$success, msat_nonan$failure) ~ msat_nonan$lat, family = "binomial")
summary(msat_nonan_lat_glm)

msat_nonan_lat_glm_quad <- glm(cbind(msat_nonan$success, msat_nonan$failure) ~ msat_nonan$lat + I((msat_nonan$lat)^2), family = "binomial")
summary(msat_nonan_lat_glm_quad)

msat_nonan_abslon_glm <- glm(cbind(msat_nonan$success, msat_nonan$failure) ~ msat_nonan$abslon, family = "binomial")
summary(msat_nonan_abslon_glm)

msat_nonan_abslon_glm_quad <- glm(cbind(msat_nonan$success, msat_nonan$failure) ~ msat_nonan$abslon + I((msat_nonan$abslon)^2), family = "binomial")
summary(msat_nonan_abslon_glm_quad)

msat_nonan_lon_glm <- glm(cbind(msat_nonan$success, msat_nonan$failure) ~ msat_nonan$lon, family = "binomial")
summary(msat_nonan_lon_glm)

msat_nonan_lon_glm_quad <- glm(cbind(msat_nonan$success, msat_nonan$failure) ~ msat_nonan$lon + I((msat_nonan$lon)^2), family = "binomial")
summary(msat_nonan_lon_glm_quad)

#plots
msat_nonan_abslat_plot_CP <- ggplot(data = msat_nonan, aes(x = abslat, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("absolute latitude") + ylab("msat he") + scale_color_manual(values = c("#999999", "#3333CC"))
msat_nonan_abslat_plot_annotated_CP <- msat_nonan_abslat_plot_CP + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_nonan_abslat_plot_annotated_CP

msat_nonan_abslat_plot <- ggplot(data = msat_nonan, aes(x = abslat, y = He, succ = success, fail = failure)) + 
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2) + 
  xlab("absolute latitude") + ylab("msat he")
msat_nonan_abslat_plot_annotated <- msat_nonan_abslat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_nonan_abslat_plot_annotated

msat_nonan_lat_plot_CP <- ggplot(data = msat_nonan, aes(x = lat, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("latitude") + ylab("msat he") + scale_color_manual(values = c("#999999", "#3333CC"))
msat_nonan_lat_plot_annotated_CP <- msat_nonan_lat_plot_CP + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_nonan_lat_plot_annotated_CP

msat_nonan_lat_plot <- ggplot(data = msat_nonan, aes(x = lat, y = He, succ = success, fail = failure)) + 
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2) + 
  xlab("latitude") + ylab("msat he")
msat_nonan_lat_plot_annotated <- msat_nonan_lat_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_nonan_lat_plot_annotated

msat_nonan_abslon_plot_CP <- ggplot(data = msat_nonan, aes(x = abslon, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("absolute longitude") + ylab("msat he") + scale_color_manual(values = c("#999999", "#3333CC"))
msat_nonan_abslon_plot_annotated_CP <- msat_nonan_abslon_plot_CP + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_nonan_abslon_plot_annotated_CP

msat_nonan_lon_plot_CP <- ggplot(data = msat_nonan, aes(x = lon, y = He, succ = success, fail = failure, color = Pelagic_Coastal)) + 
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2, aes(linetype = Pelagic_Coastal)) + 
  xlab("longitude") + ylab("msat he") + scale_color_manual(values = c("#999999", "#3333CC"))
msat_nonan_lon_plot_annotated_CP <- msat_nonan_lon_plot_CP + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_nonan_lon_plot_annotated_CP

msat_nonan_lon_plot <- ggplot(data = msat_nonan, aes(x = lon, y = He, succ = success, fail = failure)) + 
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), formula = cbind(succ, fail) ~ x + I(x^2), size = 2) + 
  xlab("longitude") + ylab("msat he")
msat_nonan_lon_plot_annotated <- msat_nonan_lon_plot + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 20), legend.title = element_blank())
msat_nonan_lon_plot_annotated

abslat_all_quad_plot_CP <- grid.arrange(mtdna_nonapi_small_abslat_plot_annotated_CP, mtdna_nonahe_abslat_plot_annotated_CP, 
                                     msat_nonan_abslat_plot_annotated_CP, ncol = 3)
abslat_all_quad_plot_CP

lat_all_quad_plot <- grid.arrange(mtdna_nonapi_small_lat_plot_annotated_CP, mtdna_nonahe_lat_plot_annotated_CP, 
                                  msat_nonan_lat_plot_annotated_CP, ncol = 3)
lat_all_quad_plot

abslon_all_quad_plot <- grid.arrange(mtdna_nonapi_small_abslon_plot_annotated, mtdna_nonahe_abslon_plot_annotated_CP, 
                                     msat_nonan_abslon_plot_annotated_CP, ncol = 3)
abslon_all_quad_plot

lon_all_quad_plot <- grid.arrange(mtdna_nonapi_small_lon_plot_annotated_CP, mtdna_nonahe_lon_plot_annotated_CP, 
                                     msat_nonan_lon_plot_annotated_CP, ncol = 3)
lon_all_quad_plot

mtdna_pi_all_quad_plot <- grid.arrange(mtdna_nonapi_small_lat_plot_annotated, mtdna_nonapi_small_abslat_plot_annotated,
                                       mtdna_nonapi_small_lon_plot_annotated, ncol = 3)
mtdna_pi_all_quad_plot

mtdna_hd_all_quad_plot <- grid.arrange(mtdna_nonahe_lat_plot_annotated, mtdna_nonahe_abslat_plot_annotated,
                                       mtdna_nonahe_lon_plot_annotated, ncol = 3)
mtdna_hd_all_quad_plot

msat_all_quad_plot <- grid.arrange(msat_nonan_lat_plot_annotated, msat_nonan_abslat_plot_annotated,
                                       msat_nonan_lon_plot_annotated, ncol = 3)
msat_all_quad_plot
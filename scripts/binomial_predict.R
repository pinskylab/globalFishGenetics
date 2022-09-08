#Predict w/GLMMs and plot with binned means

remove(list = ls())

library(tidyverse)
library(gridExtra)
library(plotrix)
library(lme4)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(effects)
library(sjPlot)
library(data.table)
library(DescTools)

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)
msat <- read.csv("output/msatloci_assembled.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)
mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                  'Centroid')], all.x = TRUE)
msat <- merge(msat, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                'Centroid')], all.x = TRUE)

#clean up dataframes
#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)

#calculate abslat
mtdna_small$abslat <- abs(mtdna_small$lat)
msat$abslat <- abs(msat$lat)

########## pi
#subset mtdna to remove Pi = NA columns
mtdna_small_pi <- subset(mtdna_small, mtdna_small$Pi != "NA")

#log transform and check again
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$Pi != 0) #if any zeros will screw up log transformation (log10(0) is undefined, also probably shouldn't be there anyway)
mtdna_small_pi$logpi <- log10(mtdna_small_pi$Pi)

#remove logpi = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logpi != "Inf")

######## building lat, abslat & lon model ########

#### clean up data ####
#scale geographic variables
mtdna_small_pi$lat_scale <- scale(mtdna_small_pi$lat)
mtdna_small_pi$abslat_scale <- scale(mtdna_small_pi$abslat)

#convert lon to radians
mtdna_small_pi$lon_360 <- mtdna_small_pi$lon + 180 #convert (-180,180) to (0,360)
mtdna_small_pi$lon_rad <- (2*pi*mtdna_small_pi$lon_360)/360

lat_model_pi <- lmer(logpi ~ lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                       (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

abslat_model_pi <- lmer(logpi ~ abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                       (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

lon_model_pi <- lmer(logpi ~ sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                       (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#predict for abslat
abslat <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
#need to scale this using the same scale mean and sd as in real dataset
#scale from dataset (mtdna pi) -- mean = 27.24903, SD = 15.2667 --> got from dataset with attributes(mtdna_small_pi$abslat_scale)

abslat_scale <- (abslat - 27.24903)/15.2667

pi_abslat_predict_df <- as.data.frame(abslat_scale)
  pi_abslat_predict_df$Family <- c(rep(0, times = 9))
  pi_abslat_predict_df$Genus <- c(rep(0, times = 9))
  pi_abslat_predict_df$spp <- c(rep(0, times = 9))
  pi_abslat_predict_df$Source <- c(rep(0, times = 9))
  pi_abslat_predict_df$MarkerName <- c(rep(0, times = 9))
  pi_abslat_predict_df$Site <- c(rep(0, times = 9))

pi_predict_abslat <- predict(abslat_model_pi, pi_abslat_predict_df, re.form = NA)
  pi_abslat_predict_df$predict_pi <- 10^(pi_predict_abslat)


#other variables -- just use average?
#not interested in effect of range on regression necessarily, just account for it?

#bootstrap for confidence intervals
#want to use bootMer but on predicted dataset
  
pi_abslat_boot <- bootMer(abslat_model_pi, 
                          FUN=function(x)predict(x, pi_abslat_predict_df, re.form = NA), nsim = 100)
pi_abslat_predict_df$lower_ci <- 10^(apply(pi_abslat_boot$t, 2, quantile, 0.025))
pi_abslat_predict_df$upper_ci <- 10^(apply(pi_abslat_boot$t, 2, quantile, 0.975))


#add real means

#mtdna pi lat 10 binning

mtdna_small_pi$lat_round <- NA #doesn't work quite right with negatives (BUT can specify bins as want)
mtdna_small_pi$abslat_round <- NA

for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$lat_round[i] <- DescTools::RoundTo(mtdna_small_pi$lat[i], multiple = 10, FUN = floor)}
}

for(i in 1:nrow(mtdna_small_pi)) {
  cat(paste(i, " ", sep = ''))
  {mtdna_small_pi$abslat_round[i] <- DescTools::RoundTo(mtdna_small_pi$abslat[i], multiple = 10, FUN = floor)}
}

mtdna_small_pi$abslat_round <- as.factor(mtdna_small_pi$abslat_round)
mtdna_small_pi$lat_round <- as.factor(mtdna_small_pi$lat_round)

#subset by lat bands and then get mean and plot that

#abslat bands
mtdna_small_pi <- data.table(mtdna_small_pi)
abslat_logpi_mean <- mtdna_small_pi[, mean(logpi), by = list(abslat_round)]
  colnames(abslat_logpi_mean) <- c("abslat", "logpi_mean")
  abslat_logpi_mean$pi_mean <- 10^(abslat_logpi_mean$logpi_mean)
abslat_pi_SE <- mtdna_small_pi[, std.error(Pi), by = list(abslat_round)] #not quite sure if this is correct...
  colnames(abslat_pi_SE) <- c("abslat", "pi_SE")

abslat_pi_binned_means <- merge(abslat_logpi_mean, abslat_pi_SE, by = "abslat")
  abslat_pi_binned_means <- abslat_pi_binned_means[order(abslat), ]
  abslat_pi_binned_means$abslat <- as.numeric(as.character(abslat_pi_binned_means$abslat))
  abslat_pi_binned_means <- abslat_pi_binned_means %>% add_row(abslat = 80)

#merge predicted and binned means df together
pi_abslat_full <- cbind(pi_abslat_predict_df, abslat_pi_binned_means)
pi_abslat_full$mean_lowerSE <- pi_abslat_full$pi_mean - pi_abslat_full$pi_SE
pi_abslat_full$mean_upperSE <- pi_abslat_full$pi_mean + pi_abslat_full$pi_SE

#for legend
colors <- c("Raw data" = "#87A4CA", "10-degree binned means" = "#3F6DAA", "Regression" = "black")

msat_abslat_plot_both <- ggplot(data = pi_abslat_full) + 
  geom_point(data = mtdna_small_pi, aes(x = abslat, y = Pi, color = "Raw data"), 
             size = 4, alpha = 0.25) + 
  geom_line(aes(x = abslat, y = predict_pi, color = "Regression"), size = 2) + 
  #geom_point(aes(x = abslat, y = predict_pi), size = 4, color = "black") + 
  geom_ribbon(aes(x = abslat, ymin = lower_ci, ymax = upper_ci), alpha = 0.1) + 
  geom_point(aes(x = abslat, y = pi_mean, color = "10-degree binned means"), size = 8, shape = "square") + 
  geom_errorbar(aes(x = abslat, ymin = mean_lowerSE, ymax = mean_upperSE, color = "10-degree binned means"), 
                width = 0.75, size = 1.5) + 
 # geom_smooth(aes(x = abslat, y = Pi), method = "lm", formula = y ~ splines::bs(x), size = 2) + 
 # geom_smooth(data = mtdna_pi_small, aes(x = abslat, y = Pi), 
#              method = "lm", formula = y ~ x + I(x^2), color = "black", size = 2) +
   xlab("absolute latitude") + ylab("mtdna pi") + labs(color = "Legend") + 
  scale_color_manual(values = colors)
msat_abslat_plot_annotated_both <- msat_abslat_plot_both + theme_bw() + 
  scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) + 
  scale_y_continuous(limits = c(0, 0.05)) + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 26, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 20, color = "black"), 
        legend.position = "top", legend.text = element_text(size = 18), legend.title = element_blank())
msat_abslat_plot_annotated_both

#CI --> Zoe did it this way
#used Ben Bolker's tutorial https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#predictions-andor-confidence-or-prediction-intervals-on-predictions

mm <- model.matrix(terms(abslat_model_pi), pi_predict_abslat)

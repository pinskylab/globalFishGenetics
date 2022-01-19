#Script for building models

remove(list = ls())

library(tidyverse)
library(gridExtra)
library(plotrix)
library(lme4)
library(MuMIn)
library(glmmTMB)
library(DHARMa)

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

##############################################################################################################

######## building models for mtdna Hd ########

######## check effect of independent variables ######## 
  
#subset mtdna bc are a few outliers with very high bp -- entire mtdna?
mtdna_small <- subset(mtdna, as.numeric(mtdna$bp) < 2000)

#subset mtdna to remove He = NA columns
mtdna_small_He <- subset(mtdna_small, mtdna_small$He != "NA")

#### checking effect of bp ####
mtdna_small_He$bp <- as.numeric(mtdna_small_He$bp)
He_bp <- ggplot(mtdna_small_He, aes(x = as.numeric(bp), y = He)) + 
  geom_point() + 
  geom_smooth(method = "lm") 
He_bp #slight positive relationship --> more bp = more possible haplotypes = higher div

##### checking effect of Coastal v. Pelagic ####

He_CP <- ggplot(mtdna_small_He, aes(x = as.numeric(bp), y = He, colour = Pelagic_Coastal)) + 
  geom_point() + 
  geom_smooth(method = "lm")
He_CP #looks like coastal have higher div than pelagic but is significant?

He_CP_hist <- ggplot(mtdna_small_He, aes(x = He, colour = Pelagic_Coastal)) + 
  geom_histogram(binwidth = 0.1, fill = "white", alpha = 0.5, position = "identity")
He_CP_hist #overall distribution looks similar, just much higher counts in coastal

#### checking effect of position in range ####
#fix character type
mtdna_small_He$Centroid <- as.numeric(mtdna_small_He$Centroid)
mtdna_small_He$Half_RangeSize <- as.numeric(mtdna_small_He$Half_RangeSize)

mtdna_small_He$range_position <- NA #create column to fill in

for (i in 1:nrow(mtdna_small_He)) { #calculate distance from range center as percentage (0-1 both sides of centroid, 1 = all the way at range edge)
  mtdna_small_He$range_position[i] <- abs((mtdna_small_He$lat[i] - mtdna_small_He$Centroid[i])/mtdna_small_He$Half_RangeSize[i])
}

#check those >1 and round to 1 (slight discrepancy btwn aquamaps and reality)
mtdna_check <- subset(mtdna_small_He, mtdna_small_He$range_position > 1) #often right at aquamaps limit, round to 1 and keep
mtdna_small_He$range_position[mtdna_small_He$range_position > 1] <- 1

He_range <- ggplot(mtdna_small_He, aes(x = range_position, y = He)) + 
  geom_point() + 
  geom_smooth(method = "lm")
He_range #doesn't really look like a relationship but also missing some close to range edge

#### checking effect of MarkerName ####
He_MN <- boxplot(He ~ MarkerName, data = mtdna_small_He)

MN_split_plot <- ggplot(aes(x = as.numeric(bp), y = He), data = mtdna_small_He) + 
    geom_point() + 
    facet_wrap(~ MarkerName) + 
    xlab("bp") + ylab("He")
MN_split_plot #doesn't look to be much of an effect of MarkerName on He -- diff regions of mtDNA show similar spans of diversity

MN_split_hist_plot <- ggplot(aes(x = He), data = mtdna_small_He) + 
     geom_histogram(binwidth = 0.1) + 
     facet_wrap(~ MarkerName) + 
     xlab("He")
MN_split_hist_plot #looks like more of effect but again, most MarkerName categories have very few observations

#### checking effect of spp ####
He_spp <- boxplot(He ~ spp, data = mtdna_small_He) #diff species have diff ranges of He (keep)

spp_color_plot <- ggplot(data = mtdna_small_He, aes(x = as.numeric(bp), y = He, colour = spp)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position = "none")
spp_color_plot

##### checking effect of site ####
He_site <- boxplot(He ~ Site, data = mtdna_small_He) #many sites BUT often only one observation per site so not worth keeping

site_color_plot <- ggplot(data = mtdna_small_He, aes(x = as.numeric(bp), y = He, colour = Site)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position = "none")
site_color_plot

######## build null model (only includes nuisance variables) ########
#variables to include: bp, position in spp range, (1|MarkerName), (1|Family/Genus/spp), (1|Source), (1|Site)
#dredge to see which ones matter most?

#### clean up data ####
#have to scale numerical data bc otherwise throws convergence issues
#issue with scaling is makes interpretation difficult
mtdna_small_He$bp_scale <- scale(as.numeric(mtdna_small_He$bp))

#subset to only those with range_position
mtdna_small_He <- subset(mtdna_small_He, range_position != "NA")

#transform He data to deal with 1s (no 0s in dataset as excluded monomorphic data)
mtdna_small_He$transformed_He <- NA #create column to fill in

for (i in 1:nrow(mtdna_small_He)) { #transform data to handle 1s (Douma & Weedon (2018) Methods in Ecology & Evolution)
  mtdna_small_He$transformed_He[i] <- ((mtdna_small_He$He[i]*(mtdna_small_He$n[i] - 1)) + 0.5)/mtdna_small_He$n[i]
}

#### run null model ####
#null model (no abslat/lat/lon)
beta_null_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + (1|Family/Genus/spp) + (1|Source) + 
                                  (1|Site) + (1|MarkerName), family = beta_family, data = mtdna_small_He)

beta_null_model_test <- glmmTMB(transformed_He ~ bp_scale + range_position, family = beta_family, data = mtdna_small_He)
beta_sim <- simulateResiduals(fittedModel = beta_null_model_full, n = 1000, plot = F)
plotQQunif(beta_sim)
plotResiduals(beta_sim)
countZeroes <- function(x) sum(x == 0.5)
testGeneric(sim_test, summary = countZeroes, alternative = "less")

#test binomial model again
mtdna_small_He$success <- round(mtdna_small_He$He*mtdna_small_He$n)
mtdna_small_He$failure<- round((1 - mtdna_small_He$He)*mtdna_small_He$n)

binomial_null <- glmer(cbind(success, failure) ~ bp_scale + range_position + (1|Family/Genus/spp) + (1|Source) + 
                         (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail")

binomial_null_test <- glmer(cbind(success, failure) ~ range_position + (1|Family/Genus/spp) + (1|Source) + 
                         (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail")


binomial_null_sim <- simulateResiduals(fittedModel = binomial_null, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_null_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_null_sim) #residuals against predicted value -- looking for uniformity

plotResiduals(binomial_null_sim, mtdna_small_He$bp)
plotResiduals(binomial_null_sim, mtdna_small_He$range_position)

#compare binomial w/ and w/out bp_scale
anova(binomial_null, binomial_null_test)

#add random effect for every unit
mtdna_small_He$ID <- (1:1680)
binomial_null_IDRE <- glmer(cbind(success, failure) ~ bp_scale + range_position + (1|Family/Genus/spp) + (1|Source) + 
                                   (1|Site) + (1|ID) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail")

binomial_null_IDRE_sim <- simulateResiduals(fittedModel = binomial_null_IDRE, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_null_IDRE_sim)
plotResiduals(binomial_null_IDRE_sim) #doesn't make things better -- makes them worse as far as fit w/residuals looks


#checking fit with DHARMa
null_model_mtdna_he_sim_output <- simulateResiduals(fittedModel = beta_null_model_full, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(null_model_mtdna_he_sim_output) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(null_model_mtdna_he_sim_output) #residuals against predicted value -- looking for uniformity

#testing against specific predictors
plotResiduals(null_model_mtdna_he_sim_output, mtdna_small_He$bp_scale)
plotResiduals(null_model_mtdna_he_sim_output, mtdna_small_He$range_position)

#test dispersion for binomial model
testDispersion(binomial_null_sim)

######## build lat, abslat & lon model ########
#have abslat, lat, lon or some combo of three

#### clean up data ####
#scale geographic variables
mtdna_small_He$lat_scale <- scale(mtdna_small_He$lat)
#mtdna_small_He$lon_scale <- scale(mtdna_small_He$lon)
mtdna_small_He$abslat <- abs(mtdna_small_He$lat)
mtdna_small_He$abslat_scale <- scale(mtdna_small_He$abslat)

#convert lon to radians
mtdna_small_He$lon_360 <- mtdna_small_He$lon + 180 #convert (-180,180) to (0,360)
mtdna_small_He$lon_rad <- (2*pi*mtdna_small_He$lon_360)/360

#check circular relationship
y <- sin(mtdna_small_He$lon_rad) + cos(mtdna_small_He$lon_rad)
plot(mtdna_small_He$lon_rad, y)

#lat model
beta_lat_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                  (1|Site), family = beta_family, data = mtdna_small_He) #beta_family (betabinomial is binomial)
#dredge_lat <- dredge(lat_model_full)

mtdna_he_lat_sim_output <- simulateResiduals(fittedModel = beta_lat_model_full, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(mtdna_he_lat_sim_output) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(mtdna_he_lat_sim_output) #residuals against predicted value -- looking for uniformity

#testing against specific predictors
plotResiduals(mtdna_he_lat_sim_output, mtdna_small_He$lat_scale)

binomial_lat <- glmer(cbind(success, failure) ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                         (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail")

binomial_lat_sim <- simulateResiduals(fittedModel = binomial_lat, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_lat_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_lat_sim) #residuals against predicted value -- looking for uniformity

plotResiduals(binomial_lat_sim, mtdna_small_He$lat_scale)

#lon model
beta_lon_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                     (1|Site), family = beta_family, data = mtdna_small_He) #beta_family (betabinomial is binomial)
#dredge_lon <- dredge(lon_model_full)

beta_lon_sim <- simulateResiduals(fittedModel = beta_lon_model_full, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(beta_lon_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(beta_lon_sim) #residuals against predicted value -- looking for uniformity

binomial_lon <- glmer(cbind(success, failure) ~ bp_scale + range_position + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                        (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail")

binomial_lon_sim <- simulateResiduals(fittedModel = binomial_lon, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_lon_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_lon_sim) #residuals against predicted value -- looking for uniformity

plotResiduals(binomial_lon_sim, mtdna_small_He$lon_scale)

#abslat model
beta_abslat_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), family = beta_family, data = mtdna_small_He) #beta_family (betabinomial is binomial)
#dredge_abslat <- dredge(abslat_model_full)

binomial_abslat <- glmer(cbind(success, failure) ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                        (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail")

binomial_abslat_sim <- simulateResiduals(fittedModel = binomial_abslat, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_abslat_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_abslat_sim) #residuals against predicted value -- looking for uniformity

plotResiduals(binomial_abslat_sim, mtdna_small_He$abslat_scale)

#lat & lon model
beta_lat_lon_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), family = beta_family, data = mtdna_small_He) #beta_family (betabinomial is binomial)
#dredge_lat_lon <- dredge(lat_lon_model_full)

binomial_lat_lon <- glmer(cbind(success, failure) ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                           (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

binomial_lat_lon_sim <- simulateResiduals(fittedModel = binomial_lat_lon, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_lat_lon_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_lat_lon_sim) #residuals against predicted value -- looking for uniformity
testDispersion(binomial_lat_lon_sim)

#abslat & lon model
beta_abslat_lon_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                     (1|Site), family = beta_family, data = mtdna_small_He) #beta_family (betabinomial is binomial)
#dredge_abslat_lon <- dredge(abslat_lon_model_full)

binomial_abslat_lon <- glmer(cbind(success, failure) ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                            (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

binomial_abslat_lon_sim <- simulateResiduals(fittedModel = binomial_abslat_lon, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_abslat_lon_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_abslat_lon_sim) #residuals against predicted value -- looking for uniformity

#compare models with anova
anova(binomial_null, binomial_abslat, binomial_lat, binomial_lon, binomial_abslat_lon, binomial_lat_lon)

#visualize effects
effects_range_position <- effect(term = "range_position", mod = binomial_lat_lon)
summary(effects_range_position)
x_range_position <- as.data.frame(effects_range_position)

#range_position plot
range_position_plot <- ggplot() + 
  geom_point(data = mtdna_small_He, aes(x = range_position, y = He)) + 
  geom_point(data = x_range_position, aes(x = range_position, y = fit), color = "blue") + 
  geom_line(data = x_range_position, aes(x = range_position, y = fit), color = "blue") + 
  geom_ribbon(data = x_range_position, aes(x = range_position, ymin = lower, ymax = upper), alpha = 0.3, fill = "blue") + 
  labs(x = "range_position (% from centroid)", y = "He")

effects_bp_scale <- effect(term = "bp_scale", mod = binomial_lat_lon)
summary(effects_bp_scale)
x_bp_scale <- as.data.frame(effects_bp_scale)

#bp_scale plot
bp_scale_plot <- ggplot() + 
  geom_point(data = mtdna_small_He, aes(x = bp_scale, y = He)) + 
  geom_point(data = x_bp_scale, aes(x = bp_scale, y = fit), color = "blue") + 
  geom_line(data = x_bp_scale, aes(x = bp_scale, y = fit), color = "blue") + 
  geom_ribbon(data = x_bp_scale, aes(x = bp_scale, ymin = lower, ymax = upper), alpha = 0.3, fill = "blue") + 
  labs(x = "bp scale", y = "He")

effects_abslat <- effect(term = "abslat_scale", mod = binomial_abslat_lon)
summary(effects_abslat)
x_abslat <- as.data.frame(effects_abslat)

#abslat plot
abslat_plot <- ggplot() + 
  geom_point(data = mtdna_small_He, aes(x = abslat_scale, y = He)) + 
  geom_point(data = x_abslat, aes(x = abslat_scale, y = fit), color = "blue") + 
  geom_line(data = x_abslat, aes(x = abslat_scale, y = fit), color = "blue") + 
  geom_ribbon(data = x_abslat, aes(x = abslat_scale, ymin = lower, ymax = upper), alpha = 0.3, fill = "blue") + 
  labs(x = "abslat scale", y = "He")

#bootMer

##############################################################################################################

######## building models for mtdna pi ########

######## check effect of independent variables ######## 

#subset mtdna to remove Pi = NA columns
mtdna_small_pi <- subset(mtdna_small, mtdna_small$Pi != "NA")

#check pi distribution -- is it normal?
pi_hist <- ggplot(mtdna_small_pi, aes(x = Pi)) + 
  geom_histogram(binwidth = 0.01)
pi_hist #skewed to the right (few high values, many low values)

#log transform and check again
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$Pi != 0) #if any zeros will screw up log transformation (log10(0) is undefined, also probably shouldn't be there anyway)
mtdna_small_pi$logpi <- log10(mtdna_small_pi$Pi)

logpi_hist <- ggplot(mtdna_small_pi, aes(x = logpi)) + 
  geom_histogram(binwidth = 0.01)
logpi_hist #much better

#remove logpi = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$logpi != "Inf")

#### checking effect of bp ####
Pi_bp <- ggplot(mtdna_small_pi, aes(x = as.numeric(bp), y = Pi)) + 
  geom_point() + 
  geom_smooth(method = "lm")
Pi_bp #slight negative relationship --> more bp = likely to be more similar bp (fewer average pairwise differences) = lower div

##### checking effect of Coastal v. Pelagic ####
Pi_CP <- ggplot(mtdna_small_pi, aes(x = as.numeric(bp), y = Pi, colour = Pelagic_Coastal)) + 
  geom_point() + 
  geom_smooth(method = "lm")
Pi_CP #looks like pelagic have higher div than coastal, and diff relationships as well? --> opposite of mtdna Hd??

Pi_CP_hist <- ggplot(mtdna_small_pi, aes(x = Pi, colour = Pelagic_Coastal)) + 
  geom_histogram(binwidth = 0.01, fill = "white", alpha = 0.5, position = "identity")
Pi_CP_hist #v different distributions -- coastal is lower and pelagic higher BUT pelagic much smaller sample size

#### checking effect of position in range ####
#fix character type
mtdna_small_pi$Centroid <- as.numeric(mtdna_small_pi$Centroid)
mtdna_small_pi$Half_RangeSize <- as.numeric(mtdna_small_pi$Half_RangeSize)

mtdna_small_pi$range_position <- NA #create column to fill in

for (i in 1:nrow(mtdna_small_pi)) { #gcalculate distance from range center as percentage (0-1 both sides of centroid, 1 = all the way at range edge)
  mtdna_small_pi$range_position[i] <- abs((mtdna_small_pi$lat[i] - mtdna_small_pi$Centroid[i])/mtdna_small_pi$Half_RangeSize[i])
}

#check those >1 and round to 1 (slight discrepancy btwn aquamaps and reality)
mtdna_check <- subset(mtdna_small_pi, mtdna_small_pi$range_position > 1) #often right at aquamaps limit, round to 1 and keep
mtdna_small_pi$range_position[mtdna_small_pi$range_position > 1] <- 1

Pi_range <- ggplot(mtdna_small_pi, aes(x = range_position, y = Pi)) + 
  geom_point() + 
  geom_smooth(method = "lm")
Pi_range #looks like slight negative relationship (lower div farther from range edge)

#### checking effect of MarkerName ####
Pi_MN <- boxplot(Pi ~ MarkerName, data = mtdna_small_pi)

MN_split_plot_pi <- ggplot(aes(x = as.numeric(bp), y = Pi), data = mtdna_small_pi) + 
  geom_point() + 
  facet_wrap(~ MarkerName) + 
  xlab("bp") + ylab("Pi")
MN_split_plot_pi #very different sample sizes

MN_split_hist_plot_pi <- ggplot(aes(x = Pi), data = mtdna_small_pi) + 
  geom_histogram(binwidth = 0.01) + 
  facet_wrap(~ MarkerName) + 
  xlab("Pi")
MN_split_hist_plot_pi #again, most MarkerName categories have very few observations

#### checking effect of spp ####
Pi_spp <- boxplot(Pi ~ spp, data = mtdna_small_pi) #diff species have diff ranges of pi

spp_color_plot_pi <- ggplot(data = mtdna_small_pi, aes(x = as.numeric(bp), y = Pi, colour = spp)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position = "none")
spp_color_plot_pi

######## build null model (only includes nuisance variables) ########
#variables to include: bp, position in spp range, (1|MarkerName), (1|Family/Genus/spp), (1|Source), (1|Site)
#dredge to see which ones matter most?

#have to scale numerical data bc otherwise throws convergence issues
#issue with scaling is makes interpretation difficult
mtdna_small_pi$bp_scale <- scale(as.numeric(mtdna_small_pi$bp))

#subset to only those with range_position
mtdna_small_pi <- subset(mtdna_small_pi, range_position != "NA")

#### run null model ####
#null model (no abslat/lat/lon)
null_model_full_pi <- lmer(logpi ~ range_position + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                             (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?
#dredge_null_pi <- dredge(null_model_full_pi) #not getting same scale/convergence issues as with He
#dredge_null_pi #same top models as with He

null_model_test <- lmer(logpi ~ range_position + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                             (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#checking fit with DHARMa
null_model_pi_sim_output <- simulateResiduals(null_model_test, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(null_model_pi_sim_output) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(null_model_pi_sim_output) #residuals against predicted value -- looking for uniformity
testDispersion(null_model_pi_sim_output)

#testing against specific predictors
plotResiduals(null_model_pi_sim_output, mtdna_small_pi$bp_scale)
plotResiduals(null_model_pi_sim_output, mtdna_small_pi$range_position)

#pull p-values
coefs <- data.frame(coef(summary(abslat_lon_model_full_pi)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

######## building lat, abslat & lon model ########
#have abslat, lat, lon or some combo of three

#### clean up data ####
#scale geographic variables
mtdna_small_pi$lat_scale <- scale(mtdna_small_pi$lat)
#mtdna_small_pi$lon_scale <- scale(mtdna_small_pi$lon)
mtdna_small_pi$abslat <- abs(mtdna_small_pi$lat)
mtdna_small_pi$abslat_scale <- scale(mtdna_small_pi$abslat)

#convert lon to radians
mtdna_small_pi$lon_360 <- mtdna_small_pi$lon + 180 #convert (-180,180) to (0,360)
mtdna_small_pi$lon_rad <- (2*pi*mtdna_small_pi$lon_360)/360

#lat model
lat_model_full_pi <- lmer(logpi ~ range_position + lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                          (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
#dredge_lat_pi <- dredge(lat_model_full_pi)

lat_model_pi_sim_output <- simulateResiduals(lat_model_full_pi, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(lat_model_pi_sim_output) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(lat_model_pi_sim_output) #residuals against predicted value -- looking for uniformity
#testDispersion(lat_model_pi_sim_output)

#abslat model
abslat_model_full_pi <- lmer(logpi ~ range_position + abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                            (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
#dredge_abslat_pi <- dredge(abslat_model_full_pi)

abslat_model_pi_sim_output <- simulateResiduals(abslat_model_full_pi, plot = F)
plotQQunif(abslat_model_pi_sim_output) #QQplot
plotResiduals(abslat_model_pi_sim_output) #residuals

#lon model
lon_model_full_pi <- lmer(logpi ~ range_position + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                            (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
#dredge_lon_pi <- dredge(lon_model_full_pi)

lon_model_pi_sim_output <- simulateResiduals(lon_model_full_pi, plot = F)
plotQQunif(lon_model_pi_sim_output) #QQplot
plotResiduals(lon_model_pi_sim_output) #residuals

#lat & lon model
lat_lon_model_full_pi <- lmer(logpi ~ range_position + lat_scale + I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                            (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
#dredge_lat_lon_pi <- dredge(lat_lon_model_full_pi)

lat_lon_model_pi_sim_output <- simulateResiduals(lat_lon_model_full_pi, plot = F)
plotQQunif(lat_lon_model_pi_sim_output) #QQplot
plotResiduals(lat_lon_model_pi_sim_output) #residuals

#abslat model
abslat_lon_model_full_pi <- lmer(logpi ~ range_position + abslat_scale + I(abslat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
#dredge_abslat_lon_pi <- dredge(abslat_lon_model_full_pi)

abslat_lon_model_pi_sim_output <- simulateResiduals(abslat_lon_model_full_pi, plot = F)
plotQQunif(abslat_lon_model_pi_sim_output) #QQplot
plotResiduals(abslat_lon_model_pi_sim_output) #residuals
testDispersion(abslat_lon_model_pi_sim_output)

#compare models with anova
anova(null_model_full_pi, lat_model_full_pi, lon_model_full_pi, abslat_model_full_pi, lat_lon_model_full_pi, abslat_lon_model_full_pi)

##############################################################################################################

######## building models for msat He########

######## check effect of independent variables ######## 

#remove missing n values
msat$n <- as.numeric(msat$n)
msat <- subset(msat, msat$n != "NA")

#### checking effect of CrossSpp ####
msat_crossspp <- ggplot(msat, aes(x = CrossSpp, y = He)) + 
  geom_point() + 
  geom_smooth(method = "lm")
msat_crossspp #slight negative relationship -- slightly lower He with CrossSpp markers (makes sense)

##### checking effect of Repeat ####
msat$Repeat <- as.numeric(msat$Repeat)
msat <- subset(msat, msat$Repeat != "NA")

msat_repeat <- ggplot(msat, aes(x = Repeat, y = He)) + 
  geom_point() + 
  geom_smooth(method = "lm")
msat_repeat #positive relationship, though don't know what is going on with the repeats <1?

##### checking effect of Coastal v. Pelagic ####
msat_CP_hist <- ggplot(msat, aes(x = He, colour = Pelagic_Coastal)) + 
  geom_histogram(binwidth = 0.1, fill = "white", alpha = 0.5, position = "identity")
msat_CP_hist #overall distribution looks similar, just much higher counts in Coastal

#### checking effect of position in range ####
#fix character type
msat$Centroid <- as.numeric(msat$Centroid)
msat$Half_RangeSize <- as.numeric(msat$Half_RangeSize)

msat$range_position <- NA #create column to fill in

for (i in 1:nrow(msat)) { #gcalculate distance from range center as percentage (0-1 both sides of centroid, 1 = all the way at range edge)
  msat$range_position[i] <- abs((msat$lat[i] - msat$Centroid[i])/msat$Half_RangeSize[i])
}

#check those >1 and round to 1 (slight discrepancy btwn aquamaps and reality)
msat_check <- subset(msat, msat$range_position > 1) #often right at aquamaps limit, round to 1 and keep
msat_check$latcomp <- msat_check$lat
msat$range_position[msat$range_position > 1] <- 1

msat_range <- ggplot(msat, aes(x = range_position, y = He)) + 
  geom_point() + 
  geom_smooth(method = "lm")
msat_range  #looks like slight positive relationship?

#### checking effect of Primer Note ####
msat_PN <- boxplot(He ~ PrimerNote, data = msat) #distributions look very similar

PN_split_plot <- ggplot(aes(x = Repeat, y = He), data = msat) + 
  geom_point() + 
  facet_wrap(~ PrimerNote) + 
  xlab("Repeat") + ylab("He")
PN_split_plot #doesn't look to be much of an effect of PN on He

PN_split_hist_plot <- ggplot(aes(x = He), data = msat) + 
  geom_histogram(binwidth = 0.1) + 
  facet_wrap(~ PrimerNote) + 
  xlab("He")
PN_split_hist_plot #again, distributions look to be the same (or very similar), just far more non-PrimerNote observations

#### checking effect of spp ####
msat_spp <- boxplot(He ~ spp, data = msat) #diff species have diff ranges of He (keep)

msat_spp_color_plot <- ggplot(data = msat, aes(x = Repeat, y = He, colour = spp)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position = "none")
msat_spp_color_plot

##### checking effect of site ####
msat_site <- boxplot(He ~ Site, data = msat) #many sites BUT often only one observation per site so not worth keeping

######## build null model (only includes nuisance variables) ########
#variables to include: PrimerNote, CrossSpp, Repeat, position in spp range, (1|MarkerName), (1|Family/Genus/spp), (1|Source), (1|Site)
#dredge to see which ones matter most?

#subset to only those with range_position
msat <- subset(msat, range_position != "NA")

#transform He data to deal with 1s (no 0s in dataset as excluded monomorphic data)
msat$transformed_He <- NA #create column to fill in
for (i in 1:nrow(msat)) { #transform data to handle 1s (Douma & Weedon (2018) Methods in Ecology & Evolution)
  msat$transformed_He[i] <- ((msat$He[i]*(msat$n[i] - 1)) + 0.5)/msat$n[i]
}

#### run null model ####
#null model (no abslat/lat/lon)
beta_msat_null_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + range_position + (1|Family/Genus/spp) + 
                                       (1|Source) + (1|Site), family = beta_family, data = msat)
summary(msat_null_model_full)
#msat_dredge_null <- dredge(msat_null_model_full)
#msat_dredge_null

#try randomly sampling rows from msat to better visualize diagnostics
msat_small <- sample_n(msat, 2500)

beta_msat_null_test <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + range_position + (1|Family/Genus/spp) + 
                                       (1|Source) + (1|Site), family = beta_family, data = msat_small)

beta_msat_null_test_sim <- simulateResiduals(fittedModel = beta_msat_null_test, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(beta_msat_null_test_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(beta_msat_null_test_sim) #residuals against predicted value -- looking for uniformity

plotResiduals(beta_msat_null_test_sim, msat$PrimerNote)
plotResiduals(beta_msat_null_test_sim, msat$CrossSpp)
plotResiduals(beta_msat_null_test_sim, msat$Repeat)
plotResiduals(beta_msat_null_test_sim, msat$range_position)

beta_msat_null_model_nocrossspp <- glmmTMB(transformed_He ~ PrimerNote + Repeat + range_position + (1|Family/Genus/spp) + 
                                       (1|Source) + (1|Site), family = beta_family, data = msat)
beta_msat_null_model_norepeat <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + range_position + (1|Family/Genus/spp) + 
                                       (1|Source) + (1|Site), family = beta_family, data = msat)
beta_msat_null_model_nocnor <- glmmTMB(transformed_He ~ PrimerNote + range_position + (1|Family/Genus/spp) + 
                                       (1|Source) + (1|Site), family = beta_family, data = msat)

beta_msat_null_noc_sim <- simulateResiduals(fittedModel = beta_msat_null_model_nocrossspp, n = 250, plot = F) #creates "DHARMa" residuals from simulations
beta_msat_null_nor_sim <- simulateResiduals(fittedModel = beta_msat_null_model_norepeat, n = 250, plot = F) #creates "DHARMa" residuals from simulations
beta_msat_null_nocnor_sim <- simulateResiduals(fittedModel = beta_msat_null_model_nocnor, n = 250, plot = F) #creates "DHARMa" residuals from simulations

plotQQunif(beta_msat_null_nocnor_sim)
plotResiduals(beta_msat_null_nocnor_sim, smoothScatter = T)

beta_msat_null_sim <- simulateResiduals(fittedModel = beta_msat_null_model_full, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(beta_msat_null_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(beta_msat_null_sim, smoothScatter = T) #residuals against predicted value -- looking for uniformity

plotResiduals(beta_msat_null_sim, msat$PrimerNote)
plotResiduals(beta_msat_null_sim, msat$CrossSpp)
plotResiduals(beta_msat_null_sim, msat$Repeat)
plotResiduals(beta_msat_null_sim, msat$range_position)

#test binomial model again
msat$success <- round(msat$He*msat$n)
msat$failure<- round((1 - msat$He)*msat$n)

binomial_msat_null <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + Repeat + range_position + (1|Family/Genus/spp) + (1|Source) + 
                         (1|Site), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_null_sim <- simulateResiduals(fittedModel = binomial_msat_null, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_null_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_msat_null_sim) #residuals against predicted value -- looking for uniformity
testDispersion(binomial_msat_null_sim)

#try correcting for overdispersion/heteroskedasticity
#add random effect for every unit
msat$ID <- (1:24182)
binomial_msat_null_IDRE <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + Repeat + range_position + (1|Family/Genus/spp) + (1|Source) + 
                              (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_null_IDRE_sim <- simulateResiduals(fittedModel = binomial_msat_null_IDRE, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_null_IDRE_sim)
plotResiduals(binomial_msat_null_IDRE_sim)
testDispersion(binomial_msat_null_IDRE_sim)

plotResiduals(binomial_msat_null_IDRE_sim, msat$PrimerNote)
plotResiduals(binomial_msat_null_IDRE_sim, msat$CrossSpp)
plotResiduals(binomial_msat_null_IDRE_sim, msat$Repeat)
plotResiduals(binomial_msat_null_IDRE_sim, msat$range_position)

binomial_msat_null_IDRE_nocrossspp <- glmer(cbind(success, failure) ~ PrimerNote + Repeat + range_position + (1|Family/Genus/spp) + 
                                             (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)
binomial_msat_null_IDRE_norepeat <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + (1|Family/Genus/spp) + 
                                              (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)
binomial_msat_null_IDRE_nocnor <- glmer(cbind(success, failure) ~ PrimerNote + range_position + (1|Family/Genus/spp) + 
                                              (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)

binomial_msat_null_IDRE_noc_sim <- simulateResiduals(fittedModel = binomial_msat_null_IDRE_nocrossspp, n = 250, plot = F) #creates "DHARMa" residuals from simulations
binomial_msat_null_IDRE_nor_sim <- simulateResiduals(fittedModel = binomial_msat_null_IDRE_norepeat, n = 250, plot = F) #creates "DHARMa" residuals from simulations
binomial_msat_null_IDRE_nocnor_sim <- simulateResiduals(fittedModel = binomial_msat_null_IDRE_nocnor, n = 250, plot = F) #creates "DHARMa" residuals from simulations

plotQQunif(binomial_msat_null_IDRE_nocnor_sim)
plotResiduals(binomial_msat_null_IDRE_nocnor_sim, smoothScatter = T)

anova(binomial_msat_null_IDRE, binomial_msat_null_IDRE_nocrossspp, binomial_msat_null_IDRE_norepeat, binomial_msat_null_IDRE_nocnor)

######## building lat, abslat & lon model ########
#have abslat, lat, lon or some combo of three

#### clean up data ####
#scale geographic variables
msat$lat_scale <- scale(msat$lat)
#msat$lon_scale <- scale(msat$lon)
msat$abslat_scale <- scale(msat$abslat)
  
#convert lon to radians
msat$lon_360 <- msat$lon + 180 #convert (-180,180) to (0,360)
msat$lon_rad <- (2*pi*msat$lon_360)/360

#lat model
beta_msat_lat_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + 
                                      range_position + lat_scale + I(lat_scale^2) + 
                                      (1|Family/Genus/spp) + (1|Source) + (1|Site), 
                                    family = beta_family, data = msat)
#msat_dredge_lat <- dredge(msat_lat_model_full)

binomial_msat_lat_IDRE <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + lat_scale + 
                                  I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                                   (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_lat_IDRE_sim <- simulateResiduals(fittedModel = binomial_msat_lat_IDRE, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_lat_IDRE_sim)
plotResiduals(binomial_msat_lat_IDRE_sim)

#lon model
beta_msat_lon_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + 
                                      range_position + sin(lon_rad) + cos(lon_rad) + 
                                      (1|Family/Genus/spp) + (1|Source) + (1|Site), 
                                    family = beta_family, data = msat)
#msat_dredge_lon <- dredge(msat_lon_model_full)


binomial_msat_lon_IDRE <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + sin(lon_rad) + 
                                  cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                                  (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_lon_IDRE_sim <- simulateResiduals(fittedModel = binomial_msat_lon_IDRE, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_lon_IDRE_sim)
plotResiduals(binomial_msat_lon_IDRE_sim)

#abslat model
beta_msat_abslat_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + 
                                         range_position + abslat_scale + I(abslat_scale^2) + 
                                         (1|Family/Genus/spp) + (1|Source) + (1|Site), 
                                       family = beta_family, data = msat)
#msat_dredge_abslat <- dredge(msat_abslat_model_full)

binomial_msat_abslat_IDRE <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + abslat_scale + 
                                  I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                                  (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_abslat_IDRE_sim <- simulateResiduals(fittedModel = binomial_msat_abslat_IDRE, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_abslat_IDRE_sim)
plotResiduals(binomial_msat_abslat_IDRE_sim)

#lat & lon model
beta_msat_lat_lon_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + 
                                          range_position + lat_scale + I(lat_scale^2) + 
                                          sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + 
                                          (1|Source) + (1|Site), family = beta_family, 
                                        data = msat)
# <- dredge(msat_lat_lon_model_full)

binomial_msat_lat_lon_IDRE <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + lat_scale + 
                                  I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                                  (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_lat_lon_IDRE_sim <- simulateResiduals(fittedModel = binomial_msat_lat_lon_IDRE, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_lat_lon_IDRE_sim)
plotResiduals(binomial_msat_lat_lon_IDRE_sim)

#abslat & lon model
beta_msat_abslat_lon_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + 
                                             range_position + abslat_scale + I(abslat_scale^2) + 
                                             sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + 
                                             (1|Source) + (1|Site), family = beta_family, 
                                           data = msat)
#msat_dredge_abslat_lon <- dredge(msat_abslat_lon_model_full)

binomial_msat_abslat_lon_IDRE <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + abslat_scale + 
                                  I(abslat_scale^2) + cos(lon_rad) + sin(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                                  (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_abslat_lon_IDRE_sim <- simulateResiduals(fittedModel = binomial_msat_abslat_lon_IDRE, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_abslat_lon_IDRE_sim)
plotResiduals(binomial_msat_abslat_lon_IDRE_sim)

anova(binomial_msat_null_IDRE, binomial_msat_lat_IDRE, binomial_msat_lon_IDRE, binomial_msat_abslat_IDRE, 
      binomial_msat_lat_lon_IDRE, binomial_msat_abslat_lon_IDRE)
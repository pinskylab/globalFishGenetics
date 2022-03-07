#Script for building models

remove(list = ls())

library(tidyverse)
library(gridExtra)
library(plotrix)
library(lme4)
library(MuMIn)
library(glmmTMB)
library(DHARMa)
library(effects)

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

#calculate success and failure
mtdna_small_He$success <- round(mtdna_small_He$He*mtdna_small_He$n)
mtdna_small_He$failure<- round((1 - mtdna_small_He$He)*mtdna_small_He$n)

#### run null model ####
#null model (no abslat/lat/lon)
#binomial model

binomial_null <- glmer(cbind(success, failure) ~ bp_scale + range_position + (1|Family/Genus/spp) + (1|Source) + 
                         (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                         control = glmerControl(optimizer = "bobyqa"))

binomial_null_nobp <- glmer(cbind(success, failure) ~ range_position + (1|Family/Genus/spp) + (1|Source) + 
                         (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail",
                         control = glmerControl(optimizer = "bobyqa"))

binomial_null_norp <- glmer(cbind(success, failure) ~ bp_scale + (1|Family/Genus/spp) + (1|Source) + 
                          (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, 
                          na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

binomial_null_nobp_norp <- glmer(cbind(success, failure) ~ (1|Family/Genus/spp) + (1|Source) + 
                              (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", 
                              control = glmerControl(optimizer = "bobyqa"))

#check fit with DHARMa
binomial_null_sim <- simulateResiduals(fittedModel = binomial_null, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_null_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_null_sim) #residuals against predicted value -- looking for uniformity

#against bp & range position
plotResiduals(binomial_null_sim, mtdna_small_He$bp)
plotResiduals(binomial_null_sim, mtdna_small_He$range_position)

#compare binomial w/ and w/out fixed effects
anova(binomial_null, binomial_null_nobp, binomial_null_norp, binomial_null_nobp_norp)

#test dispersion for binomial model
testDispersion(binomial_null_sim)

#look at partial residuals
bp_eff <- effect("range_position", residuals = TRUE, binomial_null) #do with squared value? seems to be the same
plot(bp_eff, smooth.residuals = TRUE)

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

#### lat model ####
binomial_lat <- glmer(cbind(success, failure) ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                         (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

binomial_lat_norp <- glmer(cbind(success, failure) ~ bp_scale + lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                        (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_lat_sim <- simulateResiduals(fittedModel = binomial_lat, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_lat_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_lat_sim) #residuals against predicted value -- looking for uniformity

plotResiduals(binomial_lat_sim, mtdna_small_He$lat_scale)

#look at partial residuals
lat_eff <- effect("I(lat_scale^2)", residuals = TRUE, binomial_lat_norp) #do with squared value? seems to be the same
plot(lat_eff, smooth.residuals = TRUE)

#### lon model ####
binomial_lon <- glmer(cbind(success, failure) ~ bp_scale + range_position + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                        (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

binomial_lon_norp <- glmer(cbind(success, failure) ~ bp_scale + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                        (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_lon_sim <- simulateResiduals(fittedModel = binomial_lon, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_lon_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_lon_sim) #residuals against predicted value -- looking for uniformity

plotResiduals(binomial_lon_sim, mtdna_small_He$lon_scale)

#look at partial residuals
lon_eff <- effect("cos(lon_rad)", residuals = TRUE, binomial_lon_norp) #not sure how to do this...
plot(lon_eff, smooth.residuals = TRUE)

#### abslat model ####
binomial_abslat <- glmer(cbind(success, failure) ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                        (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

binomial_abslat_norp <- glmer(cbind(success, failure) ~ bp_scale + abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                           (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_abslat_sim <- simulateResiduals(fittedModel = binomial_abslat, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_abslat_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_abslat_sim) #residuals against predicted value -- looking for uniformity

plotResiduals(binomial_abslat_sim, mtdna_small_He$abslat_scale)

#look at partial residuals
abslat_eff <- effect("abslat_scale", residuals = TRUE, binomial_abslat_norp) #doesn't matter which abslat, reads together
plot(abslat_eff, smooth.residuals = TRUE)

#### lat & lon model ####
binomial_lat_lon <- glmer(cbind(success, failure) ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                           (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

binomial_lat_lon_norp <- glmer(cbind(success, failure) ~ bp_scale + lat_scale + I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                            (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_lat_lon_sim <- simulateResiduals(fittedModel = binomial_lat_lon, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_lat_lon_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_lat_lon_sim) #residuals against predicted value -- looking for uniformity
testDispersion(binomial_lat_lon_sim)

#look at partial residuals
lat_eff <- effect("lat_scale", residuals = TRUE, binomial_lat_lon_norp)
lon_eff <- effect("sin(lon_rad)", residuals = TRUE, binomial_lat_lon_norp)
plot(lat_eff, smooth.residuals = TRUE)

#### abslat & lon model ###
binomial_abslat_lon <- glmer(cbind(success, failure) ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                            (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

binomial_abslat_lon_norp <- glmer(cbind(success, failure) ~ bp_scale + abslat_scale + I(abslat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                               (1|Site) + (1|MarkerName), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
binomial_abslat_lon_sim <- simulateResiduals(fittedModel = binomial_abslat_lon, n = 1000, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_abslat_lon_sim) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(binomial_abslat_lon_sim) #residuals against predicted value -- looking for uniformity

#look at partial residuals
abslat_eff <- effect("abslat_scale", residuals = TRUE, binomial_abslat_lon_norp)
lon_eff <- effect("sin(lon_rad)", residuals = TRUE, binomial_abslat_lon_norp)
plot(lat_eff, smooth.residuals = TRUE)
#x_range_position <- as.data.frame(effects_range_position) #if want to turn these into ggplot

#compare models with anova
anova(binomial_null, binomial_abslat, binomial_lat, binomial_lon, binomial_abslat_lon, binomial_lat_lon)

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
null_model_pi <- lmer(logpi ~ bp_scale + range_position + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                             (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

null_model_pi_nobp <- lmer(logpi ~ range_position + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                             (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

null_model_pi_norp <- lmer(logpi ~ bp_scale + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                             (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

null_model_pi_nobp_norp <- lmer(logpi ~ (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                             (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?

#checking fit with DHARMa
null_model_pi_sim_output <- simulateResiduals(null_model_pi, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(null_model_pi_sim_output) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(null_model_pi_sim_output) #residuals against predicted value -- looking for uniformity
testDispersion(null_model_pi_sim_output)

#testing against specific predictors
plotResiduals(null_model_pi_sim_output, mtdna_small_pi$bp_scale)
plotResiduals(null_model_pi_sim_output, mtdna_small_pi$range_position)

#pull p-values
coefs <- data.frame(coef(summary(abslat_lon_model_pi_norp)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#look at partial residuals
rp_eff <- effect("range_position", residuals = TRUE, null_model_pi)
plot(rp_eff, smooth.residuals = TRUE)

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

#### lat model ###
lat_model_pi <- lmer(logpi ~ range_position + lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                          (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

lat_model_pi_norp <- lmer(logpi ~ lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                       (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
lat_model_pi_sim_output <- simulateResiduals(lat_model_pi, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(lat_model_pi_sim_output) #QQplot --> looks like underdispersion? more residuals around 0.5, fewer in tail
plotResiduals(lat_model_pi_sim_output) #residuals against predicted value -- looking for uniformity
#testDispersion(lat_model_pi_sim_output)

#look at partial residuals
lat_eff <- effect("lat_scale", residuals = TRUE, lat_model_pi_norp)
plot(lat_eff, smooth.residuals = TRUE)

#### abslat model ###
abslat_model_pi <- lmer(logpi ~ range_position + abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                            (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

abslat_model_pi_norp <- lmer(logpi ~ abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                          (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
abslat_model_pi_sim_output <- simulateResiduals(abslat_model_full_pi, plot = F)
plotQQunif(abslat_model_pi_sim_output) #QQplot
plotResiduals(abslat_model_pi_sim_output) #residuals

#look at partial residuals
abslat_eff <- effect("abslat_scale", residuals = TRUE, abslat_model_pi_norp)
plot(abslat_eff, smooth.residuals = TRUE)

#l### on model ###
lon_model_pi <- lmer(logpi ~ range_position + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                            (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

lon_model_pi_norp <- lmer(logpi ~ sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                       (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
lon_model_pi_sim_output <- simulateResiduals(lon_model_full_pi, plot = F)
plotQQunif(lon_model_pi_sim_output) #QQplot
plotResiduals(lon_model_pi_sim_output) #residuals

#look at partial residuals
lon_eff <- effect("sin(lon_rad)", residuals = TRUE, lon_model_pi_norp)
plot(lon_eff, smooth.residuals = TRUE)

#### lat & lon model ###
lat_lon_model_pi <- lmer(logpi ~ range_position + lat_scale + I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                            (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

lat_lon_model_pi_norp <- lmer(logpi ~ lat_scale + I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                           (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
lat_lon_model_pi_sim_output <- simulateResiduals(lat_lon_model_full_pi, plot = F)
plotQQunif(lat_lon_model_pi_sim_output) #QQplot
plotResiduals(lat_lon_model_pi_sim_output) #residuals

#look at partial residuals
lat_eff <- effect("lat_scale", residuals = TRUE, lat_lon_model_pi_norp)
lon_eff <- effect("sin(lon_rad)", residuals = TRUE, lat_lon_model_pi_norp)
plot(lon_eff, smooth.residuals = TRUE)

#### abslat & lon model ####
abslat_lon_model_pi <- lmer(logpi ~ range_position + abslat_scale + I(abslat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

abslat_lon_model_pi_norp <- lmer(logpi ~ abslat_scale + I(abslat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                              (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))

#checking fit with DHARMa
abslat_lon_model_pi_sim_output <- simulateResiduals(abslat_lon_model_full_pi, plot = F)
plotQQunif(abslat_lon_model_pi_sim_output) #QQplot
plotResiduals(abslat_lon_model_pi_sim_output) #residuals
testDispersion(abslat_lon_model_pi_sim_output)

#look at partial residuals
abslat_eff <- effect("abslat_scale", residuals = TRUE, abslat_lon_model_pi_norp)
lon_eff <- effect("sin(lon_rad)", residuals = TRUE, abslat_lon_model_pi_norp)
plot(lon_eff, smooth.residuals = TRUE)

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

#calculate successes & failures
msat$success <- round(msat$He*msat$n)
msat$failure<- round((1 - msat$He)*msat$n)

#add random effect for every unit
msat$ID <- (1:24182)

#### run null model ####
#null model (no abslat/lat/lon)
#binomial model

binomial_msat_null <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + Repeat + range_position + (1|Family/Genus/spp) + (1|Source) + 
                              (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_null_norp <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + Repeat + (1|Family/Genus/spp) + (1|Source) + 
                              (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_null_nocrossspp <- glmer(cbind(success, failure) ~ PrimerNote + Repeat + range_position + (1|Family/Genus/spp) + 
                                              (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)

binomial_msat_null_norepeat <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + (1|Family/Genus/spp) + 
                                       (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)

binomial_msat_null_nopn <- glmer(cbind(success, failure) ~ CrossSpp + Repeat + range_position + (1|Family/Genus/spp) + 
                                         (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)

binomial_msat_null_nopn_norepeat <- glmer(cbind(success, failure) ~ CrossSpp + range_position + (1|Family/Genus/spp) + 
                                   (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)

binomial_msat_null_nopn_nocrossspp <- glmer(cbind(success, failure) ~ Repeat + range_position + (1|Family/Genus/spp) + 
                                   (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)

binomial_msat_null_nopn_norp <- glmer(cbind(success, failure) ~ CrossSpp + range_position + (1|Family/Genus/spp) + 
                                   (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)

binomial_msat_null_nocrossspp_norp <- glmer(cbind(success, failure) ~ PrimerNote + Repeat + (1|Family/Genus/spp) + 
                                         (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)

binomial_msat_null_norepeat_norp <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + (1|Family/Genus/spp) + 
                                       (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)

binomial_msat_null_nocrosspp_norepeat <- glmer(cbind(success, failure) ~ PrimerNote + range_position + (1|Family/Genus/spp) + 
                                          (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)

binomial_msat_null_nocrosspp_norepeat_norp <- glmer(cbind(success, failure) ~ PrimerNote + (1|Family/Genus/spp) + 
                                    (1|Source) + (1|Site) + (1|ID), family = binomial, data = msat)

binomial_msat_null_nocrosspp_norp_nopn <- glmer(cbind(success, failure) ~ Repeat + (1|Family/Genus/spp) + (1|Source) + 
                              (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_null_nocrosspp_norepeat_nopn <- glmer(cbind(success, failure) ~ range_position + (1|Family/Genus/spp) + (1|Source) + 
                              (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_null_nopn_norepeat_norp <- glmer(cbind(success, failure) ~ CrossSpp + (1|Family/Genus/spp) + (1|Source) + 
                              (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_null_nopn_norp_nocrosspp_norepeat <- glmer(cbind(success, failure) ~ (1|Family/Genus/spp) + (1|Source) + 
                              (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

#checking fit with DHARMa
binomial_msat_null_sim <- simulateResiduals(fittedModel = binomial_msat_null, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_null_sim)
plotResiduals(binomial_msat_null_sim)
testDispersion(binomial_msat_null_sim)

plotResiduals(binomial_msat_null_IDRE_sim, msat$PrimerNote)
plotResiduals(binomial_msat_null_IDRE_sim, msat$CrossSpp)
plotResiduals(binomial_msat_null_IDRE_sim, msat$Repeat)
plotResiduals(binomial_msat_null_IDRE_sim, msat$range_position)

#comparing models
anova(binomial_msat_null_IDRE, binomial_msat_null_IDRE_nocrossspp, binomial_msat_null_IDRE_norepeat, binomial_msat_null_IDRE_nocnor)

#look at partial residuals
primernote_eff <- effect("PrimerNote", residuals = TRUE, binomial_msat_null)
crossspp_eff <- effect("CrossSpp", residuals = TRUE, binomial_msat_null)
repeat_eff <- effect("Repeat", residuals = TRUE, binomial_msat_null)
rangepos_eff <- effect("range_position", residuals = TRUE, binomial_msat_null)
plot(rangepos_eff, smooth.residuals = TRUE)

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

#### lat model ####
binomial_msat_lat <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + lat_scale + 
                                  I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                                   (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_lat_norp <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + lat_scale + 
                             I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                             (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

#checking fit with DHARMa
binomial_msat_lat_sim <- simulateResiduals(fittedModel = binomial_msat_lat, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_lat_sim)
plotResiduals(binomial_msat_lat_sim)

#look at partial residuals
lat_eff <- effect("lat_scale", residuals = TRUE, binomial_msat_lat_norp)
plot(lat_eff, smooth.residuals = TRUE)

#### lon model ####
binomial_msat_lon <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + sin(lon_rad) + 
                                  cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                                  (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_lon_norp <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + sin(lon_rad) + 
                             cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                             (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_lon_sim <- simulateResiduals(fittedModel = binomial_msat_lon, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_lon_sim)
plotResiduals(binomial_msat_lon_sim)

#look at partial residuals
lon_eff <- effect("sin(lon_rad)", residuals = TRUE, binomial_msat_lon_norp)
plot(lon_eff, smooth.residuals = TRUE)

#a### bslat model ####
binomial_msat_abslat <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + abslat_scale + 
                                  I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                                  (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_abslat_norp <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + abslat_scale + 
                                     I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + 
                                     (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

#checking fit with DHARMa
binomial_msat_abslat_sim <- simulateResiduals(fittedModel = binomial_msat_abslat, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_abslat_sim)
plotResiduals(binomial_msat_abslat_sim)

#look at partial residuals
abslat_eff <- effect("abslat_scale", residuals = TRUE, binomial_msat_abslat_norp)
plot(abslat_eff, smooth.residuals = TRUE)

#### lat & lon model ####
binomial_msat_lat_lon <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + lat_scale + 
                                  I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                                  (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_lat_lon_norp <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + lat_scale + 
                                 I(lat_scale^2) + sin(lon_rad) + cos(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                                 (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

#checking fit with DHARMa
binomial_msat_lat_lon_sim <- simulateResiduals(fittedModel = binomial_msat_lat_lon, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_lat_lon_sim)
plotResiduals(binomial_msat_lat_lon_sim)

#look at partial residuals
lat_eff <- effect("lat_scale", residuals = TRUE, binomial_msat_lat_lon_norp)
lon_eff <- effect("sin(lon_rad)", residuals = TRUE, binomial_msat_lat_lon_norp)
plot(lat_eff, smooth.residuals = TRUE)

#### abslat & lon model ####
binomial_msat_abslat_lon <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + range_position + abslat_scale + 
                                  I(abslat_scale^2) + cos(lon_rad) + sin(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                                  (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

binomial_msat_abslat_lon_norp <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + abslat_scale + 
                                    I(abslat_scale^2) + cos(lon_rad) + sin(lon_rad) + (1|Family/Genus/spp) + (1|Source) + 
                                    (1|Site) + (1|ID), family = binomial, data = msat, na.action = "na.fail")

#checking fit with DHARMa
binomial_msat_abslat_lon_im <- simulateResiduals(fittedModel = binomial_msat_abslat_lon, n = 250, plot = F) #creates "DHARMa" residuals from simulations
plotQQunif(binomial_msat_abslat_lon_sim)
plotResiduals(binomial_msat_abslat_lon_sim)

#look at partial residuals
abslat_eff <- effect("abslat_scale", residuals = TRUE, binomial_msat_abslat_lon_norp)
lon_eff <- effect("sin(lon_rad)", residuals = TRUE, binomial_msat_abslat_lon_norp)
plot(abslat_eff, smooth.residuals = TRUE)

#comparing models
anova(binomial_msat_null, binomial_msat_lat, binomial_msat_lon, binomial_msat_abslat, 
      binomial_msat_lat_lon, binomial_msat_abslat_lon)
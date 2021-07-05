#Script for building models

remove(list = ls())

library(tidyverse)
library(gridExtra)
library(plotrix)
library(lme4)
library(MuMIn)
library(glmmTMB)

#read in data
mtdna <- read.csv("output/mtdna_assembled_2.csv", stringsAsFactors = FALSE)
msat <- read.csv("output/msatloci.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)
  mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                    'Centroid')], all.x = TRUE)
  msat <- merge(msat, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 'Half_RangeSize', 
                                  'Centroid')], all.x = TRUE)

##############################################################################################################

######## building models for mtdna He ########
  
#build success and failure columns
mtdna$success <- round(mtdna$He*mtdna$n)
mtdna$failure <- round((1 - mtdna$He)*mtdna$n)

mtdna$abslat <- abs(mtdna$lat)

#subset mtdna bc are a few outliers with very high bp -- entire mtdna?
mtdna_small <- subset(mtdna, as.numeric(mtdna$bp) < 2000)

#subset mtdna to remove He = NA columns
mtdna_small_He <- subset(mtdna_small, mtdna_small$He != "NA")

#### checking effect of bp ####
mtdna_small_He$bp <- as.numeric(mtdna_small_He$bp)
He_bp <- ggplot(mtdna_small_He, aes(x = as.numeric(bp), y = He)) + geom_point() + geom_smooth(method = "lm") #slight positive relationship --> more bp = more possible haplotypes
He_bp

mtdna_bp_glm <- glm(cbind(success, failure) ~ as.numeric(bp), family = "binomial", data = mtdna_small_He)
summary(mtdna_bp_glm)

##### checking effect of Coastal v. Pelagic ####
He_CP <- ggplot(mtdna_small_He, aes(x = as.numeric(bp), y = He, colour = Pelagic_Coastal)) + geom_point() + geom_smooth(method = "lm") #slight positive relationship --> more bp = more possible haplotypes
He_CP #looks like have higher div coastal than pelagic but is significant?

He_CP_hist <- ggplot(mtdna_small_He, aes(x = He, colour = Pelagic_Coastal)) + 
  geom_histogram(binwidth = 0.1, fill = "white", alpha = 0.5, position = "identity")
He_CP_hist #overall distribution looks similar, just much higher counts in Coastal

mtdna_CP_glm <- glm(cbind(success, failure) ~ Pelagic_Coastal, family = "binomial", data = mtdna_small_He)
summary(mtdna_CP_glm)

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

He_range <- ggplot(mtdna_small_He, aes(x = range_position, y = He)) + geom_point() + geom_smooth(method = "lm") #doesn't really look like a relationship but also missing some close to range edge
He_range

mtdna_range_glm <- glm(cbind(success, failure) ~ range_position, family = "binomial", data = mtdna_small_He)
summary(mtdna_range_glm)

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

mtdna_MN_glmer <- glmer(cbind(success,failure) ~ (1|MarkerName), family = "binomial", data = mtdna_small_He)
summary(mtdna_MN_glmer)

#### checking effect of spp ####
He_spp <- boxplot(He ~ spp, data = mtdna_small_He) #diff species have diff ranges of He (keep)

spp_color_plot <- ggplot(data = mtdna_small_He, aes(x = as.numeric(bp), y = He, colour = spp)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position = "none")
spp_color_plot

mtdna_spp_glmer <- glmer(cbind(success,failure) ~ (1|spp), family = "binomial", data = mtdna_small_He)
summary(mtdna_spp_glmer)

##### checking effect of site ####
He_site <- boxplot(He ~ Site, data = mtdna_small_He) #many sites BUT often only one observation per site so not worth keeping

site_color_plot <- ggplot(data = mtdna_small_He, aes(x = as.numeric(bp), y = He, colour = Site)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position = "none")
site_color_plot

##### build null model (only includes nuisance variables) ####
#variables to include: bp, lat, lon, (1|spp), position in spp range, Pelagic v. Coastal, body.length?
#dredge to see which ones matter most

#have to scale numerical data bc otherwise throws convergence issues
#issue with scaling is makes interpretation difficult
mtdna_small_He$bp_scale <- scale(as.numeric(mtdna_small_He$bp))
mtdna_small_He <- subset(mtdna_small_He, range_position != "NA")

#transform He data to deal with 1s (no 0s in dataset as excluded monomorphic data)
mtdna_small_He$transformed_He <- NA #create column to fill in
for (i in 1:nrow(mtdna_small_He)) { #transform data to handle 1s (Douma & Weedon (2018) Methods in Ecology & Evolution)
  mtdna_small_He$transformed_He[i] <- ((mtdna_small_He$He[i]*(mtdna_small_He$n[i] - 1)) + 0.5)/mtdna_small_He$n[i]
}

#mtdna_small_He$He[mtdna_small_He$He == 1] <- 0.99999 #problem with beta_family -- y can't equal 1 --> set 1 to 0.999999
#mtdna_small_He <- subset(mtdna_small_He, He != "NA")
#mtdna_small_He <- subset(mtdna_small_He, n != "NA")

#null_model_full <- glmer(cbind(success, failure) ~ bp_scale + range_position + (1|spp) + (1|Source) + (1|MarkerName) + 
               #          (1|Site), family = binomial, data = mtdna_small_He, na.action = "na.fail")
dredge_null <- dredge(null_model_full)
dredge_null
beta_null_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + (1|Family/Genus/spp) + (1|Source) + 
                                  (1|Site) + (1|MarkerName), family = beta_family, data = mtdna_small_He)
#not a ton of difference between the two in terms of coefficient, relationships, variance explained by random effects, etc.

######## building lat, abslat & long model for mtdna Hd ########

#scale lat
mtdna_small_He$lat_scale <- scale(mtdna_small_He$lat)
mtdna_small_He$lon_scale <- scale(mtdna_small_He$lon)
mtdna_small_He$abslat_scale <- scale(mtdna_small_He$abslat)

#lat_model_full <- glmer(cbind(success, failure) ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + (1|spp) + (1|Source) + (1|MarkerName) + 
       #                    (1|Site), family = binomial, data = mtdna_small_He, na.action = "na.fail")
beta_lat_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                  (1|Site), family = beta_family, data = mtdna_small_He) #beta_family (betabinomial is binomial)
#dredge_lat <- dredge(lat_model_full)

#lon_model_full <- glmer(cbind(success, failure) ~ bp_scale + range_position + lon_scale + I(lon_scale^2) + (1|spp) + (1|Source) + (1|MarkerName) + 
                    #      (1|Site), family = binomial, data = mtdna_small_He, na.action = "na.fail")
beta_lon_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + lon_scale + I(lon_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), family = beta_family, data = mtdna_small_He) #beta_family (betabinomial is binomial)
#dredge_lon <- dredge(lon_model_full)

#abslat_model_full <- glmer(cbind(success, failure) ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + (1|spp) + (1|Source) + (1|MarkerName) + 
                        #  (1|Site), family = binomial, data = mtdna_small_He, na.action = "na.fail")
beta_abslat_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), family = beta_family, data = mtdna_small_He) #beta_family (betabinomial is binomial)
#dredge_abslat <- dredge(abslat_model_full)

#lat_lon_model_full <- glmer(cbind(success, failure) ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + lon_scale + I(lon_scale^2) + (1|spp) + (1|Source) + (1|MarkerName) + 
                   #       (1|Site), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa")) #switching to faster optimizer --> can check with allFit to see if optimizers converge
#ss <- getME(lat_lon_model_full, c("theta", "fixef")) #trying to restart from previous fit (more iterations to converge in)
#ll_model_2 <- update(lat_lon_model_full, start = ss, control = glmerControl(optCtrl= list(maxfun = 2e4)))
beta_lat_lon_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + lon_scale + I(lon_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                 (1|Site), family = beta_family, data = mtdna_small_He) #beta_family (betabinomial is binomial)
#dredge_lat_lon <- dredge(lat_lon_model_full)

#abslat_lon_model_full <- glmer(cbind(success, failure) ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + lon_scale + I(lon_scale^2) + (1|spp) + (1|Source) + (1|MarkerName) + 
                         #     (1|Site), family = binomial, data = mtdna_small_He, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa")) #switching to faster optimizer --> can check with allFit to see if optimizers converge
beta_abslat_lon_model_full <- glmmTMB(transformed_He ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + lon_scale + I(lon_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                     (1|Site), family = beta_family, data = mtdna_small_He) #beta_family (betabinomial is binomial)
#dredge_abslat_lon <- dredge(abslat_lon_model_full)

##############################################################################################################

######## building null model for mtdna pi ########

#subset mtdna bc are a few outliers with very high bp -- entire mtdna?
mtdna_small_pi <- subset(mtdna, as.numeric(mtdna$bp) < 2000)

#subset mtdna to remove Pi = NA columns
mtdna_small_pi <- subset(mtdna_small_pi, mtdna_small_pi$Pi != "NA")

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
Pi_bp <- ggplot(mtdna_small_pi, aes(x = as.numeric(bp), y = Pi)) + geom_point() + geom_smooth(method = "lm") #slight negative relationship --> more bp = likely to be more similar bp (fewer average pairwise differences)
Pi_bp

mtdna_bp_glm_pi <- glm(Pi ~ as.numeric(bp), data = mtdna_small_pi)
summary(mtdna_bp_glm_pi)

##### checking effect of Coastal v. Pelagic ####
Pi_CP <- ggplot(mtdna_small_pi, aes(x = as.numeric(bp), y = Pi, colour = Pelagic_Coastal)) + geom_point() + geom_smooth(method = "lm") #slight positive relationship --> more bp = more possible haplotypes
Pi_CP #looks like have higher div pelagic than coastal, and diff relationships as well?

Pi_CP_hist <- ggplot(mtdna_small_pi, aes(x = Pi, colour = Pelagic_Coastal)) + 
  geom_histogram(binwidth = 0.01, fill = "white", alpha = 0.5, position = "identity")
Pi_CP_hist #v different distributions -- coastal is lower and pelagic higher BUT pelagic much smaller sample size

mtdna_CP_glm_pi <- glm(Pi ~ Pelagic_Coastal, data = mtdna_small_pi)
summary(mtdna_CP_glm_pi)

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

Pi_range <- ggplot(mtdna_small_pi, aes(x = range_position, y = Pi)) + geom_point() + geom_smooth(method = "lm") #looks like slight negative relationship
Pi_range

mtdna_range_glm_pi <- glm(Pi ~ range_position, data = mtdna_small_pi)
summary(mtdna_range_glm_pi)

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

mtdna_MN_lmer_pi <- lmer(Pi ~ (1|MarkerName), data = mtdna_small_pi)
summary(mtdna_MN_lmer_pi)

#### checking effect of spp ####
Pi_spp <- boxplot(Pi ~ spp, data = mtdna_small_pi) #diff species have diff ranges of pi

spp_color_plot_pi <- ggplot(data = mtdna_small_pi, aes(x = as.numeric(bp), y = Pi, colour = spp)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position = "none")
spp_color_plot_pi

mtdna_spp_lmer_pi <- lmer(Pi ~ (1|spp), data = mtdna_small_pi)
summary(mtdna_spp_lmer_pi)

##### build null model (only includes nuisance variables) ####
#variables to include: bp, lat, lon, (1|spp), position in spp range, Pelagic v. Coastal, body.length?, (1|MarkerName)?
#dredge to see which ones matter most

#have to scale numerical data bc otherwise throws convergence issues
#issue with scaling is makes interpretation difficult
mtdna_small_pi$bp_scale <- scale(as.numeric(mtdna_small_pi$bp))
mtdna_small_pi <- subset(mtdna_small_pi, range_position != "NA")

null_model_full_pi <- lmer(logpi ~ bp_scale + range_position + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                             (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa")) #want REML = FALSE as want maximum-likelihood, bp doesn't have to be scaled here --> should scale it?
dredge_null_pi <- dredge(null_model_full_pi) #not getting same scale/convergence issues as with He
dredge_null_pi #same top models as with He

#pull p-values
coefs <- data.frame(coef(summary(abslat_lon_model_full_pi)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

######## building lat, abslat & long model for mtdna Hd ########

#scale lat
mtdna_small_pi$lat_scale <- scale(mtdna_small_pi$lat)
mtdna_small_pi$lon_scale <- scale(mtdna_small_pi$lon)
mtdna_small_pi$abslat_scale <- scale(mtdna_small_pi$abslat)

lat_model_full_pi <- lmer(logpi ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                          (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
dredge_lat_pi <- dredge(lat_model_full_pi)

abslat_model_full_pi <- lmer(logpi ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                            (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
dredge_abslat_pi <- dredge(abslat_model_full_pi)

lon_model_full_pi <- lmer(logpi ~ bp_scale + range_position + lon_scale + I(lon_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                            (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
dredge_lon_pi <- dredge(lon_model_full_pi)

lat_lon_model_full_pi <- lmer(logpi ~ bp_scale + range_position + lat_scale + I(lat_scale^2) + lon_scale + I(lon_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                            (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
dredge_lat_lon_pi <- dredge(lat_lon_model_full_pi)

abslat_lon_model_full_pi <- lmer(logpi ~ bp_scale + range_position + abslat_scale + I(abslat_scale^2) + lon_scale + I(lon_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|MarkerName) + 
                                (1|Site), REML = FALSE, data = mtdna_small_pi, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa"))
dredge_abslat_lon_pi <- dredge(abslat_lon_model_full_pi)

##############################################################################################################

######## building null model for msat He ########

#remove missing n values
msat$n <- as.numeric(msat$n)
msat <- subset(msat, msat$n != "NA")

#build success and failure columns
msat$success <- round(msat$He*msat$n)
msat$failure <- round((1 - msat$He)*msat$n)

#### checking effect of CrossSpp ####
msat_crossspp <- ggplot(msat, aes(x = CrossSpp, y = He)) + geom_point() + geom_smooth(method = "lm") #slight negative relationship -- slightly lower He with CrossSpp markers (makes sense)
msat_crossspp

msat_crossspp_glm <- glm(cbind(success, failure) ~ CrossSpp, family = "binomial", data = msat)
summary(msat_crossspp_glm) #significant relationship, keep for now

##### checking effect of Repeat ####
msat$Repeat <- as.numeric(msat$Repeat)
msat <- subset(msat, msat$Repeat != "NA")

msat_repeat <- ggplot(msat, aes(x = Repeat, y = He)) + geom_point() + geom_smooth(method = "lm") #positive relationship, though don't know what is going on with the repeats <1?
msat_repeat

msat_repeat_glm <- glm(cbind(success, failure) ~ Repeat, family = "binomial", data = msat)
summary(msat_repeat_glm) #significant relationship, keep for now

##### checking effect of Coastal v. Pelagic ####
msat_CP_hist <- ggplot(msat, aes(x = He, colour = Pelagic_Coastal)) + 
  geom_histogram(binwidth = 0.1, fill = "white", alpha = 0.5, position = "identity")
msat_CP_hist #overall distribution looks similar, just much higher counts in Coastal

msat_CP_glm <- glm(cbind(success, failure) ~ Pelagic_Coastal, family = "binomial", data = msat)
summary(msat_CP_glm)

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

msat_range <- ggplot(msat, aes(x = range_position, y = He)) + geom_point() + geom_smooth(method = "lm") #looks like slight positive relationship?
msat_range

msat_range_glm <- glm(cbind(success, failure) ~ range_position, family = "binomial", data = msat)
summary(msat_range_glm)

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

msat_PN_glm <- glm(cbind(success,failure) ~ PrimerNote, family = "binomial", data = msat)
summary(msat_PN_glm)

#### checking effect of spp ####
msat_spp <- boxplot(He ~ spp, data = msat) #diff species have diff ranges of He (keep)

msat_spp_color_plot <- ggplot(data = msat, aes(x = Repeat, y = He, colour = spp)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position = "none")
msat_spp_color_plot

msat_spp_glmer <- glmer(cbind(success,failure) ~ (1|spp), family = "binomial", data = msat)
summary(msat_spp_glmer)

##### checking effect of site ####
msat_site <- boxplot(He ~ Site, data = msat) #many sites BUT often only one observation per site so not worth keeping

##### build null model (only includes nuisance variables) ####
#variables to include: PrinerNote, CrossSpp, Repeat, lat, lon, (1|spp), position in spp range, Pelagic v. Coastal, body.length?
#dredge to see which ones matter most

#have to scale numerical data bc otherwise throws convergence issues
#issue with scaling is makes interpretation difficult
msat <- subset(msat, range_position != "NA")
#msat$He[msat$He == 1] <- 0.99999 #problem with beta_family -- y can't equal 1 --> set 1 to 0.999999

#transform He data to deal with 1s (no 0s in dataset as excluded monomorphic data)
msat$transformed_He <- NA #create column to fill in
for (i in 1:nrow(msat)) { #transform data to handle 1s (Douma & Weedon (2018) Methods in Ecology & Evolution)
  msat$transformed_He[i] <- ((msat$He[i]*(msat$n[i] - 1)) + 0.5)/msat$n[i]
}

#msat_null_model_full <- glmer(cbind(success,failure) ~ PrimerNote + CrossSpp + Repeat + range_position + (1|spp) + (1|Source) + (1|Site), family = "binomial", data = msat, na.action = "na.fail")
beta_msat_null_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + range_position + (1|Family/Genus/spp) + 
                                       (1|Source) + (1|Site), family = beta_family, data = msat)
msat_dredge_null <- dredge(msat_null_model_full)
msat_dredge_null
summary(msat_null_model_full)
#confint(msat_null_model_full)

######## building lat, abslat & long model for mtdna Hd ########

#scale lat
msat$lat_scale <- scale(msat$lat)
msat$lon_scale <- scale(msat$lon)
msat$abslat <- abs(msat$lat)
  msat$abslat_scale <- scale(msat$abslat)

#msat_lat_model_full <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + Repeat + range_position + lat_scale + I(lat_scale^2) + (1|spp) + (1|Source) + 
 #                         (1|Site), family = binomial, data = msat, na.action = "na.fail")
beta_msat_lat_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + range_position + lat_scale + I(lat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|Site), family = beta_family, data = msat)
#msat_dredge_lat <- dredge(msat_lat_model_full)

#msat_lon_model_full <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + Repeat + range_position + lon_scale + I(lon_scale^2) + (1|spp) + (1|Source) + 
 #                              (1|Site), family = binomial, data = msat, na.action = "na.fail")
beta_msat_lon_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + range_position + lon_scale + I(lon_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|Site), family = beta_family, data = msat)
#msat_dredge_lon <- dredge(msat_lon_model_full)

#msat_abslat_model_full <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + Repeat + range_position + abslat_scale + I(abslat_scale^2) + (1|spp) + (1|Source) + 
 #                              (1|Site), family = binomial, data = msat, na.action = "na.fail")
beta_msat_abslat_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + range_position + abslat_scale + I(abslat_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|Site), family = beta_family, data = msat)
#msat_dredge_abslat <- dredge(msat_abslat_model_full)

#msat_lat_lon_model_full <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + Repeat + range_position + lat_scale + I(lat_scale^2) + lon_scale + I(lon_scale^2) + (1|spp) + (1|Source) + 
 #                              (1|Site), family = binomial, data = msat, na.action = "na.fail")
beta_msat_lat_lon_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + range_position + lat_scale + I(lat_scale^2) + lon_scale + I(lon_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|Site), family = beta_family, data = msat)
# <- dredge(msat_lat_lon_model_full)

#msat_abslat_lon_model_full <- glmer(cbind(success, failure) ~ PrimerNote + CrossSpp + Repeat + range_position + abslat_scale + I(abslat_scale^2) + lon_scale + I(lon_scale^2) + (1|spp) + (1|Source) + 
 #                                  (1|Site), family = binomial, data = msat, na.action = "na.fail")
beta_msat_abslat_lon_model_full <- glmmTMB(transformed_He ~ PrimerNote + CrossSpp + Repeat + range_position + abslat_scale + I(abslat_scale^2) + lon_scale + I(lon_scale^2) + (1|Family/Genus/spp) + (1|Source) + (1|Site), family = beta_family, data = msat)
#msat_dredge_abslat_lon <- dredge(msat_abslat_lon_model_full)
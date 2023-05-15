#supplemental figures

library(tidyverse)
library(scales)

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)
cp_info <- read.csv("output/spp_combined_info.csv", stringsAsFactors = FALSE)
shared_studies <- read.csv("output/sharedstudies.csv", stringsAsFactors = FALSE)

#merge dataframes
mtdna <- merge(mtdna, cp_info[, c('spp', 'Pelagic_Coastal', 'Genus', 'Family', 'Northernmost', 'Southernmost', 
                                  'Half_RangeSize', 'Centroid')], all.x = TRUE)

#subset mtdna bc are a few outliers with very high bp --- entire mtdna (mitogenome?)
mtdna_small <-subset(mtdna, as.numeric(mtdna$bp) < 2000)

######################################################################################################################

######## Clean up dataframe ########

#add position in range
mtdna_small$Centroid <- as.numeric(mtdna_small$Centroid) #change to numeric for calculations
mtdna_small$Half_RangeSize <- as.numeric(mtdna_small$Half_RangeSize)

mtdna_small$range_position <- NA #create column to fill in

for (i in 1:nrow(mtdna_small)) { #calculate distance from range center as percentage (0-1 both sides of centroid, 1 = all the way at range edge)
  mtdna_small$range_position[i] <- abs((mtdna_small$lat[i] - mtdna_small$Centroid[i])/mtdna_small$Half_RangeSize[i])
}

#check those >1 and round to 1 (slight discrepancy btwn aquamaps and reality)
mtdna_check <- subset(mtdna_small, mtdna_small$range_position > 1) #often right at aquamaps limit, round to 1 and keep
mtdna_small$range_position[mtdna_small$range_position > 1] <- 1

#subset to only those with range_position
mtdna_small <- subset(mtdna_small, range_position != "NA")

#subset mtdna to remove Hd = NA or pi = NA columns
mtdna_small_hd <- subset(mtdna_small, mtdna_small$He != "NA")
mtdna_small_pi <- subset(mtdna_small, mtdna_small$Pi != "NA")

#####################################################################################################################

######## diversity correlation plots ########
#1000 x 1000

#### Hd v pi plot ####

hd_pi_corr <- cor(x = mtdna_small$Pi, y = mtdna_small$He, use = "complete.obs", method = "spearman") #r= 0.8179

hd_pi_plot <- ggplot() + 
  geom_point(data = mtdna_small, aes(x = Pi, y = He), 
             alpha = 0.5, size = 8, color = "#0E3C45") + 
  annotate("text", x = 0.00011, y = 0.98, label = "A", size = 20) + 
  annotate("text", x = 0.03, y = 0.02, label = "r = 0.9178", size = 12) +
  scale_x_continuous(trans = "log10", limits = c(0.0001, 0.1), labels = scales::label_log()) +
  ylim(0, 1) +
  xlab("π (mtDNA)") + ylab(bquote(H[d]~"(mtDNA)")) + 
  theme_minimal() + 
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 34),
        axis.title.y = element_text(size = 34),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"))
hd_pi_plot

#### He v pi plot ####
#uses shared study dataframe bc combined mtdna, msat

he_pi_corr <- cor(x = shared_studies$Pi, y = shared_studies$He, use = "complete.obs", method = "spearman") #r= 0.2416

he_pi_plot <- ggplot() + 
  geom_point(data = shared_studies, aes(x = Pi, y = He), 
             alpha = 0.5, size = 8, color = "#0E3C45") + 
  annotate("text", x = 0.00011, y = 0.98, label = "B", size = 20) + 
  annotate("text", x = 0.03, y = 0.02, label = "r = 0.2416", size = 12) +
  scale_x_continuous(trans = "log10", limits = c(0.0001, 0.1), labels = scales::label_log()) +
  ylim(0, 1) +
  xlab("π (mtDNA)") + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme_minimal() + 
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 34),
        axis.title.y = element_text(size = 34),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"))
he_pi_plot

#### He v Hd plot ####
#uses shared study dataframe bc combined mtdna, msat

he_hd_corr <- cor(x = shared_studies$Hd, y = shared_studies$He, use = "complete.obs", method = "spearman") #r= 0.3492

he_hd_plot <- ggplot() + 
  geom_point(data = shared_studies, aes(x = Hd, y = He), 
             alpha = 0.5, size = 8, color = "#0E3C45") + 
  annotate("text", x = 0.02, y = 0.98, label = "C", size = 20) + 
  annotate("text", x = 0.85, y = 0.02, label = "r = 0.3492", size = 12) +
  xlim(0, 1) + 
  ylim(0, 1) +
  xlab(bquote(H[d]~"(mtDNA)")) + ylab(bquote(H[e]~"(nucDNA)")) + 
  theme_minimal() + 
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 4),
        axis.title.x = element_text(size = 34),
        axis.title.y = element_text(size = 34),
        axis.ticks = element_line(color = "black", linewidth = 2),
        axis.text.x = element_text(size = 34, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"))
he_hd_plot

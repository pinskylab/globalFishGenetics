rm(list=ls())

mtdna_spp <- read.csv('output/fbdat_mtdna_2.csv', stringsAsFactors=FALSE)
mtdna_studies <- read.csv('output/mtdna_assembled_2.csv', stringsAsFactors=FALSE)

spp <- unique(mtdna_spp$spp)
studies <- unique(mtdna_studies$Source)

spp_data <- data.frame(spp)
spp_data$species_name <- spp_data$spp
spp_data <- separate(spp_data, spp, c("genus", "species"), sep = " ")
genus <- unique(spp_data$genus)
spp_genus <- data.frame(genus)

#install.packages("rfishbase")
library(rfishbase)
library(tidyverse)

spp_ecology <- ecology(spp)

write.csv(spp_ecology, file='output/ecology_mtdna.csv', row.names=FALSE)
write.csv(spp_genus, file='output/genus_mtdna.csv', row.names=FALSE)
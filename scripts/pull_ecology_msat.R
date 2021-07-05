rm(list=ls())

msat_spp <- read.csv('output/fbdat_msat_2.csv', stringsAsFactors=FALSE)
msat_studies <- read.csv('output/msat_assembled_2.csv', stringsAsFactors=FALSE)

spp <- unique(msat_spp$spp)
studies <- unique(msat_studies$Source)

spp_data <- data.frame(spp)
spp_data$species_name <- spp_data$spp
spp_data <- separate(spp_data, spp, c("genus", "species"), sep = " ")
genus <- unique(spp_data$genus)
spp_genus <- data.frame(genus)

#install.packages("rfishbase")
library(rfishbase)
library(tidyverse)

spp_ecology <- ecology(spp)

write.csv(spp_ecology, file='output/ecology_msat.csv', row.names=FALSE)
write.csv(spp_genus, file='output/genus_msat.csv', row.names=FALSE)
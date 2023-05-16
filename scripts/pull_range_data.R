#################################### Script to Pull Range Data from FishBase #######################################################

#Records Range extent for each species in datasets
#From Fishbase (really, AquaMaps data)

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(rfishbase)
library(tidyverse)

#read in data
msat <- read.csv("output/msat_assembled.csv")
mtdna <- read.csv("output/mtdna_assembled.csv")

##############################################################################################################################################

######## Create shared database ########
#get one list of species (all species in either msat & mtDNA datasets, or both)

#identify unique species
spps_msat <- msat[!duplicated(msat$spp), c('spp', 'CommonName')] #grab unique species and their common name
spps_mtdna <- mtdna[!duplicated(mtdna$spp), c('spp','CommonName')]

#combine species list
spps_list <- rbind(spps_msat, spps_mtdna)
spps_list <- spps_list[order(spps_list$spp, spps_list$CommonName), ]
  spps_list <- spps_list[!duplicated(spps_list$spp), c('spp', 'CommonName')] #remove duplicates (ones that show up in both mtdna and msat)

#################################################################################################################################

######## Pull range data ########

#create empty dataframe to populate
spps_distinfo <- spps_list
  spps_distinfo$Northernmost <- NA
  spps_distinfo$Southernmost <- NA
  spps_distinfo$Easternmost <- NA
  spps_distinfo$Westernmost <- NA

#for loop to pull Northernmost range extent
for (i in 1:nrow(spps_distinfo)){
  cat(paste(i, " ", sep = '')) 
  spps_distinfo$Northernmost[i] <- stocks(spps_distinfo$spp[i], fields = "Northernmost")
  }

#for loop to pull Southernmost range extent
for (i in 1:nrow(spps_distinfo)){
  cat(paste(i, " ", sep = '')) 
  spps_distinfo$Southernmost[i] <- stocks(spps_distinfo$spp[i], fields = "Southermost")
}

#for loop to pull Easternmost range extent
for (i in 1:nrow(spps_distinfo)){
  cat(paste(i, " ", sep = '')) 
  spps_distinfo$Easternmost[i] <- stocks(spps_distinfo$spp[i], fields = "Easternmost")
}

#for loop to pull Westernmost range extent
for (i in 1:nrow(spps_distinfo)){
  cat(paste(i, " ", sep = '')) 
  spps_distinfo$Westernmost[i] <- stocks(spps_distinfo$spp[i], fields = "Westernmost")
}

## identify any species with missing data (really only care about North/South extent)
inds_NS <- is.na(spps_distinfo$Northernmost) | is.na(spps_distinfo$Southernmost)
  sum(inds_NS) #number of missing occurrences
spps_missinginfo_NS <- spps_distinfo[inds_NS, ]

#create dataframe of species with missing information to try and gather by hand
#sometimes this is on FishBase/AquaMaps but not caught in for loop
spps_toget <- data.frame(spps_missinginfo$spp)
  spps_commonname <- data.frame(spps_missinginfo$CommonName)
spps_toget <- cbind(spps_toget, spps_commonname)

## write out dataframe ##
spps_distinfo <- apply(spps_distinfo, 2, as.character) #turn information into characters
write.csv(spps_distinfo, "Output/spps_distinfo.csv")

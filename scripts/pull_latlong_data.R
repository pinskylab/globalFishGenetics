library(rfishbase)
library(tidyverse)

msat <- read.csv("output/msat_assembled.csv")
mtdna <- read.csv("output/mtdna_assembled.csv")

spps_msat <- msat[!duplicated(msat$spp), c('spp', 'CommonName')] #grab unique species and their common name
spps_mtdna <- mtdna[!duplicated(mtdna$spp), c('spp','CommonName')]

spps_list <- rbind(spps_msat, spps_mtdna)
spps_list <- spps_list[order(spps_list$spp, spps_list$CommonName), ]
spps_list <- spps_list[!duplicated(spps_list$spp), c('spp', 'CommonName')]

spps_distinfo <- spps_list
spps_distinfo$Northernmost <- NA
spps_distinfo$Southernmost <- NA
spps_distinfo$Easternmost <- NA
spps_distinfo$Westernmost <- NA

for (i in 1:nrow(spps_distinfo)){
  cat(paste(i, " ", sep = '')) 
  spps_distinfo$Northernmost[i] <- stocks(spps_distinfo$spp[i], fields = "Northernmost")
  }

for (i in 1:nrow(spps_distinfo)){
  cat(paste(i, " ", sep = '')) 
  spps_distinfo$Southernmost[i] <- stocks(spps_distinfo$spp[i], fields = "Southermost")
}

for (i in 1:nrow(spps_distinfo)){
  cat(paste(i, " ", sep = '')) 
  spps_distinfo$Easternmost[i] <- stocks(spps_distinfo$spp[i], fields = "Easternmost")
}

for (i in 1:nrow(spps_distinfo)){
  cat(paste(i, " ", sep = '')) 
  spps_distinfo$Westernmost[i] <- stocks(spps_distinfo$spp[i], fields = "Westernmost")
}

inds <- is.na(spps_distinfo$Northernmost) | is.na(spps_distinfo$Southernmost) | is.na(spps_distinfo$Easternmost) | is.na(spps_distinfo$Westernmost)
sum(inds)
spps_missinginfo <- spps_distinfo[inds, ]

inds_NS <- is.na(spps_distinfo$Northernmost) | is.na(spps_distinfo$Southernmost)
sum(inds_NS)
spps_missinginfo_NS <- spps_distinfo[inds_NS, ]
spps_toget <- data.frame(spps_missinginfo$spp)
spps_commonname <- data.frame(spps_missinginfo$CommonName)
spps_toget <- cbind(spps_toget, spps_commonname)

spps_all <- spps_distinfo[!inds, ]
spps_all <- data.frame(spps_all)

spps_distinfo <- apply(spps_distinfo, 2, as.character)
write.csv(spps_distinfo, "Output/spps_distinfo.csv")

#Script to pull env data from bio oracle

remove(list = ls())

#load libraries
library(sdmpredictors)
library(leaflet)
library(raster)

#read in data
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)
msat <- read.csv("output/msatloci_assembled.csv", stringsAsFactors = FALSE)

##############################################################################################################

######## Gathering data from Bio-ORACLE ########

#list layers in Bio-CLIM database (could also do MARSPEC)
#only want annual (no monthly)
list_layers(datasets = "Bio-ORACLE", monthly = FALSE)
#for SST: Bo_sstmean, BO_sstmax, BO_sstmin, BO_sstrange
#lat lon cell size of 0.0833 degrees (5 arcmin or 9.2 km)
#Satellite (Aqua-MODIS), monthly climatologies

#for oxygen: BO_dissox (80)
#lat lon cell size of 0.0833 degrees (5 arcmin or 9.2 km)
#in situ measurement

#for chlorophyll (primary productivity): BO_chlomax, BO_chlomean (71), BO_chlomin, BO_chlorange
#lat lon cell size of 0.0833 degrees (5 arcmin or 9.2 km)
#Satellite (Aqua-MODIS), monthly climatologies

#download specific layers
sst <- load_layers(c("BO_sstmean", "BO_sstrange", "BO_sstmax", "BO_sstmin"))
ox <- load_layers(c("BO_dissox"))
chloroA <- load_layers(c("BO_chlomean", "BO_chlorange", "BO_chlomax", "BO_chlomin"))

###################################################################################################################

######## Pull environmental data for each site ########

#### mtdna #### 
mtdna_env <- data.frame(X = mtdna$X, Site = mtdna$Site, 
                        sst = extract(sst, mtdna[, c("lon", "lat")]), 
                        ox = extract(ox, mtdna[, c("lon", "lat")]), 
                        chloroA = extract(chloroA, mtdna[, c("lon", "lat")]))

write.csv(mtdna_env, "output/mtdna_env.csv")

#### msat #### 
msat_env <- data.frame(X = msat$X, Site = msat$Site, 
                        sst = extract(sst, msat[, c("lon", "lat")]), 
                        ox = extract(ox, msat[, c("lon", "lat")]), 
                        chloroA = extract(chloroA, msat[, c("lon", "lat")]))

write.csv(msat_env, "output/msat_env.csv")

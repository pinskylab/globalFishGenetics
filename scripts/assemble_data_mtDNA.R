################################################### Script for Assembling mtDNA Data  #######################################################

#assemble the mtDNA data and do basic QA/QC
#starts from the new individual mtDNA files as well as master table Malin put together from earlier work (prior to 2013)
#skipping matching with fishing data for now as don't have access to those datasets

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(tidyverse)

#read in data
mtdna1 <- read.csv("data/mtDNA/Fishery lat mtDNA Complete Database.csv", stringsAsFactors = FALSE) #csv from previous msat assembly
mtdna2 <- read.csv("data/mtDNA/mtdna_2013-2020_data.csv", stringsAsFactors = FALSE)
  mtdna2 <- mtdna2[, !(names(mtdna2) == 'Notes')] #remove notes column from df
  mtdna2 <- subset(mtdna2, lat_deg != "NA" & lon_deg != "NA") #remove extra bottom 3 rows with empty data
mtdna3 <- read.csv("data/mtDNA/mtdna_2013-2020_data_1site.csv", stringsAsFactors = FALSE)
  mtdna3 <- mtdna3[, !(names(mtdna3) == 'Notes')]
srdbmatch <- read.csv("data/srdb_matching/mtdna_to_match.csv", stringsAsFactors=FALSE) #to match genetic data to SRDB stocks

##########################################################################################################################################

######## Check for duplicates ########

#add an indicator for source file
mtdna1$file <- 'mtdna101'
mtdna2$file <- 'mtdna102'
mtdna3$file <- 'mtdna103'

#check for duplicate studies between new and old mtdna dataframes
overlaps <- intersect(c(mtdna2$spp, mtdna3$spp), mtdna1$spp) #species included in old and new datasets
oldmtdna_possdup_studies <- sort(unique(mtdna1$Source[mtdna1$spp %in% overlaps])) #list of studies from previous mtdna df on these species
mtdna_possdup_studies <- sort(unique(c(mtdna2$Source[mtdna2$spp %in% overlaps], 
                                       mtdna3$Source[mtdna3$spp %in% overlaps]))) #list of studies from new mtdna df on these species
#no duplicates, all good

#exclude duplicated data in mtdna1
#Lemair et al. 2005 Journal of Evolutionary Biology 18:70-80 also recorded as Lemair, Versini & Bonhomme 2005
dim(mtdna1) #1332 rows
mtdna1 <- subset(mtdna1, Source != "Lemair et al. 2005 Journal of Evolutionary Biology 18:70-80") #remove duplicate study
dim(mtdna1) #1318 rows

######## Merge dataframes ########

#verify column names match
all(names(mtdna1) == names(mtdna2))
all(names(mtdna1) == names(mtdna3)) #all match

#merge new mtdna dataframes together
mtdna <- rbind(mtdna1, mtdna2, mtdna3)
dim(mtdna) #2117 x 20

######## Merge SRDB stock info ########

#match on spp, Source, & site
#read in fbsci, stockid, lat & lon (lat & lon) for error-checking
names(srdbmatch)[names(srdbmatch) == 'lat'] <- 'lat_srdb'
names(srdbmatch)[names(srdbmatch) == 'lon'] <- 'lon_srdb'

#trim to unique spp, Source, site
inds <- duplicated(srdbmatch[, c('spp', 'Source', 'Country', 'Site', 'CollectionYear')])
inds2 <- duplicated(srdbmatch[, c('spp', 'Source', 'Country', 'Site', 'CollectionYear', 'lat_srdb', 'lon_srdb')])
sum(inds)  #9
sum(inds2) #9: matches at 9, means lat/lon is same for duplicated rows (checked --> duplicated bc diff markers at same location/species, OK to trim)
nrow(srdbmatch) #320
srdbmatch <- srdbmatch[!inds, ] #trimming to unique rows (not trimming anything)
nrow(srdbmatch) #320

#merge SRDB stock info
sort(setdiff(srdbmatch$Source, mtdna$Source)) #papers in srdbmatch that aren't in mtdna (none)
nrow(mtdna) #2131
mtdna <- merge(mtdna, srdbmatch[, c('spp', 'Source', 'Country', 'Site', 'CollectionYear', 'fbsci', 'stockid', 'lat_srdb', 'lon_srdb')], 
               all.x = TRUE, by = c('spp', 'Source', 'Country', 'Site', 'CollectionYear'))
nrow(mtdna) #2117

##########################################################################################################################################

######## Clean newly merged msat dataframe ########

######## Calc lat and lon in decimal degrees ########

#remove lines without lat/lon
inds <- is.na(mtdna$lat_deg) | is.na(mtdna$lon_deg)
sum(inds) #0: good, all have lat/lon

#change minutes & seconds to negative if lat/long is negative
mtdna$lat_min <- ifelse(mtdna$lat_deg < 0, -mtdna$lat_min, mtdna$lat_min) #imp so things sum properly in next step
mtdna$lat_sec <- ifelse(mtdna$lat_deg < 0, -mtdna$lat_sec, mtdna$lat_sec)
mtdna$lon_min <- ifelse(mtdna$lon_deg < 0, -mtdna$lon_min, mtdna$lon_min)
mtdna$lon_sec <- ifelse(mtdna$lon_deg < 0, -mtdna$lon_sec, mtdna$lon_sec)

#turn minutes & seconds into decimals
mtdna$lat <- rowSums(cbind(mtdna$lat_deg, mtdna$lat_min/60, mtdna$lat_sec/3600), na.rm = TRUE)
mtdna$lon <- rowSums(cbind(mtdna$lon_deg, mtdna$lon_min/60, mtdna$lon_sec/3600), na.rm = TRUE)

#remove studies averaged across too broad a lat/lon range
dim(mtdna) #2131 rows
mtdna <- subset(mtdna, Source != "Keskin et al. 2012 Mitochondrial DNA 23(2):62-69") #remove bc all sites averaged together
mtdna <- subset(mtdna, Source != "Le Port and Lavery. 2012 Journal of Heredity 103(2):174-185") #remove bc all sites averaged by country
mtdna <- subset(mtdna, Site != "Western North Atlantic, US mid-Atlantic") #remove bc averaging across mid-Atlantic states
dim(mtdna) #2116 rows

#remove sites in black sea
dim(mtdna) #2116 rows
mtdna <- subset(mtdna, Site != "BS1" & Site != "BS2" & Site != "BS3" & Site != "BS4" & Site != "BL" & Site != "Black Sea (west)")
dim(mtdna) #2110 rows

######## Fix species name mistakes ########

#get list of species names
mtdna$spp <- gsub(' $', '', mtdna$spp) #remove trailing spaces
spps <- sort(unique(mtdna$spp)) #check list against FB

#fix species name mistakes
mtdna$spp[mtdna$spp == "Acanthopagrus schlegeli"] <- "Acanthopagrus schlegelii" #fix spelling error
mtdna$spp[mtdna$spp == "Amphiprion akallopsisos"] <- "Amphiprion akallopisos" #fixed spelling error
mtdna$spp[mtdna$spp == "Carcharinus altimus"] <- "Carcharhinus altimus" #fix spelling error
mtdna$spp[mtdna$spp == "Centrophorus zeehani"] <- "Centrophorus zeehaani" #fix spelling error
mtdna$spp[mtdna$spp == "Amodytes personatus"] <- "Ammodytes personatus" #fix spelling error
mtdna$spp[mtdna$spp == "Crystallichthys matsushimae"] <- "Crystallias matsushimae" #FB calls this Crystallias
mtdna$spp[mtdna$spp == "Dasyatis brevicaudata"] <- "Bathytoshia brevicaudata" #FB calls this Bathytoshia
mtdna$spp[mtdna$spp == "Diplodus sargus sargus"] <- "Diplodus sargus" #get rid of extra sargus
mtdna$spp[mtdna$spp == "Grahamina capito"] <- "Forsterygion capito" #FB calls this Forsterygion
mtdna$spp[mtdna$spp == "Grahamina gymnota"] <- "Forsterygion gymnotum" #FB calls this Forsterygion
mtdna$spp[mtdna$spp == "Grahamina nigripenne"] <- "Forsterygion nigripenne" #FB calls this Forsterygion
mtdna$spp[mtdna$spp == "Heptranchis perlo"] <- "Heptranchias perlo" #fixed spelling error
mtdna$spp[mtdna$spp == "Holocentrus ascensionis"] <- "Holocentrus adscensionis" #fix spelling error
mtdna$spp[mtdna$spp == "Paralichthys olivaceus `"] <- "Paralichthys olivaceus" #get rid of '
mtdna$spp[mtdna$spp == "Pleuragramma antarcticum"] <- "Pleuragramma antarctica" #FB calls this antarctica
mtdna$spp[mtdna$spp == "Pleuronectes herzensteini"] <- "Pseudopleuronectes herzensteini" #FB calls this Pseudopleuronectes
mtdna$spp[mtdna$spp == "Pleuronectes yokohamae"] <- "Pseudopleuronectes yokohamae" #FB calls this Pseudopleuronectes
mtdna$spp[mtdna$spp == "Saurida elongate"] <- "Saurida elongata" #FB calls this elongata
mtdna$spp[mtdna$spp == "Sparus aurata L."] <- "Sparus aurata" #get rid of L.

######## Fix common name mistakes ########

#get list of common names
mtdna$CommonName <- gsub(' $', '', mtdna$CommonName) #remove trailing spaces
com_spps <- mtdna[!duplicated(mtdna[, c('spp', 'CommonName')]), c('spp', 'CommonName')] #grab unique spp/CommonName combos
com_spps <- com_spps[order(com_spps$spp, com_spps$CommonName), ] #order by spp then by COmmonName, check against FB

#fix common name mistakes based on FB suggestion
mtdna$CommonName[mtdna$spp == "Abudefduf saxatilis"] <- "Seargeant-major"
mtdna$CommonName[mtdna$spp == "Acanthopagrus schlegelii"] <- "Blackhead seabream"
mtdna$CommonName[mtdna$spp == "Acanthurus nigroris"] <- "Bluelined surgeonfish"
mtdna$CommonName[mtdna$spp == "Alopias vulpinus"] <- "Thresher"
mtdna$CommonName[mtdna$spp == "Ammodytes personatus"] <- "Pacific sandlance"
mtdna$CommonName[mtdna$spp == "Atherinella brasiliensis"] <- "Brazilian silverside"
mtdna$CommonName[mtdna$spp == "Balistes capriscus"] <- "Grey triggerfish"
mtdna$CommonName[mtdna$spp == "Bathygobius cocosensis"] <- "Cocos frill-goby"
mtdna$CommonName[mtdna$spp == "Bathytoshia brevicaudata"] <- "Short-tail stingray"
mtdna$CommonName[mtdna$spp == "Boreogadus saida"] <- "Polar cod"
mtdna$CommonName[mtdna$spp == "Branchiostegus japonicus"] <- "Horsehead tilefish"
mtdna$CommonName[mtdna$spp == "Caffrogobius caffer"] <- "Banded goby"
mtdna$CommonName[mtdna$spp == "Carcharodon carcharias"] <- "Great white shark"
mtdna$CommonName[mtdna$spp == "Centrophorus harrissoni"] <- "Dumb gulper shark"
mtdna$CommonName[mtdna$spp == "Centrophorus moluccensis"] <- "Smallfin gulper shark"
mtdna$CommonName[mtdna$spp == "Centrophorus zeehaani"] <- "Southern dogfish"
mtdna$CommonName[mtdna$spp == "Centroscymnus crepidater"] <- "Longnose velvet dogfish"
mtdna$CommonName[mtdna$spp == "Chromis viridis"] <- "Blue green damselfish"
mtdna$CommonName[mtdna$spp == "Cleisthenes herzensteini"] <- "Pointhead flounder"
mtdna$CommonName[mtdna$spp == "Collichthys lucidus"] <- "Big head croaker"
mtdna$CommonName[mtdna$spp == "Conger conger"] <- "European conger"
mtdna$CommonName[mtdna$spp == "Dicentrarchus labrax"] <- "European seabass"
mtdna$CommonName[mtdna$spp == "Diplodus vulgaris"] <- "Common two-banded seabream"
mtdna$CommonName[mtdna$spp == "Engraulis encrasicolus"] <- "European anchovy"
mtdna$CommonName[mtdna$spp == "Epinephelus fasciatus"] <- "Blacktip grouper"
mtdna$CommonName[mtdna$spp == "Etelis carbunculus"] <- "Deep-water red snapper"
mtdna$CommonName[mtdna$spp == "Etelis coruscans"] <- "Deepwater longtail red snapper"
mtdna$CommonName[mtdna$spp == "Eviota albolineata"] <- "White-line dwarfgoby"
mtdna$CommonName[mtdna$spp == "Girella punctata"] <- "Largescale blackfish"
mtdna$CommonName[mtdna$spp == "Hippocampus trimaculatus"] <- "Longnose seahorse"
mtdna$CommonName[mtdna$spp == "Larimichthys polyactis"] <- "Yellow croaker"
mtdna$CommonName[mtdna$spp == "Lutjanus fulvus"] <- "Blacktail snapper"
mtdna$CommonName[mtdna$spp == "Lutjanus kasmira"] <- "Common bluestripe snapper"
mtdna$CommonName[mtdna$spp == "Lutjanus synagris"] <- "Lane snapper"
mtdna$CommonName[mtdna$spp == "Merluccius albidus"] <- "Offshore silver hake"
mtdna$CommonName[mtdna$spp == "Merluccius bilinearis"] <- "Silver hake"
mtdna$CommonName[mtdna$spp == "Miichthys miiuy"] <- "Mi-iuy croaker"
mtdna$CommonName[mtdna$spp == "Myripristis berndti"] <- "Blotcheye soldierfish"
mtdna$CommonName[mtdna$spp == "Odonus niger"] <- "Red-toothed triggerfish"
mtdna$CommonName[mtdna$spp == "Oplegnathus fasciatus"] <- "Barred knifejaw"
mtdna$CommonName[mtdna$spp == "Pampus chinensis"] <- "Chinese silver pomfret"
mtdna$CommonName[mtdna$spp == "Paralichthys californicus"] <- "California flounder"
mtdna$CommonName[mtdna$spp == "Paralichthys olivaceus"] <- "Bastard halibut"
mtdna$CommonName[mtdna$spp == "Pomacentrus coelestis"] <- "Neon damselfish"
mtdna$CommonName[mtdna$spp == "Pseudopleuronectes herzensteini"] <- "Yellow striped flounder"
mtdna$CommonName[mtdna$spp == "Pseudopleuronectes yokohamae"] <- "Marbeled flounder"
mtdna$CommonName[mtdna$spp == "Pterois miles"] <- "Devil firefish"
mtdna$CommonName[mtdna$spp == "Raja miraletus"] <- "Brown ray"
mtdna$CommonName[mtdna$spp == "Rhizoprionodon acutus"] <- "Milk shark"
mtdna$CommonName[mtdna$spp == "Sardinella lemuru"] <- "Bali sardinella"
mtdna$CommonName[mtdna$spp == "Saurida elongata"] <- "Slender lizardfish"
mtdna$CommonName[mtdna$spp == "Scomberomorus commerson"] <- "Narrow-barred Spanish mackerel"
mtdna$CommonName[mtdna$spp == "Scomberomorus sierra"] <- "Pacific sierra"
mtdna$CommonName[mtdna$spp == "Sebastes emphaeus"] <- "Puget Sound rockfish"
mtdna$CommonName[mtdna$spp == "Sebastes schlegelii"] <- "Korean rockfish"
mtdna$CommonName[mtdna$spp == "Sebastiscus marmoratus"] <- "False kelpfish"
mtdna$CommonName[mtdna$spp == "Seriola lalandi"] <- "Yellowtail amberjack"
mtdna$CommonName[mtdna$spp == "Siganus guttatus"] <- "Orange-spotted spinefoot"
mtdna$CommonName[mtdna$spp == "Sillago japonica"] <- "Japanese sillago"
mtdna$CommonName[mtdna$spp == "Solea solea"] <- "Common sole"
mtdna$CommonName[mtdna$spp == "Sparus aurata"] <- "Gilthead seabream"
mtdna$CommonName[mtdna$spp == "Sphyraena barracuda"] <- "Great barracuda"
mtdna$CommonName[mtdna$spp == "Sphyrna lewini"] <- "Scalloped hammerhead"
mtdna$CommonName[mtdna$spp == "Thalassoma bifasciatum"] <- "Bluehead"
mtdna$CommonName[mtdna$spp == "Xiphias gladius"] <- "Swordfish"
mtdna$CommonName[mtdna$spp == "Zebrasoma flavescens"] <- "Yellow tang"

######## Fix locus name mistakes ########

#get list of locus names
mtdna$MarkerName <- gsub(' $', '', mtdna$MarkerName) #remove trailing spaces
mtdna$MarkerName <- gsub('^ ', '', mtdna$MarkerName) #remove leading space
MarkerName <- sort(unique(mtdna$MarkerName)) #check list for duplicates

#fix locus names (standardize across datasheets & papers)
mtdna$MarkerName[mtdna$MarkerName == "5' end of D-loop region, 3' end of D-region, tRNAPhe gene, and 5' end of 12S rRNA gene"] <- "5' D-loop, 3' D-loop, tRNAPhe, 5' 12S rRNAPhe"
mtdna$MarkerName[mtdna$MarkerName %in% c("ATPase 6 and 8", "ATPase6/8")] <- "ATPase 6, 8"
mtdna$MarkerName[mtdna$MarkerName == "Atpase, cytb"] <- "ATPase, cytb"
mtdna$MarkerName[mtdna$MarkerName %in% c("COI", "cytochrome-c oxidase I", "cytochrome c oxidase subunit 1", "cytochrome oxidase c subunit I")] <- "cytochrome c oxidase subunit I"
mtdna$MarkerName[mtdna$MarkerName == "COIII-ND3-ND4L"] <- "COIII, ND3, ND4L"
mtdna$MarkerName[mtdna$MarkerName %in% c("Control region", "CR")] <- "control region"
mtdna$MarkerName[mtdna$MarkerName %in% c("control region (D-loop)", "D-Loop", "D-loop region", "D-loop sequence")] <- "D-loop"
mtdna$MarkerName[mtdna$MarkerName %in% c("control region and COI", "control region and mitochondrial cytochrome c oxidase I")] <- "control region, cytochrome c oxidase subunit I"
mtdna$MarkerName[mtdna$MarkerName %in% c("Cox1")] <- "cox1"
mtdna$MarkerName[mtdna$MarkerName %in% c("cyt b", "Cyt b", "Cytb", "cytochrome b", "Cytochrome b")] <- "cytb"
mtdna$MarkerName[mtdna$MarkerName == "D-loop partial"] <- "D-loop (partial)"
mtdna$MarkerName[mtdna$MarkerName %in% c("HVR-1", "HVR-I", "hyper-variable region I")] <- "HVR I"
mtdna$MarkerName[mtdna$MarkerName %in% c("NADH-dehydrogenase subunit 4 (ND-4)", "ND4 region")] <- "ND4"
mtdna$MarkerName[mtdna$MarkerName == "ND2 and ND5"] <- "ND2, ND5"
mtdna$MarkerName[mtdna$MarkerName == "ND2/Cytb"] <- "ND2, cytb"
mtdna$MarkerName[mtdna$MarkerName == "transfer RNA Pro and D-loop loci"] <- "tRNAPro, D-loop"
mtdna$MarkerName[mtdna$MarkerName == "tRNA Thr to control region"] <- "tRNAThr to control region"
mtdna$MarkerName[mtdna$MarkerName == "tRNAG, ND6"] <- "major non-coding region, tRNAG, ND6 (partial)" #double-checked paper

##########################################################################################################################################

######## QA/QC ########

######## Character check ########

#were numeric fields read properly?
summary(mtdna) #bp no but bc sometimes a range

#double-check species names visually
t(t(sort(unique(mtdna$spp)))) #print in one column (267 spp)

#double-check locus names visually
t(t(sort(unique(mtdna$MarkerName)))) #print in one column (36 diff markers)

#double-check source names
t(t(sort(unique(mtdna$Source)))) #print in one column (241 studies)

#fix CollectionYear for studies where mis-reported
mtdna$CollectionYear[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "CYS" & mtdna$MarkerName == "control region"] <- "2004-2006"
  mtdna$CollectionYear[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "ISB" & mtdna$MarkerName == "control region"] <- "2007-2008"
  mtdna$CollectionYear[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "SIS" & mtdna$MarkerName == "control region"] <- "2001-2002"
  mtdna$CollectionYear[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "SRK" & mtdna$MarkerName == "control region"] <- "2003-2007"

#fix N for studies where mis-reported
mtdna$n[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "FKS" & mtdna$MarkerName == "control region"] <- 77
  mtdna$n[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "SRK" & mtdna$MarkerName == "control region"] <- 115
  mtdna$n[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "TBB" & mtdna$MarkerName == "control region"] <- 76

######### Check latitude ########

#make sure no instance >90 or <-90
mtdna[mtdna$lat > 90 & !is.na(mtdna$lat), c('spp', 'Source', 'Country', 'Site', 'lat_deg', 'lat_min', 'lat_sec', 'lat')] #0
mtdna[mtdna$lat < -90 & !is.na(mtdna$lat), c('spp', 'Source', 'Country', 'Site', 'lat_deg', 'lat_min', 'lat_sec', 'lat')] #0
  
#check for missing lat data
sum(is.na(mtdna$lat) & !is.na(mtdna$lat_deg)) #anywhere lat is NA but lat_deg is not: 0
inds <- is.na(mtdna$lat)
sum(inds) #0

#histogram of latitude
hist(mtdna$lat) # mostly northern hemisphere, but spans globe

#compare lat to srdbmatch lat, fill in where missing
plot(mtdna$lat, mtdna$lat_srdb); abline(0, 1) ##OK with Dammannagoda et al. 2011 and Bremer et al. 1999, since Malin entered lat/lon in the mtdna sheet
mtdna[is.na(mtdna$lat) & !is.na(mtdna$lat_srdb), ] #0

######## Check longitude ########

#make sure no instance >180, <-180 or == 0
mtdna[mtdna$lon > 180 & !is.na(mtdna$lon), c('spp', 'Source', 'Country', 'Site', 'lon_deg', 'lon_min', 'lon_sec', 'lon')] #0
mtdna[mtdna$lon < -180 & !is.na(mtdna$lon), c('spp', 'Source', 'Country', 'Site', 'lon_deg', 'lon_min', 'lon_sec', 'lon')] #0 
mtdna[mtdna$lon == 0 & !is.na(mtdna$lon), c('spp', 'Source', 'Country', 'Site', 'lon_deg', 'lon_min', 'lon_sec', 'lon')] #1, but OK (off coast of Ghana)

#check for missing lon data
sum(is.na(mtdna$lon) & !is.na(mtdna$lon_deg)) #anywhere lon is NA but lon_deg is not: 0
inds <- is.na(mtdna$lon)
sum(inds) #0

#histogram of longitude
hist(mtdna$lon) #mostly Eastern hemisphere, but spans globe

#compare lon to srdbmatch lon, fill in where missing
plot(mtdna$lon, mtdna$lon_srdb); abline(0, 1) #OK with Tripp-Valdez, Stepien, Scoles et al., Dammannagoda et al., and Bremer et al, since Malin entered lat/lon on mtdna sheet
mtdna[is.na(mtdna$lon) & !is.na(mtdna$lon_srdb), ] #0

######## Fix lat/long based on mtdna map ########

#correct lat/long sites where map shows datapoint not in ocean
mtdna$lat[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Oregon"] <- 46.227506 #based on google maps
  mtdna$lon[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Oregon"] <- -124.230677 #based on google maps
  mtdna$lat[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "California"] <- 38.384647 #based on google maps
  mtdna$lon[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "California"] <- -123.751836 #based on google maps
  mtdna$lat[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Fairweather Island, Gulf of Alaska"] <- 58.175171 #based on google maps
  mtdna$lon[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Fairweather Island, Gulf of Alaska"] <- -136.557556 #based on google maps
  mtdna$lat[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Sitka, Gulf of Alaska"] <- 57.037355 #based on google maps
  mtdna$lon[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Sitka, Gulf of Alaska"] <- -135.462987 #based on google maps
  mtdna$lat[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Vancouver Island, British Columbia"] <- 49.624768 #based on google maps
  mtdna$lon[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Vancouver Island, British Columbia"] <- -124.849541 #based on google maps
mtdna$lat[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "San Diego"] <- 32.74804 #based on google maps
  mtdna$lon[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "San Diego"] <- -117.305393 #based on google maps
  mtdna$lat[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Los Angeles"] <- 33.854146 #based on google maps
  mtdna$lon[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Los Angeles"] <- -118.604253 #based on google maps
  mtdna$lat[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "San Dieguito"] <- 33.007198 #based on google maps
  mtdna$lon[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "San Dieguito"] <- -117.287726 #based on google maps
  mtdna$lat[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Mission Bay"] <- 32.777789 #based on google maps
  mtdna$lon[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Mission Bay"] <- -117.264859 #based on google maps
  mtdna$lat[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Oceanside Harbor"] <- 33.192299 #based on google maps
  mtdna$lon[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Oceanside Harbor"] <- -117.458968 #based on google maps
  mtdna$lat[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "San Diego Bay"] <- 32.6339229 #based on google maps
  mtdna$lon[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "San Diego Bay"] <- -117.190135 #based on google maps
  mtdna$lat[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Tijuana Estuary"] <- 32.529308 #based on google maps
  mtdna$lon[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Tijuana Estuary"] <- -117.132024 #based on google maps
  mtdna$lat[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Agua Hedionda"] <- 33.141301 #based on google maps
  mtdna$lon[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Agua Hedionda"] <- -117.343668 #based on google maps
  mtdna$lat[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Punta Banda"] <- 27.863113 #based on google maps
  mtdna$lon[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Punta Banda"] <- -114.473020 #based on google maps
  mtdna$lat[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Bahia Asunction"] <- 26.987588 #based on google maps
  mtdna$lon[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Bahia Asunction"] <- -114.397554 #based on google maps
  mtdna$Site[mtdna$Site == "Bahia Asunction"] <- "Bahia Asuncion" #fix typo based on paper spelling
  mtdna$lat[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Bahia Magdelena"] <- 24.560852 #based on google maps
  mtdna$lon[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151" & mtdna$Site == "Bahia Magdelena"] <- -112.012381 #based on google maps
mtdna$lat[mtdna$Source == "Sotka et al. 2005 Marine Biotechnology 7:223-230" & mtdna$Site == "Seattle"] <- 47.699 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Sotka et al. 2005 Marine Biotechnology 7:223-230" & mtdna$Site == "Seattle"] <- -122.386667 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Sotka et al. 2005 Marine Biotechnology 7:223-230" & mtdna$Site == "Pt. George"] <- 48.620768 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Sotka et al. 2005 Marine Biotechnology 7:223-230" & mtdna$Site == "Pt. George"] <- -122.585878 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Sotka et al. 2005 Marine Biotechnology 7:223-230" & mtdna$Site == "Shaw Island"] <- 48.654406 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Sotka et al. 2005 Marine Biotechnology 7:223-230" & mtdna$Site == "Shaw Island"] <- -123.075304 #based on paper coordinates
mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Puget Sound"] <- 47.583333 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Puget Sound"] <- -122.5 #based on paper coordinates
mtdna$lat[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "SIN"] <- 23.055571 #based on google maps
  mtdna$lon[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "SIN"] <- -107.275269 #based on google maps
  mtdna$lat[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "MICH"] <- 18.01491 #based on google maps
  mtdna$lon[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "MICH"] <- -103.061221 #based on google maps
  mtdna$lat[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "OAX"] <- 15.987714 #based on google maps
  mtdna$lon[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "OAX"] <- -95.04588 #based on google maps
  mtdna$lat[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "CH"] <- 14.553582 #based on google maps
  mtdna$lon[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "CH"] <- -92.739448 #based on google maps
mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Magdalen Island, QUE"] <- -61.51 #based on paper coordinates (flipping sign)
  mtdna$lat[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Mill River, PEI"] <- 46.475675 #based on google maps
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Mill River, PEI"] <- -64.117766 #based on  google maps
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Joggins, NS"] <- -64.29 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "St. John, NB"] <- -66.08 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "St. Andrews, NB"] <- -67.04 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Mt. Desert Island, ME"] <- -68.2 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Broad Cove, ME"] <- -69.24 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Kittery Point, ME"] <- -70.4 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Narragansett Bay, RI"] <- -71.24 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Waquoit Bay, MA"] <- -70.31 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Patchogue, NJ"] <- -73 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Sandy Hook, NJ"] <- -73.59 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Tuckerton, NJ"] <- -74.2 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Sandy Pt Park, MD"] <- -76.24 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Patterson Park, MD"] <- -76.31 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Silver Beach, VA"] <- -75.57 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Oregon Inlet, NC"] <- -75.31 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Hatteras Inlet, NC"] <- -75.42 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Morehead City, NC"] <- -76.41 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Little River, SC"] <- -78.36 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Pawleys Island, SC"] <- -79.08 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Folly Beach, SC"] <- -79.05 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Tybee Island, GA"] <- -80.57 #based on paper coordinates (flipping sign)
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "Jekyll Island, GA"] <- -81.26 #based on paper coordinates (flipping sign)
  mtdna$lat[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "St. Augustine Beach, FL"] <- 29.430933 #based on google maps
  mtdna$lon[mtdna$Source == "Mach et al. 2011 Marine Biology 158:515-530 " & mtdna$Site == "St. Augustine Beach, FL"] <- -81.073096 #based on google maps
mtdna$lon[mtdna$Source == "Carr et al. 1995" & mtdna$Site == "Heart's Ease Ledge, near Random Island"] <- -53.463483 #based on google maps
  mtdna$lat[mtdna$Source == "Carr et al. 1995" & mtdna$Site == "Flatrock , east coast of Avalon Peninsula north of St. John's"] <- 47.713069 #based on google maps
  mtdna$lon[mtdna$Source == "Carr et al. 1995" & mtdna$Site == "Flatrock , east coast of Avalon Peninsula north of St. John's"] <- -52.696191 #based on google maps
mtdna$lon[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Para"] <- -47.682662 #based on google maps
mtdna$lon[mtdna$Source == "Rocha 2004" & mtdna$Site == "NE Brazil"] <- -34.752947 #based on google maps
  mtdna$lat[mtdna$Source == "Rocha 2004" & mtdna$Site == "Trindade"] <- -20.51854 #based on google maps 
  mtdna$lon[mtdna$Source == "Rocha 2004" & mtdna$Site == "Trindade"] <- -29.334836 #based on google maps 
  mtdna$Country[mtdna$Source == "Rocha 2004" & mtdna$Site == "Trindade"] <- "Brazil"
  mtdna$lat[mtdna$Source == "Rocha 2004" & mtdna$Site == "Bahamas"] <- 27.099024 #based on google maps 
  mtdna$lon[mtdna$Source == "Rocha 2004" & mtdna$Site == "Bahamas"] <- -78.078476 #based on google maps 
  mtdna$lat[mtdna$Source == "Rocha 2004" & mtdna$Site == "St. Croix" | mtdna$Site == "St.Croix"] <- 17.778058 #based on google maps 
  mtdna$lon[mtdna$Source == "Rocha 2004" & mtdna$Site == "St. Croix" | mtdna$Site == "St.Croix"] <- -64.828443 #based on google maps 
mtdna$lat[mtdna$Source == "Cardenas et al. 2009 Fisheries Research 100:109-115" & mtdna$Site == "San Antonio"]<- -33.583333 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Cardenas et al. 2009 Fisheries Research 100:109-115" & mtdna$Site == "San Antonio"] <- -71.616667 #based on paper coordinates
mtdna$lat[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "East Atlantic"] <- 3.833166 #based on google maps
  mtdna$lon[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "East Atlantic"] <- 7.926534 #based on google maps
  mtdna$lat[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "San Blas"] <- 9.568288 #based on google maps
  mtdna$lon[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "San Blas"] <- -78.710984 #based on google maps
mtdna$lat[mtdna$Source == "Shui et. al 2009 Fish Science 75:593-600" & mtdna$Site == "NB"] <- 29.487855 #based on google maps
  mtdna$lon[mtdna$Source == "Shui et. al 2009 Fish Science 75:593-600" & mtdna$Site == "NB"] <- 122.258372 #based on google maps
mtdna$lon[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "South Portugal"] <- -9.25 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Mallorca"] <- 2.333333 #based on paper coordinates & google maps
mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "TGOU"] <- 36.784108 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "TGOU"] <- 10.426040 #based on google maps
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "BANV"] <- 51.291558 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "BANV"] <- 2.798765 #based on google maps
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "FBRE"] <- 48.579197 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "FBRE"] <- -2.689803 #based on google maps
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "FYEU"] <- 46.731184 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "FYEU"] <- -2.328105 #based on google maps
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "FPRE"] <- 43.45 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "FPRE"] <- 3.783333 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "PAVR"] <- 40.659132 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "PAVR"] <- -8.776491 #based on google maps
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "EMAM"] <- 37.763914 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "EMAM"] <- -0.780593 #based on google maps
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "TBIZ"] <- 37.416073 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "TBIZ"] <- 9.815259 #based on google maps
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "AGLA"] <- 36.930153 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "AGLA"] <- 7.777423 #based on google maps
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "MMAR"] <- 35.626061 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "MMAR"] <- -5.256336 #based on google maps
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "MNAD"] <- 35.175243 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "MNAD"] <- -2.906444 #based on google maps
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "MRBT"] <- 34.003488 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "MRBT"] <- -6.921824 #based on google maps
  mtdna$lat[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "TISK"] <- 37.209644 #based on google maps
  mtdna$lon[mtdna$Source == "Lemaire, Versini & Bonhomme 2005" & mtdna$Site == "TISK"] <- 9.025964 #based on google maps
  mtdna$lon[mtdna$Source == "Kochzius Blohm 2005" & mtdna$Site == "Northern Red Sea"] <- 34.179298 #based on google maps
  mtdna$lat[mtdna$Source == "Kochzius Blohm 2005" & mtdna$Site == "Gulf of Aqaba"] <- 29.103332 #based on google maps
  mtdna$lon[mtdna$Source == "Kochzius Blohm 2005" & mtdna$Site == "Gulf of Aqaba"] <- 34.802618 #based on google maps
mtdna$lat[mtdna$Source == "Coscia et al. 2012 Heredity 108:537-546" & mtdna$Site == "Spain"] <- 36.508689 #based on google maps
  mtdna$lon[mtdna$Source == "Coscia et al. 2012 Heredity 108:537-546" & mtdna$Site == "Spain"] <- -6.309368 #based on google maps
  mtdna$lat[mtdna$Source == "Coscia et al. 2012 Heredity 108:537-546" & mtdna$Site == "France"] <- 45.800263 #based on google maps
  mtdna$lon[mtdna$Source == "Coscia et al. 2012 Heredity 108:537-546" & mtdna$Site == "France"] <- -1.441997 #based on google maps
  mtdna$lat[mtdna$Source == "Coscia et al. 2012 Heredity 108:537-546" & mtdna$Site == "Ireland"] <- 52.120047 #based on google maps
  mtdna$lon[mtdna$Source == "Coscia et al. 2012 Heredity 108:537-546" & mtdna$Site == "Ireland"] <- -6.45999 #based on google maps
mtdna$lon[mtdna$Source == "Francisco, Faria, Lengkeek et al. 2011" & mtdna$Site == "Cabo-do-Mundo"] <- -8.72211 #based on google maps
  mtdna$lat[mtdna$Source == "Francisco, Faria, Lengkeek et al. 2011" & mtdna$Site == "Cadiz"] <- 36.579781 #based on google maps
  mtdna$lon[mtdna$Source == "Francisco, Faria, Lengkeek et al. 2011" & mtdna$Site == "Cadiz"] <- -6.354163 #based on google maps
  mtdna$lat[mtdna$Source == "Francisco, Faria, Lengkeek et al. 2011" & mtdna$Site == "Carantec"] <- 48.711922 #based on google maps
  mtdna$lon[mtdna$Source == "Francisco, Faria, Lengkeek et al. 2011" & mtdna$Site == "Carantec"] <- -3.910381 #based on google maps
mtdna$lat[mtdna$Source == "Kumar et al. 2012 Marine Biology Research 8:992-1002" | mtdna$Source == "Kumar et al. 2012 Conservation Genetics 13:1119-1131" & mtdna$Site == "Pondicherry"] <- 11.927044 #based on google maps
  mtdna$lon[mtdna$Source == "Kumar et al. 2012 Marine Biology Research 8:992-1002" | mtdna$Source == "Kumar et al. 2012 Conservation Genetics 13:1119-1131" & mtdna$Site == "Pondicherry"] <- 79.853967 #based on google maps
mtdna$lat[mtdna$Source == "Chen et al. 2010 Zoological Studies 49: 270-282 " & mtdna$Site == "Autumn Taiwan"] <- 24.81293 #based on google maps
  mtdna$lon[mtdna$Source == "Chen et al. 2010 Zoological Studies 49: 270-282 " & mtdna$Site == "Autumn Taiwan"] <- 121.899067 #based on google maps
mtdna$lon[mtdna$Source == "Johnson et al. 2009 Copeia 3: 465-474" & mtdna$Site == "Malibu Lagoon State Beach, CA (LA)"] <- -118.6667 #based on google maps
  mtdna$lon[mtdna$Source == "Johnson et al. 2009 Copeia 3: 465-474" & mtdna$Site == "Doheny State Beach, CA (OC)"] <- -117.6833 #based on google maps
  mtdna$lon[mtdna$Source == "Johnson et al. 2009 Copeia 3: 465-474" & mtdna$Site == "Crown Memorial State Beach, CA (SF)"] <- -122.25 #based on google maps
mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "East China Sea"] <- 38.983333 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "East China Sea"] <- 128.7 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Kodiak Island (AK)"] <- 57.735545 #based on google maps
  mtdna$lon[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Kodiak Island (AK)"] <- -152.373848 #based on google maps
  mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Unimak Pass (AK)"] <- 54.633333 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Unimak Pass (AK)"] <- -168.166667 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Hecate Strait"] <- 53.216667 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Hecate Strait"] <- -130.95 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Near Islands (AK)"] <- 52.516667 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Near Islands (AK)"] <- 173.866667 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Adak Island (AK)"] <- 51.666667 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Adak Island (AK)"] <- -176.6 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Strait of Georgia"] <- 48.9 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Strait of Georgia"] <- -123.1 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Coastal Washington Sea"] <- 47.916667 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Coastal Washington Sea"] <- -125.55 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Sea of Okhotsk"] <- 44.333333 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Sea of Okhotsk"] <- 135.866667 #based on paper coordinates
  mtdna$lat[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Puget Sound"] <- 47.583333 #based on paper coordinates
  mtdna$lon[mtdna$Source == "Canino et al. 2010 Molecular Ecology 19:4339-4351" & mtdna$Site == "Puget Sound"] <- -122.5 #based on paper coordinates
mtdna$lat[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Gangneung"] <- 37.822835 #based on google maps
  mtdna$lon[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Gangneung"] <- 128.97092 #based on google maps
  mtdna$lat[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Teshio"] <- 44.834785 #based on google maps
  mtdna$lon[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Teshio"] <- 141.551989 #based on google maps
  mtdna$lat[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Tomamae"] <- 44.333148 #based on google maps
  mtdna$lon[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Tomamae"] <- 141.465883 #based on google maps
  mtdna$lat[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Onishika"] <- 44.159328 #based on google maps
  mtdna$lon[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Onishika"] <- 141.626534 #based on google maps
  mtdna$lat[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Yoichi"] <- 43.266507 #based on google maps
  mtdna$lon[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Yoichi"] <- 140.805261 #based on google maps
  mtdna$lat[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Tomakomai"] <- 42.489643 #based on google maps
  mtdna$lon[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Tomakomai"] <- 141.607896 #based on google maps
  mtdna$lat[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Erimo"] <- 41.994359 #based on google maps
  mtdna$lon[mtdna$Source == "Kim et al. 2010 Animal Cells and Systems 14:197-206" & mtdna$Site == "Erimo"] <- 143.378365 #based on google maps
mtdna$lat[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Kodial Island"] <- 56.833301 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Kodial Island"] <- -153.096023 #based on google maps
  mtdna$Site[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Kodial Island"] <- "Kodiak Island"
  mtdna$lat[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Hecate Strait"] <- 53.179008 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Hecate Strait"] <- -130.843189 #based on google maps
  mtdna$lat[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Unimak Pass"] <- 54.295662 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Unimak Pass"] <- -164.041092 #based on google maps
  mtdna$lat[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Aleutian Island"] <- 52.467578 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Aleutian Island"] <- -177.445597 #based on google maps
  mtdna$lat[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Gulf of Georgia"] <- 49.17593 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Gulf of Georgia"] <- -123.620235 #based on google maps
  mtdna$lat[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Okhotsk Sea"] <- 49.010155 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Okhotsk Sea"] <- 148.336611 #based on google maps
  mtdna$lat[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Sea of Japan"] <- 42.74775 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Sea of Japan"] <- 137.889463 #based on google maps
  mtdna$lat[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Eastern Hokkaido"] <- 41.96959 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Eastern Hokkaido"] <- 145.265709 #based on google maps
  mtdna$lat[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Yellow Sea"] <- 38.343053 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Yellow Sea"] <- 123.282168 #based on google maps
mtdna$lat[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "YL"] <- 25.22693 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "YL"] <- 121.678158 #based on google maps
  mtdna$lat[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "MA"] <- 25.025705 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "MA"] <- 121.99622 #based on google maps
  mtdna$lat[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "SK"] <- 26.657566 #based on google maps
  mtdna$lon[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "SK"] <- 127.860392 #based on google maps
mtdna$lat[mtdna$Source == "Muki et al. 2009 Ichythology Res 56:380-387" & mtdna$Site == "Wakayama"] <- 34.045771 #based on google maps
  mtdna$lon[mtdna$Source == "Muki et al. 2009 Ichythology Res 56:380-387" & mtdna$Site == "Wakayama"] <- 134.984552 #based on google maps
  mtdna$lat[mtdna$Source == "Muki et al. 2009 Ichythology Res 56:380-387" & mtdna$Site == "Okinawa-jima"] <- 26.902202 #based on google maps
  mtdna$lon[mtdna$Source == "Muki et al. 2009 Ichythology Res 56:380-387" & mtdna$Site == "Okinawa-jima"] <- 128.033732 #based on google maps
  mtdna$lat[mtdna$Source == "Muki et al. 2009 Ichythology Res 56:380-387" & mtdna$Site == "Bonin Islands"] <- 26.914747 #based on google maps
  mtdna$lon[mtdna$Source == "Muki et al. 2009 Ichythology Res 56:380-387" & mtdna$Site == "Bonin Islands"] <- 142.171266 #based on google maps
  mtdna$lat[mtdna$Source == "Muki et al. 2009 Ichythology Res 56:380-387" & mtdna$Site == "Iriomote-jima"] <- 24.348439 #based on google maps
  mtdna$lon[mtdna$Source == "Muki et al. 2009 Ichythology Res 56:380-387" & mtdna$Site == "Iriomote-jima"] <- 123.652159 #based on google maps
  mtdna$Source[mtdna$Source == "Muki et al. 2009 Ichthyology Res 56:380-387"] <- "Mukai et al. 2009 Ichthyology Res 56:380-387"
mtdna$lat[mtdna$Source == "Xiao et. al 2011 Genetics 139:187-198" & mtdna$Site == "Sea of Okhotsk"] <- 44.5 #based on paper coordinates
mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Western Australia"] <- "P01"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P01"] <- -19.7739 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P01"] <- 119.0553 #based on google maps
  mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Northern Territory" & mtdna$n == 57] <- "P02"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P02"] <- -12.164 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P02"] <- 130.6386 #based on google maps
  mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Northern Territory" & mtdna$n == 26] <- "P03"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P03"] <- -13.9478 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P03"] <- 136.2193 #based on google maps
  mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Northern Territory" & mtdna$n == 21] <- "P04"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P04"] <- -15.4743 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P04"] <- 136.9343 #based on google maps
  mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Eastern Gulf of Carpentaria (Queensland Coast)" & mtdna$n == 19] <- "P05"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P05"] <- -16.2962 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P05"] <- 139.2712 #based on google maps
  mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Eastern Gulf of Carpentaria (Queensland Coast)" & mtdna$n == 21] <- "P06"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P06"] <- -16.8227 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P06"] <- 140.7536 #based on google maps
  mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Eastern Gulf of Carpentaria (Queensland Coast)" & mtdna$n == 45] <- "P07"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P07"] <- -15.323 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P07"] <- 141.3204 #based on google maps
  mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Eastern Gulf of Carpentaria (Queensland Coast)" & mtdna$n == 22] <- "P08"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P08"] <- -10.8396 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P08"] <- 141.9656 #based on google maps
  mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Eastern Coast of Queensland" & mtdna$n == 43] <- "P09"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P09"] <- -16.1371 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P09"] <- 145.5669 #based on google maps
  mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Eastern Coast of Queensland" & mtdna$n == 17] <- "P10"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P10"] <- -16.7766 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P10"] <- 145.9205 #based on google maps
  mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Eastern Coast of Queensland" & mtdna$n == 36] <- "P11"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P11"] <- -19.0822 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P11"] <- 146.7426 #based on google maps
  mtdna$Site[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "Eastern Coast of Queensland" & mtdna$n == 27] <- "P12"
  mtdna$lat[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P12"] <- -20.8971 #based on google maps
  mtdna$lon[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P12"] <- 149.0813 #based on google maps
mtdna$lat[mtdna$Source == "Farnsworth et al. 2010 Marine Biology 157:945-953 " & mtdna$Site == "Hicks Reef & Day Reef "] <- -14.478984 #based on google maps
  mtdna$lon[mtdna$Source == "Farnsworth et al. 2010 Marine Biology 157:945-953 " & mtdna$Site == "Hicks Reef & Day Reef "] <- 145.508362 #based on google maps
  mtdna$lat[mtdna$Source == "Farnsworth et al. 2010 Marine Biology 157:945-953 " & mtdna$Site == "North Direction Island & Lizard Island "] <- -14.712794 #based on google maps
  mtdna$lon[mtdna$Source == "Farnsworth et al. 2010 Marine Biology 157:945-953 " & mtdna$Site == "North Direction Island & Lizard Island "] <- 145.483486 #based on google maps
  mtdna$lat[mtdna$Source == "Farnsworth et al. 2010 Marine Biology 157:945-953 " & mtdna$Site == "Martin Reef & Linnet Reef "] <- -14.777422 #based on google maps
  mtdna$lon[mtdna$Source == "Farnsworth et al. 2010 Marine Biology 157:945-953 " & mtdna$Site == "Martin Reef & Linnet Reef "] <- 145.347378 #based on google maps
mtdna$lat[mtdna$Source == "Evans et al. 2010 Fisheries Research 102:16-25" & mtdna$Site == "Palm Islands "] <- -18.710198 #based on google maps
  mtdna$lon[mtdna$Source == "Evans et al. 2010 Fisheries Research 102:16-25" & mtdna$Site == "Palm Islands "] <- 146.550978 #based on google maps
  mtdna$lat[mtdna$Source == "Evans et al. 2010 Fisheries Research 102:16-25" & mtdna$Site == "Whitsunday Islands "] <- -20.237411 #based on google maps
  mtdna$lon[mtdna$Source == "Evans et al. 2010 Fisheries Research 102:16-25" & mtdna$Site == "Whitsunday Islands "] <- 148.925811 #based on google maps
  mtdna$lat[mtdna$Source == "Evans et al. 2010 Fisheries Research 102:16-25" & mtdna$Site == "Capricorn Bunker Islands "] <- -23.279422 #based on google maps
  mtdna$lon[mtdna$Source == "Evans et al. 2010 Fisheries Research 102:16-25" & mtdna$Site == "Capricorn Bunker Islands "] <- 151.83499 #based on google maps
  mtdna$lat[mtdna$Source == "Evans et al. 2010 Fisheries Research 102:16-25" & mtdna$Site == "Keppel Islands "] <- -23.120186 #based on google maps
  mtdna$lon[mtdna$Source == "Evans et al. 2010 Fisheries Research 102:16-25" & mtdna$Site == "Keppel Islands "] <- 150.922916 #based on google maps
mtdna$lat[mtdna$Source == "Hickey et. al 2009 Molecular Ecology 18:680-697" & mtdna$Site == "New Plysmouth"] <-  -39.004456 #based on google maps
  mtdna$lon[mtdna$Source == "Hickey et. al 2009 Molecular Ecology 18:680-697" & mtdna$Site == "New Plysmouth"] <- 174.035392 #based on google maps
  mtdna$Site[mtdna$Source == "Hickey et. al 2009 Molecular Ecology 18:680-697" & mtdna$Site == "New Plysmouth"] <- "New Plymouth"
mtdna$lat[mtdna$Source == "Anderson and Karel 2009 Marine and Coastal Fisheries: Dynamics, Management, and Ecosystem Science 1:121-132" & mtdna$Site == "SL"] <- 29.686301 #based on google maps
  mtdna$lon[mtdna$Source == "Anderson and Karel 2009 Marine and Coastal Fisheries: Dynamics, Management, and Ecosystem Science 1:121-132" & mtdna$Site == "SL"] <- -93.832946 #based on google maps
mtdna$lat[mtdna$Source == "Bowen & Grant 1997" & mtdna$Site == "Iquique"] <- -20.285195 #based on google maps
  mtdna$lon[mtdna$Source == "Bowen & Grant 1997" & mtdna$Site == "Iquique"] <- -70.240529 #based on google maps
  mtdna$lat[mtdna$Source == "Bowen & Grant 1997" & mtdna$Site == "San Diego, Southern California"] <- 32.648754 #based on google maps
  mtdna$lon[mtdna$Source == "Bowen & Grant 1997" & mtdna$Site == "San Diego, Southern California"] <- -117.213337 #based on google maps
mtdna$lat[mtdna$Source == "Cardenas et al. 2009 Fisheries Research 100:109-115" & mtdna$Site == "Iquique"] <- -20.285195 #based on google maps
  mtdna$lon[mtdna$Source == "Cardenas et al. 2009 Fisheries Research 100:109-115" & mtdna$Site == "Iquique"] <- -70.240529 #based on google maps
mtdna$lat[mtdna$Source == "Castillo-Olguin et al. 2012 Ciencias Marinas 38(4):635-652" & mtdna$Site == "Nayarit"] <- 21.16842 #based on google maps
  mtdna$lon[mtdna$Source == "Castillo-Olguin et al. 2012 Ciencias Marinas 38(4):635-652" & mtdna$Site == "Nayarit"] <- -105.46868 #based on google maps
  mtdna$lat[mtdna$Source == "Castillo-Olguin et al. 2012 Ciencias Marinas 38(4):635-652" & mtdna$Site == "Sinaloa"] <- 22.949148 #based on google maps
  mtdna$lon[mtdna$Source == "Castillo-Olguin et al. 2012 Ciencias Marinas 38(4):635-652" & mtdna$Site == "Sinaloa"] <- -106.259761 #based on google maps
mtdna$lat[mtdna$Source == "Cheng et al. 2011 Biochemical Systematics and Ecology 39:718-724" & mtdna$Site == "Ruian "] <- 27.031518 #based on google maps
  mtdna$lon[mtdna$Source == "Cheng et al. 2011 Biochemical Systematics and Ecology 39:718-724" & mtdna$Site == "Ruian "] <- 120.330339 #baased on google maps
mtdna$lat[mtdna$Source == "Chevolot et al. 2007 Mar Biol 151:1275-1286" & mtdna$Site == "Kattegat"] <- 57.600212 #based on google maps
  mtdna$lon[mtdna$Source == "Chevolot et al. 2007 Mar Biol 151:1275-1286" & mtdna$Site == "Kattegat"] <- 10.682336 #based on google maps
  mtdna$lat[mtdna$Source == "Chevolot et al. 2007 Mar Biol 151:1275-1286" & mtdna$Site == "Oxafjordur"] <- 66.335182 #based on google maps
  mtdna$lon[mtdna$Source == "Chevolot et al. 2007 Mar Biol 151:1275-1286" & mtdna$Site == "Oxafjordur"] <- -16.779878 #based on google maps
  mtdna$lat[mtdna$Source == "Chevolot et al. 2007 Mar Biol 151:1275-1286" & mtdna$Site == "Skagafjordur"] <- 66.291407 #based on google maps
  mtdna$lon[mtdna$Source == "Chevolot et al. 2007 Mar Biol 151:1275-1286" & mtdna$Site == "Skagafjordur"] <- -17.229453 #based on google maps
mtdna$lat[mtdna$Source == "Daley et al. 2012 Marine and Freshwater Research 63:708-714" & mtdna$Site == "Victoria"] <- -40.104311 #based on google maps
  mtdna$lon[mtdna$Source == "Daley et al. 2012 Marine and Freshwater Research 63:708-714" & mtdna$Site == "Victoria"] <- 147.931332 #basd on google maps
mtdna$lat[mtdna$Source == "Garcia-Rodriguez et al. 2011 Fisheries Research 107:169-176" & mtdna$Site == "Ensenada "] <- 31.821834 #based on google maps
  mtdna$lon[mtdna$Source == "Garcia-Rodriguez et al. 2011 Fisheries Research 107:169-176" & mtdna$Site == "Ensenada "] <- -116.74386 #based on google maps
mtdna$lat[mtdna$Source == "Goswami et. al 2009 Hydrobiologia 621:213-221" & mtdna$Site == "Kollam"] <- 8.868616 #based on google maps
  mtdna$lon[mtdna$Source == "Goswami et. al 2009 Hydrobiologia 621:213-221" & mtdna$Site == "Kollam"] <- 76.577859 #based on google maps
mtdna$lat[mtdna$Source == "Hickey et. al 2009 Molecular Ecology 18:680-697" & mtdna$Site == "Chatham Islands"] <- -43.796497 #based on google maps
  mtdna$lon[mtdna$Source == "Hickey et. al 2009 Molecular Ecology 18:680-697" & mtdna$Site == "Chatham Islands"] <- -176.285205 #based on google maps
mtdna$lat[mtdna$Source == "Hyde & Vetter 2009 Canadian Journal of Fisheries and Aquatic Sciences 66:1569-1581" & mtdna$Site == "Depoe Bay, Oregon"] <- 44.816667 #based on google maps
  mtdna$lon[mtdna$Source == "Hyde & Vetter 2009 Canadian Journal of Fisheries and Aquatic Sciences 66:1569-1581" & mtdna$Site == "Depoe Bay, Oregon"] <- -124.090539 #based on google maps
  mtdna$lat[mtdna$Source == "Hyde & Vetter 2009 Canadian Journal of Fisheries and Aquatic Sciences 66:1569-1581" & mtdna$Site == "Shelter Cove, Oregon"] <- 40.028397 #based on google maps
  mtdna$lon[mtdna$Source == "Hyde & Vetter 2009 Canadian Journal of Fisheries and Aquatic Sciences 66:1569-1581" & mtdna$Site == "Shelter Cove, Oregon"] <- -124.083687 #based on google maps
mtdna$lat[mtdna$Source == "Kumar et al. 2012 Conservation Genetics 13:1119-1131" & mtdna$Site == "Tuticorin"] <- 8.514448 #based on google maps
  mtdna$lon[mtdna$Source == "Kumar et al. 2012 Conservation Genetics 13:1119-1131" & mtdna$Site == "Tuticorin"] <- 78.148665 #based on google maps
mtdna$lat[mtdna$Source == "Menezez et al. 2012 Journal of Fish Biology 80:2198-2212" & mtdna$Site == "Pondicherry"] <- 11.594979 #based on google maps
  mtdna$lon[mtdna$Source == "Menezez et al. 2012 Journal of Fish Biology 80:2198-2212" & mtdna$Site == "Pondicherry"] <- 79.812898 #based on google maps
mtdna$lat[mtdna$Source == "Nance et. al 2011 PLoS One 6(7):1-12" & mtdna$Site == "MAZ"] <- 23.216408 #based on google maps
  mtdna$lon[mtdna$Source == "Nance et. al 2011 PLoS One 6(7):1-12" & mtdna$Site == "MAZ"] <- -106.441885 #based on google maps
mtdna$lat[mtdna$Source == "Neethling et. al 2008 BMC Evolutionary Biology 8:325" & mtdna$Site == "Haga Haga"] <- -32.764742 #bsed on google maps
  mtdna$lon[mtdna$Source == "Neethling et. al 2008 BMC Evolutionary Biology 8:325" & mtdna$Site == "Haga Haga"] <- 28.254598 #based on google maps
  mtdna$lat[mtdna$Source == "Neethling et. al 2008 BMC Evolutionary Biology 8:325" & mtdna$Site == "Jongensfontein"] <- -34.432198 #bsed on google maps
  mtdna$lon[mtdna$Source == "Neethling et. al 2008 BMC Evolutionary Biology 8:325" & mtdna$Site == "Jongensfontein"] <- 21.344819 #based on google maps
  mtdna$lat[mtdna$Source == "Neethling et. al 2008 BMC Evolutionary Biology 8:325" & mtdna$Site == "Port Alfred"] <- -33.615675 #bsed on google maps
  mtdna$lon[mtdna$Source == "Neethling et. al 2008 BMC Evolutionary Biology 8:325" & mtdna$Site == "Port Alfred"] <- 26.899218 #based on google maps
mtdna$lat[mtdna$Source == "Nunez et. al 2010 Revista de Biologia Marina y Oceanografia 45(1):565-573" & mtdna$Site == "Aysen"] <- -45.374055 #based on google maps
  mtdna$lon[mtdna$Source == "Nunez et. al 2010 Revista de Biologia Marina y Oceanografia 45(1):565-573" & mtdna$Site == "Aysen"] <- -73.627466 #based on google maps
mtdna$lat[mtdna$Source == "Palsson et. al 2009 Polar Biology 32:471-479" & mtdna$Site == "9"] <- 74.144419 #based on google maps
  mtdna$lon[mtdna$Source == "Palsson et. al 2009 Polar Biology 32:471-479" & mtdna$Site == "9"] <- -19.791641 #based on google maps
  mtdna$lat[mtdna$Source == "Palsson et. al 2009 Polar Biology 32:471-479" & mtdna$Site == "10"] <- 75.014057 #based on google maps
  mtdna$lon[mtdna$Source == "Palsson et. al 2009 Polar Biology 32:471-479" & mtdna$Site == "10"] <- -19.537832 #based on google maps
mtdna$lat[mtdna$Source == "Saito et. al 2008 Journal of Fish Biology 73:1937-1945" & mtdna$Site == "Joetsu"] <- 37.230622 #based on google maps
  mtdna$lon[mtdna$Source == "Saito et. al 2008 Journal of Fish Biology 73:1937-1945" & mtdna$Site == "Joetsu"] <- 138.19126 #based on google maps
mtdna$lat[mtdna$Source == "Santa Brigida et. al 2007 Brazilian Journal of Biology 67(4):919-924" & mtdna$Site == "Macapa"] <- 0.53635 #based on google maps
  mtdna$lon[mtdna$Source == "Santa Brigida et. al 2007 Brazilian Journal of Biology 67(4):919-924" & mtdna$Site == "Macapa"] <- -49.824463 #based on google maps
mtdna$lat[mtdna$Source == "Shigenobu et. al 2007 Fisheries Science 73:1104-1112" & mtdna$Site == "AKT"] <- 40.05258 #based on google maps
  mtdna$lon[mtdna$Source == "Shigenobu et. al 2007 Fisheries Science 73:1104-1112" & mtdna$Site == "AKT"] <- 139.831335 #based on google maps
  mtdna$lat[mtdna$Source == "Shigenobu et. al 2007 Fisheries Science 73:1104-1112" & mtdna$Site == "TTR"] <- 35.572204 #based on google maps
  mtdna$lon[mtdna$Source == "Shigenobu et. al 2007 Fisheries Science 73:1104-1112" & mtdna$Site == "TTR"] <- 133.667041 #based on google maps
mtdna$lat[mtdna$Source == "Teske et al. 2005" & mtdna$Site == "Pulai Estuary, Johor, Peninsular Malaysia"] <- 1.790253 #based on google maps
  mtdna$lon[mtdna$Source == "Teske et al. 2005" & mtdna$Site == "Pulai Estuary, Johor, Peninsular Malaysia"] <- 102.80579 #based on google maps
  mtdna$lat[mtdna$Source == "Teske et al. 2005" & mtdna$Site == "Tayabas Bay, Quezon"] <- 13.799611 #based on google maps
  mtdna$lon[mtdna$Source == "Teske et al. 2005" & mtdna$Site == "Tayabas Bay, Quezon"] <- 121.649335 #based on google maps
  mtdna$lat[mtdna$Source == "Teske et. al 2010 Marine Biology 157:2029-2042" & mtdna$Site == "False Bay"] <- -34.12184 #based on google maps
  mtdna$lon[mtdna$Source == "Teske et. al 2010 Marine Biology 157:2029-2042" & mtdna$Site == "False Bay"] <- 18.664569 #based on google maps
mtdna$lat[mtdna$Source == "Vis et al. 1997" & mtdna$Site == "Estuary of Gulf of St. Lawrence"] <- 49.00478 #based on google maps
  mtdna$lon[mtdna$Source == "Vis et al. 1997" & mtdna$Site == "Estuary of Gulf of St. Lawrence"] <- -67.030212 #based on google maps
  mtdna$lat[mtdna$Source == "Vis et al. 1997" & mtdna$Site == "Western Iceland"] <- 65.361062 #based on google maps
  mtdna$lon[mtdna$Source == "Vis et al. 1997" & mtdna$Site == "Western Iceland"] <- -22.480626 #based on google maps
mtdna$lat[mtdna$Source == "Whitney et. al 2012 Journal of Biogeography 39:1144-1156" & mtdna$Site == "Central GBR"] <- -20.065821 #based on google maps
  mtdna$lon[mtdna$Source == "Whitney et. al 2012 Journal of Biogeography 39:1144-1156" & mtdna$Site == "Central GBR"] <- 148.5103 #based on google maps
mtdna$lat[mtdna$Source == "Woodall et. al 2011 Journal of Fish Biology 78:1738-1756" & mtdna$Site == "Alicante"] <- 38.334933 #based on google maps
  mtdna$lon[mtdna$Source == "Woodall et. al 2011 Journal of Fish Biology 78:1738-1756" & mtdna$Site == "Alicante"] <- -0.456393 #based on google maps
  mtdna$lat[mtdna$Source == "Woodall et. al 2011 Journal of Fish Biology 78:1738-1756" & mtdna$Site == "Napoli"] <- 40.722512 #based on google maps
  mtdna$lon[mtdna$Source == "Woodall et. al 2011 Journal of Fish Biology 78:1738-1756" & mtdna$Site == "Napoli"] <- 14.196295 #based on google maps
mtdna$lat[mtdna$Source == "Xiao et. al 2010 Biochemical Genetics 48:402-417" | mtdna$Source == "Xiao et. al 2011 Genetics 139:187-198" & mtdna$Site == "Sea of Okhotsk"] <- 44.354565 #based on google maps
  mtdna$lon[mtdna$Source == "Xiao et. al 2010 Biochemical Genetics 48:402-417" | mtdna$Source == "Xiao et. al 2011 Genetics 139:187-198" & mtdna$Site == "Sea of Okhotsk"] <- 143.812988 #based on google maps
  
######## Check He ########

#make sure no He <0 or >1 (all percentages)
mtdna[(mtdna$He < 0 | mtdna$He > 1) & !is.na(mtdna$He), c('spp', 'Source', 'Country', 'Site', 'He')] #0
  
#make sure not including any monomorphic loci
mtdna <- subset(mtdna, He != 0 | Pi != 0) #both Pi and He are 0
  
#make sure all instances have pi
inds <- is.na(mtdna$He)
sum(inds) #60
mtdna[inds, c('spp', 'Source', 'He', 'Pi')] #OK, all that are NA for He have Pi recorded

#correct He & Hese when reported incorrectly in raw data (often paper has SD and not converted correctly)  
mtdna$Hese[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Mediterranean/Crete"] <- 0.005
  mtdna$Hese[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Madeira"] <- 0.007
  mtdna$Hese[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Azores"] <- 0.008
  mtdna$Hese[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Gulf of Mexico/Louisiana"] <- 0.010 
  mtdna$Hese[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Gulf of Mexico/Panama City, FL"] <- 0.024 
  mtdna$Hese[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Georgia"] <- 0.016
  mtdna$Hese[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "North Carolina"] <- 0.007 
  mtdna$Hese[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "South Carolina"] <- 0.004
mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "East Cantabrian Sea (Basque Country)" & mtdna$MarkerName == "cytb"] <- 0.00783
  mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "East Cantabrian Sea (Basque Country)" & mtdna$MarkerName == "16S"] <- 0.185
  mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "East Cantabrian Sea (Getaria Coast)" & mtdna$MarkerName == "cytb"] <- 0.00367
  mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "East Cantabrian Sea (Getaria Coast)" & mtdna$MarkerName == "16S"] <- 0.0269
  mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "West Cantabrian Sea" & mtdna$MarkerName == "cytb"] <- 0.0113
  mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "West Cantabrian Sea" & mtdna$MarkerName == "16S"] <- 0.0141
  mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "West Cantabrian Sea (Valdearenas)" & mtdna$MarkerName == "cytb"] <- 0.0158
  mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "West Cantabrian Sea (Valdearenas)" & mtdna$MarkerName == "16S"] <- 0.00077
  mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "North Nantes" & mtdna$MarkerName == "cytb"] <- 0.0103
  mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "North Nantes" & mtdna$MarkerName == "16S"] <- 0.00442
  mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "The Adriatic Sea" & mtdna$MarkerName == "cytb"] <- 0.00458
  mtdna$Hese[mtdna$Source == "Borrell et al. 2012 ICES Journal of Marine Science 69:1357-1371" & mtdna$Site == "The Adriatic Sea" & mtdna$MarkerName == "16S"] <- 0.0181
mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P09" & mtdna$MarkerName == "control region"] <- 0.0127
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P10" & mtdna$MarkerName == "control region"] <- 0.0197
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P11" & mtdna$MarkerName == "control region"] <- 0.0139
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P12" & mtdna$MarkerName == "control region"] <- 0.0139
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P09" & mtdna$MarkerName == "ATPase"] <- 0.0132
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P10" & mtdna$MarkerName == "ATPase"] <- 0.0249
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P11" & mtdna$MarkerName == "ATPase"] <- 0.0156
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P12" & mtdna$MarkerName == "ATPase"] <- 0.019
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P05" & mtdna$MarkerName == "control region"] <- 0.0112
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P06" & mtdna$MarkerName == "control region"] <- 0.0175
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P07" & mtdna$MarkerName == "control region"] <- 0.0112
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P08" & mtdna$MarkerName == "control region"] <- 0.0187
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P05" & mtdna$MarkerName == "ATPase"] <- 0.0243
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P06" & mtdna$MarkerName == "ATPase"] <- 0.0234
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P07" & mtdna$MarkerName == "ATPase"] <- 0.0105
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P08" & mtdna$MarkerName == "ATPase"] <- 0.0262
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P02" & mtdna$MarkerName == "control region"] <- 0.0047
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P03" & mtdna$MarkerName == "control region"] <- 0.0175
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P04" & mtdna$MarkerName == "control region"] <- 0.0106
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P02" & mtdna$MarkerName == "ATPase"] <- 0.0075
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P03" & mtdna$MarkerName == "ATPase"] <- 0.0205
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P04" & mtdna$MarkerName == "ATPase"] <- 0.0144
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P01" & mtdna$MarkerName == "control region"] <- 0.0087 
  mtdna$Hese[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P01" & mtdna$MarkerName == "ATPase"] <- 0.0128
mtdna$Hese[mtdna$Source == "Campo & Garcia-Vazquez. 2010 Journal of Sea Research 64:360-368" & mtdna$Site == "Western Coast"] <- 0.0119
  mtdna$Hese[mtdna$Source == "Campo & Garcia-Vazquez. 2010 Journal of Sea Research 64:360-368" & mtdna$Site == "Southern Ireland"] <- 0.0101
  mtdna$Hese[mtdna$Source == "Campo & Garcia-Vazquez. 2010 Journal of Sea Research 64:360-368" & mtdna$Site == "Aegean Sea"] <- 0.0125
  mtdna$Hese[mtdna$Source == "Campo & Garcia-Vazquez. 2010 Journal of Sea Research 64:360-368" & mtdna$Site == "Atlantic"] <- 0.0071
  mtdna$Hese[mtdna$Source == "Campo & Garcia-Vazquez. 2010 Journal of Sea Research 64:360-368" & mtdna$Site == "Cantabric Sea"] <- 0.0109
mtdna$Hese[mtdna$Source == "Correia et al. 2006. Fisheries Science. 72: 20-27" & mtdna$Site == "Azores"] <- 0.0171
  mtdna$Hese[mtdna$Source == "Correia et al. 2006. Fisheries Science. 72: 20-27" & mtdna$Site == "North Portuguese Continental Slope"] <- 0.0219
mtdna$Hese[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Azores"] <- 0.00209
  mtdna$Hese[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Ireland"] <- 0.00272
  mtdna$Hese[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Madeira"] <- 0.00516
  mtdna$Hese[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "North Portugal"] <- 0.00387
  mtdna$Hese[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "South Portugal"] <- 0.00588
  mtdna$Hese[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Mallorca"] <- 0.00809
mtdna$Hese[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151"] <- NA
mtdna$Hese[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Mid-Atlantic Ridge"] <- 0.0047
  mtdna$Hese[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Chatham Rise"] <- 0.072
  mtdna$Hese[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Azores"] <- 0.01
  mtdna$Hese[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Madeira"] <- 0.072
  mtdna$Hese[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Tasman Sea"] <- 0.0264
  mtdna$Hese[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Rockall Trough"] <- 0.0116
mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Bay of Biscay BISB07"] <- 0.072
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Bay of Biscay BISC07"] <- 0.0177
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Celtic Sea CEL08"] <- 0.00799
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "English Channel ENG08"] <- 0.00866
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Irish Sea"] <- 0.00965
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Kagerrak"] <- 0.00336
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Kattegat KATB07"] <- 0.00817
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea BEL08f"] <- 0.0061
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea BEL08j"] <- 0.0114
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea BEL08s"] <- 0.00969
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea GER07"] <- 0.0119
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea NOR08"] <- 0.0165
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea THA07j"] <- 0.0218
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea THA08"] <- 0.00795
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Scheldt estuary ZAN07"] <- 0.00931
  mtdna$Hese[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Wadden Sea TEX07"] <- 0.00808
mtdna$He[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Broome"] <- 0.62
  mtdna$Hese[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Broome"] <- 0.0494
  mtdna$Hese[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Great Barrier Reef"] <- 0.0702
  mtdna$Hese[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Bimini"] <- 0.0047
  mtdna$Hese[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Mariana Islands"] <- 0.0224 
  mtdna$Hese[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Seychelles"] <- 0.0105
  mtdna$Hese[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "San Blas"] <- 0.0119
  mtdna$Hese[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Hawaii"] <- 0.0093
  mtdna$Hese[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Virgin Islands"] <- 0.0878
mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "American Samoa"] <- 0.00709
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Line Islands"] <- 0.00592
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Marshall Islands"] <- 0.0039
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Society Islands"] <- 0.00212
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "French Frigate Shoals"] <- 0.0174
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Gardner Pinnacles"] <- 0.0117
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Hawaii, Big Island"] <- 0.0150
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Johnston Atoll"] <- 0.016
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Kauai"] <- 0.0158
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Kure Atoll"] <- 0.0145
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Laysan  Island"] <- 0.0101
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Lisianski Island"] <- 0.0131
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Maro Reef"] <- 0.108
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Midway Island"] <- 0.0146
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Nihoa"] <- 0.0161
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Oahu"] <- 0.0173
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Pearl & Hermes Reef"] <- 0.0106
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Christmas Island" & mtdna$spp == "Chaetodon meyeri"] <- 0.0126
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Cocos-Keeling Islands" & mtdna$spp == "Chaetodon meyeri"] <- 0.0102 
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Diego Garcia" & mtdna$spp == "Chaetodon meyeri"] <- 0.0017
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Republic of Kiribati" & mtdna$spp == "Chaetodon meyeri"] <- 0.0055 
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Republic of Palau" & mtdna$spp == "Chaetodon meyeri"] <- 0.0756
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Republic of Seychelles" & mtdna$spp == "Chaetodon meyeri"] <- 0.049
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Moorea" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.007
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Christmas Island" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.0053
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Cocos-Keeling Islands" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.0088
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Republic of Kiribati" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.0085
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Nuku Hiva" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.0082
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Palmyra Atoll" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.016
  mtdna$Hese[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Republic of Palau" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.008
mtdna$Hese[mtdna$Source == "Duncan et al. 2006" & mtdna$Site == "Pac Panama"] <- 0.0491
  mtdna$Hese[mtdna$Source == "Duncan et al. 2006" & mtdna$Site == "Seychelles"] <- 0.0387
  mtdna$Hese[mtdna$Source == "Duncan et al. 2006" & mtdna$Site == "South Africa"] <- 0.0234
  mtdna$Hese[mtdna$Source == "Duncan et al. 2006" & mtdna$Site == "Hawaii"] <- 0.0089
mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Cocos Islands"] <- 0.00707
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Moorea"] <- 0.0144
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Christmas Island"] <- 0.0042
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Diego Garcia"] <- 0.00884
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Seychelles"] <- 0.0124
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Kiritimati"] <- 0.00676
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Fiji"] <- 0.00548
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "French Frigate Shoals"] <- 0.0139
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Gardner Pinnacles"] <- 0.016
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Hawaii Island"] <- 0.0157
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Kauai"] <- 0.0118
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Maui"] <- 0.01
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Molokai"] <- 0.0149
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Nihoa"] <- 0.00949
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Oahu"] <- 0.0112
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Pearl & Hermes Reef"] <- 0.0134
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Palau"] <- 0.005
mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Chichi-jima"] <- 0.0106
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "French Frigate Atoll"] <- 0.0047
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Hilo"] <- 0.0062
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Johnston Atoll"] <- 0.0053 
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Kauai"] <- 0.0062
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Kealakekua"] <- 0.0083
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Kure Atoll"] <- 0.0068
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Lanai"] <- 0.0101
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Maro Reef"] <- 0.0078
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Maui"] <- 0.0005
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Midway Atoll"] <- 0.007
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Molokai"] <- 0.0058
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Nihoa"] <- 0.0112
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Oahu"] <- 0.0068
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Peral & Hermes Atoll"] <- 0.0074
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Pohnpei"] <- 0.0186
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Saipan"] <- 0.0087
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "South Point"] <- 0.0066
  mtdna$Hese[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Wawaloli"] <- 0.009
mtdna$Hese[mtdna$Source == "Froukh Koschzius 2007" & mtdna$Site == "Tor"] <- 0.00433
  mtdna$Hese[mtdna$Source == "Froukh Koschzius 2007" & mtdna$Site == "Aqaba"] <- 0.00693
  mtdna$Hese[mtdna$Source == "Froukh Koschzius 2007" & mtdna$Site == "Lith"] <- 0.0032
  mtdna$Hese[mtdna$Source == "Froukh Koschzius 2007" & mtdna$Site == "Thuwal"] <- 0.00442
  mtdna$Hese[mtdna$Source == "Froukh Koschzius 2007" & mtdna$Site == "Hodeidah"] <- 0.00277
mtdna$Hese[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "French Frigate Shoals"] <- 0.00095
  mtdna$Hese[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Hilo"] <- 0.00112
  mtdna$Hese[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Kauai"] <- 0.00167
  mtdna$Hese[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Kona"] <- 0.001
  mtdna$Hese[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Kure"] <- 0.0173
  mtdna$Hese[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Maro"] <- 0.00502
  mtdna$Hese[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Maui Nui"] <- 0.00176
  mtdna$Hese[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Midway"] <- 0.00142
  mtdna$Hese[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Necker"] <- 0.00114
  mtdna$Hese[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Oahu"] <- 0.00085
mtdna$Hese[mtdna$Source == "Gaither et al. 2011 BMC Evolutionary Biology 11:1-16"] <- NA
mtdna$Hese[mtdna$Source == "Gnaither et al. 2011 PLoS ONE 6:1-13"] <- NA
mtdna$Hese[mtdna$Source == "Gaither et al. 2012 Proceedings of the Royal Society B 279:3948-3957"] <- NA
mtdna$Hese[mtdna$Source == "Garcia-Rodriguez et al. 2011 Fisheries Research 107:169-176" & mtdna$Site == "Bahia Magdalena" & mtdna$CollectionYear == 2006] <- 0.000269
  mtdna$Hese[mtdna$Source == "Garcia-Rodriguez et al. 2011 Fisheries Research 107:169-176" & mtdna$Site == "Bahia Magdalena" & mtdna$CollectionYear == 2007] <- 0.00181
  mtdna$Hese[mtdna$Source == "Garcia-Rodriguez et al. 2011 Fisheries Research 107:169-176" & mtdna$Site == "Ensenada "] <- 0.000934
mtdna$Hese[mtdna$Source == "Gotoh et. al 2009 Genes Genetic Systems 84:287-295"] <- NA
mtdna$Hese[mtdna$Source == "Hoolihan Anand Herwherden 2006" & mtdna$Site == "Bahrain"] <- 0.00195
  mtdna$Hese[mtdna$Source == "Hoolihan Anand Herwherden 2006" & mtdna$Site == "Bushehr"] <- 0.00558
  mtdna$Hese[mtdna$Source == "Hoolihan Anand Herwherden 2006" & mtdna$Site == "Kuwait"] <- 0.00327
  mtdna$Hese[mtdna$Source == "Hoolihan Anand Herwherden 2006" & mtdna$Site == "Masirah"] <- 0.00655
  mtdna$Hese[mtdna$Source == "Hoolihan Anand Herwherden 2006" & mtdna$Site == "Abu Dhabi"] <- 0.00062
  mtdna$Hese[mtdna$Source == "Hoolihan Anand Herwherden 2006" & mtdna$Site == "Dibba"] <- 0.0058
  mtdna$Hese[mtdna$Source == "Hoolihan Anand Herwherden 2006" & mtdna$Site == "Ras Al Khaimah"] <- 0.00294
mtdna$Hese[mtdna$Source == "Itoi et al. 2011 Fisheries Science 77: 975-981" & mtdna$Site == "Iwate Prefecture"] <- 0.00221
  mtdna$Hese[mtdna$Source == "Itoi et al. 2011 Fisheries Science 77: 975-981" & mtdna$Site == "Izu Islands"] <- 0.00312
mtdna$Hese[mtdna$Source == "Iwamoto et al. 2012 Fisheries Science 78: 249-257" & mtdna$Site == "Pingtung, Taiwan"] <- 0.00742
  mtdna$Hese[mtdna$Source == "Iwamoto et al. 2012 Fisheries Science 78: 249-257" & mtdna$Site == "Ishigaki"] <- 0.00316
  mtdna$Hese[mtdna$Source == "Iwamoto et al. 2012 Fisheries Science 78: 249-257" & mtdna$Site == "Miyako"] <- 0.00185
  mtdna$Hese[mtdna$Source == "Iwamoto et al. 2012 Fisheries Science 78: 249-257" & mtdna$Site == "Okinawa"] <- 0.00428
  mtdna$Hese[mtdna$Source == "Iwamoto et al. 2012 Fisheries Science 78: 249-257" & mtdna$Site == "Cebu"] <- 0.00153
mtdna$Hese[mtdna$Source == "Johnson et al. 2009 Copeia 3: 465-474" & mtdna$Site == "Crown Memorial State Beach, CA (SF)"] <- 0.0133
  mtdna$Hese[mtdna$Source == "Johnson et al. 2009 Copeia 3: 465-474" & mtdna$Site == "Doheny State Beach, CA (OC)"] <- 0.0167
  mtdna$Hese[mtdna$Source == "Johnson et al. 2009 Copeia 3: 465-474" & mtdna$Site == "Malibu Lagoon State Beach, CA (LA)"] <- 0.0216 
mtdna$Hese[mtdna$Source == "Karl et al. 2012 Marine Biology 159:489-498" & mtdna$Site == "BMN"] <- 0.0284
  mtdna$Hese[mtdna$Source == "Karl et al. 2012 Marine Biology 159:489-498" & mtdna$Site == "BLZ"] <- 0.007
  mtdna$Hese[mtdna$Source == "Karl et al. 2012 Marine Biology 159:489-498" & mtdna$Site == "BRI"] <- 0.0099
  mtdna$Hese[mtdna$Source == "Karl et al. 2012 Marine Biology 159:489-498" & mtdna$Site == "NTL"] <- 0.0232
  mtdna$Hese[mtdna$Source == "Karl et al. 2012 Marine Biology 159:489-498" & mtdna$Site == "FLD"] <- 0.0132
mtdna$Hese[mtdna$Source == "Kim Park Kim 2006" & mtdna$Site == "Jumunjin, East Sea" & mtdna$MarkerName == "control region"] <- 0.0245
  mtdna$Hese[mtdna$Source == "Kim Park Kim 2006" & mtdna$Site == "Sangju, South Sea" & mtdna$MarkerName == "control region"] <- 0.0028
  mtdna$Hese[mtdna$Source == "Kim Park Kim 2006" & mtdna$Site == "Taean, West Sea" & mtdna$MarkerName == "control region"] <- 0.0245
mtdna$Hese[mtdna$Source == "Kochzius Blohm 2005" & mtdna$Site == "Northern Red Sea"] <- 0.00892
  mtdna$Hese[mtdna$Source == "Kochzius Blohm 2005" & mtdna$Site == "Gulf of Aqaba"] <- 0.0064
mtdna$Hese[mtdna$Source == "Kojima et. al 2009 Journal of Oceanography 65:187-193"] <- NA
mtdna$Hese[mtdna$Source == "Limborg et al. 2012 Heredity 109:96-107" & mtdna$Site == "Bay of Biscay"] <- 0.00292
  mtdna$Hese[mtdna$Source == "Limborg et al. 2012 Heredity 109:96-107" & mtdna$Site == "North Sea"] <- 0.00237
  mtdna$Hese[mtdna$Source == "Limborg et al. 2012 Heredity 109:96-107" & mtdna$Site == "Eastern Baltic"] <- 0.00603
  mtdna$Hese[mtdna$Source == "Limborg et al. 2012 Heredity 109:96-107" & mtdna$Site == "Adriatic Sea"] <- 0.0142
  mtdna$Hese[mtdna$Source == "Limborg et al. 2012 Heredity 109:96-107" & mtdna$Site == "Bosporus"] <- 0.00493
  mtdna$Hese[mtdna$Source == "Limborg et al. 2012 Heredity 109:96-107" & mtdna$Site == "Gulf of Lion"] <- 0.00822
mtdna$Hese[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Aleutian Island"] <- 0.12
  mtdna$Hese[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Eastern Hokkaido"] <- 0.07
  mtdna$Hese[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Gulf of Georgia"] <- 0.12
  mtdna$Hese[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Hecate Strait"] <- 0.1
  mtdna$Hese[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Kodiak Island"] <- 0.13
  mtdna$Hese[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Okhotsk Sea"] <- 0.04
  mtdna$Hese[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Sea of Japan"] <- 0.07
  mtdna$Hese[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Unimak Pass"] <- 0.13
  mtdna$Hese[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Yellow Sea"] <- 0.08
  mtdna$Hese[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "MA"] <- 0.0064
  mtdna$Hese[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "SK"] <- 0.0158
  mtdna$Hese[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "YL"] <- 0.0067
mtdna$Hese[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "AU"] <- 0.0156
  mtdna$Hese[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "CE"] <- 0.035
  mtdna$Hese[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "CH"] <- 0.0295
  mtdna$Hese[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "FI"] <- 0.157
  mtdna$Hese[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "KA"] <- 0.0364
  mtdna$Hese[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "KO"] <- 0.0516
  mtdna$Hese[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "KOR"] <- 0.0352
  mtdna$Hese[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "KW"] <- 0.0516
  mtdna$Hese[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "MO"] <- 0.0479
  mtdna$Hese[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "TW"] <- 0.0221
  mtdna$Hese[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "WA"] <- 0.0634
mtdna$Hese[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "CH"] <- 0.000583
  mtdna$Hese[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "MICH"] <- 0.0018
  mtdna$Hese[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "OAX"] <- 0.000566
  mtdna$Hese[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "SIN"] <- 0.000762
mtdna$Hese[mtdna$Source == "Machado-Schiaffino et al. 2010 Molecular Phylogenetics and Evolution 55:552-558" & mtdna$Site == "Ma-C"] <- 0.114
  mtdna$Hese[mtdna$Source == "Machado-Schiaffino et al. 2010 Molecular Phylogenetics and Evolution 55:552-558" & mtdna$Site == "Ma-N"] <- 0.105
  mtdna$Hese[mtdna$Source == "Machado-Schiaffino et al. 2010 Molecular Phylogenetics and Evolution 55:552-558" & mtdna$Site == "Ma-S"] <- 0.081
  mtdna$Hese[mtdna$Source == "Machado-Schiaffino et al. 2010 Molecular Phylogenetics and Evolution 55:552-558" & mtdna$Site == "Mb-CN"] <- 0.022
  mtdna$Hese[mtdna$Source == "Machado-Schiaffino et al. 2010 Molecular Phylogenetics and Evolution 55:552-558" & mtdna$Site == "Mb-CS"] <- 0.055
  mtdna$Hese[mtdna$Source == "Machado-Schiaffino et al. 2010 Molecular Phylogenetics and Evolution 55:552-558" & mtdna$Site == "Mb-N"] <- 0.051
  mtdna$Hese[mtdna$Source == "Machado-Schiaffino et al. 2010 Molecular Phylogenetics and Evolution 55:552-558" & mtdna$Site == "Mb-S"] <- 0.054
mtdna$Hese[mtdna$Source == "McDowell et al. 2007 Gulf and Caribbean Research 19:75-82" & mtdna$Site == "Eastern Atlantic"] <-  0.00495 
  mtdna$Hese[mtdna$Source == "McDowell et al. 2007 Gulf and Caribbean Research 19:75-82" & mtdna$Site == "Caribeean Sea"] <- 0.0118
  mtdna$Hese[mtdna$Source == "McDowell et al. 2007 Gulf and Caribbean Research 19:75-82" & mtdna$Site == "Western North Atlantic, US mid-Atlantic"] <- 0.0062
mtdna$He[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Europa "] <- 0.817
  mtdna$He[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Geyser "] <- 0.869
  mtdna$He[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Glorieuses "] <- 0.892
  mtdna$He[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Juan de Nova "] <- 0.907
  mtdna$He[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Kenya"] <- 0.933
  mtdna$He[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Madagascar"] <- 0.709
  mtdna$He[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Mayotte "] <- 0.886
  mtdna$He[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Moheli "] <- 0.817
  mtdna$He[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Moroni "] <- 0.861
  mtdna$He[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Reunion "] <- 0.864
mtdna$Hese[mtdna$Source == "Rocha et. al 2008 BMC Evolutionary Biology 8:157"] <- NA
mtdna$Hese[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "California"] <- 0.0128
  mtdna$Hese[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Fairweather Island, Gulf of Alaska"] <- 0.0125 
  mtdna$Hese[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Oregon"] <- 0.0051
  mtdna$Hese[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Sitka, Gulf of Alaska"] <- 0.00053
  mtdna$Hese[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Vancouver Island, British Columbia"] <- 0.0307
mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Argentina"] <- 0.0458
  mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Amapa"] <- 0.0481
  mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Bahia"] <- 0.0702
  mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Espirito Santo"] <- 0.0352 
  mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Maranhao"] <- 0.0271
  mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Para"] <- 0.0494
  mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Parana"] <- 0.028
  mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Pernambuco"] <- 0.0326
  mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Rio Grande do Sul"] <- 0.056
  mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Santa Catarina"] <- 0.04
  mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Sao Paulo"] <- 0.0188
  mtdna$Hese[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Venezuela"] <- 0.0481
mtdna$Hese[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "FKS" & mtdna$MarkerName == "control region"] <- 0.52
  mtdna$Hese[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "SRK" & mtdna$MarkerName == "control region"] <- 0.53
  mtdna$Hese[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "TBB" & mtdna$MarkerName == "control region"] <- 0.1
mtdna$Hese[mtdna$Source == "Teske et al. 2005" & mtdna$Site == "Pulai Estuary, Johor, Peninsular Malaysia"] <- 0.00507
  mtdna$Hese[mtdna$Source == "Teske et al. 2005" & mtdna$Site == "Tayabas Bay, Quezon"] <- 0.00118
  mtdna$Hese[mtdna$Source == "Teske et al. 2005" & mtdna$Site == "Tamil Nadu"] <- 0.0122
mtdna$He[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "Kenya"] <- 0.984
  mtdna$He[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "Mauritius"] <- 0.978
  mtdna$He[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "KSeychelles"] <- 0.993
  mtdna$Hese[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "Kenya"] <- 0.00396
  mtdna$Hese[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "Mauritius"] <- 0.00359
  mtdna$Hese[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "KSeychelles"] <- 0.00154
mtdna$Hese[mtdna$Source == "Wilson 2006" & mtdna$Site == "Newport"] <- 0.0205
  mtdna$Hese[mtdna$Source == "Wilson 2006" & mtdna$Site == "Padilla Bay"] <- 0.0301
  mtdna$Hese[mtdna$Source == "Wilson 2006" & mtdna$Site == "San Diego"] <- 0.0278
  mtdna$Hese[mtdna$Source == "Wilson 2006" & mtdna$Site == "Sitka, Alaska"] <- 0.0253
mtdna$Hese[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q1 - center-south of Yellow Sea" & mtdna$MarkerName == "cytb"] <- 0.04
  mtdna$Hese[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q1 - center-south of Yellow Sea" & mtdna$MarkerName == "cytochrome c oxidase subunit I"] <- 0.061
  mtdna$Hese[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q2 - northeast East China Sea" & mtdna$MarkerName == "cytb"] <- 0.041
  mtdna$Hese[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q2 - northeast East China Sea" & mtdna$MarkerName == "cytochrome c oxidase subunit I"] <- 0.031
  mtdna$Hese[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q3 - southwest East China Sea" & mtdna$MarkerName == "cytochrome c oxidase subunit I"] <- 0.039
  mtdna$Hese[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q3 - southwest East China Sea" & mtdna$MarkerName == "cytb"] <- 0.031
mtdna$Hese[mtdna$Source == "Zane et al. 2006" & mtdna$Site == "Elephant Island, South Shetland Islands"] <- 0.011
  mtdna$Hese[mtdna$Source == "Zane et al. 2006" & mtdna$Site == "Halley Bay, Weddell Sea"] <- 0.0064
  mtdna$Hese[mtdna$Source == "Zane et al. 2006" & mtdna$Site == "King George Island, Southern Scotia Arc"] <- 0.0064
  mtdna$Hese[mtdna$Source == "Zane et al. 2006" & mtdna$Site == "Terra Nova Bay, Ross Sea"] <- 0.0038

######## Check Pi ########
  
#make sure no Pi <0 or >1 (all percentages)
mtdna[(mtdna$Pi < 0 | mtdna$Pi > 1) & !is.na(mtdna$Pi), c('spp', 'Source', 'Country', 'Site', 'He')] #0
  
#make sure all instances have Pi
inds <- is.na(mtdna$Pi)
sum(inds) #168
mtdna[inds, c('spp', 'Source', 'He', 'Pi')] #OK, all that are NA for Pi have He recorded
  
#make sure no NA for both He and Pi
sum(is.na(mtdna$Pi) & is.na(mtdna$He)) #0: good

#correct Pi & Pise when reported incorrectly in raw data (often paper has as % or SD and not converted correctly)  
mtdna$Pise[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Mediterranean/Crete"] <- 0.007
  mtdna$Pise[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Madeira"] <- 0.007
  mtdna$Pise[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Azores"] <- 0.007
  mtdna$Pise[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Gulf of Mexico/Louisiana"] <- 0.006 
  mtdna$Pise[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Gulf of Mexico/Panama City, FL"] <- 0.005 
  mtdna$Pise[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "Georgia"] <- 0.005
  mtdna$Pise[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "North Carolina"] <- 0.006 
  mtdna$Pise[mtdna$Source == "Ball et al. 2007 Mar Biol 150:1321-1332" & mtdna$Site == "South Carolina"] <- 0.005
mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P09" & mtdna$MarkerName == "control region"] <- 0.000518
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P10" & mtdna$MarkerName == "control region"] <- 0.00136
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P11" & mtdna$MarkerName == "control region"] <- 0.000667
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P12" & mtdna$MarkerName == "control region"] <- 0.000885
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P09" & mtdna$MarkerName == "ATPase"] <- 0.000107
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P10" & mtdna$MarkerName == "ATPase"] <- 0.000315
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P11" & mtdna$MarkerName == "ATPase"] <- 0.000167
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P12" & mtdna$MarkerName == "ATPase"] <- 0.000231
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P05" & mtdna$MarkerName == "control region"] <- 0.00133
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P06" & mtdna$MarkerName == "control region"] <- 0.00111
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P07" & mtdna$MarkerName == "control region"] <- 0.000671
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P08" & mtdna$MarkerName == "control region"] <- 0.000831
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P05" & mtdna$MarkerName == "ATPase"] <- 0.000275
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P06" & mtdna$MarkerName == "ATPase"] <- 0.000262
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P07" & mtdna$MarkerName == "ATPase"] <- 0.000164
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P08" & mtdna$MarkerName == "ATPase"] <- 0.000192
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P02" & mtdna$MarkerName == "control region"] <- 0.000834
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P03" & mtdna$MarkerName == "control region"] <- 0.00102
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P04" & mtdna$MarkerName == "control region"] <- 0.00153
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P02" & mtdna$MarkerName == "ATPase"] <- 0.0000132
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P03" & mtdna$MarkerName == "ATPase"] <- 0.000216
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P04" & mtdna$MarkerName == "ATPase"] <- 0.00024
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P01" & mtdna$MarkerName == "control region"] <- 0.00184
  mtdna$Pise[mtdna$Source == "Broderick et al. 2011 Journal of Fish Biology 79:633-661" & mtdna$Site == "P01" & mtdna$MarkerName == "ATPase"] <- 0.000321
mtdna$Pise[mtdna$Source == "Correia et al. 2006. Fisheries Science. 72: 20-27" & mtdna$Site == "Azores"] <- 0.00158
  mtdna$Pise[mtdna$Source == "Correia et al. 2006. Fisheries Science. 72: 20-27" & mtdna$Site == "North Portuguese Continental Slope"] <- 0.00212
mtdna$Pi[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Azores"] <- 0.0045
  mtdna$Pi[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Ireland"] <- 0.0050
  mtdna$Pi[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Madeira"] <- 0.0038
  mtdna$Pi[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "North Portugal"] <- 0.0043
  mtdna$Pi[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "South Portugal"] <- 0.0040
  mtdna$Pi[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Mallorca"] <- 0.0021
  mtdna$Pise[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Azores"] <- 0.000281
  mtdna$Pise[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Ireland"] <- 0.000408
  mtdna$Pise[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Madeira"] <- 0.000297
  mtdna$Pise[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "North Portugal"] <- 0.000336
  mtdna$Pise[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "South Portugal"] <- 0.00049
  mtdna$Pise[mtdna$Source == "Correia et al. 2012 Marine Biology 159:1509-1525" & mtdna$Site == "Mallorca"] <- 0.000202
mtdna$Pise[mtdna$Source == "Craig et al. 2011 Bulletin, Southern California Academy of Sciences 110:141-151"] <- NA
mtdna$Pise[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Mid-Atlantic Ridge"] <- 0.00157
  mtdna$Pise[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Chatham Rise"] <- 0.0076
  mtdna$Pise[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Azores"] <- 0.00188
  mtdna$Pise[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Madeira"] <- 0.000894
  mtdna$Pise[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Tasman Sea"] <- 0.000218
  mtdna$Pise[mtdna$Source == "Cunha et al. 2012 PLoS ONE 7:e49196" & mtdna$Site == "Rockall Trough"] <- 0.00129
mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Bay of Biscay BISB07"] <- 0.000523
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Bay of Biscay BISC07"] <- 0.000143
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Celtic Sea CEL08"] <- 0.000069
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "English Channel ENG08"] <- 0.000092
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Irish Sea"] <- 0.000095
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Kagerrak"] <- 0.000099
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Kattegat KATB07"] <- 0.000095
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea BEL08f"] <- 0.000071
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea BEL08j"] <- 0.000144
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea BEL08s"] <- 0.0001
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea GER07"] <- 0.000145
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea NOR08"] <- 0.000167
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea THA07j"] <- 0.000206
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "North Sea THA08"] <- 0.000082
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Scheldt estuary ZAN07"] <- 0.000141
  mtdna$Pise[mtdna$Source == "Cuveliers et al. 2012 Marine Biology 159:1239-1253" & mtdna$Site == "Wadden Sea TEX07"] <- 0.00009
mtdna$Pise[mtdna$Source == "Daley et al. 2012 Marine and Freshwater Research 63:708-714" & mtdna$Site == "North-eastern Tasmania"] <- 0.00834
  mtdna$Pise[mtdna$Source == "Daley et al. 2012 Marine and Freshwater Research 63:708-714" & mtdna$Site == "Western Australia"] <- 0.000447
  mtdna$Pise[mtdna$Source == "Daley et al. 2012 Marine and Freshwater Research 63:708-714" & mtdna$Site == "Victoria"] <- 0.000471
mtdna$Pi[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Broome"] <- 0.0031
  mtdna$Pise[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Broome"] <- 0.000633
  mtdna$Pise[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Great Barrier Reef"] <- 0.000816
  mtdna$Pise[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Bimini"] <- 0.000471
  mtdna$Pise[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Mariana Islands"] <- 0.000353 
  mtdna$Pise[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Seychelles"] <- 0.000111
  mtdna$Pise[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "San Blas"] <- 0.000291
  mtdna$Pise[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Hawaii"] <- 0.000213
  mtdna$Pise[mtdna$Source == "Daly-Engel et al. 2012 Marine Biology 159:975-985" & mtdna$Site == "Virgin Islands"] <- 0.000531
mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "American Samoa"] <- 0.000646
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Line Islands"] <- 0.000633
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Marshall Islands"] <- 0.000849
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Society Islands"] <- 0.000707
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "French Frigate Shoals"] <- 0.000108
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Gardner Pinnacles"] <- 0.000139
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Hawaii, Big Island"] <- 0.000166
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Johnston Atoll"] <- 0.000076
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Kauai"] <- 0.00008
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Kure Atoll"] <- 0.000088
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Laysan  Island"] <- 0.000109
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Lisianski Island"] <- 0.000133
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Maro Reef"] <- 0.000135
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Midway Island"] <- 0.000148
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Nihoa"] <- 0.000133
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Oahu"] <- 0.000136
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2011 Journal of Marine Biology 2011:839134" & mtdna$Site == "Pearl & Hermes Reef"] <- 0.00015
mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Christmas Island" & mtdna$spp == "Chaetodon meyeri"] <- 0.000269
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Cocos-Keeling Islands" & mtdna$spp == "Chaetodon meyeri"] <- 0.000367
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Diego Garcia" & mtdna$spp == "Chaetodon meyeri"] <- 0.0002
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Republic of Kiribati" & mtdna$spp == "Chaetodon meyeri"] <- 0.000347
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Republic of Palau" & mtdna$spp == "Chaetodon meyeri"] <- 0.000416
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Republic of Seychelles" & mtdna$spp == "Chaetodon meyeri"] <- 0.000694
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Moorea" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.000383
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Christmas Island" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.000513
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Cocos-Keeling Islands" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.000477
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Republic of Kiribati" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.000283
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Nuku Hiva" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.000329
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Palmyra Atoll" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.000695
  mtdna$Pise[mtdna$Source == "DiBattista et al. 2012 Journal of Heredity 103:617-629" & mtdna$Site == "Republic of Palau" & mtdna$spp == "Chaetodon ornatissimus"] <- 0.000411
mtdna$Pise[mtdna$Source == "Duncan et al. 2006" & mtdna$Site == "Pac Panama"] <- 0.000566
  mtdna$Pise[mtdna$Source == "Duncan et al. 2006" & mtdna$Site == "Seychelles"] <- 0.000318
  mtdna$Pise[mtdna$Source == "Duncan et al. 2006" & mtdna$Site == "South Africa"] <- 0.00058
  mtdna$Pise[mtdna$Source == "Duncan et al. 2006" & mtdna$Site == "Hawaii"] <- 0.000045
mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Cocos Islands"] <- 0.000336
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Moorea"] <- 0.000341
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Christmas Island"] <- 0.000224
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Diego Garcia"] <- 0.000301
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Seychelles"] <- 0.000212
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Kiritimati"] <- 0.000338
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Fiji"] <- 0.000365
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "French Frigate Shoals"] <- 0.000226
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Gardner Pinnacles"] <- 0.00028
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Hawaii Island"] <- 0.000261
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Kauai"] <- 0.000294
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Maui"] <- 0.0003
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Molokai"] <- 0.000297
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Nihoa"] <- 0.000253
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Oahu"] <- 0.000208
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Pearl & Hermes Reef"] <- 0.000604
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Journal of Marine Biology 2011:518516" & mtdna$Site == "Palau"] <- 0.0003
mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Chichi-jima"] <- 0.000354
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "French Frigate Atoll"] <- 0.000316
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Hilo"] <- 0.000309
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Johnston Atoll"] <- 0.000354 
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Kauai"] <- 0.000309
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Kealakekua"] <- 0.0005
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Kure Atoll"] <- 0.000338
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Lanai"] <- 0.000338
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Maro Reef"] <- 0.000392
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Maui"] <- 0.00032
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Midway Atoll"] <- 0.000348
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Molokai"] <- 0.000433
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Nihoa"] <- 0.000447
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Oahu"] <- 0.000338
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Peral & Hermes Atoll"] <- 0.000371
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Pohnpei"] <- 0.000371
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Saipan"] <- 0.000348
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "South Point"] <- 0.000329
  mtdna$Pise[mtdna$Source == "Eble et al. 2011 Marine Ecology Progress Series 428:245-258" & mtdna$Site == "Wawaloli"] <- 0.000359
mtdna$Pise[mtdna$Source == "Farnsworth et al. 2010 Marine Biology 157:945-953 "] <- NA
mtdna$Pise[mtdna$Source == "Froukh Koschzius 2007" & mtdna$Site == "Tor"] <- 0.000577
  mtdna$Pise[mtdna$Source == "Froukh Koschzius 2007" & mtdna$Site == "Aqaba"] <- 0.000374
  mtdna$Pise[mtdna$Source == "Froukh Koschzius 2007" & mtdna$Site == "Lith"] <- 0.000769
  mtdna$Pise[mtdna$Source == "Froukh Koschzius 2007" & mtdna$Site == "Thuwal"] <- 0.000516
  mtdna$Pise[mtdna$Source == "Froukh Koschzius 2007" & mtdna$Site == "Hodeidah"] <- 0.000652
mtdna$Pise[mtdna$Source == "Gaither et al. 2010 Journal of Biogeography 37:133-147 "] <- NA
  mtdna$Pise[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "French Frigate Shoals"] <- 0.00269
  mtdna$Pise[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Hilo"] <- 0.00252
  mtdna$Pise[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Kauai"] <- 0.003
  mtdna$Pise[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Kona"] <- 0.00269
  mtdna$Pise[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Kure"] <- 0.00633
  mtdna$Pise[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Maro"] <- 0.00393
  mtdna$Pise[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Maui Nui"] <- 0.00272
  mtdna$Pise[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Midway"] <- 0.00364
  mtdna$Pise[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Necker"] <- 0.00257
  mtdna$Pise[mtdna$Source == "Gaither et al. 2010 Molecular Ecology 19:1107-1121" & mtdna$Site == "Oahu"] <- 0.0024
mtdna$Pise[mtdna$Source == "Gaither et al. 2011 BMC Evolutionary Biology 11:1-16"] <- NA
mtdna$Pise[mtdna$Source == "Gnaither et al. 2011 PLoS ONE 6:1-13"] <- NA
mtdna$Pise[mtdna$Source == "Gaither et al. 2012 Proceedings of the Royal Society B 279:3948-3957"] <- NA
mtdna$Pise[mtdna$Source == "Garcia-Rodriguez et al. 2011 Fisheries Research 107:169-176" & mtdna$Site == "Bahia Magdalena" & mtdna$CollectionYear == 2006] <- 0.00103
  mtdna$Pise[mtdna$Source == "Garcia-Rodriguez et al. 2011 Fisheries Research 107:169-176" & mtdna$Site == "Bahia Magdalena" & mtdna$CollectionYear == 2007] <- 0.00165
  mtdna$Pise[mtdna$Source == "Garcia-Rodriguez et al. 2011 Fisheries Research 107:169-176" & mtdna$Site == "Ensenada "] <- 0.00142
mtdna$Pi[mtdna$Source == "Gotoh et. al 2009 Genes Genetic Systems 84:287-295" & mtdna$Site == "CLM"] <- 0.00008
  mtdna$Pi[mtdna$Source == "Gotoh et. al 2009 Genes Genetic Systems 84:287-295" & mtdna$Site == "CMR"] <- 0.00198
  mtdna$Pi[mtdna$Source == "Gotoh et. al 2009 Genes Genetic Systems 84:287-295" & mtdna$Site == "JFL"] <- 0.00063
  mtdna$Pi[mtdna$Source == "Gotoh et. al 2009 Genes Genetic Systems 84:287-295" & mtdna$Site == "JFOS"] <- 0.00103
  mtdna$Pi[mtdna$Source == "Gotoh et. al 2009 Genes Genetic Systems 84:287-295" & mtdna$Site == "NLK"] <- 0.00113
  mtdna$Pi[mtdna$Source == "Gotoh et. al 2009 Genes Genetic Systems 84:287-295" & mtdna$Site == "PPE"] <- 0.00147
  mtdna$Pise[mtdna$Source == "Gotoh et. al 2009 Genes Genetic Systems 84:287-295"] <- NA
mtdna$Pise[mtdna$Source == "Hempelmann et al. 2004" & mtdna$Site == "East side of San Juan Island, near Friday Harbor Laboratories"] <- 0.000378
  mtdna$Pise[mtdna$Source == "Hempelmann et al. 2004" & mtdna$Site == "North of Shaw Island, Broken Point"] <- 0.0004
  mtdna$Pise[mtdna$Source == "Hempelmann et al. 2004" & mtdna$Site == "Off Boeing Creek, Seattle"] <- 0.000626
  mtdna$Pise[mtdna$Source == "Hempelmann et al. 2004" & mtdna$Site == "Point George, South of Shaw Island"] <- 0.00053 
  mtdna$Pise[mtdna$Source == "Hempelmann et al. 2004" & mtdna$Site == "West side of San Juan Island, off Eagle Cove"] <- 0.000459
mtdna$Pise[mtdna$Source == "Itoi et al. 2011 Fisheries Science 77: 975-981" & mtdna$Site == "Iwate Prefecture"] <- 0.000247
  mtdna$Pise[mtdna$Source == "Itoi et al. 2011 Fisheries Science 77: 975-981" & mtdna$Site == "Izu Islands"] <- 0.000345
mtdna$Pise[mtdna$Source == "Iwamoto et al. 2012 Fisheries Science 78: 249-257" & mtdna$Site == "Pingtung, Taiwan"] <- 0.00126
  mtdna$Pise[mtdna$Source == "Iwamoto et al. 2012 Fisheries Science 78: 249-257" & mtdna$Site == "Ishigaki"] <- 0.00151
  mtdna$Pise[mtdna$Source == "Iwamoto et al. 2012 Fisheries Science 78: 249-257" & mtdna$Site == "Miyako"] <- 0.0017
  mtdna$Pise[mtdna$Source == "Iwamoto et al. 2012 Fisheries Science 78: 249-257" & mtdna$Site == "Okinawa"] <- 0.0012
  mtdna$Pise[mtdna$Source == "Iwamoto et al. 2012 Fisheries Science 78: 249-257" & mtdna$Site == "Cebu"] <- 0.00125
mtdna$Pise[mtdna$Source == "Johnson et al. 2009 Copeia 3: 465-474" & mtdna$Site == "Crown Memorial State Beach, CA (SF)"] <- 0.000218
  mtdna$Pise[mtdna$Source == "Johnson et al. 2009 Copeia 3: 465-474" & mtdna$Site == "Doheny State Beach, CA (OC)"] <- 0.000289
  mtdna$Pise[mtdna$Source == "Johnson et al. 2009 Copeia 3: 465-474" & mtdna$Site == "Malibu Lagoon State Beach, CA (LA)"] <- 0.000243
mtdna$Pise[mtdna$Source == "Karl et al. 2012 Marine Biology 159:489-498" & mtdna$Site == "BMN"] <- 0.000087
  mtdna$Pise[mtdna$Source == "Karl et al. 2012 Marine Biology 159:489-498" & mtdna$Site == "BLZ"] <- 0.00007
  mtdna$Pise[mtdna$Source == "Karl et al. 2012 Marine Biology 159:489-498" & mtdna$Site == "BRI"] <- 0.000066
  mtdna$Pise[mtdna$Source == "Karl et al. 2012 Marine Biology 159:489-498" & mtdna$Site == "NTL"] <- 0.000013
  mtdna$Pise[mtdna$Source == "Karl et al. 2012 Marine Biology 159:489-498" & mtdna$Site == "FLD"] <- 0.000033
mtdna$Pise[mtdna$Source == "Kim Park Kim 2006" & mtdna$Site == "Jumunjin, East Sea" & mtdna$MarkerName == "control region"] <- 0.00267
  mtdna$Pise[mtdna$Source == "Kim Park Kim 2006" & mtdna$Site == "Sangju, South Sea" & mtdna$MarkerName == "control region"] <- 0.00106
  mtdna$Pise[mtdna$Source == "Kim Park Kim 2006" & mtdna$Site == "Taean, West Sea" & mtdna$MarkerName == "control region"] <- 0.00106
mtdna$Pise[mtdna$Source == "Kochzius Blohm 2005" & mtdna$Site == "Northern Red Sea"] <- 0.00177
  mtdna$Pise[mtdna$Source == "Kochzius Blohm 2005" & mtdna$Site == "Gulf of Aqaba"] <- 0.00155
mtdna$Pise[mtdna$Source == "Kojima et. al 2009 Journal of Oceanography 65:187-193"] <- NA
mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "East China Sea, E8-3"] <- 0.02
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "F4"] <- 0.018
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "J1"] <- 0.018
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "Yellow Sea, L3"] <- 0.017
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "Yellow Sea, PT4"] <- 0.015
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "Anamizu, Ishikawa"] <- 0.018
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "Kagoshima Bay"] <- 0.018
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "Kosasa, Nagasaki"] <- 0.019
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "Miura Kanagawa"] <- 0.016
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "Mugi, Tokushima"] <- 0.02
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "Mutsu Bay"] <- 0.018
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "Teuri Island"] <- 0.018
  mtdna$Pi[mtdna$Source == "Liu et al. 2006" & mtdna$Site == "Utatsu, Miyagi"] <- 0.02
mtdna$Pi[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Aleutian Island"] <- 0.00099
  mtdna$Pi[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Eastern Hokkaido"] <- 0.00069
  mtdna$Pi[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Gulf of Georgia"] <- 0.00022
  mtdna$Pi[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Hecate Strait"] <- 0.00074
  mtdna$Pi[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Kodiak Island"] <- 0.00081
  mtdna$Pi[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Okhotsk Sea"] <- 0.0004
  mtdna$Pi[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Sea of Japan"] <- 0.0012
  mtdna$Pi[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Unimak Pass"] <- 0.00123
  mtdna$Pi[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Yellow Sea"] <- 0.00065
  mtdna$Pise[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Aleutian Island"] <- 0.00102
  mtdna$Pise[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Eastern Hokkaido"] <- 0.0008
  mtdna$Pise[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Gulf of Georgia"] <- 0.00044
  mtdna$Pise[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Hecate Strait"] <- 0.00085
  mtdna$Pise[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Kodiak Island"] <- 0.0009
  mtdna$Pise[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Okhotsk Sea"] <- 0.00059
  mtdna$Pise[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Sea of Japan"] <- 0.00114
  mtdna$Pise[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Unimak Pass"] <- 0.00118
  mtdna$Pise[mtdna$Source == "Liu et al. 2010 Journal of Fish Biology 77:1071-1082 " & mtdna$Site == "Yellow Sea"] <- 0.00078
mtdna$Pise[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "MA"] <- 0.0012
  mtdna$Pise[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "SK"] <- 0.0015
  mtdna$Pise[mtdna$Source == "Liu et al. 2010 Marine and Freshwater Research 61:1416-1424" & mtdna$Site == "YL"] <- 0.0014
mtdna$Pise[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "AU"] <- 0.00155
  mtdna$Pise[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "CE"] <- 0.00022
  mtdna$Pise[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "CH"] <- 0.00081
  mtdna$Pise[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "FI"] <- 0.0117
  mtdna$Pise[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "KA"] <- 0.0111
  mtdna$Pise[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "KO"] <- 0.00471
  mtdna$Pise[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "KOR"] <- 0.000695
  mtdna$Pise[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "KW"] <- 0.00471
  mtdna$Pise[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "MO"] <- 0.000869
  mtdna$Pise[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "TW"] <- 0.00325
  mtdna$Pise[mtdna$Source == "Liu et al. 2012 Marine Ecology Progress Series 458:155-167" & mtdna$Site == "WA"] <- 0.00318
mtdna$Pise[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "CH"] <- 0.0012
  mtdna$Pise[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "MICH"] <- 0.00086
  mtdna$Pise[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "OAX"] <- 0.00102
  mtdna$Pise[mtdna$Source == "Lopez et al. 2010 BMC Genetics 11:34" & mtdna$Site == "SIN"] <- 0.00148
mtdna$Pise[mtdna$Source == "McDowell et al. 2007 Gulf and Caribbean Research 19:75-82" & mtdna$Site == "Eastern Atlantic"] <-  0.0144 
  mtdna$Pise[mtdna$Source == "McDowell et al. 2007 Gulf and Caribbean Research 19:75-82" & mtdna$Site == "Caribeean Sea"] <- 0.0205
  mtdna$Pise[mtdna$Source == "McDowell et al. 2007 Gulf and Caribbean Research 19:75-82" & mtdna$Site == "Western North Atlantic, US mid-Atlantic"] <- 0.0152
mtdna$Pi[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Europa "] <- 0.0031
  mtdna$Pi[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Geyser "] <- 0.0042
  mtdna$Pi[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Glorieuses "] <- 0.0032
  mtdna$Pi[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Juan de Nova "] <- 0.0036
  mtdna$Pi[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Kenya"] <- 0.0045
  mtdna$Pi[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Madagascar"] <- 0.0027
  mtdna$Pi[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Mayotte "] <- 0.0036
  mtdna$Pi[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Moheli "] <- 0.003
  mtdna$Pi[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Moroni "] <- 0.0039
  mtdna$Pi[mtdna$Source == "Muths et al. 2011 Marine Ecology Progress Series 443: 167-180" & mtdna$Site == "Reunion "] <- 0.0028
mtdna$Pise[mtdna$Source == "Rocha et. al 2008 BMC Evolutionary Biology 8:157"] <- NA
mtdna$Pise[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "California"] <- 0.00112
  mtdna$Pise[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Fairweather Island, Gulf of Alaska"] <- 0.00136 
  mtdna$Pise[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Oregon"] <- 0.00129
  mtdna$Pise[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Sitka, Gulf of Alaska"] <- 0.00091
  mtdna$Pise[mtdna$Source == "Rocha-Olivares Vetter 1999" & mtdna$Site == "Vancouver Island, British Columbia"] <- 0.00197
mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Argentina"] <- 0.000452
  mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Amapa"] <- 0.000411
  mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Bahia"] <- 0.000653
  mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Espirito Santo"] <- 0.000194 
  mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Maranhao"] <- 0.00031
  mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Para"] <- 0.000452
  mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Parana"] <- 0.000362
  mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Pernambuco"] <- 0.000632
  mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Rio Grande do Sul"] <- 0.00053
  mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Santa Catarina"] <- 0.000567
  mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Sao Paulo"] <- 0.00052
  mtdna$Pise[mtdna$Source == "Santos et al. 2006" & mtdna$Site == "Venezuela"] <- 0.000474
mtdna$Pise[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "FKS" & mtdna$MarkerName == "control region"] <- 0.0012
  mtdna$Pise[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "SRK" & mtdna$MarkerName == "control region"] <- 0.0012
  mtdna$Pise[mtdna$Source == "Sekino et. al 2011 Conservation Genetics 12:139-159" & mtdna$Site == "TBB" & mtdna$MarkerName == "control region"] <- 0.0001
mtdna$Pise[mtdna$Source == "Sotka et al. 2005 Marine Biotechnology 7:223-230" & mtdna$Site == "Friday Harbor"] <- 0.000378
  mtdna$Pise[mtdna$Source == "Sotka et al. 2005 Marine Biotechnology 7:223-230" & mtdna$Site == "Pt. George"] <- 0.00053
  mtdna$Pise[mtdna$Source == "Sotka et al. 2005 Marine Biotechnology 7:223-230" & mtdna$Site == "Seattle"] <- 0.000417
  mtdna$Pise[mtdna$Source == "Sotka et al. 2005 Marine Biotechnology 7:223-230" & mtdna$Site == "Shaw Island"] <- 0.0004
  mtdna$Pise[mtdna$Source == "Sotka et al. 2005 Marine Biotechnology 7:223-230" & mtdna$Site == "Westside"] <- 0.000459
mtdna$Pise[mtdna$Source == "Teske et al. 2005" & mtdna$Site == "Pulai Estuary, Johor, Peninsular Malaysia"] <- 0.000068
  mtdna$Pise[mtdna$Source == "Teske et al. 2005" & mtdna$Site == "Tayabas Bay, Quezon"] <- 0.000068
  mtdna$Pise[mtdna$Source == "Teske et al. 2005" & mtdna$Site == "Tamil Nadu"] <- 0.000122
mtdna$Pi[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "Kenya"] <- 0.049
  mtdna$Pi[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "Mauritius"] <- 0.048
  mtdna$Pi[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "KSeychelles"] <- 0.048
  mtdna$Pise[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "Kenya"] <- 0.000626
  mtdna$Pise[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "Mauritius"] <- 0.000359
  mtdna$Pise[mtdna$Source == "Visram et. al 2010 Marine Biology 157:1475-1487" & mtdna$Site == "KSeychelles"] <- 0.000343
mtdna$Pise[mtdna$Source == "Wilson 2006" & mtdna$Site == "Newport"] <- 0.000492
  mtdna$Pise[mtdna$Source == "Wilson 2006" & mtdna$Site == "Padilla Bay"] <- 0.000291
  mtdna$Pise[mtdna$Source == "Wilson 2006" & mtdna$Site == "San Diego"] <- 0.000783
  mtdna$Pise[mtdna$Source == "Wilson 2006" & mtdna$Site == "Sitka, Alaska"] <- 0.000179
mtdna$Pi[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q1 - center-south of Yellow Sea" & mtdna$MarkerName == "cytb"] <- 0.00744
  mtdna$Pi[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q1 - center-south of Yellow Sea" & mtdna$MarkerName == "cytochrome c oxidase subunit I"] <- 0.00519
  mtdna$Pi[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q2 - northeast East China Sea" & mtdna$MarkerName == "cytb"] <- 0.00589
  mtdna$Pi[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q2 - northeast East China Sea" & mtdna$MarkerName == "cytochrome c oxidase subunit I"] <- 0.00532
  mtdna$Pi[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q3 - southwest East China Sea" & mtdna$MarkerName == "cytochrome c oxidase subunit I"] <- 0.00832
  mtdna$Pi[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q3 - southwest East China Sea" & mtdna$MarkerName == "cytb"] <- 0.00603
  mtdna$Pise[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q1 - center-south of Yellow Sea" & mtdna$MarkerName == "cytb"] <- 0.00458
  mtdna$Pise[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q1 - center-south of Yellow Sea" & mtdna$MarkerName == "cytochrome c oxidase subunit I"] <- 0.00320
  mtdna$Pise[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q2 - northeast East China Sea" & mtdna$MarkerName == "cytb"] <- 0.00376
  mtdna$Pise[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q2 - northeast East China Sea" & mtdna$MarkerName == "cytochrome c oxidase subunit I"] <- 0.00325
  mtdna$Pise[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q3 - southwest East China Sea" & mtdna$MarkerName == "cytochrome c oxidase subunit I"] <- 0.00479
  mtdna$Pise[mtdna$Source == "Yu et al. 2005" & mtdna$Site == "Q3 - southwest East China Sea" & mtdna$MarkerName == "cytb"] <- 0.00385
mtdna$Pise[mtdna$Source == "Zane et al. 2006" & mtdna$Site == "Elephant Island, South Shetland Islands"] <- 0.00131
  mtdna$Pise[mtdna$Source == "Zane et al. 2006" & mtdna$Site == "Halley Bay, Weddell Sea"] <- 0.00067
  mtdna$Pise[mtdna$Source == "Zane et al. 2006" & mtdna$Site == "King George Island, Southern Scotia Arc"] <- 0.00071
  mtdna$Pise[mtdna$Source == "Zane et al. 2006" & mtdna$Site == "Terra Nova Bay, Ross Sea"] <- 0.00073

#note: 2 instances where He is non-zero but Pi is zero

##########################################################################################################################################
  
######## Trim & write out data ########
  
#check to make sure all data included
inds <- (is.na(mtdna$He) & is.na(mtdna$Pi) | is.na(mtdna$spp) | is.na(mtdna$lat) | is.na(mtdna$lon) | is.na(mtdna$n))
sum(inds) #0
mtdna[inds,]

#trim to just relevant columns
mtdna <- mtdna[ ,c('spp', 'CommonName', 'Source', 'Country', 'Site', 'lat', 'lon', 'stockid', 'CollectionYear', 
                   'MarkerName', 'n', 'bp', 'He', 'Hese', 'Pi', 'Pise', 'file')]
dim(mtdna) #2052 x 17

#write out mtdna data (allow multiple loci per line)
write.csv(mtdna, file='output/mtdna_assembled.csv')
################################################### Script for Assembling msat Data  #######################################################

#assemble the msat data and do basic QA/QC
#starts from the new individual msat files as well as old files Malin & undergrads put together from earlier work (prior to 2011)
#skipping matching with fishing data for now as don't have access to those datasets

##########################################################################################################################################

######## Set-up ########

remove(list = ls())

#load libraries
library(tidyverse)

#write basic functions
calcHe <- function(x) { # calculate expected he from allele frequencies (does NOT check that freqs sum to 1)
	x <- x[!is.na(x)]
	h <- 1 - sum(x^2)
	return(h)
}

sumna <- function(x){ # remove NAs, unless all are NAs, then return NA
	if(all(is.na(x)))
	  return(NA)
	else
	  return(sum(x, na.rm=TRUE))
}

#read in data
msat1 <- read.csv("data/msat/msat_2011-2020_data_301.csv", stringsAsFactors=FALSE)
  msat1 <- msat1[, !(names(msat1) == 'Notes')] #remove notes column from df
msat2 <- read.csv("data/msat/msat_2011-2020_data_2_302.csv", stringsAsFactors = FALSE)
  msat2 <- msat2[, !(names(msat2) == 'Notes')]
msat3 <- read.csv("data/msat/msat_2011-2020_data_3_303.csv", stringsAsFactors = FALSE)
  msat3 <- msat3[, !(names(msat3) == 'Notes')]
msat4 <- read.csv("data/msat/msat_2011-2020_data_4_304.csv", stringsAsFactors = FALSE)
  msat4 <- msat4[, !(names(msat4) == 'Notes')]
msat5 <- read.csv("data/msat/msat_2011-2020_data_marial.csv", stringsAsFactors = FALSE) #csv with Marial's data
  msat5 <- msat5[, !names(msat5) == 'Notes']
oldmsat1 <- read.csv("data/msat/Fishery lat msats 000 2015-05-20 SSP.csv", stringsAsFactors = FALSE) #csv from previous msat assembly
oldmsat2 <- read.csv("data/msat/Fishery lat msats 001 2015-07-04 MLP.csv", stringsAsFactors = FALSE)
oldmsat3 <- read.csv("data/msat/Fishery lat msats 002 2015-08-08 SSP.csv", stringsAsFactors = FALSE)
  oldmsat3 <- oldmsat3[oldmsat3$spp != '', ] #trim empty rows from df
oldmsat4 <- read.csv("data/msat/Fishery lat msats 100 2015-08-20 MLP.csv", stringsAsFactors = FALSE)
  oldmsat4 <- oldmsat4[, (names(oldmsat4) != 'X')] #drop extra column from df
  oldmsat4 <- oldmsat4[oldmsat4$spp != '', ] #trim empty rows from df
oldmsat5 <- read.csv("data/msat/Fishery lat msats 101 2015-07-17 SSP.csv", stringsAsFactors = FALSE)
  oldmsat5 <- oldmsat5[, (names(oldmsat5) != 'X')] #drop extra column from df
  oldmsat5 <- oldmsat5[oldmsat5$spp != '', ] #trim empty rows from df
oldmsat6 <- read.csv("data/msat/Fishery lat msats 200 2015-02-10 MLP.csv", stringsAsFactors = FALSE)
  oldmsat6 <- oldmsat6[, (names(oldmsat6) != 'X')] #drop extra column from df
  oldmsat6 <- oldmsat6[oldmsat6$spp != '', ] #trim empty rows from df
oldmsat7 <- read.csv("data/msat/Fishery lat msats 201 2015-10-13 MLP.csv", stringsAsFactors = FALSE)
ppdat <- read.csv("data/msat/ppdat_2016-03-04wLL.csv", stringsAsFactors=FALSE) #Pinsky & Palumbi 2014 data
  names(ppdat)[names(ppdat) == "NumMarker"] <- "NumMarkers" #for some reason, the contents weren't copied over before
srdbmatch <- read.csv("data/srdb_matching/msat_to_match.csv", stringsAsFactors = FALSE) #to match genetic data to SRDB stocks
  
##########################################################################################################################################

######## Check for duplicates ########

#add an indicator for source file
oldmsat1$file <- 'msats000'
oldmsat2$file <- 'msats001'
oldmsat3$file <- 'msats002'
oldmsat4$file <- 'msats100'
oldmsat5$file <- 'msats101'
oldmsat6$file <- 'msats200'
oldmsat7$file <- 'msats201'
ppdat$file <- 'ppdat'
msat1$file <- 'msats301'
msat2$file <- 'msats302'
msat3$file <- 'msats303'
msat4$file <- 'msats304'
msat5$file <- 'msats305'

#check for duplicate studies between oldmsat studies & ppdat
overlaps <- intersect(c(oldmsat1$spp, oldmsat2$spp, oldmsat3$spp, oldmsat4$spp, 
                        oldmsat5$spp, oldmsat6$spp, oldmsat7$spp), ppdat$spp) #species included in all old datasets
ppdat_possdup_studies <- sort(unique(ppdat$Source[ppdat$spp %in% overlaps])) #list of studies from ppdat on these species
oldmsat_possdup_studies <- sort(unique(c(oldmsat1$Source[oldmsat1$spp %in% overlaps], oldmsat2$Source[oldmsat2$spp %in% overlaps], 
                                         oldmsat3$Source[oldmsat3$spp %in% overlaps], oldmsat4$Source[oldmsat4$spp %in% overlaps], 
                                         oldmsat5$Source[oldmsat5$spp %in% overlaps], oldmsat6$Source[oldmsat6$spp %in% overlaps], 
                                         oldmsat7$Source[oldmsat7 %in% overlaps]))) #list of studies from oldmsat dataframes on these species
dups <- c("Bekkevold et al. 2005 Evolution 59:2656-2668", "Canino et al. 2005 Mol Ecol Notes 5:908-910", 
          "Carlsson et al. 2004. Mol. Ecol.13:3345-3356.", "Roques et al. 2007 Mol. Ecol. Notes 7:661-663", 
          "Selkoe et al. 2005 Mol Ecol Notes", "Wang et al. 2010 Gen & Mol Res 9:931-934", 
          "Yagishita & Kobayashi 2008 Mol. Ecol. Res. 8:302-304") #studies in oldmsat dfs that appear to have been previously recorded in ppdat (except Selkoe --> that pulled out of ppdat bc oldmsat has lat/long coordinates)

#remove studies from oldmsat dataframes (and Selkoe from ppdat)
oldmsat1 <- oldmsat1[!(oldmsat1$Source %in% dups), ]
oldmsat2 <- oldmsat2[!(oldmsat2$Source %in% dups), ]
oldmsat3 <- oldmsat3[!(oldmsat3$Source %in% dups), ]
oldmsat4 <- oldmsat4[!(oldmsat4$Source %in% dups), ]
oldmsat5 <- oldmsat5[!(oldmsat5$Source %in% dups), ]
oldmsat6 <- oldmsat6[!(oldmsat6$Source %in% dups), ]
oldmsat7 <- oldmsat7[!(oldmsat7$Source %in% dups), ]
ppdat <- ppdat[!(ppdat$Source %in% dups), ]

#check for duplicate studies between new msat and oldmsat studies
#NOTE: some dups will be okay (come from double-checking old studies, and adding new sites/markers/species when previously missed)
overlaps <- intersect(c(msat1$spp, msat2$spp, msat3$spp, msat4$spp, msat5$spp), 
                      c(oldmsat1$spp, oldmsat2$spp, oldmsat3$spp, oldmsat4$spp, oldmsat5$spp, 
                        oldmsat6$spp, oldmsat7$spp, ppdat)) #species included in old and new datasets
oldmsat_possdup_studies <- sort(unique(c(oldmsat1$Source[oldmsat1$spp %in% overlaps], oldmsat2$Source[oldmsat2$spp %in% overlaps], 
                                         oldmsat3$Source[oldmsat3$spp %in% overlaps], oldmsat4$Source[oldmsat4$spp %in% overlaps], 
                                         oldmsat5$Source[oldmsat5$spp %in% overlaps], oldmsat6$Source[oldmsat6$spp %in% overlaps], 
                                         oldmsat7$Source[oldmsat7 %in% overlaps], ppdat$Source[ppdat$spp %in% overlaps]))) #list of studies from oldmsat dataframes on these species
msat_possdup_studies <- sort(unique(c(msat1$Source[msat1$spp %in% overlaps], msat2$Source[msat2$spp %in% overlaps], 
                                      msat3$Source[msat3$spp %in% overlaps], msat4$Source[msat4$spp %in% overlaps], 
                                      msat5$Source[msat5$spp %in% overlaps]))) #list of studies from new msat dataframes on these species
dups <- c("D'Anatro et al. 2011 Env Biol Fish", "Dammannagoda et al. 2011 CJFAS 68:210-223", "Helyar et al. 2011 Cons Gen Res 3:173-176", 
          "Lawton et al. 2011 Mol Ecol 20:3584-3598", "Machado-Schiaffino et al. 2011 Biol Cons 144:330-338", 
          "Muths & Bourjea 2011 Cons Gen Res 3:629-631", "Nance et al. 2011 PLoS ONE 6:e21459", "Palof et al. 2011 Mar Biol 158:779-792", 
          "Palof et al. 2011 Mar Biol 158:779-792", "Qiu & Miyamoto 2011 Copeia 2:264-269", "Schunter et al. 2011 J Exp Mar Biol & Ecol 401:126-133", 
          "Sekino et al. 2011 Cons Gen 12:139-159", "White et al. 2010 Heredity 106:690-699") #studies in oldmsat dataframes that were re-recorded in new msat dataframes (again, ignoring studies purposefully added for additional sites)

#remove studies from oldmsat dataframes
oldmsat1 <- oldmsat1[!(oldmsat1$Source %in% dups), ]
oldmsat2 <- oldmsat2[!(oldmsat2$Source %in% dups), ]
oldmsat3 <- oldmsat3[!(oldmsat3$Source %in% dups), ]
oldmsat4 <- oldmsat4[!(oldmsat4$Source %in% dups), ]
oldmsat5 <- oldmsat5[!(oldmsat5$Source %in% dups), ]
oldmsat6 <- oldmsat6[!(oldmsat6$Source %in% dups), ]
oldmsat7 <- oldmsat7[!(oldmsat7$Source %in% dups), ]
ppdat <- ppdat[!(ppdat$Source %in% dups), ]

#fuzzy matching to find other potential duplicates
studs <- vector('list', length = 12)
studs[[1]] <- sort(unique(oldmsat1$Source))
studs[[2]] <- sort(unique(oldmsat2$Source))
studs[[3]] <- sort(unique(oldmsat3$Source))
studs[[4]] <- sort(unique(oldmsat4$Source))
studs[[5]] <- sort(unique(oldmsat5$Source))
studs[[6]] <- sort(unique(oldmsat6$Source))
studs[[7]] <- sort(unique(oldmsat7$Source))
studs[[8]] <- sort(unique(ppdat$Source))
studs[[9]] <- sort(unique(msat1$Source))
studs[[10]] <- sort(unique(msat2$Source))
studs[[11]] <- sort(unique(msat3$Source))
studs[[12]] <- sort(unique(msat4$Source))
studs[[13]] <- sort(unique(msat5$Source))

poss_matches <- vector('list', length = length(studs))
for(j in 1:(length(studs) - 1)) {
  poss_matches[[j]] <- vector('list', length = length(studs) - j)
  for(i in (j + 1):length(studs)) {
    poss_matches[[j]][[i - j]] <- lapply(studs[[j]], FUN = agrep, x = studs[[i]], fixed = TRUE, value = TRUE, max.distance = list(all = 0.2), useBytes = TRUE)
    lens <- sapply(poss_matches[[j]][[i - j]], FUN = length)
    if(any(lens > 0)) {
      print(paste('j=', j, ', i=', i, sep = ''))
      for(k in 1:sum(lens > 0)) {
        print(cbind(studs[[j]][which(lens > 0)[k]], paste(poss_matches[[j]][[i - j]][[which(lens > 0)[k]]], collapse = ',')))
      }
    }
  }
}

#examine and trim out matching studies found by fuzzy matching
#Bagley & Geller 1998
oldmsat6[oldmsat6$Source == 'Bagley & Geller 1998 Molecular Ecology 7:1083-1090', !grepl('^p', names(oldmsat6))] #not printing allele freq columns
oldmsat7[oldmsat7$Source == 'Bagley&Geller 1998 Molecular Ecology 7:1083-1090', !grepl('^p', names(oldmsat7))] 
oldmsat7 <- oldmsat7[oldmsat7$Source !='Bagley&Geller 1998 Molecular Ecology 7:1083-1090',] #keep oldmsat6 version, since it appropriately dropped Ra11 for signs of a null allele

#Bagley et al. 1999
oldmsat6[oldmsat6$Source == 'Bagley et al. 1999 Marine Biology 134:609-620', !grepl('^p', names(oldmsat6))]
oldmsat7[oldmsat7$Source == 'Bagley et al. 1999 Marine Biology 134:609-620', !grepl('^p', names(oldmsat7))] 
oldmsat6 <- oldmsat6[oldmsat6$Source != 'Bagley et al. 1999 Marine Biology 134:609-620', ] #keep oldmsat7 version, since it appropriately dropped Ra1, Ra5, and Ra6 for signs of HWP departure

#Garcia de Leon et al. 1997 Molecular Ecology 6:51-62
oldmsat6[oldmsat6$Source == 'Garcia de Leon et al. 1997 Molecular Ecology 6:51-62', !grepl('^p', names(oldmsat6))]
oldmsat7[oldmsat7$Source == 'Garcia de Leon et al. 1997 Molecular Ecology 6:51-62', !grepl('^p', names(oldmsat7))] 
oldmsat7 <- oldmsat7[oldmsat7$Source !='Garcia de Leon et al. 1997 Molecular Ecology 6:51-62', ] #keep oldmsat6 version arbitrarily, since they look the same

#Naciri et al. 1999 Journal of Heredity 90:592-596		
oldmsat6[oldmsat6$Source == 'Naciri et al. 1999 Journal of Heredity 90:592-596', !grepl('^p', names(oldmsat6))]
oldmsat7[oldmsat7$Source == 'Naciri et al. 1999 The Journal of Heredity 90: 591-596', !grepl('^p', names(oldmsat7))] 
oldmsat7 <- oldmsat7[oldmsat7$Source !='Naciri et al. 1999 The Journal of Heredity 90: 591-596', ] #keep oldmsat6 version since Malin was more careful about which sites to keep or drop

#examine and trim out matching studies found by accident
#Bahri-Sfar et al. 2000 Proc. R. Soc. Lond. B 267:929-935
oldmsat6[oldmsat6$Source == 'Bahri-Sfar et al. 2000 Proceedings of the Royal Society B 267:929-935', !grepl('^p', names(oldmsat6))]
oldmsat7[oldmsat7$Source == 'Bahri-Sfar et al. 2000 Proc. R. Soc. Lond. B 267:929-935', !grepl('^p', names(oldmsat7))] 
oldmsat7 <- oldmsat7[oldmsat7$Source !='Bahri-Sfar et al. 2000 Proc. R. Soc. Lond. B 267:929-935',] #keep msat6 version since Malin was more careful about which sites to keep or drop

######## Merge dataframes ########

#verify oldmsat column names match
all(names(oldmsat1) == names(oldmsat2))
all(names(oldmsat1) == names(oldmsat3))
all(names(oldmsat1) == names(oldmsat4))
all(names(oldmsat1) == names(oldmsat5))
all(names(oldmsat1) == names(oldmsat6))
all(names(oldmsat1) == names(oldmsat7)) #fails

#fix oldmsat7 column names
nn <- setdiff(names(oldmsat1), names(oldmsat7)) #check to see diff
for(i in 1:length(nn)) {
  oldmsat7[[nn[i]]] <- NA #pad extra columns in msat
}
setdiff(names(oldmsat1), names(oldmsat7)) #none: good
setdiff(names(oldmsat7), names(oldmsat1)) #none: good

#merge oldmsat dataframes together
oldmsat <- rbind(oldmsat1, oldmsat2, oldmsat3, oldmsat4, oldmsat5, oldmsat6, oldmsat7)
dim(oldmsat) # 3817 rows of data

#verify new msat column names match
all(names(msat1) == names(msat2))
all(names(msat1) == names(msat3))
all(names(msat1) == names(msat4))
all(names(msat1) == names(msat5)) #fails

#fix msat5 column names
nn <- setdiff(names(msat1), names(msat5)) #check to see diff
names(msat5)[names(msat5) == "NumMarker"] <- "NumMarkers" #change to make sure names match
all(names(msat1) == names(msat5)) #all match

#merge new msat dataframes together
msat <- rbind(msat1, msat2, msat3, msat4, msat5)
dim(msat) #11748 rows of data

#verify Pinsky & Palumbi column names
nms <- setdiff(names(oldmsat), names(ppdat)) #column names in oldmsat that aren't in ppdat
if(length(names) > 0) for (i in 1:length(nms)) {
  ppdat[[nms[i]]] <- NA #add these columns to ppdat
}
ppdatmsat <- ppdat[, names(oldmsat)] #trim P & P to only columns that match oldmsat
dim(ppdatmsat) #8204 x 89
dim(oldmsat) #3817 x 89

#merge Pinsky & Palumbi and oldmsat dataframes together
oldmsat <- rbind(oldmsat, ppdatmsat)
dim(oldmsat) #12021 x 89

#check old and new msat dataframes have same columns & col names
new_cols <- setdiff(names(oldmsat), names(msat)) #column names in oldmsat that aren't in msat
new_cols
if(length(names)>0) for(i in 1:length(new_cols)) {
  msat[[new_cols[i]]] <- NA #add these cols to msat
}
setdiff(names(oldmsat), names(msat)) #none: good
setdiff(names(msat), names(oldmsat)) #none: good
dim(oldmsat) #12021 x 89
dim(msat) #11748 x 89

#merge new and old msat dataframes together
msat <- rbind(msat, oldmsat)
dim(msat) #23769 x 89
		
######## Merge SRDB stock info ########

#match on spp, Source, & site
#read in fbsci, stockid, lat & lon (lat & lon for error-checking)
names(srdbmatch)[names(srdbmatch) == 'lat'] <- 'lat_srdb'
names(srdbmatch)[names(srdbmatch) == 'lon'] <- 'lon_srdb'

#trim to unique spp, Source, site
inds <- duplicated(srdbmatch[, c('spp', 'Source', 'Country', 'Site', 'CollectionYear')])
inds2 <- duplicated(srdbmatch[, c('spp', 'Source', 'Country', 'Site', 'CollectionYear', 'lat_srdb', 'lon_srdb')])
sum(inds) #6478
sum(inds2) #6478: matches which means lat/lon will be the same for each row kept (same number of duplicates whether include lat/lon or not)
nrow(srdbmatch) #7884
srdbmatch <- srdbmatch[!inds, ] #trimming to unique rows
nrow(srdbmatch) #1406

#merge SRDB stock info
#sort(setdiff(srdbmatch$Source, msat$Source)) #papers in srdbmatch that aren't in msat --> double-checked to make sure excluded bc He = NA
nrow(msat) #23769
msat <- merge(msat, srdbmatch[, c('spp', 'Source', 'Country', 'Site', 'CollectionYear', 'fbsci', 'stockid', 'lat_srdb', 'lon_srdb')], 
              all.x = TRUE, by = c('spp', 'Source', 'Country', 'Site', 'CollectionYear'))
nrow(msat) #23769

##########################################################################################################################################

######## Clean newly merged msat dataframe ########

#remove temporal replicates
dim(msat) #23769 rows
msat <- subset(msat, Source != "Larsson et al. 2010 Heredity 104:40-51" | Site != "Himmerfjarden" | CollectionYear != "1979-1980")
  msat <- subset(msat, Source != "Larsson et al. 2010 Heredity 104:40-51" | Site != "Kalix" | CollectionYear != "1980" & CollectionYear != "2002")
  msat <- subset(msat, Source != "Larsson et al. 2010 Heredity 104:40-51" | Site != "Vaxholm" | CollectionYear != "1979")
msat <- subset(msat, Source != "Pampoulie et al. 2006 Can. J. Fish. Aquat. Sci. 63:2660-2675" | Site != "311b" & Site != "511b" & Site != "812b")
msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Kattegat 02")
  msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Kolding Fjord 03")
  msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Rugen April 02")
  msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Rugen April 03")
  msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Rugen March 02")
  msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Rugen May 02")
  msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Flatbrotten 02")
  msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Maseskar 02")
  msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Tjome 02")
msat <- subset(msat, Source != "Chevolot et al. 2006 J Sea Research 56:305-316" | Site != "North Thames Estuary" | CollectionYear != 2003)
  msat <- subset(msat, Source != "Chevolot et al. 2006 J Sea Research 56:305-316" | Site != "Thames Estuary" | CollectionYear != 2004)
  msat <- subset(msat, Source != "Chevolot et al. 2006 J Sea Research 56:305-316" | Site != "East English Channel" | CollectionYear != 2004)
  msat <- subset(msat, Source != "Chevolot et al. 2006 J Sea Research 56:305-316" | Site != "Carmarthen Bay" | CollectionYear != 2003)
  msat <- subset(msat, Source != "Chevolot et al. 2006 J Sea Research 56:305-316" | Site != "Caernarfon Bay" | CollectionYear != 2004)
msat <- subset(msat, Source != "McCusker & Bentzen 2010 Mol Ecol 19:4228-4241" | Site != "Scotian Shelf" | CollectionYear != 2002)
  msat <- subset(msat, Source != "McCusker & Bentzen 2010 Mol Ecol 19:4228-4241" | Site != "Iceland" | CollectionYear != 2004)
  msat <- subset(msat, Source != "McCusker & Bentzen 2010 Mol Ecol 19:4228-4241" | Site != "Rockall Bank" | CollectionYear != 2005)
msat <- subset(msat, Source != "Saillant et al. 2010 ICES J Mar Sci 67:1240-1250" | Site != "Mississippi/Alabama" | CollectionYear != 2004)
  msat <- subset(msat, Source != "Saillant et al. 2010 ICES J Mar Sci 67:1240-1250" | Site != "Brownsville, TX" | CollectionYear != 2005)
  msat <- subset(msat, Source != "Saillant et al. 2010 ICES J Mar Sci 67:1240-1250" | Site != "Port Aransas, TX" | CollectionYear != 2005)
  msat <- subset(msat, Source != "Saillant et al. 2010 ICES J Mar Sci 67:1240-1250" | Site != "Freeport, TX" | CollectionYear != 2005)
  msat <- subset(msat, Source != "Saillant et al. 2010 ICES J Mar Sci 67:1240-1250" | Site != "Louisiana" | CollectionYear != 2005)
msat <- subset(msat, Source != "Tripp-Valdez et al. 2010 Fish Res 105:172-177" | Site != "Cabo San Lucas" | CollectionYear != 2005)
  msat <- subset(msat, Source != "Tripp-Valdez et al. 2010 Fish Res 105:172-177" | Site != "Loreto" | CollectionYear != 2006)
  msat <- subset(msat, Source != "Tripp-Valdez et al. 2010 Fish Res 105:172-177" | Site != "Guaymas" | CollectionYear != 2006)
msat <- subset(msat, Source != "Cunningham et al. 2009 Can. J. Fish. Aquat. Sci 66:153-166" | Site != "Unimak Pass" | CollectionYear != 2005)
  msat <- subset(msat, Source != "Cunningham et al. 2009 Can. J. Fish. Aquat. Sci 66:153-166" | Site != "Kodiak Island" | CollectionYear != 2003)
msat <- subset(msat, Source != "Appleyard, Williams & Ward. 2004. CCAMLR Science 11:21-32.") #all are temporal replicates of same site reported in msat304
msat <- subset(msat, Source != "D'Amato 2006 Mar Biotech" | Site != 1 & Site != 2)
msat <- subset(msat, Source != "Durand et al. 2005 Mar Bio" | Site != "Cape-2 (Off Cape of Good Hope)" & Site != "Cape-3 (Off Cape of Good Hope)" & Site != "Cape-4 (Off Cape of Good Hope)" & Site != "Cape-5 (Off Cape of Good Hope)")
msat <- subset(msat, Source != "Gomez-Uchida & Banks 2005 CJFAS" | Site != "O1C" & Site != "O1D" & Site != "O3B" & Site != "O3C" & Site != "O5A" & Site != "O5B")
msat <- subset(msat, Source != "Gomez-Uchida 2006 Dissertation" | Site != "4CA01" & Site != "4OR01" & Site != "8WA03" & Site != "CB73" & Site != "CB92" & Site != "NP86" & Site != "NP92")
msat <- subset(msat, Source != "Gonzalez & Zardoya 2007 BMC Evol Biol") #all are temporal replicates of same site reported in msat304
msat <- subset(msat, Source != "Gonzalez et al. 2008 BMC Evol" | Site != "Canada2" & Site != "Canary2" & Site != "Guinea2")
msat <- subset(msat, Source != "Hauser et al. 2007 Book" | Site != "Point Heyer, Vashon Island, Puget Sound, WA (Adults)")
msat <- subset(msat, Source != "Jorgensen et al. 2005 ICES J Mar Sci" | Site != "GD02-1" & Site != "GD02-2" & Site != "GD02-3" & Site != "GD03-2" & Site != "GD03-3" & Site != "RU02-1" & Site != "RU02-2" & Site != "RU02-3" & Site != "RU03-1" & Site != "RU03-3")
msat <- subset(msat, Source != "Limborg et al. 2009 MEPS" | Site != "Bornholm Basin (BOR06)" & Site != "German Bight (GER04)")
msat <- subset(msat, Source != "Lundy et al. 2000 Mol Ecol" | Site != "Bay of Biscay north 1999") #temporal replicate of same site reported in msat304
msat <- subset(msat, Source != "Mariani et al. 2005 MEPS" | Site != "Cape Wrath 02" & Site != "Flamborough 03" & Site != "Orkney/Shetland 03" & Site != "Whiten Head 03")
msat <- subset(msat, Source != "Mitchell 2006 Thesis" | Site != "Cherry Point, WA 1999 non-spawning")
msat <- subset(msat, Source != "O'Reilly et al. 2004 Mol Ecol" | Site != "Prince William Sound, AK A")
msat <- subset(msat, Source != "Pumitinsee et al. 2009 Aq Res" | Site != "Trang (April)" & Site != "Trang (January)" & Site != "Trang (November)")
msat <- subset(msat, Source != "Riccioni et al. 2010 PNAS" | Site != "Messina shore, Southern Tyrrhenian Sea" & Site != "North Adriatic Sea")
msat <- subset(msat, Source != "Ruzzante et al. 1996 inshore-offshore CJFAS" | Site != "North Cape, Grand Bank" | CollectionYear != 1993 & CollectionYear != 1994)
  msat <- subset(msat, Source != "Ruzzante et al. 1996 inshore-offshore CJFAS" | Site != "Smith Sound, Trinity Bay" | CollectionYear != 1993)
  msat <- subset(msat, Source != "Ruzzante et al. 1996 inshore-offshore CJFAS" | Site != "SW Arm, Trinity Bay" | CollectionYear != 1992 & CollectionYear != 1994)
  msat <- subset(msat, Source != "Ruzzante et al. 1996 inshore-offshore CJFAS" | Site != "SW Arm, Trinity Bay" | He != 0.875)
msat <- subset(msat, Source != "Ruzzante et al. 1998 Mol Ecol" | Site != "Placentia Bay")
msat <- subset(msat, Source != "Saillant & Gold 2006 Fish Bull" | Site != "Dauphin Island, AL 1997" & Site != "Port Fourchon, LA 1997") #both are temporal replicates of same sites reported in msat304
msat <- subset(msat, Source != "Saillant et al. 2006 ICES J Mar Sci" | Site != "Brownsville reference" & Site != "Dauphin Island reference")
msat <- subset(msat, Source != "Sala-Bozano et al. 2009 Mol Ecol" | Site != "Cadiz 07" & Site != "Foce Verde 06" & Site != "L'Estartit 07")
msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MASEIc01 (SE-Iceland)" & Site != "MASWIc01 (SW-Iceland)")
msat <- subset(msat, Source != "Stefansson et al. 2009 Heredity" | Site != "Irminger Sea 2")
msat <- subset(msat, Source != "Stefansson et al. 2009 ICES J Mar Sci" | Site != "Faroe Islands east (9)" & Site != "Norwegian international waters (13)")
msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "3RB(03)" & Site != "9PB(03)")
msat <- subset(msat, Source != "Watts et al. 2009 ICES J Mar Sci" | Site != "KEN" & Site != "TRA")
msat <- subset(msat, Source != "Yamanaka et al. 2000 Fisheries and Oceans Canada Research Report" | Site != "Barber Point (B) Sept 1999" & Site != "Bowie Seamount (B) Aug 1999")
dim(msat) #23078 rows

#remove sites reported twice
msat <- subset(msat, Source != "Larsson et al. 2010 Heredity 104:40-51" | Site != "Flatbrotten") #reported in Bekkevold et al. 2005 Evolution
  msat <- subset(msat, Source != "Larsson et al. 2010 Heredity 104:40-51" | Site != "Maseskar") #reported in Bekkevold et al. 2005 Evolution
msat <- subset(msat, Source != "Chevolot et al. 2006 Mol Ecol 15:3693-3705" | Site != "Outer Wash") #reported in Chevolot et al. 2006 J Sea Research 56:305-316
  msat <- subset(msat, Source != "Chevolot et al. 2006 Mol Ecol 15:3693-3705" | Site != "North Thames Estuary") #reported in Chevolot et al. 2006 J Sea Research 56:305-316
  msat <- subset(msat, Source != "Chevolot et al. 2006 Mol Ecol 15:3693-3705" | Site != "English Channel") #reported in Chevolot et al. 2006 J Sea Research 56:305-316
  msat <- subset(msat, Source != "Chevolot et al. 2006 Mol Ecol 15:3693-3705" | Site != "East English Channel") #reported in Chevolot et al. 2006 J Sea Research 56:305-316
  msat <- subset(msat, Source != "Chevolot et al. 2006 Mol Ecol 15:3693-3705" | Site != "Tremadog Bay") #reported in Chevolot et al. 2006 J Sea Research 56:305-316
  msat <- subset(msat, Source != "Chevolot et al. 2006 Mol Ecol 15:3693-3705" | Site != "Liverpool Pay") #reported in Chevolot et al. 2006 J Sea Research 56:305-316
msat <- subset(msat, Source != "Maggio et al. 2009 ICES J Mar Sci 66:1883-1891" | Site != "Rimini, Adriatic Sea") #reported in Garoia et al. 2004. Marine Biotechnology. 6:446-452.
  msat <- subset(msat, Source != "Maggio et al. 2009 ICES J Mar Sci 66:1883-1891" | Site != "Fano, Adriatic Sea") #reported in Garoia et al. 2004. Marine Biotechnology. 6:446-452.
  msat <- subset(msat, Source != "Maggio et al. 2009 ICES J Mar Sci 66:1883-1891" | Site != "Bari, Adriatic Sea") #reported in Garoia et al. 2004. Marine Biotechnology. 6:446-452.
  msat <- subset(msat, Source != "Maggio et al. 2009 ICES J Mar Sci 66:1883-1891" | Site != "Albania, Adriatic Sea") #reported in Garoia et al. 2004. Marine Biotechnology. 6:446-452.
msat <- subset(msat, Source != "Withler et al. 2004. Env Bio of Fishes 69:345-357.") #entire study reported twice in same dataset
msat <- subset(msat, Source != "Reilly&Ward 1999 Molecular Ecology 8:1000-1015" | Site != "Macquarie Island (Site 2)") #Macquarie Island sampled 2X in paper and don't have specific lat/lon coordinates for either site (keeping one with larger N)
msat <- subset(msat, Source != "Buonaccorsi et al. 2005 Cons Gen" | Site != "Puget Sound-Bainbridge Is." | n != 19) #reported twice in same dataset
msat <- subset(msat, Source != "Burford & Larson 2007" | Site != "Fort Bragg, CA") #previously reported in Burford & Bernardi 2008 Mar Biol 154:701-717
msat <- subset(msat, Source != "Burford 2009 J Evol Biol" | Site != "Big Creek, CA" & Site != "Monterey, CA") #previously reported in Burford & Larson 2007
  msat <- subset(msat, Source != "Burford 2009 J Evol Biol" | Site != "Gaviota, CA" & Site != "Neah Bay, WA") #previously reported in Burford & Bernardi 2008 Mar Biol 154:701-717
msat <- subset(msat, Source != "Burridge & Smolenski 2000 Mol Ecol") #previously reported in Burridge & Smolenski 2003 NZ J Mar Fres Res
msat <- subset(msat, Source != "Chapman et al. 2002 Mar Biotech" | Site != "Atlantic") #reported twice in same study (Atlantic is sum of all other sampling sites)
msat <- subset(msat, Source != "Jorgensen et al. 2005 Mol Ecol" | Site != "Gdansk Bay (GD03-3)" & Site != "Rugen (RU0503)") #previously reported in Jorgensen et al. 2005 ICES J Mar Sci
msat <- subset(msat, Source != "Lynch 2008 Thesis") #reported in msat304 (Lynch 2008 Thesis (W&M))
msat <-subset(msat, Source != "Puebla et al. 2009 Ecol" | Site != "Barb." & Site != "Bocas 2004" & Site != "Carrie B.") #previously reported in Puebla et al. 2008 Mol Ecol
msat <- subset(msat, Source != "Rhodes et al. 2003 Mar Bio" | Site != "Pohnpei 1998") #reported in msat304 (Pohnpei combined)
msat <- subset(msat, Source != "Ruzzante et al. 1998 Mol Ecol" | Site != "North Cape") #previously reported in Ruzzante et al. 1996 inshore-offshore CJFAS 
msat <- subset(msat, Source != "Ruzzante et al. 2000 J Fish Bio" | Site != "Newfoundland and Labrador") #previously reported in earlier Ruzzante studies (1996 & 1998)
msat <- subset(msat, Source != "Saillant et al. 2006 ICES J Mar Sci" | Site != "Port Aransas bycatch" & Site != "Port Aransas reference") #previously reported in Saillant & Gold 2006 Fish Bull
msat <- subset(msat, Source != "Shubina et al. 2004 Env Biol Fish") #reported in Shubina et al. 2009 Mol Biol (keeping 2009 study bc reports more markers)
msat <- subset(msat, Source != "Stefansson et al. 2009 ICES J Mar Sci" | Site != "Irminger Sea archive (1)" & Site != "Irminger Sea deep northeast (3)" & Site != "Irminger Sea deep northeast (4)" & Site != "Irminger Sea shallow southwest (5)" & Site != "Irminge rSea shallow southwest(6)") #previously reported in Stefansson et al. 2009 Heredity
msat <- subset(msat, Source != "Watts et al. 2009 ICES J Mar Sci" | Site != "ZO1" & Site != "ZO2" & Site != "ZO3" & Site != "ZO4" & Site != "Z05" & Site != "Z06") #reported in Watts et al. 2004 J Sea Res
dim(msat) #22778 rows

#remove EST-based markers
msat <- subset(msat, Source != "Kim et al. 2009 Cons Gen") #entire study uses EST markers
dim(msat) #22697 rows

######## Process allele frequencies into He where needed ########

#check that allele freqs sum to 1
x <- apply(msat[, grep('p[[:digit:]]+', names(msat))], MARGIN = 1, FUN = sumna)
msat[which(x < 0.99), ] #6, all Em-07 in Antoro et al. 2006 Marine Biotechnology 8:17-26 --> checked paper, and allele freqs don't sum to 1 there either
dim(msat) #22697 rows

msat <- subset(msat, Source != "Antoro et al. 2006 Marine Biotechnology 8:17-26" | MarkerName != "Em-07") #remove Em-07 from dataset
dim(msat) #22691 rows
msat[which(x < 0.99 & is.na(msat$He)), ] #0

#calculate He
inds <- is.na(msat$He) & !is.na(msat$p1)
sum(inds) #54
msat[inds, ] #all LeClair et al. 2006 TAFS
msat$He[inds] <- apply(msat[inds, grep('p[[:digit:]]+', names(msat))], MARGIN = 1, FUN = calcHe)

######## Calc lat and lon in decimal degrees ########

#remove lines without lat/lon
inds <- is.na(msat$lat_deg) | is.na(msat$lon_deg)
sum(inds) #442
msat[inds, c('spp', 'Source', 'Site', 'lat_deg', 'lon_deg', 'file')] #site description is too vague or no location described in the paper
msat <- msat[!inds, ]
nrow(msat) #22249 rows

#remove studies where lat/lon provided but averaged across too broad a range
dim(msat) #22249 rows
msat <- subset(msat, Source != "Zhang et al. 2001 Mol Ecol Notes") #remove bc sites span too broad a range
msat <- subset(msat, Source != "Nielsen et al. 2010 Cons Gen 11:999-1012" | Site != "Aleutian Islands, AK") #remove bc site span too broad a range
msat <- subset(msat, Source != "Chevolot et al. 2006 Mol Ecol 15:3693-3705" | Site != "Adriatic Sea" & Site != "Flores") #remove bc sites span too broad a range
msat <- subset(msat, Source != "Verissimo et al. 2010 Mol Ecol 8:1651-1662") #removing entire study because can't independently verify data
msat <- subset(msat, Source != "Vagelli, Burford & Bernardi 2009 Marine Genomics 1:129-134") #removing entire study because can't independently verify data
msat <- subset(msat, Source != "O'Leary et al. 2007 Journal of Fish Biology 70:310-335" | Site != "Scotian Shelf" & Site != "Faeroes Ridge" & Site != "Barents Sea" & Site != "West Greenland" & Site != "East Greenland")
msat <- subset(msat, Source != "Nielsen et al. 2003 Mol Ecol 12:1497-1508" | Site != "North Sea") #remove bc sites span too broad a range
msat <- subset(msat, Source != "Arigoni & Largiader 2000 Molecular Ecology 9:2155-2234") #remove entire study bc He is reported as average across populations
msat <- subset(msat, Source != "Feldheim et al. 2001 Molecular Ecology 10:295-303") #remove entire study bc all sites pooled to calculate global He
msat <- subset(msat, Source != "Nugruho et al. 2001 Fisheries Science 67:843-850" | Site != "NSW coast") #remove bc gives latitude but not longitude in paper
msat <- subset(msat, Source != "Anderson 2007 Fish Bull") #remove entire study bc all sites span too broad a range
msat <- subset(msat, Source != "Carlsson et al. 2004 Mol Ecol") #remove entire study bc all sites span too broad a range
msat <- subset(msat, Source != "Catanese et al. 2007 Mol Ecol Notes") #remove entire study bc no coords provided (just ocean basins)
msat <- subset(msat, Source != "Clark et al. 2004 Mol Ecol Notes") #remove entire study bc no coords provided
msat <- subset(msat, Source != "Durand et al. 2005 Mar Bio" | Site != "Peru" & Site != "Sri Lanka") #remove bc sites span too broad a range
msat <- subset(msat, Source != "Garoia et al. 2007 Mol Ecol") #remove entire study bc all sites span too broad a range
msat <- subset(msat, Source != "Gonzalez & Zardoya 2007 Mol Ecol Notes") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Hutchinson et al. 2001 MEPS") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Kim et al. 2010 Fish Sci" | Site != "Russia" & Site != "USA") #remove sites bc no coords provided (just vague locations)
msat <- subset(msat, Source != "Lage et al. 2001 CJFAS") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Landi et al. 2005 Mol Ecol Notes") #remove entire study bc no coords provided
msat <- subset(msat, Source != "Li et al. 2006 Mol Ecol Notes") #remove entire study bc no coords provided
msat <- subset(msat, Source != "Machado-Schiaffino & Garcia-Vazquez 2009") #remove entire study bc no coords provided
msat <- subset(msat, Source != "Miller et al. 2000 Mol Ecol") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Miller et al. 2001 Mol Ecol Notes") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Moran et al. 1999 Mol Ecol") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Narum et al. 2004 Copeia") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Nugroho et al. 1998 Fish Sci") #remove entire study bc no coords provided
msat <- subset(msat, Source != "Pereyra et al. 2004 Mol Ecol Notes") #remove entire study bc no coords provided
msat <- subset(msat, Source != "Perez-Enriquez & Taniguchi 1999 Fish Sci" | Site != "Australia" & Site != "New Zealand") #remove sites bc span too broad a range
msat <- subset(msat, Source != "Ramsak et al. 2003  Mol Ecol Notes") #remove entire study bc no coords provided
msat <- subset(msat, Source != "Riccioni et al. 2010 PNAS" | Site != "Algerian coast") #remove site bc spans too broad a range
msat <- subset(msat, Source != "Ruzzante et al. 2000 Ecol Appl") #remove entire study bc spans too broad a range (also some inds previously reported in earlier studies)
msat <- subset(msat, Source != "Saxton 2009 Thesis") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "Norway and Iceland") #remove site bc spans too broad a range
msat <- subset(msat, Source != "Seibert & Ruzzante et al. 2006 Mol Ecol Notes") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Sekino & Hara 2000 Mol Ecol") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Sekino et al. 2000 Mol Ecol") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Takagi et al. 2001 Fish Bull" | Site != "NW Pacific" & Site != "SE Pacific" & Site != "SW Atlantic" & Site != "SW Pacific") #remove sites bc no coords provided
msat <- subset(msat, Source != "Tzeng et al. 2009 ICES J Mar Sci") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "White et al. 2009 Cons Gen") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Zhang et al. 2001 Mol Ecol Notes") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Zhang et al. 2006 Mol Ecol Notes") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Zhang et al. 2010 Cons Gen Res") #remove entire study bc spans too broad a range
msat <- subset(msat, Source != "Zhu et al. 2005 Mol Ecol Notes") #remove entire study bc spans too broad a range
dim(msat) #21757 rows

#remove sites in black sea, Caspian sea, lakes/rivers or of potentially farmed sample
dim(msat) #21757 rows
msat <- subset(msat, Source != "Cemal Turan (2015) Biochemical Systematics and Ecology 63: 174e182" | Site != "Bulgarian Varna Coast" 
               & Site != "Black Sea Igneada" & Site != "the Black Sea Samsun" & Site != "Black Sea Duzce" & Site != "Istanbul Bosporus " 
               & Site != "the Black Sea Trabzon " & Site != "Marmara Sea BandÄ±rma") #excluded bc in Black Sea
msat <- subset(msat, Source != "Boissin et al. (2016) Molecular Ecology 25:2195-2209" | Country != "Ukraine" & Country != "Georgia" 
               & Country != "Bulgaria" & Country != "Romania") #excluded bc in Black Sea
msat <- subset(msat, Source != "Boissin et al. (2016) Molecular Ecology 25:2195-2209" | Site != "Sinop" & Site != "Sile") #excluded bc in Black Sea
msat <- subset(msat, Source != "L. C. Woodall et al. (2015) Conserv Genet 16:1139â€“1153" | Site != "VBU: Varna") #excluded bc in Black Sea
msat <- subset(msat, Source != "Miralles, Juanes & Garcia-Vazquez (2014) Transactions of the American Fisheries Society 143:1308-1315" | Site != "Istanbul") #excluded bc in Black Sea
msat <- subset(msat, Source != "Sala-Bozano et al. 2009 Mol Ecol" | Site != "Foce Verde 07" & Site != "Foce Verde 06") #excluded bc in Black Sea
msat <- subset(msat, Source != "Feldheim et al. 2009 Mol. Ecol. Res. 9:639-644" | Site != "Caspian Sea, Niyazabad") #excluded bc in Caspian Sea
msat <- subset(msat, Source != "Ma et al. 2011 Env Fish Bio" | Site != "Wuhan region, Yangtze River") #excluded bc in freshwater river
msat <- subset(msat, Source != "Ding et al. 2009 Cons Gen") #exclude whole study bc is of Atlantic halibut but sampled off China (farmed?)
dim(msat) #21701 rows	

#change minutes & seconds to negative if lat/long is negative
msat$lat_min <- ifelse(msat$lat_deg < 0, -msat$lat_min, msat$lat_min) #imp so things sum properly in next step
msat$lat_sec <- ifelse(msat$lat_deg < 0, -msat$lat_sec, msat$lat_sec)
msat$lon_min <- ifelse(msat$lon_deg < 0, -msat$lon_min, msat$lon_min)
msat$lon_sec <- ifelse(msat$lon_deg < 0, -msat$lon_sec, msat$lon_sec)

#turn minutes & seconds into decimals
msat$lat <- rowSums(cbind(msat$lat_deg, msat$lat_min/60, msat$lat_sec/3600), na.rm = TRUE)
msat$lon <- rowSums(cbind(msat$lon_deg, msat$lon_min/60, msat$lon_sec/3600), na.rm = TRUE)

######## Fix species name mistakes ########

#get list of species names
msat$spp <- gsub(' $', '', msat$spp) #remove trailing spaces
spps <- sort(unique(msat$spp)) #check list against FB

#fix species name mistakes
msat$spp[msat$spp == "Apogon doederleini"] <- "Ostorhinchus doederleini" #FB calls this Ostorhinchus
msat$spp[msat$spp == "Aphanius fascitus"] <- "Aphanius fasciatus" #Add missing a
msat$spp[msat$spp == "Centroselachus crepidater"] <- "Centroscymnus crepidater" #FB calls this Centroscymnus
msat$spp[msat$spp == "Chrysophrys auratus"] <- "Pagrus auratus" #FB calls this Pagrus
msat$spp[msat$spp == "Clupea pallasii"] <- "Clupea pallasii pallasii" #Add extra pallasii
msat$spp[msat$spp == "Diplodus sargus sargus"] <- "Diplodus sargus" #Remove extra sargus
msat$spp[msat$spp == "Epinephelus acanthistius"] <- "Hyporthodus acanthistius" #FB calls this Hyporthodus
msat$spp[msat$spp == "Epinephelus septemfasciatus"] <- "Hyporthodus septemfasciatus" #FB calls this Hypothordus
msat$spp[msat$spp == "Gnatholepsis cauerensis"] <- "Gnatholepis cauerensis" #fix extra s
msat$spp[msat$spp == "Hapalogenys nitens"] <- "Hapalogenys nigripinnis" #FB calls this nigripinnis
msat$spp[msat$spp == "Lepidonotothen larseni"] <- "Nototheniops larseni" #FB calls this Norotheniops
msat$spp[msat$spp == "Liza affinis"] <- "Planiliza affinis" #FB calls this Planiliza
msat$spp[msat$spp == "Lophius piscatorious"] <- "Lophius piscatorius" #Remove extra o
msat$spp[msat$spp == "Macrourus breglax"] <- "Macrourus berglax" #Fixed misspelling
msat$spp[msat$spp == "Mullus barbatus"] <- "Mullus barbatus barbatus" #Add extra barbatus
msat$spp[msat$spp == "Navodon septentrionalis"] <- "Thamnaconus septentrionalis" #FB calls this Navodon
msat$spp[msat$spp == "Pleuragramma antarcticum"] <- "Pleuragramma antarctica" #FB calls this antarctica
msat$spp[msat$spp == "Pleuronectes herzensteini"] <- "Pseudopleuronectes herzensteini" #FB calls this Pseudopleuronectes
msat$spp[msat$spp == "Pleuronectes yokohamae"] <- "Pseudopleuronectes yokohamae" #FB calls this Pseudopleuronectes
msat$spp[msat$spp == "Psetta maxima"] <- "Scophthalmus maximus" #FB calls this Scophthalmus
msat$spp[msat$spp == "Pseudosciaena crocea"] <- "Larimichthys crocea" #FB calls this Larimichthys
msat$spp[msat$spp == "Pseudosciaena polyactis"] <- "Larimichthys polyactis" #FB calls this Larimichthys
msat$spp[msat$spp == "Sarda Sarda"] <- "Sarda sarda" #Change to lowercase species name
msat$spp[msat$spp == "Sebastes marinus"] <- "Sebastes norvegicus" #FB calls this norvegicus
msat$spp[msat$spp == "Sebastes schlegeli"] <- "Sebastes schlegelii" #Fixed misspelling
msat$spp[msat$spp == "Tetrapturus albidus"] <- "Kajikia albida" #FB calls this Kajikia albida
msat$spp[msat$spp == "Theragra chalcogramma"] <- "Gadus chalcogrammus" #FB calls this Gadus chalcogrammus

#fix species that were mis-named in dataset
msat$spp[msat$Source == "Takagi et al. 1999 Thunnus orientalis Fish Sci" & msat$CommonName == "Yellowfin tuna"] <- "Thunnus albacares" #originally labeled Thunnus alalunga but paper says yellowfin tuna

#remove species that aren't truly marine
dim(msat) #21695 rows
msat <- subset(msat, spp != "Alosa alosa" & spp != "Alosa fallax" & spp != "Alosa sapidissima" & 
                 spp != "Aphanius fasciatus" & spp!= "Coregonus lavaretus" & spp != "Cynoscion acoupa" & 
                 spp!= "Esox lucius" & spp != "Fundulus heteroclitus" & spp!= "Gasterosteus aculeatus" & 
                 spp != "Hypomesus transpacificus" & spp != "Mallotus villosus" & spp != "Osmerus eperlanus" & 
                 spp != "Osmerus mordax" & spp != "Odontesthes argentinensis" & spp != "Platichthys flesus" & 
                 spp != "Platichthys stellatus" & spp != "Pomatoschistus microps" & spp != "Pomatoschistus minutus" &
                 spp != "Pungitius pungitius")
dim(msat) #20861 rows

######## Fix common name mistakes ########

#get list of common names
msat$CommonName <- gsub(' $', '', msat$CommonName) #remove trailing spaces
com_spps <- msat[!duplicated(msat[, c('spp', 'CommonName')]), c('spp', 'CommonName')] #grab unique spp/CommonName combos
com_spps <- com_spps[order(com_spps$spp, com_spps$CommonName), ] #order by spp then by CommonName, check against FB

#fix common name mistakes based on FB suggestion
msat$CommonName[msat$spp == "Acanthopagrus australis"] <- "Yellowfin bream"
msat$CommonName[msat$spp == "Acanthopagrus taiwanensis"] <- "Taiwan picnic seabream"
msat$CommonName[msat$spp == "Acanthurus nigricans"] <- "Whitecheek surgeonfish"
msat$CommonName[msat$spp == "Aetobatus narinari"] <- "Whitespotted eagle ray"
msat$CommonName[msat$spp == "Amphiprion melanopus"] <- "Fire clownfish"
msat$CommonName[msat$spp == "Brevoortia gunteri"] <- "Finescale menhaden"
msat$CommonName[msat$spp == "Brosme brosme"] <- "Cusk"
msat$CommonName[msat$spp == "Chaetodon tricinctus"] <- "Three-striped butterflyfish"
msat$CommonName[msat$spp == "Chaetodon trifasciatus"] <- "Melon butterflyfish"
msat$CommonName[msat$spp == "Coris julis"] <- "Mediterranean rainbow wrasse"
msat$CommonName[msat$spp == "Coryphaenoides brevibarbis"] <- "Shortbeard grenadier"
msat$CommonName[msat$spp == "Ctenolabrus rupestris"] <- "Goldsinny-wrasse"
msat$CommonName[msat$spp == "Cynoscion arenarius"] <- "Sand weakfish"
msat$CommonName[msat$spp == "Cynoscion nothus"] <- "Silver seatrout"
msat$CommonName[msat$spp == "Diplodus hottentotus"] <- "Zebra"
msat$CommonName[msat$spp == "Electrona antarctica"] <- "Antarctic lanternfish" 
msat$CommonName[msat$spp == "Eleginops maclovinus"] <- "Patagonian blennie"
msat$CommonName[msat$spp == "Epinephelus bruneus"] <- "Kelp grouper"
msat$CommonName[msat$spp == "Epinephelus caninus"] <- "Dogtooth grouper"
msat$CommonName[msat$spp == "Epinephelus costae"] <- "Goldblotch grouper"
msat$CommonName[msat$spp == "Epinephelus marginatus"] <- "Dusky grouper"
msat$CommonName[msat$spp == "Epinephelus morio"] <- "Red grouper"
msat$CommonName[msat$spp == "Etmopterus spinax"] <- "Velvet belly"
msat$CommonName[msat$spp == "Girella laevifrons"] <- "" #no common name
msat$CommonName[msat$spp == "Glaucosoma hebraicum"] <- "West Australian dhufish"
msat$CommonName[msat$spp == "Gobionotothens gibberifrons"] <- "Humped rockcod"
msat$CommonName[msat$spp == "Gymnothorax chilospilus"] <- "Lipspot moray"
msat$CommonName[msat$spp == "Hapalogenys nigripinnis"] <- "Short barbeled velvetchin"
msat$CommonName[msat$spp == "Hippocampus angustus"] <- "Western spiny seahorse"
msat$CommonName[msat$spp == "Hippocampus hippocampus"] <- "Short-snouted seahorse"
msat$CommonName[msat$spp == "Hyporthodus acanthistius"] <- "Rooster hind"
msat$CommonName[msat$spp == "Kajikia albida"] <- "Atlantic white marlin"
msat$CommonName[msat$spp == "Larimichthys polyactis"] <- "Yellow croaker"
msat$CommonName[msat$spp == "Lepidorhombus boscii"] <- "Four-spot megrim"
msat$CommonName[msat$spp == "Lutjanus campechanus"] <- "Northern red snapper"
msat$CommonName[msat$spp == "Macrourus berglax"] <- "Roughhead grenadier"
msat$CommonName[msat$spp == "Menidia menidia"] <- "Atlantic silverside"
msat$CommonName[msat$spp == "Microchirus azevia"] <- "Bastard sole"
msat$CommonName[msat$spp == "Molva molva"] <- "Ling"
msat$CommonName[msat$spp == "Mullus barbatus barbatus"] <- "Red mullet"
msat$CommonName[msat$spp == "Mullus surmuletus"] <- "Surmullet"
msat$CommonName[msat$spp == "Mustelus mustelus"] <- "Smooth-hound"
msat$CommonName[msat$spp == "Mycteroperca phenax"] <- "Scamp"
msat$CommonName[msat$spp == "Mycteroperca rubra"] <- "Mottled grouper"
msat$CommonName[msat$spp == "Oblada melanura"] <- "Saddled seabream"
msat$CommonName[msat$spp == "Pagrus auratus"] <- "Silver seabream"
msat$CommonName[msat$spp == "Palabrax albomaculatus"] <- "Camotillo"
msat$CommonName[msat$spp == "Paralabrax nebulifer"] <- "Barred sand bass"
msat$CommonName[msat$spp == "Polyprion americanus"] <- "Wreckfish"
msat$CommonName[msat$spp == "Pseudopleuronectes herzensteini"] <- "Yellow striped flounder"
msat$CommonName[msat$spp == "Pseudopleuronectes yokohamae"] <- "Marbled flounder"
msat$CommonName[msat$spp == "Pterapogon kauderni"] <- "Banggai cardinalfish"
msat$CommonName[msat$spp == "Rhomboplites aurorubens"] <- "Vermilion snapper"
msat$CommonName[msat$spp == "Sardina pilchardus"] <- "European pilchard"
msat$CommonName[msat$spp == "Sciaenops ocellatus"] <- "Red drum"
msat$CommonName[msat$spp == "Scomber australasicus"] <- "Blue mackerel"
msat$CommonName[msat$spp == "Scomber japonicus"] <- "Chub mackerel"
msat$CommonName[msat$spp == "Scomberomorus commerson"] <- "Narrow-barred Spanish mackerel"
msat$CommonName[msat$spp == "Scyliorhinus canicula"]<- "Lesser spotted dogfish"
msat$CommonName[msat$spp == "Sebastes auriculatus"] <- "Brown rockfish"
msat$CommonName[msat$spp == "Sebastes borealis"] <- "Shortraker rockfish"
msat$CommonName[msat$spp == "Sebastes caurinus"] <- "Copper rockfish"
msat$CommonName[msat$spp == "Sebastes fasciatus"] <- "Acadian redfish"
msat$CommonName[msat$spp == "Sebastes flavidus"] <- "Yellowtail rockfish"
msat$CommonName[msat$spp == "Sebastes inermis"] <- "Dark-banded rockfish"
msat$CommonName[msat$spp == "Sebastes melanops"] <- "Black rockfish"
msat$CommonName[msat$spp == "Sebastes mentella"] <- "Beaked redfish"
msat$CommonName[msat$spp == "Sebastes norvegicus"] <- "Golden redfish"
msat$CommonName[msat$spp == "Sebastes paucispinis"] <- "Bocaccio rockfish"
msat$CommonName[msat$spp == "Sebastes pinniger"] <- "Canary rockfish"
msat$CommonName[msat$spp == "Sebastes ruberrimus"] <- "Yelloweye rockfish"
msat$CommonName[msat$spp == "Sebastes schlegelii"] <- "Korean rockfish"
msat$CommonName[msat$spp == "Sebastes thompsoni"] <- "Goldeye rockfish"
msat$CommonName[msat$spp == "Sebastes viviparus"] <- "Norway redfish"
msat$CommonName[msat$spp == "Sebastiscus marmoratus"] <- "False kelpfish"
msat$CommonName[msat$spp == "Seriola lalandi"] <- "Yellowtail amberjack"
msat$CommonName[msat$spp == "Serranus atricauda"] <- "Blacktail comber"
msat$CommonName[msat$spp == "Serranus cabrilla"] <- "Comber"
msat$CommonName[msat$spp == "Serranus hepatus"] <- "Brown comber"
msat$CommonName[msat$spp == "Serranus scriba"] <- "Painted comber"
msat$CommonName[msat$spp == "Siganus fuscescens"] <- "Mottled spinefoot"
msat$CommonName[msat$spp == "Siganus spinus"] <- "Little spinefoot"
msat$CommonName[msat$spp == "Solea senegalensis"] <- "Senegalese sole"
msat$CommonName[msat$spp == "Solea solea"] <- "Common sole"
msat$CommonName[msat$spp == "Sparus aurata"] <- "Gilthead seabream"
msat$CommonName[msat$spp == "Sphyrna lewini"] <- "Scalloped hammerhead"
msat$CommonName[msat$spp == "Spinachia spinachia"] <- "Sea stickleback"
msat$CommonName[msat$spp == "Sprattus fuegensis"] <- "Falkland sprat"
msat$CommonName[msat$spp == "Thalassoma bifasciatum"] <- "Bluehead"
msat$CommonName[msat$spp == "Thunnus alalunga"] <- "Albacore"
msat$CommonName[msat$spp == "Thunnus thynnus"] <- "Atlantic bluefin tuna"
msat$CommonName[msat$spp == "Trachinotus goodei"] <- "Great pompano"
msat$CommonName[msat$spp == "Trachurus japonicus"] <- "Japanese jack mackerel"
msat$CommonName[msat$spp == "Trachurus murphyi"] <- "Chilean jack mackerel"

##########################################################################################################################################

######## QA/QC ########

######## Character check ########

#were numeric fields read properly?
summary(msat) #repeats are NOT all numeric

#double-check species names visually
t(t(sort(unique(msat$spp)))) #print in one column (351 species)

#fix study names where mis-reported
msat$Source[msat$Source == "Galarza et al. 2007 Conserv Genet 8:1251-1253"] <- "Galarza et al. 2007 Molecular Ecology Notes 7:230-232"

#fix CollectionYear for studies where mis-reported
msat$CollectionYear[msat$Source == "Bisol et al. 2007 Hydrobiologia 577:151-159" & msat$Site == "Chioggia"] <- 2001
  msat$CollectionYear[msat$Source == "Bisol et al. 2007 Hydrobiologia 577:151-159" & msat$Site == "Lido"] <- 2001
  msat$CollectionYear[msat$Source == "Bisol et al. 2007 Hydrobiologia 577:151-159" & msat$Site == "Lago dei Teneri"] <- 2001
msat$CollectionYear[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Atlantic Iberia"] <- 2000
  msat$CollectionYear[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Bay of Biscay"] <- 2001

#fix N for studies where mis-reported
msat$n[msat$Source == "Bisol et al. 2007 Hydrobiologia 577:151-159" & msat$Site == "Chioggia"] <- 26
  msat$n[msat$Source == "Bisol et al. 2007 Hydrobiologia 577:151-159" & msat$Site == "Lido"] <- 25
  msat$n[msat$Source == "Bisol et al. 2007 Hydrobiologia 577:151-159" & msat$Site == "Lago dei Teneri"] <- 36
msat$n[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Atlantic Iberia"] <- 204
  msat$n[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Bay of Biscay"] <- 58
msat$n[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Turneffe Atoll, Belize"] <- 44
  msat$n[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Turneffe Atoll, Belize" & msat$MarkerName == "SpAAT40-1"] <- 37
  msat$n[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Roatan Island"] <- 35
  msat$n[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Roatan Island" & msat$MarkerName == "SpAAC42-1"] <- 34
  msat$n[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Roatan Island" & msat$MarkerName == "SpAAT40-1"] <- 32
  msat$n[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Roatan Island" & msat$MarkerName == "SpAAC41-1"] <- 33

#check repeat definitions
msat$Repeat[msat$Repeat %in% c('?', '') & !is.na(msat$Repeat)] <- NA #fix ? in repeat definitions
x <- as.numeric(msat$Repeat) #if Repeat not numeric then coerced to NA in X
msat$Repeat[is.na(x) & !is.na(msat$Repeat)] #pulls instances where Repeat col is not numeric
#only one case where repeat is "9 or 2." Leave for now.
  
#add Repeat info based on papers
msat$Repeat[msat$Source == "Appleyard, Williams & Ward. 2004. CCAMLR Science 11:21-32." & msat$MarkerName == "To2"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Appleyard, Williams & Ward. 2004. CCAMLR Science 11:21-32." & msat$MarkerName == "To5"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$MarkerName == "Gmo1"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$MarkerName == "Gmo2"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$MarkerName == "Gmo120"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$MarkerName == "Gmo132"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$MarkerName == "Gmo141"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$MarkerName == "Gmo8"] <- 4 #based on ref paper
  msat$Repeat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$MarkerName == "Gmo19"] <- 4 #based on ref paper
  msat$Repeat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$MarkerName == "Gmo34"] <- 4 #based on ref paper
  msat$Repeat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$MarkerName == "Gmo37"] <- 4 #based on ref paper
msat$Repeat[msat$Source == "Smith&McVeagh 2000 Journal of Fish Biology 57:72-83" & msat$MarkerName == "To2"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Smith&McVeagh 2000 Journal of Fish Biology 57:72-83" & msat$MarkerName == "To3"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Smith&McVeagh 2000 Journal of Fish Biology 57:72-83" & msat$MarkerName == "To5"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Shaw et al. 1999 Heredity 83:490-499" & msat$MarkerName == "Cha17"] <- 2 #based on genbank accession #
  msat$Repeat[msat$Source == "Shaw et al. 1999 Heredity 83:490-499" & msat$MarkerName == "Cha20"] <- 2 #based on genbank accession #
  msat$Repeat[msat$Source == "Shaw et al. 1999 Heredity 83:490-499" & msat$MarkerName == "Cha63"] <- 2 #based on genbank accession #
  msat$Repeat[msat$Source == "Shaw et al. 1999 Heredity 83:490-499" & msat$MarkerName == "Cha113"] <- 2 #based on genbank accession #
msat$Repeat[msat$Source == "McGowan&Reith 1999 Molecular Ecology 8:1761-1763" & msat$MarkerName == "HhiA44"] <- 2 #based on genbank accession #
  msat$Repeat[msat$Source == "McGowan&Reith 1999 Molecular Ecology 8:1761-1763" & msat$MarkerName == "HhiC17"] <- 2 #based on genbank accession #
  msat$Repeat[msat$Source == "McGowan&Reith 1999 Molecular Ecology 8:1761-1763" & msat$MarkerName == "HhiD34"] <- 2 #based on genbank accession #
  msat$Repeat[msat$Source == "McGowan&Reith 1999 Molecular Ecology 8:1761-1763" & msat$MarkerName == "HhiI29"] <- 2 #based on genbank accession #
  msat$Repeat[msat$Source == "McGowan&Reith 1999 Molecular Ecology 8:1761-1763" & msat$MarkerName == "HhiJ42"] <- 2 #based on genbank accession #
msat$Repeat[msat$Source == "McPherson et al. 2001 Journal of Fish Biology 59:356-370"] <- 4 #based on ref papers
msat$Repeat[msat$Source == "Nugruho et al. 2001 Fisheries Science 67:843-850"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "McPherson et al. 2001 Molecular Ecology Notes 1: 31-32"] <- 4 #provided in paper (primer note)
msat$Repeat[msat$Source == "Anderson & Karel 2007 J Fish Bio"] <- 3.4 #based on ref paper
msat$Repeat[msat$Source == "Anderson & Karel 2010 Fish Aq J" & msat$MarkerName == "SOC12"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Anderson & Karel 2010 Fish Aq J" & msat$MarkerName == "SOC243"] <- 3 #based on ref papers
  msat$Repeat[msat$Source == "Anderson & Karel 2010 Fish Aq J" & msat$MarkerName == "SOC412"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Anderson & Karel 2010 Fish Aq J" & msat$MarkerName == "SOC415"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Anderson & Karel 2010 Fish Aq J" & msat$MarkerName == "SOC432"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Anderson & Karel 2010 Fish Aq J" & msat$MarkerName == "SOC50"] <- 2 #based on ref papers
msat$Repeat[msat$Source == "Anderson et al. 2009 Fish Bull"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Anderson et al. 2009 Fish Bull" & msat$MarkerName == "SOC243"] <- 3 #based on ref papers
msat$Repeat[msat$Source == "Bentzen et al. 1996 CJFAS"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Bernal-Ramirez et al. 2003 Mar Bio"] <- 2 #based on ref papers
msat$Repeat[msat$Source == "Bradbury et al. 2009 J Fish Bio"] <- 4 #based on ref papers
  msat$Repeat[msat$Source == "Bradbury et al. 2009 J Fish Bio" & msat$MarkerName == "Gmo2"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Bradbury et al. 2009 J Fish Bio" & msat$MarkerName == "Gmo35"] <- 3 #based on ref papers
msat$Repeat[msat$Source == "Beldade et al. 2009 Mol Ecol Res/Abercrombie eta l. 2009 Mol Ecol Res"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Burford & Larson 2007"] <- 2.43 #based on paper
msat$Repeat[msat$Source == "Burford 2009 J Evol Biol"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Burford 2009 J Evol Biol" & msat$MarkerName == "Sra.15-8"] <- 4 #based on ref paper
  msat$Repeat[msat$Source == "Burford 2009 J Evol Biol" & msat$MarkerName == "Sra.16-5"] <- 4 #based on ref paper
msat$Repeat[msat$Source == "Chapman et al. 2002 Mar Biotech" & msat$MarkerName == "Soc014"] <- 2 #based on paper
  msat$Repeat[msat$Source == "Chapman et al. 2002 Mar Biotech" & msat$MarkerName == "Soc017"] <- 2 #based on paper
  msat$Repeat[msat$Source == "Chapman et al. 2002 Mar Biotech" & msat$MarkerName == "Soc029"] <- 2 #based on paper
msat$Repeat[msat$Source == "Charrier et al. 2007 MEPS" & msat$MarkerName == "Gmo02"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Charrier et al. 2007 MEPS" & msat$MarkerName == "MpourBW13"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Charrier et al. 2007 MEPS" & msat$MarkerName == "Tch10"] <- 4 #based on ref paper
  msat$Repeat[msat$Source == "Charrier et al. 2007 MEPS" & msat$MarkerName == "Tch20"] <- 4 #based on ref paper
  msat$Repeat[msat$Source == "Charrier et al. 2007 MEPS" & msat$MarkerName == "Tch8"] <- 4 #based on ref paper
msat$Repeat[msat$Source == "Danancher & Garcia-Vazquez 2009 Mar Bio"] <- 2 #based on paper
msat$Repeat[msat$Source == "Durand et al. 2005 Mar Bio" & msat$MarkerName == "cmrTA-208"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Florin & Hoglund 2007 Mol Ecol"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Florin & Hoglund 2007 Mol Ecol" & msat$MarkerName == "Sma1-125"] <- 4 #based on ref paper
msat$Repeat[msat$Source == "Galarza et al. 2009 PNAS"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Galarza et al. 2009 PNAS" & msat$MarkerName == "St245"] <- 4 #based on ref paper
msat$Repeat[msat$Source == "Gold et al. 2010 Fish Res"] <- 4 #based on ref paper
  msat$Repeat[msat$Source == "Gold et al. 2010 Fish Res" & msat$MarkerName == "Sbr19"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Gold et al. 2010 Fish Res" & msat$MarkerName == "Sbr26"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Gold et al. 2010 Fish Res" & msat$MarkerName == "Sbr35"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Gold et al. 2010 Fish Res" & msat$MarkerName == "Sca49"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Gold et al. 2010 Fish Res" & msat$MarkerName == "Sca61"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Gold et al. 2010 Fish Res" & msat$MarkerName == "Sni13"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Gold et al. 2010 Fish Res" & msat$MarkerName == "Sni26"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Gold et al. 2010 Fish Res" & msat$MarkerName == "Sni29"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Gomez-Uchida & Banks 2005 CJFAS"] <- 3.29 #based on ref papers
msat$Repeat[msat$Source == "Gomez-Uchida 2006 Dissertation"] <- 3.3 #based on ref papers
msat$Repeat[msat$Source == "Gonzalez et al. 2008 BMC Evol"] <- 2 #based on ref papers
msat$Repeat[msat$Source == "Gonzalez-Wanguemert et al. 2010 JEMBE"] <- 2 #based on paper
msat$Repeat[msat$Source == "Hauser et al. 2002 PNAS"] <- 2 #based on ref papers
msat$Repeat[msat$Source == "Hoarau et al. 2002 Mol Ecol"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Hoarau et al. 2002 Mol Ecol" & msat$MarkerName == "List1003"] <- 3 #based on ref paper
msat$Repeat[msat$Source == "Jorgensen et al. 2005 Mol Ecol"] <- 4 #based on ref papers
msat$Repeat[msat$Source == "Jue 2010 Thesis"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Jue 2010 Thesis" & msat$MarkerName == "MBO048"] <- 4 #based on ref paper
msat$Repeat[msat$Source == "Karlsson et al. 2009 Mar Bio"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Kasapidis & Magoulas 2008 Fish Res"] <- 2 #based on paper
msat$Repeat[msat$Source == "Lage & Kornfield 1999 Mol Ecol"] <- 2 #based on paper
msat$Repeat[msat$Source == "Limborg et al. 2009 MEPS"] <- 2 #based on paper
msat$Repeat[msat$Source == "McClelland et al. 2005 J Fish Bio"] <- 2 #based on paper
msat$Repeat[msat$Source == "Mitchell 2006 Thesis"] <- 4 #based on ref papers
  msat$Repeat[msat$Source == "Mitchell 2006 Thesis" & msat$MarkerName == "Cha113"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Mitchell 2006 Thesis" & msat$MarkerName == "Cha123"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Mitchell 2006 Thesis" & msat$MarkerName == "Cha17"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Mitchell 2006 Thesis" & msat$MarkerName == "Cha20"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Mitchell 2006 Thesis" & msat$MarkerName == "Cha63"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Morishima et al. 2009 Mol Ecol Res"] <- 2 #based on paper
  msat$Repeat[msat$Source == "Morishima et al. 2009 Mol Ecol Res" & msat$MarkerName == "Kto32"] <- 4 #based on paper
msat$Repeat[msat$Source == "Nielsen et al. 2004 Mol Ecol"] <- 2 #based on paper
  msat$Repeat[msat$Source == "Nielsen et al. 2004 Mol Ecol" & msat$MarkerName == "Sma-125"] <- 4 #based on paper
msat$Repeat[msat$Source == "Olsen et al. 2002 Mol Ecol Notes"] <- 4 #based on paper
msat$Repeat[msat$Source == "Pampoulie & Danielsdottir 2008 Fish Res"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Pampoulie & Danielsdottir 2008 Fish Res" & msat$MarkerName == "Sal1"] <- 4 #based on ref paper
  msat$Repeat[msat$Source == "Pampoulie & Danielsdottir 2008 Fish Res" & msat$MarkerName == "Sal3"] <- 5 #based on ref paper
msat$Repeat[msat$Source == "Perez-Enriquez & Taniguchi 1999 Fish Sci"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Puebla et al. 2008 Mol Ecol" & msat$MarkerName == "Pam013"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Puebla et al. 2009 Ecol" & msat$MarkerName == "pam013"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Roques et al. 2001 Mol Ecol"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Roques et al. 2002 Mar Bio"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Ruzzante et al. 1996 inshore-offshore CJFAS"] <- 2 #based on ref papers
msat$Repeat[msat$Source == "Saillant & Gold 2006 Fish Bull"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Saillant & Gold 2006 Fish Bull" & msat$MarkerName == "Prs275"] <- 3 #based on ref paper
msat$Repeat[msat$Source == "Saillant et al. 2006 ICES J Mar Sci"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Saillant et al. 2006 ICES J Mar Sci" & msat$MarkerName == "Prs275"] <- 3 #based on ref paper
msat$Repeat[msat$Source == "Sala-Bozano et al. 2009 Mol Ecol"] <- 2 #based on ref papers
msat$Repeat[msat$Source == "Shishidou et al. 2008 Nipp Suis Gakk"] <- 2 #based on ref paper
msat$Repeat[msat$Source == "Small et al. 2005 TAFS"] <- 3.5 #based on paper
msat$Repeat[msat$Source == "Stefansson et al. 2009 Heredity"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Stefansson et al. 2009 Heredity" & msat$MarkerName == "Sal1"] <- 4 #based on ref paper
  msat$Repeat[msat$Source == "Stefansson et al. 2009 Heredity" & msat$MarkerName == "Sal3"] <- 5 #based on ref paper
  msat$Repeat[msat$Source == "Stefansson et al. 2009 Heredity" & msat$MarkerName == "Smen10"] <- 4 #based on paper
  msat$Repeat[msat$Source == "Stefansson et al. 2009 Heredity" & msat$MarkerName == "Smen5"] <- 4 #based on paper
msat$Repeat[msat$Source == "Was et al. 2008 ICES J Mar Sci"] <- 2 #based on ref paper
  msat$Repeat[msat$Source == "Was et al. 2008 ICES J Mar Sci" & msat$MarkerName == "MmerUEAw01"] <- 4 #based on ref papers
  msat$Repeat[msat$Source == "Was et al. 2008 ICES J Mar Sci" & msat$MarkerName == "Tch10"] <- 4 #based on ref paper
msat$Repeat[msat$Source == "Was et al. 2010 Mar Bio"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Was et al. 2010 Mar Bio" & msat$MarkerName == "LIST1003"] <- 3 #based on ref paper
msat$Repeat[msat$Source == "Watts et al. 2004 J Sea Res"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Watts et al. 2004 J Sea Res" & msat$MarkerName == "LIST1003"] <- 3 #based on ref paper
  msat$Repeat[msat$Source == "Watts et al. 2004 J Sea Res" & msat$MarkerName == "LIST1007"] <- 3 #based on ref paper
  msat$Repeat[msat$Source == "Watts et al. 2004 J Sea Res" & msat$MarkerName == "MmerUAE01"] <- 3 #based on ref paper
msat$Repeat[msat$Source == "Watts et al. 2009 ICES J Mar Sci"] <- 2 #based on ref papers
  msat$Repeat[msat$Source == "Watts et al. 2009 ICES J Mar Sci" & msat$MarkerName == "LIST1-003"] <- 3 #based on ref paper
  msat$Repeat[msat$Source == "Watts et al. 2009 ICES J Mar Sci" & msat$MarkerName == "LIST1-007"] <- 3 #based on ref paper
msat$Repeat[msat$Source == "White et al. 2009 Mol Ecol"] <- 2 #based on ref papers

######## Check latitude ########

#make sure no instances >90 or <-90
msat[msat$lat > 90 & !is.na(msat$lat), c('spp', 'Source', 'Country', 'Site', 'file', 'lat_deg', 'lat_min', 'lat_sec', 'lat')]	#0
msat[msat$lat < -90 & !is.na(msat$lat), c('spp', 'Source', 'Country', 'Site', 'file', 'lat_deg', 'lat_min', 'lat_sec', 'lat')] #0

#check for missing lat data
sum(is.na(msat$lat) & !is.na(msat$lat_deg)) #anywhere lat is NA but lat_deg is not: 0
inds <- is.na(msat$lat)
sum(inds) # 0

#histogram of latitude
hist(msat$lat) #mostly northern hemisphere, but spans globe

#compare lat to srdbmatch lat, fill in where missing
plot(msat$lat, msat$lat_srdb); abline(0, 1) #falls right on 1:1 line
msat[is.na(msat$lat) & !is.na(msat$lat_srdb), ] #0

####### Check longitude ########

#make sure no instances >180, <-180 or == 0
msat[msat$lon > 180 & !is.na(msat$lon), c('spp', 'Source', 'Country', 'Site', 'lon_deg', 'lon_min', 'lon_sec', 'lon')] #0
msat[msat$lon < -180 & !is.na(msat$lon), c('spp', 'Source', 'Country', 'Site', 'lon_deg', 'lon_min', 'lon_sec', 'lon')] #0
msat[msat$lon == 0 & !is.na(msat$lon), c('spp', 'Source', 'Country', 'Site', 'lon_deg', 'lon_min', 'lon_sec', 'lon')] #1, but OK (off coast of Ghana)

#check for missing lon data
sum(is.na(msat$lon) & !is.na(msat$lon_deg)) #anywhere lon is NA but lon_deg is not: 0
inds <- is.na(msat$lon)
sum(inds) #0

#histogram of longitude
hist(msat$lon) #mostly western hemisphere, but also large chunk in Europe

#compare lon to srdbmatch lon, fill in where missing
plot(msat$lon, msat$lon_srdb); abline(0, 1) #falls right on 1:1 line (ones that do not have improper signs --> fixed below)
msat[is.na(msat$lon) & !is.na(msat$lon_srdb), ] #0

######## Fix lat/long based on msat map ########

#correct lat/long sites where map shows datapoint not in ocean
msat$lon[msat$Source == "Shubina et al. 2004 Env Biol Fish" & msat$Site == "Olyutor"] <- 169.2293 #based on google maps
msat$lat[msat$Source == "Schmidt 2005 Dissertation" & msat$Site == "MAG96 (Irminger Sea)"] <- 60.44 #based on other coordinates at this site
  msat$lon[msat$Source == "Schmidt 2005 Dissertation" & msat$Site == "MAG96 (Irminger Sea)"] <- -28.2 #based on other coordinates at this site
msat$lon[msat$Source == "Zarraonaindia et al. 2009 ICES J Mar Sci" & msat$Site == "LIONS"] <- 4.508017 #based on google maps
msat$lon[msat$Source == "Hoarau et al. 2002 Mol Ecol" & msat$Site == "Bay of Vilaine"] <- -2.782727 #based on google maps
msat$lat[msat$Source == "Larson et al. 2011 Cons Gen 12:679-690" & msat$Site == "Puget Sound, WA"] <- 47.606949 #based on google maps
  msat$lon[msat$Source == "Larson et al. 2011 Cons Gen 12:679-690" & msat$Site == "Puget Sound, WA"] <- -122.407869 #based on google maps
msat$lat[msat$Source == "Anderson & Karel 2010 Fish Aq J" & msat$Site == "Sabine Lake"] <- 29.686301 #based on google maps
  msat$lon[msat$Source == "Anderson & Karel 2010 Fish Aq J" & msat$Site == "Sabine Lake"] <- -93.832946 #based on google maps
msat$lat[msat$Source == "Antoro et al. 2006 Marine Biotechnology 8:17-26" & msat$Site == "Trang"] <- 7.468683 #based on google maps
  msat$lon[msat$Source == "Antoro et al. 2006 Marine Biotechnology 8:17-26" & msat$Site == "Trang"] <- 99.29126 #based on google maps
msat$lat[msat$Source == "Bahri-Sfar et al. 2000 Proceedings of the Royal Society B 267:929-935" & msat$Site == "Northern Tunis lagoon"] <- 36.803542 #based on google maps
  msat$lon[msat$Source == "Bahri-Sfar et al. 2000 Proceedings of the Royal Society B 267:929-935" & msat$Site == "Northern Tunis lagoon"] <- 10.326354 #based on google maps
msat$lat[msat$Source == "Beacham et al. 2001 CSAS; Beacham et al. 2002 CSAS" & msat$Site == "Hunt's Inlet, BC"] <- 54.148063 #based on google maps
  msat$lon[msat$Source == "Beacham et al. 2001 CSAS; Beacham et al. 2002 CSAS" & msat$Site == "Hunt's Inlet, BC"] <- -130.523557 #based on google maps
  msat$lat[msat$Source == "Beacham et al. 2001 CSAS; Beacham et al. 2002 CSAS" & msat$Site == "Louscoone Inlet, BC"] <- 52.12305 #based on google maps
  msat$lon[msat$Source == "Beacham et al. 2001 CSAS; Beacham et al. 2002 CSAS" & msat$Site == "Louscoone Inlet, BC"] <- -131.214599 #based on google maps
  msat$lat[msat$Source == "Beacham et al. 2001 CSAS; Beacham et al. 2002 CSAS" & msat$Site == "Skaat Harbour, BC"] <- 52.40954 #based on google maps
  msat$lon[msat$Source == "Beacham et al. 2001 CSAS; Beacham et al. 2002 CSAS" & msat$Site == "Skaat Harbour, BC"] <- -131.38896 #based on google maps
  msat$lat[msat$Source == "Beacham et al. 2001 CSAS; Beacham et al. 2002 CSAS" & msat$Site == "Tomales Bay, CA"] <- 38.27215 #based on google maps
  msat$lon[msat$Source == "Beacham et al. 2001 CSAS; Beacham et al. 2002 CSAS" & msat$Site == "Tomales Bay, CA"] <- -123.01157 #based on google maps
msat$lat[msat$Source == "Bekkevold et al. 2005 Evolution" & msat$Site == "Kolding Fjord 02"] <- 55.502447 #based on google maps
  msat$lon[msat$Source == "Bekkevold et al. 2005 Evolution" & msat$Site == "Kolding Fjord 02"] <- 9.651237 #based on 
  msat$lat[msat$Source == "Bekkevold et al. 2005 Evolution" & msat$Site == "Limfjord 03"] <- 56.983754 #based on google maps
  msat$lon[msat$Source == "Bekkevold et al. 2005 Evolution" & msat$Site == "Limfjord 03"] <- 10.359377 #based on 
msat$lat[msat$Source == "Bouza et al. 2002 CJFAS" & msat$Site == "Vilagarcia"] <- 42.607834 #based on google maps
  msat$lon[msat$Source == "Bouza et al. 2002 CJFAS" & msat$Site == "Vilagarcia"] <- -8.814009 #based on google maps
msat$lat[msat$Source == "Bradbury et al. 2009 J Fish Bio" & msat$Site == "Holyrood Pond, St. Mary's Bay, and Foxtrap, Conception Bay (Newfoundland)"] <- 47.499221 #based on google maps
  msat$lon[msat$Source == "Bradbury et al. 2009 J Fish Bio" & msat$Site == "Holyrood Pond, St. Mary's Bay, and Foxtrap, Conception Bay (Newfoundland)"] <- -53.055702 #based on google maps 
  msat$lat[msat$Source == "Bradbury et al. 2009 J Fish Bio" & msat$Site == "Holyrood Pond, St. Mary's Bay, and Newman Sound, Bonavista Bay (Newfoundland)"] <- 48.819264 #based on google maps
  msat$lon[msat$Source == "Bradbury et al. 2009 J Fish Bio" & msat$Site == "Holyrood Pond, St. Mary's Bay, and Newman Sound, Bonavista Bay (Newfoundland)"] <-  -53.414021 #based on google maps
msat$lat[msat$Source == "Buonaccorsi et al. 2004 Mar Bio" & msat$Site == "Brookings, OR"] <- 42.100612 #based on google maps
  msat$lon[msat$Source == "Buonaccorsi et al. 2004 Mar Bio" & msat$Site == "Brookings, OR"] <- -124.363406 #based on google maps
  msat$lat[msat$Source == "Buonaccorsi et al. 2004 Mar Bio" & msat$Site == "Hardy Creek, CA"] <- 39.700243 #based on google maps
  msat$lon[msat$Source == "Buonaccorsi et al. 2004 Mar Bio" & msat$Site == "Hardy Creek, CA"] <- -123.849132 #based on google maps
  msat$lat[msat$Source == "Buonaccorsi et al. 2004 Mar Bio" & msat$Site == "San Diego, CA"] <-  32.796148 #based on google maps
  msat$lon[msat$Source == "Buonaccorsi et al. 2004 Mar Bio" & msat$Site == "San Diego, CA"] <- -117.274361 #based on google maps
msat$lat[msat$Source == "Buonaccorsi et al. 2005 Cons Gen" & msat$Site == "Ensenada, MX 2001"] <- 31.813247 #based on google maps
  msat$lon[msat$Source == "Buonaccorsi et al. 2005 Cons Gen" & msat$Site == "Ensenada, MX 2001"] <- -116.73849 #based on google maps
msat$lat[msat$Source == "Burford 2009 J Evol Biol" & msat$Site == "Avila, CA"] <- 35.160262 #based on google maps
  msat$lon[msat$Source == "Burford 2009 J Evol Biol" & msat$Site == "Avila, CA"] <- -120.725973 #based on google maps
  msat$lat[msat$Source == "Burford 2009 J Evol Biol" & msat$Site == "Ft. Ross, CA"] <- 38.495275 #based on google maps
  msat$lon[msat$Source == "Burford 2009 J Evol Biol" & msat$Site == "Ft. Ross, CA"] <- -123.254294 #based on google maps
  msat$lat[msat$Source == "Burford 2009 J Evol Biol" & msat$Site == "Pacific City, OR"] <-  45.220207 #based on google maps
  msat$lon[msat$Source == "Burford 2009 J Evol Biol" & msat$Site == "Pacific City, OR"] <- -124.025876 #based on google maps
msat$lat[msat$Source == "Chapman et al. 2002 Mar Biotech" & msat$Site == "Ashley River, SC"] <- 32.766778 #based on google maps
  msat$lon[msat$Source == "Chapman et al. 2002 Mar Biotech" & msat$Site == "Ashley River, SC"] <- -79.93664 #based on google maps
  msat$lat[msat$Source == "Chapman et al. 2002 Mar Biotech" & msat$Site == "Wando River, SC"] <- 32.807213 #based on google maps
  msat$lon[msat$Source == "Chapman et al. 2002 Mar Biotech" & msat$Site == "Wando River, SC"] <- -79.91383 #based on google maps
msat$lat[msat$Source == "Crivello et al. 2004 J Fish Bio" & msat$Site == "Niantic River, CT"] <- 41.319463 #based on google maps
  msat$lon[msat$Source == "Crivello et al. 2004 J Fish Bio" & msat$Site == "Niantic River, CT"] <- -72.185354 #based on google maps
msat$lat[msat$Source == "D'Anatro, Pereira & Lessa (2011) Environ. Biol. Fish 91:407-420" & msat$Site == "Laguna de Rocha"] <- -34.688445 #based on google maps
  msat$lon[msat$Source == "D'Anatro, Pereira & Lessa (2011) Environ. Biol. Fish 91:407-420" & msat$Site == "Laguna de Rocha"] <- -54.245872 #based on google maps
msat$lat[msat$Source == "Florin & Hoglund 2007 Mol Ecol" & msat$Site == "Gdynia"] <- 54.534548 #based on google maps
  msat$lon[msat$Source == "Florin & Hoglund 2007 Mol Ecol" & msat$Site == "Gdynia"] <- 18.62135 #based on google maps
  msat$lat[msat$Source == "Florin & Hoglund 2007 Mol Ecol" & msat$Site == "Gotland 2004"] <- 57.357056 #based on google maps
  msat$lon[msat$Source == "Florin & Hoglund 2007 Mol Ecol" & msat$Site == "Gotland 2004"] <- 18.763634 #based on google maps
  msat$lat[msat$Source == "Florin & Hoglund 2007 Mol Ecol" & msat$Site == "Latvia"] <- 56.424177 #based on google maps
  msat$lon[msat$Source == "Florin & Hoglund 2007 Mol Ecol" & msat$Site == "Latvia"] <- 20.985281 #based on google maps
msat$lat[msat$Source == "Galarza et al. 2009 CJFAS 66:1478-1490" & msat$Site == "Porticello"] <- 38.132188 #based on google maps
  msat$lon[msat$Source == "Galarza et al. 2009 CJFAS 66:1478-1490" & msat$Site == "Porticello"] <- 13.571875 #based on google maps
msat$lat[msat$Source == "Gold et al. 2010 Fish Res" & msat$Site == "Isla de Margarita"] <- 11.116422 #based on google maps
  msat$lon[msat$Source == "Gold et al. 2010 Fish Res" & msat$Site == "Isla de Margarita"] <- -64.064255 #based on google maps
msat$lat[msat$Source == "Gonzalez-Wanguemert et al. 2010 JEMBE" & msat$Site == "Faro"] <- 36.977872 #based on google maps
  msat$lon[msat$Source == "Gonzalez-Wanguemert et al. 2010 JEMBE" & msat$Site == "Faro"] <- -7.961712 #based on google maps
msat$lat[msat$Source == "Graves & McDowell 2006 Bull Mar Sci 79:469-482" & msat$Site == "Santos"] <- -24.033708 #based on google maps
  msat$lon[msat$Source == "Graves & McDowell 2006 Bull Mar Sci 79:469-482" & msat$Site == "Santos"] <- -46.353149 #based on google maps
msat$lat[msat$Source == "Hess et al. 2001 CJFAS 68:89-104" & msat$Site == "Crescent City, CA"] <- 41.739754 #based on google maps
  msat$lon[msat$Source == "Hess et al. 2001 CJFAS 68:89-104" & msat$Site == "Crescent City, CA"] <- -124.214166 #based on google maps
msat$lat[msat$Source == "Jorgensen et al. 2005 Mol Ecol" & msat$Site == "Bothnian Bay, Finnish side (FS02)"] <- 60.977692 #based on google maps
  msat$lon[msat$Source == "Jorgensen et al. 2005 Mol Ecol" & msat$Site == "Bothnian Bay, Finnish side (FS02)"] <- 20.968159 #based on google maps
  msat$lat[msat$Source == "Jorgensen et al. 2005 Mol Ecol" & msat$Site == "Bothnian Bay, Swedish side (SS02)"] <- 60.592956 #based on google maps
  msat$lon[msat$Source == "Jorgensen et al. 2005 Mol Ecol" & msat$Site == "Bothnian Bay, Swedish side (SS02)"] <-  17.729939 #based on google maps
msat$lat[msat$Source == "Karlsson et al. 2009 Mar Bio" & msat$Site == "Aransas"] <- 27.934264 #based on google maps
  msat$lon[msat$Source == "Karlsson et al. 2009 Mar Bio" & msat$Site == "Aransas"] <- -96.921631 #based on google maps
  msat$lat[msat$Source == "Karlsson et al. 2009 Mar Bio" & msat$Site == "Port Isabel"] <- 26.164982 #based on google maps
  msat$lon[msat$Source == "Karlsson et al. 2009 Mar Bio" & msat$Site == "Port Isabel"] <- -97.237822 #based on google maps
msat$lat[msat$Source == "LeClair et al. 2006 Trans of the Am Fish Soc 135:1631-1643" & msat$Site == "South Puget Sound, WA"] <- 47.177815 #based on google maps
  msat$lon[msat$Source == "LeClair et al. 2006 Trans of the Am Fish Soc 135:1631-1643" & msat$Site == "South Puget Sound, WA"] <- -122.734525 #based on google maps
msat$lat[msat$Source == "Limborg et al. 2009 MEPS" & msat$Site == "Belt Sea (BEL)"] <- 55.569365 #based on google maps
  msat$lon[msat$Source == "Limborg et al. 2009 MEPS" & msat$Site == "Belt Sea (BEL)"] <- 10.541138 #based on google maps
  msat$lat[msat$Source == "Limborg et al. 2009 MEPS" & msat$Site == "Northern Kattegat (KAT)"] <- 57.421479 #based on google maps
  msat$lon[msat$Source == "Limborg et al. 2009 MEPS" & msat$Site == "Northern Kattegat (KAT)"] <- 10.559651 #based on google maps
msat$lat[msat$Source == "Lundy et al. 2000 Mol Ecol" & msat$Site == "Bay of Biscay south 1999"] <- 43.717779 #based on google maps
  msat$lon[msat$Source == "Lundy et al. 2000 Mol Ecol" & msat$Site == "Bay of Biscay south 1999"] <- -7.130774 #based on google maps
msat$lat[msat$Source == "McCraney et al. 2010 Mol Ecol 19:3315-3327" & msat$Site == "Gannon Pond, CA"] <- 40.86991 #based on google maps
  msat$lon[msat$Source == "McCraney et al. 2010 Mol Ecol 19:3315-3327" & msat$Site == "Gannon Pond, CA"] <- -124.175753 #based on google maps
  msat$lat[msat$Source == "McCraney et al. 2010 Mol Ecol 19:3315-3327" & msat$Site == "Gannon Slough, CA"] <- 40.852669 #based on google maps
  msat$lon[msat$Source == "McCraney et al. 2010 Mol Ecol 19:3315-3327" & msat$Site == "Gannon Slough, CA"] <- -124.182241 #based on google maps
  msat$lat[msat$Source == "McCraney et al. 2010 Mol Ecol 19:3315-3327" & msat$Site == "Jacoby Creek, CA"] <- 40.862252 #based on google maps
  msat$lon[msat$Source == "McCraney et al. 2010 Mol Ecol 19:3315-3327" & msat$Site == "Jacoby Creek, CA"] <- -124.183478 #based on google maps
  msat$lat[msat$Source == "McCraney et al. 2010 Mol Ecol 19:3315-3327" & msat$Site == "McDaniel Slough, CA"] <- 40.868353 #based on google maps
  msat$lon[msat$Source == "McCraney et al. 2010 Mol Ecol 19:3315-3327" & msat$Site == "McDaniel Slough, CA"] <- -124.170088 #based on google maps
  msat$lat[msat$Source == "McCraney et al. 2010 Mol Ecol 19:3315-3327" & msat$Site == "Wood Creek, CA"] <- 40.810171 #based on google maps
  msat$lon[msat$Source == "McCraney et al. 2010 Mol Ecol 19:3315-3327" & msat$Site == "Wood Creek, CA"] <- -124.212244 #based on google maps
msat$lat[msat$Source == "McPherson et al. 2001 Journal of Fish Biology 59:356-370" & msat$Site == "Bras d'Or Lakes"] <- 46.371304 #based on google maps
  msat$lon[msat$Source == "McPherson et al. 2001 Journal of Fish Biology 59:356-370" & msat$Site == "Bras d'Or Lakes"] <- -60.276712 #based on google maps
msat$lat[msat$Source == "Mitchell 2006 Thesis" & msat$Site == "Case Inlet, WA 1999"] <- 47.239152 #based on google maps
  msat$lon[msat$Source == "Mitchell 2006 Thesis" & msat$Site == "Case Inlet, WA 1999"] <- -122.663842 #based on google maps
msat$lat[msat$Source == "Mobley et al. 2010 J Biogeography 37:1363-1377" & msat$Site == "Aransas Pass, TX"] <- 28.108717 #based on google mas
  msat$lon[msat$Source == "Mobley et al. 2010 J Biogeography 37:1363-1377" & msat$Site == "Aransas Pass, TX"] <- -96.63383 #based on google maps
  msat$lat[msat$Source == "Mobley et al. 2010 J Biogeography 37:1363-1377" & msat$Site == "Morehead City, NC"] <- 34.61275 #based on google mas
  msat$lon[msat$Source == "Mobley et al. 2010 J Biogeography 37:1363-1377" & msat$Site == "Morehead City, NC"] <- -76.700871 #based on google maps
  msat$lat[msat$Source == "Mobley et al. 2010 J Biogeography 37:1363-1377" & msat$Site == "St. Joseph Bay, FL"] <- 29.785717 #based on google mas
  msat$lon[msat$Source == "Mobley et al. 2010 J Biogeography 37:1363-1377" & msat$Site == "St. Joseph Bay, FL"] <- -85.438702 #based on google maps
msat$lat[msat$Source == "Nielsen et al. 2004 Mol Ecol" & msat$Site == "Belt Sea 2001b"] <- 55.06027 #based on google maps
  msat$lon[msat$Source == "Nielsen et al. 2004 Mol Ecol" & msat$Site == "Belt Sea 2001b"] <- 10.798179 #based on google maps
  msat$lat[msat$Source == "Nielsen et al. 2004 Mol Ecol" & msat$Site == "Central Kattegat 2001"] <- 56.583897 #based on google maps
  msat$lon[msat$Source == "Nielsen et al. 2004 Mol Ecol" & msat$Site == "Central Kattegat 2001"] <- 10.448179 #based on google maps
  msat$lat[msat$Source == "Nielsen et al. 2004 Mol Ecol" & msat$Site == "Northern Kattegat 2001"] <- 57.261885 #based on google maps
  msat$lon[msat$Source == "Nielsen et al. 2004 Mol Ecol" & msat$Site == "Northern Kattegat 2001"] <-  10.607192 #based on google maps
msat$lat[msat$Source == "O'Connell et al. 1998 J Fish Bio" & msat$Site == "Kodiak Island, AK"] <- 58.191321 #based on google maps
  msat$lon[msat$Source == "O'Connell et al. 1998 J Fish Bio" & msat$Site == "Kodiak Island, AK"] <- -153.278154 #based on google maps
msat$lat[msat$Source == "Plank et al. 2010 J Fish Bio 77:329-340" & msat$Site == "Gulf of California, CA"] <- 29.304995 #based on google maps
  msat$lon[msat$Source == "Plank et al. 2010 J Fish Bio 77:329-340" & msat$Site == "Gulf of California, CA"] <- -122.655152 #based on google maps
msat$lat[msat$Source == "Putte et al. 2009 Polar Biol 32:1731-1741" & msat$Site == "Larsen B Ice Shelf"] <- -66.539349 #based on google maps
  msat$lon[msat$Source == "Putte et al. 2009 Polar Biol 32:1731-1741" & msat$Site == "Larsen B Ice Shelf"] <- -66.975285 #based on google maps
msat$lat[msat$Source == "Romo et al. 2006 Fish Sci" & msat$Site == "Nagasaki"] <- 32.726931 #based on google maps
  msat$lon[msat$Source == "Romo et al. 2006 Fish Sci" & msat$Site == "Nagasaki"] <- 129.792108 #based on google maps
msat$lat[msat$Source == "Salas et al. 2009 Mar Biol 157:437-445" & msat$Site == "Cahuita"] <- 9.771318 #based on google maps
  msat$lon[msat$Source == "Salas et al. 2009 Mar Biol 157:437-445" & msat$Site == "Cahuita"] <- -82.83332 #based on google maps
msat$lat[msat$Source == "Sekino & Hara 2001 Mar Biotech" & msat$Site == "Hokkaido"] <- 42.995274 #based on google maps
  msat$lon[msat$Source == "Sekino & Hara 2001 Mar Biotech" & msat$Site == "Hokkaido"] <- 140.268679 #based on google maps
  msat$lat[msat$Source == "Sekino & Hara 2001 Mar Biotech" & msat$Site == "Niigata"] <- 37.5758 #based on google maps
  msat$lon[msat$Source == "Sekino & Hara 2001 Mar Biotech" & msat$Site == "Niigata"] <- 138.631949 #based on google maps
  msat$lat[msat$Source == "Sekino & Hara 2001 Mar Biotech" & msat$Site == "Tottori"] <- 35.549353 #based on google maps
  msat$lon[msat$Source == "Sekino & Hara 2001 Mar Biotech" & msat$Site == "Tottori"] <- 13.77606 #based on google maps
msat$lat[msat$Source == "Sellas et al. 2011 Cons Gen Res 3:609-611" & msat$Site == "Sarasota, FL"] <- 27.341885 #based on google maps
  msat$lon[msat$Source == "Sellas et al. 2011 Cons Gen Res 3:609-611" & msat$Site == "Sarasota, FL"] <- -82.622105 #based on google maps
msat$lat[msat$Source == "Small et al. 2005 TAFS" & msat$Site == "Cherry Point, WA 2003"] <- 48.858482 #based on google maps
  msat$lon[msat$Source == "Small et al. 2005 TAFS" & msat$Site == "Cherry Point, WA 2003"] <- -122.77157 #based on google maps
  msat$lat[msat$Source == "Small et al. 2005 TAFS" & msat$Site == "Squaxin Pass, WA 2002"] <- 47.170636 #based on google maps
  msat$lon[msat$Source == "Small et al. 2005 TAFS" & msat$Site == "Squaxin Pass, WA 2002"] <- -122.637305 #based on google maps
msat$lat[msat$Source == "Teske et al. 2010 Mar Biol 157:2029-2042" & msat$Site == "False Bay"] <- -34.221826 #based on google maps
  msat$lon[msat$Source == "Teske et al. 2010 Mar Biol 157:2029-2042" & msat$Site == "False Bay"] <- 18.650836 #based on google maps
msat$lat[msat$Source == "Tripp-Valdez et al. 2010 Fish Res 105:172-177" & msat$Site == "Topolobampo"] <- 25.551579 #based on google maps
  msat$lon[msat$Source == "Tripp-Valdez et al. 2010 Fish Res 105:172-177" & msat$Site == "Topolobampo"] <- -109.150244 #based on google maps
msat$lat[msat$Source == "Umino et al. 2009 Fish Sci 4:909-919" & msat$Site == "Owase Bay, Mie Prefecture"] <- 34.074843 #based on google maps
  msat$lon[msat$Source == "Umino et al. 2009 Fish Sci 4:909-919" & msat$Site == "Owase Bay, Mie Prefecture"] <- -9.144507 #based on google maps
msat$lat[msat$Source == "Was et al. 2010 Mar Bio" & msat$Site == "Galway Bay"] <- 53.204764 #based on google maps
  msat$lon[msat$Source == "Was et al. 2010 Mar Bio" & msat$Site == "Galway Bay"] <- -9.144507 #based on google maps
msat$lat[msat$Source == "Wilson 2006 Mol Ecol 15:809-824" & msat$Site == "San Diego Bay, CA"] <- 32.676035 #based on google maps
  msat$lon[msat$Source == "Wilson 2006 Mol Ecol 15:809-824" & msat$Site == "San Diego Bay, CA"] <- -117.197796 #based on google maps
   
#correct lat/long sites where not calculated properly originally
msat$lat[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "North Thames Estuary"] <- 52.095 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "North Thames Estuary"] <- 2.04 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "East English Channel"] <- 1.155 #needed to be averaged across collection points
  msat$lat[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "Lyme Bay"] <- 50.44286 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "Lyme Bay"] <- -3.202857 #needed to be averaged across collection points
  msat$lat[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "Carmarthen Bay"] <- 51.634 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "Carmarthen Bay"] <- -4.484 #needed to be averaged across collection points
  msat$lat[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "Tremadog Bay"] <- 52.56 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "Tremadog Bay"] <- -4.383333#needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "Caernarfon Bay"] <- -4.625 #needed to be averaged across collection points
  msat$lat[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "Liverpool Bay"] <- 53.61333 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 J Sea Research 56:305-316" & msat$Site == "Liverpool Bay"] <- -3.426667 #needed to be averaged across collection points
msat$lat[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Ajo"] <- 43.73 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Ajo"] <- -5.29 #needed to be averaged across collection points
  msat$lat[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Penas"] <- 43.47 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Penas"] <- -3.355 #needed to be averaged across collection points
  msat$lat[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Estaca"] <-43.7225  #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Estaca"] <- -6.8925 #needed to be averaged across collection points
  msat$lat[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Corsica, Mediterraean Basin"] <- 41.8825 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Corsica, Mediterraean Basin"] <- 9.5375 #needed to be averaged across collection points
  msat$lat[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Graciosa"] <- 39.10333 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Graciosa"] <- -28.00667 #needed to be averaged across collection points
  msat$lat[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Sao Miguel"] <- 37.7475 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Sao Miguel"] <- -25.6975 #needed to be averaged across collection points
  msat$lat[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Terceira"] <- 38.670711 #needed to be averaged across collection points
  msat$lon[msat$Source == "Chevolot et al. 2006 Mol Ecol 15:3693-3705" & msat$Site == "Terceira"] <- -27.355072 #needed to be averaged across collection points
msat$lat[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Atlantic Iberia"] <- 42.89583 #needed to be averaged across collection points
  msat$lon[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Atlantic Iberia"] <- -9.187495 #needed to be averaged across collection points
  msat$lat[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Bay of Biscay"] <- 43.511111 #needed to be averaged across collection points
  msat$lon[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Bay of Biscay"] <- -2.0888889#needed to be averaged across collection points
  msat$lat[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Porcupine Bank"] <- 50.363333 #needed to be averaged across collection points
  msat$lon[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Porcupine Bank"] <- -13.46 #needed to be averaged across collection points
msat$lat[msat$Source == "Stockley et al. 2005 Mar Bio 146:793-804" & msat$Site == "Azores east"] <- 37.311 #needed to be averaged across collection points
  msat$lon[msat$Source == "Stockley et al. 2005 Mar Bio 146:793-804" & msat$Site == "Azores east"] <- -25.42233 #needed to be averaged across collection points
msat$lat[msat$Source == "Stockley et al. 2005 Mar Bio 146:793-804" & msat$Site == "Azores central"] <- 38.60834 #needed to be averaged across collection points
  msat$lon[msat$Source == "Stockley et al. 2005 Mar Bio 146:793-804" & msat$Site == "Azores central"] <- -28.43334 #needed to be averaged across collection points
  msat$lat[msat$Source == "Stockley et al. 2005 Mar Bio 146:793-804" & msat$Site == "Princess Alice Bank"] <- 38.06665 #needed to be averaged across collection points
  msat$lon[msat$Source == "Stockley et al. 2005 Mar Bio 146:793-804" & msat$Site == "Princess Alice Bank"] <- -29.19167 #needed to be averaged across collection points
msat$lat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$Site == "Eastern Baltic, Bornholm basin"] <- 54.85 #needed to be averaged across collection points
  msat$lon[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$Site == "Eastern Baltic, Bornholm basin"] <- 15.41 #needed to be averaged across collection points 
  msat$lat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$Site == "Eastern Baltic, Gdansk basin"] <- 54.565 #needed to be averaged across collection points
  msat$lon[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$Site == "Eastern Baltic, Gdansk basin"] <- 19.75 #needed to be averaged across collection points 
  msat$lat[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$Site == "Eastern Baltic, Gotland basin"] <- 56.35 #needed to be averaged across collection points
  msat$lon[msat$Source == "Nielsen et al. 2003 Mol Ecol 12:1497-1508" & msat$Site == "Eastern Baltic, Gotland basin"] <- 19.295 #needed to be averaged across collection points 
msat$lat[msat$Source == "Ruzzante et al. 1996 inshore-offshore CJFAS" & msat$Site == "Smith Sound, Trinity Bay"] <- 48.22 #needed to be changed to specific sampling date lat
  msat$lon[msat$Source == "Ruzzante et al. 1996 inshore-offshore CJFAS" & msat$Site == "Smith Sound, Trinity Bay"] <- -53.56 #needed to be changed to specific sampling date lon
  msat$lon[msat$Source == "Ruzzante et al. 1996 inshore-offshore CJFAS" & msat$Site == "SW Arm, Trinity Bay"] <- -53.93 #needed to be changed to specific sampling date lon
msat$lat[msat$Source == "Ruzzante et al. 1998 Mol Ecol" & msat$Site == "Gilbert Bay, Labrador"] <- 52.584 #needed to be averaged across collection points
  msat$lon[msat$Source == "Ruzzante et al. 1998 Mol Ecol" & msat$Site == "Gilbert Bay, Labrador"] <- -55.978 #needed to be averaged across collection points   
msat$lat[msat$Source == "Zatcoff et al. 2004 Mar Bio" & msat$spp == "Epinephelus morio" & msat$Site == "Eastern Gulf of Mexico"] <- 29.40721 #needed to be averaged across collection points
  msat$lon[msat$Source == "Zatcoff et al. 2004 Mar Bio" & msat$spp == "Epinephelus morio" & msat$Site == "Eastern Gulf of Mexico"] <- -84.99995 #needed to be averaged across collection points
  msat$lat[msat$Source == "Zatcoff et al. 2004 Mar Bio" & msat$spp == "Epinephelus morio" & msat$Site == "Mexico"] <- 21.37031 #needed to be averaged across collection points
  msat$lon[msat$Source == "Zatcoff et al. 2004 Mar Bio" & msat$spp == "Epinephelus morio" & msat$Site == "Mexico"] <- -89.44977 #needed to be averaged across collection points
  msat$lat[msat$Source == "Zatcoff et al. 2004 Mar Bio" & msat$spp == "Epinephelus morio" & msat$Site == "North and South Carolina"] <- 33.6043 #needed to be averaged across collection points
  msat$lon[msat$Source == "Zatcoff et al. 2004 Mar Bio" & msat$spp == "Epinephelus morio" & msat$Site == "North and South Carolina"] <- -77.39228 #needed to be averaged across collection points
  msat$lat[msat$Source == "Zatcoff et al. 2004 Mar Bio" & msat$spp == "Mycteroperca phenax" & msat$Site == "Florida east coast"] <- 29.23332 #needed to be changed to specific sampling pt lat
  msat$lon[msat$Source == "Zatcoff et al. 2004 Mar Bio" & msat$spp == "Mycteroperca phenax" & msat$Site == "Florida east coast"] <- -80.92255 #needed to be changed to specific sampling pt lon
  msat$lat[msat$Source == "Zatcoff et al. 2004 Mar Bio" & msat$spp == "Mycteroperca phenax" & msat$Site == "Florida panhandle"] <- 29.6597 #needed to be changed to specific sampling pt lat
  msat$lon[msat$Source == "Zatcoff et al. 2004 Mar Bio" & msat$spp == "Mycteroperca phenax" & msat$Site == "Florida panhandle"] <- -84.70982 #needed to be changed to specific sampling pt lon 
msat$lon[msat$Source == "Knutsen et al. 2007 Mol. Ecol. Notes 7:851-853"] <- 5.58 #typo in paper

######## Check He ########

#make sure no He <0 or >1 (all percentages)
msat[(msat$He<0 | msat$He>1) & !is.na(msat$He),c('spp', 'Source', 'Country', 'Site', 'He')]	#0

#make sure not including any monomorphic loci
dim(msat) #20861 rows
msat <- subset(msat, He != 0)
dim(msat) #20824 rows

#make sure all instances have He
inds <- is.na(msat$He)
sum(inds) #0
msat[inds, ]

#correct He & Hese when reported incorrectly in raw data (often paper has SD and not converted correctly)  
msat$Hese[msat$Source == "White et al. 2010 Heredity 106:690-699" & msat$Site == "NW Mid-Atlantic Ridge" & msat$CollectionYear == 2007] <- 0.0229
  msat$Hese[msat$Source == "White et al. 2010 Heredity 106:690-699" & msat$Site == "NE Mid-Atlantic Ridge" & msat$CollectionYear == 2007] <- 0.0204
  msat$Hese[msat$Source == "White et al. 2010 Heredity 106:690-699" & msat$Site == "NW Mid-Atlantic Ridge" & msat$CollectionYear == 2005] <- 0.0268
  msat$Hese[msat$Source == "White et al. 2010 Heredity 106:690-699" & msat$Site == "SW Mid-Atlantic Ridge" & msat$CollectionYear == 2005] <- 0.0346
  msat$Hese[msat$Source == "White et al. 2010 Heredity 106:690-699" & msat$Site == "NE Mid-Atlantic Ridge" & msat$CollectionYear == 2005] <- 0.0291
  msat$Hese[msat$Source == "White et al. 2010 Heredity 106:690-699" & msat$Site == "SE Mid-Atlantic Ridge" & msat$CollectionYear == 2005] <- 0.0127
  msat$Hese[msat$Source == "White et al. 2010 Heredity 106:690-699" & msat$Site == "Porcupine" & msat$CollectionYear == 2001] <- 0.0209
  msat$Hese[msat$Source == "White et al. 2010 Heredity 106:690-699" & msat$Site == "North of the Azores" & msat$CollectionYear == 2005] <- 0.0225
  msat$Hese[msat$Source == "White et al. 2010 Heredity 106:690-699" & msat$Site == "Greenland" & msat$CollectionYear == 2008] <- 0.0104
  msat$Hese[msat$Source == "White et al. 2010 Heredity 106:690-699" & msat$Site == "SE Mid-Atlantic Ridge" & msat$CollectionYear == 2009] <- 0.0191
  msat$Hese[msat$Source == "White et al. 2010 Heredity 106:690-699" & msat$Site == "NW Mid-Atlantic Ridge" & msat$CollectionYear == 2009] <- 0.0126
msat$Hese[msat$Source == "Teske et al. 2010 Mar Biol 157:2029-2042" & msat$Site == "False Bay"] <- 0.0212
  msat$Hese[msat$Source == "Teske et al. 2010 Mar Biol 157:2029-2042" & msat$Site == "Alphard Banks"] <- 0.0267
  msat$Hese[msat$Source == "Teske et al. 2010 Mar Biol 157:2029-2042" & msat$Site == "Plettenberg Bay"] <- 0.0119
  msat$Hese[msat$Source == "Teske et al. 2010 Mar Biol 157:2029-2042" & msat$Site == "Tsitsikamma NP"] <- 0.0069
  msat$Hese[msat$Source == "Teske et al. 2010 Mar Biol 157:2029-2042" & msat$Site == "Bird Island"] <- 0.0134
msat$He[msat$Source == "Bisol et al. 2007 Hydrobiologia 577:151-159" & msat$Site == "Chioggia"] <- 0.322
  msat$He[msat$Source == "Bisol et al. 2007 Hydrobiologia 577:151-159" & msat$Site == "Lido"] <- 0.307
  msat$He[msat$Source == "Bisol et al. 2007 Hydrobiologia 577:151-159" & msat$Site == "Lago dei Teneri"] <- 0.289
msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Antwerpen"] <- 0.0113
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Saint-Brieuc"] <- 0.0096
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Aveiro"] <- 0.0057
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Rabat"] <- 0.009
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Tahadarte 2"] <-0.0101 
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Tanger 2"] <- 0.0074
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Ksar-Sghir"] <- 0.0074
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Martil"] <- 0.0134
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Nador"] <- 0.0167
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Marseille"] <- 0.0099
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Sete"] <- 0.0098
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Annaba"] <- 0.0137
  msat$Hese[msat$Source == "Naciri et al. 1999 Journal of Heredity 90:592-596" & msat$Site == "Marsala"] <- 0.0157
msat$He[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Atlantic Iberia" & msat$MarkerName == "Mmer-hk3b"] <- 0.844 
  msat$He[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Atlantic Iberia" & msat$MarkerName == "Mmer-hk20b"] <- 0.924
  msat$He[msat$Source == "Pita et al. 2011 Continental Shelf Research 31:376-387" & msat$Site == "Bay of Biscay" & msat$MarkerName == "Mmer-hk3b"] <- 0.819
msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Turneffe Atoll, Belize" & msat$MarkerName == "SpAAC42-1"] <- 0.92
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Turneffe Atoll, Belize" & msat$MarkerName == "SpAAT40-1"] <- 0.88
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Turneffe Atoll, Belize" & msat$MarkerName == "SpAAC44-1"] <- 0.36
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Turneffe Atoll, Belize" & msat$MarkerName == "SpGATA40"] <- 0.95
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Turneffe Atoll, Belize" & msat$MarkerName == "SpAAC33-1"] <- 0.88
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Turneffe Atoll, Belize" & msat$MarkerName == "SpAAC41-1"] <- 0.94
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Banco Chinchorro" & msat$MarkerName == "SpAAC41-1"] <-  0.94 
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Roatan Island" & msat$MarkerName == "SpAAC42-1"] <- 0.92
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Roatan Island" & msat$MarkerName == "SpAAT40-1"] <- 0.86
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Roatan Island" & msat$MarkerName == "SpAAC44-1"] <- 0.37
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Roatan Island" & msat$MarkerName == "SpGATA40"] <- 0.95
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Roatan Island" & msat$MarkerName == "SpAAC33-1"] <- 0.87
  msat$He[msat$Source == "Hepburn et al. 2009 Coral Reefs 28:277-288" & msat$Site == "Roatan Island" & msat$MarkerName == "SpAAC41-1"] <- 0.94
msat <- subset(msat, Source != "Ball et al. 2007 Mar Biol") #removing study and adding per-marker stats (in 304)
msat <- subset(msat, Source != "Burford & Bernardi 2008 Mar Biol") #removing study and adding per-marker stats (in 304)
msat$Hese[msat$Source == "Roques et al. 2002 Mar Bio" & msat$Site == "Saguenay River"] <- 0.008
  msat$Hese[msat$Source == "Roques et al. 2002 Mar Bio" & msat$Site == "St. Lawrence East"] <- 0.009
  msat$Hese[msat$Source == "Roques et al. 2002 Mar Bio" & msat$Site == "St. Lawrence South"] <- 0.009
  msat$Hese[msat$Source == "Roques et al. 2002 Mar Bio" & msat$Site == "Traenaegga"] <- 0.012

#remove sites and/or loci not in HWE
msat <- subset(msat, Source != "Antoro et al. 2006 Marine Biotechnology 8:17-26" | Site != "Nakornsrithammarat" | MarkerName != "Em-07")
  msat <- subset(msat, Source != "Antoro et al. 2006 Marine Biotechnology 8:17-26" | Site != "Trang" | MarkerName != "Em-07")
  msat <- subset(msat, Source != "Antoro et al. 2006 Marine Biotechnology 8:17-26" | Site != "Lampung" | MarkerName != "Em-07")
  msat <- subset(msat, Source != "Antoro et al. 2006 Marine Biotechnology 8:17-26" | Site != "Jepara" | MarkerName != "Em-07" & MarkerName != "Em-08" & MarkerName != "Em-10")
  msat <- subset(msat, Source != "Antoro et al. 2006 Marine Biotechnology 8:17-26" | Site != "Flores" | MarkerName != "Em-01" & MarkerName != "Em-07" & MarkerName != "Em-08")
msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Kattegat 03" | MarkerName != "Cha1020")
  msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Lillebaelt 03" | MarkerName != "Cpa101")
  msat <- subset(msat, Source != "Bekkevold et al. 2005 Evolution" | Site != "Tjome 03" | MarkerName != "Cha1202")
msat <- subset(msat, Source != "Chevolot et al. 2006 J Sea Research 56:305-316" | Site != "Outer Wash" | MarkerName != "Rc-B4")
  msat <- subset(msat, Source != "Chevolot et al. 2006 J Sea Research 56:305-316" | Site != "North Thames Estuary" | MarkerName != "Rc-B4" & MarkerName != "Rc-E9")
  msat <- subset(msat, Source != "Chevolot et al. 2006 J Sea Research 56:305-316" | Site != "Lyme Bay" | MarkerName != "Rc-B4")
  msat <- subset(msat, Source != "Chevolot et al. 2006 J Sea Research 56:305-316" | Site != "Carmarthen Bay" | MarkerName != "Rc-B3")
  msat <- subset(msat, Source != "Chevolot et al. 2006 J Sea Research 56:305-316" | Site != "Tremadog Bay" | MarkerName != "Rc-B4")
  msat <- subset(msat, Source != "Chevolot et al. 2006 J Sea Research 56:305-316" | Site != "Liverpool Bay" | MarkerName != "Rc-B3" & MarkerName != "Rc-B4")
msat <- subset(msat, Source != "Chevolot et al. 2006 Mol Ecol 15:3693-3705" | Site != "Ajo" | MarkerName != "Rc-E9" & MarkerName != "Rc-G2")
  msat <- subset(msat, Source != "Chevolot et al. 2006 Mol Ecol 15:3693-3705" | Site != "Penas" | MarkerName != "Rc-B3")
  msat <- subset(msat, Source != "Chevolot et al. 2006 Mol Ecol 15:3693-3705" | Site != "Corsica, Mediterraean Basin" | MarkerName != "Rc-B6")
  msat <- subset(msat, Source != "Chevolot et al. 2006 Mol Ecol 15:3693-3705" | Site != "Sao Miguel" | MarkerName != "Rc-B4")
msat <- subset(msat, Source != "Hauser et al. 2006 Intl Pacific Halibut Commission Sci Report 81:1-28" | MarkerName != "Hhi59" & MarkerName != "HhiJ42" & MarkerName != "HhiD34")
msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Haemulon flavolineatum" | Site != "Cayos Cochinos" | MarkerName != "AAC3" & MarkerName != "AAC54")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Haemulon flavolineatum" | Site != "Glovers Reef" | MarkerName != "AAC46")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Haemulon flavolineatum" | Site != "Lower Keys, FL" | MarkerName != "AAC43")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Haemulon flavolineatum" | Site != "Lee Stocking Island, FL" | MarkerName != "AAT3")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Haemulon flavolineatum" | Site != "Mahahual" | MarkerName != "AAC37")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Haemulon flavolineatum" | Site != "Puerto Morelos" | MarkerName != "AAC41" & MarkerName != "AAC46")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Haemulon flavolineatum" | Site != "St. Lucia" | MarkerName != "AAC3" & MarkerName != "AAC54")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Haemulon flavolineatum" | Site != "Tobago" | MarkerName != "AAC10" & MarkerName != "AAC41" & MarkerName != "AAC43")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Haemulon flavolineatum" | Site != "Upper Keys, FL" | MarkerName != "AAC41" & MarkerName != "AAC46" & MarkerName != "AAT3")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Thalassoma bifasciatum" | Site != "Anguilla" | MarkerName != "AAT42")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Thalassoma bifasciatum" | Site != "Aruba" | MarkerName != "Cap")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Thalassoma bifasciatum" | Site != "Barbados" | MarkerName != "AAT4")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Thalassoma bifasciatum" | Site != "Cayos Cochinos" | MarkerName != "AAT42")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Thalassoma bifasciatum" | Site != "Dominica" | MarkerName != "AAC34")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Thalassoma bifasciatum" | Site != "Glovers Reef" | MarkerName != "AAT4" & MarkerName != "AAT42")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Thalassoma bifasciatum" | Site != "Lower Keys, FL" | MarkerName != "AAT42")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Thalassoma bifasciatum" | Site != "Puerto Rico" | MarkerName != "AAC34" & MarkerName != "AAT42")
  msat <- subset(msat, Source != "Purcell et al. 2006 Proc R Soc B 273:1483-90" | spp != "Thalassoma bifasciatum" | Site != "Anguilla" | MarkerName != "AAT4" & MarkerName != "AAT42")
msat <- subset(msat, Source != "Wilson 2006 Mol Ecol 15:809-824" | spp != "Syngnathus auliscus" | MarkerName != "Slep3")
msat <- subset(msat, Source != "Landinez-Garcia et al. 2009 Ciencias Marinas 35:321-331") #all loci out of HWE
msat <- subset(msat, Source != "Leray et al. 2010 Evolution 64:1218-1230" | Site != "French Polynesia" | MarkerName != "A114" & MarkerName != "C12")
  msat <- subset(msat, Source != "Leray et al. 2010 Evolution 64:1218-1230" | Site != "Marquesas" | MarkerName != "A7" & MarkerName != "C12" & MarkerName != "A101")
msat <- subset(msat, Source != "Liu et al. 2011 Acta Oceanologica Sinica 30:76-83" | Site != "Laizhou Bay" | MarkerName != "HSTS_2" & MarkerName != "HSTS_4")
  msat <- subset(msat, Source != "Liu et al. 2011 Acta Oceanologica Sinica 30:76-83" | Site != "Weihai" | MarkerName != "HSTS_4" & MarkerName != "HSTS_7")
  msat <- subset(msat, Source != "Liu et al. 2011 Acta Oceanologica Sinica 30:76-83" | Site != "Qingdao" | MarkerName != "HSTS_c" & MarkerName != "HSTS_4" & MarkerName != "HSTS_7")
  msat <- subset(msat, Source != "Liu et al. 2011 Acta Oceanologica Sinica 30:76-83" | Site != "Rizhao" | MarkerName != "HSTS_h" & MarkerName != "HSTS_i")
msat <- subset(msat, Source != "Ma et al. 2011 Gen & Mol Res 10:1455-1460" | MarkerName != "Lap2" & MarkerName != "Lap5" & MarkerName != "Lap10")
msat <- subset(msat, Source != "Ma et al. 2011 Fish Sci 77:707-711" | MarkerName != "Niba9" & MarkerName != "Niba10" & MarkerName != "Niba13")
msat <- subset(msat, Source != "Dammannagoda et al. 2011 CJFAS 68:210-223" | Site != "Negombo" | MarkerName != "UTD328")
  msat <- subset(msat, Source != "Dammannagoda et al. 2011 CJFAS 68:210-223" | Site != "Weligama" | MarkerName != "UTD535" & MarkerName != "UTD523" & MarkerName != "UTD172" & MarkerName != "UTD328")
  msat <- subset(msat, Source != "Dammannagoda et al. 2011 CJFAS 68:210-223" | Site != "Tangalle" | MarkerName != "UTD172" & MarkerName != "UTD328" & MarkerName != "UTD73")
  msat <- subset(msat, Source != "Dammannagoda et al. 2011 CJFAS 68:210-223" | Site != "Kalmunai" | MarkerName != "UTD172" & MarkerName != "UTD328" & MarkerName != "UTD73")
  msat <- subset(msat, Source != "Dammannagoda et al. 2011 CJFAS 68:210-223" | Site != "Trincomalee" | MarkerName != "UTD523" & MarkerName != "UTD328" & MarkerName != "UTD73")
  msat <- subset(msat, Source != "Dammannagoda et al. 2011 CJFAS 68:210-223" | Site != "Laccadive Islands" | MarkerName != "UTD535" & MarkerName != "UTD523" & MarkerName != "UTD172" & MarkerName != "UTD328")
  msat <- subset(msat, Source != "Dammannagoda et al. 2011 CJFAS 68:210-223" | Site != "Maldives" | MarkerName != "UTD523" & MarkerName != "UTD172" & MarkerName != "UTD328")
msat <- subset(msat, Source != "Muths & Bourjea 2011 Cons Gen Res 3:629-631" | Site != "Moheli" | MarkerName != "EPI-04" & MarkerName != "EPI-48")
msat <- subset(msat, Source != "Nance et al. 2011 PLoS ONE 6:e21459" | Site != "Mazatlan" | MarkerName != "SLE053")
msat <- subset(msat, Source != "Nance et al. 2011 PLoS ONE 6:e21459" | Site != "Santa Catalina" | MarkerName != "SLE053" & MarkerName != "SLE071" & MarkerName != "SLE086")
  msat <- subset(msat, Source != "Pita et al. 2011 Continental Shelf Research 31:376-387" | Site != "Bay of Biscay" | MarkerName != "Mmer-hk20b")
msat <- subset(msat, Source != "Putte et al. 2009 Polar Biol 32:1731-1741" | spp != "Trematomus hansoni" | MarkerName != "Trne055")
msat <- subset(msat, Source != "Saillant et al. 2010 ICES J Mar Sci 67:1240-1250" | Site != "Northwest Florida" | MarkerName != "Prs221")
  msat <- subset(msat, Source != "Saillant et al. 2010 ICES J Mar Sci 67:1240-1250" | Site != "Mississippi/Alabama" | MarkerName != "Lca91")
msat <- subset(msat, Source != "Salas et al. 2009 Mar Biol 157:437-445" | Site != "Barrier North" | MarkerName != "SpGATA40"& MarkerName != "SpTG16")
  msat <- subset(msat, Source != "Salas et al. 2009 Mar Biol 157:437-445" | Site != "Barrier Central" | MarkerName != "SpGATA40")
  msat <- subset(msat, Source != "Salas et al. 2009 Mar Biol 157:437-445" | Site != "Uvita Island" | MarkerName != "SpGATA40")
msat <- subset(msat, Source != "Tripp-Valdez et al. 2010 Fish Res 105:172-177" | Site != "Loreto" | MarkerName != "chi002")
  msat <- subset(msat, Source != "Tripp-Valdez et al. 2010 Fish Res 105:172-177" | Site != "Guaymas" | MarkerName != "chi008a")
  msat <- subset(msat, Source != "Tripp-Valdez et al. 2010 Fish Res 105:172-177" | Site != "Puerto Libertad" | MarkerName != "chi008")
  msat <- subset(msat, Source != "Tripp-Valdez et al. 2010 Fish Res 105:172-177" | Site != "Cabo San Lucas" | MarkerName != "chi023")
  msat <- subset(msat, Source != "Tripp-Valdez et al. 2010 Fish Res 105:172-177" | Site != "Mazatlan" | MarkerName != "chi002")
  msat <- subset(msat, Source != "Tripp-Valdez et al. 2010 Fish Res 105:172-177" | Site != "Nayarit" | MarkerName != "chi002")
  msat <- subset(msat, Source != "Tripp-Valdez et al. 2010 Fish Res 105:172-177" | Site != "Punta Lobos" | MarkerName != "chi037")
  msat <- subset(msat, Source != "Tripp-Valdez et al. 2010 Fish Res 105:172-177" | Site != "Pacific Ocean" | MarkerName != "chi023" & MarkerName != "chi037")
msat <- subset(msat, Source != "Maggio et al. 2009 ICES J Mar Sci 66:1883-1891" | Site != "Sete, Gulf of Lions" | MarkerName != "Mb7" & MarkerName != "Mb39")
  msat <- subset(msat, Source != "Maggio et al. 2009 ICES J Mar Sci 66:1883-1891" | Site != "Genova, Tyrrhenian Sea" | MarkerName != "Mb26" & MarkerName != "Mb39")
  msat <- subset(msat, Source != "Maggio et al. 2009 ICES J Mar Sci 66:1883-1891" | Site != "S. Agata di Militello, Tyrrhenian Sea" | MarkerName != "Mb15" & MarkerName != "Mb26" & MarkerName != "Mb31" & MarkerName != "Mb39")
  msat <- subset(msat, Source != "Maggio et al. 2009 ICES J Mar Sci 66:1883-1891" | Site != "Porticello, Tyrrhenian Sea" | MarkerName != "Mb39")
  msat <- subset(msat, Source != "Maggio et al. 2009 ICES J Mar Sci 66:1883-1891" | Site != "Castellamare del Golfo, Tyrrhenian Sea" | MarkerName != "Mb26" & MarkerName != "Mb39")
  msat <- subset(msat, Source != "Maggio et al. 2009 ICES J Mar Sci 66:1883-1891" | Site != "Licata, Strait of Sicily" | MarkerName != "Mb15" & MarkerName != "Mb39")
  msat <- subset(msat, Source != "Maggio et al. 2009 ICES J Mar Sci 66:1883-1891" | Site != "Siarcusa, Strait of Sicily" | MarkerName != "Mb7" & MarkerName != "Mb15" & MarkerName != "Mb39")
msat <- subset(msat, Source != "Christie & Eble 2009 Mol. Ecol. Res. 9:544-546" | MarkerName != "Zefl03" & MarkerName != "Zefl04" & MarkerName != "Zefl05" & MarkerName != "Zefl11" & MarkerName != "Zefl16")
msat <- subset(msat, Source != "Hepburn et al. 2009 Coral Reefs 28:277-288" | Site != "Turneffe Atoll, Belize"| MarkerName != "SpAAT9-2" & MarkerName != "SpAAT39")
  msat <- subset(msat, Source != "Hepburn et al. 2009 Coral Reefs 28:277-288" | Site != "Banco Chinchorro"| MarkerName != "SpGATA40" & MarkerName != "SpAAT9-2" & MarkerName != "SpAAT39")
  msat <- subset(msat, Source != "Hepburn et al. 2009 Coral Reefs 28:277-288" | Site != "Roatan Island"| MarkerName != "SpAAT9-2" & MarkerName != "SpAAT39")
msat <- subset(msat, Source != "Matschiner, Hanel & Salzburger 2009 Mol. Ecol. Res. 18:2574-2587" | MarkerName != "Trne35")
msat <- subset(msat, Source != "Miller-Sims et al. 2008: Mol. Ecol. 17:5036-5048" | MarkerName != "AC33?" & MarkerName != "POM3?" & MarkerName != "POM15?")
  msat <- subset(msat, Source != "Miller-Sims et al. 2008: Mol. Ecol. 17:5036-5048" | Site != "Sykes reef Capricorn-Bunker group of Southern Great Barrier Reef" | MarkerName != "POM3")
msat <- subset(msat, Source != "Buston et al. 2007 Mol. Ecol. Res. 16:3671-3678" | MarkerName != "Cf4")
msat <- subset(msat, Source != "Knutsen et al. 2007 Mol. Ecol. Notes 7:851-853" | MarkerName != "Bbrom 8" & MarkerName != "Bbrom 11")
msat <- subset(msat, Source != "Papetti et al. 2007 Polar Biol 30:1605-1613" | Site != "Elephant Island" | MarkerName != "Ca40")
  msat <- subset(msat, Source != "Papetti et al. 2007 Polar Biol 30:1605-1613" | Site != "South Shetlands" | MarkerName != "Ca55")
msat <- subset(msat, Source != "Galarza et al. 2007 Molecular Ecology Notes 7:230-232" | spp != "Mullus surmuletus" | MarkerName != "Mbar14" & MarkerName != "Mbar46" & MarkerName != "Mbar55")
msat <- subset(msat, Source != "Miller-Sims et al. 2005. Mol Ecol Notes 5:424-426 " | MarkerName != "POM 6")
msat <- subset(msat, Source != "Stockley et al. 2005 Mar Bio 146:793-804" | Site != "Azores west" | MarkerName != "PbMS1" & MarkerName != "PbMS2")
  msat <- subset(msat, Source != "Stockley et al. 2005 Mar Bio 146:793-804" | Site != "Azores east" | MarkerName != "PbMS1" & MarkerName != "PbMS2")
  msat <- subset(msat, Source != "Stockley et al. 2005 Mar Bio 146:793-804" | Site != "Azores central" | MarkerName != "PbMS1" & MarkerName != "PbMS2")
  msat <- subset(msat, Source != "Stockley et al. 2005 Mar Bio 146:793-804" | Site != "Princess Alice Bank" | MarkerName != "PbMS1" & MarkerName != "PbMS2" & MarkerName != "PbMS6" & MarkerName != "PbMS20")
  msat <- subset(msat, Source != "Stockley et al. 2005 Mar Bio 146:793-804" | Site != "Portuguese slope" | MarkerName != "PbMS1" & MarkerName != "PbMS2")
msat <- subset(msat, Source != "Garoia et al. 2004. Marine Biotechnology. 6:446-452." | Site != "Ortona" | MarkerName != "Mb39")
  msat <- subset(msat, Source != "Garoia et al. 2004. Marine Biotechnology. 6:446-452." | Site != "Bari" | MarkerName != "Mb15" & MarkerName != "Mb26b")
msat <- subset(msat, Source != "Hoffman et al. 2004. Mol Ecol Notes 4:342-344." | MarkerName != "Pka16")
msat <- subset(msat, Source != "Williams et al. 2004. Mol Ecol Notes 4:525-527" | MarkerName != "TbAAT4-B" & MarkerName != "TbAAT18-B" & MarkerName != "TbAAC50-C" & MarkerName != "TbAAT42-D" & MarkerName != "T323modt-A")
msat <- subset(msat, Source != "Van Herwerden, Benzie & Davies. 2003. Journal of Fish Biology. 62:987-999." | Site != "Sandshoe, Swains" | MarkerName != "58rte" & MarkerName != "67rte")
  msat <- subset(msat, Source != "Van Herwerden, Benzie & Davies. 2003. Journal of Fish Biology. 62:987-999." | Site != "20-137, Mackay" | MarkerName != "58rte" & MarkerName != "67rte" & MarkerName != "80rte" & MarkerName != "90rte")
  msat <- subset(msat, Source != "Van Herwerden, Benzie & Davies. 2003. Journal of Fish Biology. 62:987-999." | Site != "Glow, Townsville" | MarkerName != "67rte" & MarkerName != "80rte")
msat <- subset(msat, Source != "Aboim et al. 2003. Mol Ecol Notes 3:18-20" | MarkerName != "Hd 008" & MarkerName != "Hd 044" & MarkerName != "Hd 095")
msat <- subset(msat, Source != "Garcia de Leon et al. 1997 Molecular Ecology 6:51-62" | Site != "Valencia" | MarkerName != "Labrax-8" & MarkerName != "Labrax-29")
  msat <- subset(msat, Source != "Garcia de Leon et al. 1997 Molecular Ecology 6:51-62" | Site != "Eastern Rhone" | MarkerName != "Labrax-6" & MarkerName != "Labrax-8")
msat <- subset(msat, Source != "De Innocentiis et al. 2001 Molecular Ecology 10:2163-2175" | Site != "Thyrrenian Sea" | MarkerName != "GAG045" & MarkerName != "GAG049")
  msat <- subset(msat, Source != "De Innocentiis et al. 2001 Molecular Ecology 10:2163-2175" | Site != "Ionian Sea" | MarkerName != "GAG038")
  msat <- subset(msat, Source != "De Innocentiis et al. 2001 Molecular Ecology 10:2163-2175" | Site != "Lampedusa Island" | MarkerName != "GAG007")
  msat <- subset(msat, Source != "De Innocentiis et al. 2001 Molecular Ecology 10:2163-2175" | Site != "South Mediterranean Sea" | MarkerName != "GAG010" & MarkerName != "GAG049")
msat <- subset(msat, Source != "Smith&McVeagh 2000 Journal of Fish Biology 57:72-83" | Site != "Macquarie Island" | MarkerName != "To3")
msat <- subset(msat, Source != "Shaw et al. 1999 Heredity 83:490-499" | Site != "Barents Sea (NS1)" | MarkerName != "Cha20" & MarkerName != "Cha113")
  msat <- subset(msat, Source != "Shaw et al. 1999 Heredity 83:490-499" | Site != "Norwegian Sea (NS2)" | MarkerName != "Cha113")
msat <- subset(msat, Source != "McGowan&Reith 1999 Molecular Ecology 8:1761-1763" | MarkerName != "HhiJ42")
msat <- subset(msat, Source != "McPherson et al. 2001 Journal of Fish Biology 59:356-370" | Site != "Eastern Passage" | MarkerName != "Cha1045")
msat <- subset(msat, Source != "Nugruho&Taniguchi 1999 Fisheries Science 65:353-357") #excluding bc didn't test for HWE
msat <- subset(msat, Source != "McPherson et al. 2001 Molecular Ecology Notes 1: 31-32" | MarkerName != "1014" & MarkerName != "1235")
msat <- subset(msat, Source != "An et al. 2009 Cons Gen" | MarkerName != "KSs19B" & MarkerName != "KSs5")
msat <- subset(msat, Source != "An et al. 2009 Genes Genomics" | MarkerName != "KSi221B")
msat <- subset(msat, Source != "Anderson & Karel 2007 J Fish Bio" | Site != "Galveston, TX" & Site != "New Jersey")
msat <- subset(msat, Source != "Anderson & McDonald 2007 J Fish Bio" | MarkerName != "AF039660")
msat <- subset(msat, Source != "Bentzen et al. 1996 CJFAS" | Site != "Flemish Cap" | MarkerName != "Gmo141")
  msat <- subset(msat, Source != "Bentzen et al. 1996 CJFAS" | Site != "Northeast Spur" | MarkerName != "Gmo141")
  msat <- subset(msat, Source != "Bentzen et al. 1996 CJFAS" | Site != "Scotian Shelf" | MarkerName != "Gmo141")
  msat <- subset(msat, Source != "Bentzen et al. 1996 CJFAS" | Site != "SOUTH (North Cape, Grand Bank, Nose of the Bank)" | MarkerName != "Gmo141")
msat <- subset(msat, Source != "Bernal-Ramirez et al. 2003 Mar Bio" | Site != "Hawke Bay" | MarkerName != "GA2B")
msat <- subset(msat, Source != "Bradbury et al. 2009 J Fish Bio" | MarkerName != "Ute35")
msat <- subset(msat, Source != "Buonaccorsi et al. 2005 Cons Gen" | Site != "San Diego, CA")
msat <- subset(msat, Source != "Burford 2009 J Evol Biol" | Site != "Newport/Depoe, OR" | MarkerName != "Sra.7-25")
  msat <- subset(msat, Source != "Burford 2009 J Evol Biol" | Site != "Ocean Cove, CA" | MarkerName != "Sra.7-7")
  msat <- subset(msat, Source != "Burford 2009 J Evol Biol" | Site != "Pacific City, OR" | MarkerName != "Sra.7-25")
  msat <- subset(msat, Source != "Burford 2009 J Evol Biol" | Site != "Pt. Sur, CA" | MarkerName != "Sra.7-7")
  msat <- subset(msat, Source != "Burford 2009 J Evol Biol" | Site != "Santa Rosa Is., CA" | MarkerName != "Sra.7-2")
msat <- subset(msat, Source != "Burridge & Smolenski 2003 NZ J Mar Fresh Res" | Site != "Eden" | MarkerName != "Nma305" & MarkerName != "Nma311")
  msat <- subset(msat, Source != "Burridge & Smolenski 2003 NZ J Mar Fresh Res" | Site != "New Zealand" | MarkerName != "Nma106")
  msat <- subset(msat, Source != "Burridge & Smolenski 2003 NZ J Mar Fresh Res" | Site != "Tasman" | MarkerName != "Nma106" & MarkerName != "Nma311")
  msat <- subset(msat, Source != "Burridge & Smolenski 2003 NZ J Mar Fresh Res" | Site != "Western Australia" | MarkerName != "Nma106" & MarkerName != "Nma305")
msat <- subset(msat, Source != "Canales-Aguirre et al. 2010 Cons Gen" | MarkerName != "TmurA115" & MarkerName != "TmurB104")
msat <- subset(msat, Source != "Cardenas et al. 2009 Fish Res" | Site != "Iquique" | MarkerName != "Tt74")
  msat <- subset(msat, Source != "Cardenas et al. 2009 Fish Res" | Site != "New Zealand" | MarkerName != "Tt74")
msat <- subset(msat, Source != "Carreras-Carbonell et al. 2006 Mol Ecol Notes" | MarkerName != "Sc05" & MarkerName != "Sc08")
msat <- subset(msat, Source != "Carreras-Carbonell et al. 2008 Cons Gen") #excluding bc didn't test for HWE (and only reported He)
msat <- subset(msat, Source != "Castillo et al. 2005 ICES J Mar Sci" | Site != "IX-C" & Site != "IX-S" & Site != "Med" & Site != "VIIIc-E" & Site != "VIIIc-W")
msat <- subset(msat, Source != "Catanese et al. 2008 Mol Ecol Res" | Site != "Atlantic (Gulf of Cadiz)" | MarkerName != "Maz1" & MarkerName != "Maz4" & MarkerName != "Maz7" & MarkerName != "Maz9")
  msat <- subset(msat, Source != "Catanese et al. 2008 Mol Ecol Res" | Site != "Mediterranean (Malaga and Estepona)" | MarkerName != "Maz2" & MarkerName != "Maz7" & MarkerName != "Maz9")
msat <- subset(msat, Source != "Cha et al. 2009 Cons Gen Res" | MarkerName != "KSj18")
msat <- subset(msat, Source != "Chang et al. 2009 Gen Genom" | MarkerName != "KTj13" & MarkerName != "KTj20" & MarkerName != "KTj28")
msat <- subset(msat, Source != "Chapman et al. 2002 Mar Biotech" | MarkerName != "RD201")
  msat <- subset(msat, Source != "Chapman et al. 2002 Mar Biotech" | Site != "Ace Basin, SC" | MarkerName != "Soc014")
  msat <- subset(msat, Source != "Chapman et al. 2002 Mar Biotech" | Site != "Ashley River, SC" | MarkerName != "Soc014")
  msat <- subset(msat, Source != "Chapman et al. 2002 Mar Biotech" | Site != "Broad River, SC" | MarkerName != "Soc014")
  msat <- subset(msat, Source != "Chapman et al. 2002 Mar Biotech" | Site != "Cape Romain, SC" | MarkerName != "Soc029")
  msat <- subset(msat, Source != "Chapman et al. 2002 Mar Biotech" | Site != "Charleston Harbor, SC" | MarkerName != "Soc014" & MarkerName != "Soc029")
  msat <- subset(msat, Source != "Chapman et al. 2002 Mar Biotech" | Site != "Colleton River, SC" | MarkerName != "Soc014")
  msat <- subset(msat, Source != "Chapman et al. 2002 Mar Biotech" | Site != "Sapelo Sound, SC" | MarkerName != "Soc014")
  msat <- subset(msat, Source != "Chapman et al. 2002 Mar Biotech" | Site != "Wando River, SC" | MarkerName != "Soc014")
msat <- subset(msat, Source != "Charrier et al. 2007 MEPS" | Site != "Bay of Audierne" | MarkerName != "Tch8")
  msat <- subset(msat, Source != "Charrier et al. 2007 MEPS" | Site != "Cape Breton" | MarkerName != "Gmo02" & MarkerName != "Tch10")
  msat <- subset(msat, Source != "Charrier et al. 2007 MEPS" | Site != "Celtic Sea" | MarkerName != "Gmo02")
  msat <- subset(msat, Source != "Charrier et al. 2007 MEPS" | Site != "Dingle" | MarkerName != "Gmo02")
  msat <- subset(msat, Source != "Charrier et al. 2007 MEPS" | Site != "Dogger Bank" | MarkerName != "Gmo02" & MarkerName != "Tch10")
  msat <- subset(msat, Source != "Charrier et al. 2007 MEPS" | Site != "Flamborough Head" | MarkerName != "Pop14c")
  msat <- subset(msat, Source != "Charrier et al. 2007 MEPS" | Site != "Irish Sea" | MarkerName != "Gmo02" & MarkerName != "Tch10")
  msat <- subset(msat, Source != "Charrier et al. 2007 MEPS" | Site != "Skarvoyni" | MarkerName != "Gmo02" & MarkerName != "Tch10")
  msat <- subset(msat, Source != "Charrier et al. 2007 MEPS" | Site != "Southern Bight" | MarkerName != "Gmo02" & MarkerName != "Tch8")
  msat <- subset(msat, Source != "Charrier et al. 2007 MEPS" | Site != "Western Hebridges" | MarkerName != "Gmo02" & MarkerName != "Tch10")
msat <- subset(msat, Source != "Crivello et al. 2004 J Fish Bio" | MarkerName != "A441" & MarkerName != "D34" & MarkerName != "J42" & MarkerName != "P157")
  msat <- subset(msat, Source != "Crivello et al. 2004 J Fish Bio" | Site != "Thames River, CT" | MarkerName != "I29")
msat <- subset(msat, Source != "Cushman et al. 2009 CJFAS" | Site != "East coast Florida" | MarkerName != "MBO088")
msat <- subset(msat, Source != "D'Amato 2006 Mar Biotech" | Site != 3 | MarkerName != "Mm110-13" & MarkerName != "Mm110-8" & MarkerName != "19-1g2" & MarkerName != "Mm422" & MarkerName != "Mm5-4")
  msat <- subset(msat, Source != "D'Amato 2006 Mar Biotech" | Site != 4 | MarkerName != "19-1g2" & MarkerName != "Mm5-4")
  msat <- subset(msat, Source != "D'Amato 2006 Mar Biotech" | Site != 5 | MarkerName != "19-1g2")
  msat <- subset(msat, Source != "D'Amato 2006 Mar Biotech" | Site != 6 | MarkerName != "19-1g2" & MarkerName != "Mm422" & MarkerName != "Mm5-4")
  msat <- subset(msat, Source != "D'Amato 2006 Mar Biotech" | Site != 7 | MarkerName != "19-1g2" & MarkerName != "Mm422" & MarkerName != "Mm5-4")
  msat <- subset(msat, Source != "D'Amato 2006 Mar Biotech" | Site != 8 | MarkerName != "19-1g2" & MarkerName != "Mm5-4")
  msat <- subset(msat, Source != "D'Amato 2006 Mar Biotech" | Site != 9 | MarkerName != "Mm110-8" & MarkerName != "19-1g2")
msat <- subset(msat, Source != "D'Amato et al. 1999 Mol Ecol" | MarkerName != "Mm 19-1 g-1" & MarkerName != "Mm 19-1 g-2")
msat <- subset(msat, Source != "Dammannagoda 2007 Thesis" | Site != "Kalmunei" | MarkerName != "UTD328" & MarkerName != "UTD73")
  msat <- subset(msat, Source != "Dammannagoda 2007 Thesis" | Site != "Maldives" | MarkerName != "UTD328")
  msat <- subset(msat, Source != "Dammannagoda 2007 Thesis" | Site != "Negombo 2001" | MarkerName != "UTD328")
  msat <- subset(msat, Source != "Dammannagoda 2007 Thesis" | Site != "Weligama 2001" | MarkerName != "UTD328")
msat <- subset(msat, Source != "Dammannagoda et al. 2008 Fish Res" | Site != "Kirinda" | MarkerName != "UTD402")
  msat <- subset(msat, Source != "Dammannagoda et al. 2008 Fish Res" | Site != "Negombo" | MarkerName != "UTD402" & MarkerName != "UTD494")
msat <- subset(msat, Source != "Danancher & Garcia-Vazquez 2009 Mar Bio" | spp != "Lepidorhombus whiffiagonis" | Site != "VIIIabd" | MarkerName != "Loc-34")
  msat <- subset(msat, Source != "Danancher & Garcia-Vazquez 2009 Mar Bio" | spp != "Lepidorhombus whiffiagonis" | Site != "VIIj"| MarkerName != "Loc-34")
msat <- subset(msat, Source != "De Innocentiis et al. 2004 Fish Sci" | Site != "At-A" | MarkerName != "SaGT01" & MarkerName != "SaGT32")
  msat <- subset(msat, Source != "De Innocentiis et al. 2004 Fish Sci" | Site != "Ty-B" | MarkerName != "SaGT31")
msat <- subset(msat, Source != "Diaz-Jaimes & Uribe-Alcocer 2006 Fish Sci" | Site != "Clipperton Islands" | MarkerName != "Tth10")
  msat <- subset(msat, Source != "Diaz-Jaimes & Uribe-Alcocer 2006 Fish Sci" | Site != "Gulf of California" | MarkerName != "Tth5")
  msat <- subset(msat, Source != "Diaz-Jaimes & Uribe-Alcocer 2006 Fish Sci" | Site != "Mexico" | MarkerName != "Tth21" & MarkerName != "Tth5")
  msat <- subset(msat, Source != "Diaz-Jaimes & Uribe-Alcocer 2006 Fish Sci" | Site != "South-west Revillagigedo Islands" | MarkerName != "Tth8")
msat <- subset(msat, Source != "Ding et al. 2009 Mol Ecol Res" | MarkerName != "PLL4")
msat <- subset(msat, Source != "Dos Santos et al. 2008 Mol Ecol Res" | MarkerName != "Elf19" & MarkerName != "Elf37" & MarkerName != "Elf39" & MarkerName != "Elf44" & MarkerName != "Elf46" & MarkerName != "Elf50")
msat <- subset(msat, Source != "Funes et al. 2004 Mol Ecol Notes" | MarkerName != "SseCA47" & MarkerName != "SseGATA3" & MarkerName != "SseGATA9")
msat <- subset(msat, Source != "Galarza et al. 2009 PNAS" | Site != "Blanes" | MarkerName != "Omel27" & MarkerName != "Omel38")
  msat <- subset(msat, Source != "Galarza et al. 2009 PNAS" | Site != "Cabo de Gata" | MarkerName != "Omel54")
  msat <- subset(msat, Source != "Galarza et al. 2009 PNAS" | Site != "Herradura" | MarkerName != "Omel54")
msat <- subset(msat, Source != "Garoia et al. 2006 Mol Ecol Notes" | MarkerName != "Sos(AC)3" & MarkerName != "Sos(AC)45")
msat <- subset(msat, Source != "Gilbert-Horvath et al. 2006 Mol Ecol" | Site != "Big Creek")
msat <- subset(msat, Source != "Gold et al. 2010 Fish Res" | Site != "Cumana" | MarkerName != "Sbr18")
  msat <- subset(msat, Source != "Gold et al. 2010 Fish Res" | Site != "Isla de Margarita" | MarkerName != "Sbr19")
  msat <- subset(msat, Source != "Gold et al. 2010 Fish Res" | Site != "North Coast" | MarkerName != "Sbr16")
  msat <- subset(msat, Source != "Gold et al. 2010 Fish Res" | Site != "South Coast" | MarkerName != "Sbr28")
  msat <- subset(msat, Source != "Gold et al. 2010 Fish Res" | Site != "West Coast" | MarkerName != "Sbr14" & MarkerName != "Sbr28" & MarkerName != "Sbr36")
msat <- subset(msat, Source != "Gomez-Uchida & Banks 2005 CJFAS"  | Site != "C2" & Site != "C3" & Site != "C5" & Site != "O1B" & Site != "W1")
msat <- subset(msat, Source != "Gomez-Uchida 2006 Dissertation" | Site != "10OR01" & Site != "6WA03" & Site != "9OR01")
msat <- subset(msat, Source != "Gonzalez et al. 2008 BMC Evol" | Site != "Guinea" | MarkerName != "TA102" & MarkerName != "TA113" & MarkerName != "TA117")
  msat <- subset(msat, Source != "Gonzalez et al. 2008 BMC Evol" | Site != "Indian" | MarkerName != "TTH208")
  msat <- subset(msat, Source != "Gonzalez et al. 2008 BMC Evol" | Site != "Pacific" | MarkerName != "TTH208")
msat <- subset(msat, Source != "Gonzalez-Wanguemert et al. 2010 JEMBE" | MarkerName != "Dvul6" & MarkerName != "Dvul61")
  msat <- subset(msat, Source != "Gonzalez-Wanguemert et al. 2010 JEMBE" | Site != "Azores" | MarkerName != "DVUL4")
  msat <- subset(msat, Source != "Gonzalez-Wanguemert et al. 2010 JEMBE" | Site != "Banyuls" | MarkerName != "DVUL4")
  msat <- subset(msat, Source != "Gonzalez-Wanguemert et al. 2010 JEMBE" | Site != "Canarian Islands" | MarkerName != "DVUL4")
  msat <- subset(msat, Source != "Gonzalez-Wanguemert et al. 2010 JEMBE" | Site != "Faro" | MarkerName != "DVUL4")
  msat <- subset(msat, Source != "Gonzalez-Wanguemert et al. 2010 JEMBE" | Site != "Mallorca" | MarkerName != "DVUL4")
  msat <- subset(msat, Source != "Gonzalez-Wanguemert et al. 2010 JEMBE" | Site != "Murcia" | MarkerName != "DVUL4")
  msat <- subset(msat, Source != "Gonzalez-Wanguemert et al. 2010 JEMBE" | Site != "Tunisia" | MarkerName != "DVUL4")
msat <- subset(msat, Source != "Guo et al. 2007 Mol Ecol Notes" | MarkerName != "Lru005" & MarkerName != "Lru008" & MarkerName != "Lru013" & MarkerName != "Lru019" & MarkerName != "Lru022" & MarkerName != "Lru026" & MarkerName != "Lru028" & MarkerName != "Lru029" & MarkerName != "Lru031" & MarkerName != "Lru039")
msat <- subset(msat, Source != "Hauser et al. 2007 Book" | MarkerName != "Sma10" & MarkerName != "Sra16-5")
msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | MarkerName != "PL09")
  msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | Site != "Alftanes" | MarkerName != "List1003" & MarkerName != "PL09" & MarkerName != "PL142" & MarkerName != "PL92")
  msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | Site != "Amrum" | MarkerName != "PL09" & MarkerName != "PL92")
  msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | Site != "Balgzand" | MarkerName != "PL09" & MarkerName != "PL115" & MarkerName != "PL92")
  msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | Site != "Bay of Vilaine" | MarkerName != "PL09" & MarkerName != "PL115" & MarkerName != "PL92")
  msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | Site != "Faeroe's Plateau" | MarkerName != "PL09")
  msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | Site != "Femer Baelt" | MarkerName != "List1003" & MarkerName != "PL09" & MarkerName != "PL142" & MarkerName != "PL92")
  msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | Site != "Irish Sea" | MarkerName != "PL09" & MarkerName != "PL142")
  msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | Site != "Lofoten" | MarkerName != "PL09" & MarkerName != "PL142" & MarkerName != "PL92")
  msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | Site != "Oban" | MarkerName != "PL09")
  msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | Site != "Terschellinger Bank" | MarkerName != "PL09" & MarkerName != "PL115" & MarkerName != "PL92")
  msat <- subset(msat, Source != "Hoarau et al. 2002 Mol Ecol" | Site != "Trondheim fjord" | MarkerName != "PL09" & MarkerName != "PL142" & MarkerName != "PL92")
msat <- subset(msat, Source != "Hyde et al. 2008 Mol Ecol" | MarkerName != "Sra15-8" & MarkerName != "Sra7-2")
msat <- subset(msat, Source != "Jeong et al. 2003 Fish Sci" | Site != "Hiroshima (pre-stocking)" | MarkerName != "Acs1" & MarkerName != "Acs3")
  msat <- subset(msat, Source != "Jeong et al. 2003 Fish Sci" | Site != "Yosu" | MarkerName != "Acs1")
msat <- subset(msat, Source != "Jorgensen et al. 2005 ICES J Mar Sci" | Site != "RU03-2" | MarkerName != "Cpa107")
msat <- subset(msat, Source != "Jorgensen et al. 2005 Mol Ecol" | Site != "Baltic Proper (BP02)" | MarkerName != "Cpa101" & MarkerName != "Cpa107" & MarkerName != "Cpa111")
  msat <- subset(msat, Source != "Jorgensen et al. 2005 Mol Ecol" | Site != "Gulf of Riga (R03)" | MarkerName != "Cpa107")
  msat <- subset(msat, Source != "Jorgensen et al. 2005 Mol Ecol" | Site != "Hano Bay (HB02)" | MarkerName != "Cha1017")
msat <- subset(msat, Source != "Jue 2010 Thesis" | Site != "Campeche Bank" | MarkerName != "CA-2" & MarkerName != "GAG045")
  msat <- subset(msat, Source != "Jue 2010 Thesis" | Site != "Northwest Florida Shelf" | MarkerName != "GAG045")
msat <- subset(msat, Source != "Karaiskou et al. 2009 J Fish Bio"| Site != "Mesaloggi" | MarkerName != "SauNINRA")
  msat <- subset(msat, Source != "Karaiskou et al. 2009 J Fish Bio"| Site != "Thermaikos" | MarkerName != "Saul47" & MarkerName != "SauNINRA")
msat <- subset(msat, Source != "Karlsson et al. 2008 Fish Bull" | MarkerName != "Soc637" & MarkerName != "Soc645" & MarkerName != "Soc649" & MarkerName != "Soc650" & MarkerName != "Soc653" & MarkerName != "Soc661" & MarkerName != "Soc666" & MarkerName != "Soc675" & MarkerName != "Soc727")
msat <- subset(msat, Source != "Karlsson et al. 2008 Mol Ecol Res" | MarkerName != "Soc501" & MarkerName != "Soc513" & MarkerName != "Soc525" & MarkerName != "Soc549" & MarkerName != "Soc556" & MarkerName != "Soc567" & MarkerName != "Soc570" & MarkerName != "Soc610" & MarkerName != "Soc625" & MarkerName != "Soc626")
msat <- subset(msat, Source != "Karlsson et al. 2009 Mar Bio" | Site != "Galveston" | MarkerName != "Lca20")
  msat <- subset(msat, Source != "Karlsson et al. 2009 Mar Bio" | Site != "Port Isabel" | MarkerName != "Lca20")
  msat <- subset(msat, Source != "Karlsson et al. 2009 Mar Bio" | Site != "Port Lavaca" | MarkerName != "Lca20")
msat <- subset(msat, Source != "Kim et al. 2003" | MarkerName != "KOP10" & MarkerName != "KOP11" & MarkerName != "KOP14" & MarkerName != "KOP20" & MarkerName != "KOP5" & MarkerName != "KOP6")
msat <- subset(msat, Source != "Kim et al. 2007 Mol Ecol Notes" | MarkerName != "Phz12" & MarkerName != "Phz3" & MarkerName != "Phz8")
msat <- subset(msat, Source != "Kim et al. 2009 Fish Sci" | spp != "Kareius bicoloratus" | MarkerName != "Phz3")
  msat <- subset(msat, Source != "Kim et al. 2009 Fish Sci" | spp != "Limanda punctatissima" | MarkerName != "Phz2" & MarkerName != "Phz3")
  msat <- subset(msat, Source != "Kim et al. 2009 Fish Sci" | spp != "Pseudopleuronectes yokohamae" | MarkerName != "Phz12" & MarkerName != "Phz3")
msat <- subset(msat, Source != "Kim et al. 2009 Mol Ecol Notes" | MarkerName != "KOP28" & MarkerName != "KOP29" & MarkerName != "KOP37" & MarkerName != "KOP47" & MarkerName != "KOP49" & MarkerName != "KOP52" & MarkerName != "KOP54" & MarkerName != "KOP56" & MarkerName != "KOP59" & MarkerName != "KOP64" & MarkerName != "KOP70" & MarkerName != "KOP74" & MarkerName != "KOP90")
msat <- subset(msat, Source != "Kim et al. 2010 Fish Sci" | Site != "Eastern Coast of Korean Peninsula" | MarkerName != "Gmo-C83" & MarkerName != "Gmo19" & MarkerName != "KGM5" & MarkerName != "KGM9")
  msat <- subset(msat, Source != "Kim et al. 2010 Fish Sci" | Site != "Japan" | MarkerName != "Gmo-C82" & MarkerName != "Gmo19" & MarkerName != "Tch18" & MarkerName != "Tch20")
  msat <- subset(msat, Source != "Kim et al. 2010 Fish Sci" | Site != "Southern Coast of Korean Peninsula" | MarkerName != "Gmo-C83" & MarkerName != "Gmo19" & MarkerName != "KGM5" & MarkerName != "Tch18")
  msat <- subset(msat, Source != "Kim et al. 2010 Fish Sci" | Site != "Western Coast of Korean Peninsula" | MarkerName != "Gmo19" & MarkerName != "KGM5" & MarkerName != "Tch18")
msat <- subset(msat, Source != "Kim et al. 2010 J Fish Bio" | Site != "Boryeong (BR)" | MarkerName != "Kop16" & MarkerName != "Kop21")
  msat <- subset(msat, Source != "Kim et al. 2010 J Fish Bio" | Site != "Wando (WD)" | MarkerName != "Kop26")
msat <- subset(msat, Source != "Koedprang et al. 2007 Fish Sci" | spp != "Epinephelus akaara" | MarkerName != "Em-07" & MarkerName != "Em-10")
  msat <- subset(msat, Source != "Koedprang et al. 2007 Fish Sci" | spp != "Epinephelus bleekeri" | MarkerName != "Em-07")
  msat <- subset(msat, Source != "Koedprang et al. 2007 Fish Sci" | spp != "Epinephelus coioides" | MarkerName != "Em-07")
  msat <- subset(msat, Source != "Koedprang et al. 2007 Fish Sci" | spp != "Epinephelus maculatus" | MarkerName != "Em-10")
  msat <- subset(msat, Source != "Koedprang et al. 2007 Fish Sci" | spp != "Epinephelus malabaricus" | MarkerName != "Em-10")
  msat <- subset(msat, Source != "Koedprang et al. 2007 Fish Sci" | spp != "Epinephelus merra" | MarkerName != "CA-07")
  msat <- subset(msat, Source != "Koedprang et al. 2007 Fish Sci" | spp != "Epinephelus ongus" | MarkerName != "Em-07" & MarkerName != "Em-10")
msat <- subset(msat, Source != "Lage & Kornfield 1999 Mol Ecol" | MarkerName != "Mae110")
msat <- subset(msat, Source != "Limborg et al. 2009 MEPS" | Site != "Adriatic Sea" | MarkerName != "Spsp133")
  msat <- subset(msat, Source != "Limborg et al. 2009 MEPS" | Site != "Belt Sea (BEL)" | MarkerName != "Spsp133")
  msat <- subset(msat, Source != "Limborg et al. 2009 MEPS" | Site != "Bornholm Basin (BOR06)" | MarkerName != "Spsp133")
  msat <- subset(msat, Source != "Limborg et al. 2009 MEPS" | Site != "Celtic Sea" | MarkerName != "Spsp170")
  msat <- subset(msat, Source != "Limborg et al. 2009 MEPS" | Site != "Gdansk Deep (GDA)" | MarkerName != "Spsp154")
  msat <- subset(msat, Source != "Limborg et al. 2009 MEPS" | Site != "Northern Kattegat (KAT)" | MarkerName != "Spsp170")
msat <- subset(msat, Source != "Liu et al. 2007 Mol Ecol Notes 7:1178-1180 Acanthopagrus schlegelii" | MarkerName != "Bsb_4" & MarkerName != "Bsb_6")
msat <- subset(msat, Source != "Lundy et al. 1999 Mol Ecol") #all loci out of HWE
msat <- subset(msat, Source != "Lundy et al. 2000 Mol Ecol" | Site != "Bay of Biscay south 1999" | MarkerName != "Hk20" & MarkerName != "Hk34b")
msat <- subset(msat, Source != "Lynch et al. 2010 Fish Bull" | Site != "Chesapeake Bay (Gloucester Point, VA)" | MarkerName != "Asa16")
  msat <- subset(msat, Source != "Lynch et al. 2010 Fish Bull" | Site != "New England (Gloucester, MA)" | MarkerName != "Asa16")
msat <- subset(msat, Source != "Ma & Chen 2009 Cons Gen" | MarkerName != "BaFd11" & MarkerName != "Vem103" & MarkerName != "Vem14")
msat <- subset(msat, Source != "McClelland et al. 2005 J Fish Bio" | Site != "NE Georges Bank" | MarkerName != "Pam27DU")
  msat <- subset(msat, Source != "McClelland et al. 2005 J Fish Bio" | Site != "St. Mary's Bay (SW Nova Scotia)" | MarkerName != "Pam4DU")
msat <- subset(msat, Source != "Menezes et al. 2008 J Fish Bio" | Site != "Japan" | MarkerName != "SK-5")
msat <- subset(msat, Source != "Miao et al. 2009 Cons Gen Verasper moseri" | MarkerName != "Vemo12")
msat <- subset(msat, Source != "Morishima et al. 2009 Mol Ecol Res" | Site != "Kushimoto Wakayama prefecture, Pacific Ocean" | MarkerName != "Kto11" & MarkerName != "Kto42")
  msat <- subset(msat, Source != "Morishima et al. 2009 Mol Ecol Res" | Site != "Shimane prefecture, Sea of Japan" | MarkerName != "Kto09" & MarkerName != "Kto42")
msat <- subset(msat, Source != "O'Connell et al. 1998 J Fish Bio" | Site != "Port Chalmers late, AK" | MarkerName != "Cha123" & MarkerName != "Cha63")
  msat <- subset(msat, Source != "O'Connell et al. 1998 J Fish Bio" | Site != "Rocky Bay late, AK" | MarkerName != "Cha123")
msat <- subset(msat, Source != "Olsen et al. 2002 Mol Ecol Notes" | Site != "Togiak Bay, AK" | MarkerName != "Cpa101" & MarkerName != "Cpa109")
  msat <- subset(msat, Source != "Olsen et al. 2002 Mol Ecol Notes" | Site != "Uganik Bay, Kodiak Island, AK" | MarkerName != "Cpa102" & MarkerName != "Cpa109" & MarkerName != "Cpa112")
msat <- subset(msat, Source != "O'Reilly et al. 2004 Mol Ecol" | Site != "Funka Bay" | MarkerName != "Tch11" & MarkerName != "Tch13" & MarkerName != "Tch14" & MarkerName != "Tch19" & MarkerName != "Tch20" & MarkerName != "Tch6")
  msat <- subset(msat, Source != "O'Reilly et al. 2004 Mol Ecol" | Site != "North Central Bering Sea" | MarkerName != "Tch11" & MarkerName != "Tch14" & MarkerName != "Tch15" & MarkerName != "Tch19" & MarkerName != "Tch6")
  msat <- subset(msat, Source != "O'Reilly et al. 2004 Mol Ecol" | Site != "Prince William Sound, AK B" | MarkerName != "Tch11" & MarkerName != "Tch14" & MarkerName != "Tch15" & MarkerName != "Tch19" & MarkerName != "Tch6")
  msat <- subset(msat, Source != "O'Reilly et al. 2004 Mol Ecol" | Site != "Puget Sound, WA" | MarkerName != "Tch11" & MarkerName != "Tch14" & MarkerName != "Tch19" & MarkerName != "Tch3" & MarkerName != "Tch6")
  msat <- subset(msat, Source != "O'Reilly et al. 2004 Mol Ecol" | Site != "Shelikof Strait, AK A" | MarkerName != "Tch10" & MarkerName != "Tch11" & MarkerName != "Tch14" & MarkerName != "Tch15" & MarkerName != "Tch19" & MarkerName != "Tch6" & MarkerName != "Tch8")
  msat <- subset(msat, Source != "O'Reilly et al. 2004 Mol Ecol" | Site != "Shelikof Strait, AK B" | MarkerName != "Tch11" & MarkerName != "Tch14" & MarkerName != "Tch19" & MarkerName != "Tch6")
  msat <- subset(msat, Source != "O'Reilly et al. 2004 Mol Ecol" | Site != "Unimak, AK A" | MarkerName != "Tch11" & MarkerName != "Tch14" & MarkerName != "Tch15" & MarkerName != "Tch19" & MarkerName != "Tch6" & MarkerName != "Tch8")
  msat <- subset(msat, Source != "O'Reilly et al. 2004 Mol Ecol" | Site != "Unimak, AK B" | MarkerName != "Tch11" & MarkerName != "Tch14" & MarkerName != "Tch19" & MarkerName != "Tch6" & MarkerName != "Tch8")
  msat <- subset(msat, Source != "O'Reilly et al. 2004 Mol Ecol" | Site != "Unimak, AK C" | MarkerName != "Tch11" & MarkerName != "Tch14" & MarkerName != "Tch19" & MarkerName != "Tch3" & MarkerName != "Tch6")
msat <- subset(msat, Source != "Pakaki et al. 2009 Mol Ecol Res" | Site != "Ionian Sea" | MarkerName != "Ee2-376" & MarkerName != "Ee2-507")
  msat <- subset(msat, Source != "Pakaki et al. 2009 Mol Ecol Res" | Site != "South" | MarkerName != "Ee2-477")
msat <- subset(msat, Source != "Pampoulie & Danielsdottir 2008 Fish Res" | Site != "Reykjanes Ridge Iceland (MaG96)" | MarkerName != "SEB33")
  msat <- subset(msat, Source != "Pampoulie & Danielsdottir 2008 Fish Res" | Site != "Southwest Iceland (MaIc97)" | MarkerName != "SEB31")
  msat <- subset(msat, Source != "Pampoulie & Danielsdottir 2008 Fish Res" | Site != "Icelandic Shelf (MeIc97)" | MarkerName != "SEB31")
  msat <- subset(msat, Source != "Pampoulie & Danielsdottir 2008 Fish Res" | Site != "Norway (ViNo01)" | MarkerName != "SEB31")
  msat <- subset(msat, Source != "Pampoulie & Danielsdottir 2008 Fish Res" | Site != "Southwest Iceland (ViIc96)" | MarkerName != "SEB31" & MarkerName != "SEB33" & MarkerName != "SEB45")
msat <- subset(msat, Source != "Pampoulie et al. 2009 ICES J Mar Sci" | Site != 5)
msat <- subset(msat, Source != "Pardo et al. 2005 Mol Ecol Notes" | MarkerName != "Sma-USC1" & MarkerName != "Sma-USC2")
msat <- subset(msat, Source != "Perez-Enriquez & Taniguchi 1999 Fish Sci" | MarkerName != "Pma4")
  msat <- subset(msat, Source != "Perez-Enriquez & Taniguchi 1999 Fish Sci" | Site != "Kochi")
msat <- subset(msat, Source != "Ponce et al. 2006 Mol Ecol Notes" | MarkerName != "PauGAH3")
msat <- subset(msat, Source != "Pumitinsee et al. 2009 Aq Res" | Site != "Trang (July)" | MarkerName != "CA06" & MarkerName != "EM10")
msat <- subset(msat, Source != "Renshaw et al. 2007 Mol Ecol Notes" | spp != "Lutjanus analis" | MarkerName != "Lsy10" & MarkerName != "Lsy14" & MarkerName != "Och11")
  msat <- subset(msat, Source != "Renshaw et al. 2007 Mol Ecol Notes" | spp != "Lutjanus synagris" | MarkerName != "Och6")
  msat <- subset(msat, Source != "Renshaw et al. 2007 Mol Ecol Notes" | spp != "Ocyurus chrysurus" | MarkerName != "Och14")
msat <- subset(msat, Source != "Renshaw et al. 2009 Mol Ecol Res" | MarkerName != "Sbr20" & MarkerName != "Sbr35")
msat <- subset(msat, Source != "Rhodes et al. 2003 Mar Bio" | Site != "Great Barrier Reef" | MarkerName != "Epo1")
  msat <- subset(msat, Source != "Rhodes et al. 2003 Mar Bio" | Site != "Marshall Islands" | MarkerName != "Epo1")
  msat <- subset(msat, Source != "Rhodes et al. 2003 Mar Bio" | Site != "Palau" | MarkerName != "Epo1")
msat <- subset(msat, Source != "Riccioni et al. 2010 PNAS" | Site != "Camogli Shore, Ligurian Sea" | MarkerName != "Tth157" & MarkerName != "Tth208" & MarkerName != "Tth62")
  msat <- subset(msat, Source != "Riccioni et al. 2010 PNAS" | Site != "Carlo Forte Island, Southwestern Sardinia" | MarkerName != "Tth208" & MarkerName != "Tth62")
  msat <- subset(msat, Source != "Riccioni et al. 2010 PNAS" | Site != "Central Adriatic Sea" | MarkerName != "Tth10" & MarkerName != "Tth157" & MarkerName != "Tth208" & MarkerName != "Tth34")
  msat <- subset(msat, Source != "Riccioni et al. 2010 PNAS" | Site != "Puerto Mazarron, Alboran Sea" | MarkerName != "Tth1-31" & MarkerName != "Tth157" & MarkerName != "Tth5" & MarkerName != "Tth62")
  msat <- subset(msat, Source != "Riccioni et al. 2010 PNAS" | Site != "Southern Tyrrhenian Sea" | MarkerName != "Tth208" & MarkerName != "Tth34")
msat <- subset(msat, Source != "Romo et al. 2006 Fish Sci" | Site != "Hakatajima, Ehime" | MarkerName != "Vmo4" & MarkerName != "Vva11")
  msat <- subset(msat, Source != "Romo et al. 2006 Fish Sci" | Site != "Souma, Fukushima" | MarkerName != "Vva9")
msat <- subset(msat, Source != "Roques et al. 1999 Mol Ecol 8:1703" | Site != "Newfoundland")
msat <- subset(msat, Source != "Roques et al. 2001 Mol Ecol" | Site != "Grand Banks 2" & Site != "Labrador U2H")
msat <- subset(msat, Source != "Roques et al. 2002 Mar Bio" | Site != "East Greenland" & Site != "Faroe Islands" & Site != "Irminger Sea" & Site != "West" & Site != "West Greenland")
msat <- subset(msat, Source != "Roques et al. 2006 Mol Ecol Notes" | MarkerName != "Dvul2" & MarkerName != "Dvul63")
msat <- subset(msat, Source != "Saillant & Gold 2006 Fish Bull" | Site != "Port Aransas, TX 1997" | MarkerName != "Lca22")
msat <- subset(msat, Source != "Saillant et al. 2006 ICES J Mar Sci" | Site != "Brownsville bycatch" | MarkerName != "Prs303")
  msat <- subset(msat, Source != "Saillant et al. 2006 ICES J Mar Sci" | Site != "Port Mansfield bycatch" | MarkerName != "Prs229")
msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MAEGr01 (Greenland East)" | MarkerName != "SEB30")
  msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MASWIc00 (SW-Iceland)" | MarkerName != "SEB31" & MarkerName != "SEB33")
  msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MECIrm01 (Central Irminger Sea)" | MarkerName != "SEB37" & MarkerName != "SEB46")
  msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MECIrm03 (Central Irminger Sea)" | MarkerName != "SEB37")
  msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MEEGr00 (Greenland East)" | MarkerName != "SEB46")
  msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MEEGr01 (Greenland East)" | MarkerName != "SEB37" & MarkerName != "SEB46")
  msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MENAF02J01" | MarkerName != "SEB46")
  msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MENAFO1F01" | MarkerName != "SEB46")
  msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MESEIc01 (SE-Iceland)" | MarkerName != "SEB46")
  msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MESIrm01 (Southern Irminger Sea)" | MarkerName != "SEB46")
  msat <- subset(msat, Source != "Schmidt 2005 Dissertation" | Site != "MEWGr01 (Greenland West)" | MarkerName != "SEB37")
msat <- subset(msat, Source != "Sekino & Hara 2001 Mar Biotech" | Site != "Iwate" | MarkerName != "Po42")
msat <- subset(msat, Source != "Sekino et al. 2001 Mar Biotech" | Site != "Akita" | MarkerName != "Sth3B")
  msat <- subset(msat, Source != "Sekino et al. 2001 Mar Biotech" | Site != "Aomori" | MarkerName != "Sth24")
msat <- subset(msat, Source != "Shubina et al. 2009 Mol Biol" | MarkerName != "Tch5" & MarkerName != "Tch14" & MarkerName != "Tch17" & MarkerName != "Tch19")
msat <- subset(msat, Source != "Shulzitski 2005 Thesis" | Site != "Dry Tortugas, FL" | MarkerName != "La39")
  msat <- subset(msat, Source != "Shulzitski 2005 Thesis" | Site != "Gladden Spit" | MarkerName != "La27a" & MarkerName != "LaC-16")
  msat <- subset(msat, Source != "Shulzitski 2005 Thesis" | Site != "Jupiter, FL" | MarkerName != "La39")
  msat <- subset(msat, Source != "Shulzitski 2005 Thesis" | Site != "Mayaguez" | MarkerName != "La39" & MarkerName != "La49a")
  msat <- subset(msat, Source != "Shulzitski 2005 Thesis" | Site != "Roatan" | MarkerName != "La39")
msat <- subset(msat, Source != "Small et al. 2005 TAFS" | Site != "Port Gamble, WA 2002")
msat <- subset(msat, Source != "Takagi et al. 1999 Thunnus orientalis Fish Sci" | spp != "Thunnus alalunga" | Site != "Celebes Sea" | MarkerName != "Ttho-1" & MarkerName != "Ttho-4" & MarkerName != "Ttho-7")
  msat <- subset(msat, Source != "Takagi et al. 1999 Thunnus orientalis Fish Sci" | Site != "E Australia" | MarkerName != "Ttho-1" & MarkerName != "Ttho-4")
msat <- subset(msat, Source != "Takagi et al. 2001 Fish Bull" | Site != "Biscay Bay, NE Atlantic" | MarkerName != "Ttho-7")
msat <- subset(msat, Source != "Tang et al. 2009 Mol Ecol Res" | MarkerName != "Sa2683" & MarkerName != "Sa2873")
msat <- subset(msat, Source != "Tian et al. 2009 Cons Gen" | MarkerName != "Kabi06" & MarkerName != "Kabi22")
msat <- subset(msat, Source != "Tysklind et al. 2009 Mol Ecol Res Limanda limanda" | Site != "Irish Sea" | MarkerName != "DAC2-15" & MarkerName != "DAC4-34" & MarkerName != "DAG2-22" & MarkerName != "DAG4-91")
  msat <- subset(msat, Source != "Tysklind et al. 2009 Mol Ecol Res Limanda limanda" | Site != "North Sea" | MarkerName != "DAC2-15" & MarkerName != "DAC4-34" & MarkerName != "DAC5-70" & MarkerName != "DAG2-22" & MarkerName != "DAG4-91")
msat <- subset(msat, Source != "van Herwerden et al. 2006 Fish Res" | Site != "Aden" | MarkerName != "L42Sc")
  msat <- subset(msat, Source != "van Herwerden et al. 2006 Fish Res" | Site != "Bandar Abbas" | MarkerName != "H96Sc")
  msat <- subset(msat, Source != "van Herwerden et al. 2006 Fish Res" | Site != "Dhofar" | MarkerName != "L42Sc")
  msat <- subset(msat, Source != "van Herwerden et al. 2006 Fish Res" | Site != "Muscat" | MarkerName != "L42Sc")
msat <- subset(msat, Source != "Wang et al. 2010 Gen Mol Res" | MarkerName != "Mimi-12" & MarkerName != "Mimi-3" & MarkerName != "Mimi-6")
msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "10PB(04)" | MarkerName != "MmerUEAW01" & MarkerName != "MpouBW7" & MarkerName != "Tch10")
  msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "11PB(05)" | MarkerName != "MmerUEAW01")
  msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "12PB(03)" | MarkerName != "MpouBW13")
  msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "13PB(03)" | MarkerName != "Tch10")
  msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "14CS(03)" | MarkerName != "MmerUEAW01" & MarkerName != "MpouBW8")
  msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "15CS(04)" | MarkerName != "MmerUEAW01" & MarkerName != "MpouBW8")
  msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "16BB" | MarkerName != "Tch10")
  msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "1PpB" | MarkerName != "MpouBW7" & MarkerName != "Tch10")
  msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "2SB" | MarkerName != "Tch10")
  msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "5RB(05)" | MarkerName != "MmerUEAW01")
  msat <- subset(msat, Source != "Was et al. 2008 ICES J Mar Sci" | Site != "8HS(05)" | MarkerName != "Tch10")
msat <- subset(msat, Source != "Was et al. 2010 Mar Bio" | Site != "Baltic Sea 2006" | MarkerName != "PL06" & MarkerName != "PL09" & MarkerName != "PL92")
  msat <- subset(msat, Source != "Was et al. 2010 Mar Bio" | Site != "Bantry Bay" | MarkerName != "PL06" & MarkerName != "PL09" & MarkerName != "PL115" & MarkerName != "PL167")
  msat <- subset(msat, Source != "Was et al. 2010 Mar Bio" | Site != "Cork Harbour" | MarkerName != "PL06" & MarkerName != "PL09" & MarkerName != "PL115" & MarkerName != "PL142")
  msat <- subset(msat, Source != "Was et al. 2010 Mar Bio" | Site != "Galway Bay" | MarkerName != "PL06" & MarkerName != "PL09" & MarkerName != "PL167")
  msat <- subset(msat, Source != "Was et al. 2010 Mar Bio" | Site != "Iceland 2006" | MarkerName != "PL06" & MarkerName != "PL09" & MarkerName != "PL142" & MarkerName != "PL167")
  msat <- subset(msat, Source != "Was et al. 2010 Mar Bio" | Site != "Irish Sea 2006" | MarkerName != "LIST1003" & MarkerName != "PL06" & MarkerName != "PL09" & MarkerName != "PL167")
msat <- subset(msat, Source != "Watts et al. 1999 Mol Ecol" | MarkerName != "LIST1004")
msat <- subset(msat, Source != "Watts et al. 2004 J Sea Res" | Site != "Balgzand" | MarkerName != "LIST1002")
  msat <- subset(msat, Source != "Watts et al. 2004 J Sea Res" | Site != "Zone 1" | MarkerName != "LIST1002")
msat <- subset(msat, Source != "Watts et al. 2009 ICES J Mar Sci" | Site != "AUC" | MarkerName != "PI92")
  msat <- subset(msat, Source != "Watts et al. 2009 ICES J Mar Sci" | Site != "BLA" | MarkerName != "PI92")
  msat <- subset(msat, Source != "Watts et al. 2009 ICES J Mar Sci" | Site != "FIR" | MarkerName != "PI92")
  msat <- subset(msat, Source != "Watts et al. 2009 ICES J Mar Sci" | Site != "TOR" | MarkerName != "LIST1-001" & MarkerName != "LIST1-003")
msat <- subset(msat, Source != "Xing et al. 2009 Cons Gen" | MarkerName != "Spma20" & MarkerName != "Spma21")
msat <- subset(msat, Source != "Yagashita & Kobayashi 2008 Mol Ecol Res" | Site !="Tanabe Bay, Wakayama" | MarkerName != "Scja01" & MarkerName != "Scja02" & MarkerName != "Scja03" & MarkerName != "Scja07" & MarkerName != "Scja09")
  msat <- subset(msat, Source != "Yagashita & Kobayashi 2008 Mol Ecol Res" | Site != "western Wakasa Bay, Kyoto" | MarkerName != "Scja02" & MarkerName != "Scja04" & MarkerName != "Scja07")
msat <- subset(msat, Source != "Yu et al. 2002 Mar Biotech" | Site != "I-Lan 1999" | MarkerName != "EJ2" & MarkerName != "EJ27.1" & MarkerName != "EJ27.2" & MarkerName != "EJ35" & MarkerName != "EJ9")
  msat <- subset(msat, Source != "Yu et al. 2002 Mar Biotech" | Site != "Peng-Hu 2000")
msat <- subset(msat, Source != "Zarraonaindia et al. 2009 ICES J Mar Sci")
msat <- subset(msat, Source != "Zatcoff et al. 2004 Mar Bio" | spp != "Epinephelus morio" | Site != "Eastern Gulf of Mexico" | MarkerName != "Gag23")
  msat <- subset(msat, Source != "Zatcoff et al. 2004 Mar Bio" | spp != "Epinephelus morio" | Site != "Florida Keys" | MarkerName != "Gag23" & MarkerName != "Mbo48")
  msat <- subset(msat, Source != "Zatcoff et al. 2004 Mar Bio" | spp != "Epinephelus morio" | Site != "Mexico" | MarkerName != "Gag23")
  msat <- subset(msat, Source != "Zatcoff et al. 2004 Mar Bio" | spp != "Epinephelus morio" | Site != "North and South Carolina" | MarkerName != "Gag23" & MarkerName != "Mbo48")
msat <- subset(msat, Source != "Zhao et al. 2009 Cons Gen: Epinephelus awoara" | MarkerName != "Epaw17" & MarkerName != "Epaw19" & MarkerName != "Epaw25" & MarkerName != "Epaw34" & MarkerName != "Epaw6")
msat <- subset(msat, Source != "Zhao et al. 2009 Cons Gen: Epinephelus septemfasciatus" | MarkerName != "Ese33" & MarkerName != "Ese36" & MarkerName != "Ese43")
dim(msat) #19900 rows

######## Check CrossSpp ########

#make sure all numerical
sort(unique(msat$CrossSpp)) #looks like all numbers

#check fractional cases
inds <- which(!(msat$CrossSpp %in% c(0, 1)) & !is.na(msat$CrossSpp)) #grab instances where not 0 or 1 (and not NA)
length(inds) #144 cases
t(sort(unique(msat$CrossSpp[inds]))) #look to be fractional
msat[inds, c('CrossSpp', 'NumMarkers')] #all CrossSpp not 0 or 1 have NumMarkers >1

#check places where CrossSpp isn't known
inds <- is.na(msat$CrossSpp)
sum(inds) #15 --> fixed below
msat[inds, c('spp', 'Source', 'CrossSpp', 'MarkerName', 'file')]

#add CrossSpp info based on papers (in Durand etal. 2005 Mar. Bio. 147:313-322 some just unknown)
msat$CrossSpp[msat$Source == "Appleyard, Williams & Ward. 2004. CCAMLR Science 11:21-32." & msat$MarkerName == "To2"] <- 0 #based on ref paper
  msat$CrossSpp[msat$Source == "Appleyard, Williams & Ward. 2004. CCAMLR Science 11:21-32." & msat$MarkerName == "To5"] <- 0 #based on ref paper
msat$CrossSpp[msat$Source == "Diaz-Jaimes & Uribe-Alcocer 2006 Fish Sci"] <- 1 #ALL markers developed in cross species based on ref paper
msat$CrossSpp[msat$Source == "McPherson et al. 2001 Journal of Fish Biology 59:356-370" & msat$MarkerName == "Cpa113"] <- 1 #based on ref paper
  msat$CrossSpp[msat$Source == "McPherson et al. 2001 Journal of Fish Biology 59:356-370" & msat$MarkerName == "Cpa102"] <- 1 #based on ref paper
  msat$CrossSpp[msat$Source == "McPherson et al. 2001 Journal of Fish Biology 59:356-370" & msat$MarkerName == "Cpa108"] <- 1 #based on ref paper
msat$CrossSpp[msat$Source == "Gomez-Uchida 2006 Dissertation"] <- 0.6 #based on ref papers
msat$CrossSpp[msat$Source == "Hyde et al. 2008 Mol Ecol" & msat$MarkerName == "Spi4"] <- 0 #based on ref paper
msat$CrossSpp[msat$Source == "Mitchell 2006 Thesis" & msat$MarkerName == "Cha113"] <- 1 #based on ref paper
  msat$CrossSpp[msat$Source == "Mitchell 2006 Thesis" & msat$MarkerName == "Cha123"] <- 1 #based on ref paper
  msat$CrossSpp[msat$Source == "Mitchell 2006 Thesis" & msat$MarkerName == "Cha17"] <- 1 #based on ref paper
  msat$CrossSpp[msat$Source == "Mitchell 2006 Thesis" & msat$MarkerName == "Cha20"] <- 1 #based on ref paper
  msat$CrossSpp[msat$Source == "Mitchell 2006 Thesis" & msat$MarkerName == "Cha63"] <- 1 #based on ref paper
msat$CrossSpp[msat$Source == "Pampoulie & Danielsdottir 2008 Fish Res" & msat$MarkerName == "Smen10"] <- 1 #based on ref paper
  msat$CrossSpp[msat$Source == "Pampoulie & Danielsdottir 2008 Fish Res" & msat$MarkerName == "Smen5"] <- 1 #based on ref paper
  msat$CrossSpp[msat$Source == "Pampoulie & Danielsdottir 2008 Fish Res" & msat$spp == "Sebastes mentella" & msat$MarkerName == "Smen10"] <- 0 #based on ref paper
  msat$CrossSpp[msat$Source == "Pampoulie & Danielsdottir 2008 Fish Res" & msat$spp == "Sebastes mentella" & msat$MarkerName == "Smen5"] <- 0 #based on ref paper
msat$CrossSpp[msat$Source == "Sala-Bozano et al. 2009 Mol Ecol"] <- 0.44 #based on ref papers

######## Check NumMarkers ########

#make sure all numerical
sort(unique(msat$NumMarkers)) #looks like all numbers < 0

##########################################################################################################################################

######## Trim & write out data ########
	
#check to make sure all data included
inds <- is.na(msat$He) | is.na(msat$spp) | is.na(msat$lat) | is.na(msat$lon) | is.na(msat$PrimerNote) | is.na(msat$n) #columns need to have data for, already checked most explicitly
sum(inds) #40: all from Wilson et al. (2019), don't have n for specific sites, leave for now
msat[inds, ]
	
#trim to just relevant columns
msat <- msat[, c('spp', 'CommonName', 'Source', 'PrimerNote', 'Country', 'Site', 'lat', 'lon', 'stockid', 'CollectionYear', 
                 'NumMarkers', 'MarkerName', 'CrossSpp', 'n', 'Repeat', 'He', 'Hese', 'file')]
dim(msat) #19900 x 18

#write out msat data (allow multiple loci per line)
write.csv(msat, file = "output/msat_assembled.csv")

#########################################################################################################################################

######## Make duplicates of datapoints that are averages across multiple loci (so each row is one locus) ########

#check to make sure that know NumMarkers for all data points
dim(msat) #19900
inds <- !is.na(msat$NumMarkers) #only want instances where NumMarkers is known
sum(inds) #19900
sum(!inds) #0

#replicate rows so that each is one locus
msatloci <- msat[inds, ][rep(row.names(msat[inds, ]), msat$NumMarkers[inds]), ]
dim(msatloci) #28551 rows

#keep original He
msatloci$He_orig <- msatloci$He
summary(msatloci$He_orig) #0.73
summary(msat$He) #0.74

#re-number so no duplicate row names
rownames(msatloci) <- 1:nrow(msatloci)
nrow(msatloci) #28551 rows

######## Add SE to He where Hese is known ########

#subset to instances where Hese is known
inds <- which(!is.na(msatloci$Hese) & msatloci$NumMarkers > 1) #only when NumMarkers > 1
length(inds) #1723

#add SD row
msatloci$SD <- NA

#generate random deviations to add to the mean
for(i in inds) { #based on study-specific SE
  msatloci$SD[i] <- msatloci$Hese[i]*sqrt(msat$NumMarkers[i])
  e <- rnorm(n = 1, mean = 0, sd = msatloci$SD[i])
  msatloci$He[i] <- msatloci$He_orig[i] + e
}

#check He distribution
summary(msatloci$He_orig[inds]) #mean is 0.70
summary(msatloci$He[inds]) #-0.29 to 2.01, mean = 0.70
msatloci$He[msatloci$He > 1] = 1 #enforce He <= 1
msatloci$He[msatloci$He < 0] = 0 #enforce He >= 0 -->  then exclude monomorphic loci
summary(msatloci$He[inds]) #mean is 0.68
msatloci <- subset(msatloci, He != 0)
summary(msatloci$He[inds]) #mean is 0.69

######## Add SE to He where Hese is not known ########

#grab Hese's currently present in database
msat[!is.na(msat$Hese), c('spp', 'Source', 'Hese')]

#calculate SD (based on known SE values)
hesdmean <- mean(msat$Hese*sqrt(msat$NumMarkers), na.rm = TRUE)
hesdmean #SD = 0.24

#grab values to add SD to (only instances where NumMarkers > 1)
inds <- which(msatloci$NumMarkers > 1 & is.na(msatloci$Hese) & !is.na(msatloci$He))
length(inds) #7879

#generate random deviations to add to the mean
e <- rnorm(length(inds), mean = 0, sd = hesdmean)
msatloci$He[inds] <- msatloci$He_orig[inds] + e #add sd to the mean

#check He distribution
summary(msatloci$He_orig[inds]) #mean is 0.72
summary(msatloci$He[inds]) #-0.42 to 1.95, mean = 0.72
msatloci$He[msatloci$He > 1] = 1 #enforce He <= 1
msatloci$He[msatloci$He < 0] = 0 #enforce He >= 0 -->  then exclude monomorphic loci
summary(msatloci$He[inds]) #mean is 0.70
msatloci <- subset(msatloci, He != 0)
summary(msatloci$He[inds]) #mean is 0.71

######## Format and write out ########

#write out msat data (one locus per line)
write.csv(msatloci, file = "output/msatloci_assembled.csv")

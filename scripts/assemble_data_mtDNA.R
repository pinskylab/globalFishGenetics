setwd('/Users/mpinsky/Documents/Princeton/Latitudinal gradients of diversity/Analysis/')


######################
# read in mtDNA files and assemble
# if needed
######################
#mtDNA1.1 = read.csv('../Data/mtDNA Data/csvs/Fishery lat mtDNA 1.1 2015-10-22 MW.csv', stringsAsFactors=FALSE)
#	head(mtDNA1.1)
#mtDNA1.2 = read.csv('../Data/mtDNA Data/csvs/Fishery lat mtDNA 1.2 2015-10-01 MW.csv', stringsAsFactors=FALSE)
#	head(mtDNA1.2)
#mtDNA1.3 = read.csv('../Data/mtDNA Data/csvs/Fishery lat mtDNA 1.3 2015-10-26 MW.csv', stringsAsFactors=FALSE)
#	head(mtDNA1.3)
#mtDNA1.4 = read.csv('../Data/mtDNA Data/csvs/Fishery lat mtDNA 1.4 2015-10-02 MW.csv', stringsAsFactors=FALSE)
#	head(mtDNA1.4)
#mtDNA1.5 = read.csv('../Data/mtDNA Data/csvs/Fishery lat mtDNA 1.5 2015-10-03 MW.csv', stringsAsFactors=FALSE)
#	head(mtDNA1.5)
#mtDNA2.1 = read.csv('../Data/mtDNA Data/csvs/Fishery lat mtDNA 2.1 2015-09-20 MLP.csv', stringsAsFactors=FALSE)
#	head(mtDNA2.1)
#mtDNA2.2 = read.csv('../Data/mtDNA Data/csvs/Fishery lat mtDNA 2.2 2015-09-20 MLP.csv', stringsAsFactors=FALSE)
#	head(mtDNA2.2)
#mtDNA3 = read.csv('../Data/mtDNA Data/csvs/Fishery lat mtDNA 3 2015-10-26 MW.csv', stringsAsFactors=FALSE)
#	head(mtDNA3)
#mtDNA4 = read.csv('../Data/mtDNA Data/csvs/Fishery lat mtDNA 4 2015-10-26 MW.csv', stringsAsFactors=FALSE)
#	head(mtDNA4)
#mtDNA5 = read.csv('../Data/mtDNA Data/csvs/Fishery lat mtDNA 5 2015-10-28.2 MLP.csv', stringsAsFactors=FALSE)
#	head(mtDNA5)
#
## Trim off rows with no data
#	nrow(mtDNA1.1); inds = is.na(mtDNA1.1$spp) | mtDNA1.1$spp == ''; sum(inds)
#	mtDNA1.1 <- mtDNA1.1[!inds,]; nrow(mtDNA1.1)
#
#	nrow(mtDNA1.2); inds = is.na(mtDNA1.2$spp) | mtDNA1.2$spp == ''; sum(inds)
#	mtDNA1.2 <- mtDNA1.2[!inds,]; nrow(mtDNA1.2)
#
#	nrow(mtDNA1.3); inds = is.na(mtDNA1.3$spp) | mtDNA1.3$spp == ''; sum(inds)
#	mtDNA1.3 <- mtDNA1.3[!inds,]; nrow(mtDNA1.3)
#
#	nrow(mtDNA1.4); inds = is.na(mtDNA1.4$spp) | mtDNA1.4$spp == ''; sum(inds)
#	mtDNA1.4 <- mtDNA1.4[!inds,]; nrow(mtDNA1.4)
#
#	nrow(mtDNA1.5); inds = is.na(mtDNA1.5$spp) | mtDNA1.5$spp == ''; sum(inds)
#	mtDNA1.5 <- mtDNA1.5[!inds,]; nrow(mtDNA1.5)
#
#	nrow(mtDNA2.1); inds = is.na(mtDNA2.1$spp) | mtDNA2.1$spp == ''; sum(inds)
#	mtDNA2.1 <- mtDNA2.1[!inds,]; nrow(mtDNA2.1)
#
#	nrow(mtDNA2.2); inds = is.na(mtDNA2.2$spp) | mtDNA2.2$spp == ''; sum(inds)
#	mtDNA2.2 <- mtDNA2.2[!inds,]; nrow(mtDNA2.2)
#
#	nrow(mtDNA3); inds = is.na(mtDNA3$spp) | mtDNA3$spp == ''; sum(inds)
#	mtDNA3 <- mtDNA3[!inds,]; nrow(mtDNA3)
#
#	nrow(mtDNA4); inds = is.na(mtDNA4$spp) | mtDNA4$spp == ''; sum(inds)
#	mtDNA4 <- mtDNA4[!inds,]; nrow(mtDNA4)
#
#	nrow(mtDNA5); inds = is.na(mtDNA5$spp) | mtDNA5$spp == ''; sum(inds)
#	mtDNA5 <- mtDNA5[!inds,]; nrow(mtDNA5)
#
## Check for duplicate studies
#	# TO DO
#
## verify column names match
#	all(names(mtDNA1.1) == names(mtDNA1.2))
#	all(names(mtDNA1.1) == names(mtDNA1.3))
#	all(names(mtDNA1.1) == names(mtDNA1.4))
#	all(names(mtDNA1.1) == names(mtDNA1.5))
#	all(names(mtDNA1.1) == names(mtDNA2.1))
#	all(names(mtDNA1.1) == names(mtDNA2.2))
#	all(names(mtDNA1.1) == names(mtDNA3))
#	all(names(mtDNA1.1) == names(mtDNA4))
#	all(names(mtDNA1.1) == names(mtDNA4))
#
## merge files
#	mtDNA = rbind(mtDNA1.1, mtDNA1.2, mtDNA1.3, mtDNA1.4, mtDNA1.5, mtDNA2.1, mtDNA2.2, mtDNA3, mtDNA4, mtDNA5) # merge files
#		dim(mtDNA) # 1363

#######################################
# Alternative: read in assembled file
#######################################
mtDNA <- read.csv('../Data/mtDNA Data/Compiled/Fishery lat mtDNA Complete Database 2016-03-02 MLP.csv')


# Calc lat and lon in decimal degrees
	mtDNA$lat_min[is.na(mtDNA$lat_min)] = 0
	mtDNA$lat_sec[is.na(mtDNA$lat_sec)] = 0

	mtDNA$lon_min[is.na(mtDNA$lon_min)] = 0
	mtDNA$lon_sec[is.na(mtDNA$lon_sec)] = 0

	mtDNA$lat = rowSums(cbind(mtDNA$lat_deg, mtDNA$lat_min/60, mtDNA$lat_sec/3600))
		range(mtDNA$lat)	
	mtDNA$lon = rowSums(cbind(mtDNA$lon_deg, mtDNA$lon_min/60, mtDNA$lon_sec/3600))
		range(mtDNA$lon)	

# Fix species names
	mtDNA$spp <- gsub(' $', '', mtDNA$spp) # remove trailing space
	mtDNA$spp[mtDNA$spp == 'Hippocampus Kuda'] <- 'Hippocampus kuda' # fix capitalization
	mtDNA <- mtDNA[mtDNA$spp != 'Merluccius albidus and Merluccius bilinearis',] # remove double spp
	mtDNA$spp[mtDNA$spp ==  "Paralichthys olivaceus `"] <- "Paralichthys olivaceus"
	mtDNA$spp[mtDNA$spp ==  "Sparus aurata L."] <- "Sparus aurata"
	mtDNA$spp[mtDNA$spp ==  "Diplodus sargus sargus"] <- "Diplodus sargus"
	mtDNA$spp[mtDNA$spp ==  "Clupea pallasii pallasii"] <- "Clupea pallasii"
	mtDNA$spp[mtDNA$spp ==  "Centrophorus zeehani"] <- "Centrophorus zeehaani"

# Trim just to relevant columns
mtDNA <- mtDNA[,c('spp', 'CommonName', 'Source', 'Country', 'Site', 'CollectionYear', 'MarkerName', 'n', 'bp', 'He', 'Hese', 'Pi', 'Pise', 'lat', 'lon')]

# Process collection year (to do later)
	# need to clean this first



###################
# QA/QC
###################
# were numeric fields read as numerics?
summary(mtDNA)

# check species names
t(t(sort(unique(mtDNA$spp)))) # to print as one long list. 162 species.

# check locus names
t(t(sort(unique(mtDNA$MarkerName)))) # to print as one long list

# check lat
mtDNA[(mtDNA$lat>90 | mtDNA$lat< -90) & !is.na(mtDNA$lat),c('spp', 'Source', 'Country', 'Site', 'lat')]	

# check lon
mtDNA[(mtDNA$lon>180 | mtDNA$lon< -180) & !is.na(mtDNA$lon),c('spp', 'Source', 'Country', 'Site', 'lon')]	

# check He
mtDNA[(mtDNA$He<0 | mtDNA$He>1) & !is.na(mtDNA$He),c('spp', 'Source', 'Country', 'Site', 'He')]	

# check Pi
mtDNA[(mtDNA$Pi<0 | mtDNA$Pi>1) & !is.na(mtDNA$Pi),c('spp', 'Source', 'Country', 'Site', 'He')]	

# Check species names and make a translation table to FishBase
	require(rfishbase)
	fbdatmtdna <- data.frame(spp=sort(unique(mtDNA$spp)), fbsci=NA) # translation table from my scientific names to those in fishbase

	options(nwarnings=300)
	for(i in 1:nrow(fbdatmtdna)){ # check sci names
		cat(i)
		fbdatmtdna$fbsci[i] <- validate_names(fbdatmtdna$spp[i]) # check that sci names are correct. returned warnings, all about potential synonyms.
	}
		warnings()
		# for(i in 1:nrow(fbdatmtdna)) print(synonyms(fbdatmtdna$fbsci[i])) # there is only ever 1 accepted name

		inds <- which(fbdatmtdna$spp != fbdatmtdna$fbsci)
		fbdatmtdna[inds,] # examine rows where fishbase has a different name
	
	# write out
	write.csv(fbdatmtdna, file='output/fbdat_mtdna.csv', row.names=FALSE)



# Remove problematic datapoints
	# not used right now



# write out mtDNA data (allow multiple loci per line)
	write.csv(mtDNA, file=paste('output/mtDNA_', Sys.Date(), '.csv', sep=''))



#####################################################################################
# write out lat/lon and species data for appending environmental and trait data
#####################################################################################

# Write out list of stocks, if useful
#write.csv(mtDNA[,c('spp', 'Source', 'Country', 'Site')], file=paste('output/stocks_', Sys.Date(), '.csv', sep=''))


	latlonmtdna <- mtDNA[!duplicated(mtDNA[,c('lat', 'lon')]), c('lat', 'lon')]
		latlonmtdna <- latlon[order(latlonmtdna$lat, latlonmtdna$lon),]
	sppsmtdna <- mtDNA[!duplicated(mtDNA$spp), c('spp', 'CommonName')]
	
	# Add FishBase sci name
	sppsmtdna2 <- merge(sppsmtdna, fbdatmtdna)
		dim(sppsmtdna)
		dim(sppsmtdna2)
	
	# Add FB SpecCode
	sppsmtdna2$SpecCode <- NA
	for(i in 1:nrow(sppsmtdna2)){
		cat(i)
		sppsmtdna2$SpecCode[i] <- as.numeric(species(sppsmtdna2$fbsci[i], fields='SpecCode')$SpecCode)
	}
	summary(sppsmtdna2) # no NAs
	
	dim(latlonmtdna) # 1156 locations
	dim(sppsmtdna2) # 162 species

	write.csv(latlonmtdna, file=paste('output/latlon_mtDNA_', Sys.Date(), '.csv', sep=''), row.names=FALSE)
	write.csv(sppsmtdna2, file=paste('output/spps_mtDNA_', Sys.Date(), '.csv', sep=''), row.names=FALSE)



#############################
# plot sampling locations
#############################

require(maps)
require(mapdata)

pdf('figures/map_mtDNA.pdf')

plot(latlon$lon, latlon$lat, col='red', cex=0.5, xlab='Longitude', ylab='Latitude', main='mtDNA data', type='n')
map(database='worldHires', add=TRUE, fill=TRUE, col='grey', boundary=FALSE, interior=FALSE)
points(latlon$lon, latlon$lat, col='red', cex=0.5)

dev.off()




#################################
## Add species and stock data
#################################

# Add overfished status
	tsview <- read.csv('Data/srdb/timeseries_values_views.csv')
	stocks <- read.csv('Data/srdb/stock.csv')

	# find whether overfished or not (for each srdb stock)
	tsmin <- aggregate(list(minB=tsview$B.Bmsytouse), by=list(stockid=tsview$stockid, stocklong=tsview$stocklong), FUN=min) # find minimum B/Bmsy
	tsmin$srdb_overfished <- tsmin$minB<1 # overfished if B/Bmsy ever less than 1

	# add species name to srdb data
	tsminstock <- merge(tsmin, stocks)
		dim(tsmin)
		dim(tsminstock)
	head(tsminstock)

	# summarize by spp
	# would prefer to match by stock, but would have to do that partially by hand
	tsminsp <- aggregate(list(propOF = tsminstock$srdb_overfished), by=list(spp=tsminstock$scientificname), FUN=mean)
		head(tsminsp)
		dim(tsminsp)

	# examine match species from mtDNA to tsminsp
	nodataspp <- sort(setdiff(mtDNA$spp, tsminsp$spp)) # species with no stock status information
	nodataspp
	dataspp <- sort(intersect(mtDNA$spp, tsminsp$spp)) # species with no stock status information
	dataspp
	
		# search for slight mismatches in names
		grep('Clupea', tsminsp$spp, value=TRUE)
		grep('radiata', tsminsp$spp, value=TRUE)
		grep('Ammodytes', tsminsp$spp, value=TRUE)
		grep('thazard', tsminsp$spp, value=TRUE)
		grep('Beryx', tsminsp$spp, value=TRUE)
		grep('argus', tsminsp$spp, value=TRUE)
		grep('herz', tsminsp$spp, value=TRUE)
		grep('labrax', tsminsp$spp, value=TRUE)
		grep('sargus', tsminsp$spp, value=TRUE)
		grep('hynnus', tsminsp$spp, value=TRUE)
		grep('Lutjanus', tsminsp$spp, value=TRUE)
		grep('Merluccius', tsminsp$spp, value=TRUE)
		grep('californicus', tsminsp$spp, value=TRUE)
		grep('Scomber', tsminsp$spp, value=TRUE)
		grep('Trachurus', tsminsp$spp, value=TRUE)
		grep('variegat', tsminsp$spp, value=TRUE)
		
	# fix some names to match mtDNA
	#not needed
	
	# set NA to 0
	tsminsp$propOF[is.na(tsminsp$propOF)] <- 0
	
	# merge in overfished status (at species level)
	mtDNA2 <- merge(mtDNA, tsminsp, all.x=TRUE)
		dim(mtDNA)
		dim(mtDNA2)

	# set missing data to 0
	mtDNA2$propOF[is.na(mtDNA2$propOF)] <- 0

# Add average global fisheries catch
	fao=read.csv("Data/FAO FIGIS/figis_guestnull.csv") # original fao data
	ccols = grep("^X", names(fao)) # columns that have catch data (e.g. "X1950")
	faospp <- aggregate(list(fao[,ccols]), by=list(spp=fao$Scientific.name), FUN=sum) # sum across countries
	ccols2 = grep("^X", names(faospp)) # columns that have catch data (e.g. "sumcatch.X1950")
	faospp$sumcatch = apply(faospp[,ccols2], 1, sum, na.rm=T)
	faospp$nyears = apply(faospp[,ccols2]>1, 1, sum, na.rm=T)
	faospp$avecatch = faospp$sumcatch/faospp$nyears

	# examine match species from mtDNA to tsminsp
	sort(setdiff(mtDNA2$spp, faospp$spp)) # species with no stock status information
	sort(intersect(mtDNA2$spp, faospp$spp)) # species with no stock status information

		# search for slight mismatches in names
		grep('radiata', faospp$spp, value=TRUE)
		grep('alistes', faospp$spp, value=TRUE)
		grep('herz', faospp$spp, value=TRUE)
		grep('pineph', faospp$spp, value=TRUE)
		grep('Glyptocephalus', faospp$spp, value=TRUE)
		grep('Lutjanus', faospp$spp, value=TRUE)
		grep('Merluccius', faospp$spp, value=TRUE)
		grep('major', faospp$spp, value=TRUE)
		grep('maculatus', faospp$spp, value=TRUE)
		grep('Trachurus', faospp$spp, value=TRUE)
		grep('variegat', faospp$spp, value=TRUE)

	# fix names to match my names
	faospp$spp <- as.character(faospp$spp)
	faospp$spp[faospp$spp=='Raja radiata'] <- 'Amblyraja radiata'
	faospp$spp[faospp$spp=='Pseudopleuronectes herzenst.'] <- 'Pleuronectes herzensteini'

	# merge in ave global catch
	mtDNA3 <- merge(mtDNA2, faospp[,c('spp', 'avecatch')], all.x=TRUE)
		dim(mtDNA2)
		dim(mtDNA3)

	# set missing data to 0
	mtDNA3$avecatch[is.na(mtDNA3$avecatch)] <- 0

# Add body size
	fbdatmtdna$maxlength<-NA
	for(i in 1:nrow(fbdatmtdna)){ # get length data
		fbdatmtdna$maxlength[i] <- as.numeric(species(fbdatmtdna$fbsci[i], fields='Length')$Length)
		if(is.na(fbdatmtdna$maxlength[i])){ # if male maxlength is NA, use female maxlength
			fbdatmtdna$maxlength[i] <- as.numeric(species(fbdatmtdna$fbsci[i], fields='LengthFemale')$LengthFemale)
		}
	}
	
	summary(fbdatmtdna)

	mtDNA4 <- merge(mtDNA3, fbdatmtdna[,c('spp', 'maxlength')], all.x=TRUE)
		dim(mtDNA3)
		dim(mtDNA4)

# write out mtDNA data with species traits
	write.csv(mtDNA4, file=paste('output/mtDNAlh_', Sys.Date(), '.csv', sep=''))

# Script to assemble the mtDNA dataset
# Starts from the mtDNA file compiled by hand

#######################################
# Read in assembled file
#######################################
#mtdna <- read.csv('data/mtDNA/Fishery lat mtDNA Complete Database.csv', stringsAsFactors=FALSE)
#srdbmatch <- read.csv("data/srdb_matching/mtdna_to_match.csv", stringsAsFactors=FALSE) # to match genetic data to SRDB stocks
mtdna <- read.csv('data/mtDNA/mtDNA_merged.csv', stringsAsFactors=FALSE)

# Calc lat and lon in decimal degrees
	mtdna$lat_min[is.na(mtdna$lat_min)] = 0
	mtdna$lat_sec[is.na(mtdna$lat_sec)] = 0

	mtdna$lon_min[is.na(mtdna$lon_min)] = 0
	mtdna$lon_sec[is.na(mtdna$lon_sec)] = 0

	mtdna$lat = rowSums(cbind(mtdna$lat_deg, mtdna$lat_min/60, mtdna$lat_sec/3600))
		range(mtdna$lat)	
	mtdna$lon = rowSums(cbind(mtdna$lon_deg, mtdna$lon_min/60, mtdna$lon_sec/3600))
		range(mtdna$lon)	

# remove lines without lat/lon
	inds <- is.na(mtdna$lat) | is.na(mtdna$lon)
	sum(inds) # 0


# merge in srdb stock information	
	# bring in fbsci, stockid, lat, lon (latter two for error-checking)
	names(srdbmatch)[names(srdbmatch)=='lat'] <- 'lat_srdb'
	names(srdbmatch)[names(srdbmatch)=='lon'] <- 'lon_srdb'

	# trim to unique spp, Source, Site
	inds <- duplicated(srdbmatch[,c('spp', 'Source', 'Country', 'Site', 'CollectionYear')])
	inds2 <- duplicated(srdbmatch[,c('spp', 'Source', 'Country', 'Site', 'CollectionYear', 'lat_srdb', 'lon_srdb')])
	sum(inds) # 9
	sum(inds2) # 9: matches: good, means that lat/lon the same for each row we'll keep
		# srdbmatch[inds &  !inds2,]
	nrow(srdbmatch) # 329

	srdbmatch <- srdbmatch[!inds,] # trim to unique sites
	nrow(srdbmatch) # 320

	# sort(setdiff(srdbmatch$Source, msat$Source)) # papers in srdbmatch that aren't in msat. checked by eye the mismatch wasn't a typo.

		nrow(mtdna) # 1332
	mtdna <- merge(mtdna, srdbmatch[,c('spp', 'Source', 'Country', 'Site', 'CollectionYear', 'fbsci', 'stockid', 'lat_srdb', 'lon_srdb')], all.x=TRUE, by=c('spp', 'Source', 'Country', 'Site', 'CollectionYear'))
		nrow(mtdna) # 1332


# Fix species names
	# sort(unique(mtdna$spp))

	mtdna$spp <- gsub(' $', '', mtdna$spp) # remove trailing space
	mtdna$spp[mtdna$spp == 'Hippocampus Kuda'] <- 'Hippocampus kuda' # fix capitalization
	mtdna <- mtdna[mtdna$spp != 'Merluccius albidus and Merluccius bilinearis',] # remove double spp entry
	mtdna$spp[mtdna$spp ==  "Paralichthys olivaceus `"] <- "Paralichthys olivaceus"
	mtdna$spp[mtdna$spp ==  "Sparus aurata L."] <- "Sparus aurata"
	mtdna$spp[mtdna$spp ==  "Diplodus sargus sargus"] <- "Diplodus sargus"
	mtdna$spp[mtdna$spp ==  "Clupea pallasii pallasii"] <- "Clupea pallasii"
	mtdna$spp[mtdna$spp ==  "Centrophorus zeehani"] <- "Centrophorus zeehaani"
	mtdna$spp[mtdna$spp ==  "Amodytes personatus"] <- "Ammodytes personatus"
	mtdna$spp[mtdna$spp ==  "Acanthopagrus schlegeli"] <- "Acanthopagrus schlegelii"
	mtdna$spp[mtdna$spp ==  "Amphiprion akallopsisos"] <- "Amphiprion akallopisos"
	mtdna$spp[mtdna$spp ==  "Saurida elongate"] <- "Saurida elongata"
	mtdna$spp[mtdna$spp ==  "Holocentrus ascensionis"] <- "Holocentrus adscensionis"
	mtdna$spp[mtdna$spp ==  "Pleuragramma antarcticum"] <- "Pleuragramma antarctica"

# Fix locus names
	mtdna$MarkerName <- gsub(' $', '', mtdna$MarkerName) # remove trailing space
	mtdna$MarkerName <- gsub('^ ', '', mtdna$MarkerName) # remove leading space
	mtdna$MarkerName[mtdna$MarkerName %in% c('cyt b', 'Cytb', 'cytochrome b', 'Cytochrome b', 'Cyt b')] <- 'cytb' # standardize names
	mtdna$MarkerName[mtdna$MarkerName %in% c('cytochrome c oxidase subunit 1', 'cytochrome oxidase c subunit I', 'cytochrome-c oxidase I', 'COI')] <- 'cytochrome c oxidase subunit I' # standardize names
	mtdna$MarkerName[mtdna$MarkerName %in% c('D-Loop', 'D-loop region', 'D-loop sequence')] <- 'D-loop' # standardize names
	mtdna$MarkerName[mtdna$MarkerName %in% c('Control region', 'control region (D-loop)', 'CR')] <- 'control region' # standardize names
	mtdna$MarkerName[mtdna$MarkerName %in% c('HVR-1', 'HVR-I', 'hyper-variable region I')] <- 'HVR I' # standardize names
	mtdna$MarkerName[mtdna$MarkerName %in% c('ATPase 6 and 8', 'ATPase6/8')] <- 'ATPase 6, 8' # standardize names
	mtdna$MarkerName[mtdna$MarkerName %in% c('NADH-dehydrogenase subunit 4 (ND-4)', 'ND4 region')] <- 'ND4' # standardize names
	mtdna$MarkerName[mtdna$MarkerName %in% c('Cox1')] <- 'cox1' # standardize names
	

# Process collection year (to do later)
	# need to clean this first



###################
# QA/QC
###################
# were numeric fields read as numerics?
	summary(mtdna)

# check species names
	spp_list <- t(t(sort(unique(mtdna$spp)))) # to print as one long list. 245 species.

# check locus names
	locus_list <- t(t(sort(unique(mtdna$MarkerName)))) # to print as one long list

# check source names
	source_list <- t(t(sort(unique(mtdna$Source)))) # to print as one long list


# check lat
	mtdna[(mtdna$lat>90 | mtdna$lat< -90) & !is.na(mtdna$lat),c('spp', 'Source', 'Country', 'Site', 'lat')]	

	hist(mtdna$lat) # mostly northern hemisphere, but not only


# check lon
	mtdna[(mtdna$lon>180 | mtdna$lon< -180) & !is.na(mtdna$lon),c('spp', 'Source', 'Country', 'Site', 'lon')]	

	hist(mtdna$lon) # # most around 100degW, but quite evenly spread


# compare lat and lon to srdbmatch lat/lon, fill in where missing
	plot(mtdna$lat, mtdna$lat_srdb); abline(0,1) # most fall right on 1:1 line
		mtdna[abs(mtdna$lat - mtdna$lat_srdb)>0.5 & !is.na(mtdna$lat_srdb), c('spp', 'Source', 'Site', 'lat', 'lat_srdb')] # I'm OK with Dammannagoda et al. 2011 and Bremer et al. 1999, since I entered lat/lon in the mtdna sheet
	plot(mtdna$lon, mtdna$lon_srdb); abline(0,1) # most fall right on 1:1 line
		mtdna[abs(mtdna$lon - mtdna$lon_srdb)>0.5 & !is.na(mtdna$lon_srdb), c('spp', 'Source', 'Site', 'lon', 'lon_srdb')] # I'm OK with Tripp-Valdez, Stepien, Scoles et al., Dammannagoda et al., and Bremer et al, since I entered lat/lon on mtdna sheet
	
	mtdna[is.na(mtdna$lat) & !is.na(mtdna$lat_srdb),] # 0
	mtdna[is.na(mtdna$lon) & !is.na(mtdna$lon_srdb),] # 0



# check He
	mtdna[(mtdna$He<0 | mtdna$He>1) & !is.na(mtdna$He),c('spp', 'Source', 'Country', 'Site', 'He')]	# 0

	inds <- is.na(mtdna$He)
	sum(inds) # 48

# check Pi
	mtdna[(mtdna$Pi<0 | mtdna$Pi>1) & !is.na(mtdna$Pi),c('spp', 'Source', 'Country', 'Site', 'He')]	 # 0

	inds <- is.na(mtdna$Pi)
	sum(inds) # 42
	
	sum(is.na(mtdna$Pi) & is.na(mtdna$He)) # 0: good

# All data included?
	inds <- (is.na(mtdna$He) & is.na(mtdna$Pi)) | is.na(mtdna$spp) | is.na(mtdna$lat) | is.na(mtdna$lon) | is.na(mtdna$n)
	sum(inds) # 0
	mtdna[inds,]


# Remove problematic datapoints
	# not used right now


# Trim just to relevant columns

mtdna <- mtdna[,c('spp', 'CommonName', 'Source', 'Country', 'Site', 'lat', 'lon', 'stockid', 'CollectionYear', 'MarkerName', 'n', 'bp', 'He', 'Hese', 'Pi', 'Pise')]
	dim(mtdna) # 1332


# write out mtdna data (allow multiple loci per line)
	#write.csv(mtdna, file='output/mtdna.csv')
	write.csv(mtdna, file='output/mtdna_merged_cleaned.csv')



####################################################################
# Check species names and make a translation table to FishBase
####################################################################
require(rfishbase)
fbdatmtdna <- data.frame(spp=sort(unique(mtdna$spp)), fbsci=NA) # translation table from my scientific names to those in fishbase
fbdatmtdna <- subset(fbdatmtdna, spp!="Branchiostoma belcheri" & spp!="Branchiostoma japonicum") #not in fishbase??

options(nwarnings=300)
nrow(fbdatmtdna)
for(i in 1:nrow(fbdatmtdna)){ # check sci names
	cat(i)
	fbdatmtdna$fbsci[i] <- validate_names(fbdatmtdna$spp[i]) # check that sci names are correct. returned warnings, all about potential synonyms.
}
	warnings() # all related to X "can also be misapplied to other species"
	# for(i in 1:nrow(fbdatmtdna)) print(synonyms(fbdatmtdna$fbsci[i])) # there is only ever 1 accepted name

	inds <- which(fbdatmtdna$spp != fbdatmtdna$fbsci)
	fbdatmtdna[inds,] # examine rows where fishbase has a different name. all look good

# write out
write.csv(fbdatmtdna, file='output/fbdat_mtdna_merged.csv', row.names=FALSE)


#####################################################################################
# write out lat/lon and species data for appending environmental and trait data
#####################################################################################

# Write out list of stocks, if useful
#write.csv(mtdna[,c('spp', 'Source', 'Country', 'Site')], file=paste('output/stocks_', Sys.Date(), '.csv', sep=''))


	latlonmtdna <- mtdna[!duplicated(mtdna[,c('lat', 'lon')]), c('lat', 'lon')]
		latlonmtdna <- latlon[order(latlonmtdna$lat, latlonmtdna$lon),]
	sppsmtdna <- mtdna[!duplicated(mtdna$spp), c('spp', 'CommonName')]
	
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

	write.csv(latlonmtdna, file=paste('output/latlon_mtdna_', Sys.Date(), '.csv', sep=''), row.names=FALSE)
	write.csv(sppsmtdna2, file=paste('output/spps_mtdna_', Sys.Date(), '.csv', sep=''), row.names=FALSE)



#############################
# plot sampling locations
#############################

require(maps)
require(mapdata)

pdf('figures/map_mtdna_merged.pdf')

plot(mtdna$lon, mtdna$lat, col='red', cex=0.5, xlab='Longitude', ylab='Latitude', main='mtdna data', type='n')
map(database='world', add=TRUE, fill=TRUE, col='grey', boundary=FALSE, interior=FALSE)
points(mtdna$lon, mtdna$lat, col='red', cex=0.5)

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
	nodataspp <- sort(setdiff(mtdna$spp, tsminsp$spp)) # species with no stock status information
	nodataspp
	dataspp <- sort(intersect(mtdna$spp, tsminsp$spp)) # species with no stock status information
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
		
	# fix some names to match mtdna
	#not needed
	
	# set NA to 0
	tsminsp$propOF[is.na(tsminsp$propOF)] <- 0
	
	# merge in overfished status (at species level)
	mtdna2 <- merge(mtdna, tsminsp, all.x=TRUE)
		dim(mtdna)
		dim(mtdna2)

	# set missing data to 0
	mtdna2$propOF[is.na(mtdna2$propOF)] <- 0

# Add average global fisheries catch
	fao=read.csv("Data/FAO FIGIS/figis_guestnull.csv") # original fao data
	ccols = grep("^X", names(fao)) # columns that have catch data (e.g. "X1950")
	faospp <- aggregate(list(fao[,ccols]), by=list(spp=fao$Scientific.name), FUN=sum) # sum across countries
	ccols2 = grep("^X", names(faospp)) # columns that have catch data (e.g. "sumcatch.X1950")
	faospp$sumcatch = apply(faospp[,ccols2], 1, sum, na.rm=T)
	faospp$nyears = apply(faospp[,ccols2]>1, 1, sum, na.rm=T)
	faospp$avecatch = faospp$sumcatch/faospp$nyears

	# examine match species from mtdna to tsminsp
	sort(setdiff(mtdna2$spp, faospp$spp)) # species with no stock status information
	sort(intersect(mtdna2$spp, faospp$spp)) # species with no stock status information

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
	mtdna3 <- merge(mtdna2, faospp[,c('spp', 'avecatch')], all.x=TRUE)
		dim(mtdna2)
		dim(mtdna3)

	# set missing data to 0
	mtdna3$avecatch[is.na(mtdna3$avecatch)] <- 0

# Add body size
	fbdatmtdna$maxlength<-NA
	for(i in 1:nrow(fbdatmtdna)){ # get length data
		fbdatmtdna$maxlength[i] <- as.numeric(species(fbdatmtdna$fbsci[i], fields='Length')$Length)
		if(is.na(fbdatmtdna$maxlength[i])){ # if male maxlength is NA, use female maxlength
			fbdatmtdna$maxlength[i] <- as.numeric(species(fbdatmtdna$fbsci[i], fields='LengthFemale')$LengthFemale)
		}
	}
	
	summary(fbdatmtdna)

	mtdna4 <- merge(mtdna3, fbdatmtdna[,c('spp', 'maxlength')], all.x=TRUE)
		dim(mtdna3)
		dim(mtdna4)

# write out mtdna data with species traits
	write.csv(mtdna4, file=paste('output/mtdnalh_', Sys.Date(), '.csv', sep=''))

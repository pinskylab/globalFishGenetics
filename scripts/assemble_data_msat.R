# assemble the msat data and do basic QA/QC

# basic functions
calcHe = function(x){ # calculate expected heterozygosit from allele frequencies. does not check that they sum to 1
	x = x[!is.na(x)]
	h = 1-sum(x^2)
	return(h)
}

sumna = function(x){ # remove NAs, unless all are in NAs, then return NA
	if(all(is.na(x))) return(NA)
	else return(sum(x, na.rm=TRUE))
}


# main script
setwd('/Users/mpinsky/Documents/Princeton/Latitudinal gradients of diversity/Analysis/')


###########
## msats ##
###########
# read in msat files
msat1 = read.csv('../Data/msat Data/csv/Fishery lat msats 000 2015-05-20 SSP.csv', stringsAsFactors=FALSE)
#	head(msat1)
msat2 = read.csv('../Data/msat Data/csv/Fishery lat msats 001 2015-07-04 MLP.csv', stringsAsFactors=FALSE)
#	head(msat2)
msat3 = read.csv('../Data/msat Data/csv/Fishery lat msats 002 2015-08-08 SSP.csv', stringsAsFactors=FALSE)
#	head(msat3)
	msat3 <- msat3[msat3$spp != '', ] # trim empty rows
msat4 = read.csv('../Data/msat Data/csv/Fishery lat msats 100 2015-08-20 MLP.csv', stringsAsFactors=FALSE)
#	head(msat4)
	msat4 <- msat4[,!(names(msat4)=='X')] # drop an extra column
	msat4 <- msat4[msat4$spp != '', ] # trim empty rows
msat5 = read.csv('../Data/msat Data/csv/Fishery lat msats 101 2015-07-17 SSP.csv', stringsAsFactors=FALSE)
#	head(msat5)
	msat5 = msat5[,!(names(msat5)=='X')] # drop an extra column
	msat5 <- msat5[msat5$spp != '', ] # trim empty rows
msat6 = read.csv('../Data/msat Data/csv/Fishery lat msats 200 2015-02-10 MLP.csv', stringsAsFactors=FALSE)
#	head(msat6)
	msat6 = msat6[,!(names(msat6)=='X')] # drop an extra column
	msat6 <- msat6[msat6$spp != '', ] # trim empty rows
msat7 = read.csv('../Data/msat Data/csv/Fishery lat msats 201 2015-10-13 MLP.csv', stringsAsFactors=FALSE)
#	head(msat7)

ppdat = read.csv("../Data/msat Data/Pinsky&Palumbi/ppdat_2016-03-04wLL.csv", stringsAsFactors=FALSE) # Pinsky & Palumbi 2014 data
	ppdat$NumMarkers <- ppdat$NumMarker # for some reason, the contents weren't copied over before
#	head(ppdat)

# add an indicator for source file
msat1$file <- 'msats000'
msat2$file <- 'msats001'
msat3$file <- 'msats002'
msat4$file <- 'msats100'
msat5$file <- 'msats101'
msat6$file <- 'msats200'
msat7$file <- 'msats201'
ppdat$file <- 'ppdat'

# Check for duplicate studies between msat and ppdat
#overlaps = intersect(c(msat1$spp, msat2$spp, msat3$spp, msat4$spp, msat5$spp), ppdat$spp) # species included in both datasets
#sort(unique(ppdat$Source[ppdat$spp %in% overlaps])) # list of studies from ppdat on these species
#sort(unique(c(msat1$Source[msat1$spp %in% overlaps], msat2$Source[msat2$spp %in% overlaps], msat3$Source[msat3$spp %in% overlaps], msat4$Source[msat4$spp %in% overlaps], msat5$Source[msat5$spp %in% overlaps]))) # studies from msat dataframe on these species

dups = c('Antoro et al. 2006 Marine Biotechnology 8:17-26', 'Bekkevold et al. 2005 Evolution 59:2656-2668', 'Canino et al. 2005 Mol Ecol Notes 5:908-910', 'Carlsson et al. 2004. Mol. Ecol. 13, 3345-3356.', 'Roques et al. 2007 Mol. Ecol. Notes 7:661-663', 'Selkoe et al. 2006 Ecology 87:3082-3094', 'Wang et al. 2010 Gen & Mol Res 9:931-934', 'Yagishita & Kobayashi 2008 Mol. Ecol. Res. 8:302-304') # studies in msat data.frames that appear to overlap those in ppdat
msat1 = msat1[!(msat1$Source %in% dups),]
msat2 = msat2[!(msat2$Source %in% dups),]
msat3 = msat3[!(msat3$Source %in% dups),]
msat4 = msat4[!(msat4$Source %in% dups),]
msat5 = msat5[!(msat5$Source %in% dups),]
msat6 = msat6[!(msat6$Source %in% dups),]
msat7 = msat7[!(msat7$Source %in% dups),]

	# use approximate matching (commented out since decisions already made)
#studs <- vector('list', length=7)
#studs[[1]] <- sort(unique(msat1$Source))
#studs[[2]] <- sort(unique(msat2$Source))
#studs[[3]] <- sort(unique(msat3$Source))
#studs[[4]] <- sort(unique(msat4$Source))
#studs[[5]] <- sort(unique(msat5$Source))
#studs[[6]] <- sort(unique(msat6$Source))
#studs[[7]] <- sort(unique(msat7$Source))
#studs[[8]] <- sort(unique(ppdat$Source))
#
#poss_matches <- vector('list', length=length(studs))
#for(j in 1:(length(studs)-1)){
#	poss_matches[[j]] <- vector('list', length=length(studs)-j)
#	for(i in (j+1):length(studs)){
#		poss_matches[[j]][[i-j]] <- lapply(studs[[j]], FUN=agrep, x=studs[[i]], fixed=TRUE, value=TRUE, max.distance=list(all=0.2), useBytes=TRUE)
#		lens <- sapply(poss_matches[[j]][[i-j]], FUN=length)
#		if(any(lens>0)){
#			print(paste('j=', j, ', i=', i, sep=''))
#			for(k in 1:sum(lens>0)){
##				print(cbind(studs[[j]][which(lens>0)[k]], paste(poss_matches[[j]][[i-j]][[which(lens>0)[k]]], collapse=',')))
#			}
#		}
#	}
#}

	# examine and trim out 4 matching studies between msat6 and msat7
		# Bagley & Geller 1998
#	msat6[msat6$Source=='Bagley & Geller 1998 Molecular Ecology 7:1083-1090',]
#	msat7[msat7$Source=='Bagley&Geller 1998 Molecular Ecology 7:1083-1090',] 
	msat7 <- msat7[msat7$Source !='Bagley&Geller 1998 Molecular Ecology 7:1083-1090',] # keep msat6 version, since it appropriately dropped Ra11 for signs of a null allele

		# Bagley et al. 1999
#	msat6[msat6$Source=='Bagley et al. 1999 Marine Biology 134:609-620',!grepl('^p', names(msat6))]
#	msat7[msat7$Source=='Bagley et al. 1999 Marine Biology 134:609-620',!grepl('^p', names(msat7))] 
	msat6 <- msat6[msat6$Source !='Bagley et al. 1999 Marine Biology 134:609-620',] # keep msat7 version, since it appropriately dropped Ra1, Ra5, and Ra6 for signs of HWP departure
	
		# Garcia de Leon et al. 1997 Molecular Ecology 6:51-62
#	msat6[msat6$Source=='Garcia de Leon et al. 1997 Molecular Ecology 6:51-62',!grepl('^p', names(msat6))]
#	msat7[msat7$Source=='Garcia de Leon et al. 1997 Molecular Ecology 6:51-62',!grepl('^p', names(msat7))] 
	msat7 <- msat7[msat7$Source !='Garcia de Leon et al. 1997 Molecular Ecology 6:51-62',] # keep msat6 version arbitrarily, since they look the same

		# Naciri et al. 1999 Journal of Heredity 90:592-596		
#	msat6[msat6$Source=='Naciri et al. 1999 Journal of Heredity 90:592-596',!grepl('^p', names(msat6))]
#	msat7[msat7$Source=='Naciri et al. 1999 The Journal of Heredity 90: 591-596',!grepl('^p', names(msat7))] 
	msat7 <- msat7[msat7$Source !='Naciri et al. 1999 The Journal of Heredity 90: 591-596',] # keep msat6 version since I was more careful about which sites to keep or drop

		# Bahri-Sfar et al. 2000 Proc. R. Soc. Lond. B 267:929-935 (msats201) and Bahri-Sfar et al. 2000 Proceedings of the Royal Society B 267:929-935 (msats200): found by accident, not by fuzzy matching above
	msat6[msat6$Source=='Bahri-Sfar et al. 2000 Proceedings of the Royal Society B 267:929-935',!grepl('^p', names(msat6))]
	msat7[msat7$Source=='Bahri-Sfar et al. 2000 Proc. R. Soc. Lond. B 267:929-935',!grepl('^p', names(msat7))] 
	msat7 <- msat7[msat7$Source !='Bahri-Sfar et al. 2000 Proc. R. Soc. Lond. B 267:929-935',] # keep msat6 version since I was more careful about which sites to keep or drop

# verify column names match
	all(names(msat1) == names(msat2))
	all(names(msat1) == names(msat3))
	all(names(msat1) == names(msat4))
	all(names(msat1) == names(msat5))
	all(names(msat1) == names(msat6))
	all(names(msat1) == names(msat7)) # fails
		nn <- setdiff(names(msat1), names(msat7))
		for(i in 1:length(nn)) msat7[[nn[i]]] <- NA # pad extra columns in msat7
		setdiff(names(msat1), names(msat7))
		setdiff(names(msat7), names(msat1))

# merge files from students
	msat = rbind(msat1, msat2, msat3, msat4, msat5, msat6, msat7) # merge student files
		dim(msat) # 4585
			
# merge in file from Pinsky & Palumbi
	msat$lat = NA # need this to merge with ppdat
	nms = setdiff(names(msat), names(ppdat)) # column names in msat that aren't in ppdat
	nms
	for(i in 1:length(nms)) ppdat[[nms[i]]] = NA # add these columns to ppdat
	ppdat$lon <- ppdat$long # move to new column name to match msat
	ppdatmsat = ppdat[,names(msat)] # trim Pinsky & Palumbi to only msats and only columns that match the student data
	dim(ppdatmsat) # 11888
	dim(msat) # 4585
	
	msat = rbind(msat, ppdatmsat) # merge in Pinsky & Palumbi 2014 data
		dim(msat) # 16473 x 90
		
# Process allele frequencies into He where needed
	inds = is.na(msat$He) & !is.na(msat$p1)
	sum(inds) # 54
#	msat[inds,] # all LeClair et al. 2006 TAFS
	msat$He[inds] = apply(msat[inds,grep('p[[:digit:]]+', names(msat))], MARGIN=1, FUN=calcHe)
	 
# Calc lat and lon in decimal degrees
	inds = is.na(msat$lat) & !is.na(msat$lat_deg)
	sum(inds) # 15067
#	msat[inds,c('lat', 'lat_deg', 'lat_min')]
	msat$lat[inds] = rowSums(cbind(msat$lat_deg[inds], msat$lat_min[inds]/60, msat$lat_sec[inds]/3600), na.rm=TRUE)

	inds = !is.na(msat$lon_deg)
	sum(inds) # 15067
	msat$lon[inds] = rowSums(cbind(msat$lon_deg[inds], msat$lon_min[inds]/60, msat$lon_sec[inds]/3600), na.rm=TRUE)


# fix species name mistakes
	# msat[msat$spp %in% c("Dicentrarchus labrax", "Dicentrarchus labrax "),!grepl('^p', names(msat6))]
	msat$spp[msat$spp == "Dicentrarchus labrax "] <- "Dicentrarchus labrax" # fix the extra space

	#msat[msat$spp %in% c("Pleuronectes herzensteini", "Pseudopleuronectes herzensteini"),!grepl('^p', names(msat6))]
	msat$spp[msat$spp == "Pleuronectes herzensteini"] <- "Pseudopleuronectes herzensteini" # Pseudopleuronectes is correct according to FishBase
	
	# msat[msat$spp %in% c("Rhomboplites aurorubens", "Rhomboplites aurorubens "),!grepl('^p', names(msat6))]
	msat$spp[msat$spp == "Rhomboplites aurorubens "] <- "Rhomboplites aurorubens" # fix the extra space

	msat$spp[msat$spp == "Aphanius fascitus"] <- "Aphanius fasciatus" # fix the missing a

	msat$spp[msat$spp == "Gnatholepsis cauerensis"] <- "Gnatholepis cauerensis" # fix the missing a

	msat$spp[msat$spp == "Lophius piscatorious"] <- "Lophius piscatorius" # fix the missing a

	# msat[msat$spp %in% c("Pagrus auratus"),!grepl('^p', names(msat6))] # all samples from New Zealand, therefore Pagrus auratus, not Sparus aurata

	# msat[msat$spp %in% c("Sebastes marinus"),!grepl('^p', names(msat6))] # FB calls this norvegicus

	# msat[msat$spp %in% c("Chaetodon lunulatus"),!grepl('^p', names(msat6))] # no reason to suspect this is actually C. lunula (from the paper)
	
	# msat[msat$spp %in% c("Trachinotus falcatus"),c('spp', 'CommonName', 'Source', 'PrimerNote', 'Country', 'Site')] # no reason to suspect this is actually C. lunula (from the paper)


###########
# QA/QC
###########

# were numeric fields read as numerics?
summary(msat)
	# repeats are NOT all numeric

# check species names visually
#t(t(sort(unique(msat$spp)))) # to print as one long list. 363 species.

# check for slightly misspelled species names with fuzzy matching
#spps <- sort(unique(msat$spp))
#
#poss_matches <- vector('list', length=length(spps))
#for(i in 1:(length(spps)-1)){
#	poss_matches[[i]] <- agrep(pattern=spps[i], x=spps[(i+1):length(spps)], fixed=TRUE, value=TRUE, max.distance=list(all=0.2), useBytes=TRUE)
#	len <- length(poss_matches[[i]])
#	if(len>0){
#		print(paste(', i=', i, sep=''))
#		print(cbind(spps[i], paste(poss_matches[[i]][which(len>0)], collapse=',')))
#	}
#}



# find Fishbase species names
	fbdat <- data.frame(spp=sort(unique(msat$spp)), fbsci=NA) # translation table from my scientific names to those in fishbase

	options(nwarnings=300)
	for(i in 1:nrow(fbdat)){ # check sci names
		cat(i)
		if(!(fbdat$spp[i] %in% c('Sebastes marinus', 'Pagrus auratus', 'Chaetodon lunulatus', 'Trachinotus falcatus'))){ # chokes for no apparent reason on S. marinus. Returns 3 answers for P. auratus and C. lunulatus.
			fbdat$fbsci[i] <- validate_names(fbdat$spp[i]) # check that sci names are correct. returned warnings, all about potential synonyms.
		}
		if(fbdat$spp[i] == 'Sebastes marinus'){
			fbdat$fbsci[i] <- 'Sebastes norvegicus' # what FB calls it (from a web search)
		}
		if(fbdat$spp[i] == 'Pagrus auratus'){
			fbdat$fbsci[i] <- 'Pagrus auratus' # correct since our samples are from New Zealand
		}
		if(fbdat$spp[i] == 'Chaetodon lunulatus'){
			fbdat$fbsci[i] <- 'Chaetodon lunulatus' # assume this isn't C. lunula (other option FB gives)
		}
		if(fbdat$spp[i] == 'Trachinotus falcatus'){
			fbdat$fbsci[i] <- 'Trachinotus falcatus' # correct since our samples are from Puerto Rico (not Indo-pacific)
		}
	}
		warnings()
		# for(i in 1:nrow(fbdat)) print(synonyms(fbdat$fbsci[i])) # there is only ever 1 accepted name

		inds <- which(fbdat$spp != fbdat$fbsci)
		fbdat[inds,] # examine rows where fishbase has a different name
		
		
	# write out
	write.csv(fbdat, file='output/fbdat_msat.csv', row.names=FALSE)

# check repeat definitions: STILL need to do
	x=as.numeric(msat$Repeat)
	msat$Repeat[is.na(x)]

# check that allele freqs sum to 1
	x = apply(msat[,grep('p[[:digit:]]+', names(msat))], MARGIN=1, FUN=sumna) 
	which(x<0.9)
	msat[which(x<0.99 & is.na(msat$He)),] # none

# check lat
	msat[msat$lat>90 & !is.na(msat$lat),c('spp', 'Source', 'Country', 'Site', 'lat_deg', 'lat_min', 'lat_sec', 'lat')]	
	msat[msat$lat < -90 & !is.na(msat$lat),c('spp', 'Source', 'Country', 'Site', 'lat_deg', 'lat_min', 'lat_sec', 'lat')]	

	sum(is.na(msat$lat) & !is.na(msat$lat_deg)) # missing lat but not lat_deg?

	hist(msat$lat)

# check lon
	msat[msat$lon > 180 & !is.na(msat$lon),c('spp', 'Source', 'Country', 'Site', 'lon_deg', 'lon_min', 'lon_sec', 'lon')]	
	msat[msat$lon < -180 & !is.na(msat$lon),c('spp', 'Source', 'Country', 'Site', 'lon_deg', 'lon_min', 'lon_sec', 'lon')]	
	msat[msat$lon == 0 & !is.na(msat$lon),c('spp', 'Source', 'Country', 'Site', 'lon_deg', 'lon_min', 'lon_sec', 'lon')]	

	sum(is.na(msat$lon) & !is.na(msat$lon_deg)) # missing lon but not lon_deg?

	hist(msat$lon) # 

# check He
	msat[(msat$He<0 | msat$He>1) & !is.na(msat$He),c('spp', 'Source', 'Country', 'Site', 'He')]	

# check CrossSpp
	sort(unique(msat$CrossSpp))	
	inds = which(!(msat$CrossSpp %in% c(0, 1))) # good that all are numbers
	length(inds)
	t(sort(unique(msat$CrossSpp[inds])))

# check places where NumMarkers isn't known (not yet done)

# Process collection year (to do later)
	# need to clean this first
	
# Trim just to relevant columns
msat <- msat[,c('spp', 'CommonName', 'Source', 'PrimerNote', 'Country', 'Site', 'CollectionYear', 'NumMarkers', 'MarkerName', 'CrossSpp', 'n', 'Repeat', 'He', 'Hese', 'file', 'lat', 'lon')]


# write out msat data (allow multiple loci per line)
	write.csv(msat, file=paste('output/msat_', Sys.Date(), '.csv', sep=''))


#####################################################################################
# write out lat/lon and species data for appending environmental and trait data
#####################################################################################
	latlon <- msat[!duplicated(msat[,c('lat', 'lon')]), c('lat', 'lon')]
		latlon <- latlon[order(latlon$lat, latlon$lon),]
	spps <- msat[!duplicated(msat$spp), c('spp', 'CommonName')]
	
	# Add FishBase sci name
	spps2 <- merge(spps, fbdat)
		dim(spps)
		dim(spps2)
	
	# Add FB SpecCode
	spps2$SpecCode <- NA
	for(i in 1:nrow(spps2)){
		cat(i)
		spps2$SpecCode[i] <- as.numeric(species(spps2$fbsci[i], fields='SpecCode')$SpecCode)
	}

	summary(spps2) # no NAs

	dim(latlon) # 2124 locations
	dim(spps2) # 360 species

	write.csv(latlon, file=paste('output/latlon_msat_', Sys.Date(), '.csv', sep=''), row.names=FALSE)
	write.csv(spps2, file=paste('output/spps_msat_', Sys.Date(), '.csv', sep=''), row.names=FALSE)



#############################
# plot sampling locations
#############################

require(maps)
require(mapdata)
pdf('figures/map_msats.pdf')

plot(latlon$lon, latlon$lat, col='red', cex=0.5, xlab='Longitude', ylab='Latitude', main='Microsatellite data')
map(database='worldHires', add=TRUE, fill=TRUE, col='grey', boundary=FALSE, interior=FALSE)
points(latlon$lon, latlon$lat, col='red', cex=0.5)

dev.off()



#####################################################################################################
# Make duplicates of datapoints that are averages across multiple loci (so each row is one locus)
#####################################################################################################

dim(msat) # 16473
inds = !is.na(msat$NumMarkers) # only where NumMarkers is known
sum(inds) # 16406
sum(!inds) # only 67 missing
msatloci = msat[inds,][rep(row.names(msat[inds,]), msat$NumMarkers[inds]),] # replicate rows so that each is one locus
dim(msatloci) # 22239
msatloci$He_orig = msatloci$He # save original He

# Add error to He where Hese is known
inds = which(!is.na(msatloci$Hese) & msatloci$NumMarkers>1)
	length(inds) # 706
sds <- msatloci$Hese[inds]*sqrt(msatloci$NumMarkers[inds]) # convert SE to SD
	summary(sds)
e = rnorm(length(inds), mean = 0, sd = sds) # generate random deviations to add to the mean
msatloci$He[inds] = msatloci$He_orig[inds] + e # add to the mean
	summary(msatloci$He[inds])
	msatloci$He[msatloci$He>1] = 1 # enforce He <= 1
	msatloci$He[msatloci$He<0] = 0 # enforce He >=0
	summary(msatloci$He[inds])

# Add error to He where Hese is not known
	msat[!is.na(msat$Hese), c('spp', 'Source', 'Hese')] # many Hese values available
hesdmean = mean(msat$Hese*sqrt(msat$NumMarkers), na.rm=T) # convert SE to SD. SD=0.1897
	hesdmean
	# summary(lm(I(Hese*sqrt(NumMarkers)) ~ He, data = msat)) # p=0.06 lack definitive evidence that SD increases with mean
inds = which(msatloci$NumMarkers > 1 & is.na(msatloci$Hese) & !is.na(msatloci$He)) # values to add error to
	length(inds) # 4513
e = rnorm(length(inds), mean = 0, sd = hesdmean)
	summary(e)
msatloci$He[inds] = msatloci$He_orig[inds] + e
	summary(msatloci$He[inds])
	msatloci$He[msatloci$He>1] = 1 # enforce He <= 1
	msatloci$He[msatloci$He<0] = 0 # enforce He >=0
	summary(msatloci$He[inds])

# new row names so that none are duplicated
rownames(msatloci) = 1:nrow(msatloci)

# write out msat data (one locus per line)
	write.csv(msatloci, file=paste('output/msatloci_', Sys.Date(), '.csv', sep=''))



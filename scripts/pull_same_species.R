remove(list = ls())

library(data.table)
library(tidyverse)

# put data from same studies of mtdna & msat into one csv file

msat <- read.csv("output/msatloci_assembled.csv", stringsAsFactors = FALSE)
mtdna <- read.csv("output/mtdna_assembled.csv", stringsAsFactors = FALSE)

 #check for duplicate studies between oldmsat studies & ppdat
overlaps <- intersect(msat$spp, mtdna$spp) #species included in both mtdna & msat

msat_studies <- sort(unique(msat$Source[msat$spp %in% overlaps])) #msat studies with those species
mtdna_studies <- sort(unique(mtdna$Source[mtdna$spp %in% overlaps]))# mtdna studies with those species

all_studies <- sort(c(msat_studies, mtdna_studies))

#matching studies
matches <- c("Ball et al. 2007 Mar Biol 150:1321-1332", "Ball et al. 2007 Mar. Biol. 150:1321-1332", 
             "Cardenas et al. 2009 Fish Res", "Cardenas et al. 2009 Fisheries Research 100:109-115", 
             "Carson et al. (2011) Fishery Bulletin 109:416-428", "Carson et al. 2011 Fishery Bulletin 109:416-428",
             "Castillo-Olguin et al. 2012 Ciencias Marinas 38(4):635-652", "Castillo-Olguín,Uribe-Alcocer & Díaz-Jaimes (2012) Ciencias Marinas 38:635-652", 
             "Catarino et al. 2017 PLoS ONE 12:e0174988", "Catrino et al. (2017) PLoS ONE 12:e0174988",
             "Coscia et al. (2012) Heredity 108:537-546", "Coscia et al. 2012 Heredity 108:537-546",
             "Cuveliers et al. (2012) Mar. Biol. 159:1239-1253", "Cuveliers et al. 2012 Marine Biology 159:1239-1253",
             "Dammannagoda et al. 2008 Fish Res", "Dammannagoda et. al 2008 Fisheries Research 90:147-157",
             "Duncan et al. (2015) Fisheries Research 164:64-72", "Duncan et al. 2015 Fisheries Research 164:64-72", 
             "Durand et al. 2005 Mar Bio", "Durand et al. 2005 Mar. Bio. 147:313-322",
             "Fauvelot & Borsa (2011) Biological Journal of the Linnean Society 104:886-902", "Fauvelot & Borso 2001 Biological Journal of the Linnean Society 104: 886-902",
             "Fernandez-Silva 2015 J. Biogeogr. 42:2402-2413", "Fernandez-Silva et al. (2015) J. Biogeogr. 42:2402-2413",
             "Gaither et al. (2011) PLoS ONE 6:e28913", "Gaither et al. 2011 PLoS ONE 6:1-13",
             "Galarza et al. 2009 PNAS", "Galerza et al. 2009 PNAS 106:1473-1478",
             "Gold et al. 2010 Fish Res", "Gold et al. 2010 Fisheries Research 103:30-39",
             "Gold et al. (2011) North American Journal of Fisheries Management 31:209-223", "Gold et al. 2011 North American Journal of Fisheries Management 31:209-223", 
             "González-Wangüemert & Pérez-Ruzafa (2012) Marine Ecology 33:337-349", "Gonzalez-Wanguemert & Perez-Ruzafa 2012 Marine Ecology 33:337-349",
             "Haney Silliman Rand 2007", "Haney, Silliman, & Rand 2007 Heredity 98:294-302",
             "Healey et al. (2018) Ecology and Evolution 8:2182-2195", "Healey et al. 2017 Ecology & Evolution 10.1002/ece3.3775", 
             "Healey et al. (2018) ICES Journal of Marine Science 75:1465-1472", "Healey et al. 2018 ICES Journal of Marine Science 75:1465-1472",
             "Henriques et al. (2017) Fisheries Research 187:86-95", "Henriques et al. 2017 Fisheries Research 187:86-97",
             "Hess et al. 2001 CJFAS 68:89-104", "Hess et al. 2011 Canadian Journal of Fisheries and Aquatic Sciences 68:89-104", 
             "Horne et al. (2013) Journal of Heredity 104:532-546", "Horne et al. 2013 Journal of Heredity 104:532-546",
             "Jackson et al. (2014) PLoS ONE 9:e97508", "Jackson et al. 2014 PLoS ONE 9:e97508",
             "Karl et al. 2012 Marine Biology 159:489-498", "Karl, Castro & Garla (2012) Mar. Biol. 159:489-498",
             "Kim et al. 2007 Mol Ecol Notes", "Kim, Morishima & Arai 2007 Mol Ecol Notes 7:79-81",
             "Lane, Symonds & Ritchie (2016) Fisheries Research 174:19-29", "Lane, Symonds & Ritchie 2016 Fisheries Research 174:19-29",
             "Limborg et al. (2012) Heredity 109:96-107", "Limborg et al. 2012 Heredity 109:96-107",
             "Tripp-Valdez et. al 2012 Journal of Applied Ichthyology 28:516-523", "M. A. Tripp-Valdez et al. (2012) J. Appl. Ichthyol. 28: 516â€“523",
             "Magallón-Gayón, Diaz-Jaimes & Uribe-Alcocer (2016) PeerJ 4:e2583", "Magallón-Gayón, Diaz-Jaimes & Uribe-Alcocer 2016 PeerJ 4:e2583",
             "Mzingirwa et al. (2019) Molecular Biology Reports 46:5079-5088", "Mzingirwa et al. 2019 Molecular Biology Reports 46:5079-5088",
             "Nance et al. (2011) PLoS ONE 6:e21459", "Nance et. al 2011 PLoS One 6(7):1-12", 
             "Nielsen et al. 2010 Cons Gen 11:999-1012", "Nielsen et. al 2010 Conservation Genetics 11:999-1012", 
             "Oretga-Villaizan et al. 2006 Fisheries Sci 72:556-567", "Ortega Villaizan Romo et al. 2006", 
             "Piñeros et al. (2015) Evol. Biol. 42:235-249", "Piñeros et al. 2015 Evol. Biol. 42:235-249", 
             "Renshaw et al. 2006 Mol Ecol Notes 6:1162-1164", "Renshaw et al. 2007 Mol Ecol Notes", 
             "Ruggeri et al. (2016) PLoS ONE 11:e0151507", "Ruggeri et al. (2016) PLoS ONE 11:e0153061", 
             "Sala-Bozano et al. 2009 Mol Ecol", "Sala-Bozano et. al 2009 Molecular Ecology 18:4811-4826", 
             "Šegvi?-Bubi? et al. (2016) Fisheries Research 179:271-279", "Šegvi?-Bubi? et al. 2016 Fisheries Research 179:271-279", 
             "Sekino et al. (2011) Conserv. Genet. 12:139-159", "Sekino et. al 2011 Conservation Genetics 12:139-159", 
             "Silva et al. (2017) Scientific Reports 7:2983", "Silva et al. 2017 Scientific Reports 7:2893", 
             "Silva, Horne & Castilho (2014) J. Biogeogr. 41:1171-1182", "Silva, Horne & Castilho 2014 J. Biogeogr. 41:1171-1182", 
             "Stefanni et al. 2015 Heredity 115:527-537", "Steffani et al. (2015) Heredity 115:527-537",
             "Steinberg et al. (2016) Coral Reefs 35:959-970", "Steinberg et al. 2016 Coral Reefs 35:959-970", 
             "Teske et al. 2010 Mar Biol 157:2029-2042", "Teske et. al 2010 Marine Biology 157:2029-2042", 
             "Vignaud et al. (2014) Molecular Ecology 23: 2590â€“2601", "Vignaud et al. 2014 Molecular Ecology 23:2590-2601", 
             "Waldrop et al. (2016) Journal of Biogeography (J. Biogeogr.) 43, 1116â€“1129", "Waldrop et al. 2016 J. Biogeogr. 43:1116-1129", 
             "Wilson 2006", "Wilson 2006 Mol Ecol 15:809-824")

#pull matching studies from mtdna & msat
msat_matches <- as.data.table(msat[(msat$Source %in% matches), ])
mtdna_matches <- as.data.table(mtdna[(mtdna$Source %in% matches), ])

#take mean of msat markers
msat_means <- msat_matches[, mean(He, na.rm = TRUE), by = .(spp, Source, Site)]
  msat_source_sum <- msat_means[, .N, by = Source] #count # sites by sources

#take mean of mtdna markers (He)
mtdna_means <- mtdna_matches[, mean(He, na.rm = TRUE), by = .(spp, Source, Site, MarkerName)]
  mtdna_source_sum <- mtdna_means[, .N, by = Source] #count # sites by sources

#check where sitenames may differ
check_sitenames_mtdna <- c(mtdna_means[spp == "Chaetodon lunulatus", 3 ])
  check_sitenames_mtdna
check_sitenames_msat <- c(msat_means[spp == "Verasper variegatus", 3 ])
  check_sitenames_msat

#fix site names (make match across dataframes)
mtdna_means$Site[mtdna_means$spp == "Anoplopoma fimbria" & mtdna_means$Site == "SQ"] <- "San Quintin"
  mtdna_means$Site[mtdna_means$spp == "Anoplopoma fimbria" & mtdna_means$Site == "BS"] <- "Bering Sea"
  mtdna_means$Site[mtdna_means$spp == "Anoplopoma fimbria" & mtdna_means$Site == "GA"] <- "Gulf of Alaska"
  mtdna_means$Site[mtdna_means$spp == "Anoplopoma fimbria" & mtdna_means$Site == "OR"] <- "Oregon"
msat_means$Site[msat_means$spp == "Chaetodon lunulatus" & msat_means$Site == "Moâ€™orea"] <- "Mo'orea"
  msat_means$Site[msat_means$spp == "Chaetodon lunulatus" & msat_means$Site == "Okinawa"] <- "Okinawa Island"
mtdna_means$Site[mtdna_means$spp == "Diplodus sargus" & mtdna_means$Site == "Quarteira "] <- "Quarteira"
  mtdna_means$Site[mtdna_means$spp == "Diplodus sargus" & mtdna_means$Site == "Galicia "] <- "Galicia"
  mtdna_means$Site[mtdna_means$spp == "Diplodus sargus" & mtdna_means$Site == "Murica "] <- "Murcia"
mtdna_means$Site[mtdna_means$spp == "Diplodus vulgaris" & mtdna_means$Site == "Portugal S."] <- "Portugal South"
  mtdna_means$Site[mtdna_means$spp == "Diplodus vulgaris" & mtdna_means$Site == "Madeira Island"] <- "Madeira Islands"
  mtdna_means$Site[mtdna_means$spp == "Diplodus vulgaris" & mtdna_means$Site == "Balearic Island"] <- "Balearic Islands"
mtdna_means$Site[mtdna_means$spp == "Epinephelus striatus" & mtdna_means$Site == " Dog Rocks"] <- "Dog Rocks"
mtdna_means$Site[mtdna_means$spp == "Ginglymostoma cirratum" & mtdna_means$Site == "BMN"] <- "Bimini"
  mtdna_means$Site[mtdna_means$spp == "Ginglymostoma cirratum" & mtdna_means$Site == "BLZ"] <- "Glover's Reef"
  mtdna_means$Site[mtdna_means$spp == "Ginglymostoma cirratum" & mtdna_means$Site == "BRI"] <- "Brazilian Islands"
  mtdna_means$Site[mtdna_means$spp == "Ginglymostoma cirratum" & mtdna_means$Site == "NTL"] <- "Natal"
  mtdna_means$Site[mtdna_means$spp == "Ginglymostoma cirratum" & mtdna_means$Site == "FLD"] <- "Florida"
mtdna_means$Site[mtdna_means$spp == "Hippoglossus stenolepis" & mtdna_means$Site == "Gulf of Alaska"] <- "Gulf of Alaska, AK"
mtdna_means$Site[mtdna_means$spp == "Lithognathus mormyrus" & mtdna_means$Site == "BA 06"] <- "L'Estartit 06"
  mtdna_means$Site[mtdna_means$spp == "Lithognathus mormyrus" & mtdna_means$Site == "AEG 08"] <- "Iraklio 08"
  mtdna_means$Site[mtdna_means$spp == "Lithognathus mormyrus" & mtdna_means$Site == "IND 08"] <- "Durban 08"
  mtdna_means$Site[mtdna_means$spp == "Lithognathus mormyrus" & mtdna_means$Site == "ADR 07"] <- "Duce 07"
  mtdna_means$Site[mtdna_means$spp == "Lithognathus mormyrus" & mtdna_means$Site == "ALB 07"] <- "Malaga 07"
mtdna_means$Site[mtdna_means$spp == "Lutjanus analis" & mtdna_means$Site == "Puerto Rico-East"] <- "Puerto Rico - East"
  mtdna_means$Site[mtdna_means$spp == "Lutjanus analis" & mtdna_means$Site == "Puerto Rico-West"] <- "Puerto Rico - West"
mtdna_means$Site[mtdna_means$spp == "Lutjanus synagris" & mtdna_means$Site == "Puerto Rico-west "] <- "Puerto Rico - West"
  mtdna_means$Site[mtdna_means$spp == "Lutjanus synagris" & mtdna_means$Site == "Puerto Rico-east "] <- "Puerto Rico - East"
  mtdna_means$Site[mtdna_means$spp == "Lutjanus synagris" & mtdna_means$Site == "Florida Keys (Marathon)"] <- "Florida Keys"
  mtdna_means$Site[mtdna_means$spp == "Lutjanus synagris" & mtdna_means$Site == "St. Croix "] <- "St. Croix"
  mtdna_means$Site[mtdna_means$spp == "Lutjanus synagris" & mtdna_means$Site == "St. Thomas "] <- "St. Thomas"
mtdna_means$Site[mtdna_means$spp == "Pagrus pagrus" & mtdna_means$Site == "Mediterranean/Crete"] <- "Crete"
mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Christmas Island "] <- "Christmas Island"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Rowley Shoals "] <- "Rowley Shoals"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Scott Reef "] <- "Scott Reef"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Tonga "] <- "Tonga"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Seychelles "] <- "Seychelles"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Brooks Banks "] <- "Brooks Banks"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Gardner "] <- "Gardner"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Guam "] <- "Guam"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Hawai'i Island "] <- "Hawai'i Island"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Kaua'i "] <- "Kaua'i"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Lana'i "] <- "Lana'i"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Maro "] <- "Maro"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Maui "] <- "Maui"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Moloka'i "] <- "Moloka'i"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "North Hampton "] <- "North Hampton"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Penguin Banks "] <- "Penguin Banks"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "Pioneer "] <- "Pioneer"
  mtdna_means$Site[mtdna_means$spp == "Pristipomoides filamentosus" & mtdna_means$Site == "St. Rogatien "] <- "St. Rogatien"
mtdna_means$Site[mtdna_means$spp == "Scomberomorus brasiliensis" & mtdna_means$Site == "North Coast Trinidad "] <- "North Coast"
  mtdna_means$Site[mtdna_means$spp == "Scomberomorus brasiliensis" & mtdna_means$Site == "South Coast Trinidad "] <- "South Coast"
  mtdna_means$Site[mtdna_means$spp == "Scomberomorus brasiliensis" & mtdna_means$Site == "West Coast Trinidad "] <- "West Coast"
  mtdna_means$Site[mtdna_means$spp == "Scomberomorus brasiliensis" & mtdna_means$Site == "Cumana "] <- "Cumana"
  mtdna_means$Site[mtdna_means$spp == "Scomberomorus brasiliensis" & mtdna_means$Site == "Isla de Margarita "] <- "Isla de Margarita"
mtdna_means$Site[mtdna_means$spp == "Scomberomorus commerson" & mtdna_means$Site == "Coral Sea"] <- "New Caledonia"
mtdna_means$Site[mtdna_means$spp == "Sebasates flavidus" & mtdna_means$Site == "N. Vancouver Is., BC"] <- "N. Vancouver Island, BC"
  mtdna_means$Site[mtdna_means$spp == "Sebasates flavidus" & mtdna_means$Site == "S. Vancouver Is., BC"] <- "S. Vancouver Island, BC"
  mtdna_means$Site[mtdna_means$spp == "Sebasates flavidus" & mtdna_means$Site == "San Miguel Is., CA"] <- "San Miguel Island, CA"
mtdna_means$Site[mtdna_means$spp == "Solea solea" & mtdna_means$Site == "Celtic Sea CEL08"] <- "Celtic Sea"
  mtdna_means$Site[mtdna_means$spp == "Solea solea" & mtdna_means$Site == "English Channel ENG08"] <- "Eng_English Channel"
  mtdna_means$Site[mtdna_means$spp == "Solea solea" & mtdna_means$Site == "North Sea NOR08"] <- "Ger_North Sea"
  mtdna_means$Site[mtdna_means$spp == "Solea solea" & mtdna_means$Site == "Scheldt estuary ZAN07"] <- "Scheldt estuary"
  mtdna_means$Site[mtdna_means$spp == "Solea solea" & mtdna_means$Site == "Wadden Sea TEX07"] <- "Wadden Sea"
mtdna_means$Site[mtdna_means$spp == "Sparus aurata" & mtdna_means$Site == "France"] <- "Ille d'Oleron"
  mtdna_means$Site[mtdna_means$spp == "Sparus aurata" & mtdna_means$Site == "Ireland"] <- "Wexford"
  mtdna_means$Site[mtdna_means$spp == "Sparus aurata" & mtdna_means$Site == "Portugal"] <- "Aveiro"
  mtdna_means$Site[mtdna_means$spp == "Sparus aurata" & mtdna_means$Site == "Spain"] <- "Cadiz"
mtdna_means$Site[mtdna_means$spp == "Sphyrna lewini" & mtdna_means$Site == "Michoac’n"] <- "Michoacán"
  mtdna_means$Site[mtdna_means$spp == "Sphyrna lewini" & mtdna_means$Site == "LAP"] <- "La Paz"
  mtdna_means$Site[mtdna_means$spp == "Sphyrna lewini" & mtdna_means$Site == "TAR"] <- "Tarcoles"
  mtdna_means$Site[mtdna_means$spp == "Sphyrna lewini" & mtdna_means$Site == "MAZ"] <- "Mazatlan"
  mtdna_means$Site[mtdna_means$spp == "Sphyrna lewini" & mtdna_means$Site == "CEB"] <- "Cebaco Island"
  mtdna_means$Site[mtdna_means$spp == "Sphyrna lewini" & mtdna_means$Site == "GPA"] <- "Gulf of Panama"
  mtdna_means$Site[mtdna_means$spp == "Sphyrna lewini" & mtdna_means$Site == "SCA"] <- "Santa Catalina"
mtdna_means$Site[mtdna_means$spp == "Syngnathus leptorhynchus" & mtdna_means$Site == "San Diego"] <- "San Diego Bay, CA"
mtdna_means$Site[mtdna_means$spp == "Thalassoma bifasciatum" & mtdna_means$Site == "Butler Bay"] <- "Butler Bay, St Croix"
  mtdna_means$Site[mtdna_means$spp == "Thalassoma bifasciatum" & mtdna_means$Site == "Grapetree Bay"] <- "Grapetree Bay, St Croix"
  mtdna_means$Site[mtdna_means$spp == "Thalassoma bifasciatum" & mtdna_means$Site == "Jack Bay"] <- "Jack Bay, St Croix"
  mtdna_means$Site[mtdna_means$spp == "Thalassoma bifasciatum" & mtdna_means$Site == "Northstar Bay"] <- "Northstar Bay, St Croix"
mtdna_means$Site[mtdna_means$spp == "Thunnus albacares" & mtdna_means$Site == "KK"] <- "Kandakuliya"
  mtdna_means$Site[mtdna_means$spp == "Thunnus albacares" & mtdna_means$Site == "KR"] <- "Kirinda"
  mtdna_means$Site[mtdna_means$spp == "Thunnus albacares" & mtdna_means$Site == "MD"] <- "Maldives"
  mtdna_means$Site[mtdna_means$spp == "Thunnus albacares" & mtdna_means$Site == "NE"] <- "Negombo"
  mtdna_means$Site[mtdna_means$spp == "Thunnus albacares" & mtdna_means$Site == "TA"] <- "Tangalle"
  mtdna_means$Site[mtdna_means$spp == "Thunnus albacares" & mtdna_means$Site == "TR"] <- "Trincomalee"
  mtdna_means$Site[mtdna_means$spp == "Thunnus albacares" & mtdna_means$Site == "WE"] <- "Weligama"
mtdna_means$Site[mtdna_means$spp == "Verasper variegatus" & mtdna_means$Site == "CYS"] <- "Yellow Sea"
  mtdna_means$Site[mtdna_means$spp == "Verasper variegatus" & mtdna_means$Site == "ISB"] <- "Ise Bay"
  
#merge by site/species
match_list <- list(mtdna_means, msat_means)
match_all <- match_list %>% reduce(full_join, by = c("spp", "Site"), all = TRUE)
  colnames(match_all) <- c("spp", "Source_mtdna", "Site", "MarkerName", "Hd", "Source_msat", "He")  
  
  
match_all_noNA <- na.omit(match_all)
cor(match_all_noNA$He, match_all_noNA$Hd) #0.459 --> in general He (msat) higher than Hd
  
#plot Hd v. He
Hd_He_plot <- ggplot(data = match_all, aes(x = He, y = Hd)) + 
  geom_point(aes (color = spp)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#666666") + 
 geom_smooth() + 
  theme_bw() + 
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1))
Hd_He_plot

#plot Hd v. pi
Hd_pi_plot <- ggplot(data = mtdna, aes(x = He, y = Pi)) + 
  geom_point(aes (color = spp), show.legend = FALSE) + 
  geom_smooth() + 
  theme_bw()
Hd_pi_plot
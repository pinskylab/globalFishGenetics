# Microsatellite Data Entry
---

Guidelines for how to enter data in microsatellite data spreadsheets (including filtering criteria).

## Screening Instructions

We only want to enter data from population genetic studies on marine fish that use microsatellite loci, report genetic diversity, and report where the samples were taken.

**Do NOT enter data from:**
1. Captive, farmed, or stocked populations, or anything else that isn't wild.
2. Hybrid populations or samples.
    * If the study doesn't mention this explicitly, it is probably OK.
4. EST (Expressed Sequence Tag)-linked microsatellites (keywords indicating this: linkage groups, family groups).
    * If a study uses EST-linked microsatellites, it will probably say so in teh abstract. If it doesn't mention anything, the study is probably OK. ESTs are bits of genetic material copied from an organisms' transcriptome (the mRNA) and so is likely to be closely associated with genes that matter to the organisms' fitness and more likely to be experiencing natural selection.
5. Loci that are not in Hardy-Weinberg equilibrium (they may be under selection).
6. Loci that are monomorphic.
7. Anything that isn't a microsatellite (RFLPs and minisatellites are not microsatellites).
8. Species that are anadromous, catadromous, amphidromous, estuarine, or freshwater (anything that isn't marine).
    * If you're not sure if the species is marine, enter the species name in www.fishbase.org or www.fishbase.us and read about the species.
9. Data on diversity averaged across a wide area (> 3° latitude or longitude) or where the site is only vaguely identified (accuracy worse than 3° latitude or longitude).
10. Multiple samples from the same site (for example, in different years if each year is reported separately). To decide which year to enter, pick the one with the highest sample size, and to break a tie, use the most recent sample.
11. Samples with a sample size smaller than 4 individuals.
12. Don't enter the same data twice. For example, a paper may erport samples and data that were first analyzed and presented in another paper that you have or will enter. If the data are the same, enter the data from the original paper.

## Entry Instructions

Once you are certain the study passes the above criteria, we can record data. In the spreadsheet for entering data, each line is the diversity for one microsatellite at one site (population) in one species. If the paper only reports averages across multiple microsatellite loci, then the line will have that average.

If you have to choose between entering averaged data for many microsatellites at one site, or averaged data for many sites at one microsatellite, choose the average across microsatellites at one site. We care more about the spatial variation in diversity.

**Data to enter:**
1. **Spp:** Scientific name
    * Check on FishBase to make sure you have the most up-to-date scientific name, as it may differ from what is on the paper.
2. **CommonName:** Also from FishBase
3. **Source:** Study data comes from - enter in the form *Adams et al. 2010 Molecular Ecology 8:1000-1015*, where 8 is the volume and 1000-1015 are the page numbers
    * If only one or two authors, use a citation similar to Adams, or Adams & Benner, respectively.
4. **PrimerNote:** 1 for studies which are Primer Notes, 0 for all others
    * Some studies only exist to report the discovery of new microsatellites; these are called Primer Notes. Often they're labeled as such, but not always. They tend to be short and only evaluate a small number of individuals. Many of them are in the journals *Molecular Ecology Notes* or *Conservation Genetics Resources*.
5. **Country:** Country where the sample was taken
6. **Site:** Site where sample was taken, as named by the author
7. **Lat_Deg:** Degrees latitude (can be decimal degrees)
    * Use negative numbers for the Southern Hemisphere.
    * Find the site on Google Maps if needed. Right click and select "What's Here?" to get Google to give you the lat and lon coordinates.
8. **Lat_Min:** Minutes latitude, use only if you don't have decimal degrees for Lat_Deg
9. **Lat_Sec:** Seconds latitude, use only if you don't have decimal degrees for Lat_Deg
10. **Lon_Deg:** Degrees longitude (can be decimal degrees)
    * Use negative numbers for the Western Hemisphere.
    * Find the site on Google Maps if needed. Right click and select "What's Here?" to get Google to give you the lat and lon coordinates.
11. **Lon_Min:** Minutes longitude, use only if you don't have decimal degrees for Lon_Deg
12. **Lon_Sec:** Seconds longitude, use only if you don't have decimal degrees for Lon_Deg
13. **CollectionYear:** Year in which the samples were taken
14. **NumMarkers:** Number of microsatellite markers whose data are entered on this line
    * Will often be 1 (will only be >1 if the paper only reports an average across multiple loci).
15. **MarkerName:** Name of the microsatellite marker, as listed in the paper
    * Leave blank if there are multiple markers reported together on this line.
16. **CrossSpp:** Was the microsatellite originally developed in a different specie? 1 if yes, 0 if no
    * You may have to find the original reference for a microsatellite locus to learn whether it was developed in a different species.
    * If the line of data is for a number of markers averaged together and some were cross-species and some were not, input the proportion that were cross-species (e.g.: 0.2 if 1 of 5 was cross-species).
17. **n:** Number of individuals sampled
    * If NumMarkers >1, then this will be an average across all of the markers reported on this line.
18. **Repeat:** The length of the microsatellite repeat
    * Will often be 2-6. A 2-nucleotide repeat is called dinucleotide, 3 is called trinucleotide, etc.
    * If the repeat is complex (has di's, tri's, tetra's or some mix), then use the most common repeat motif length. For example, (GA7)(CACA)9(CT)6 would be 2, since 7+6 > 9.
    * If NumMarkers >1, then this will be an average across all the repeat sizes of those markers.
19. **He:** Expected heterozgosity
    * This will often be reported in a table or sometimes an appendix/supplemental material for the paper.
20. **Hese:** Standard error of the He measurement
    * Only enter this if NumMarkers >1 and the paper reports it. You can convert standard deviation to standard error by dividing by sqrt(n).
21. **P1 to P68:** Extra columns that can be added to hold allele frequencies, if He is not reported
    * One column for each allele. If allele numbers are reported rather than frequencies, divide by 2n to get frequencies.

**Each line MUST have Spp, Source, NumMarkers, n, He, and LatDeg!!** The other cells are important to have, but are not reasons for skipping the line of data.

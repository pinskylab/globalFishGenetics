# Microsatellite data entry
2/11/2015

Do not enter data from:
1. Captive, farmed, or stocked populations, or anything else that isn’t wild
2. EST-linked microsatellites (keywords: linkage groups, family groups)
3. Loci that are not in Hardy-Weinberg equilibrium (they may be under selection)
4. Anything that isn’t a microsatellite (minisatellites are not microsatellites)
5. Species that are anadromous, catadromous, or freshwater (anything that isn’t marine)
6. Microsatellites that have only one allele (monomorphic).
7. Data on diversity averaged across a wide area (> 3° latitude or longitude) or where the site is only vaguely identified (accuracy worse than 3° latitude).
8. Papers that I’ve already entered (see list in SourcesEntered_2012-02-05.csv). Note, however, that there are a few large primer notes where I’ve only entered data from some of the species listed in the primer note. 
9. Multiple samples from the same site (for example, in different years if each year is reported separately). To decide which year to enter, pick the one with the highest sample size, and to break a tie, use the most recent sample. 
10. Don’t enter the same data twice. For example, a paper may report samples and data that were first analyzed and presented in another paper that you have or will enter. If the data are the same, enter the data from the original paper.

## Step 3. Enter data
We have a separate spreadsheet for entering data. Each line is the diversity for one microsatellite at one site. If the paper only reports averages across multiple microsatellite loci, then the line will have that average.

If you have to choose between entering averaged data for many microsatellites at one site, or averaged data for many sites at one microsatellite, choose the average across microsatellites at one site. We care more about the spatial variation in diversity.

Data to enter:
1. Spp
a. Scientific name. Check on Fishbase to make sure you have the most up-to-date scientific name, as it may differ from what is in the paper.
2. CommonName
a. Also from Fishbase
3. Source
a. Enter in the form “Adams et al. 2010 Molecular Ecology 8:1000-1015,” where 8 is the volume and 1000-1015 are the page numbers.
b. If only one or two authors, use a citation similar to Adams, or Adams & Benner, respectively.
4. PrimerNote
a. 1 for Primer Note, 0 for not.
b. Some studies only exist to report the discovery of new microsatellites, and these are called Primer Notes. Often they’re labeled as such, but not always. They tend to be short and only evaluate a small number of individuals. Many of them are in the journals Molecular Ecology Resources, Molecular Ecology Notes, or Conservation Genetics Resources.
5. Country
a. Country where the sample was taken
6. Site
a. Site name, as named by the author
7. Lat_Deg
a. Degrees latitude. This can be decimal degrees. Use negative numbers for the Southern Hemisphere
b. Find the site on Google Maps if needed. Right click and select ‘What’s Here?’ to get Google to give you the lat and long.
8. Lat_Min
a. Minutes latitude. Use this field only if you don’t have decimal degrees for LatDeg
9. Lat_Sec
a. Seconds latitude. Use this field only if you have neither decimal degrees nor decimal minutes for Lat_Deg and Lat_Min
10. Lon_Deg
a. Degrees longitude. This can be decimal degrees. Use negative numbers for the Western Hemisphere
b. Find the site on Google Maps if needed. Right click and select ‘What’s Here?’ to get Google to give you the lat and long.
11. Lon_Min
a. Minutes longitude. Use this field only if you don’t have decimal degrees for LonDeg
12. Lon_Sec
a. Seconds longitude. Use this field only if you have neither decimal degrees nor decimal minutes for Lon_Deg and Lon_Min
13. CollectionYear
a. Year in which the samples were taken. Leave blank if not stated.
14. NumMarkers
a. Number of microsatellite markers whose data are entered on this line. Will often be 1. It will be >1 if the paper only reports an average across multiple loci.
15. MarkerName
a. Name of the microsatellite marker, as listed in the paper.
b. Leave blank if there are multiple markers reported together on this line.
16. CrossSpp
a. Was the microsatellite originally developed in a different species?
b. 1 if yes, 0 if no
c. You may have to find the original reference for a microsatellite locus to learn whether it was developed in a different species.
d. If the line of data is for a number of  markers averaged together and some were cross-species and some were not, input the proportion that were cross-species (e.g., 0.2 if 1 of 5 was cross-species).
17. n
a. Number of individuals sampled
b. If NumMarkers>1, then this might have to be an average across all of the markers reported on this line.
18. Repeat
a. The length of the microsatellite repeat. Will be 2-5. A 2-nucleotide repeat is a called dinucleotide, 3 is called trinucleotide, etc. If the repeat is complex (has di’s, tri’s, and tetra’s, or some mix), then use the most common repeat motif length. For example, (GA)7(GACA)9(CT)6 would be 2, since 7+6 > 9. 
b. If NumMarkers>1, then this will be an average across all the repeat sizes of those markers.
19. He
a. Expected heterozygosity
b. This will often be reported in a table, or sometimes an appendix for the paper.
20. Hese
a. Standard error of the He measurement
b. Only enter this if NumMarkers > 1 and the paper reports it
c. You can convert standard deviation to standard error by dividing by sqrt(number of values averaged together).
21. P1 to P68
a. Extra columns that can be added to hold allele frequencies, if He is not reported. One column for each allele
b. If allele numbers are reported rather that frequencies, divide by 2n to get frequency.

Each line MUST have Spp, Source, NumMarkers, n, He, and LatDeg. The other cells are important to have, but are not reasons for skipping the line of data.


## Keep notes!
Keep a dated logbook (a text file is good) where you make notes

## Saving your file
Save your file often. 

When it’s time to save your file, title it:

Fishery lat msats A III YYYY-MM-DD.xlsx

where A is the number for the batch of files, III are your initials, YYYY is the year, MM is the month, and DD is the day.

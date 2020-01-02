# mtDNA data entry
7/15/2015

Our goal is to extract key data from a large number of published scientific papers and enter it into a large database. This database contains information on patterns of genetic diversity for marine fish globally. Nothing like this has been put together before.

There are three key steps when looking at a paper: 1) Determine whether we can use the paper, 2) Extract the data, and 3) Document what you did in a logbook. Sometimes, you will be able to use some but not all of the data from a paper.

I fully expect that some of the papers you read or aspects of these instructions will be confusing or foreign. PLEASE ask me when you have questions (malin.pinsky@rutgers.edu). I would much prefer that you ask me questions than enter data incorrectly.

## Step 1. Screening: Do not enter data from:
1. Captive, farmed, or stocked populations, or anything else that isn’t wild
2. Hybrids between different species
3. Populations with a sample size of less than 4 individuals
4. Anything that isn’t mitochondrial DNA of at least 200 bp
5. Any loci that are RFLPs rather than sequencing
6. Species that are anadromous, catadromous, or freshwater (anything that isn’t marine)
7. Papers with only one sampled site
8. Data from populations that have no lat/long, and for which you can’t figure out the lat/long.
9. Data on diversity averaged across a wide area (> 3° latitude or longitude) or where the site is only vaguely identified (accuracy worse than 3° latitude).
10. Multiple samples from the same site (for example, in different years if each year is reported separately). To decide which year to enter, pick the one with the highest sample size, and to break a tie, use the most recent sample. 
11. Don’t enter the same data twice. For example, a paper may report samples and data that were first analyzed and presented in another paper that you have or will enter. If the data are the same, enter the data from the original paper.


## Step 2. Enter data
We have a separate spreadsheet for entering data. Each line is the diversity for one mtDNA locus at one site.

Data to enter:
1. Spp
a. Scientific name. Check on Fishbase (http://www.fishbase.org) to make sure you have the most up-to-date scientific name, as it may differ from what is in the paper.
2. CommonName
a. Also from Fishbase
b. If nothing in Fishbase, then use a common name mentioned in the paper.
3. Source
a. Enter in the form “Adams et al. 2010 Molecular Ecology 8:1000-1015,” where 8 is the volume and 1000-1015 are the page numbers.
b. If only one or two authors, use a citation similar to “Adams”, or “Adams & Benner”, respectively.
4. Country
a. Country where the sample was taken
b. If the location is far from shore and can’t be tagged to a country, please use the ocean name (e.g., Pacific Ocean). 
5. Site
a. Site name, as named by the author
6. LatDeg
a. Degrees latitude. This can be decimal degrees. Use negative numbers for the Southern Hemisphere
b. If exact latitude is not given, use any maps in the paper to narrow in on the site location. Place names may also be helpful. However, if you don’t think that you can nail down the location more precisely than 3 degrees, we will have to skip the data (see criterion for exclusion above).
c. Find the site on Google Maps if needed. Right click and select ‘What’s Here?’ to get Google to give you the lat and long.
d. If a range of latitudes are reported by the paper, use the average (as long as <= 3 degrees, see criterion for exclusion above).
7. LatMin
a. Minutes latitude. Use this field only if you don’t have decimal degrees for LatDeg
8. LatSec
a. Seconds latitude. Use this field only if you don’t have decimal degrees for LatDeg
9. LonDeg
a. Degrees longitude. This can be decimal degrees. Use negative numbers for the Western Hemisphere
b. Find the site on Google Maps if needed. Right click and select ‘What’s Here?’ to get Google to give you the lat and long.
c. If a range of longitude are reported by the paper, use the average (as long as <= 3 degrees, see criterion for exclusion above).
10. LonMin
a. Minutes longitude. Use this field only if you don’t have decimal degrees for LonDeg
11. LonSec
a. Seconds longitude. Use this field only if you don’t have decimal degrees for LonDeg
12. CollectionYear
a. Year in which the samples were taken. Leave blank if not stated. Average if samples from multiple years were pooled.
13. MarkerName
a. Name of the mtDNA locus, as listed in the paper, e.g., “control region” or “cytochrome b”
14. n
a. Number of individuals sampled
15. bp
a. The length of the mtDNA locus, in basepairs.
16. He
a. Expected heterozygosity, aka the haplotype diversity, and perhaps other names often denoted by H or h (check first with me).
b. This will often be reported in a table, or sometimes an appendix for the paper.
17. Hese
a. Standard error of the He measurement
b. You can convert standard deviation to standard error by dividing by sqrt(number of values averaged together).
18. Pi
a. The nucleotide diversity
19. Pise
a. The standard error of pi
b. You can convert standard deviation to standard error by dividing by sqrt(number of values averaged together).

Each line MUST have Spp, Source, (He or pi), and LatDeg. The other cells are important to have, but are not reasons for skipping the line of data.


## Step 3. Keep a detailed logbook
Keep a dated logbook (a Google Docs file shared with me is good) where you make notes.  Each day that you do work, you should enter the day’s date on a new line, and then list the papers that you worked on that day. Under each paper, keep notes on whether or not you entered any data from paper, why not (if not entered), what was confusing in the paper, and how you resolved your confusion. If you have questions on any piece, please email me and then make notes in your logbook how you resolved your confusion. 

## Step 4. Saving your file
Save your file often. Start over with a new file whenever you hit 2000 lines.

When it’s time to save your file, title it:

Fishery lat mtDNA A YYYY-MM-DD III.xlsx

where A is the paper set number, III are your initials, YYYY is the year, MM is the month, and DD is the day.



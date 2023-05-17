# mtDNA Data Entry
---

Guidelines for how to enter data in mitochondrial DNA (mtDNA) data spreadsheets (including filtering criteria).

## Screening Instructions

We only want to enter data from population genetic studies on marine fish that use mtDNA sequences, report haplotype sequences, report genetic diversity, and report where the samples were taken.

**Do NOT enter data from:**
1. Captive, farmed, or stocked populations, or anything else that isn't wild.
2. Hybrid populations or samples.
    * If the study doesn't mention this explicitly, it is probably OK.
3. Anything that isn't from a mitochondrial DNA sequence of at least 200 bp in length.
4. Loci that are monomorphic.
5. Anything that isn't mitochondrial DNA sequence data (RFLPs do not count).
6. Species that are anadromous, catadromous, amphidromous, estuarine, or freshwater (anything that isn't marine).
    * If you're not sure if the species is marine, enter the species name in www.fishbase.org or www.fishbase.us and read about the species.
7. Data on diversity averaged across a wide area (> 3° latitude or longitude) or where the site is only vaguely identified (accuracy worse than 3° latitude or longitude).
8. Multiple samples from the same site (for example, in different years if each year is reported separately). To decide which year to enter, pick the one with the highest sample size, and to break a tie, use the most recent sample.
9. Samples with a sample size smaller than 4 individuals.
10. Don't enter the same data twice. For example, a paper may erport samples and data that were first analyzed and presented in another paper that you have or will enter. If the data are the same, enter the data from the original paper.

## Entry Instructions

Once you are certain the study passes the above criteria, we can record data. In the spreadsheet for entering data, each line is the diversity for one mtDNA locus at one site (population) in one species.

**Data to enter:**
1. **Spp:** Scientific name
    * Check on FishBase to make sure you have the most up-to-date scientific name, as it may differ from what is on the paper.
2. **CommonName:** Also from FishBase
3. **Source:** Study data comes from - enter in the form *Adams et al. 2010 Molecular Ecology 8:1000-1015*, where 8 is the volume and 1000-1015 are the page numbers
    * If only one or two authors, use a citation similar to Adams, or Adams & Benner, respectively.
4. **Country:** Country where the sample was taken
5. **Site:** Site where sample was taken, as named by the author
6. **Lat_Deg:** Degrees latitude (can be decimal degrees)
    * Use negative numbers for the Southern Hemisphere.
    * Find the site on Google Maps if needed. Right click and select "What's Here?" to get Google to give you the lat and lon coordinates.
7. **Lat_Min:** Minutes latitude, use only if you don't have decimal degrees for Lat_Deg
8. **Lat_Sec:** Seconds latitude, use only if you don't have decimal degrees for Lat_Deg
9. **Lon_Deg:** Degrees longitude (can be decimal degrees)
    * Use negative numbers for the Western Hemisphere.
    * Find the site on Google Maps if needed. Right click and select "What's Here?" to get Google to give you the lat and lon coordinates.
10. **Lon_Min:** Minutes longitude, use only if you don't have decimal degrees for Lon_Deg
11. **Lon_Sec:** Seconds longitude, use only if you don't have decimal degrees for Lon_Deg
12. **CollectionYear:** Year in which the samples were taken
13. **MarkerName:** Name of the mtDNA locus, as listed in the paper (e.g. "control region" or "cytochrome b")
    * Leave blank if not stated.
14. **n:** Number of individuals sampled
15. **bp:** The length of the mtDNA locus, in basepairs
16. **He:** Expected haplotype diversity
    * May be denoted by other names (e.g. H, h, or Hs).
    * This will often be reported in a table or sometimes an appendix/supplemental material for the paper.
17. **Hese:** Standard error of the He measurement
    * You can convert standard deviation to standard error by dividing by sqrt(n).
18. **Pi:** The nucleotide diversity
19. **Pise:** Standard error of the pi measurement
    * You can convert standard deviation to standard error by dividing by sqrt(n).

**Each line MUST have Spp, Source, n, He OR pi, and LatDeg!!** The other cells are important to have, but are not reasons for skipping the line of data.

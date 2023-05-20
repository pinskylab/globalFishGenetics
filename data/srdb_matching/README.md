CSV files with information on matching stock IDs (RAM Legacy fishery stocks) for some of the microsatellite & mtDNA data.

Both files contain columns with following information:
1. **fbsci:** Species scientific name (of matching stock from RAM database)
2. **Country:** Country of the matching stock from RAM database
3. **Site:** Name of the location of the matching stock from RAM database
4. **stockid:** Name of the stock from the RAM database
1. **Spp:** Species scientific name
2. **CommonName:** Species common name
3. **Source:** Study the data comes from
13. **CollectionYear:** Year in which the samples were taken
15. **MarkerName:** Name of the marker, as listed in the paper
17. **n:** Number of individuals sampled
19. **He:** Expected heterozygosity
20. **Hese:** Standard error of the He measurement
21. **file:** name of the data file that contains data from the matching population
22. **lat:** latitude 
23. **lon:** longitude

`msat_to_match.csv` also contains *PrimerNote*, *NumMarkers*, *CrossSpp*, and *Repeat* columns.

`mtdna_to_match.csv` also contains *bp*, *Pi*, and *Pise* columns. It does not contain a *files* column, because it only matches to the `Fishery lat mtDNA Complete Database.csv` file.

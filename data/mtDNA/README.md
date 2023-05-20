CSV files with raw mitochondrial DNA data (eg, the spreadsheets where data were initially recorded). These files are read into assemble_data_mtdna.R to create one unified, merged dataset with all mtDNA data.

All files contain the following information in each column (for more information, look at `data/data_entry_instructions/mtDNA_data_entry_instructions.md`):

1. **Spp:** Species scientific name
2. **CommonName:** Species common name
3. **Source:** Study the data comes from
4. **Country:** Country where the sample was taken
5. **Site:** Site where sample was taken, as named by author
6. **Lat_Deg:** Degrees latitude (may be decimal degrees)
7. **Lat_Min:** Minutes latitude
8. **Lat_Sec:** Seconds latitude
9. **Lon_Deg:** Degrees longitude (may be decimal degrees)
10. **Lon_Min:** Minutes longitude
11. **Lon_Sec:** Seconds longitude
12. **CollectionYear:** Year in which the samples were taken
13. **MarkerName:** Name of the mtDNA locus, as listed in the paper
14. **n:** Number of individuals sampled
15. **bp:** Length of the mtDNA locus, in basepairs
16. **He:** Expected haplotype diversity
17. **Hese:** Standard error of the He measurement
18. **Pi:** Nucleotide diversity
19. **Pise:** Standard error of the pi measurement
20. **Notes:** Notes made by the recorder during the data collection process

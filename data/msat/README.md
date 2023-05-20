CSV files with raw microsatellite data (eg, the spreadsheets where data were initially recorded). These files are read into assemble_data_msat.R to create one unified, merged dataset with all microsatellite data.
*  `ppdat_2016-03-04wLL.csv` contains data that were originally recorded for use in [Pinsky & Palumbi (2014)](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.12509). The rest of the files contain data initially recorded for this analysis.

All files contain the following information in each column (for more information, look at `data/data_entry_instructions/msat_data_entry_instructions.md`):
1. **Spp:** Species scientific name
2. **CommonName:** Species common name
3. **Source:** Study the data comes from
4. **PrimerNote:** Indicates whether the study is a primer note or not (1 for yes, 0 for no)
5. **Country:** Country where the sample was taken
6. **Site:** Site where sample was taken, as named by author
7. **Lat_Deg:** Degrees latitude (may be decimal degrees)
8. **Lat_Min:** Minutes latitude
9. **Lat_Sec:** Seconds latitude
10. **Lon_Deg:** Degrees longitude (may be decimal degrees)
11. **Lon_Min:** Minutes longitude
12. **Lon_Sec:** Seconds longitude
13. **CollectionYear:** Year in which the samples were taken
14. **NumMarkers:** Number of microsatellite markers
15. **MarkerName:** Name of the microsatellite marker, as listed in the paper
16. **CrossSpp:** Indicates whether the microsatellite was originally developed in a different species (1 if yes, 0 if no)
17. **n:** Number of individuals sampled
18. **Repeat:** Length of the microsatellite repeat
19. **He:** Expected heterozygosity
20. **Hese:** Standard error of the He measurement
21. **Columns 21-89:** Extra columnsto hold allele frequencies (for calculating He), if He is not reported

*Note:* For the `msat_2011-2020_data*csvs`, Column 21 is instead a "Notes" column, that contains any notes made by the recorder during the data collection process.

Contains files read into R scripts for analyses.
* *Note: The main exception to this are the csvs that contain the raw data read into either of the assemble R scripts. These can be found in the `data` repository.*

Data for R scripts:

*  **msat_assembled.csv:** Compiled microsatellite database. Read into `pull_range_data.R`.
*  **msat_env.csv:** Environmental data (SST, chlorophyll, & dissolved oxygen) that corresponds with populations in the nuclear (msat) dataset. Read into `bootstrap_he.R`, `maps.R`, `msat_he_models.R`, and `msat_he_predict.R`.
*  **msatloci_assembled.csv:** Compiled microsatellite database, where averaged He values have been split out into separate rows by locus. Read into `ID_shared_species.R`, `bootstrap_he.R`, `maps.R`, `msat_he_models.R`, `msat_he_predict.R`, and `pull_env_data.R`.
*  **mtdna_assembled.csv:** Compiled mtDNA database. Read into `ID_shared_species.R`, `bootstrap_hd.R`, `bootstrap_pi.R`, `maps.R`, `mtdna_hd_models.R`, `mtdna_hd_predict.R`, `mtdna_pi_models.R`, `mtdna_pi_predict.R`, `pull_env_data.R`, and `pull_range_data.R`.
*  **mtdna_env.csv:** Environmental data (SST, chlorophyll, & dissolved oxygen) that corresponds with populations in the mtDNA dataset. Read into `bootstrap_hd.R`, `bootstrap_pi.R`, `mtdna_hd_models.R`, `mtdna_hd_predict.R`, `mtdna_pi_models.R`, and `mtdna_pi_predict.R`.
*  **sharedstudies.csv:** Combined nuclear (msat) He, mtdna Hd & pi diversity observations for populations where both nuclear & mitochondrial genetic diversity were recorded. Read into (and created by) `ID_shared_species.R`.
*  **spp_combined_info.csv:** Taxonomic information (order, family, genus), IUCN status, and range extent information for all species in either (or both) the nuclear (msat) and mitochondrial datasets. Read into `bootstrap_hd.R`, `bootstrap_he.R`, `bootstrap_pi.R`, `msat_he_models.R`, `msat_he_predicts.R`, `mtdna_hd_models.R`, `mtdna_hd_predict.R`, `mtdna_pi_models.R`, and `mtdna_pi_predict.R`.

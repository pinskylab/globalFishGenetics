**Scripts:**

* **assemble_data_msat.R:** Cleans and assembles all csv files in the `data/msat` directory into one cohesive dataframe.
* **assemble_data_mtDNA.R:** Cleans and assembles all csv files in the `data/mtDNA` directory into one cohesive dataframe.
* **bootstrap_hd.R:** Bootstraps mtDNA Hd models to calculate 95% confidence intervals for model coefficients.
* **bootstrap_he.R:** Bootstraps nuclear msat He models to calculate 95% confidence intervals for model coefficients.
* **bootstrap_pi.R:** Bootstraps mtDNA pi models to calculate 95% confidence intervals for model coefficients.
* **coefficient_bootstrap_cis.R:** Reads in output from the `bootstrap*.R` files  to create supplemental figures summarizing results.
* **ID_shared_species.R:** Identifies populations/species where both mtDNA and nuclear (microsatellite) genetic diversity was measured.
* **maps.R:** Creates maps of the global distribution of all genetic diversity estimates, as well as mean chlorophyll and SST.
* **msat_he_models.R:** Code for all nuclear (microsatellite) He models.
* **msat_he_predict.R:** Plots marginal effects (relationship between msat He and a given predictor variable) for all msat He models.
* **mtdna_hd_models.R:** Code for all mitochondrial Hd models.
* **mtdna_hd_predict.R:** Plots marginal effects (relationship between mtDNA Hd and a given predictor variable) for all mtDNA Hd models.
* **mtdna_pi_models.R:** Code for all mitochondrial pi models.
* **mtdna_pi_predict.R:** Plots marginal effects (relationship between mtDNA pi and a given predictor variable) for all mtDNA pi models.
* **pull_env_data.R:** Pulls corresponding environmental data from Bio-ORACLE for populations in both nuclear and mtDNA datasets.
* **pull_range_data.R:** Pulls range extent data from FishBase (e.g. AquaMaps) for all species in either (or both) the nuclear and mtDNA datasets.

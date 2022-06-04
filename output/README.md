# Output Readme

This directory has files that have been compiled by scripts.

- msat_assembled.csv: Compiled microsatellite database
  - spp: scientific name
  - CommonName: common name for the species
  - Source: citation for the source paper
  - PrimerNote: 1 if the paper is a Primer Note paper published to report on microsatellite primers
  - Country: country
  - Site: site name, as reported by the authors
  - lat: latitude
  - lon: longitude
  - stockid: RAM Legacy stockid, if known
  - CollectionYear: year samples were collected
  - NumMarkers: number of markers averaged together on this line of data
  - MarkerName: name of the marker, as reported by the authors
  - CrossSpp: 1 if the primers were developed in a different species, 0 if not
  - n: number of individuals
  - Repeat: length of the microsatellite repeat (dinucleotide, tri-, or tetra-)
  - He: heterozygosity
  - Hese: standard error of He, if known
  - file: file in data/ that this line of data came from
- mtdna_assembled.csv: Compiled mtDNA database
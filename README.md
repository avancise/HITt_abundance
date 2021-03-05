# HITt_abundance
These scripts support the POPAN abundance estimation described in Van Cise et al. 2021 (Endangered Species Research, DOI: https://doi.org/10.3354/esr01117). Descriptions of each script can be found below. 

## Ttru_abundance.Rmd
This Rmarkdown file compiles all analysis scripts in the project. It also includes manuscript text and code for figures and tables included in the manuscript.

## RMark_abundance.R
This script initiates data input for the program MARK and uses RMark to run POPAN model and generate parameter estimates in MARK.

## RMark_abundance_datainput.R
This script reads in photo identification data from the Cascadia Research Collective database, cleans and filters the data, and stratifies indviduals by sighting location.

## test.rand.sample.R
This script tests whether sampling effort in each stock differs significantly from random. Analysis modified from original code provided by K. Alexandra Curtis (alex.curtis@noaa.gov), shared March 6, 2020.

## distance.permutation.test.greatcircle.R
This script tests whether individual inter-annual movements are significantly smaller than stock inter-annual movements. Analysis modified from original code provided by K. Alexandra Curtis (alex.curtis@noaa.gov), shared March 6, 2020.

## Effort testing.R
This script generates summaries of effort, e.g. number of observations, number of individuals observed, median observations per individual.

## Prop_distinctive.R
This script estimates the proportion of highly distinctive animals in the population, used for downstream correction of abundance estimates.

## subarea resights.R
This script estimates the proportion of times an individual is resighted within the same subarea where as its original sighting location.

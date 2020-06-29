# HITt_abundance
These scripts support the POPAN abundance estimation described in Van Cise et al. (in prep). Descriptions of each script can be found below. 

## RMark_abundance_datainput.R
This script reads in photo identification data from the Cascadia Research Collective database, cleans and filters the data, and stratifies indviduals by sighting location.

## RMark_abundance.R
This script initiates data input for the program MARK and uses RMark to run POPAN estimates in MARK.

## test.rand.sample.R
This script tests whether sampling effort in each stock differs significantly from random.

## distance.permutation.test.greatcircle.R
This script tests whether individual inter-annual movements are significantly smaller than stock inter-annual movements.

## Effort testing.R
This script generates summaries of effort, e.g. number of observations, number of individuals observed, median observations per individual.

## Prop_distinctive.R
This script estimates the proportion of highly distinctive animals in the population, used for downstream correction of abundance estimates.


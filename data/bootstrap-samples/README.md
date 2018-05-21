# Bootstrap Standard Errors

Need to resample the microdata with replacement (same sample size), then run it through the smoothing and estimation scripts.

## Resampling scheme

Resampling should proxy what a new sample would look like.
As the ACS is sampled at the household level, need to resample from a list of unique household IDs.
The sampling is likely stratified, so can resample within relevant groupings to mimic repeated sampling:
* By region: to get enough observations in each area; so block resample within each MSA
* By year: same sample size targeted each year

## Computational method

1. Create minimal table of microdata as base table for resampling
	* Top 20 MSAs
	* Only the variables used in estimation
2. Generate bootstrap resamples in R for each MSA and YEAR:
	* Get list of unique household IDs
	* Resample within YEAR, MSA
	* Inner join back to base data on YEAR, MSA, SERIAL
3. Verify: the fraction of unique SERIAL in resample should be approximately `1 - exp(-1)` 
	* Success!

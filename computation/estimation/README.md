# Estimation

Use `csv2sqlite.py` to save the raw csv.gz data file into an sqlite database.
Then use SQL queries to get aggregated values to avoid loading the entire dataset into memory.
Queries apply categorizations (race, edu) on-the-fly, so no need to pre-clean the data.

## Estimation of rate parameters

* Using `R` and `data.table`
* Need marriage and divorce rates for each couple-type (globally)
	* Marriage rate (directly observable): SQL queries for flows and stocks to compute rates
	* Divorce rate (infer from non-divorce rate and death rate)
* Weighted OLS (by stocks of couples)

## Estimation of model objects

* Using `R` and `data.table`: query aggregated population counts in each desired metro
	* Smooth raw population counts appropriately with kernel density estimator (Nadaraya-Watson with Sheather-Jones bandwidth)
	* Save smoothed data for loading into `julia`
* Using `julia`: convert DataFrames to multidimensional arrays (per metro)
	* populations: single men, single women, couples
	* compute estimates: alpha, surplus, production function
	* average the production function estimates to get a global result

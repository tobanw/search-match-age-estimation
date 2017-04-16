# Estimation

## Prepare ACS data

Use `csv2sqlite.py` to save the raw csv.gz data file into an sqlite database.
Then use SQL queries to get aggregated values to avoid loading the entire dataset into memory.
Queries apply categorizations (race, edu) on-the-fly, so no need to pre-clean the data.

## Estimation of rate parameters by OLS

* `estimate-rates-full.R` and `estimate-rates.R`
	* very poor accuracy due to noisy inference on divorce flows
* Need marriage and divorce rates for each couple-type (globally)
	* Marriage rate (directly observable): SQL queries for flows and stocks to compute rates
	* Divorce rate (infer from non-divorce rate and death rate)
* Weighted OLS (by stocks of couples)

## Estimation of model objects

* `smooth-pops.R`: query aggregated population counts in each desired metro
	* Smooth raw population counts: non-parametric regression (local-linear)
	* Save smoothed data to csv for loading into `julia`
* `mort-rates.R`: interpolates and saves death rates
* `prepare-pops.jl`: convert DataFrames to multidimensional arrays (per metro)
	* loads in populations and death rates
* `main-estim.jl`: calls `estim-functions.jl` to
	* compute estimates: alpha, surplus, production function
	* average the production function estimates to get a global result

## Estimation of rate parameters by SMM

* moments: marriage and divorce flows by individual types (not couple pairs)
* `smooth-flows.R`: queries and smooths data moments
	* `prepare-pops.jl`: loads up data moments to arrays per MSA
* `estim-functions.jl`: provides a function to compute the model moments given alpha
* `main-estim.jl`: performs the MD estimation
	* loads data moments, computes model moments
	* runs an optimizer to minimize weighted (by population size) sum of squared differences

## 10 largest metro areas

* 35620: New York-Newark-Jersey City, NY-NJ-PA
* 31080: Los Angeles-Long Beach-Anaheim, CA
* 16980: Chicago-Naperville-Elgin, IL-IN-WI
* 19100: Dallas-Fort Worth-Arlington, TX
* 37980: Philadelphia-Camden-Wilmington, PA-NJ-DE-MD
* 26420: Houston-The Woodlands-Sugar Land, TX
* 47900: Washington-Arlington-Alexandria, DC-VA-MD-WV
* 33100: Miami-Fort Lauderdale-West Palm Beach, FL
* 12060: Atlanta-Sandy Springs-Roswell, GA
* 14460: Boston-Cambridge-Newton, MA-NH

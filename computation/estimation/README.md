# Estimation

## Prepare ACS data

Use `csv2sqlite.py` to save the raw csv.gz data file into an sqlite database.
* `acs`: table from ACS microdata
* `migtopuma`: table to convert migration state/puma to puma
* `pumatomet`: table to convert puma to metro
Then use SQL queries to get aggregated values to avoid loading the entire dataset into memory.
Queries apply categorizations (race, edu) on-the-fly, so no need to pre-clean the data.

## Estimation of model objects

* set up data:
	* `smooth-pops.R`: query and smooth aggregated population counts in each desired metro
	* `smooth-flows.R`: query and smooth data moments
	* Smoothing by non-parametric regression (local-linear)
	* Saves smoothed data to csv for loading into `julia`
	* `mort-rates.R`: interpolates and saves death rates
	* `prepare-pops.jl`: load above csv files, then convert DataFrames to multidimensional arrays (per metro) and saves as JLD files
* `main-estim.jl`: runs the show, but need to set options first
	* loads populations from saved JLD files, or calls `prepare-pops.jl` to generate them anew
	* estimate arrival rates and then non-parametric objects using `estim-functions.jl` and `compute-npobj.jl`
	* can also do a parameter grid search or monte carlo estimation

## Largest metro areas by adult population (millions)

1. 35620: 14.5m - New York-Newark-Jersey City, NY-NJ-PA
2. 31080: 9.4m - Los Angeles-Long Beach-Anaheim, CA
3. 16980: 6.8m - Chicago-Naperville-Elgin, IL-IN-WI
4. 19100: 4.6m - Dallas-Fort Worth-Arlington, TX
5. 37980: 4.4m - Philadelphia-Camden-Wilmington, PA-NJ-DE-MD
6. 26420: 4.2m - Houston-The Woodlands-Sugar Land, TX
7. 47900: 4.1m - Washington-Arlington-Alexandria, DC-VA-MD-WV
8. 33100: 4.1m - Miami-Fort Lauderdale-West Palm Beach, FL
9. 12060: 3.8m - Atlanta-Sandy Springs-Roswell, GA
10. 14460: 3.5m - Boston-Cambridge-Newton, MA-NH
11. 41860: 3.3m - San Francisco-Oakland-Hayward, CA
12. 19820: 3.1m - Detroit-Warren-Dearborn, MI
13. 38060: 3.1m - Phoenix-Mesa-Scottsdale, AZ
14. 40140: 3.0m - Riverside-San Bernardino-Ontario, CA
15. 42660: 2.6m - Seattle-Tacoma-Bellevue, WA
16. 33460: 2.4m - Minneapolis-St. Paul-Bloomington, MN-WI
17. 41740: 2.3m - San Diego-Carlsbad, CA
18. 45300: 2.1m - Tampa-St. Petersburg-Clearwater, FL
19. 41180: 2.0m - St. Louis, MO-IL
20. 12580: 2.0m - Baltimore-Columbia-Towson, MD

30th metro is 1.4m, 40th is 0.9m.

## (Deprecated) Estimation of rate parameters by OLS

* `estimate-rates-full.R` and `estimate-rates.R`
	* very poor accuracy due to noisy inference on divorce flows
* Need marriage and divorce rates for each couple-type (globally)
	* Marriage rate (directly observable): SQL queries for flows and stocks to compute rates
	* Divorce rate (infer from non-divorce rate and death rate)
* Weighted OLS (by stocks of couples)

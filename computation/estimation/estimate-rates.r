library(DBI) # RSQLite database functions
library(data.table)
library(zoo) # for interpolating mortality rates with na.approx

# Order of traits: husband then wife; age, edu, race
# Notation: husband gets _SP suffix, wife is default
# Min age is 25 because of endogeneity of college

# load up CDC mortality table: by sex, age, and minority
mort.dt <- fread('mort_08-15.csv')

# connect to sqlite database
# table name: acs
db = dbConnect(RSQLite::SQLite(), 'acs_08-14.db')


### Categorization ###

case_minority = ' case when "RACESING" in (1,4) and "HISPAN" = 0 then 0 else 1 end '
case_minority_sp = ' case when "RACESING_SP" in (1,4) and "HISPAN_SP" = 0 then 0 else 1 end '
case_college = ' case when "EDUC" >= 10 then 1 else 0 end '
case_college_sp = ' case when "EDUC_SP" >= 10 then 1 else 0 end '


### Queries ###

# flows: counts of new marriages by type pair
qry_marr_flow = paste0('select
		"AGE_SP" as AGE_M, ', case_college_sp, ' as COLLEGE_M, ', case_minority_sp, ' as MINORITY_M,
		"AGE" as AGE_F, ', case_college, ' as COLLEGE_F, ', case_minority, ' as MINORITY_F,
		sum("HHWT") as MARFLOW
	from acs
	where "MARRINYR" = 2 and "MARST" <= 2
		and "SEX" = 2
		and AGE_F >= 25 and AGE_M >= 25
	group by AGE_M, COLLEGE_M, MINORITY_M, AGE_F, COLLEGE_F, MINORITY_F')

# stocks (lagged): counts of all marriages by type pair (by year for separating cohorts)
qry_marr_stock = paste0('select "YEAR" + 1 as YEAR,
		"AGE_SP" + 1 as AGE_M, ', case_college_sp, ' as COLLEGE_M, ', case_minority_sp, ' as MINORITY_M,
		"AGE" + 1 as AGE_F, ', case_college, ' as COLLEGE_F, ', case_minority, ' as MINORITY_F,
		sum("HHWT") as MARSTOCK
	from acs
	where "MARST" <= 2
		and "SEX" = 2
		and AGE_F >= 25 and AGE_M >= 25
	group by YEAR, AGE_M, COLLEGE_M, MINORITY_M, AGE_F, COLLEGE_F, MINORITY_F')

# stocks of surviving marriages (by year for separating cohorts)
qry_m1_stock = paste0('select "YEAR" as YEAR,
		"AGE_SP" as AGE_M, ', case_college_sp, ' as COLLEGE_M, ', case_minority_sp, ' as MINORITY_M,
		"AGE" as AGE_F, ', case_college, ' as COLLEGE_F, ', case_minority, ' as MINORITY_F,
		sum("HHWT") as M1STOCK
	from acs
	where "MARST" <= 2 and "MARRINYR" != 2
		and "SEX" = 2
		and AGE_F >= 25 and AGE_M >= 25
	group by YEAR, AGE_M, COLLEGE_M, MINORITY_M, AGE_F, COLLEGE_F, MINORITY_F')

# counts of singles eligible to marry within the past year
qry_men_elig = paste0('select
		"AGE" as AGE_M, ', case_college, ' as COLLEGE_M, ', case_minority, ' as MINORITY_M,
		sum("PERWT") as SNG_M
	from acs
	where "SEX" = 1 and AGE_M >= 25
		and ("MARRINYR" = 2 or "MARST" >= 3)
		and "DIVINYR" != 2
	group by AGE_M, COLLEGE_M, MINORITY_M')
qry_wom_elig = paste0('select
		"AGE" as AGE_F, ', case_college, ' as COLLEGE_F, ', case_minority, ' as MINORITY_F,
		sum("PERWT") as SNG_F
	from acs
	where "SEX" = 2 and AGE_F >= 25
		and ("MARRINYR" = 2 or "MARST" >= 3)
		and "DIVINYR" != 2
	group by AGE_F, COLLEGE_F, MINORITY_F')

# queries return dataframes, convert to data.table
master.dt = data.table(dbGetQuery(db, qry_marr_flow))
mar.stock = data.table(dbGetQuery(db, qry_marr_stock))
m1.stock = data.table(dbGetQuery(db, qry_m1_stock))

men.elig = data.table(dbGetQuery(db, qry_men_elig))
wom.elig = data.table(dbGetQuery(db, qry_wom_elig))


### Divorce rates ###

div.dt <- merge(mar.stock, m1.stock,
				by=c("YEAR", "AGE_M", "COLLEGE_M", "MINORITY_M", "AGE_F", "COLLEGE_F", "MINORITY_F"),
				all=TRUE) # merge the MARSTOCK columns

# master.dt is aggregate, div.dt is by cohort-year
# need to compute DIVRATE by cohort and then average and merge into master.dt

## Death rates

# prep mortality table for merging
mort.dt[, DTHRT := DR100 / 100000] # convert individual probabilities
# lag by one year (5+1) so that death rate is for past year
mort.dt[, AGE := AGE_L + 6] # midpoints for interpolation: 31,41,51,61,71
mort.dt[AGE_H==24, AGE := 25] # boundaries for interpolation
mort.dt[AGE_H==84, AGE := 79] # boundaries for interpolation

# merge sparse mortality rates into full age table
men.elig <- merge(men.elig, mort.dt[SEX==1, .(AGE_M = AGE, MINORITY_M = MINORITY, DTHRT_M = DTHRT)],
				   by=c("AGE_M", "MINORITY_M"), all=TRUE) # male 
wom.elig <- merge(wom.elig, mort.dt[SEX==2, .(AGE_F = AGE, MINORITY_F = MINORITY, DTHRT_F = DTHRT)],
				   by=c("AGE_F", "MINORITY_F"), all=TRUE) # female 

# interpolate death rates over ages within minority group
men.elig[, DTHRT_M := na.approx(DTHRT_M, AGE_M), MINORITY_M] # male
wom.elig[, DTHRT_F := na.approx(DTHRT_F, AGE_F), MINORITY_F] # female

# merge in death rates (and stocks of singles)
div.dt <- merge(div.dt, men.elig, by=c("AGE_M", "COLLEGE_M", "MINORITY_M"))
div.dt <- merge(div.dt, wom.elig, by=c("AGE_F", "COLLEGE_F", "MINORITY_F"))

# separation rate due to deaths: assuming independence, which overestimates because of
#  correlation between couple
div.dt[, CPLDR := DTHRT_F + DTHRT_M - DTHRT_F * DTHRT_M]

# divorce flows: total attrition less death flows
#div.dt[, DIVFLOW := MARSTOCK * (1 - CPLDR) - M1STOCK]

# divorce rates calculated per separate cohort
div.dt[, DIVRTC := 1 - CPLDR - M1STOCK / MARSTOCK] # 1 - death rate - non-div rate

# merge in the average divorce rate and total marriage stock (by cohort)
#   as well as stocks of singles
master.dt <- merge(master.dt,
				   div.dt[, .(DIVRATE = mean(DIVRTC), WEIGHT = sum(MARSTOCK),
							  SNG_M = SNG_M[1], SNG_F = SNG_F[1]),
						  .(AGE_M, COLLEGE_M, MINORITY_M, AGE_F, COLLEGE_F, MINORITY_F)],
				   by=c("AGE_M", "COLLEGE_M", "MINORITY_M", "AGE_F", "COLLEGE_F", "MINORITY_F"))


### Marriage rates ###

master.dt[, MARRATE := MARFLOW / SNG_M / SNG_F] # sequential division to prevent integer overflow error


### Weighted OLS ###

# drop NA rates, which may really bias the estimates
reg.dt <- master.dt[!is.na(MARRATE) & !is.na(DIVRATE),
					.(Y=1, MARRATE, DIVRATE, WEIGHT)] # add in Y for regression

# weighted regression without an intercept
ols.reg <- lm(Y ~ 0 + MARRATE + DIVRATE, data=reg.dt, weights=WEIGHT)

#summary(ols.reg) to view results

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

case_minority = 'case
when "RACESING" in (1,4) and "HISPAN" = 0 then 0
else 1
end as MINORITY'
case_minority_sp = 'case
when "RACESING_SP" in (1,4) and "HISPAN_SP" = 0 then 0
else 1
end as MINORITY_SP'
case_college = 'case
when "EDUC" >= 10 then 1
else 0
end as COLLEGE'
case_college_sp = 'case
when "EDUC_SP" >= 10 then 1
else 0
end as COLLEGE_SP'


### Marriage rates ###

# counts of new marriages by type pair
qry_marr_flow = paste0('select
		"AGE_SP", ', case_college_sp, ', ', case_minority_sp,
		', "AGE", ', case_college, ', ', case_minority,
		', sum("HHWT") as MARFLOW
	from acs
	where ("MARRINYR" = 2 and "MARST" <= 2)
		and "SEX" = 2
		and "AGE" >= 25
	group by "AGE_SP", COLLEGE_SP, MINORITY_SP, "AGE", COLLEGE, MINORITY')

# counts of all marriages by type pair
qry_marr_stock = paste0('select
		"AGE_SP", ', case_college_sp, ', ', case_minority_sp,
		', "AGE", ', case_college, ', ', case_minority,
		', sum("HHWT") as MARSTOCK
	from acs
	where "MARST" <= 2
		and "SEX" = 2
		and "AGE" >= 25
	group by "AGE_SP", COLLEGE_SP, MINORITY_SP, "AGE", COLLEGE, MINORITY')

# counts of singles eligible to marry
qry_men_elig = paste0('select
		"AGE", ', case_college, ', ', case_minority,
		', sum("PERWT") as SINGMEN
	from acs
	where "SEX" = 1
		and ("MARRINYR" = 2 or "MARST" >= 3)
		and "AGE" >= 25
	group by "AGE", COLLEGE, MINORITY')
qry_wom_elig = paste0('select
		"AGE", ', case_college, ', ', case_minority,
		', sum("PERWT") as SINGWOM
	from acs
	where "SEX" = 2
		and ("MARRINYR" = 2 or "MARST" >= 3)
		and "AGE" >= 18
	group by "AGE", COLLEGE, MINORITY')

# population counts of men and women for weighting cells in regression
qry_men_pop = paste0('select
		"AGE", ', case_college, ', ', case_minority,
		', sum("PERWT") as MENPOP
	from acs
	where "SEX" = 1 and "AGE" >= 25
	group by "AGE", COLLEGE, MINORITY')
qry_wom_pop = paste0('select
		"AGE", ', case_college, ', ', case_minority,
		', sum("PERWT") as WOMPOP
	from acs
	where "SEX" = 2 and "AGE" >= 25
	group by "AGE", COLLEGE, MINORITY')

# queries return dataframes, convert to data.table
master.dt = data.table(dbGetQuery(db, qry_marr_flow))
mar.stock = data.table(dbGetQuery(db, qry_marr_stock))

men.elig = data.table(dbGetQuery(db, qry_men_elig))
wom.elig = data.table(dbGetQuery(db, qry_wom_elig))

men.pop = data.table(dbGetQuery(db, qry_men_pop))
wom.pop = data.table(dbGetQuery(db, qry_wom_pop))

# set keys for easy merging
setkeyv(master.dt, colnames(master.dt)[1:6])
setkeyv(mar.stock, colnames(mar.stock)[1:6])

master.dt <- merge(master.dt, mar.stock, all=TRUE) # merge in the MARSTOCK column

# merge in stocks of eligible singles (for computing rates from flows)
master.dt <- merge(master.dt, men.elig,
				   by.x=c("AGE_SP", "COLLEGE_SP", "MINORITY_SP"),
				   by.y=c("AGE", "COLLEGE", "MINORITY")) # all=TRUE)
master.dt <- merge(master.dt, wom.elig,
				   by=c("AGE", "COLLEGE", "MINORITY")) # all=TRUE)

# compute marriage rates
master.dt[, MARRATE := MARFLOW / SINGMEN / SINGWOM] # sequential division to prevent integer overflow error


### Divorce rates ###

master.dt[, MARSTOCK1 := MARSTOCK - MARFLOW] # marriages at least one year old

## Death rates

# prep mortality table for merging
mort.dt[, DTHRT := DR100 / 100000] # convert individual probabilities
mort.dt[, AGE := AGE_L + 5] # midpoints for interpolation: 30,40,50,60,70
mort.dt[AGE_H==24, AGE := 25] # boundaries for interpolation
mort.dt[AGE_H==84, AGE := 79] # boundaries for interpolation

# merge sparse mortality rates into full age table
wom.pop <- merge(wom.pop, mort.dt[SEX==2, .(AGE, MINORITY, FDTHRT = DTHRT)],
				   by=c("AGE", "MINORITY"), all=TRUE) # female 
men.pop <- merge(men.pop, mort.dt[SEX==1, .(AGE, MINORITY, MDTHRT = DTHRT)],
				   by=c("AGE", "MINORITY"), all=TRUE) # male 

# interpolate death rates over ages within minority group
men.pop[, MDTHRT := na.approx(MDTHRT, AGE), MINORITY] # male
wom.pop[, FDTHRT := na.approx(FDTHRT, AGE), MINORITY] # female

# merge in population stocks and death rates
master.dt <- merge(master.dt, men.pop,
				   by.x=c("AGE_SP", "COLLEGE_SP", "MINORITY_SP"),
				   by.y=c("AGE", "COLLEGE", "MINORITY"))
master.dt <- merge(master.dt, wom.pop,
				   by=c("AGE", "COLLEGE", "MINORITY"))

# separation rate due to deaths: assuming independence, which overestimates because of
#  correlation between couple
master.dt[, CPLDR := FDTHRT + MDTHRT - FDTHRT * MDTHRT]

#FIXME: shift within couple types for BOTH ages!
# death flows: product of younger cohort's stock and rate 
master.dt[, DTHFLOW := shift(CPLDR) * shift(MARSTOCK)] # lag: initial year NA

#FIXME: shift within couple types for BOTH ages!
# divorce flows: total attrition less death flows
master.dt[, DIVFLOW := MARSTOCK - shift(MARSTOCK1, type="lead") - DTHFLOW] # lead: final year NA

# compute divorce rates
master.dt[, DIVRATE := DIVFLOW / MARSTOCK]
# FIXME: this is way too big, almost all values are larger than 0.1 in abs!


### Weighted OLS ###

# TODO: weight by MARSTOCK and drop obs with small values that give noisy DIVRATE
# weight by pop size ℓm*ℓf, sample representativeness (as in GJR)
master.dt[, WEIGHT := 1e-9 * MENPOP * WOMPOP] # shrink to prevent integer overflow error

# drop NA rates, which may really bias the estimates
reg.dt <- master.dt[!is.na(MARRATE) & !is.na(DIVRATE),
					.(Y=1, MARRATE, DIVRATE, WEIGHT)] # add in Y for regression

# weighted regression without an intercept
ols.reg <- lm(Y ~ 0 + MARRATE + DIVRATE, data=reg.dt, weights=WEIGHT)

# print results
summary(ols.reg)

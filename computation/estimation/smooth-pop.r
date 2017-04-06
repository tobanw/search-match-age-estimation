# Non-parametric (local-linear) regression to smooth population stocks

# Order of traits: husband then wife; age, edu, race
# Notation: husband gets _SP suffix, wife is default
# Min age is 25 because of endogeneity of college
# ACS data: 2008-2014, ages 18-79 (note: AGE_SP is not limited).

# TODO: try quarterly age resolution for npreg; save smoothed values on chosen (coarser) age grid

library(DBI) # RSQLite database functions
library(data.table)
library(np) # non-parametric regression

# connect to sqlite database
# table name: acs
db = dbConnect(RSQLite::SQLite(), 'data/acs_08-14.db')


### Categorization ###

case_minority = ' case when "RACESING" in (1,4) and "HISPAN" = 0 then 0 else 1 end '
case_minority_sp = ' case when "RACESING_SP" in (1,4) and "HISPAN_SP" = 0 then 0 else 1 end '
case_college = ' case when "EDUC" >= 10 then 1 else 0 end '
case_college_sp = ' case when "EDUC_SP" >= 10 then 1 else 0 end '


### Queries ###

# total population counts
qry_men_tot = paste0('select "MET2013" as MSA,
		"AGE" as AGE_M, ', case_college, ' as COLLEGE_M, ', case_minority, ' as MINORITY_M,
		sum("PERWT") as MEN
	from acs
	where "SEX" = 1 and AGE_M >= 25
	group by MSA, AGE_M, COLLEGE_M, MINORITY_M')
qry_wom_tot = paste0('select "MET2013" as MSA,
		"AGE" as AGE_F, ', case_college, ' as COLLEGE_F, ', case_minority, ' as MINORITY_F,
		sum("PERWT") as WOM
	from acs
	where "SEX" = 2 and AGE_F >= 25
	group by MSA, AGE_F, COLLEGE_F, MINORITY_F')

# singles counts
qry_men_sng = paste0('select "MET2013" as MSA,
		"AGE" as AGE_M, ', case_college, ' as COLLEGE_M, ', case_minority, ' as MINORITY_M,
		sum("PERWT") as SNG_M
	from acs
	where "SEX" = 1 and AGE_M >= 25 and "MARST" >= 3
	group by MSA, AGE_M, COLLEGE_M, MINORITY_M')
qry_wom_sng = paste0('select "MET2013" as MSA,
		"AGE" as AGE_F, ', case_college, ' as COLLEGE_F, ', case_minority, ' as MINORITY_F,
		sum("PERWT") as SNG_F
	from acs
	where "SEX" = 2 and AGE_F >= 25 and "MARST" >= 3
	group by MSA, AGE_F, COLLEGE_F, MINORITY_F')

# married counts
qry_marr = paste0('select "MET2013" as MSA,
		"AGE_SP" as AGE_M, ', case_college_sp, ' as COLLEGE_M, ', case_minority_sp, ' as MINORITY_M,
		"AGE" as AGE_F, ', case_college, ' as COLLEGE_F, ', case_minority, ' as MINORITY_F,
		sum("HHWT") as MARRIAGES
	from acs
	where "MARST" <= 2
		and "SEX" = 2
		and AGE_F >= 25 and AGE_M >= 25
	group by MSA, AGE_M, COLLEGE_M, MINORITY_M, AGE_F, COLLEGE_F, MINORITY_F')

# queries return dataframes, convert to data.table
men.tot = data.table(dbGetQuery(db, qry_men_tot))
wom.tot = data.table(dbGetQuery(db, qry_wom_tot))
men.sng = data.table(dbGetQuery(db, qry_men_sng))
wom.sng = data.table(dbGetQuery(db, qry_wom_sng))
marriages = data.table(dbGetQuery(db, qry_marr))

# TODO: master grids to ensure the entire domain is computed for conversion to array
# for each table, merge in the grids and fill NA with zeros

### Smoothing ###

# convert categorical variables to factors for npreg
for (col in c('COLLEGE_M', 'MINORITY_M')) set(men.tot, j = col, value = factor(men.tot[[col]]))
for (col in c('COLLEGE_M', 'MINORITY_M')) set(men.sng, j = col, value = factor(men.sng[[col]]))
for (col in c('COLLEGE_F', 'MINORITY_F')) set(wom.tot, j = col, value = factor(wom.tot[[col]]))
for (col in c('COLLEGE_F', 'MINORITY_F')) set(wom.sng, j = col, value = factor(wom.sng[[col]]))
for (col in c('COLLEGE_M', 'MINORITY_M', 'COLLEGE_F', 'MINORITY_F')) set(marriages, j = col, value = factor(marriages[[col]]))

# batch smoothing by MSA: add MASS column to each data table
men.sng[, MASS := predict(npreg(bws=npregbw(formula = SNG_M ~ AGE_M + COLLEGE_M + MINORITY_M,
									regtype="ll",
									bwmethod="cv.aic",
									data=.SD))
                          newdata=.SD),
				  by = MSA]

men.tot[, MASS := predict(npreg(bws=npregbw(formula = MEN ~ AGE_M + COLLEGE_M + MINORITY_M,
									regtype="ll",
									bwmethod="cv.aic",
									data=.SD))
                          newdata=.SD),
				  by = MSA]

wom.sng[, MASS := predict(npreg(bws=npregbw(formula = SNG_F ~ AGE_F + COLLEGE_F + MINORITY_F,
									regtype="ll",
									bwmethod="cv.aic",
									data=.SD))
                          newdata=.SD),
				  by = MSA]

wom.tot[, MASS := predict(npreg(bws=npregbw(formula = WOM ~ AGE_F + COLLEGE_F + MINORITY_F,
									regtype="ll",
									bwmethod="cv.aic",
									data=.SD))
                          newdata=.SD),
				  by = MSA]

# TODO: batch smoothing for marriages
# extremely slow with product kernel: on the order of 10h per MSA (feasible for ~10 MSAs)
# try sample splitting, 16 categorical combos
# -> formula = MARRIAGES ~ AGE_M + AGE_F,... by = .(MSA, COLLEGE_M, MINORITY_M, COLLEGE_F, MINORITY_F)
# if still slow, just limit to top MSAs


### Save output ###

# TODO: make sure to save the smoothed tables
fwrite(marriages, file="data/marriages.csv")
fwrite(men.tot, file="data/men-total.csv")
fwrite(wom.tot, file="data/wom-total.csv")
fwrite(men.sng, file="data/men-single.csv")
fwrite(wom.sng, file="data/wom-single.csv")


### Plots ###

# men.msa: table for men in one msa, with a col of smoothed pop
# to plot masses for each (college, minority) combo:
#ggplot(men.msa, aes(x=AGE_M, y=SMTH)) + geom_line() + facet_grid(COLLEGE_M ~ MINORITY_M)

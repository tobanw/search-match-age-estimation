# Non-parametric (local-linear) regression to smooth population stocks

# Order of traits: husband then wife; age, edu, race
# Notation: husband gets _SP suffix, wife is default
# Min age is 25 because of endogeneity of college
# ACS data: 2008-2014, ages 18-79 (note: AGE_SP is not limited).

# TODO: try quarterly age resolution for npreg; save smoothed values on chosen (coarser) age grid
# TODO: merge death rates

library(DBI) # RSQLite database functions
library(data.table)
library(np) # non-parametric regression

# connect to sqlite database
# table name: acs
db <- dbConnect(RSQLite::SQLite(), 'data/acs_08-14.db')


### Categorization ###

case_minority <- ' case when "RACESING" in (1,4) and "HISPAN" = 0 then 0 else 1 end '
case_minority_sp <- ' case when "RACESING_SP" in (1,4) and "HISPAN_SP" = 0 then 0 else 1 end '
case_college <- ' case when "EDUC" >= 10 then 1 else 0 end '
case_college_sp <- ' case when "EDUC_SP" >= 10 then 1 else 0 end '

# top 10 largest MSAs
top.msa <- c(35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460)
#top_msa <- ' MSA in (35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460) '
top_msa <- ' MSA == 35620 ' # TEST

# terminal age
max.age <- 65

# complete grids to merge in
men.grid <- CJ(AGE_M = 25:79, COLLEGE_M = 0:1, MINORITY_M = 0:1)
wom.grid <- CJ(AGE_F = 25:79, COLLEGE_F = 0:1, MINORITY_F = 0:1)
mar.grid <- CJ(AGE_M = 25:79, COLLEGE_M = 0:1, MINORITY_M = 0:1, AGE_F = 25:79, COLLEGE_F = 0:1, MINORITY_F = 0:1)


### Queries ###

# total population counts
qry_men_tot <- paste0('select "MET2013" as MSA,
		"AGE" as AGE_M, ', case_college, ' as COLLEGE_M, ', case_minority, ' as MINORITY_M,
		sum("PERWT") as MEN
	from acs
	where "SEX" = 1 and AGE_M >= 25 and ', top_msa,
	'group by MSA, AGE_M, COLLEGE_M, MINORITY_M')
qry_wom_tot <- paste0('select "MET2013" as MSA,
		"AGE" as AGE_F, ', case_college, ' as COLLEGE_F, ', case_minority, ' as MINORITY_F,
		sum("PERWT") as WOM
	from acs
	where "SEX" = 2 and AGE_F >= 25 and ', top_msa,
	'group by MSA, AGE_F, COLLEGE_F, MINORITY_F')

# singles counts
qry_men_sng <- paste0('select "MET2013" as MSA,
		"AGE" as AGE_M, ', case_college, ' as COLLEGE_M, ', case_minority, ' as MINORITY_M,
		sum("PERWT") as SNG_M
	from acs
	where "SEX" = 1 and AGE_M >= 25 and "MARST" >= 3 and ', top_msa,
	'group by MSA, AGE_M, COLLEGE_M, MINORITY_M')
qry_wom_sng <- paste0('select "MET2013" as MSA,
		"AGE" as AGE_F, ', case_college, ' as COLLEGE_F, ', case_minority, ' as MINORITY_F,
		sum("PERWT") as SNG_F
	from acs
	where "SEX" = 2 and AGE_F >= 25 and "MARST" >= 3 and ', top_msa,
	'group by MSA, AGE_F, COLLEGE_F, MINORITY_F')

# married counts
qry_marr <- paste0('select "MET2013" as MSA,
		"AGE_SP" as AGE_M, ', case_college_sp, ' as COLLEGE_M, ', case_minority_sp, ' as MINORITY_M,
		"AGE" as AGE_F, ', case_college, ' as COLLEGE_F, ', case_minority, ' as MINORITY_F,
		sum("HHWT") as MARRIAGES
	from acs
	where "MARST" <= 2
		and "SEX" = 2
		and AGE_F >= 25 and AGE_M >= 25
		and ', top_msa,
	'group by MSA, AGE_M, COLLEGE_M, MINORITY_M, AGE_F, COLLEGE_F, MINORITY_F')

# TODO: merging doesn't add MSA, so all the merged rows have NA
# queries return dataframes, convert to data.table and merge into the complete grid
men.tot <- merge(data.table(dbGetQuery(db, qry_men_tot)), men.grid, all=TRUE, by=c("AGE_M", "COLLEGE_M", "MINORITY_M"))
wom.tot <- merge(data.table(dbGetQuery(db, qry_wom_tot)), wom.grid, all=TRUE, by=c("AGE_F", "COLLEGE_F", "MINORITY_F"))
men.sng <- merge(data.table(dbGetQuery(db, qry_men_sng)), men.grid, all=TRUE, by=c("AGE_M", "COLLEGE_M", "MINORITY_M"))
wom.sng <- merge(data.table(dbGetQuery(db, qry_wom_sng)), wom.grid, all=TRUE, by=c("AGE_F", "COLLEGE_F", "MINORITY_F"))
marriages <- merge(data.table(dbGetQuery(db, qry_marr)), mar.grid, all=TRUE, by=c("AGE_M", "COLLEGE_M", "MINORITY_M", "AGE_F", "COLLEGE_F", "MINORITY_F"))

# drop men older than 79
marriages <- marriages[AGE_M <= 79]

# fill NA with zeros before smoothing
men.tot[is.na(MEN), MEN := 0]
men.sng[is.na(SNG_M), SNG_M := 0]
wom.tot[is.na(WOM), WOM := 0]
wom.sng[is.na(SNG_F), SNG_F := 0]
marriages[is.na(MARRIAGES), MARRIAGES := 0]


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
											data=.SD)),
						  newdata=.SD),
        by = MSA]

men.tot[, MASS := predict(npreg(bws=npregbw(formula = MEN ~ AGE_M + COLLEGE_M + MINORITY_M,
											regtype="ll",
											bwmethod="cv.aic",
											data=.SD)),
						  newdata=.SD),
        by = MSA]

wom.sng[, MASS := predict(npreg(bws=npregbw(formula = SNG_F ~ AGE_F + COLLEGE_F + MINORITY_F,
											regtype="ll",
											bwmethod="cv.aic",
											data=.SD)),
						  newdata=.SD),
        by = MSA]

wom.tot[, MASS := predict(npreg(bws=npregbw(formula = WOM ~ AGE_F + COLLEGE_F + MINORITY_F,
											regtype="ll",
											bwmethod="cv.aic",
											data=.SD)),
						  newdata=.SD),
        by = MSA]

## combine groups with age >= max.age (truncation at T)

# lump together truncated masses to be merged back in
trm.men.sng <- men.sng[AGE_M >= max.age, .(AGE_M = max.age, MASS = sum(MASS)), by = .(MSA, COLLEGE_M, MINORITY_M)]
trm.men.tot <- men.tot[AGE_M >= max.age, .(AGE_M = max.age, MASS = sum(MASS)), by = .(MSA, COLLEGE_M, MINORITY_M)]
trm.wom.sng <- wom.sng[AGE_F >= max.age, .(AGE_F = max.age, MASS = sum(MASS)), by = .(MSA, COLLEGE_F, MINORITY_F)]
trm.wom.tot <- wom.tot[AGE_F >= max.age, .(AGE_F = max.age, MASS = sum(MASS)), by = .(MSA, COLLEGE_F, MINORITY_F)]

# drop and merge truncated masses
men.sng <- merge(men.sng[AGE_M < max.age], trm.men.sng, all=TRUE, by=c("MSA", "AGE_M", "COLLEGE_M", "MINORITY_M", "MASS"))
men.tot <- merge(men.tot[AGE_M < max.age], trm.men.tot, all=TRUE, by=c("MSA", "AGE_M", "COLLEGE_M", "MINORITY_M", "MASS"))
wom.sng <- merge(wom.sng[AGE_F < max.age], trm.wom.sng, all=TRUE, by=c("MSA", "AGE_F", "COLLEGE_F", "MINORITY_F", "MASS"))
wom.tot <- merge(wom.tot[AGE_F < max.age], trm.wom.tot, all=TRUE, by=c("MSA", "AGE_F", "COLLEGE_F", "MINORITY_F", "MASS"))


# marriages
# TODO: batch smoothing for marriages
# extremely slow with product kernel: on the order of 15h per MSA (feasible for ~10 MSAs)
# try sample splitting, 16 categorical combos
# -> formula = MARRIAGES ~ AGE_M + AGE_F,... by = .(MSA, COLLEGE_M, MINORITY_M, COLLEGE_F, MINORITY_F)
# if still slow, just limit to top MSAs

marriages[, MASS := predict(npreg(bws=npregbw(formula = MARRIAGES ~ AGE_M + COLLEGE_M + MINORITY_M + AGE_F + COLLEGE_F + MINORITY_F,
											  regtype="ll",
											  bwmethod="cv.aic",
											  data=.SD)),
							newdata=.SD),
          by = MSA]

# both >= max.age
trm.both.mar <- marriages[AGE_M >= max.age & AGE_F >= max.age,
						  .(AGE_M = max.age, AGE_F = max.age, MASS := sum(MASS)),
						  by = .(MSA, COLLEGE_M, MINORITY_M, COLLEGE_F, MINORITY_F)]
# husband >= max.age, wife < max.age: group within each wife age
trm.husb <- marriages[AGE_M >= max.age & AGE_F < max.age,
					  .(AGE_M = max.age, MASS := sum(MASS)),
					  by = .(MSA, COLLEGE_M, MINORITY_M, AGE_F, COLLEGE_F, MINORITY_F)]
# husband < max.age, wife >= max.age: group within each husband age
trm.wife <- marriages[AGE_M < max.age & AGE_F >= max.age,
					  .(AGE_F = max.age, MASS := sum(MASS)),
					  by = .(MSA, AGE_M, COLLEGE_M, MINORITY_M, COLLEGE_F, MINORITY_F)]

# drop and merge each segment
marr <- merge(marriages[AGE_M < max.age & AGE_F < max.age], trm.both.mar, all=TRUE,
			  by=c("MSA", "AGE_M", "COLLEGE_M", "MINORITY_M", "AGE_F", "COLLEGE_F", "MINORITY_F", "MASS"))
marr <- merge(marr, trm.husb, all=TRUE,
			  by=c("MSA", "AGE_M", "COLLEGE_M", "MINORITY_M", "AGE_F", "COLLEGE_F", "MINORITY_F", "MASS"))
marr <- merge(marr, trm.wife, all=TRUE,
			  by=c("MSA", "AGE_M", "COLLEGE_M", "MINORITY_M", "AGE_F", "COLLEGE_F", "MINORITY_F", "MASS"))


#TODO: verify no zeros, fill with ones if necessary


### Save output ###

fwrite(marr[, .(MSA, AGE_M, COLLEGE_M, MINORITY_M, AGE_F, COLLEGE_F, MINORITY_F, MASS)], file="data/marriages.csv")
fwrite(men.tot[, .(MSA, AGE_M, COLLEGE_M, MINORITY_M, MASS)], file="data/men-total.csv")
fwrite(wom.tot[, .(MSA, AGE_F, COLLEGE_F, MINORITY_F, MASS)], file="data/wom-total.csv")
fwrite(men.sng[, .(MSA, AGE_M, COLLEGE_M, MINORITY_M, MASS)], file="data/men-single.csv")
fwrite(wom.sng[, .(MSA, AGE_F, COLLEGE_F, MINORITY_F, MASS)], file="data/wom-single.csv")


### Plots ###

# men.msa: table for men in one msa, with a col of smoothed pop
# to plot masses for each (college, minority) combo:
#ggplot(men.msa, aes(x=AGE_M, y=MASS)) + geom_line() + facet_grid(COLLEGE_M ~ MINORITY_M)
# couples:
#ggplot(mar.msa, aes(x=AGE_M, y=AGE_F, z=MASS)) + geom_contour() + facet_grid(COLLEGE_M + MINORITY_M ~ COLLEGE_F + MINORITY_F)

# Non-parametric (local-linear) regression to smooth flows

# Order of traits: husband then wife; age, edu, race
# Notation: husband gets _SP suffix, wife is default
# Min age is 25 because of endogeneity of college
# ACS data: 2008-2014, ages 18-79 (note: AGE_SP is not limited).

library(DBI) # RSQLite database functions
library(data.table)
library(np) # non-parametric regression

indiv.bw <- TRUE # set to FALSE to use bw=cv.aic
age.bw <- 6 # bw to use for manual age smoothing

# connect to sqlite database
# table name: acs
db <- dbConnect(RSQLite::SQLite(), 'data/acs_08-14.db')


### Categorization ###

case_minority <- ' case when "RACESING" in (1,4) and "HISPAN" = 0 then 1 else 2 end '
case_college <- ' case when "EDUC" >= 10 then 2 else 1 end '

# top 20 largest MSAs
top.msa <- c(35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460, 41860, 19820, 38060, 40140, 42660, 33460, 41740, 45300, 41180, 12580)
top_msa <- ' MSA in (35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460, 41860, 19820, 38060, 40140, 42660, 33460, 41740, 45300, 41180, 12580) '

# terminal age (must be same in all files)
max.age <- 65
min.age <- 25

# complete grids to merge in
type.grid <- CJ(MSA = top.msa, SEX = 1:2, AGE = min.age:79, COLLEGE = 1:2, MINORITY = 1:2)


### Queries ###

# flows: counts of new marriages by MSA, SEX, and type
# MF(x), MF(y)
qry_mar_flow <- paste0('select "MET2013" as MSA, "SEX" as SEX,
	"AGE" as AGE, ', case_college, ' as COLLEGE,', case_minority, ' as MINORITY,
	sum("PERWT") as MF
	from acs
	where "MARRINYR" = 2 and AGE >= ', min.age, ' and ', top_msa,
	'group by MSA, SEX, AGE, COLLEGE, MINORITY')

# DF(x), DF(y)
qry_div_flow <- paste0('select "MET2013" as MSA, "SEX" as SEX,
	"AGE" as AGE, ', case_college, ' as COLLEGE,', case_minority, ' as MINORITY,
	sum("PERWT") as DF
	from acs
	where "DIVINYR" = 2 and AGE >= ', min.age, ' and ', top_msa,
	'group by MSA, SEX, AGE, COLLEGE, MINORITY')

# queries return dataframes, convert to data.table and merge into the complete grid
mar.flow <- merge(data.table(dbGetQuery(db, qry_mar_flow)), type.grid,
				all=TRUE, by=c("MSA", "SEX", "AGE", "COLLEGE", "MINORITY"))
div.flow <- merge(data.table(dbGetQuery(db, qry_div_flow)), type.grid,
				all=TRUE, by=c("MSA", "SEX", "AGE", "COLLEGE", "MINORITY"))

# fill NA with zeros before smoothing
mar.flow[is.na(MF), MF := 0]
div.flow[is.na(DF), DF := 0]


### Smoothing by Non-parametric Regression ###

# convert categorical variables to factors for npreg
for (col in c('COLLEGE', 'MINORITY')) set(mar.flow, j = col, value = factor(mar.flow[[col]]))
for (col in c('COLLEGE', 'MINORITY')) set(div.flow, j = col, value = factor(div.flow[[col]]))

# batch smoothing by MSA and SEX: add FLOW column to each data table
if (indiv.bw) {
	mar.flow[, FLOW := predict(npreg(bws=age.bw, txdat=AGE, tydat=MF, regtype="ll")),
			   by = .(MSA, SEX, COLLEGE, MINORITY)]
	div.flow[, FLOW := predict(npreg(bws=age.bw, txdat=AGE, tydat=DF, regtype="ll")),
			   by = .(MSA, SEX, COLLEGE, MINORITY)]
} else {
	mar.flow[, FLOW := predict(npreg(bws=npregbw(formula = MF ~ AGE + COLLEGE + MINORITY,
												 regtype="ll",
												 bwmethod="cv.aic",
												 data=.SD)),
							   newdata=.SD),
			   by = .(MSA, SEX)]
	div.flow[, FLOW := predict(npreg(bws=npregbw(formula = DF ~ AGE + COLLEGE + MINORITY,
												 regtype="ll",
												 bwmethod="cv.aic",
												 data=.SD)),
							   newdata=.SD),
			   by = .(MSA, SEX)]
}


### Trim and Clean ###

# trim max.age (NOTE: max.age will be dropped in the estimation because of truncation)
mar.flow <- mar.flow[AGE <= max.age]
div.flow <- div.flow[AGE <= max.age]

# llr may give negative or tiny values: truncate to min of 1
mar.flow[FLOW < 1.0, FLOW := 1.0]
div.flow[FLOW < 1.0, FLOW := 1.0]


### Save flows and stocks ###

# by MSA and type
fwrite(mar.flow[SEX == 1, .(MSA, AGE, COLLEGE, MINORITY, FLOW)], file="data/men-MF.csv")
fwrite(mar.flow[SEX == 2, .(MSA, AGE, COLLEGE, MINORITY, FLOW)], file="data/wom-MF.csv")
fwrite(div.flow[SEX == 1, .(MSA, AGE, COLLEGE, MINORITY, FLOW)], file="data/men-DF.csv")
fwrite(div.flow[SEX == 2, .(MSA, AGE, COLLEGE, MINORITY, FLOW)], file="data/wom-DF.csv")

# Non-parametric (local-linear) regression to smooth population stocks

# Order of traits: husband then wife; age, edu, race
# Notation: husband gets _SP suffix, wife is default
# Min age is 25 because of endogeneity of college
# ACS data: 2008-2015, ages 18-79 (note: AGE_SP is not limited).

library(DBI) # RSQLite database functions
library(data.table)
library(np) # non-parametric regression


### Usage Options ###

# use edu and race as types? modify output filename at bottom
#	strategy: keep redundant types for compat
static.types <- FALSE

# choose what bandwidth to use for ages in marriages and for individuals
marr.bw <- FALSE # set to FALSE to use bw=cv.aic
indiv.bw <- FALSE # set to FALSE to use bw=cv.aic

age.bw <- 3 # bw to use for manual age smoothing

# connect to sqlite database
# table names: acs, mig2met
# NOTE: edit queries to adjust number of years to divide masses by
db <- dbConnect(RSQLite::SQLite(), 'data/acs_08-15.db')

# top 20 largest MSAs
top.msa <- c(35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460, 41860, 19820, 38060, 40140, 42660, 33460, 41740, 45300, 41180, 12580)
top_msa <- ' (35620, 31080, 16980, 19100, 37980, 26420, 47900, 33100, 12060, 14460, 41860, 19820, 38060, 40140, 42660, 33460, 41740, 45300, 41180, 12580) '

# terminal age
max.age <- 65
min.age <- 18


### Categorization ###

case_minority <- ' case when "RACESING" in (1,4) and "HISPAN" = 0 then 1 else 2 end '
case_minority_sp <- ' case when "RACESING_SP" in (1,4) and "HISPAN_SP" = 0 then 1 else 2 end '
case_college <- ' case when "EDUC" >= 10 then 2 else 1 end '
case_college_sp <- ' case when "EDUC_SP" >= 10 then 2 else 1 end '

# customize queries for desired types
if (static.types) {
	ind.types <- paste0(case_college, 'as COLLEGE,', case_minority, 'as MINORITY,')
	husb.types <- paste0(case_college_sp, 'as COLLEGE_M,', case_minority_sp, 'as MINORITY_M,')
	wife.types <- paste0(case_college, 'as COLLEGE_F,', case_minority, 'as MINORITY_F,')
} else {
	ind.types <- '1 as COLLEGE, 1 as MINORITY,'
	husb.types <- '1 as COLLEGE_M, 1 as MINORITY_M,'
	wife.types <- '1 as COLLEGE_F, 1 as MINORITY_F,'
}

ind.grp <- ', COLLEGE, MINORITY'
husb.grp <- ', COLLEGE_M, MINORITY_M'
wife.grp <- ', COLLEGE_F, MINORITY_F'

ind.mrg <- c("COLLEGE", "MINORITY")
husb.mrg <- c("COLLEGE_M", "MINORITY_M")
wife.mrg <- c("COLLEGE_F", "MINORITY_F")

# complete grids to merge in
if (static.types) {
	ind.grid <- CJ(MSA = top.msa, SEX = 1:2, AGE = min.age:79, COLLEGE = 1:2, MINORITY = 1:2)
	mar.grid <- CJ(MSA = top.msa, AGE_M = min.age:79, COLLEGE_M = 1:2, MINORITY_M = 1:2,
				   AGE_F = min.age:79, COLLEGE_F = 1:2, MINORITY_F = 1:2)
} else {
	ind.grid <- CJ(MSA = top.msa, SEX = 1:2, AGE = min.age:79, COLLEGE = 1, MINORITY = 1)
	mar.grid <- CJ(MSA = top.msa, AGE_M = min.age:79, COLLEGE_M = 1, MINORITY_M = 1,
				   AGE_F = min.age:79, COLLEGE_F = 1, MINORITY_F = 1)
}


### Queries ###

# total population counts
qry_tot <- paste0('select "MET2013" as MSA, "SEX" as SEX,
					  "AGE" as AGE, ', ind.types,
					  ' sum("PERWT")/8 as RAW_POP
				  from acs
				  where AGE >= ', min.age, ' and MSA in ', top_msa, ' and "GQ" < 3
				  group by MSA, SEX, AGE', ind.grp)

# singles counts
qry_sng <- paste0('select "MET2013" as MSA, "SEX" as SEX,
					  "AGE" as AGE, ', ind.types,
					  ' sum("PERWT")/8 as RAW_SNG
				  from acs
				  where AGE >= ', min.age, ' and "MARST" >= 3 and MSA in ', top_msa, ' and "GQ" < 3
				  group by MSA, SEX, AGE', ind.grp)

# married counts
qry_marr <- paste0('select "MET2013" as MSA,
					   "AGE_SP" as AGE_M, ', husb.types,
					   ' "AGE" as AGE_F, ', wife.types,
					   ' sum("HHWT")/8 as RAW_MASS
				   from acs
				   where "MARST" <= 2 and "SEX" = 2
					   and AGE_F >=', min.age, '  and AGE_M >= ', min.age,
					   ' and MSA in ', top_msa,
				   'group by MSA, AGE_M', husb.grp, ', AGE_F', wife.grp)

# MF(x,y): counts of new marriages by MSA and type pair
qry_marr_flow <- paste0('select "MET2013" as MSA,
							"AGE_SP" as AGE_M, ', husb.types,
							' "AGE" as AGE_F, ', wife.types,
							' sum("HHWT")/8 as RAW_MF
						from acs
						where "MARRINYR" = 2 and "MARST" <= 2
							and AGE_F >= ', min.age, ' and AGE_M >= ', min.age,
							' and MSA in ', top_msa, ' and "SEX" = 2
						group by MSA, AGE_M', husb.grp,', AGE_F', wife.grp)

# MF(x), MF(y)
#qry_mar_flow <- paste0('select "MET2013" as MSA, "SEX" as SEX,
#	"AGE" as AGE, ', case_college, ' as COLLEGE,', case_minority, ' as MINORITY,
#	sum("PERWT") as MF
#	from acs
#	where "MARRINYR" = 2 and AGE >= ', min.age, ' and MSA in ', top_msa,
#	'group by MSA, SEX, AGE, COLLEGE, MINORITY')

# DF(x), DF(y)
qry_div_flow <- paste0('select "MET2013" as MSA, "SEX" as SEX,
						   "AGE" as AGE, ', ind.types,
						   ' sum("PERWT")/8 as RAW_DF
					   from acs
					   where "DIVINYR" = 2 and AGE >= ', min.age, ' and MSA in ', top_msa,
					   ' group by MSA, SEX, AGE', ind.grp)

# migration inflow: counts of people who weren't in current MSA last year
qry_inflow <- paste0('select acs."MET2013" as MSA, acs."SEX" as SEX,
						 acs."AGE" as AGE, ', ind.types, 
						 ' sum(acs."PERWT")/4 as INFLOW_RAW
                     from acs
                     left join mig2met as m on acs."MIGPUMA1"=m."MIGPUMA1" and acs."MIGPLAC1"=m."MIGPLAC1"
                     where acs."AGE" >= ', min.age, ' and acs."MET2013" in ', top_msa, ' and m."MSA" != acs."MET2013"
                         and acs."MIGRATE1D" >= 24 and acs."GQ" < 3 and acs."YEAR" >= 2012
                     group by acs."MET2013", acs."SEX", acs."AGE"', ind.grp)

# migration outflow: counts of people who were in MSA last year but no longer
qry_outflow <- paste0('select m."MSA" as MSA, acs."SEX" as SEX,
						  acs."AGE" as AGE, ', ind.types,
						  ' sum(acs."PERWT")/4 as OUTFLOW_RAW
                      from acs
                      left join mig2met as m on acs."MIGPUMA1"=m."MIGPUMA1" and acs."MIGPLAC1"=m."MIGPLAC1"
                      where acs."AGE" >= ', min.age, ' and m."MSA" in ', top_msa, ' and m."MSA" != acs."MET2013"
						  and acs."MIGRATE1D" >= 24 and acs."GQ" < 3 and acs."YEAR" >= 2012
                      group by m."MSA", acs."SEX", acs."AGE"', ind.grp)

# couple migration inflows
qry_inmar <- paste0('select acs."MET2013" as MSA,
						acs."AGE_SP" as AGE_M, ', husb.types,
						' acs."AGE" as AGE_F, ', wife.types,
						' sum(acs."HHWT")/4 as INFLOW_RAW
                     from acs
                     left join mig2met as m on acs."MIGPUMA1"=m."MIGPUMA1" and acs."MIGPLAC1"=m."MIGPLAC1"
                     where acs."MET2013" in ', top_msa, ' and m."MSA" != acs."MET2013"
                        and "MARST" <= 2 and "SEX" = 2
                        and AGE_F >= ', min.age, ' and AGE_M >= ', min.age,
						' and acs."MIGRATE1D" >= 24 and acs."YEAR" >= 2012
                     group by acs."MET2013", acs."AGE_SP"', husb.grp,', acs."AGE"', wife.grp)

# couple migration outflows
qry_outmar <- paste0('select m."MSA" as MSA,
						acs."AGE_SP" as AGE_M, ', husb.types,
						' acs."AGE" as AGE_F, ', wife.types,
						' sum(acs."HHWT")/4 as OUTFLOW_RAW
                     from acs
                     left join mig2met as m on acs."MIGPUMA1"=m."MIGPUMA1" and acs."MIGPLAC1"=m."MIGPLAC1"
                     where m."MSA" in ', top_msa, '  and m."MSA" != acs."MET2013"
                        and "MARST" <= 2 and "SEX" = 2
                        and AGE_F >= ', min.age, ' and AGE_M >= ', min.age,
                        ' and acs."MIGRATE1D" >= 24 and acs."YEAR" >= 2012
                     group by m."MSA", acs."AGE_SP"', husb.grp,', acs."AGE"', wife.grp)

# queries return dataframes, convert to data.table and merge into the complete grid
pop.dt <- merge(merge(data.table(dbGetQuery(db, qry_tot)), data.table(dbGetQuery(db, qry_sng)),
					  all=TRUE, by=c("MSA", "SEX", "AGE", ind.mrg)),
				ind.grid,
				all.y=TRUE, by=c("MSA", "SEX", "AGE", ind.mrg))

# right join to force to grid (drops husbands > 79)
marriages <- merge(data.table(dbGetQuery(db, qry_marr)), mar.grid,
				   all.y=TRUE, by=c("MSA", "AGE_M", husb.mrg, "AGE_F", wife.mrg))

mar.flow <- merge(data.table(dbGetQuery(db, qry_marr_flow)), mar.grid,
				  all.y=TRUE, by=c("MSA", "AGE_M", husb.mrg, "AGE_F", wife.mrg))
div.flow <- merge(data.table(dbGetQuery(db, qry_div_flow)), ind.grid,
				  all.y=TRUE, by=c("MSA", "SEX", "AGE", ind.mrg))

ind.mig <- merge(merge(data.table(dbGetQuery(db, qry_inflow)), data.table(dbGetQuery(db, qry_outflow)),
					   all=TRUE, by=c("MSA", "SEX", "AGE", ind.mrg)),
				 ind.grid,
				 all.y=TRUE, by=c("MSA", "SEX", "AGE", ind.mrg))

mar.mig <- merge(merge(data.table(dbGetQuery(db, qry_inmar)), data.table(dbGetQuery(db, qry_outmar)),
					   all=TRUE, by=c("MSA", "AGE_M", husb.mrg, "AGE_F", wife.mrg)),
				 mar.grid,
				 all.y=TRUE, by=c("MSA", "AGE_M", husb.mrg, "AGE_F", wife.mrg))

# fill NA with zeros before smoothing
pop.dt[is.na(pop.dt)] <- 0
marriages[is.na(RAW_MASS), RAW_MASS := 0]
mar.flow[is.na(RAW_MF), RAW_MF := 0]
div.flow[is.na(RAW_DF), RAW_DF := 0]
ind.mig[is.na(ind.mig)] <- 0
mar.mig[is.na(mar.mig)] <- 0

ind.mig[, NET_OUTFLOW_RAW := OUTFLOW_RAW - INFLOW_RAW]
mar.mig[, NET_OUTFLOW_RAW := OUTFLOW_RAW - INFLOW_RAW]


### Smoothing by Non-parametric Regression ###

# npregbw cv methods give following bandwidths on age per msa
#	* wom.sng: 2.2
#	* marriages (split sample): 0.5 for each age (3 hour runtime)
#	* marriages (product kernel): takes >24h for 1/5 starts... aborted

# convert categorical variables to factors for npreg
for (col in ind.mrg) set(pop.dt, j = col, value = factor(pop.dt[[col]]))
for (col in ind.mrg) set(div.flow, j = col, value = factor(div.flow[[col]]))
for (col in ind.mrg) set(ind.mig, j = col, value = factor(ind.mig[[col]]))
for (col in c(husb.mrg, wife.mrg)) set(marriages, j = col, value = factor(marriages[[col]]))
for (col in c(husb.mrg, wife.mrg)) set(mar.flow, j = col, value = factor(mar.flow[[col]]))
for (col in c(husb.mrg, wife.mrg)) set(mar.mig, j = col, value = factor(mar.mig[[col]]))

if (indiv.bw) {
	# manual smoothing with sample splitting
	pop.dt[, `:=`(SNG = predict(npreg(bws=age.bw, txdat=AGE, tydat=RAW_SNG, regtype="ll")),
				  POP = predict(npreg(bws=age.bw, txdat=AGE, tydat=RAW_POP, regtype="ll"))),
			by = c("MSA", "SEX", ind.mrg)] # typically use list for `by`, but need vector here to join

	div.flow[, FLOW := predict(npreg(bws=age.bw,
									 txdat=AGE, tydat=RAW_DF, regtype="ll")),
			   by = c("MSA", "SEX", ind.mrg)]

	ind.mig[, NET_OUTFLOW := predict(npreg(bws=age.bw, txdat=AGE, tydat=NET_OUTFLOW_RAW, regtype="ll")),
			by = c("MSA", "SEX", ind.mrg)]
} else {
	# cross-validated bandwidth (product kernel)
	pop.dt[, `:=`(SNG = predict(npreg(bws=npregbw(formula = RAW_SNG ~ AGE,
												  regtype="ll",
												  bwmethod="cv.aic",
												  data=.SD))),
			      POP = predict(npreg(bws=npregbw(formula = RAW_POP ~ AGE,
												  regtype="ll",
												  bwmethod="cv.aic",
												  data=.SD)))),
			by = c("MSA", "SEX", ind.mrg)]

	div.flow[, FLOW := predict(npreg(bws=npregbw(formula = RAW_DF ~ AGE,
												 regtype="ll",
												 bwmethod="cv.aic",
												 data=.SD))),
			  by = c("MSA", "SEX", ind.mrg)]

	# Note: untested
	ind.mig[, NET_OUTFLOW := predict(npreg(bws=npregbw(formula = NET_OUTFLOW_RAW ~ AGE,
													   regtype="ll",
													   bwmethod="cv.aic",
													   data=.SD))),
			 by = c("MSA", "SEX", ind.mrg)]
}

## combine groups with age >= max.age (truncation at T)

# lump together truncated masses to be merged back in
trm.pop <- pop.dt[AGE >= max.age,
				   .(AGE = max.age, SNG = sum(SNG), POP = sum(POP)),
				   by = c("MSA", "SEX", ind.mrg)]

# drop and merge truncated masses
pop.dt <- merge(pop.dt[AGE < max.age], trm.pop,
				 all=TRUE, by=c("MSA", "SEX", "AGE", ind.mrg, "SNG", "POP"))

# lump together truncated masses to be merged back in
trm.ind.mig <- ind.mig[AGE >= max.age,
					   .(AGE = max.age, NET_OUTFLOW = sum(NET_OUTFLOW)),
					   by = c("MSA", "SEX", ind.mrg)]

# drop and merge truncated masses
ind.mig <- merge(ind.mig[AGE < max.age], trm.ind.mig,
				 all=TRUE, by=c("MSA", "SEX", "AGE", ind.mrg, "NET_OUTFLOW"))


#	But npreg didn't get past 1/5 in 15h. The full grid (zero-filled) is that much larger.
#	Idea: first do a fast search with generous tolerances, then do a precise search from that starting value 

# marriages

# sample splitting, 16 categorical combos: tolerable (on the order of 10m per group but 16 groups per MSA)
#	* using all 6 variables (product kernels) is extremely slow due to curse of dimensionality (multiple days per MSA)
#	* try: use a shell script to launch 10 different R instances with one MSA each

print("Starting non-parametric regression on couple masses.")

if (marr.bw) {
	# manual smoothing with sample splitting
	marriages[, MASS := predict(npreg(bws=c(age.bw, age.bw),
									  txdat=.(AGE_M, AGE_F), tydat=RAW_MASS,
									  regtype="ll")),
			   by = c("MSA", husb.mrg, wife.mrg)]

	mar.flow[, FLOW := predict(npreg(bws=c(age.bw, age.bw),
									 txdat=.(AGE_M, AGE_F), tydat=RAW_MF,
									 regtype="ll")),
			  by = c("MSA", husb.mrg, wife.mrg)]

	mar.mig[, NET_OUTFLOW := predict(npreg(bws=c(age.bw, age.bw),
										   txdat=.(AGE_M, AGE_F), tydat=NET_OUTFLOW_RAW,
										   regtype="ll")),
		     by = c("MSA", husb.mrg, wife.mrg)]
} else {
	# cross-validated bandwidth with sample splitting
	print("CV bandwidth will take a while...")

	marriages[, MASS := predict(npreg(bws=npregbw(formula = RAW_MASS ~ AGE_M + AGE_F,
												  regtype="ll",
												  bwmethod="cv.aic",
												  data=.SD))),
				by = c("MSA", husb.mrg, wife.mrg)]

	mar.flow[, FLOW := predict(npreg(bws=npregbw(formula = RAW_MF ~ AGE_M + AGE_F,
												 regtype="ll",
												 bwmethod="cv.aic",
												 data=.SD))),
			   by = c("MSA", husb.mrg, wife.mrg)]

	mar.mig[, NET_OUTFLOW := predict(npreg(bws=npregbw(formula = NET_OUTFLOW_RAW ~ AGE_M + AGE_F,
													   regtype="ll",
													   bwmethod="cv.aic",
													   data=.SD))),
			 by = c("MSA", husb.mrg, wife.mrg)]
}
# full product kernel on all 6 variables
#marriages[, MASS := predict(npreg(bws=npregbw(formula = MARRIAGES ~ AGE_M + COLLEGE_M + MINORITY_M + AGE_F + COLLEGE_F + MINORITY_F,
#											  regtype="ll",
#											  bwmethod="cv.aic",
#											  data=.SD)),
#							newdata=.SD),
#          by = MSA]


### Trim, Truncate, and Clean ###

# both >= max.age
trm.both.mar <- marriages[AGE_M >= max.age & AGE_F >= max.age,
						  .(AGE_M = max.age, AGE_F = max.age, MASS = sum(MASS)),
						  by = c("MSA", husb.mrg, wife.mrg)]
# husband >= max.age, wife < max.age: group within each wife age
trm.husb.mar <- marriages[AGE_M >= max.age & AGE_F < max.age,
					  .(AGE_M = max.age, MASS = sum(MASS)),
					  by = c("MSA", husb.mrg, "AGE_F", wife.mrg)]
# husband < max.age, wife >= max.age: group within each husband age
trm.wife.mar <- marriages[AGE_M < max.age & AGE_F >= max.age,
					  .(AGE_F = max.age, MASS = sum(MASS)),
					  by = c("MSA", "AGE_M", husb.mrg, wife.mrg)]

# both >= max.age
trm.both.mig <- mar.mig[AGE_M >= max.age & AGE_F >= max.age,
						.(AGE_M = max.age, AGE_F = max.age, NET_OUTFLOW = sum(NET_OUTFLOW)),
						by = c("MSA", husb.mrg, wife.mrg)]
# husband >= max.age, wife < max.age: group within each wife age
trm.husb.mig <- mar.mig[AGE_M >= max.age & AGE_F < max.age,
					.(AGE_M = max.age, NET_OUTFLOW = sum(NET_OUTFLOW)),
					by = c("MSA", husb.mrg, "AGE_F", wife.mrg)]
# husband < max.age, wife >= max.age: group within each husband age
trm.wife.mig <- mar.mig[AGE_M < max.age & AGE_F >= max.age,
					.(AGE_F = max.age, NET_OUTFLOW = sum(NET_OUTFLOW)),
					by = c("MSA", "AGE_M", husb.mrg, wife.mrg)]

# drop and merge each segment
marriages <- merge(marriages[AGE_M < max.age & AGE_F < max.age], trm.both.mar, all=TRUE,
				   by=c("MSA", "AGE_M", husb.mrg, "AGE_F", wife.mrg, "MASS"))
marriages <- merge(marriages, trm.husb.mar, all=TRUE,
			  by=c("MSA", "AGE_M", husb.mrg, "AGE_F", wife.mrg, "MASS"))
marriages <- merge(marriages, trm.wife.mar, all=TRUE,
			  by=c("MSA", "AGE_M", husb.mrg, "AGE_F", wife.mrg, "MASS"))

mar.mig <- merge(mar.mig[AGE_M < max.age & AGE_F < max.age], trm.both.mig, all=TRUE,
				 by=c("MSA", "AGE_M", husb.mrg, "AGE_F", wife.mrg, "NET_OUTFLOW"))
mar.mig <- merge(mar.mig, trm.husb.mig, all=TRUE,
				 by=c("MSA", "AGE_M", husb.mrg, "AGE_F", wife.mrg, "NET_OUTFLOW"))
mar.mig <- merge(mar.mig, trm.wife.mig, all=TRUE,
				 by=c("MSA", "AGE_M", husb.mrg, "AGE_F", wife.mrg, "NET_OUTFLOW"))

# trim max.age (NOTE: max.age will be dropped in the estimation because of truncation)
mar.flow <- mar.flow[AGE_M <= max.age & AGE_F <= max.age]
div.flow <- div.flow[AGE <= max.age]

# llr may give negative or tiny values: truncate to min of 0
marriages[MASS < 0, MASS := 0.0]
mar.flow[FLOW < 0, FLOW := 0.0]
div.flow[FLOW < 0, FLOW := 0.0]


### Save flows and stocks ###

fwrite(marriages, file="data/ageonly/marriages.csv")
fwrite(pop.dt, file="data/ageonly/pop.csv")
fwrite(mar.flow, file="data/ageonly/pair-MF.csv")
fwrite(div.flow, file="data/ageonly/ind-DF.csv")
fwrite(ind.mig, file="data/ageonly/indiv-migration.csv")
fwrite(mar.mig, file="data/ageonly/mar-migration.csv")

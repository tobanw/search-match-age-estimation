library(data.table)
library(zoo) # for interpolating mortality rates with na.approx

# use edu and race as types? modify output filename at bottom
#	strategy: keep redundant types for compat
static.types <- TRUE

# terminal age: must match `smooth_pop.r`
max.age <- 65
min.age <- 25 # exclusive (inflow group)

# load up CDC mortality table
if (static.types) { # by sex, age, and minority
	mort.dt <- fread('data/mort-full_08-15.csv')

	# full grid to merge and interpolate mortality rates
	ind.grid <- CJ(AGE = min.age:90, COLLEGE = 1:2, MINORITY = 1:2)
} else { # by sex, age
	mort.dt <- fread('data/mort-ageonly_08-15.csv')
	mort.dt[, MINORITY := 1] # for compat

	ind.grid <- CJ(AGE = min.age:90, COLLEGE = 1, MINORITY = 1)
}


### Annual Death Rates ###

# prep mortality table for merging
if (static.types) mort.dt[, MINORITY := MINORITY + 1] # convert from 0:1 to 1:2 
mort.dt[, DTHRT := DR100 / 100000] # convert individual probabilities
mort.dt[, AGE := AGE_L + 5] # midpoints for interpolation: 30,40,50,60,70
#mort.dt[AGE_H==24, AGE := 25] # boundaries for interpolation
#mort.dt[AGE_H==94, AGE := 79] # boundaries for interpolation

# merge sparse mortality rates into full age table
men.death <- merge(ind.grid, mort.dt[SEX==1, .(AGE, MINORITY, DTHRT)],
				   by=c("AGE", "MINORITY"), all=TRUE)
wom.death <- merge(ind.grid, mort.dt[SEX==2, .(AGE, MINORITY, DTHRT)],
				   by=c("AGE", "MINORITY"), all=TRUE)

# fill youngest age boundaries with youngest mort rates (otherwise interpolation misbehaves)
men.death[AGE == min.age, DTHRT := mort.dt[SEX == 1 & AGE_L == 20, DTHRT]]
wom.death[AGE == min.age, DTHRT := mort.dt[SEX == 2 & AGE_L == 20, DTHRT]]
# oldest age 90 also coincides with 85-99 group in data

# interpolate death rates over ages within minority group
men.death[, DTHRT := na.approx(DTHRT, AGE), MINORITY]
wom.death[, DTHRT := na.approx(DTHRT, AGE), MINORITY]

# death arrival rates: ψ = -log(1-mort) ≈ mort (when mort ≈ 0)
men.death[, PSI := -log(1 - DTHRT)]
wom.death[, PSI := -log(1 - DTHRT)]

# truncation at max.age: use mean
trm.men.death <- men.death[AGE >= max.age,
						   .(AGE = max.age, PSI = mean(PSI)),
						   by = .(COLLEGE, MINORITY)]
men.death <- merge(men.death[AGE < max.age], trm.men.death,
				 all=TRUE, by=c("AGE", "COLLEGE", "MINORITY", "PSI"))

trm.wom.death <- wom.death[AGE >= max.age,
						   .(AGE = max.age, PSI = mean(PSI)),
						   by = .(COLLEGE, MINORITY)]
wom.death <- merge(wom.death[AGE < max.age], trm.wom.death,
				 all=TRUE, by=c("AGE", "COLLEGE", "MINORITY", "PSI"))

# save tables to csv
if (static.types) { # by sex, age, and minority
	fwrite(men.death[, .(AGE, COLLEGE, MINORITY, PSI)], file="data/racedu/men-psi.csv")
	fwrite(wom.death[, .(AGE, COLLEGE, MINORITY, PSI)], file="data/racedu/wom-psi.csv")
} else { # by sex, age
	fwrite(men.death[, .(AGE, COLLEGE, MINORITY, PSI)], file="data/ageonly/men-psi.csv")
	fwrite(wom.death[, .(AGE, COLLEGE, MINORITY, PSI)], file="data/ageonly/wom-psi.csv")
}

library(data.table)
library(ggplot2)

# SELECT model
age.only <- TRUE
data.path <- "data/ageonly16-13/" # input dir: smoothed data
estim.path <- "results/estimates-csv/dynamic-const/ageonly16-13/" # input dir: model estimates
plot.path <- "results/result-plots/dynamic-const/ageonly/" # output dir: plots
img.ext <- ".pdf"

# Cities for faceted plots
top.cities <- c(35620, 31080, 16980, 19100)

if(age.only) {
	min.age <- 18
	edu.types <- c(A = 1)
	rac.types <- c(A = 1)
} else {
	min.age <- 25
	edu.types <- c(N = 1, C = 2)
	rac.types <- c(W = 1, M = 2)
}

max.age <- 65


##### SETUP DATA #####

msa.dt <- fread('data/top20msa.csv') # for getting names from codes
setnames(msa.dt, old = "METCODE", new = "MSA") # rename

# joint couple objects
prod.dt <- merge(fread(file.path(estim.path, "prod.csv")), msa.dt, by="MSA")
alpha.dt <- merge(fread(file.path(estim.path, "alpha.csv")), msa.dt, by="MSA")
m.stock.dt <- merge(fread(file.path(data.path, "marriages.csv")), msa.dt, by="MSA")
m.mig.dt <- merge(fread(file.path(data.path, "mar-migration.csv")), msa.dt, by="MSA")
mf.dt <- merge(fread(file.path(estim.path, "mMF.csv")), fread(file.path(data.path, "pair-MF.csv")),
               by = c("MSA", "AGE_M", "COLLEGE_M", "MINORITY_M", "AGE_F", "COLLEGE_F", "MINORITY_F"))
setnames(mf.dt, old = c("FLOW", "VALUE"), new = c("DATA", "MODEL"))

# individual objects
val.dt <- merge(fread(file.path(estim.path, "men_val.csv")), fread(file.path(estim.path, "wom_val.csv")),
                suffixes = c("_M", "_F"), by = c("MSA", "AGE", "COLLEGE", "MINORITY"))

df.model.dt <- merge(fread(file.path(estim.path, "mDF_m.csv")), fread(file.path(estim.path, "mDF_f.csv")),
                     suffixes = c("_M", "_F"), by = c("MSA", "AGE", "COLLEGE", "MINORITY"))
setnames(df.model.dt, old = c("VALUE_M", "VALUE_F"), new = c("FLOW_M", "FLOW_F")) # rename

# cast long to wide; sex from indiv files: ind-DF.csv, pop.csv
df.data.dt <- dcast(fread(file.path(data.path, "ind-DF.csv")), MSA + AGE + COLLEGE + MINORITY ~ SEX, value.var = "FLOW")
setnames(df.data.dt, old = c("1", "2"), new = c("FLOW_M", "FLOW_F"))
df.dt <- merge(df.model.dt, df.data.dt, suffixes = c("_MODEL", "_DATA"), by = c("MSA", "AGE", "COLLEGE", "MINORITY"))

pop.dt <- dcast(fread(file.path(data.path, "pop.csv")), MSA + AGE + COLLEGE + MINORITY ~ SEX,
                value.var = c("SNG", "POP"))
setnames(pop.dt, old = c("SNG_1", "SNG_2", "POP_1", "POP_2"), new = c("SNG_M", "SNG_F", "POP_M", "POP_F"))

# merge in city names
mf.dt <- merge(mf.dt, msa.dt, by="MSA")
val.dt <- merge(val.dt, msa.dt, by="MSA")
pop.dt <- merge(pop.dt, msa.dt, by="MSA")
df.dt <- merge(df.dt, msa.dt, by="MSA")

# convert categoricals to factors
val.dt[COLLEGE==1, EDU:="No college"]
val.dt[COLLEGE==2, EDU:="College"]
val.dt[MINORITY==1, RAC:="Asian or White"]
val.dt[MINORITY==2, RAC:="Minority"]
    
pop.dt[COLLEGE==1, EDU:="No college"]
pop.dt[COLLEGE==2, EDU:="College"]
pop.dt[MINORITY==1, RAC:="Asian or White"]
pop.dt[MINORITY==2, RAC:="Minority"]

df.dt[COLLEGE==1, EDU:="No college"]
df.dt[COLLEGE==2, EDU:="College"]
df.dt[MINORITY==1, RAC:="Asian or White"]
df.dt[MINORITY==2, RAC:="Minority"]

for (col in c('MSA', 'CITY', 'EDU', 'RAC')) set(val.dt, j = col, value = factor(val.dt[[col]]))
for (col in c('MSA', 'CITY', 'EDU', 'RAC')) set(pop.dt, j = col, value = factor(pop.dt[[col]]))
for (col in c('MSA', 'CITY', 'EDU', 'RAC')) set(df.dt, j = col, value = factor(df.dt[[col]]))

for (col in c('METCODE', 'CITY')) set(msa.dt, j = col, value = factor(msa.dt[[col]]))
set(prod.dt, j = 'CITY', value = factor(prod.dt[['CITY']]))
set(alpha.dt, j = 'CITY', value = factor(alpha.dt[['CITY']]))
set(mf.dt, j = 'CITY', value = factor(mf.dt[['CITY']]))
set(m.mig.dt, j = 'CITY', value = factor(m.mig.dt[['CITY']]))
set(m.stock.dt, j = 'CITY', value = factor(m.stock.dt[['CITY']]))



##### PLOTTING #####

# plotting layers
xy.line <- geom_segment(x=min.age, y=min.age, xend=max.age, yend=max.age, colour="white", size=1.2)
age.labs <- labs(x="Husband age", y="Wife age")

flow.color.grad <- scale_colour_gradient2(low="black", mid="blue", high="red", midpoint = 0)
alpha.color.grad <- scale_colour_gradient2(low="navy", mid="red", high="white", midpoint=0.5)
prod.color.grad <- scale_colour_gradient2(low="navy", mid="red", high="white", midpoint=0)

# clipping function for sensible production plots
truncator <- function(x, low, high) {
    x[x>high] <- high
    x[x<low] <- low
    return(x)
}


for (r in names(rac.types)) for (e in names(edu.types)) {

	### MF fit: data and model ###

	# data: top 4
	ggplot(mf.dt[AGE_M < max.age-1 & AGE_F < max.age-1 & MSA %in% top.cities &
		         COLLEGE_M==edu.types[[e]] & MINORITY_M==rac.types[[r]] &
		         COLLEGE_F==edu.types[[e]] & MINORITY_F==rac.types[[r]] ],
	       aes(x=AGE_M, y=AGE_F, z=DATA)) +
		stat_contour(aes(colour=..level..), bins=15, size=0.5) +
		flow.color.grad +
		xy.line +
		age.labs +
		facet_wrap(~CITY)
	ggsave(paste0("mf-data-", r, e, img.ext), path = plot.path)

	# model: top 4
	ggplot(mf.dt[MODEL >= 0 & AGE_M < max.age-1 & AGE_F < max.age-1 & MSA %in% top.cities &
		         COLLEGE_M==edu.types[[e]] & MINORITY_M==rac.types[[r]] &
		         COLLEGE_F==edu.types[[e]] & MINORITY_F==rac.types[[r]] ],
	       aes(x=AGE_M, y=AGE_F, z=MODEL)) +
		stat_contour(aes(colour=..level..), bins=15, size=0.5) +
		flow.color.grad +
		xy.line +
		age.labs +
		facet_wrap(~CITY)
	ggsave(paste0("mf-model-", r, e, img.ext), path = plot.path)

	# fit gap, top 4: (model - data)
	ggplot(mf.dt[AGE_M < max.age-1 & AGE_F < max.age-1 & MSA %in% top.cities &
		         COLLEGE_M==edu.types[[e]] & MINORITY_M==rac.types[[r]] &
		         COLLEGE_F==edu.types[[e]] & MINORITY_F==rac.types[[r]],
		         .(AGE_M, AGE_F, FIT = MODEL - DATA, CITY)],
		   aes(x=AGE_M, y=AGE_F, z=FIT)) +
		stat_contour(aes(colour=..level..), size=0.5) +
		flow.color.grad +
		xy.line +
		age.labs +
		facet_wrap(~CITY)
	ggsave(paste0("mf-fit-", r, e, img.ext), path = plot.path)


	### DF fit: data and model ###

	ggplot(df.dt[MSA %in% top.cities & AGE<max.age &
		         COLLEGE==edu.types[[e]] & MINORITY==rac.types[[r]] ],
		   aes(x=AGE)) +
		geom_line(aes(y=FLOW_M_DATA, color = "Male"), size=1) +
		geom_line(aes(y=FLOW_F_DATA, color = "Female"), size=1) +
		geom_line(aes(y=FLOW_M_MODEL, color = "Male"), size=1, linetype=2) +
		geom_line(aes(y=FLOW_F_MODEL, color = "Female"), size=1, linetype=2) +
		labs(x="Age", y="Divorces", color="Sex") +
		facet_wrap(~CITY, scales = "free_y")
	ggsave(paste0("df-", r, e, img.ext), path = plot.path)


	### Alpha ###

	ggplot(alpha.dt[MSA %in% top.cities &
		            COLLEGE_M==edu.types[[e]] & MINORITY_M==rac.types[[r]] &
		            COLLEGE_F==edu.types[[e]] & MINORITY_F==rac.types[[r]] ],
		   aes(x=AGE_M, y=AGE_F, z=VALUE)) +
		stat_contour(aes(colour=..level..), size=0.5) +
		alpha.color.grad +
		xy.line +
		age.labs +
		facet_wrap(~CITY)
	ggsave(paste0("alpha-", r, e, img.ext), path = plot.path)


	### Production ###

	ggplot(prod.dt[AGE_M < max.age-1 & AGE_F < max.age-1 & MSA %in% top.cities &
		           COLLEGE_M==edu.types[[e]] & MINORITY_M==rac.types[[r]] &
		           COLLEGE_F==edu.types[[e]] & MINORITY_F==rac.types[[r]],
				   .(CITY, AGE_M, AGE_F, TVAL = truncator(VALUE, -4, 2.4))],
		   aes(x=AGE_M, y=AGE_F, z=TVAL)) +
		#geom_tile(aes(fill = TVAL)) +
		stat_contour(aes(colour=..level..), size=0.5) +
		prod.color.grad +
		xy.line +
		age.labs +
		facet_wrap(~CITY)
	ggsave(paste0("prod-", r, e, img.ext), path = plot.path)


	### Value functions ###

	ggplot(val.dt[MSA %in% top.cities & AGE<max.age &
		          COLLEGE==edu.types[[e]] & MINORITY==rac.types[[r]] ],
		   aes(x=AGE)) +
		geom_line(aes(y=VALUE_M, color = "Male"), size=1) +
		geom_line(aes(y=VALUE_F, color = "Female"), size=1) +
		labs(x="Age", y="Value of search", color="Sex") +
		facet_wrap(~CITY)
	ggsave(paste0("val-", r, e, img.ext), path = plot.path)


	### Marriage Stocks ###
		
	ggplot(m.stock.dt[AGE_M < max.age-1 & AGE_F < max.age-1 & MSA %in% top.cities &
		              COLLEGE_M==edu.types[[e]] & MINORITY_M==rac.types[[r]] &
		              COLLEGE_F==edu.types[[e]] & MINORITY_F==rac.types[[r]] ],
		   aes(x=AGE_M, y=AGE_F, z=MASS)) +
		stat_contour(aes(colour=..level..), bins=15, size=0.5) +
		flow.color.grad +
		xy.line +
		age.labs +
		facet_wrap(~CITY)
	ggsave(paste0("mar-stock-", r, e, img.ext), path = plot.path)


	### Marriage Emigration Flows ###

	ggplot(m.mig.dt[AGE_M < max.age-1 & AGE_F < max.age-1 & MSA %in% top.cities &
		            COLLEGE_M==edu.types[[e]] & MINORITY_M==rac.types[[r]] &
		            COLLEGE_F==edu.types[[e]] & MINORITY_F==rac.types[[r]] ],
		   aes(x=AGE_M, y=AGE_F, z=NET_OUTFLOW)) +
		stat_contour(aes(colour=..level..), bins=15, size=0.5) +
		flow.color.grad +
		xy.line +
		age.labs +
		facet_wrap(~CITY)
	ggsave(paste0("mar-mig-", r, e, img.ext), path = plot.path)


	### Population Stocks ###

	ggplot(pop.dt[MSA %in% top.cities & AGE<max.age &
		          COLLEGE==edu.types[[e]] & MINORITY==rac.types[[r]] ],
		   aes(x=AGE)) +
		geom_line(aes(y=POP_M, color = "Male"), size=1) +
		geom_line(aes(y=POP_F, color = "Female"), size=1) +
		geom_line(aes(y=SNG_M, color = "Male"), size=1, linetype=2) +
		geom_line(aes(y=SNG_F, color = "Female"), size=1, linetype=2) +
		labs(x="Age", y="Counts", color="Sex") +
		facet_wrap(~CITY, scales = "free_y")
	ggsave(paste0("pop-", r, e, img.ext), path = plot.path)

}

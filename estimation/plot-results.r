# Command line usage: `Rscript plot-results.r model`
#	model: ageonly or racedu

# command line args
model = commandArgs(trailingOnly=TRUE)

# SELECT model: ageonly16 or racedu24
if (model == "ageonly") {
	model.dir <- "ageonly16"
	age.only <- TRUE
} else if (model == "racedu") {
	model.dir <- "racedu24"
	age.only <- FALSE
} else {
	model.dir <- "ageonly16"
	age.only <- TRUE
	warning(paste("Model", model, "not available, defaulting to 'ageonly'."))
}

data.path <- file.path("data/", model.dir) # input dir: smoothed pop data
estim.path <- file.path("results/estimates-csv/dynamic-const/", model.dir) # input dir: model estimates
plot.path <- file.path("results/result-plots/dynamic-const/", model.dir) # output dir: plots
img.ext <- ".png"


message("Setting up...")

library(data.table)
library(ggplot2)

source("local-regression.r") # load local-polynomial regression function

# bandwidth matrix for smoothed prod plots
prod.bw.cor <- 0.5
prod.bw <- matrix(8 * c(1, prod.bw.cor, prod.bw.cor, 1), nrow = 2, ncol = 2)
# order of local polynomial for smoothing
prod.smooth.order <- 1

# Cities for faceted plots
top.cities <- c(35620, 31080, 16980, 19100)

if (age.only) {
	min.age <- 18
	edu.types <- c(A = 1)
	rac.types <- c(A = 1)
} else {
	min.age <- 25
	edu.types <- c(N = 1, C = 2)
	rac.types <- c(W = 1, M = 2)
}

max.age <- 65
max.age.estim <- 50 # max age for plotting estimates (some strange results at later ages)


#' Clipping function for sensible production plots
truncator <- function(x, low, high) {
	x[x>high] <- high
	x[x<low] <- low
	return(x)
}


##### SETUP DATA #####

msa.dt <- fread('data/top20msa.csv') # for getting names from codes
setnames(msa.dt, old = "METCODE", new = "MSA") # rename

# joint couple objects
prod.dt <- merge(fread(file.path(estim.path, "prod.csv")), msa.dt, by = "MSA")
alpha.dt <- merge(fread(file.path(estim.path, "alpha.csv")), msa.dt, by = "MSA")
m.stock.dt <- merge(fread(file.path(data.path, "marriages.csv")), msa.dt, by = "MSA")
m.mig.dt <- merge(fread(file.path(data.path, "mar-migration.csv")), msa.dt, by = "MSA")
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
mf.dt <- merge(mf.dt, msa.dt, by = "MSA")
val.dt <- merge(val.dt, msa.dt, by = "MSA")
pop.dt <- merge(pop.dt, msa.dt, by = "MSA")
df.dt <- merge(df.dt, msa.dt, by = "MSA")

# convert categoricals to factors
val.dt[COLLEGE == 1, EDU := "No college"]
val.dt[COLLEGE == 2, EDU := "College"]
val.dt[MINORITY == 1, RAC := "Asian or White"]
val.dt[MINORITY == 2, RAC := "Minority"]

pop.dt[COLLEGE == 1, EDU := "No college"]
pop.dt[COLLEGE == 2, EDU := "College"]
pop.dt[MINORITY == 1, RAC := "Asian or White"]
pop.dt[MINORITY == 2, RAC := "Minority"]

df.dt[COLLEGE == 1, EDU := "No college"]
df.dt[COLLEGE == 2, EDU := "College"]
df.dt[MINORITY == 1, RAC := "Asian or White"]
df.dt[MINORITY == 2, RAC := "Minority"]

for (col in c('MSA', 'CITY', 'EDU', 'RAC')) set(val.dt, j = col, value = factor(val.dt[[col]]))
for (col in c('MSA', 'CITY', 'EDU', 'RAC')) set(pop.dt, j = col, value = factor(pop.dt[[col]]))
for (col in c('MSA', 'CITY', 'EDU', 'RAC')) set(df.dt, j = col, value = factor(df.dt[[col]]))

for (col in c('METCODE', 'CITY')) set(msa.dt, j = col, value = factor(msa.dt[[col]]))
set(prod.dt, j = 'CITY', value = factor(prod.dt[['CITY']]))
set(alpha.dt, j = 'CITY', value = factor(alpha.dt[['CITY']]))
set(mf.dt, j = 'CITY', value = factor(mf.dt[['CITY']]))
set(m.mig.dt, j = 'CITY', value = factor(m.mig.dt[['CITY']]))
set(m.stock.dt, j = 'CITY', value = factor(m.stock.dt[['CITY']]))

# smooth production estimates: 20 (ageonly), 20 x 2^4 = 320 (rac-edu)
prod.dt[, SMOOTH := loc.poly.reg(AGE_M, AGE_F, VALUE, prod.bw, order = prod.smooth.order),
        by = .(MSA, COLLEGE_M, MINORITY_M, COLLEGE_F, MINORITY_F)]


##### PLOTTING #####

# lists of plots (for each size) to convert to tikz
plots.full <- list() # full page square: 5.5 x 6.35
plots.rect <- list() # large rectangle: 5.5 x 3.5
plots.std <- list() # 5.5 x 2.5 standard size

# plotting layers
xy.line <- geom_segment(x = min.age, y = min.age, xend = max.age, yend = max.age, colour = "white", size = 1.2)
age.labs <- labs(x = "Husband age", y = "Wife age")

error.fill.grad <- scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0, name = "Error")
flow.color.grad <- scale_colour_gradient(low = "navy", high = "red", name = "Flow")
stock.color.grad <- scale_colour_gradient(low = "navy", high = "red", name = "Stock")
alpha.color.grad <- scale_colour_gradient2(low = "navy", mid = "red", high = "white", midpoint = 0.5, name = "Probability")
prod.color.grad <- scale_colour_gradient2(low = "navy", mid = "red", high = "white", midpoint = 0, name = "Output")
theme.colorbar <- theme(legend.position = "bottom", legend.key.width = unit(2, "cm"), legend.key.height = unit(0.3, "cm"))


message("Start plotting...")

# lifecycle prod for homogamous couples with 2 year age gap (from smoothed global average f)
if (age.only) {
	p <- ggplot(
			prod.dt[MSA %in% top.cities & AGE_F <= max.age.estim & AGE_F > min.age + 5 &
			        ((AGE_M - AGE_F) %in% c(-2, 0, 2, 4)),
			        .(AGE_F, SMOOTH, CITY, GAP = AGE_M - AGE_F)],
			aes(x = AGE_F, y = SMOOTH, color = as.factor(GAP))) +
		geom_line(size = 1) +
		labs(x = "Wife age", y = "Output", color = "Age gap") +
		facet_wrap(~CITY, scales = "free_y")

	plots.rect[["prod_smooth_lifecycle"]] <- p
	ggsave(paste0("prod-smooth-lifecycle", img.ext), path = plot.path)

	p <- ggplot(
			prod.dt[AGE_F <= max.age.estim & AGE_F > min.age + 5 &
			        ((AGE_M - AGE_F) %in% c(-2, 0, 2, 4)),
			        .(AVG = mean(SMOOTH)), # average across MSAs
			        by = .(AGE_M, AGE_F)],
			aes(x = AGE_F, y = AVG, color = as.factor(AGE_M - AGE_F))) +
		geom_line(size = 1) +
		labs(x = "Wife age", y = "Output", color = "Age gap")

	plots.std[["prod_smooth_global_lifecycle"]] <- p
	ggsave(paste0("prod-smooth-global-lifecycle", img.ext), path = plot.path)
} else {
	p <- ggplot(
			prod.dt[MSA %in% top.cities & AGE_M <= max.age.estim & AGE_F == AGE_M - 2 &
			        COLLEGE_M == COLLEGE_F & MINORITY_M == MINORITY_F],
			aes(x = AGE_F, y = SMOOTH, color = as.factor(paste0(COLLEGE_M, MINORITY_M)))) +
		geom_line(size = 1) +
		scale_color_discrete(name = "Types",
		                     breaks = c('11',    '12',   '21',   '22'),
		                     labels = c("NC,NM", "NC,M", "C,NM", "C,M")) +
		labs(x = "Wife age", y = "Output") +
		facet_wrap(~CITY, scales = "free_y")

	plots.rect[["prod_smooth_lifecycle_racedu"]] <- p
	ggsave(paste0("prod-smooth-lifecycle_racedu", img.ext), path = plot.path)

	p <- ggplot(
			prod.dt[AGE_M <= max.age.estim & AGE_F == AGE_M - 2 & COLLEGE_M == COLLEGE_F & MINORITY_M == MINORITY_F,
			        .(AVG = mean(SMOOTH)), # average across MSAs
			        by = .(AGE_M, COLLEGE_M, MINORITY_M, AGE_F, COLLEGE_F, MINORITY_F)],
			aes(x = AGE_F, y = AVG, color = as.factor(paste0(COLLEGE_M, MINORITY_M)))) +
		geom_line(size = 1) +
		scale_color_discrete(name = "Types",
		                     breaks = c('11',    '12',   '21',   '22'),
		                     labels = c("NC,NM", "NC,M", "C,NM", "C,M")) +
		labs(x = "Wife age", y = "Output")

	plots.std[["prod_smooth_global_lifecycle_racedu"]] <- p
	ggsave(paste0("prod-smooth-global-lifecycle-racedu", img.ext), path = plot.path)
}


message("Pairwise plots...")

for (r in names(rac.types)) for (e in names(edu.types)) { # r,e: indiv/wife types

	##### Individual-level plots first #####

	### Population Stocks ###

	p <- ggplot(
			pop.dt[MSA %in% top.cities & AGE < max.age &
			       COLLEGE == edu.types[[e]] & MINORITY == rac.types[[r]] ],
			aes(x = AGE)) +
		geom_line(aes(y = POP_M / 1000, color = "Male"), size = 1) +
		geom_line(aes(y = POP_F / 1000, color = "Female"), size = 1) +
		geom_line(aes(y = SNG_M / 1000, color = "Male"), size = 1, linetype = 2) +
		geom_line(aes(y = SNG_F / 1000, color = "Female"), size = 1, linetype = 2) +
		labs(x = "Age", y = "Counts (thousands)", color = "Sex") +
		theme(legend.position = "bottom") +
		facet_wrap(~CITY, scales = "free_y")

	plots.full[[paste0("pop_", r, e)]] <- p
	ggsave(paste0("pop-", r, e, img.ext), path = plot.path)


	### DF fit: data and model ###

	p <- ggplot(
			df.dt[MSA %in% top.cities & AGE < max.age &
			      COLLEGE == edu.types[[e]] & MINORITY == rac.types[[r]] ],
			aes(x = AGE)) +
		geom_line(aes(y = FLOW_M_DATA, color = "Male"), size = 1) +
		geom_line(aes(y = FLOW_F_DATA, color = "Female"), size = 1) +
		geom_line(aes(y = FLOW_M_MODEL, color = "Male"), size = 1, linetype = 2) +
		geom_line(aes(y = FLOW_F_MODEL, color = "Female"), size = 1, linetype = 2) +
		labs(x = "Age", y = "Divorces", color = "Sex") +
		theme(legend.position = "bottom") +
		facet_wrap(~CITY, scales = "free_y")

	plots.rect[[paste0("DF_", r, e)]] <- p
	ggsave(paste0("df-", r, e, img.ext), path = plot.path)


	### Value functions ###

	p <- ggplot(
			val.dt[MSA %in% top.cities & AGE <= max.age.estim &
			       COLLEGE == edu.types[[e]] & MINORITY == rac.types[[r]] ],
			aes(x = AGE)) +
		geom_line(aes(y = VALUE_M, color = "Male"), size = 1) +
		geom_line(aes(y = VALUE_F, color = "Female"), size = 1) +
		labs(x = "Age", y = "Value of search", color = "Sex") +
		theme(legend.position = "bottom") +
		facet_wrap(~CITY)

	plots.rect[[paste0("val_", r, e)]] <- p
	ggsave(paste0("val-", r, e, img.ext), path = plot.path)



	##### Pairwise plot inner loop #####
	for (r_h in names(rac.types)) for (e_h in names(edu.types)) { # r_h, e_h: husband types

		### MF fit: data and model ###

		# data: top 4
		p <- ggplot(
				mf.dt[AGE_M < max.age-1 & AGE_F < max.age-1 & MSA %in% top.cities &
				      COLLEGE_M == edu.types[[e_h]] & MINORITY_M == rac.types[[r_h]] &
				      COLLEGE_F == edu.types[[e]] & MINORITY_F == rac.types[[r]] ],
				aes(x = AGE_M, y = AGE_F, z = DATA)) +
			stat_contour(aes(colour = ..level..), bins = 15, size = 0.5) +
			flow.color.grad +
			xy.line +
			age.labs +
			theme.colorbar +
			facet_wrap(~CITY)

		plots.full[[paste0("dMF_", r_h, e_h, r, e)]] <- p
		ggsave(paste0("mf-data-", r_h, e_h, r, e, img.ext), path = plot.path)

		# model: top 4
		p <- ggplot(
				mf.dt[MODEL >= 0 & AGE_M < max.age-1 & AGE_F < max.age-1 & MSA %in% top.cities &
				      COLLEGE_M == edu.types[[e_h]] & MINORITY_M == rac.types[[r_h]] &
				      COLLEGE_F == edu.types[[e]] & MINORITY_F == rac.types[[r]] ],
				aes(x = AGE_M, y = AGE_F, z = MODEL)) +
			stat_contour(aes(colour = ..level..), bins = 15, size = 0.5) +
			flow.color.grad +
			xy.line +
			age.labs +
			theme.colorbar +
			facet_wrap(~CITY)

		plots.full[[paste0("MF_", r_h, e_h, r, e)]] <- p
		ggsave(paste0("mf-model-", r_h, e_h, r, e, img.ext), path = plot.path)

		# fit pct error, top 4: (model - data) / data, with data smoothed
		p <- ggplot(
				mf.dt[DATA > 2 & AGE_M < max.age-1 & AGE_F < max.age-1 & MSA %in% top.cities &
				      COLLEGE_M == edu.types[[e_h]] & MINORITY_M == rac.types[[r_h]] &
				      COLLEGE_F == edu.types[[e]] & MINORITY_F == rac.types[[r]],
				      .(AGE_M, AGE_F, CITY,
				        FIT_ERROR = truncator((MODEL - DATA) / DATA, -0.6, 0.6))],
				aes(x = AGE_M, y = AGE_F, z = FIT_ERROR)) +
			geom_raster(aes(fill = FIT_ERROR), interpolate = TRUE) +
			#stat_contour(aes(colour = ..level..), size = 0.5) +
			error.fill.grad +
			xy.line +
			age.labs +
			theme.colorbar +
			facet_wrap(~CITY)

		plots.full[[paste0("fitMF_", r_h, e_h, r, e)]] <- p
		ggsave(paste0("mf-fit-", r_h, e_h, r, e, img.ext), path = plot.path)


		### Alpha ###

		p <- ggplot(
				alpha.dt[AGE_M <= max.age.estim & AGE_F <= max.age.estim & MSA %in% top.cities &
				         COLLEGE_M == edu.types[[e_h]] & MINORITY_M == rac.types[[r_h]] &
				         COLLEGE_F == edu.types[[e]] & MINORITY_F == rac.types[[r]] ],
				aes(x = AGE_M, y = AGE_F, z = VALUE)) +
			stat_contour(aes(colour = ..level..), size = 0.5) +
			alpha.color.grad +
			xy.line +
			age.labs +
			theme.colorbar +
			facet_wrap(~CITY)

		plots.full[[paste0("alpha_", r_h, e_h, r, e)]] <- p
		ggsave(paste0("alpha-", r_h, e_h, r, e, img.ext), path = plot.path)


		### Production ###

		# re-usable dt, especially for partial plots
		prod.smooth.dt <- prod.dt[AGE_M <= max.age.estim & AGE_F <= max.age.estim &
		                          COLLEGE_M == edu.types[[e_h]] & MINORITY_M == rac.types[[r_h]] &
		                          COLLEGE_F == edu.types[[e]] & MINORITY_F == rac.types[[r]],
		                          .(TVAL = truncator(mean(SMOOTH), -5, 4)), # average across MSAs
		                          by = .(AGE_M, COLLEGE_M, MINORITY_M, AGE_F, COLLEGE_F, MINORITY_F)]

		# global average f
		p <- ggplot(
				prod.dt[AGE_M <= max.age.estim & AGE_F <= max.age.estim &
				        COLLEGE_M == edu.types[[e_h]] & MINORITY_M == rac.types[[r_h]] &
				        COLLEGE_F == edu.types[[e]] & MINORITY_F == rac.types[[r]],
				        .(TVAL = truncator(mean(VALUE), -4, 4)), # average across MSAs
				        by = .(AGE_M, COLLEGE_M, MINORITY_M, AGE_F, COLLEGE_F, MINORITY_F)],
				aes(x = AGE_M, y = AGE_F, z = TVAL)) +
			stat_contour(aes(colour = ..level..), size = 0.5) +
			expand_limits(x = min.age, y = min.age) + # prevent plot from cropping off
			prod.color.grad +
			xy.line +
			age.labs +
			theme.colorbar

		plots.full[[paste0("prod_global_", r_h, e_h, r, e)]] <- p
		ggsave(paste0("prod-global-", r_h, e_h, r, e, img.ext), path = plot.path)

		# smoothed global average f
		p <- ggplot(prod.smooth.dt, aes(x = AGE_M, y = AGE_F, z = TVAL)) +
			stat_contour(aes(colour = ..level..), bins = 12, size = 0.5) +
			expand_limits(x = min.age+3, y = min.age+3) + # prevent plot from cropping off
			prod.color.grad +
			xy.line +
			age.labs +
			theme.colorbar

		plots.full[[paste0("prod_smooth_global_", r_h, e_h, r, e)]] <- p
		ggsave(paste0("prod-smooth-global-", r_h, e_h, r, e, img.ext), path = plot.path)

		# partial lines plot from smoothed global average f
		p <- ggplot() +
			geom_line(data = prod.smooth.dt[AGE_M == 38 & AGE_F > min.age + 5],
			          aes(x = AGE_F, y = TVAL, color = "Female"), size = 1) +
			geom_line(data = prod.smooth.dt[AGE_F == 36 & AGE_M > min.age + 5],
			          aes(x = AGE_M, y = TVAL, color = "Male"), size = 1) +
			labs(x = "Age", y = "Output", color = "Sex") +
			theme(legend.position = "bottom")

		plots.std[[paste0("prod_smooth_global_partials_", r_h, e_h, r, e)]] <- p
		ggsave(paste0("prod-smooth-global-partials-", r_h, e_h, r, e, img.ext), path = plot.path)

		# local estimates
		p <- ggplot(
				prod.dt[AGE_M <= max.age.estim & AGE_F <= max.age.estim & MSA %in% top.cities &
				        COLLEGE_M == edu.types[[e_h]] & MINORITY_M == rac.types[[r_h]] &
				        COLLEGE_F == edu.types[[e]] & MINORITY_F == rac.types[[r]],
				        .(CITY, AGE_M, AGE_F, TVAL = truncator(VALUE, -4, 2.4))],
				aes(x = AGE_M, y = AGE_F, z = TVAL)) +
			stat_contour(aes(colour=..level..), size=0.5) +
			prod.color.grad +
			xy.line +
			age.labs +
			theme.colorbar +
			facet_wrap(~CITY)

		plots.full[[paste0("prod_", r_h, e_h, r, e)]] <- p
		ggsave(paste0("prod-", r_h, e_h, r, e, img.ext), path = plot.path)

		# smoothed local estimates
		p <- ggplot(
				prod.dt[AGE_M <= max.age.estim & AGE_F <= max.age.estim & MSA %in% top.cities &
				        COLLEGE_M == edu.types[[e_h]] & MINORITY_M == rac.types[[r_h]] &
				        COLLEGE_F == edu.types[[e]] & MINORITY_F == rac.types[[r]],
				        .(CITY, AGE_M, AGE_F, TVAL = truncator(SMOOTH, -4, 2.4))],
				aes(x = AGE_M, y = AGE_F, z = TVAL)) +
			stat_contour(aes(colour = ..level..), size = 0.5) +
			prod.color.grad +
			xy.line +
			age.labs +
			theme.colorbar +
			facet_wrap(~CITY)

		plots.full[[paste0("prod_smooth_", r_h, e_h, r, e)]] <- p
		ggsave(paste0("prod-smooth-", r_h, e_h, r, e, img.ext), path = plot.path)


		### Marriage Stocks ###

		p <- ggplot(
				m.stock.dt[AGE_M < max.age-1 & AGE_F < max.age-1 & MSA %in% top.cities &
				           COLLEGE_M == edu.types[[e_h]] & MINORITY_M == rac.types[[r_h]] &
				           COLLEGE_F == edu.types[[e]] & MINORITY_F == rac.types[[r]] ],
				aes(x = AGE_M, y = AGE_F, z = MASS)) +
			stat_contour(aes(colour = ..level..), bins = 15, size = 0.5) +
			stock.color.grad +
			xy.line +
			age.labs +
			theme.colorbar +
			facet_wrap(~CITY)

		plots.full[[paste0("marstock_", r_h, e_h, r, e)]] <- p
		ggsave(paste0("mar-stock-", r_h, e_h, r, e, img.ext), path = plot.path)


		### Marriage Emigration Flows ###

		p <- ggplot(
				m.mig.dt[AGE_M < max.age-1 & AGE_F < max.age-1 & MSA %in% top.cities &
				         COLLEGE_M == edu.types[[e_h]] & MINORITY_M == rac.types[[r_h]] &
				         COLLEGE_F == edu.types[[e]] & MINORITY_F == rac.types[[r]] ],
				aes(x = AGE_M, y = AGE_F, z = NET_OUTFLOW)) +
			stat_contour(aes(colour = ..level..), bins = 15, size = 0.5) +
			flow.color.grad +
			xy.line +
			age.labs +
			theme.colorbar +
			facet_wrap(~CITY)

		plots.full[[paste0("marmig_", r_h, e_h, r, e)]] <- p
		ggsave(paste0("mar-mig-", r_h, e_h, r, e, img.ext), path = plot.path)

	}
}

message("Done!")

### Save plot objects ###

# save objects for further manipulation and conversion to tex
plots <- list(full = plots.full,
              rect = plots.rect,
              std = plots.std)
saveRDS(plots, file = file.path(plot.path, "plot-objects.Rdata"))


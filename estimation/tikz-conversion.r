# Command line usage: `Rscript tikz-conversion.R model`
#	model: ageonly or racedu

# command line args
model = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggplot2)
library(tikzDevice)

### Load and set up data ###

message("Loading plot objects...")

plot.path <- "results/result-plots/dynamic-const/"
if (model == "ageonly") {
	plot.dir <- "ageonly16"
} else if (model == "racedu") {
	plot.dir <- "racedu24"
} else {
	plot.dir <- "ageonly16"
	warning(paste("Model", model, "not available, defaulting to 'ageonly'."))
}

plots <- readRDS(file = file.path(plot.path, plot.dir, "plot-objects.Rdata"))
out.dir <- file.path(plot.path, plot.dir, "tikz")


### Full page plots ###

message("Building full page tikz plots...")

for (name in names(plots$full)) {
	tikz(file = file.path(out.dir, paste(name, "tex", sep = ".")),
	     width = 5.5, height = 6.35)
	print(plots$full[[name]]) # doesn't work without print()
	dev.off()
}

### Large size plots ###

message("Building large tikz plots...")

for (name in names(plots$rect)) {
	tikz(file = file.path(out.dir, paste(name, "tex", sep = ".")),
	     width = 5.5, height = 3.5)
	print(plots$rect[[name]]) # doesn't work without print()
	dev.off()
}

### Standard size plots ###

message("Building standard tikz plots...")

for (name in names(plots$std)) {
	tikz(file = file.path(out.dir, paste(name, "tex", sep = ".")),
	     width = 5.5, height = 2.5)
	print(plots$std[[name]]) # doesn't work without print()
	dev.off()
}

message("Done!")


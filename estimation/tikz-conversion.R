# Command line usage: `Rscript tikz-conversion.R out/dir/ plots.Rdata`
# plots.Rdata: filename of a list of ggplot objects saved as RDS

# command line args
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggplot2)
library(tikzDevice)

### Load and set up data ###

out.dir <- args[1]
plots <- readRDS(file = args[2])

### Standard size plots ###

for (name in names(plots)) {
	tikz(file = file.path(out.dir, paste(name, "tex", sep = ".")),
		 width = 5.5, height = 2.5)
	print(plots[[name]]) # doesn't work without print()
	dev.off()
}

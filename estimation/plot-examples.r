library(data.table)
library(ggplot2)
library(tikzDevice)
library(grid)
library(gridExtra)

source("local-regression.r") # reuse bivariate normal function

out.dir <- "results/result-plots/examples"
img.ext <- ".png"

### Example production plots ###

grid.pts <- seq(from = 0, to = 10, length.out = 30)
dt <- CJ(a = grid.pts, b = grid.pts) # table of grid points
H <- 16 * matrix(c(1,.9,.9,1), nrow = 2, ncol = 2) # covariance matrix
num.bin <- 8

# reusable a=b line
xy.line <- geom_segment(x = 0, y = 0, xend = 10, yend = 10, colour = "red", linetype = 2, size = 1)

# a=b symmetry
p.left <- ggplot(dt[, .(a, b, Z = K.H(a, b, H))],
		aes(x=a, y=b, z=Z)) +
	stat_contour(aes(colour=..level..), bins=num.bin, size=0.5) +
	xy.line + expand_limits(x = 10, y = 10) +
	# no legend or axis ticks
	guides(colour=F) +
	theme(axis.line = element_blank(), axis.ticks = element_blank(),
	      axis.text.x = element_blank(), axis.text.y = element_blank(),
	      axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave(paste0("symmetric-prod", img.ext), path = out.dir)

# a=b+k symmetry
k <- 2
p.mid <- ggplot(dt[, .(a, b, Z = K.H(a, b + k, H))],
		aes(x=a, y=b, z=Z)) +
	stat_contour(aes(colour=..level..), bins=num.bin, size=0.5) +
	xy.line + expand_limits(x = 10, y = 10) +
	geom_segment(x = k, y = 0, xend = 10, yend = 10 - k, colour = "red", size = 0.8) +
	# no legend or axis ticks
	guides(colour=F) +
	theme(axis.line = element_blank(), axis.ticks = element_blank(),
	      axis.text.x = element_blank(), axis.text.y = element_blank(),
	      axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave(paste0("shifted-prod", img.ext), path = out.dir)

# a=b+k symmetry
q <- 1.5
p.right <- ggplot(dt[, .(a, b, Z = K.H(a, b * q, H))],
		aes(x=a, y=b, z=Z)) +
	stat_contour(aes(colour=..level..), bins=num.bin, size=0.5) +
	xy.line + expand_limits(x = 10, y = 10) +
	geom_segment(x = 0, y = 0, xend = 10, yend = 10 / q, colour = "red", size = 0.8) +
	# no legend or axis ticks
	guides(colour=F) +
	theme(axis.line = element_blank(), axis.ticks = element_blank(),
	      axis.text.x = element_blank(), axis.text.y = element_blank(),
	      axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave(paste0("scaled-prod", img.ext), path = out.dir)


tikz(file = file.path(out.dir, paste("prod_examples", "tex", sep = ".")),
     width = 6, height = 2.2)
grid.arrange(p.left, p.mid, p.right,
             ncol = 3,
             bottom = "Husband age", left = "Wife age")
dev.off()

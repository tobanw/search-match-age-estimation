## Hand-rolled bivariate local-polynomial regression functions with unconstrained bandwidth matrix
# allows for diagonal orientation of kernel, to smooth more along joint aging diagonal

#' Inverse function for positive definite matrices
matrix.inv <- function(A) {chol2inv(chol(A))}

#' Matrix sqrt function
matrix.sqrt <- function(A) {
    if (length(A) == 1)
        return(sqrt(A))
    sva <- svd(A)
    if (min(sva$d) >= 0)
        Asqrt <- sva$u %*% diag(sqrt(sva$d)) %*% t(sva$v)
    else
        stop("Matrix square root is not defined")
  return(Asqrt)
}

#' Standard bivariate normal density function
K <- function(x) {
	# x: nx2 matrix
	return(1.0/(2*pi) * exp(-0.5 * rowSums(x^2)))
}

#' Bivariate kernel with unrestricted bandwidth matrix
K.H <- function(x, y, H) {
	# x,y: vectors
    Hinv12 <- matrix.sqrt(matrix.inv(H))
    z <- cbind(x, y)
    return(1.0/sqrt(det(H)) * K(z %*% t(Hinv12)))
}

#' Local-polynomial bivariate regression
#'
#' Arguments:
#'	x,y (vector): bivariate domain vectors
#'	z (vector): outcome value vector, to be smoothed
#'	H (matrix): 2x2 bandwidth matrix
#'	order: order of polynomial (1-3) used in local regression
#'
#' Returns: smoothed value vector corresponding to z
loc.poly.reg <- function(x, y, z, H, order = 3) {
    n <- length(z)
    smth <- vector(length = n) # empty vector for smoothed z values
    
    # for each data point, computed a smoothed alternate (assuming a complete grid of data)
    for (i in 1:n) {
		xm <- x - x[i]
		ym <- y - y[i]
        X <- cbind(1, xm, ym) # local-linear
		if (order > 1) { # augment to local-quadratic
			X <- cbind(X, xm^2, ym^2, xm*ym)
		}
		if (order > 2) { # augment to local-cubic
			X <- cbind(X, xm^3, ym^3, xm^2*ym, xm*ym^2)
		}
        W <- K.H(xm, ym, H) # vector of kernel weights
        tXW <- t(X * W) # element-wise multiply each col of X by W

		# smoothed value is the intercept of local regression: kernel-weighted least squares
        theta.i <- solve(tXW %*% X, tXW %*% z) # solve the normal equation
        smth[i] <- theta.i[1]
    }
    return(smth)
}

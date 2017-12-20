rnormregmixcov <- function (n, alpha, mubeta, sigma, x = NULL) {
  # Generates mixed normal random variables with regressor x
  # Input
  #  n : number of observations
  #  alpha  : m-vector
  #  mubeta  : k by m matrix
  #  sigma  : m-vector
  #  x : (n by k-1) matrix NOT including a constant
  # Output
  #  y : n by 1 vector
  m     <- length(alpha)
  mubeta   <- matrix(mubeta, ncol=m)

  if (!is.null(x)){
    x <- as.matrix(x)
    if (nrow(x) != n) { stop("y and x must have the same number of rows.") }
    x1   <- cbind(1,x)
    ii   <- sample(m, n, replace=TRUE, prob=alpha)
    y   <- rnorm(n, mean = rowSums(x1*t(mubeta[, ii])), sd = sigma[ii])
  } else {
    ii   <- sample(m, n, replace=TRUE, prob=alpha)
    y   <- rnorm(n, mean = mubeta[, ii], sd = sigma[ii])
  }

  y

}  # end function rnormregmixcov

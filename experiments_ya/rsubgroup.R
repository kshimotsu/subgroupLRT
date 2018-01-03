rsubgroup <- function (n, tau, mubeta, sigma, x, v) {
  # Generates mixed normal random variables with regressor x and v
  # Input
  #  n : number of observations
  #  tau  : q-vector
  #  mubeta  : k by m matrix
  #  sigma  : m-vector
  #  x : (n by k-1) matrix NOT including a constant
  #  v : (n by q-1) matrix Not including a constant
  # Output
  #  y : n by 1 vector
  m <- ncol(mubeta)
  q <- length(tau)
  mubeta   <- matrix(mubeta, ncol=m)
  tau <- matrix(tau, nrow = q)

  v <- as.matrix(v)
  x <- as.matrix(x)
  v1 <- cbind(1, v)
  x1 <- cbind(1, x)

  v1tau <- v1 %*% tau
  alpha1 <- exp(v1tau) / (1 + exp(v1tau))
  alpha <- cbind(alpha1, 1 - alpha1)
  ii <- apply(alpha, 1, function(x) {sample(m, 1, replace = TRUE, prob = x)})

  x1mubeta <- rowSums(x1 * t(mubeta[, ii]))
  y <- rnorm(n = n, mean = x1mubeta, sd=sigma[ii])

  y

}  # end function rnormregmixcov

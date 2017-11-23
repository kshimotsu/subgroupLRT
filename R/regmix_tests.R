#' @description Performs the LRT statistic given the data for y and x on
#' the null hypothesis H_0: m = m_0. Using this function is equivalent to
#' calling normalmixLRT with regressors specified by x as a parameter.
#' @export
#' @title regmixLRT_homo
#' @name regmixLRT_homo
#' @param y n by 1 vector of data for y.
#' @param x n by q matrix of data for x.
#' @param m number of components in the mixture defined by a null hypothesis, m_0.
#' @param z n by p matrix of regressor associated with gamma.
#' @param tauset set of initial tau value candidates.
#' @param an tuning parameter used in the penalty function.
#' @param ninits number of randomly drawn initial values.
#' @param crit.method method used to compute the critical values, one of \code{"none"},
#' \code{"asy"}, and \code{"boot"}. The default option is \code{"asy"}. When \code{method = "asy"},
#' the p-values are computed using the asymptotic critical values. When \code{method = "boot"},
#' the p-values are computed by bootstrap.
#' @param nbtsp The number of bootstrap replicates. Default is 199.
#' @param cl cluster used for parallelization; if it is \code{NULL}, the system
#' will automatically generate a new one for computation accordingly.
#' @param parallel Determines what percentage of available cores are used,
#' represented by a double in [0,1]. Default is 1.
#' @return A list of class \code{normalregMix} with items:
#' \item{lrtstat}{the value of the LRT statistic.}
#' \item{pvals}{p-value.}
#' \item{crit}{When \code{crit.method = "asy"}, \code{crit} is vector of critical values at the 0.1, 0.05, 0.01 level.
#'  When \code{crit.method = "boot"}, \code{crit} is 3 by 3 matrix of critival values at the (0.1, 0.05, 0.01) level, jth row corresponding to k=j.}
#' \item{crit.method}{method used to compute the variance-covariance matrix.}
#' \item{parlist}{parameter estimates as a list containing alpha, mubeta, and sigma (and gam if z is included in the model).}
#' \item{call}{The matched call.}
#' \item{m}{number of components in the mixture defined by the null hypothesis, \eqn{m_0}.}
#' @examples
#' data(faithful)
#' attach(faithful)
#' \dontrun{regmixLRT_homo(y = eruptions, x = waiting, m = 1, crit.method = "asy")}
#' \dontrun{regmixLRT_homo(y = eruptions, x = waiting, m = 2, crit.method = "asy")}
regmixLRT_homo <- function (y, x, m = 2, z = NULL, tauset = c(0.1,0.3,0.5),
                           an = NULL, ninits = 100,
                           crit.method = c("none", "asy", "boot"), nbtsp = 199,
                           cl = NULL,
                           parallel = 1) {

  y <- as.vector(y)
  x <- as.matrix(x)
  q <- ncol(x)
  n <- length(y)

  if (!is.null(z))
    z <- as.matrix(z)

  regmix.pmle.result    <- regmixPMLE_homo(y=y, x=x, m=m, z=z, vcov.method="none", ninits=ninits)
  loglik0 <- regmix.pmle.result$loglik

  regmix.pmle.result1    <- regmixPMLE_homo(y=y, x=x, m=m+1, z=z, vcov.method="none", ninits=ninits)
  loglik1 <- regmix.pmle.result1$loglik

  lrtstat <- 2*(loglik1-loglik0)

  if (crit.method == "asy"){
    result  <- regmixCrit(y=y, x=x, parlist=regmix.pmle.result$parlist, z=z, values=lrtstat,
                          parallel=parallel, cl=cl, nrep=1000, ninits.crit=25)
  } else if (crit.method == "boot") {
    result  <- regmixCritBoot_homo(y=y, x=x, parlist=regmix.pmle.result$parlist, z=z, values=lrtstat,
                              ninits=ninits, nbtsp=nbtsp, parallel=parallel, cl=cl)
  } else {
    result <- list()
    result$crit <- result$pvals <- NA
  }

  a <- list(lrtstat = lrtstat, pvals = result$pvals, crit = result$crit, parlist = regmix.pmle.result$parlist,
            call = match.call(), m = m, crit.method = crit.method, nbtsp = nbtsp,
            label = "MEMtest")

  class(a) <- "normalregMix"

  a

}  # end regmixLRT_homo

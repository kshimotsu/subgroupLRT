#' @description Performs the LRT statistic given the data for y and x on
#' the null hypothesis H_0: m = m_0. Using this function is equivalent to
#' calling normalmixLRT with regressors specified by x as a parameter.
#' @export
#' @title subgroupLRT_homo_2
#' @name subgroupLRT_homo_2
#' @param y n by 1 vector of data for y.
#' @param x n by q1 matrix of data for x.
#' @param v n by q2 matrix of data for v
#' @param m number of components in the mixture defined by a null hypothesis, m_0.
#' @param z n by p matrix of regressor associated with gamma.
#' @param ninits number of randomly drawn initial values.
#' @param crit.method method used to compute the critical values, one of \code{"none"},
#' \code{"asy"}, and \code{"boot"}. The default option is \code{"none"}. When \code{method = "asy"},
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
#' \item{crit}{vector of critical values at the 0.1, 0.05, 0.01 level.}
#' \item{crit.method}{method used to compute the variance-covariance matrix.}
#' \item{parlist}{parameter estimates as a list containing alpha, mubeta, and sigma (and gam if z is included in the model).}
#' \item{call}{The matched call.}
#' \item{m}{number of components in the mixture defined by the null hypothesis, \eqn{m_0}.}
#' @examples
#' data(faithful)
#' attach(faithful)
#' \dontrun{regmixLRT_homo(y = eruptions, x = waiting, m = 1, crit.method = "none")}
#' \dontrun{regmixLRT_homo(y = eruptions, x = waiting, m = 2, crit.method = "none")}
subgroupLRT_homo_2 <- function (y, x, v, m = 1, z = NULL, ninits = 10, k = 100,
                             crit.method = c("none", "boot"), nbtsp = 199,
                             cl = NULL, parallel = 1) {

  y <- as.vector(y)
  x <- as.matrix(x)
  v <- as.matrix(v)
  n <- length(y)

  if (!is.null(z))
    z <- as.matrix(z)

  subgroup.mle.result    <- subgroupMLE_homo(y=y, x=x, v=v, m=m, z=z, vcov.method="none", ninits=ninits)
  loglik0 <- subgroup.mle.result$loglik

  ninits2 = floor(ninits/3)
  subgroup.mle.result1    <- subgroupMLE_homo_2(y=y, x=x, v=v, m=m+1, z=z, k = k,
                                               tauinit=c(1,-2), vcov.method="none", ninits=ninits2)
  loglik1 <- subgroup.mle.result1$loglik

  subgroup.mle.result2    <- subgroupMLE_homo_2(y=y, x=x, v=v, m=m+1, z=z, k = k,
                                               tauinit=c(1,2), vcov.method="none", ninits=ninits2)
  loglik2 <- subgroup.mle.result2$loglik

  subgroup.mle.result3    <- subgroupMLE_homo_2(y=y, x=x, v=v, m=m+1, z=z, k = k,
                                                tauinit=c(1,0), vcov.method="none", ninits=ninits2)
  loglik3 <- subgroup.mle.result3$loglik

    loglikem <- max(c(loglik1, loglik2, loglik3))

  lrtstat <- 2*(loglikem - loglik0)

  if (crit.method == "boot") {
    result  <-subgroupCritBoot_homo_2(y=y, x=x, v=v, parlist=subgroup.mle.result$parlist, z=z, values=lrtstat,
                                       k=k, ninits=ninits, nbtsp=nbtsp, parallel=parallel, cl=cl)
  } else {
    result <- list()
    result$crit <- result$pvals <- NA
  }

  a <- list(lrtstat = lrtstat, pvals = result$pvals, crit = result$crit, parlist = subgroup.mle.result$parlist,
            call = match.call(), m = m, crit.method = crit.method, nbtsp = nbtsp)
  # ,label = "MEMtest")

  # class(a) <- "normalregMix"

  a

}  # end subgroupLRT_homo_2


#' @description Performs the LRT statistic given the data for y and x on
#' the null hypothesis H_0: m = m_0. Using this function is equivalent to
#' calling normalmixLRT with regressors specified by x as a parameter.
#' @export
#' @title subgroupEM_homo
#' @name subgroupEM_homo
#' @param y n by 1 vector of data for y.
#' @param x n by q1 matrix of data for x.
#' @param v n by q2 matrix of data for v
#' @param m number of components in the mixture defined by a null hypothesis, m_0.
#' @param z n by p matrix of regressor associated with gamma.
#' @param ninits number of randomly drawn initial values.
#' @param crit.method method used to compute the critical values, one of \code{"none"},
#' \code{"asy"}, and \code{"boot"}. The default option is \code{"none"}. When \code{method = "asy"},
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
#' \item{crit}{vector of critical values at the 0.1, 0.05, 0.01 level.}
#' \item{crit.method}{method used to compute the variance-covariance matrix.}
#' \item{parlist}{parameter estimates as a list containing alpha, mubeta, and sigma (and gam if z is included in the model).}
#' \item{call}{The matched call.}
#' \item{m}{number of components in the mixture defined by the null hypothesis, \eqn{m_0}.}
#' @examples
#' data(faithful)
#' attach(faithful)
#' \dontrun{regmixLRT_homo(y = eruptions, x = waiting, m = 1, crit.method = "none")}
#' \dontrun{regmixLRT_homo(y = eruptions, x = waiting, m = 2, crit.method = "none")}
subgroupEM_homo <- function (y, x, v, m = 1, z = NULL, ninits = 10, k = 9,
                              crit.method = c("none", "boot"), nbtsp = 199,
                              cl = NULL, parallel = 1) {

  y <- as.vector(y)
  x <- as.matrix(x)
  v <- as.matrix(v)
  n <- length(y)

  if (!is.null(z))
    z <- as.matrix(z)

  subgroup.mle.result    <- subgroupMLE_homo(y=y, x=x, v=v, m=m, z=z, vcov.method="none", ninits=ninits)
  loglik0 <- subgroup.mle.result$loglik

  ninits2 = floor(ninits/2)
  subgroup.em.result1    <- subgroupMLE_homo_EM(y=y, x=x, v=v, m=m+1, z=z, k = k, tauinit=c(1,-2),
                                              vcov.method="none", ninits=ninits2)
  loglik1 <- subgroup.em.result1$loglik

  subgroup.em.result2    <- subgroupMLE_homo_EM(y=y, x=x, v=v, m=m+1, z=z, k = k, tauinit=c(1,2),
                                              vcov.method="none", ninits=ninits2)
  loglik2 <- subgroup.em.result2$loglik
  loglikem <- max(c(loglik1, loglik2))

  lrtstat <- 2*(loglikem - loglik0)

  if (crit.method == "boot") {
    result  <-subgroupCritBoot_homo_EM(y=y, x=x, v=v, parlist=subgroup.mle.result$parlist, z=z, values=lrtstat,
                                    k=k, ninits=ninits, nbtsp=nbtsp, parallel=parallel, cl=cl)
  } else {
    result <- list()
    result$crit <- result$pvals <- NA
  }

  a <- list(lrtstat = lrtstat, pvals = result$pvals, crit = result$crit, parlist = subgroup.mle.result$parlist,
            call = match.call(), m = m, crit.method = crit.method, nbtsp = nbtsp)
  # ,label = "MEMtest")

  # class(a) <- "normalregMix"

  a

}  # end subgroupEM_homo


#' @description Performs the LRT statistic given the data for y and x on
#' the null hypothesis H_0: m = m_0. Using this function is equivalent to
#' calling normalmixLRT with regressors specified by x as a parameter.
#' @export
#' @title subgroupLRT_homo
#' @name subgroupLRT_homo
#' @param y n by 1 vector of data for y.
#' @param x n by q1 matrix of data for x.
#' @param v n by q2 matrix of data for v
#' @param m number of components in the mixture defined by a null hypothesis, m_0.
#' @param z n by p matrix of regressor associated with gamma.
#' @param ninits number of randomly drawn initial values.
#' @param crit.method method used to compute the critical values, one of \code{"none"},
#' \code{"asy"}, and \code{"boot"}. The default option is \code{"none"}. When \code{method = "asy"},
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
#' \item{crit}{vector of critical values at the 0.1, 0.05, 0.01 level.}
#' \item{crit.method}{method used to compute the variance-covariance matrix.}
#' \item{parlist}{parameter estimates as a list containing alpha, mubeta, and sigma (and gam if z is included in the model).}
#' \item{call}{The matched call.}
#' \item{m}{number of components in the mixture defined by the null hypothesis, \eqn{m_0}.}
#' @examples
#' data(faithful)
#' attach(faithful)
#' \dontrun{regmixLRT_homo(y = eruptions, x = waiting, m = 1, crit.method = "none")}
#' \dontrun{regmixLRT_homo(y = eruptions, x = waiting, m = 2, crit.method = "none")}
subgroupLRT_homo <- function (y, x, v, m = 1, z = NULL, ninits = 10,
                           crit.method = c("none", "boot"), nbtsp = 199,
                           tauinit = c(0,0), cl = NULL, parallel = 1) {

  y <- as.vector(y)
  x <- as.matrix(x)
  v <- as.matrix(v)
  n <- length(y)

  if (!is.null(z))
    z <- as.matrix(z)

  subgroup.mle.result    <- subgroupMLE_homo(y=y, x=x, v=v, m=m, z=z, vcov.method="none", ninits=ninits)
  loglik0 <- subgroup.mle.result$loglik

  subgroup.mle.result1    <- subgroupMLE_homo(y=y, x=x, v=v, m=m+1, z=z,
                                              tauinit = tauinit,
                                              vcov.method="none", ninits=ninits)
  loglik1 <- subgroup.mle.result1$loglik

  lrtstat <- 2*(loglik1-loglik0)

  if (crit.method == "boot") {
    result  <-subgroupCritBoot_homo(y=y, x=x, v=v, parlist=subgroup.mle.result$parlist, z=z, values=lrtstat,
                              tauinit = tauinit, ninits=ninits, nbtsp=nbtsp, parallel=parallel, cl=cl)
  } else {
    result <- list()
    result$crit <- result$pvals <- NA
  }

  a <- list(lrtstat = lrtstat, pvals = result$pvals, crit = result$crit, parlist = subgroup.mle.result$parlist,
            call = match.call(), m = m, crit.method = crit.method, nbtsp = nbtsp)
  # ,label = "MEMtest")

  # class(a) <- "normalregMix"

  a

}  # end subgroupLRT_homo

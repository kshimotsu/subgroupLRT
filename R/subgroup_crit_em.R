#' @description Computes the bootstrap critical values of the LRT.
#' @export
#' @title subgroupCritBoot_homo_EM
#' @name subgroupCritBoot_homo_EM
#' @param y n by 1 vector of data for y.
#' @param x n by q vector of data for x.
#' @param parlist parameter estimates as a list containing alpha, mubeta, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mubeta = (mubeta_1', ..., mubeta_m'),
#' sigma = (sigma_1, ..., sigma_m), gam).
#' @param z n by p matrix of regressor associated with gamma.
#' @param values vector the values of the LRT statistic at which the p-values are computed.
#' @param ninits number of candidates of the initial value of the EM algorithm.
#' @param nbtsp number of bootstrap replicates. Default is 199.
#' @param parallel Determines what percentage of available cores are used, represented by a double in [0,1]. Default is 1.
#' @param cl cluster used for parallelization; if it is \code{NULL}, the system will automatically generate a cluster.
#' @return A list with the following items:
#' \item{crit}{vector of critical values at the 0.1, 0.05, 0.01 level.}
#' \item{pvals}{p-values corresponding to \code{values}.}
subgroupCritBoot_homo_EM <- function (y, x, v, parlist, z = NULL, values = NULL, ninits = 100, k,
                                      nbtsp = 199, parallel = 1, cl = NULL) {
  # if (normalregMixtest.env$normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMixtest.env$normalregMix.test.seed)

  y  <- as.vector(y)
  n  <- length(y)
  x  <- as.matrix(x)
  v <- as.matrix(v)

  tau   <- parlist$tau
  mubeta  <- parlist$mubeta
  sigma   <- parlist$sigma
  gam     <- parlist$gam
  m <- ncol(mubeta)
  sigma   <- rep(sigma,m)

  pvals <- NULL

  # Generate bootstrap observations
  ybset <- replicate(nbtsp, rsubgroup(n = n, tau = tau, mubeta = mubeta, sigma = sigma, x = x, v = v))

  if (!is.null(z)) {
    zgam <- as.matrix(z) %*% gam
    ybset <- ybset + replicate(nbtsp, as.vector(zgam))
  }

  num.cores <- max(1,floor(detectCores()*parallel))
  if (num.cores > 1) {
    if (is.null(cl)) {
      cl <- makeCluster(num.cores)
      newcluster <- TRUE
    }
    clusterSetRNGStream(cl, 8888577)
    out <- parLapplyLB(cl, 1:nbtsp, function(j) subgroupEM_homo(y=ybset[,j], x = x, v = v, k = k,
                                                                m = m, z = z, parallel = 0, ninits = ninits, crit.method = "none"))
    if (newcluster) {
      on.exit(stopCluster(cl))
    } else {
      on.exit(cl)
    }
  }
  else
    out <- apply(ybset, 2, subgroupEM_homo, x = x, v = v, k = k,
                 m = m, z = z, parallel = 0, ninits = ninits, crit.method = "none")

  lrtstat.b <- sapply(out, "[[", "lrtstat")  # 1 by nbstp vector

  lrtstat.b <- sort(lrtstat.b)

  q <- ceiling(nbtsp*c(0.90,0.95,0.99))
  crit <- lrtstat.b[q]

  if (!is.null(values)) {
    k <- length(values)
    pvals <- rowMeans(t(matrix(rep.int(lrtstat.b,k),ncol=k)) > values)
  }

  return(list(crit = crit, pvals = pvals))
}  # end function subgroupCritBoot_homo_EM


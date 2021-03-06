#' @description A variant of subgroup MLE. The objective funtion is maximized
#' while fixing tau. Once the objective function is maximized, then update
#' the parameters k times by EM step.
#' @export
#' @title subgroupMLE_homo_2
#' @name subgroupMLE_homo_2
#' @param y n by 1 vector of data for y.
#' @param x n by q1 matrix of data for x.
#' @param v n by q2 matrix of data for w
#' @param m number of components in the mixture.
#' @param z n by p matrix of regressor associated with gam.
#' @param vcov.method Method used to compute the variance-covariance matrix, one of \code{"Hessian"} and \code{"OPG"}.
#' The default option is \code{"Hessian"}. When \code{method = "Hessian"}, the variance-covarince matrix is
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009).
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @param ninits number of randomly drawn initial values.
#' @param epsilon convergence criterion. Convergence is declared when the log-likelihood increases by less than \code{epsilon}.
#' @param maxit maximum number of iterations.
#' @param epsilon.short convergence criterion in short EM. Convergence is declared when the log-likelihood increases by less than \code{epsilon.short}.
#' @param maxit.short maximum number of iterations in short EM.
#' @param binit initial value of parameter vector that is included as a candidate parameter vector.
#' @return　A list of class \code{normalregMix} with items:
#' \item{coefficients}{vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\beta_{1.1},
#' \ldots,\beta_{q.1},\ldots,\mu_m,\beta_{1.m},\ldots,\beta_{q.m},\sigma,gam}.}
#' \item{parlist}{parameter estimates as a list containing alpha, mubeta, and sigma (and gam if z is included in the model).}
#' \item{vcov}{estimated variance-covariance matrix.}
#' \item{loglik}{maximized value of the log-likelihood.}
#' \item{aic}{Akaike Information Criterion of the fitted model.}
#' \item{bic}{Bayesian Information Criterion of the fitted model.}
#' \item{postprobs}{n by m matrix of posterior probabilities for observations}
#' \item{components}{n by 1 vector of integers that indicates the indices of components
#' each observation belongs to based on computed posterior probabilities}
#' \item{call}{The matched call.}
#' \item{m}{number of components in the mixture.}
#' @note \code{subgroupMLE_homo} maximizes the log-likelihood function
#' using the EM algorithm with combining short and long runs of EM steps as in Biernacki et al. (2003).
#' \code{subgroupMLE_homo} first runs the EM algorithm from \code{ninits}\eqn{* 4m(1 + p)} initial values
#' with the convertence criterion \code{epsilon.short} and \code{maxit.short}.
#' Then, \code{subgroup_homo} uses \code{ninits} best initial values to run the EM algorithm
#' with the convertence criterion \code{epsilon} and \code{maxit}.
#' @references     Biernacki, C., Celeux, G. and Govaert, G. (2003)
#' Choosing Starting Values for the EM Algorithm for Getting the
#' Highest Likelihood in Multivariate Gaussian Mixture Models,
#' \emph{Computational Statistics and Data Analysis}, \bold{41}, 561--575.
#'
#' Boldea, O. and Magnus, J. R. (2009)
#' Maximum Likelihood Estimation of the Multivariate Normal Mixture Model,
#' \emph{Journal of the American Statistical Association},
#' \bold{104}, 1539--1549.
#'
#' Chen, J., Tan, X. and Zhang, R. (2008)
#' Inference for Normal Mixtures in Mean and Variance,
#' \emph{Statistica Sinica}, \bold{18}, 443--465.
#'
#' McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley \& Sons, Inc.
#' @examples
#' data(faithful)
#' attach(faithful)
#' regmixMLE_homo(y = eruptions, x = waiting, m = 1)
#' regmixMLE_homo(y = eruptions, x = waiting, m = 2)
#' # out <- regmixMLE_homo(y = eruptions, x = waiting, m = 2)
#' # summary(out)
subgroupMLE_homo_2 <- function (y, x, v, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", "none"),
                              ninits = 10, k = 100, tauinit, epsilon = 1e-08, maxit = 2000,
                              epsilon.short = 1e-02, maxit.short = 500, binit = NULL) {

  y <- as.vector(y)
  x <- as.matrix(x)  # n by (q1-1) matrix, x corresponds to Z in Shen and He
  v <- as.matrix(v)  # n by (q2-1) matrix, v corresponds to X in Shen and He, tau coef
  n   <- length(y)
  if (nrow(x) != n) { stop("y and x must have the same number of rows.") }
  x1 <- cbind(1, x)
  v1 <- cbind(1, v)
  q1 <- ncol(x1)
  q2 <- ncol(v1)

  p       <- 0
  gam   <- NULL
  ninits.short <- ninits*10*(q1+p)*m
  vcov.method <- match.arg(vcov.method)
  vcov.method = "none" # Currenlty, vcov matrix calculation is disabled.

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    if (nrow(z) != n) { stop("y and z must have the same number of rows.") }
  }

  npar    <- q1*m + 1 + p + q2  # number of parameters
  xz      <- cbind(x, z)
  ls.out  <- lsfit(xz, y)
  sd0     <- sqrt(mean(ls.out$residuals^2))

  if (m == 1) {
    mubeta <- as.matrix(unname(ls.out$coeff[1:q1]))
    if (!is.null(z)) {gam <- unname(ls.out$coeff[(q1+1):(q1+p)])}
    res     <- ls.out$residuals
    sigma   <- sqrt(mean(res*res))
    loglik  <- - (n/2)*(1 + log(2*pi) + 2*log(sigma))
    aic     <- -2*loglik + 2*npar
    bic     <- -2*loglik + log(n)*npar

    parlist <- list(tau = 0, mubeta = mubeta, sigma = sigma, gam = gam)
    coefficients <- c(tau = 0, mubeta = mubeta, sigma = sigma, gam = gam)
    postprobs <- rep(1, n)

  } else {  # m >= 2

    # generate initial values
    tmp <- subgroupMLEinit_homo(y = y, x = x, v = v, z = z, ninits = ninits.short, m = m)
    tmp$tau <- matrix(tauinit, nrow = q2, ncol = ninits.short)

    if (is.null(z))
      ztilde <- matrix(0) # dummy
    else
      ztilde <- z

    # short EM
    b0 <- rbind(tmp$tau, tmp$mubeta, tmp$sigma, tmp$gam)
    if (!is.null(binit)) {
      b0[ , 1] <- binit
    }
    out.short <- cppSubgroupMLE_homo(b0, y, x, v, ztilde, m, p, maxit.short,
                                     updatetau = 0, emtest = 0,
                                     ninits.short, epsilon.short)
    # long EM
    components <- order(out.short$loglikset, decreasing = TRUE)[1:ninits]
    b1 <- b0[ ,components] # b0 has been updated
    out0 <- cppSubgroupMLE_homo(b1, y, x, v, ztilde, m, p, maxit,
                               updatetau = 0, emtest = 0, ninits, epsilon)

    index <- which.max(out0$loglikset)
    b11 <- b1[,index]
    b11 <- matrix(b11, ncol=1)
    # tau <- b1[1:q2,index] # b0 has been updated
    # mubeta <- matrix(b1[(1+q2):(q2+q1*m),index],nrow=q1,ncol=m)
    # sigma <- b1[(q2+q1*m+1),index]
    # if (!is.null(z)) {
    #   gam     <- b1[((q1+1)*m+2):((q1+1)*m+p+1),index]
    # }
    # loglik    <- out$loglikset[index]
    # postprobs <- matrix(out$post[,index], nrow=n)
    # if (is.nan(loglik)) {
    #   loglik <- -Inf
    # }

    # update parameters k+1 times
    out <- cppSubgroupMLE_homo(b11, y, x, v, ztilde, m, p, maxit = k,
                               updatetau = 1, emtest = 0, ninits=1, epsilon)
    index <- which.max(out$loglikset)
    tau <- b11[1:q2] # b1 has been updated
    mubeta <- matrix(b11[(1+q2):(q2+q1*m)],nrow=q1,ncol=m)
    sigma <- b11[(q2+q1*m+1)]
    if (!is.null(z)) {
      gam     <- b11[((q1+1)*m+2):((q1+1)*m+p+1)]
    }
    loglik    <- out$loglikset
    postprobs <- matrix(out$post, nrow=n)
    if (is.nan(loglik)) {
      loglik <- -Inf
    }

    aic <- -2*loglik + 2*npar
    bic <- -2*loglik + log(n)*npar

    tau1.sign <- sign(tau[1])
    tau <- tau1.sign*tau;
    tau1.positive <- (tau1.sign+1)/2;
    mu.order <-  tau1.positive*c(1:2) + (1-tau1.positive)*c(2:1)
    # This works only when m=2

    #    mu.order  <- order(mubeta[1,])
    mubeta    <- mubeta[,mu.order]

    postprobs <- postprobs[, mu.order]

    colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))

    mubeta.name <- matrix(0,nrow = q1, ncol = m)
    mubeta.name[1,] <- paste("mu", 1:m, sep = "")

    if (q1 == 2) {
      mubeta.name[2,] <- paste("beta", 1:m,  sep = "")
    } else {
      for (i in 1:(q1-1)) {
        for (j in 1:m) {
          mubeta.name[i+1,j] <- paste("beta", i, ".", j, sep = "")
        }
      }
    }

    parlist <- list(tau = tau, mubeta = mubeta, sigma = sigma, gam = gam)
    coefficients <- unlist(parlist)
    names(coefficients)[(q2+1):(q2+q1*m)] <- c(mubeta.name)
  }  # end m >= 2

  # if (vcov.method == "none") {
  #   vcov <- NULL
  # } else {
  #   vcov <- regmixVcov(y = y, x = x, coefficients = coefficients, z = z , vcov.method = vcov.method)
  # }
  #
  # a <- list(coefficients = coefficients, parlist = parlist, vcov = vcov, loglik = loglik,
  #           aic = aic, bic = bic, postprobs = postprobs,
  #           components = getComponentcomponents(postprobs),
  #           call = match.call(), m = m, label = "PMLE")

  a <- list(coefficients = coefficients, parlist = parlist, loglik = loglik,
            aic = aic, bic = bic)
  # class(a) <- "normalregMix"

  a

}  # end function subgroupMLE_homo_EM

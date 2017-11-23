#' @description Estimates parameters of a homoscedastic finite mixture
#' of regressions by the method of maximum likelhood.
#' @export
#' @title regmixPMLE_homo
#' @name regmixPMLE_homo
#' @param y n by 1 vector of data for y.
#' @param x n by q matrix of data for x.
#' @param m number of components in the mixture.
#' @param z n by p matrix of regressor associated with gam.
#' @param vcov.method Method used to compute the variance-covariance matrix, one of \code{"Hessian"} and \code{"OPG"}.
#' The default option is \code{"Hessian"}. When \code{method = "Hessian"}, the variance-covarince matrix is
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009).
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @param ninits number of randomly drawn initial values.
#' @param epsilon convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit maximum number of iterations.
#' @param epsilon.short convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
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
#' @note \code{regmixPMLE_homo} maximizes the penalized log-likelihood function
#' using the EM algorithm with combining short and long runs of EM steps as in Biernacki et al. (2003).
#' \code{regmixPMLE_homo} first runs the EM algorithm from \code{ninits}\eqn{* 4m(1 + p)} initial values
#' with the convertence criterion \code{epsilon.short} and \code{maxit.short}.
#' Then, \code{regmixPMLE_homo} uses \code{ninits} best initial values to run the EM algorithm
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
#' regmixPMLE_homo(y = eruptions, x = waiting, m = 1)
#' regmixPMLE_homo(y = eruptions, x = waiting, m = 2)
#' out <- regmixPMLE_homo(y = eruptions, x = waiting, m = 2)
#' summary(out)
regmixPMLE_homo <- function (y, x, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", "none"),
                        ninits = 10, epsilon = 1e-08, maxit = 2000,
                        epsilon.short = 1e-02, maxit.short = 500, binit = NULL) {

  y   <- as.vector(y)
  x   <- as.matrix(x)   # n by (q1-1) matrix
  n   <- length(y)
  if (nrow(x) != n) { stop("y and x must have the same number of rows.") }
  x1  <- cbind(1, x)
  q1   <- ncol(x1)

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

  npar    <- m-1 + (q1+1)*m + p  # number of parameters
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
    penloglik <- loglik

    parlist <- list(alpha = 1, mubeta = mubeta, sigma = sigma, gam = gam)
    coefficients <- c(alpha = 1, mubeta = mubeta, sigma = sigma, gam = gam)
    postprobs <- rep(1, n)

  } else {  # m >= 2

    # generate initial values
    tmp <- regmixPMLEinit_homo(y = y, x = x, z = z, ninits = ninits.short, m = m)

    h       <- 0    # setting h=0 gives PMLE
    tau     <- 0.5  # setting tau=0.5 gives PMLE
    k <- 0 # setting k=0 gives PMLE

    sigma0  <- rep(sd0, m)
    mu0     <- double(m)    # dummy
    an      <- 1/n  # penalty term for variance

	if (is.null(z))
      ztilde <- matrix(0) # dummy
	else
      ztilde <- z

    # short EM
    b0 <- rbind( tmp$alpha, tmp$mubeta, tmp$sigma, tmp$gam )
    if (!is.null(binit)) {
      b0[ , 1] <- binit
    }
    out.short <- cppRegmixPMLE_homo(b0, y, x, ztilde, mu0, sigma0, m, p, an, maxit.short,
                                  ninits.short, epsilon.short)
    # long EM
    components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
    b1 <- b0[ ,components] # b0 has been updated
    out <- cppRegmixPMLE_homo(b1, y, x, ztilde, mu0, sigma0, m, p, an, maxit, ninits, epsilon)

    index     <- which.max(out$penloglikset)
    alpha <- b1[1:m,index] # b0 has been updated
    mubeta <- matrix(b1[(1+m):((q1+1)*m),index],nrow=q1,ncol=m)
    sigma <- b1[(1+(q1+1)*m):((q1+2)*m),index]
    if (!is.null(z)) {
      gam     <- b1[((q1+2)*m+1):((q1+2)*m+p),index]
    }
    penloglik <- out$penloglikset[index]
    loglik    <- out$loglikset[index]
    postprobs <- matrix(out$post[,index], nrow=n)

    aic <- -2*loglik + 2*npar
    bic <- -2*loglik + log(n)*npar

    mu.order  <- order(mubeta[1,])
    alpha     <- alpha[mu.order]
    mubeta    <- mubeta[,mu.order]
    sigma     <- sigma[mu.order]

    postprobs <- postprobs[, mu.order]
    colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))

    mubeta.name <- matrix(0,nrow = q1, ncol = m)
    mubeta.name[1,] <- paste("mu", 1:m, sep = "")

    if (q1 == 2) {
      mubeta.name[2,] <- paste("beta", 1:m,  sep = "")
    } else {
      for (i in 1:(q1-1)) {
        for (j in 1:m) {
          # mubeta.name[i+1,j] <- paste("beta", j, i, sep = "")
          mubeta.name[i+1,j] <- paste("beta", i, ".", j, sep = "")
        }
      }
    }

    parlist <- list(alpha = alpha, mubeta = mubeta, sigma = sigma, gam = gam)
    coefficients <- unlist(parlist)
    names(coefficients)[(m+1):((q1+1)*m)] <- c(mubeta.name)
  }  # end m >= 2

  if (vcov.method == "none") {
    vcov <- NULL
  } else {
    vcov <- regmixVcov(y = y, x = x, coefficients = coefficients, z = z , vcov.method = vcov.method)
  }

  a <- list(coefficients = coefficients, parlist = parlist, vcov = vcov, loglik = loglik,
            penloglik = penloglik, aic = aic, bic = bic, postprobs = postprobs,
            components = getComponentcomponents(postprobs),
            call = match.call(), m = m, label = "PMLE")

  class(a) <- "normalregMix"

  a

}  # end function regmixPMLE_homo

#' Generate initial values used by \code{regmixPMLE_homo}.
#' @export
#' @title regmixPMLEinit_homo
#' @name regmixPMLEinit_homo
#' @param y n by 1 vector of data for y.
#' @param x n by q matrix of data for x.
#' @param z n by p matrix of regressor associated with gamma.
#' @param ninits number of initial values to be generated.
#' @param m number of components in the mixture.
#' @return A list with the following items:
#' \item{alpha}{m by ninits matrix for alpha.}
#' \item{mubeta}{q+1 by m by ninits array for mu and beta.}
#' \item{sigma}{m by ninits matrix for sigma.}
#' \item{gam}{m by ninits matrix for gam.}
regmixPMLEinit_homo <- function (y, x, z = NULL, ninits = 1, m = 2)
{
  if (normalregMixtest.env$normalregMix.test.on) # initial values controlled by normalregMix.test.on
    set.seed(normalregMixtest.env$normalregMix.test.seed)

  n  <- length(y)
  q1  <- ncol(x)+1
  p  <- ncol(z)

  gam <- NULL
  if (!is.null(z)) {
    out     <- lsfit(cbind(x, z), y)
    gam0  <- out$coef[(q1+1):(q1+p)]
    gam   <- matrix(runif(p*ninits, min=0.5, max=1.5), nrow=p)*gam0
    mubeta_hat <- out$coef[1:q1]
    y     <- y - z %*% gam0
    r     <- out$residuals
    stdR  <- sd(r)
  } else {
    out         <- lsfit(x, y)
    mubeta_hat  <- out$coef
    r           <- out$residuals
    stdR        <- sd(r)
  }

  alpha <- matrix(runif(m*ninits), nrow=m)
  alpha <- t(t(alpha)/colSums(alpha))


  minMU <- min(y - x %*% mubeta_hat[-1])
  maxMU <- max(y - x %*% mubeta_hat[-1])
  mubeta <- matrix(0, nrow=q1*m, ncol=ninits)
  for (j in 1:m) {
    mubeta[(q1*(j-1)+1), ] <- runif(ninits, min=minMU, max=maxMU)
    for (i in 2:q1) {
      mubeta[(q1*(j-1)+i), ] <- mubeta_hat[i]*runif(ninits, min=-2, max=2)
    }
  }
  sigma <- matrix(runif(m*ninits, min=0.01, max=1), nrow=m)*stdR

  list(alpha = alpha, mubeta = mubeta, sigma = sigma, gam = gam)

}  # end function regmixPMLEinit_homo

#' Computes the variance-covariance matrix of the MLE of m-component regression mixture.
#' @export
#' @title regmixVcov_homo
#' @name regmixVcov_homo
#' @param y n by 1 vector of data for y.
#' @param x n by q matrix of data for x.
#' @param coefficients vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\beta_{1.1},
#' \ldots,\beta_{q.1},\ldots,\mu_m,\beta_{1.m},\ldots,\beta_{q.m},\sigma_1,\ldots,\sigma_m,gam}.
#' @param z n by p matrix of regressor associated with gamma.
#' @param vcov.method Method used to compute the variance-covariance matrix,
#' one of \code{"Hessian"} and \code{"OPG"}. #' The default option is \code{"Hessian"}.
#' When \code{method = "Hessian"}, the variance-covarince matrix is
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009).
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @return variance-covariance matrix of the MLE of
#' m-component normal mixture given the data and coefficients.
#' @references   Boldea, O. and Magnus, J. R. (2009)
#' Maximum Likelihood Estimation of the Multivariate Normal Mixture Model,
#' \emph{Journal of the American Statistical Association},
#' \bold{104}, 1539--1549.
regmixVcov_homo <- function(y, x, coefficients, z = NULL, vcov.method = c("Hessian", "OPG")) {
  # Computes the variance-covariance matrix of the MLE of m-component normal regression mixture
  # Input
  #  y  : n by 1 vector of dependent variable
  #  x  : n by q1-1 matrix of regressor NOT including an intercept
  #  coefficients : (alpha_1,...,alpha_m,mubeta_1^T, ...,mubeta_m^T,sigma_1, ..., sigma_m,gam^T)
  #  z  : n by p matrix of regressor associated with gam
  # Output
  #  vcov: variance-covariance matrix

  y     <- as.vector(y)
  n     <- length(y)
  len   <- length(coefficients)
  p     <- 0
  gam <- NULL
  vcov.method <- match.arg(vcov.method)

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    gam <- coefficients[(len-p+1):len]
  }

  x   <- as.matrix(x)
  x1  <- cbind(1, x)
  q1   <- ncol(x1)
  q   <- ncol(x)

  m  <- (len-p)/(3+q)
  if (round(m) != m)
    stop("The dimension of the coefficients is incompatible with x and z. Please check the data.")

  alpha   <- coefficients[1:m]  # m-vector
  mubeta  <- matrix(coefficients[(m+1):((2+q)*m)], nrow=q+1, ncol=m)  # q+1 by m
  sigma   <- coefficients[((2+q)*m+1):((3+q)*m)]  # m-vector

  if (m == 1) {
    xz1 <- cbind(x1,z)
    I <- matrix(0, nrow=q1+p+1, ncol=q1+p+1)
    I[1:(q1+p), 1:(q1+p)] <- t(xz1) %*% xz1/sigma^2
    I[(q1+p+1), (q1+p+1)] <- n*2/sigma^2
    if (p != 0){
      s.1 <- c(1:q1, (q1+p+1), (q1+1):(q1+p))
      I <- I[s.1, s.1]
    }

    vcov <- solve(I)
    # Because the variance is parameterized as sigma^2, we convert it to sigma
    c.mat.vec <- c(rep(1,q1),(1/sigma^(1/2))/2,rep(1,p))
    vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)

  } else {  # end of if (m == 1)
    # m >= 2
    # Compute posterior probabilities, and adjust y if z is present
    sigma0  <- rep(1, m)  # dummy
    mu0     <- double(m)  # dummy
    an      <- 1/n  # penalty term for variance
    h       <- 0
    tau     <- 0.5
    k <- 0
    epsilon <- 1e-08
	maxit = 2
    ninits = 1
    b <- matrix( rep( coefficients, ninits), ncol = ninits)

    if (is.null(z))
      out.p <- cppRegmixPMLE(b, y, x, matrix(0), mu0, sigma0, m, p, an, maxit, ninits, epsilon, tau, h, k)
    else
	  {
      out.p <- cppRegmixPMLE(b, y, x, z, mu0, sigma0, m, p, an, maxit, ninits, epsilon, tau, h, k)
      # Adjust y
      y <- as.vector(y - z %*% gam)
    }

    post <- matrix(out.p$post, nrow=n)

    p2 <- seq(q1+1, (q1+1)*m, by=q1+1)  # sequence of q1+1, (q1+1)*2, ... , (q1+1)*m
    p1 <- (1:((q1+1)*m))[-p2]        # other values from 1, ..., (q1+1)*m

    a <- diag(1/alpha[-m], nrow=m-1, ncol=m-1)
    a <- cbind(a, -1/alpha[m])  # m-1 by m matrix of a_i's
    abar <- a %*% t(post)  # m-1 by n

    xtheta <- x1 %*% mubeta  # n by m

    Z0 <- t(t(y-xtheta)/sigma)          # normalized data, n by m
    f <- t(t(exp(-Z0^2/2)/sqrt(2*pi))/sigma)  # pdf, n by m
    phi <- t(t(f)*alpha)            # n by m
    f0 <- rowSums(phi)              # data pdf, n by 1

    vinv <- 1/(sigma*sigma)  # m-vector

    b <- t(t(Z0)/sigma)  # n by m
    B <- t(vinv - t(b*b))  # n by m

    c0 <- array(0,dim=c(n, m, q1+1))
    c0[, , (1:q1)] <- array(tKR(x1, b), dim=c(n, m, q1))
    c0[, , (q1+1)] <- -B/2

    # Compute Hessian-based I
    if (vcov.method == "Hessian") {
      other.method = "OPG"
      C0 <- array(0, dim=c(n, m, q1+1, q1+1))
      x11 <- array(tKR(x1, x1), dim = c(n, q1, q1))
      for (i in 1:m) {
        C0[, i, (1:q1), (1:q1)] <- x11*vinv[i]
        C0[, i, (1:q1), q1+1]   <- C0[, i, q1+1, (1:q1)] <- x1*b[, i]*vinv[i] # n by q1
      }
      C0[, , q1+1, q1+1] <- t((vinv - 2*t(B))*vinv)/2      # n by m

      Q.pi <- - abar %*% t(abar)  # m-1 by m-1

      Q.pi.theta <- matrix(0,nrow=m-1,ncol=(q1+1)*m)  # m-1 by (q1+1)*m
      for (i in 1:m) {
        zi <- a[, i] - abar  # m-1 by n
        wi <- c0[, i, ]*post[, i]  # n by q1+1
        Q.i <- colSums(tKR(wi, t(zi)))  # (q1+1)*(m-1) vector
        # first q1*(m-1) elements correspond to mubeta x pi,
        # last m-1 elements correspond to sigma x pi,
        Q.pi.theta[,(q1*(i-1)+1):(q1*i)] <- matrix(Q.i[1:(q1*(m-1))],ncol=q1)  # m-1 by q1 matrix
        Q.pi.theta[, q1*m+i] <- Q.i[(q1*(m-1)+1):((q1+1)*(m-1))]  # m-1 vector
      }

      Q.theta <- matrix(0, nrow=(q1+1)*m, ncol=(q1+1)*m)
      for (i in 2:m) {  # off-diagonal blocks
        for (j in 1:(i-1)) {
          wi  <- c0[, i, ]*post[, i] # n by q1+1
          wj  <- c0[, j, ]*post[, j] # n by q1+1
          Q.ij <- - colSums(tKR(wi, wj))  # (q1+1)*(q1+1) vector
          Q.theta[((q1+1)*(i-1)+1):((q1+1)*i), ((q1+1)*(j-1)+1):((q1+1)*j)] = t(matrix(Q.ij, nrow=q1+1, ncol=q1+1))
        }
      }

      Q.theta <- Q.theta + t(Q.theta)
      for (i in 1:m) {  # diagonal blocks
        C.ii   <- array(C0[, i, , ], dim=c(n, q1+1, q1+1))
        Q.ii.1   <- apply(C.ii*post[,i], c(2, 3), sum)
        w.ii   <- tKR(c0[, i, ], c0[, i, ])*post[, i]*(1-post[, i])
        Q.ii.2   <- matrix(colSums(w.ii), nrow=q1+1, ncol=q1+1)
        Q.theta[((q1+1)*(i-1)+1):((q1+1)*i), ((q1+1)*(i-1)+1):((q1+1)*i)] <- -Q.ii.1 + Q.ii.2
      }

      # q1+1,2*(q1+1),...,m*(q1+1)th rows and columns = sigma
      # other rows and columns = mubeta
      Q.theta <- Q.theta[c(p1, p2), c(p1, p2)]  # first block = wrt mubeta, second blosk = wrt sigma

      dimI <- m-1+(q1+1)*m
      I <- matrix(0, nrow=dimI, ncol=dimI)
      I[1:(m-1), 1:(m-1)] <- - Q.pi
      I[1:(m-1), m:dimI]  <- - Q.pi.theta
      I[m:dimI, 1:(m-1)]  <- - t(Q.pi.theta)
      I[m:dimI, m:dimI]   <- - Q.theta

      if (!is.null(z)) {
        dbar <-  z*rowSums(post*b)  # n by p
        Q.gam.theta <- matrix(0, nrow=p, ncol=(q1+1)*m)  # p by (q1+1)*m matrix
        for (i in 1:m) {
          C.i <- array(C0[, i, 1, ], dim=c(n, q1+1))  # n by q1+1
          Q.i.1 <- colSums(tKR(-C.i+b[, i]*c0[, i, ], z*post[, i])) # p*(q1+1) vector
          Q.i.2 <- colSums(tKR(c0[, i, ]*post[, i], dbar))  # p*(q1+1) vector
          Q.gam.theta[, ((q1+1)*(i-1)+1):((q1+1)*i)] <- matrix(Q.i.1+Q.i.2, nrow=p, ncol=q1+1)
        }

        Q.gam.theta <- Q.gam.theta[, c(p1, p2), drop=FALSE]  # p by (q1+1)*m
        w1 <- (post*b)%*%t(a) - rowSums(post*b)*t(abar)  # n by m-1
        Q.pi.gam.0 <- colSums(tKR(w1, z))  # (m-1)*p vector
        Q.pi.gam  <- matrix(Q.pi.gam.0, nrow=m-1, ncol=p)
        Q.gam     <- - t(z)%*%(z*rowSums(post*B)) -
          matrix(colSums(tKR(dbar, dbar)), nrow=p, ncol=p)

        I <- cbind(I, -rbind(Q.pi.gam, t(Q.gam.theta)))
        I <- rbind(I, -cbind(t(Q.pi.gam), Q.gam.theta, Q.gam))
      }  # end if (!is.null(z))

    }  else {  # compute I with (method == "OPG")
      other.method = "Hessian"
      c0.a <- array(0, dim=c(n, m, 2))
      c0.a[, , 1] <- b  # n by m
      c0.a[, , 2] <- -B/2  # n by m

      score <- t(abar)

      for (j in 1:m) {
        # score.o <- cbind(score.o, c0[, j, ]*post[, j])
        score <- cbind(score, x1*c0.a[, j, 1]*post[, j], c0.a[, j, 2]*post[, j])
        # print(all.equal(score.o, score))
      }

      ind <- c(1:(m-1), p1+m-1, p2+m-1)
      score <- score[, ind]
      I <- t(score) %*% score

      if (!is.null(z))  {
        dbar <-  z*rowSums(post*b)  # n by p
        score <- cbind(score, dbar)
        I <- t(score) %*% score
      }

    }  # end if (method=="OPG")

    vcov <- try(solve(I))
    if (class(vcov) == "try-error" || any(diag(vcov) <0) ) {
      vcov <- matrix(NaN, nrow = (2+q1)*m-1+p, ncol = (2+q1)*m-1+p)
      warning("Fisher information matrix is singular and/or the
              variance is estimated to be negative. Consider using vcov.method=\"",other.method,"\".")
    }

    # Because the variance is parameterized as sigma^2, we convert it to sigma

    c.mat.vec <- c(rep(1, m-1+m*q1), (1/sigma^(1/2))/2, rep(1, p))
    vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)
    # vcov.opg <- diag(c.mat.vec) %*% vcov.opg %*% diag(c.mat.vec)

    # Add the variance of alpha_m
    M.mat <- diag(len-1)
    M.mat <- rbind(M.mat[1:(m-1),], c(rep(-1,m-1),rep(0,len-m)), M.mat[m:(len-1),])

    vcov <- M.mat %*% vcov %*% t(M.mat)
    # vcov.opg <- M.mat %*% vcov.opg %*% t(M.mat)

  }   # end else (i.e., m >= 2)

  vcov

}  # end function regmixVcov_homo


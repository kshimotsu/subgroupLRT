subgroupMLE_nloptr <- function (y, x, v, m = 2, z = NULL,
                              ninits = 10, epsilon = 1e-08, maxit = 2000) {
  y <- as.vector(y)
  x <- as.matrix(x)  # n by (q1-1) matrix, x corresponds to Z in Shen and He
  v <- as.matrix(v)  # n by (q2-1) matrix, v corresponds to X in Shen and He, tau coef
  n   <- length(y)

  x1 <- cbind(1, x)
  v1 <- cbind(1, v)
  q1 <- ncol(x1)
  q2 <- ncol(v1)

  lb <- c(-5, -5, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 0)
  ub <- c(5, 5, Inf, Inf, Inf, Inf, Inf, Inf, Inf)

  p       <- 0
  gam   <- NULL

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
    tmp <- subgroupMLEinit_homo(y = y, x = x, v = v, z = z, ninits = ninits, m = m)
    b0 <- rbind(tmp$tau, tmp$mubeta, tmp$sigma, tmp$gam)

    # long EM

    S <- list()
    for (j in 1:ninits){
      S[[j]] <- slsqp(x0 = b0[,j], fn = negLogLikelihood, lower = lb, upper = ub,
                 control = list(maxeval = maxit, ftol_abs = epsilon),
                 y = y, x = x, v = v)
    }

    loglikset <- - sapply(S,"[[","value")
    index     <- which.max(loglikset)

    b1 <-S[[index]]$par
    tau <- b1[1:q2] # b0 has been updated
    mubeta <- matrix(b1[(1+q2):(q2+q1*m)],nrow=q1,ncol=m)
    sigma <- b1[(q2+q1*m+1)]

    loglik    <- loglikset[index]
    # postprobs <- matrix(out$post[,index], nrow=n)
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

    # postprobs <- postprobs[, mu.order]

    # colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))

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

}  # end function subgroupMLE_homo

negLogLikelihood <- function (theta, y, x, v) {

  y <- as.vector(y)
  x <- as.matrix(x)  # n by (q1-1) matrix, x corresponds to Z in Shen and He
  v <- as.matrix(v)  # n by (q2-1) matrix, v corresponds to X in Shen and He, tau coef
  n   <- length(y)
  x1 <- cbind(1, x)
  v1 <- cbind(1, v)

  tau <- theta[1:2]
  mubeta <- theta[3:8]
  sigma <- theta[9]

  alpha <- exp(v1 %*% tau)/(1+exp(v1 %*% tau))
  m1 <- (y - x1 %*% mubeta[1:3])/sigma
  m2 <- (y - x1 %*% mubeta[4:6])/sigma

  f <- (2*pi)^(-1/2) * (alpha * exp(-m1*m1/2) + (1-alpha) * exp(-m2*m2/2)) / sigma

  ll <- -sum(log(f))
  ll
}

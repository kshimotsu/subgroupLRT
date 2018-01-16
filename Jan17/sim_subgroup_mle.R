# Computes the rejection frequencies of the LRT

rm(list=ls())

library(parallel)
library(normalregMix)
library(subgroupLRT)

setwd("~/Dropbox/subgroupLRT/Jan17")

# outfilename <- "sim_regmix_cov_size.RData"

ncores <- detectCores()
cl <- makePSOCKcluster(rep('localhost',ncores),master='localhost')
clusterEvalQ(cl, library(normalregMix))
clusterEvalQ(cl, library(subgroupLRT))

# sink("sim_subgroup_mle.out", append=T)

nrep <- 15
# nobset <- c(100, 200)
nobset <- c(100)
nnob <- length(nobset)
DGPset <- c(1:1)

m <- 2

a1 <- 0.5
b1 <- 0
c1 <- 1

# mubeta <- matrix(c(1,2,1,2), nrow=2, ncol=2) # one-component model in effect
# mubeta <- matrix(c(1,2,0,0), nrow=2, ncol=2) # two-component model

# mubeta <- matrix(c(1,0,2,1,0,2), nrow=3, ncol=2) # one-component model in effect
mubeta <- matrix(c(1,0,2,2,a1,2+b1), nrow=3, ncol=2) # two-component model
tautrue <- c(1, c1)
tauvec = as.matrix(tautrue, nrow=2)
sigma <- c(0.5, 0.5)


for (DGP in DGPset)
{
  for (inob in 1:nnob)
  {
    time_start <- proc.time()
    nob <- nobset[inob]
    set.seed(123456)

    x.all <- matrix(rnorm(nob*nrep), nrow = nob) - 1
    # covariate in the mean equation, corresponds to Z in the paper
    tvec.all <- matrix(as.numeric(runif(nob*nrep)>0.5), nrow = nob)
    # covariate in the mean equation, corresponds to Z in the paper
    # v.all <- matrix(rnorm(nob*nrep), nrow = nob) + 1
    # covariate in logit equation, corresponds to X in the paper
    xt.all <- array(c(x.all,tvec.all), dim = c(nob,nrep,2))
    y.all <- matrix(double(nob*nrep), nrow = nob)
    for (j in 1:nrep) {
      # xj <- cbind(tvec.all[,j],x.all[,j])
      y.all[, j] <- rsubgroup(n = nob, tau = tauvec, mubeta = mubeta,
                              sigma = sigma, x = xt.all[,j,], v = x.all[,j])
    }

    clusterExport(cl,varlist=c("y.all", "m", "x.all", "xt.all", "tautrue"))
    # clusterExport(cl,varlist=c("y.all", "m", "x.all", "v.all"))
    clusterSetRNGStream(cl, 123456)

    mleout <- parLapplyLB(cl,1:nrep, function(j) subgroupMLE_homo(y=y.all[,j],
                                          tauinit = tautrue, x = xt.all[,j,], v = x.all[,j], m=m))
    coefsum <- t(sapply(mleout,"[[","coefficients"))
    logliksum <- sapply(mleout,"[[","loglik")
    print(colMeans(coefsum))

    emout <- parLapplyLB(cl,1:nrep, function(j) subgroupMLE_homo_EM(y=y.all[,j],
                                                        x = xt.all[,j,], v = x.all[,j], m=m,
                                                        k = 9, tauinit=c(1,2)))
    coefsum.em <- t(sapply(emout,"[[","coefficients"))
    logliksum.em <- sapply(emout,"[[","loglik")

    lrtout <- parLapplyLB(cl,1:nrep, function(j) subgroupLRT_homo(y=y.all[,j],
                                        x = xt.all[,j,], v = x.all[,j], m=1, crit.method = "none"))

    lrtsum <- sapply(lrtout,"[[","lrtstat")

    emtout <- parLapplyLB(cl,1:nrep, function(j) subgroupEM_homo(y=y.all[,j],
                                        x = xt.all[,j,], v = x.all[,j], m = 1, crit.method = "none"))

    emtsum <- sapply(emtout,"[[","lrtstat")
    # lr <- cbind(loglik0sum,logliksum)

  } # end of inob loop

} # end of DGP loop

stopCluster(cl)

# sink()

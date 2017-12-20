# Computes the rejection frequencies of the LRT

rm(list=ls())

library(parallel)
library(normalregMix)
library(subgroupLRT)

setwd("~/Dropbox/yoichi/ProR/subgroupLRT/experiments")

nrep <- 10
# nobset <- c(100,200)
nobset <- c(100)
nnob <- length(nobset)
DGPset <- c(1:1)

m <- 2

mubeta <- matrix(c(1,2,1,2), nrow=2, ncol=2) # one-component model in effect
gammatrue <- c(0.3, -0.7)
gammavec = as.matrix(gammatrue, nrow=2)
sigma <- c(0.5, 0.5)

# rejfreq5all <- matrix(0, nrow=length(DGPset), ncol=3*nnob)
# rejfreq1all <- matrix(0, nrow=length(DGPset), ncol=3*nnob)

for (DGP in DGPset)
{
  for (inob in 1:nnob)
  {
    time_start <- proc.time()
    nob <- nobset[inob]
    set.seed(123456)

    x.all <- matrix(rnorm(nob*nrep), nrow = nob) - 1 # corresponds to Z in the paper
    w.all <- matrix(rnorm(nob*nrep), nrow = nob) + 1 # corresponds to X in the paper
    y.all <- matrix(double(nob*nrep), nrow = nob)
    for (j in 1:nrep) {
      w1 <- cbind(1, w.all[, j])
      w1gamma <- w1 %*% gammavec
      alpha1 <- exp(w1gamma) / (1 + exp(w1gamma))
      alpha <- cbind(alpha1, 1 - alpha1)
      ii <- apply(alpha, 1, function(x) {sample(m, 1, replace = TRUE, x)})
      x1 <- cbind(1, x.all[, j])
      x1mubeta <- rowSums(x1 * t(mubeta[, ii]))
      y.all[, j] <- rnorm(n = nob, mean = mubeta, sd=sigma[ii])
    }

    #     lrtout <- parLapply(cl,1:nrep, function(j) regmixLRT_homo(y=y.all[,j],
    #                               x=x.all[,j], m=m, crit.method="boot", parallel=0) )
    # 		pvalsum <- sapply(lrtout,"[[","pvals")
    # 		# print(pvalsum)
    # 		rejfreq10 <- 100*mean(pvalsum < 0.10)
    # 		rejfreq5 <- 100*mean(pvalsum < 0.05)
    # 		rejfreq1 <- 100*mean(pvalsum < 0.01)
    # 		# rejfreq5all[DGP,(3*inob-2):(3*inob)] <- rejfreq5
    # 		# rejfreq1all[DGP,(3*inob-2):(3*inob)] <- rejfreq1
    #
    # 		time_end <- proc.time()
    # 		runtime  <- time_end - time_start
    #
    # 		print("DGP, nob, nrep")
    # 		print(c(DGP, nob, nrep))
    # 		print(runtime)
    # 		print("rejfreq10, rejfreq5, rejfreq1")
    # 		print(c(rejfreq10, rejfreq5, rejfreq1))
    # 		rm(list=c("y.all", "x.all"))
    # 		save.image(file = outfilename)
  } # end of inob loop

} # end of DGP loop

# return(list(outall=outall,runtime=runtime))


# Computes the rejection frequencies of the LRT

rm(list=ls())

library(parallel)
library(normalregMix)
library(subgroupLRT)

setwd("~/Dropbox/subgroupLRT/Jan17")

ncores <- detectCores()
cl <- makePSOCKcluster(rep('localhost',ncores),master='localhost')
clusterEvalQ(cl, library(normalregMix))
clusterEvalQ(cl, library(subgroupLRT))

nrep <- 1000
nob <- 100
DGPset <- c(1:1)

m <- 1
k <- 100

mubeta <- matrix(c(1,0,2,1,0,2), nrow=3, ncol=2) # one-component model in effect
tautrue <- c(1, 0)
tauvec = as.matrix(tautrue, nrow=2)
sigma <- c(0.5, 0.5)

for (DGP in DGPset)
{
	time_start <- proc.time()
	set.seed(123456)

	x.all <- matrix(rnorm(nob*nrep), nrow = nob) - 1
	# covariate in the mean equation, corresponds to Z in the paper
	tvec.all <- matrix(as.numeric(runif(nob*nrep)>0.5), nrow = nob)
	# covariate in the mean equation, corresponds to Z in the paper
	# v.all <- matrix(rnorm(nob*nrep), nrow = nob) + 1
	# covariate in logit equation, corresponds to X in the paper
	xt.all <- array(c(x.all,tvec.all), dim = c(nob,nrep,2))
	y.all <- matrix(double(nob*nrep), nrow = nob)
	for (j in 1:nrep){
	  y.all[, j] <- rsubgroup(n = nob, tau = tauvec, mubeta = mubeta,
	                          sigma = sigma, x = xt.all[,j,], v = x.all[,j])
	}

	clusterExport(cl,varlist=c("y.all", "m", "xt.all", "x.all", "k"))
	clusterSetRNGStream(cl, 123456)

  lrtout <- parLapplyLB(cl,1:nrep, function(j) subgroupLRT_homo_2(y=y.all[,j],
                            x=xt.all[,j,], v = x.all[,j], m = m , k = k, crit.method="boot",
                            parallel=0) )
	pvalsum <- sapply(lrtout,"[[","pvals")
	# print(pvalsum)
	rejfreq10 <- 100*mean(pvalsum < 0.10)
	rejfreq5 <- 100*mean(pvalsum < 0.05)
	rejfreq1 <- 100*mean(pvalsum < 0.01)

	sink("sim_lrt2_size.out", append=T)

	time_end <- proc.time()
	runtime  <- time_end - time_start
	print(runtime)

	print("DGP, nob, nrep, k")
	print(c(DGP, nob, nrep, k))
	print(mubeta)
	print(tauvec)

	print("rejfreq10, rejfreq5, rejfreq1")
	print(c(rejfreq10, rejfreq5, rejfreq1))

} # end of DGP loop

stopCluster(cl)

sink()

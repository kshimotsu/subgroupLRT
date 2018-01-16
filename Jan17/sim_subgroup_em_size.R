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

nrep <- 100
nob <- 100
DGPset <- c(1:1)

m <- 1

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

	clusterExport(cl,varlist=c("y.all", "m", "xt.all", "x.all"))
	clusterSetRNGStream(cl, 123456)

	emout <- parLapplyLB(cl,1:nrep, function(j) subgroupEM_homo(y=y.all[,j],
	                          k = 9, x=xt.all[,j,], v = x.all[,j], m=m, crit.method="boot",
	                          parallel=0) )
	pvalsum.em <- sapply(emout,"[[","pvals")
	# print(pvalsum)
	rejfreq10.em <- 100*mean(pvalsum.em < 0.10)
	rejfreq5.em <- 100*mean(pvalsum.em < 0.05)
	rejfreq1.em <- 100*mean(pvalsum.em < 0.01)

	sink("sim_em_size.out", append=T)

	time_end <- proc.time()
	runtime  <- time_end - time_start
	print(runtime)

	print("DGP, nob, nrep")
	print(c(DGP, nob, nrep))
	print(mubeta)
	print(tauvec)

	print("rejfreq10.em, rejfreq5.em, rejfreq1.em")
	print(c(rejfreq10.em, rejfreq5.em, rejfreq1.em))

} # end of DGP loop

stopCluster(cl)

sink()

# Computes the MLE

rm(list=ls())

library(parallel)
library(normalregMix)
library(subgroupLRT)

setwd("~/Dropbox/subgroupLRT/experiments_ks")

outfilename <- "sim_pmle.RData"

ncores <- detectCores()
cl <- makePSOCKcluster(rep('localhost',ncores-1),master='localhost')
clusterEvalQ(cl, library(normalregMix))
clusterEvalQ(cl, library(subgroupLRT))

sink("sim_pmle.out", append=T)

nrep <- 100
# nobset <- c(100,200)
nobset <- c(400)
nnob <- length(nobset)
DGPset <- c(1:12)

m <- 2

alpha.all <- cbind(c(0.2, 0.8),c(0.4, 0.6),c(0.5, 0.5))
mubeta.all <- cbind(c(1,2,3,4),c(1,0,1,2),c(1,2,3,-1),c(-1,2,1,2))

theta.all <- list()
for (i in 1:4){
  for (j in 1:3){
    theta.all[[(i-1)*3+j]] <- list(alpha = alpha.all[,j],
                                   mubeta = mubeta.all[,i])
  }
}


# alpha <- c(0.5,0.5)
# # alpha <- c(0.2,0.8)
# mubeta <- matrix(c(1,2,3,4), nrow=2, ncol=2)
# sigma <- c(1,1)

# rejfreq5all <- matrix(0, nrow=length(DGPset), ncol=3*nnob)
# rejfreq1all <- matrix(0, nrow=length(DGPset), ncol=3*nnob)

for (DGP in DGPset)
{
  alpha <- theta.all[[DGP]]$alpha
  mubeta <- theta.all[[DGP]]$mubeta
  mubeta <- matrix(mubeta, nrow=2, ncol=2)
  sigma <- c(1,1)

	for (inob in 1:nnob)
	{
		time_start <- proc.time()
		nob <- nobset[inob]
		set.seed(123456)
		# clusterEvalQ(cl,set.seed(123456))

		x.all <- matrix(rnorm(nob*nrep),nrow=nob)-1
		y.all <- matrix(double(nob*nrep),nrow=nob)
		for (j in 1:nrep){
		  y.all[,j] <-rnormregmix(n=nob, alpha=alpha, mubeta=mubeta,
		                      sigma=sigma, x=x.all[,j])
		}

		clusterExport(cl,varlist=c("y.all","m", "x.all"))
    clusterSetRNGStream(cl, 123456)
    mleout <- parLapply(cl,1:nrep, function(j) regmixMLE_homo(y=y.all[,j],
                              m=m, x=x.all[,j]))
		coefsum <- t(sapply(mleout,"[[","coefficients"))
		# pvalsum <- t(sapply(emout,"[[","pvals"))
		# print(pvalsum)
		# rejfreq5 <- 100*colMeans(pvalsum < 0.05)
		# rejfreq1 <- 100*colMeans(pvalsum < 0.01)
		# rejfreq5all[DGP,(3*inob-2):(3*inob)] <- rejfreq5
		# rejfreq1all[DGP,(3*inob-2):(3*inob)] <- rejfreq1

		time_end <- proc.time()
		runtime  <- time_end - time_start

		print("DGP, nob, nrep")
		print(c(DGP, nob, nrep))
		# print(runtime)

		print(unlist(theta.all[[DGP]]))
		print(colMeans(coefsum), digits=3)
		# print("rejfreq5, rejfreq1")
		# print(c(rejfreq5,rejfreq1))
		# rm(list=c("y.all", "x.all"))


		# save.image(file = outfilename)
	} # end of inob loop
# system("mail -s 1v2_report kshimotsu@gmail.com < 1v2_output.out");
} # end of DGP loop

stopCluster(cl)

# return(list(outall=outall,runtime=runtime))

sink()

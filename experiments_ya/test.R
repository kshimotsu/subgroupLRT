# R script for testing

source('~/Dropbox/yoichi/ProR/subgroupLRT/experiments_ya/dgptest.R')

y <- as.matrix(y.all[,1])
x <- as.matrix(x.all[,1])
v <- as.matrix(v.all[,1])

source('~/Dropbox/yoichi/ProR/subgroupLRT/experiments_ya/rsubgroup.R')
source('~/Dropbox/yoichi/ProR/subgroupLRT/R/subgroup_funcs_homo.R')

nob <- 100
m <- 2
mubeta <- matrix(c(1,2,1,3), nrow=2, ncol=2) # one-component model in effect
tautrue <- c(0.3, -0.7)
tauvec = as.matrix(tautrue, nrow=2)
sigma <- c(0.5, 0.5)

# x <- matrix(rnorm(nob), nrow = nob) - 1 # corresponds to Z in the paper
# v <- matrix(rnorm(nob), nrow = nob) + 1 # corresponds to X in the paper
y1 <- rsubgroup(n = nob, tau = tautrue, mubeta = mubeta, sigma = sigma, x= x, v = v)

# ninits = 2;
#
# tmp <- subgroupMLEinit_homo(y,x,v, ninits = ninits)
# b0 <- rbind( tmp$tau, tmp$mubeta,  tmp$sigma, tmp$gam ) # (q1*m*ninits + q2 + 1 + p)-vector

#Rcpp::sourceCpp('~/Dropbox/yoichi/ProR/subgroupLRT/src/cppSubgroupMLE_homo.cpp')
#out.short <- cppSubgroupMLE_homo(b0, y, x, v, z = NULL, m = 2, p = 0)
#cppSubgroupMLE_homo(b0, y, x, v, z = matrix(0), m = 2, p = 0, maxit = 1)
#cppSubgroupMLE_homo(b0, y, x, v, z = matrix(0), m = 2, p = 0, maxit = 100, ninits = ninits)

out <- subgroupMLE_homo(y, x, v, m = 2, z = NULL, ninits = 10)



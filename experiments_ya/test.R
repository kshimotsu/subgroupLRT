# R script for testing

source('~/Dropbox/yoichi/ProR/subgroupLRT/experiments_ya/dgptest.R')
source('~/Dropbox/yoichi/ProR/subgroupLRT/R/subgroup_funcs_homo.R')


y <- as.matrix(y.all[,1])
x <- as.matrix(x.all[,1])
v <- as.matrix(v.all[,1])

tmp <- subgroupMLEinit_homo(y,x,v)
b0 <- rbind( tmp$tau, tmp$mubeta,  tmp$sigma, tmp$gam ) # (q1*m*ninits + q2 + 1 + p)-vector


#Rcpp::sourceCpp('~/Dropbox/yoichi/ProR/subgroupLRT/src/cppSubgroupMLE_homo.cpp')
out.short <- cppSubgroupMLE_homo(b0, y, x, v, z = NULL, m = 2, p = 0)

#subgroupMLE_homo <- function (y, x, v, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", "none"),
                              # ninits = 10, epsilon = 1e-08, maxit = 2000,
                              # epsilon.short = 1e-02, maxit.short = 500, binit = NULL)



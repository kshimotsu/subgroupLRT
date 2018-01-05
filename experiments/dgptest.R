# Computes the rejection frequencies of the LRT

rm(list=ls())

library(parallel)
library(normalregMix)
library(subgroupLRT)

nrep <- 10
m <- 2
nob <- 100

mubeta <- matrix(c(1,2,1,3), nrow=2, ncol=2) # one-component model in effect
tautrue <- c(0.3, -0.7)
tauvec = as.matrix(tautrue, nrow=2)
sigma <- c(0.5, 0.5)

x.all <- matrix(rnorm(nob*nrep), nrow = nob) - 1 # corresponds to Z in the paper
v.all <- matrix(rnorm(nob*nrep), nrow = nob) + 1 # corresponds to X in the paper

x <- as.matrix(x.all[,1])
v <- as.matrix(v.all[,1])

y1 <- rsubgroup(n = nob, tau = tautrue, mubeta = mubeta, sigma = sigma, x= x, v = v)

out <- subgroupMLE_homo(y1, x, v, m = 2, z = NULL, ninits = 10)

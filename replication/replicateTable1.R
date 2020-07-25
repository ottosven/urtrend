## ####################################################################
## ####################################################################
## Supplement for
## "Unit Root Testing with Slowly Varying Trends"
## by Sven Otto.
## This R-script allows to reproduce Table 1.
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
library(parallel)
## ##################################
## Cluster setup
## ##################################
if(is.na(strtoi(Sys.getenv(c("SLURM_NTASKS"))))){
  cl = makeCluster(detectCores()-1)
} else {
  ntasks <- strtoi(Sys.getenv(c("SLURM_NTASKS")))
  nslaves <- ntasks-1
  cl = makeCluster(nslaves, type="MPI")
}
## ##################################
## Reproducible random number generator
## ##################################
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
snow::clusterSetupRNG(cl)
## ##################################
## Simulation setting
## ##################################
MC <- 100000
T <- 50000
## ##################################
sim.dist <- function(T, bRANGE){
  innerDEN <- function(j, W, T, B) ( sum((W[(j+1):(j+B)] - W[j])^2)/T )
  tFB <- numeric(length(bRANGE))
  W <- cumsum(rnorm(T,0,sqrt(1/T)))
  for(i in 1:length(bRANGE)){
    b <- bRANGE[i]
    B <- floor(b*T)
    NUM <- sum((W[(B+1):T] - W[1:(T-B)])^2)/T - b*(1-b)
    DEN <- 2*sqrt(b*sum(sapply(1:(T-B), innerDEN, W=W, T=T, B=B))/T)
    tFB[i] <- NUM/DEN
  }
  tFB
}
##
sim.crit <- function(T){
  bRANGE <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) 
  realizations <- parSapply(cl,rep(T,MC),sim.dist, bRANGE=bRANGE)
  levels <- c(0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01, 0.001)
  quantiles <- matrix(nrow = length(levels), ncol = length(bRANGE))
  colnames(quantiles) <- bRANGE
  rownames(quantiles) <- levels 
  for(i in 1:length(bRANGE)){
    quantiles[,i] <- quantile(realizations[i,], levels)
  }
  quantiles
}
##
crits <- sim.crit(T)
crits
write.table(crits,file="./table1.csv", col.names = NA)
Sys.time()-start
stopCluster(cl)
## ####################################################################
## ####################################################################
## Supplement for
## "Unit Root Testing with Slowly Varying Trends"
## by Sven Otto.
## This R-script allows to reproduce Figure 2.
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
library(urtrend)  # install package with remotes::install_github("ottosven/urtrend")
library(parallel)
library(urca)
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
## ##################################
##
sim.statistics <- function(T, sigmaAlpha, hypothesis){
  errorprocess <- function(T) ( rnorm(T) )
  if(hypothesis == 'H0') (cumulation <- function(u) ( cumsum(u) ))
  if(hypothesis == 'H1') (cumulation <- function(u) ( filter(u, 0.9, method = "recursive") ))
  sequence <- c(rnorm(1,0,sigmaAlpha),errorprocess(T-1))
  x <- cumulation(sequence)
  y <- c(x)
  c(
    urtrend::smallb.test(y,floor(T^0.8), HC=TRUE)$teststatistic,
    urca::ur.df(y, type = "drift", lags = 0)@teststat[1],
    urca::ur.ers(y, type = "DF-GLS", model ="constant", lag.max = 0)@teststat
  )
}
##
sim.sizeadjustedpower <- function(){
  SDinitial <- (0:12)/2
  Power <- matrix(ncol = 13, nrow = 3)
  for(i in 1:13){
    realizationsH0 <- parSapply(cl,rep(100,MC),sim.statistics, sigmaAlpha=SDinitial[i], hypothesis = "H0")
    realizationsH1 <- parSapply(cl,rep(100,MC),sim.statistics, sigmaAlpha=SDinitial[i], hypothesis = "H1")
    for(j in 1:3){
      Power[j,i] <- length(which(realizationsH1[j,] < quantile(realizationsH0[j,], 0.05)))/MC
    }
  }
  rownames(Power) <- c('small-B-0.8', 'ADF', 'DFGLS')
  colnames(Power) <- SDinitial
  Power
}
##
RESULTS <- sim.sizeadjustedpower()
##
pdf('./figure2.pdf', width = 15, height = 7.5, pointsize = 14)
plot((0:12)/2, RESULTS[2,], type='l', ylim = c(0,1), lwd = 4, xlab=expression(paste("sd(", x[0], ")")), ylab = "power", cex.axis=1.3, cex.lab=1.3)
lines((0:12)/2, RESULTS[3,], lty = 5, lwd = 4)
lines((0:12)/2, RESULTS[1,], lty= 3, col = 1, lwd = 4)
names.row <- c('ADF', 'DFGLS', expression(paste(tau, "-SB, B=", T^{0.8})))
legend(x=4.8, y=0.97, legend = names.row, col = c(1,1,1), lty= c(1,5,3), lwd=4, cex = 1.3, seg.len=3)
dev.off()
##
Sys.time()-start
stopCluster(cl)
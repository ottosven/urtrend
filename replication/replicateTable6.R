## ####################################################################
## ####################################################################
## Supplement for
## "Unit Root Testing with Slowly Varying Trends"
## by Sven Otto.
## This R-script allows to reproduce Table 6.
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
sim.tablecolumn <- function(T, sigma0, rho, trendtype=1, lambda=3){
  if (trendtype == 1){
    ## SHARP BREAK
    trend <- function(T, c) ( c(rep(c,floor(2*T/3)), rep(0, T-floor(2*T/3))) )
  } else if(trendtype == 2){
    ## U-SHAPED BREAK
    trend <- function(T, c) ( c(rep(c,floor(T/4)), rep(0,T-2*floor(T/4)), rep(c,floor(T/4))) )
  } else if(trendtype == 3){
    ## CONTINUOUS BREAK
    trend <- function(T, c) ( c*c(rep(0,floor(2*T/3)) , (-8/3+4*((floor(2*T/3)+1):T)/T)) )
  } else if(trendtype == 4){
    ## ## U-SHAPED BREAK IN INTERCEPT
    trend <- function(T, c, d = 0.2) ( c*c((1:floor(T/4)/T), -1+((floor(T/4)+1):floor(3*T/4))/T, ((floor(3*T/4)+1):T)/T) )
  } else if(trendtype == 5){
    ## LSTAR BREAK
    trend <- function(T, c, phi = 20) ( c/(1+exp(phi*(1:T - 3*T/4)/T)) )
  } else if(trendtype == 6){
    ## OFFSETTING LSTAR BREAKS
    trend <- function(T, c, phi = 20) ( c/(1+exp(phi*(1:T - T/5)/T)) - 0.5*c/(1+exp(phi*(1:T - 3*T/4)/T)) )
  } else if(trendtype == 7){
    ## TRIANGULAR BREAK
    trend <- function(T, c) ( c*c(2*(1:floor(T/2))/T, (2-2*((floor(T/2)+1):T)/T)) )
  } else {
    ## FOURIER BREAK
    trend <- function(T, c) ( c/2*cos(2*pi*(1:T)/T) )
  }
  ## ##################################
  if(rho == 1){
    cumulation <- function(u) ( cumsum(u) )
  } else {
    cumulation <- function(u) ( filter(u, rho, method = "recursive") )
  }
  errorprocess <- function(T) ( arima.sim(list(order=c(1,0,1), ar=.5, ma=.0), n = T) )
  x0 <- rnorm(1,0,sigma0)
  x <- cumulation(c(x0,errorprocess(T)))
  y <- c(x[-1]) + trend(T, lambda)
  GammasSB <- c(0.5, 0.6, 0.7, 0.8)
  GammasFB <- c(0.2, 0.4, 0.6)
  Bsmallb <- floor(T^GammasSB)
  Bfixedb <- floor(T*GammasFB)
  allTauSB <- numeric(length(GammasSB))
  allTauFB <- numeric(length(GammasFB))
  lags <- urtrend::lagselection.BIC(y,4)
  y <- urtrend::get.prewhitened(y,lags)
  for(i in 1:length(GammasSB)) ( allTauSB[i] <- urtrend::smallb.test(y,Bsmallb[i], HC=TRUE)$teststatistic )
  for(i in 1:length(GammasFB)) ( allTauFB[i] <- urtrend::fixedb.test(y,Bfixedb[i], HC=TRUE)$teststatistic )
  c(allTauSB, allTauFB)
}
##
sim <- function(type){
  realizations <- list()
  realizations[[1]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=1, trendtype=type, lambda=3)
  realizations[[2]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=1, trendtype=type, lambda=6)
  realizations[[3]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=1, trendtype=type, lambda=9)
  realizations[[4]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=0.9, trendtype=type, lambda=3)
  realizations[[5]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=0.9, trendtype=type, lambda=6)
  realizations[[6]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=0.9, trendtype=type, lambda=9)
  realizations[[7]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=1, trendtype=type, lambda=3)
  realizations[[8]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=1, trendtype=type, lambda=6)
  realizations[[9]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=1, trendtype=type, lambda=9)
  realizations[[10]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=0.9, trendtype=type, lambda=3)
  realizations[[11]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=0.9, trendtype=type, lambda=6)
  realizations[[12]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=0.9, trendtype=type, lambda=9)
  CritAll <- c(
    rep(qnorm(0.05),4), 
    get.crit.FB(0.2)$crit.values["0.05"],
    get.crit.FB(0.4)$crit.values["0.05"],
    get.crit.FB(0.6)$crit.values["0.05"]
  )
  statnames <- c("SB05", "SB06", "SB07", "SB08", "FB02", "FB04", "FB06")
  rejectionrates <- matrix(ncol = 12, nrow = 7)
  rownames(rejectionrates) <- statnames
  for(j in 1:12){
    for(i in 1:7){
      rejectionrates[i,j] <- length(which(realizations[[j]][i,] < CritAll[i]))/MC
    }
  }
  rejectionrates
}
##
part1 <- sim(1)
part2 <- sim(2)
part3 <- sim(3)
part4 <- sim(5)
part5 <- sim(8)
table6 <- rbind(part1, part2, part3, part4, part5)
colnames(table6) <- c("100-H0-3","100-H0-6","100-H0-9","100-H1-3","100-H1-6","100-H1-9","300-H0-3","300-H0-6","300-H0-9","300-H1-3","300-H1-6","300-H1-9")
table6
write.table(table6,file="./table6.csv", col.names = NA)
Sys.time()-start
stopCluster(cl)
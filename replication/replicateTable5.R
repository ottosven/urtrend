## ####################################################################
## ####################################################################
## Supplement for
## "Unit Root Testing with Slowly Varying Trends"
## by Sven Otto.
## This R-script allows to reproduce Table 5.
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
  ## ##################################
  ## Conventional unit root tests:
  ## ##################################
  EndersLee<-function(y, P = 0){
    T<-length(y)
    SIN<-sin(2*pi*(1:T)/T)
    COS<-cos(2*pi*(1:T)/T)
    coeff<-lm(diff(y) ~ diff(SIN) + diff(COS))$coefficients
    S <- y - (y[1] - coeff[1] - coeff[2]*SIN[1] - coeff[3]*COS[1]) - coeff[1]*(1:T) - (coeff[2]*SIN + coeff[3]*COS)
    if(P == 0) ( reg <- lm(diff(y) ~ S[1:(T-1)] + diff(SIN) + diff(COS) ) )
    if(P > 0) ( reg <- lm(diff(y)[(P+1):(T-1)] ~ S[(P+1):(T-1)] + diff(SIN)[(P+1):(T-1)] + diff(COS)[(P+1):(T-1)] + embed(diff(S),P+1)[,-1]) )
    reg
    return(coef(summary(reg))[, "t value"][2])
  }
  ##
  DFtypeTests <- function(y, P = c(0,0,0,0)){
    results<- c(
      urca::ur.df(y, type = "drift", lags = P[1])@teststat[1],
      urca::ur.ers(y, type = "DF-GLS", model ="constant", lag.max = P[2])@teststat,
      urca::ur.ers(y, type = "DF-GLS", model ="trend", lag.max = P[3])@teststat,
      EndersLee(y, P[4])
    )
    names(results) <- c('DF-cons', 'DFGLS-cons', 'DFGLS-trend', 'EndersLee')
    results
  }
  ## ##################################
  ## ##################################
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
  ## ##################################
  if(rho == 1){
    cumulation <- function(u) ( cumsum(u) )
  } else {
    cumulation <- function(u) ( filter(u, rho, method = "recursive") )
  }
  errorprocess <- function(T) ( rnorm(T) )
  x0 <- rnorm(1,0,sigma0)
  x <- cumulation(c(x0,errorprocess(T)))
  y <- c(x[-1]) + trend(T, lambda)
  GammasSB <- c(0.5, 0.6, 0.7, 0.8)
  GammasFB <- c(0.2, 0.4, 0.6)
  Bsmallb <- floor(T^GammasSB)
  Bfixedb <- floor(T*GammasFB)
  allTauSB <- numeric(length(GammasSB))
  allTauFB <- numeric(length(GammasFB))
  for(i in 1:length(GammasSB)) ( allTauSB[i] <- urtrend::smallb.test(y,Bsmallb[i], HC=TRUE)$teststatistic )
  for(i in 1:length(GammasFB)) ( allTauFB[i] <- urtrend::fixedb.test(y,Bfixedb[i], HC=TRUE)$teststatistic )
  DF <- DFtypeTests(y)
  c(allTauSB, allTauFB, DF)
}
##
sim.5 <- function(...){
  realizations <- list()
  type <- 5
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
    get.crit.FB(0.6)$crit.values["0.05"],
    c(-2.86, -1.94, -2.89, -4.03)
  )
  statnames <- c("SB05", "SB06", "SB07", "SB08", "FB02", "FB04", "FB06", "ADF", "DFGLS", "DFGLS-t", "EL")
  rejectionrates <- matrix(ncol = 12, nrow = 11)
  rownames(rejectionrates) <- statnames
  for(j in 1:12){
    for(i in 1:11){
      rejectionrates[i,j] <- length(which(realizations[[j]][i,] < CritAll[i]))/MC
    }
  }
  rejectionrates
}
##
sim.6 <- function(...){
  realizations <- list()
  type <- 6
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
    get.crit.FB(0.6)$crit.values["0.05"],
    c(-2.86, -1.94, -2.89, -4.03)
  )
  statnames <- c("SB05", "SB06", "SB07", "SB08", "FB02", "FB04", "FB06", "ADF", "DFGLS", "DFGLS-t", "EL")
  rejectionrates <- matrix(ncol = 12, nrow = 11)
  rownames(rejectionrates) <- statnames
  for(j in 1:12){
    for(i in 1:11){
      rejectionrates[i,j] <- length(which(realizations[[j]][i,] < CritAll[i]))/MC
    }
  }
  rejectionrates
}
##
sim.7 <- function(...){
  realizations <- list()
  type <- 7
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
    get.crit.FB(0.6)$crit.values["0.05"],
    c(-2.86, -1.94, -2.89, -4.03)
  )
  statnames <- c("SB05", "SB06", "SB07", "SB08", "FB02", "FB04", "FB06", "ADF", "DFGLS", "DFGLS-t", "EL")
  rejectionrates <- matrix(ncol = 12, nrow = 11)
  rownames(rejectionrates) <- statnames
  for(j in 1:12){
    for(i in 1:11){
      rejectionrates[i,j] <- length(which(realizations[[j]][i,] < CritAll[i]))/MC
    }
  }
  rejectionrates
}
##
sim.8 <- function(...){
  realizations <- list()
  type <- 8
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
    get.crit.FB(0.6)$crit.values["0.05"],
    c(-2.86, -1.94, -2.89, -4.03)
  )
  statnames <- c("SB05", "SB06", "SB07", "SB08", "FB02", "FB04", "FB06", "ADF", "DFGLS", "DFGLS-t", "EL")
  rejectionrates <- matrix(ncol = 12, nrow = 11)
  rownames(rejectionrates) <- statnames
  for(j in 1:12){
    for(i in 1:11){
      rejectionrates[i,j] <- length(which(realizations[[j]][i,] < CritAll[i]))/MC
    }
  }
  rejectionrates
}
##
part5 <- sim.5()
part6 <- sim.6()
part7 <- sim.7()
part8 <- sim.8()
table5 <- rbind(part5, part6, part7, part8)
colnames(table5) <- c("100-H0-3","100-H0-6","100-H0-9","100-H1-3","100-H1-6","100-H1-9","300-H0-3","300-H0-6","300-H0-9","300-H1-3","300-H1-6","300-H1-9")
table5
write.table(table5,file="./table5.csv")
Sys.time()-start
stopCluster(cl)
## ####################################################################
## ####################################################################
## Supplement for
## "Unit Root Testing with Slowly Varying Trends"
## by Sven Otto.
## This R-script allows to reproduce Table 3.
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
sim.tablecolumn <- function(T, sigma0, rho, p=0){
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
  ##
  LagSelectDFGLS <- function(y, Pmax = 0, model=c("constant", "trend")){
    model <- match.arg(model)
    lagsBIC <- 0
    if(Pmax > 0){
      nobs<-length(y)
      ahat <- 1 - 7.0/nobs
      ya <- c(y[1], y[2:nobs]-ahat*y[1:(nobs-1)])
      za1 <- c(1, rep(1-ahat, nobs-1))
      if(model == "comstant"){
        yd.reg <- summary(lm(ya ~ -1 + za1))
        yd <- y - coef(yd.reg)[1]
      }else{
        trd <- 1:nobs
        za2 <- c(1, trd[2:nobs]-ahat*trd[1:(nobs-1)])
        yd.reg <- summary(lm(ya ~ -1 + za1 + za2))
        yd <- y - coef(yd.reg)[1] - coef(yd.reg)[2]*trd
      }
      yd.l <-  yd[1:(nobs-1)]
      yd.diff <- diff(yd)
      allBIC<-numeric(Pmax+1)
      allBIC[1] <- BIC(lm(yd.diff ~ -1 + yd.l))
      for( i in 1:Pmax){
        yd.dlags <- embed(diff(yd), i+1)[, -1]
        data.dfgls <- data.frame(cbind(yd.diff[-(1:i)], yd.l[-(1:i)], yd.dlags))
        colnames(data.dfgls) <- c("yd.diff", "yd.lag", paste("yd.diff.lag", 1:i, sep=""))
        dfgls.form <- formula(paste("yd.diff ~ -1 + ", paste(colnames(data.dfgls)[-1], collapse=" + ")))
        dfgls.reg <- lm(dfgls.form, data=data.dfgls)
        allBIC[i+1] <- BIC(dfgls.reg)
      }
      lagsBIC <- which.min(allBIC)-1
    }
    lagsBIC
  }
  ##
  LagSelectDF <- function(y, Pmax = 0){
    lagsBIC <- 0
    if(Pmax > 0){
      z <- diff(y)
      xx <- embed(z, Pmax+1)
      z.diff <- xx[, 1]
      z.lag.1 <- y[(Pmax+1):length(z)]
      allBIC<-numeric(Pmax+1)
      allBIC[1] <- BIC(lm(z.diff ~ z.lag.1 + 1))
      for(i in 1:Pmax){
        z.diff.lag = xx[, 2:(i+1)]
        result <- lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
        allBIC[i+1] <- BIC(result)
      }
      lagsBIC <- which.min(allBIC)-1
    }
    lagsBIC
  }
  ##
  LagSelectEL <- function(y, Pmax = 0){
    T <- length(y)
    lagsBIC <- 0
    if(Pmax > 0){
      nobs <- length(y)
      SIN<-sin(2*pi*(1:nobs)/nobs)
      COS<-cos(2*pi*(1:nobs)/nobs)
      coeff<-lm(diff(y) ~ diff(SIN) + diff(COS))$coefficients
      S <- y - (y[1] - coeff[1] - coeff[2]*SIN[1] - coeff[3]*COS[1]) - coeff[1]*(1:T) - (coeff[2]*SIN + coeff[3]*COS)
      allBIC<-numeric(Pmax+1)
      allBIC[1] <- BIC(lm(diff(y) ~ S[1:(nobs-1)] + diff(SIN) + diff(COS) ))
      for(i in 1:Pmax){
        result <- lm(diff(y)[(i+1):(nobs-1)] ~ S[(i+1):(nobs-1)] + diff(SIN)[(i+1):(nobs-1)] + diff(COS)[(i+1):(nobs-1)] + embed(diff(S),i+1)[,-1]) 
        allBIC[i+1] <- BIC(result)
      }
      lagsBIC <- which.min(allBIC)-1
    }
    lagsBIC
  }
  ##
  LagSelectionDFType <- function(y, Pmax = 0){
    lags <- c(LagSelectDF(y,Pmax), LagSelectDFGLS(y,Pmax,model="constant"), LagSelectDFGLS(y,Pmax,model="trend"), LagSelectEL(y,Pmax) )
  }
  ## ##################################
  ## ##################################
  if(rho == 1){
    cumulation <- function(u) ( cumsum(u) )
  } else {
    cumulation <- function(u) ( filter(u, rho, method = "recursive") )
  }
  if(p == 0){
    errorprocess <- function(T) ( rnorm(T) )
  } else {
    errorprocess <- function(T) ( arima.sim(list(order=c(1,0,1), ar=.5, ma=.0), n = T) )
  }
  x0 <- rnorm(1,0,sigma0)
  x <- cumulation(c(x0,errorprocess(T)))
  y <- c(x[-1])
  GammasSB <- c(0.5, 0.6, 0.7, 0.8)
  GammasFB <- c(0.2, 0.4, 0.6)
  Bsmallb <- floor(T^GammasSB)
  Bfixedb <- floor(T*GammasFB)
  allTauSB <- numeric(length(GammasSB))
  allTauFB <- numeric(length(GammasFB))
  DFlags <- rep(0,4)
  if(p == 1){
    DFlags <- rep(1,4)
    y <- urtrend::get.prewhitened(y,1)
  }
  if (p == "BIC"){
    DFlags <- LagSelectionDFType(y, 4)
    lags <- urtrend::lagselection.BIC(y,4)
    y <- urtrend::get.prewhitened(y,lags)
  }
  for(i in 1:length(GammasSB)) ( allTauSB[i] <- urtrend::smallb.test(y,Bsmallb[i], HC=TRUE)$teststatistic )
  for(i in 1:length(GammasFB)) ( allTauFB[i] <- urtrend::fixedb.test(y,Bfixedb[i], HC=TRUE)$teststatistic )
  DF <- DFtypeTests(y, DFlags)
  c(allTauSB, allTauFB, DF)
}
##
sim.1 <- function(...){
  realizations <- list()
  realizations[[1]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=1, p=0)
  realizations[[2]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=0.9, p=0)
  realizations[[3]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=1, p=0)
  realizations[[4]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=0.9, p=0)
  realizations[[5]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(5), rho=1, p=0)
  realizations[[6]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(5), rho=0.9, p=0)
  realizations[[7]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(5), rho=1, p=0)
  realizations[[8]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(5), rho=0.9, p=0)
  realizations[[9]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(10), rho=1, p=0)
  realizations[[10]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(10), rho=0.9, p=0)
  realizations[[11]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(10), rho=1, p=0)
  realizations[[12]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(10), rho=0.9, p=0)  
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
sim.2 <- function(...){
  realizations <- list()
  realizations[[1]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=1, p=1)
  realizations[[2]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=0.9, p=1)
  realizations[[3]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=1, p=1)
  realizations[[4]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=0.9, p=1)
  realizations[[5]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(5), rho=1, p=1)
  realizations[[6]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(5), rho=0.9, p=1)
  realizations[[7]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(5), rho=1, p=1)
  realizations[[8]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(5), rho=0.9, p=1)
  realizations[[9]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(10), rho=1, p=1)
  realizations[[10]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(10), rho=0.9, p=1)
  realizations[[11]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(10), rho=1, p=1)
  realizations[[12]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(10), rho=0.9, p=1)  
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
sim.3 <- function(...){
  realizations <- list()
  realizations[[1]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=1, p="BIC")
  realizations[[2]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=0, rho=0.9, p="BIC")
  realizations[[3]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=1, p="BIC")
  realizations[[4]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=0, rho=0.9, p="BIC")
  realizations[[5]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(5), rho=1, p="BIC")
  realizations[[6]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(5), rho=0.9, p="BIC")
  realizations[[7]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(5), rho=1, p="BIC")
  realizations[[8]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(5), rho=0.9, p="BIC")
  realizations[[9]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(10), rho=1, p="BIC")
  realizations[[10]] <- parSapply(cl,rep(100,MC),sim.tablecolumn, sigma0=sqrt(10), rho=0.9, p="BIC")
  realizations[[11]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(10), rho=1, p="BIC")
  realizations[[12]] <- parSapply(cl,rep(300,MC),sim.tablecolumn, sigma0=sqrt(10), rho=0.9, p="BIC")  
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
part1 <- sim.1()
part2 <- sim.2()
part3 <- sim.3()
table3 <- rbind(part1, part2, part3)
colnames(table3) <- c("100-H0-0","100-H1-0","300-H0-0","300-H1-0","100-H0-5","100-H1-5","300-H0-5","300-H1-5","100-H0-10","100-H1-10","300-H0-10","300-H1-10")
table3
write.table(table3,file="./table3.csv", col.names = NA)
Sys.time()-start
stopCluster(cl)
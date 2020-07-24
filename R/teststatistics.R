#' @import stats
NULL

#' Small-b unit root test
#'
#' Performs the small-b unit root test
#'
#' @param y Vector to be tested for a unit root.
#' @param B An integer defining the blocklength; default is \eqn{T^{0.7}}. 
#' @param HC A logical variable defining whether the heteroskedasticity robust version (\eqn{HC = TRUE}) 
#' or the homoskedasticity version (\eqn{HC = FALSE}) should be used; default is \eqn{HC = TRUE}.
#' 
#' @return A list containing the following components:
#' \item{teststatistic}{The value of the small-b test statistic.}
#' \item{p.value}{The p-value of the small-b test against stationarity.}
#' @export
#'
#' @examples
#' series1 <- rnorm(100) + c(rep(0,60), rep(5, 40))
#' smallb.test(series1)
#' series2 <- cumsum(rnorm(100)) + c(rep(0,60), rep(5, 40))
#' smallb.test(series2)
smallb.test <- function(y, B = floor(length(y)^{0.7}), HC = TRUE){
  T <- length(y)
  # get numerator and denominator of phihat
  get.phihat <- function(y,B,T){
    blocksums<-function(j){
      DeltaY<-diff(y)[j:(j+B-1)]
      DeltaYj<-y[j:(j+B-1)]-y[j]
      return(c(DeltaY %*% DeltaYj, DeltaYj %*% DeltaYj))
    }
    NumDen<-rowSums(sapply(1:(T-B),blocksums))
    return(c(NumDen[1]/(B^(3/2)*T^(1/2)), NumDen[2]/(B^2*T), NumDen[1]/NumDen[2]))
  }
  ND <- get.phihat(y,B,T)
  N <- ND[1]
  D <- ND[2]
  Rhohat <- ND[3] + 1
  vT2 <- ((T-B)*(2*B-1) - 2*(B-2))/(3*B*(T-B))
  yhat<-c(0,y[2:T]-Rhohat*y[1:(T-1)])
  sigma2block<-function(j) ( sd(yhat[(j+1):(j+B)])^2 )
  sigma2SB <- sum(sapply(1:(T-B),sigma2block))/(T-B)
  if(HC == TRUE){
    sigma4block<-function(j)( sd(yhat[(j+1):(j+B)])^2*(yhat[j]-mean(yhat))^2 )
    kappa2 <- sum(sapply(1:(T-B),sigma4block))/(T-B)
    tauSB.stat <- N/sqrt(vT2*D*kappa2/sigma2SB)
  } else {
    tauSB.stat <- N/sqrt(vT2*sigma2SB*D)
  }
  pval <- pnorm(tauSB.stat)
  return(list(teststatistic = tauSB.stat, p.value = pval))
}


#' Fixed-b unit root test
#'
#' Performs the fixed-b unit root test
#'
#' @param y Vector to be tested for a unit root.
#' @param B An integer defining the blocklength; default is \eqn{0.2*T}.
#' @param HC A logical variable defining whether the heteroskedasticity robust version (\eqn{HC = TRUE})
#' or the homoskedasticity version (\eqn{HC = FALSE}) should be used; default is \eqn{HC = TRUE}..
#'
#' @return A list containing the following components:
#' \item{teststatistic}{The value of the small-b test statistic.}
#' \item{crit.values}{A vector of critical valiues for different significance levels.}
#' \item{crit.values}{A logical vector indicating whether the unit root hypothesis is rejected in favor of stationarity for at different significance levels.}
#' \item{crit.blocklength}{The relative blocklength for which the critical values are reported.}
#' @export
#'
#' @examples
#' series1 <- rnorm(100) + c(rep(0,60), rep(5, 40))
#' fixedb.test(series1)
#' series2 <- cumsum(rnorm(100)) + c(rep(0,60), rep(5, 40))
#' fixedb.test(series2)
fixedb.test <- function(y, B = floor(0.2*length(y)), HC = TRUE){
  phihat<-function(y,B){
    T<-length(y)
    B<-min(B,T-1)
    blocksums<-function(j){
      dy<-diff(y)[j:(j+B-1)]
      dyj<-y[j:(j+B-1)]-y[j]
      return(c(dy %*% dyj ,dyj %*% dyj))
    }
    yy<-rowSums(sapply(1:(T-B),blocksums))
    return( c(yy[1]/yy[2], yy[1]/sqrt(yy[2])) )
  }
  T <- length(y)
  Phihat <- phihat(y,B)
  Rhohat <- Phihat[1] + 1
  uhat <- c(0,y[2:T]-(Phihat[1]+1)*y[1:(T-1)])
  if(HC == TRUE){
    eta<-numeric(T)
    for(s in 1:T) ( eta[s] <- sum( ( uhat[1:s] - (sum(uhat[1:s])/s) )^2 ) )
    eta <- eta/sum((uhat - mean(uhat))^2)
    eta[1]<-0
    eta[T]<-1
    Ttilde <- T*5
    Btilde <- B*5
    Tildegrid <- (1:Ttilde)/Ttilde
    i<-1
    etainv <- numeric(Ttilde)
    for(t in 1:Ttilde){
      while( eta[i] < Tildegrid[t] ) ( i <- i+1 )
      etainv[t] <-  ((i-2)+(Tildegrid[t]-eta[i-1])/(eta[i]-eta[i-1]))/T
    }
    PhihatHC<-phihat(y[pmax(floor(etainv*T),1)],Btilde)
    tauFB.stat <- PhihatHC[2]/(sd(uhat)*sqrt(B))
  } else {
    tauFB.stat <- Phihat[2]/(sd(uhat)*sqrt(B))
  }
  crit.values <- get.crit.FB(B/T)$crit.values
  return(list(
    teststatistic = tauFB.stat, 
    crit.values = crit.values, 
    rejection = tauFB.stat < crit.values,
    crit.bocklength = get.crit.FB(B/T)$crit.blocklength)
    )
}




#' Critical values for the fixed-b unit root test
#'
#' Provides critical values for the fixed-b unit root test for a given relative blocklength for different significance levels.
#'
#' @param b Relative blocklength.
#' 
#' @return A list containing the following components:
#' \item{crit.values}{A vector of critical values for different significance levels.}
#' \item{rel.blocklength}{The relative blocklength for which the critical values are reported. If the critical values for the input value of b are not available, the nearest available value is considered.}
#' @export
#'
#' @examples
#' get.crit.FB(0.2)
#' get.crit.FB(0.48)
get.crit.FB <- function(b){
  blocklengths <- 1:9/10
  b.index <- which.min(abs(b-blocklengths))
  crit.val <- FB.crit[b.index,]
  names(crit.val) <- colnames(FB.crit)
  return(list(crit.values = crit.val, crit.blocklength = blocklengths[b.index]))
}




#' Prewhitened series
#'
#' Returns the prewhitened series given some lag length.
#'
#' @param y Vector of time series values to be pre-whitened.
#' @param p Lag length of the AR model.
#' 
#' @return A vector containing the values of the pre-whitened series.
#' @export
#'
#' @examples
#' get.prewhitened(rnorm(100))
#' get.prewhitened(arima.sim(list(order=c(1,0,1), ar=.5, ma=.0), n = 100))
get.prewhitened<-function(y, p = 1){
  if(p == 0) ( y.pw <- y )
  else{
    T<-length(y)
    z <- diff(y)
    X <- matrix(data=NA,nrow=T-1,ncol=p+1)
    Y <- matrix(data=NA,nrow=T,ncol=p)
    X[,1]<-y[1:(T-1)]
    for(k in 1:p){
      X[(k+1):(T-1),k+1] <- diff(y)[1:(T-1-k)]
      Y[(k+1):T,k] <- y[1:(T-k)]
    }
    Y<-as.matrix(Y[,(1:p)])
    alphas <- lm(z ~ X - 1)$coefficients[-1]
    y.pw <- c(na.omit(y - Y %*% alphas))
  }
  return(y.pw)
}



#' BIC laglegth selection
#'
#' Returns the laglength in the prewhitening regression for which the BIC criterion is minimized.
#'
#' @param y Vector of time series values to be pre-whitened.
#' @param Pmax Maximum lag length to be considered.
#' 
#' @return Lag length determined by the BIC criterion.
#' @export
#'
#' @examples
#' lagselection.BIC(rnorm(100))
lagselection.BIC <- function(y, Pmax = 5){
  lagsBIC <- 0
  if(Pmax > 0){
    nobs <- length(y)
    z <- diff(y)
    X <- matrix(data=NA,nrow=nobs-1,ncol=Pmax+1)
    X[,1]<-y[1:(nobs-1)]
    for(k in 1:Pmax){
      X[(k+1):(nobs-1),k+1] <- diff(y)[1:(nobs-1-k)]
    }
    allBIC<-numeric(Pmax+1)
    allBIC[1] <- BIC(lm(z ~ X[,1] - 1))
    for(i in 1:Pmax){
      result <- lm(z ~ X[,1:(i+1)] - 1)
      allBIC[i+1] <- BIC(result)
    }
    lagsBIC <- which.min(allBIC)-1
  }
  return(lagsBIC)
}
## ####################################################################
## ####################################################################
## Supplement for
## "Unit Root Testing with Slowly Varying Trends"
## by Sven Otto.
## This R-script allows to reproduce Figure 1.
## ####################################################################
## ####################################################################
rm(list=ls())
start<-Sys.time()
## ##################################
##
##### SHARP BREAK #################
ELType1 <- function(T, c) ( c(rep(c,floor(2*T/3)), rep(0, T-floor(2*T/3))) )  
##### U-SHAPED BREAK ##############
ELType2 <- function(T, c) ( c(rep(c,floor(T/4)), rep(0,T-2*floor(T/4)), rep(c,floor(T/4))) )
##### CONTINUOUS BREAK ############
ELType3 <- function(T, c) ( c*c(rep(0,floor(2*T/3)) , (-8/3+4*((floor(2*T/3)+1):T)/T)) )
### U-SHAPED BREAK IN INTERCEPT ###
ELType4 <- function(T, c, d = 0.2) ( c*c((1:floor(T/4)/T), -1+((floor(T/4)+1):floor(3*T/4))/T, ((floor(3*T/4)+1):T)/T) )
##### LSTAR BREAK #################
ELType5 <- function(T, c, phi = 20) ( c/(1+exp(phi*(1:T - 3*T/4)/T)) )
##### OFFSETTING LSTAR BREAKS #####
ELType6 <- function(T, c, phi = 20) ( c/(1+exp(phi*(1:T - T/5)/T)) - 0.5*c/(1+exp(phi*(1:T - 3*T/4)/T)) )
##### TRIANGULAR BREAK ############
ELType7 <- function(T, c) ( c*c(2*(1:floor(T/2))/T, (2-2*((floor(T/2)+1):T)/T)) )
##### FOURIER BREAK ###############
ELType8 <- function(T, c) ( c/2*cos(2*pi*(1:T)/T) )
## ##################################
##
pdf(file= './figure1.pdf' ,onefile=T,width = 9, height = 5)
par(mfrow=c(2,4))
plot((1:100)/100, ELType1(100,3), type='l', main='sharp break', ylab ='d(r)', ylim = c(-0.5,3.5), xlab='r')
plot((1:100)/100, ELType2(100,3), type='l', main='u-shaped break', ylab ='d(r)', ylim = c(-0.5,3.5), xlab='r')
plot((1:100)/100, ELType3(100,3), type='l', main='continuous break', ylab ='d(r)', ylim = c(-0.5,4.5), xlab='r')
plot((1:100)/100, ELType4(100,3), type='l', main='u-shaped break in intercept', ylab ='d(r)', ylim = c(-3,3.5), xlab='r')
plot((1:100)/100, ELType5(100,3), type='l', main='LSTAR break', ylab ='d(r)', ylim = c(-0.5,3.5), xlab='r')
plot((1:100)/100, ELType6(100,3), type='l', main='offsetting LSTAR breaks', ylab ='d(r)', ylim = c(-2,2), xlab='r')
plot((1:100)/100, ELType7(100,3), type='l', main='triangular break', ylab ='d(r)', ylim = c(-0.5,3.5), xlab='r')
plot((1:100)/100, ELType8(100,3), type='l', main='Fourier break', ylab ='d(r)', ylim = c(-2,2), xlab='r')
dev.off()
Sys.time()-start
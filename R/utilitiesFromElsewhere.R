## Mon Jun 03 00:13:35 2013
## Original file Copyright 2013 A.C. Guidoum
## This file is part of the R package kedd.
## Arsalane Chouaib GUIDOUM <acguidoum@usthb.dz> and <starsalane@gmail.com>
## Department of Probabilities-Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algeris
## Algeria
##############################################################################
## Kernels functions


kernel.fun <- function(x,...) UseMethod("kernel.fun")

kernel.fun.default <- function(x=NULL,deriv.order=0,kernel=c("gaussian","epanechnikov",
                                  "uniform","triangular","triweight","tricube",
                                  "biweight","cosine","silverman"),...)
              {
     if (any(deriv.order < 0 || deriv.order != round(deriv.order)))
             stop("argument 'deriv.order' is non-negative integers")
     r <- deriv.order
     if (missing(kernel)) kernel <- "gaussian"
     if (is.null(x)){
         if (kernel == "gaussian"){x <- seq(-5,5,length=1000)}
            else if (kernel == "silverman"){x <- seq(-10,10,length=1000)}
                  else {x <- seq(-1.5,1.5,length=1000)}
        }
     kx <- kernel_fun_der(kernel,u=x,deriv.order=r)
     structure(list(kernel = kernel,deriv.order=r,x=x,kx=kx),class="kernel.fun")
}

##############
##############

kernel.conv <- function(x,...) UseMethod("kernel.conv")

kernel.conv.default <- function(x=NULL,deriv.order=0,kernel=c("gaussian","epanechnikov",
                                "uniform","triangular","triweight","tricube",
                                "biweight","cosine","silverman"),...)
              {
     if (any(deriv.order < 0 || deriv.order != round(deriv.order)))
        stop("argument 'deriv.order' is non-negative integers")
     r <- deriv.order
     if (missing(kernel)) kernel <- "gaussian"
     if (is.null(x)){
         if (kernel == "gaussian"){x <- seq(-8,8,length=1000)}
            else if (kernel == "silverman"){x <- seq(-10,10,length=1000)}
                  else {x <- seq(-2.5,2.5,length=1000)}
        }
     kx <- kernel_fun_conv(kernel,u=x,deriv.order=r)
     structure(list(kernel = kernel,deriv.order=r,x=x,kx=kx),class="kernel.conv")
}

#############
#############

plot.kernel.fun1d <- function(f,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                      type="l",las=1,lwd=1,...)
                  {
    class(f) <- "kernel.fun"
    r <- f$deriv.order
    kernel <- f$kernel
    if(is.null(xlab)) xlab <- "x"
    if(is.null(ylab)) ylab <- ""
    if(is.null(main)){
    	     if(r != 0) {main <- paste("Derivative of ",kernel,"kernel")}else{
    	                 main <- paste(kernel,"kernel")}
    	                }
    if(is.null(sub)){
    	     if(r != 0) {sub <- paste("Derivative order = ",r)}
    	                }
    plot.default(f$x,f$kx,type=type,las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		       main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)
}

plot.kernel.fun <- function(x,...) plot.kernel.fun1d (x,...)


################################
################################

plot.kernel.conv1d <- function(f,main=NULL,sub = NULL, xlab=NULL, ylab=NULL,
                             type="l",las=1,lwd=1,...)
                  {
    class(f) <- "kernel.conv"
    r <- f$deriv.order
    kernel <- f$kernel
    if(is.null(xlab)) xlab <- "x"
    if(is.null(ylab)) ylab <- ""
    if(is.null(main)){
    	     if(r != 0) {main <- paste("Convolution of derivative",kernel,"kernel")}else{
    	                 main <- paste("Convolution of",kernel,"kernel")}
    	                }
    if(is.null(sub)){
	     if(r !=0) {sub <- paste("Derivative order = ",r)}
	                }
    plot.default(f$x,f$kx,type=type,las=las,lwd=lwd,xlab=xlab,ylab=ylab,
		       main=main,sub=sub,font.main=2,cex.main=0.9,font.sub=2,cex.sub=0.7,...)
}

plot.kernel.conv <- function(x,...) plot.kernel.conv1d (x,...)

## Sat Jul 20 02:11:02 2013
## Original file Copyright 2013 A.C. Guidoum
## This file is part of the R package kedd.
## Arsalane Chouaib GUIDOUM <acguidoum@usthb.dz> and <starsalane@gmail.com>
## Department of Probabilities-Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algeris
## Algeria

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## A copy of the GNU General Public License is available at
## http://www.r-project.org/Licenses/
## Unlimited use and distribution (see LICENCE).
## kedd : Kernel estimator and bandwidth selection for density and its derivatives.
###################################################################################################

#### r(th) derivative of Kernel functions K^r(x)

kernel_fun_der <- function(kernel,u,deriv.order=0)
                 {
     if (any(deriv.order < 0 || deriv.order != round(deriv.order)))
             stop("argument 'deriv.order' is non-negative integers")
     r <- deriv.order
     if (kernel=="gaussian")          {Kr <- expression( dnorm(X) )}
     else if (kernel=="epanechnikov") {Kr <- expression( (3/4)*(1-X^2)) }
     else if (kernel=="uniform")      {Kr <- expression( 0.5 ) }
     else if (kernel=="triweight")    {Kr <- expression( (35/32)*(1-X^2)^3 )}
     else if (kernel=="biweight")     {Kr <- expression( (15/16)*(1-X^2)^2 )}
     else if (kernel=="cosine")       {Kr <- expression( (pi/4)*cos((pi*X)/2) )}
     else if (kernel=="triangular")   {Kr <- expression(1-abs(X));Kr1 <- expression(1-X); Kr2 <- expression(1+X)}
     else if (kernel=="tricube")      {Kr <- expression((70/81)*(1-abs(X)^3)^3);Kr1 <- expression((70/81)*(1-X^3)^3); Kr2 <- expression((70/81)*(1+X^3)^3)}
     else if (kernel=="silverman")    {r <- r%%8;Kr <- expression(0.5*exp(-abs(X)/sqrt(2))*sin((abs(X)/sqrt(2))+0.25*pi));Kr1 <- expression(0.5*exp(-X/sqrt(2)) * sin((X/sqrt(2)) + 0.25*pi)); Kr2 <- expression(0.5*exp(X/sqrt(2)) * sin((-X/sqrt(2)) + 0.25*pi))}

     if (kernel=="epanechnikov" && r >= 3) stop(" 'epanechnikov kernel derivative = 0' for 'order >= 3' ")
     if (kernel=="uniform" && r >= 1)      stop(" 'uniform kernel derivative = 0' for 'order >= 1' ")
     if (kernel=="triweight" && r >= 7)    stop(" 'triweight kernel derivative = 0' for 'order >= 7' ")
     if (kernel=="biweight" && r >= 5)     stop(" 'biweight kernel derivative = 0' for 'order >= 5' ")
     if (kernel=="triangular" && r >= 2)   stop(" 'triangular kernel derivative = 0' for 'order >= 2' ")
     if (kernel=="tricube" && r >= 10)     stop(" 'tricube kernel derivative = 0' for 'order >= 10' ")

        if (r == 0) {DKr <- Kr
        if (kernel=="gaussian" || kernel =="silverman"){
             K <- function(X) eval(DKr);fx <- K(u)}else{
             K <- function(X) eval(DKr)* (X >= -1 & X <= 1);fx <- K(u)}
        } else { if (kernel=="gaussian"){
		     if (r == 1){fx <- -(1/2)*u*exp(-(1/2)*u^2)*sqrt(2)/sqrt(pi)}
		     else if (r == 2){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^2-1)/sqrt(pi) }
		     else if (r == 3){fx <- -(1/2)*u*exp(-(1/2)*u^2)*sqrt(2)*(u^2-3)/sqrt(pi)}
		     else if (r == 4){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^4-6*u^2+3)/sqrt(pi)}
		     else if (r == 5){fx <- -(1/2)*u*exp(-(1/2)*u^2)*sqrt(2)*(u^4-10*u^2+15)/sqrt(pi)}
		     else if (r == 6){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^6-15*u^4+45*u^2-15)/sqrt(pi)}
		     else if (r == 7){fx <- -(1/2)*u*exp(-(1/2)*u^2)*sqrt(2)*(u^6-21*u^4+105*u^2-105)/sqrt(pi)}
		     else if (r == 8){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^8-28*u^6+210*u^4-420*u^2+105)/sqrt(pi)}
		     else if (r == 9){fx <- -(1/2)*u*exp(-(1/2)*u^2)*sqrt(2)*(u^8-36*u^6+378*u^4-1260*u^2+945)/sqrt(pi)}
		     else if (r == 10){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^10-45*u^8+630*u^6-3150*u^4+4725*u^2-945)/sqrt(pi)}
		     else if (r == 11){fx <- -(1/2)*u*exp(-(1/2)*u^2)*sqrt(2)*(u^10-55*u^8+990*u^6-6930*u^4+17325*u^2-10395)/sqrt(pi)}
		     else if (r == 12){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^12-66*u^10+1485*u^8-13860*u^6+51975*u^4-62370*u^2+10395)/sqrt(pi)}
		     else if (r == 13){fx <- -(1/2)*u*exp(-(1/2)*u^2)*sqrt(2)*(u^12-78*u^10+2145*u^8-25740*u^6+135135*u^4-270270*u^2+135135)/sqrt(pi)}
		     else if (r == 14){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^14-91*u^12+3003*u^10-45045*u^8+315315*u^6-945945*u^4+945945*u^2-135135)/sqrt(pi)}
		     else if (r == 15){fx <- -(1/2)*u*exp(-(1/2)*u^2)*sqrt(2)*(u^14-105*u^12+4095*u^10-75075*u^8+675675*u^6-2837835*u^4+4729725*u^2-2027025)/sqrt(pi)}
		     else if (r == 16){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^16-120*u^14+5460*u^12-120120*u^10+1351350*u^8-7567560*u^6+18918900*u^4-16216200*u^2+2027025)/sqrt(pi)}
		     else if (r == 17){fx <- -(1/2)*u*exp(-(1/2)*u^2)*sqrt(2)*(u^16-136*u^14+7140*u^12-185640*u^10+2552550*u^8-18378360*u^6+64324260*u^4-91891800*u^2+34459425)/sqrt(pi)}
		     else if (r == 18){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^18-153*u^16+9180*u^14-278460*u^12+4594590*u^10-41351310*u^8+192972780*u^6-413513100*u^4+310134825*u^2-34459425)/sqrt(pi)}
		     else if (r == 19){fx <- -(1/2)*u*exp(-(1/2)*u^2)*sqrt(2)*(u^18-171*u^16+11628*u^14-406980*u^12+7936110*u^10-87297210*u^8+523783260*u^6-1571349780*u^4+1964187225*u^2-654729075)/sqrt(pi)}
		     else if (r == 20){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^20-190*u^18+14535*u^16-581400*u^14+13226850*u^12-174594420*u^10+1309458150*u^8-5237832600*u^6+9820936125*u^4-6547290750*u^2+654729075)/sqrt(pi)}
		     else if (r == 21){fx <- -(1/2)*exp(-(1/2)*u^2)*sqrt(2)*u*(u^20-210*u^18+17955*u^16-813960*u^14+21366450*u^12-333316620*u^10+3055402350*u^8-15713497800*u^6+41247931725*u^4-45831035250*u^2+13749310575)/sqrt(pi)}
		     else if (r == 22){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^22-231*u^20+21945*u^18-1119195*u^16+33575850*u^14-611080470*u^12+6721885170*u^10-43212118950*u^8+151242416325*u^6-252070693875*u^4+151242416325*u^2-13749310575)/sqrt(pi)}
		     else if (r == 23){fx <- -(1/2)*exp(-(1/2)*u^2)*sqrt(2)*u*(u^22-253*u^20+26565*u^18-1514205*u^16+51482970*u^14-1081142370*u^12+14054850810*u^10-110430970650*u^8+496939367925*u^6-1159525191825*u^4+1159525191825*u^2-316234143225)/sqrt(pi)}
		     else if (r == 24){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^24-276*u^22+31878*u^20-2018940*u^18+77224455*u^16-1853386920*u^14+28109701620*u^12-265034329560*u^10+1490818103775*u^8-4638100767300*u^6+6957151150950*u^4-3794809718700*u^2+316234143225)/sqrt(pi)}
		     else if (r == 25){fx <- -(1/2)*exp(-(1/2)*u^2)*sqrt(2)*u*(u^24-300*u^22+37950*u^20-2656500*u^18+113565375*u^16-3088978200*u^14+54057118500*u^12-602350749000*u^10+4141161399375*u^8-16564645597500*u^6+34785755754750*u^4-31623414322500*u^2+7905853580625)/sqrt(pi)}
		     else if (r == 26){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^26-325*u^24+44850*u^22-3453450*u^20+164038875*u^18-5019589575*u^16+100391791500*u^14-1305093289500*u^12+10767019638375*u^10-53835098191875*u^8+150738274937250*u^6-205552193096250*u^4+102776096548125*u^2-7905853580625)/sqrt(pi)}
		     else if (r == 27){fx <- -(1/2)*exp(-(1/2)*u^2)*sqrt(2)*u*(u^26-351*u^24+52650*u^22-4440150*u^20+233107875*u^18-7972289325*u^16+180705224700*u^14-2710578370500*u^12+26428139112375*u^10-161505294575625*u^8+581419060472250*u^6-1109981842719750*u^4+924984868933125*u^2-213458046676875)/sqrt(pi)}
		     else if (r == 28){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^28-378*u^26+61425*u^24-5651100*u^22+326351025*u^20-12401338950*u^18+316234143225*u^16-5421156741000*u^14+61665657928875*u^12-452214824811750*u^10+2034966711652875*u^8-5179915266025500*u^6+6474894082531875*u^4-2988412653476250*u^2+213458046676875)/sqrt(pi)}
		     else if (r == 29){fx <- -(1/2)*exp(-(1/2)*u^2)*sqrt(2)*u*(u^28-406*u^26+71253*u^24-7125300*u^22+450675225*u^20-18928359450*u^18+539458244325*u^16-10480903032600*u^14+137561852302875*u^12-1192202719958250*u^10+6557114959770375*u^8-21459648959248500*u^6+37554385678684875*u^4-28887988983603750*u^2+6190283353629375)/sqrt(pi)}
		     else if (r == 30){fx <- (1/2)*exp(-(1/2)*u^2)*sqrt(2)*(u^30-435*u^28+82215*u^26-8906625*u^24+614557125*u^22-28392539175*u^20+899097073875*u^18-19651693186125*u^16+294775397791875*u^14-2980506799895625*u^12+19671344879311125*u^10-80473683597181875*u^8+187771928393424375*u^6-216659917377028125*u^4+92854250304440625*u^2-6190283353629375)/sqrt(pi)}
		     else {K <- function(X) (-1)^r * Hermite(X,r) * eval(Kr);
                fx <- K(u)}}
        else if (kernel=="cosine"){
		     if (r%%2==0){
			 if ((r%/%2)%%2==0){
			 K <- function(X)( ((-1)^r/2^(r+2))*pi^(r+1)*cos(0.5*pi*X))* (X >= -1 & X <= 1)}else{
			 K <- function(X)(-((-1)^r/2^(r+2))*pi^(r+1)*cos(0.5*pi*X))* (X >= -1 & X <= 1)}}else{
			 if ((r%/%2)%%2==0){
			 K <- function(X)( ((-1)^r/2^(r+2))*pi^(r+1)*sin(0.5*pi*X))* (X >= -1 & X <= 1)}else{
			 K <- function(X)(-((-1)^r/2^(r+2))*pi^(r+1)*sin(0.5*pi*X))* (X >= -1 & X <= 1)}}
             fx <- K(u)}
        else if (kernel=="tricube"){
             if (r == 1){fx <- (u < 0 & u >= -1) *((70/9)*(u^3+1)^2*u^2)+ (u >= 0 & u <= 1)*(-(70/9)*(u^3-1)^2*u^2)}
		     else if (r == 2){fx <- (u < 0 & u >= -1) *((560/9)*u^7+(700/9)*u^4+(140/9)*u)+ (u >= 0 & u <= 1)*(-(560/9)*u^7+(700/9)*u^4-(140/9)*u)}
		     else if (r == 3){fx <- (u < 0 & u >= -1) *((3920/9)*u^6+(2800/9)*u^3+(140/9))+ (u >= 0 & u <= 1)*(-(3920/9)*u^6+(2800/9)*u^3-(140/9))}
		     else if (r == 4){fx <- (u < 0 & u >= -1) *((7840/3)*u^5+(2800/3)*u^2)+ (u >= 0 & u <= 1)*(-(7840/3)*u^5+(2800/3)*u^2)}
		     else if (r == 5){fx <- (u < 0 & u >= -1) *((39200/3)*u^4+(5600/3)*u)+ (u >= 0 & u <= 1)*(-(39200/3)*u^4+(5600/3)*u)}
		     else if (r == 6){fx <- (u < 0 & u >= -1) *((156800/3)*u^3+(5600/3))+ (u >= 0 & u <= 1)*(-(156800/3)*u^3+(5600/3))}
		     else if (r == 7){fx <- (u <0 & u >= -1) *(156800*u^2)+ (u >= 0 & u <= 1)*(-156800*u^2)}
		     else if (r == 8){fx <- (u < 0 & u >= -1) *(313600*u)+ (u >= 0 & u <= 1)*( -313600*u)}
		     else if (r == 9){fx <- (u < 0 & u >= -1) *(313600)+ (u >= 0 & u <= 1)*( -313600)}}
        else if (kernel=="triweight"){
             if (r == 1){fx <- (-(105/16)*(u^2-1)^2*u)*(u >= -1 & u <= 1)}
		     else if (r == 2){fx <- (-(525/16)*u^4+(315/8)*u^2-(105/16))*(u >= -1 & u <= 1)}
		     else if (r == 3){fx <- (-(525/4)*u^3+(315/4)*u)*(u >= -1 & u <= 1)}
		     else if (r == 4){fx <- (-(1575/4)*u^2+(315/4))*(u >= -1 & u <= 1)}
		     else if (r == 5){fx <- (-(1575/2)*u)*(u >= -1 & u <= 1)}
		     else if (r == 6){fx <- (-1575/2)*(u >= -1 & u <= 1)}}
        else if (kernel=="biweight"){
             if (r == 1){fx <- (((15/4)*(u^2-1))*u)*(u >= -1 & u <= 1)}
		     else if (r == 2){fx <- ((45/4)*u^2-(15/4))*(u >= -1 & u <= 1)}
		     else if (r == 3){fx <- ((45/2)*u)*(u >= -1 & u <= 1)}
		     else if (r == 4){fx <- (45/2)*(u >= -1 & u <= 1)}}
        else if (kernel =="triangular"){
             if (r == 1){fx <- (u <= 0 & u >= -1) - (u >= 0 & u <= 1)}}
        else if (kernel=="epanechnikov"){
             if (r == 1){fx <- (-(3/2)*u)*(u >= -1 & u <= 1)}
		     else if (r == 2){fx <- (-3/2)*(u >= -1 & u <= 1)}}
        else if (kernel=="silverman"){
             if (r == 1){fx <- (u < 0)*((1/4)*sqrt(2)*exp((1/2)*u*sqrt(2))*(sin(-(1/2)*u*sqrt(2)+(1/4)*pi)-cos(-(1/2)*u*sqrt(2)+(1/4)*pi)))+ (u >= 0)* (-(1/4)*sqrt(2)*exp(-(1/2)*u*sqrt(2))*(sin((1/2)*u*sqrt(2)+(1/4)*pi)-cos((1/2)*u*sqrt(2)+(1/4)*pi)))}
		     else if (r == 2){fx <- (u < 0)*(-(1/2)*exp((1/2)*u*sqrt(2))*cos(-(1/2)*u*sqrt(2)+(1/4)*pi))+ (u >= 0)* (-(1/2)*exp(-(1/2)*u*sqrt(2))*cos((1/2)*u*sqrt(2)+(1/4)*pi))}
		     else if (r == 3){fx <- (u < 0)*(-(1/4)*exp((1/2)*u*sqrt(2))*sqrt(2)*(cos(-(1/2)*u*sqrt(2)+(1/4)*pi)+sin(-(1/2)*u*sqrt(2)+(1/4)*pi)))+ (u >= 0)* ((1/4)*exp(-(1/2)*u*sqrt(2))*sqrt(2)*(cos((1/2)*u*sqrt(2)+(1/4)*pi)+sin((1/2)*u*sqrt(2)+(1/4)*pi)))}
		     else if (r == 4){fx <- (u < 0)*(-(1/2)*exp((1/2)*u*sqrt(2))*sin(-(1/2)*u*sqrt(2)+(1/4)*pi))+ (u >= 0)* (-(1/2)*exp(-(1/2)*u*sqrt(2))*sin((1/2)*u*sqrt(2)+(1/4)*pi))}
		     else if (r == 5){fx <- (u < 0)*(-(1/4)*sqrt(2)*exp((1/2)*u*sqrt(2))*(sin(-(1/2)*u*sqrt(2)+(1/4)*pi)-cos(-(1/2)*u*sqrt(2)+(1/4)*pi)))+ (u >= 0)* ((1/4)*sqrt(2)*exp(-(1/2)*u*sqrt(2))*(sin((1/2)*u*sqrt(2)+(1/4)*pi)-cos((1/2)*u*sqrt(2)+(1/4)*pi)))}
		     else if (r == 6){fx <- (u < 0)*((1/2)*exp((1/2)*u*sqrt(2))*cos(-(1/2)*u*sqrt(2)+(1/4)*pi))+ (u >= 0)* ((1/2)*exp(-(1/2)*u*sqrt(2))*cos((1/2)*u*sqrt(2)+(1/4)*pi))}
		     else if (r == 7){fx <- (u < 0)*((1/4)*exp((1/2)*u*sqrt(2))*sqrt(2)*(cos(-(1/2)*u*sqrt(2)+(1/4)*pi)+sin(-(1/2)*u*sqrt(2)+(1/4)*pi)))+ (u >= 0)* (-(1/4)*exp(-(1/2)*u*sqrt(2))*sqrt(2)*(cos((1/2)*u*sqrt(2)+(1/4)*pi)+sin((1/2)*u*sqrt(2)+(1/4)*pi)))}}
               }
return(fx)
}

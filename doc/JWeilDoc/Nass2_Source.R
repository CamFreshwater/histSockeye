## Nass Source file
## Created by AMO May 10, 2015

library(ggplot2); library(plyr); library(reshape2);library(FSA); library(broom);library(bbmle)
library(outliers)
# add in function for standard error
se <- function(x) {sd(x)/ sqrt(length(x))}

# Round function for plotting linear model outputs into graphs

roundUp <- function(x,to=0.001)
{ to*(x%/%to + as.logical(x%%to))
}


# Climate & CONDITION ----------------------------------------

plot.climate <- function(x, npar=TRUE,print=TRUE) { 
  par(mar=c(5, 11, 4, 1) + 0.1)
  plot(x$FultonK_N5 ~ x$ocean.yr, axes=F, ylim=c(0,max(x$FultonK_N5,na.rm = TRUE)), xlab="", ylab="",
       type="p",col="grey70", main="",xlim=c(1911,1944), pch=20)
  lines(x$meanK~ x$ocean.yr,pch=20,col="black",lwd=3,lty=3)
  axis(2, ylim=c(0,max(x$FultonK_N5)),col="black",lwd=2, cex.axis=0.6)
  mtext(2,text="Condition(K)",line=1.9, cex=0.7)
  
  par(new=T)
  plot(x$meanSST ~ x$ocean.yr, axes=F, ylim=c(2,max(x$meanSST)), xlab="", ylab="", 
       type="l", lwd=3, main="",xlim=c(1911,1944), col="coral1")
  axis(2, ylim=c(0,max(x$meanSST)),lwd=2,line=3, cex.axis=0.6, col="coral1")
  #points(x$meanSST ~ x$ocean.yr,pch=17, col="coral1")
  mtext(2,text="Mean SST (C)",line=4.6,col="coral1",cex=0.7)
  
  par(new=T)
  plot(x$avgPDO ~ x$ocean.yr, axes=F, ylim=c(-6,3), xlab="", ylab="", 
       type="l", lwd=3, main="",xlim=c(1911,1944),col="dodgerblue")
  axis(2, ylim=c(0,max(x$avgPDO)),lwd=2,line=5.7,cex.axis=0.6,col="dodgerblue")
  #points(x$avgPDO ~ x$ocean.yr,pch=20)
  mtext(2,text="PDO",line=7.3, col="dodgerblue",cex=0.7)
  
  par(new=T)
  plot(x$ALPI ~ x$ocean.yr, axes=F, ylim=c(-6,20), xlab="", ylab="", 
       type="l", lwd=3, main="",xlim=c(1911,1944),col="limegreen")
  axis(2, ylim=c(0,max(x$ALPI)),lwd=2,line=8.2,cex.axis=0.6,col="limegreen")
  #points(x$avgPDO ~ x$ocean.yr,pch=20)
  mtext(2,text="ALPI",line=9.8, col="limegreen",cex=0.7)
  
  axis(1,pretty(range(x$ocean.yr),10))
  #axis(1, at=seq(1911, 1946, by=1), labels = TRUE)
  mtext("Ocean Entry Year",side=1,col="black",line=2, cex=0.8)
  #mtext(paste(substitute(x), "abundance"),side=3,col="black",line=2, cex=1.8))}
  text(1915, 20, paste(substitute(x)), cex=1.2)}


# Climate & LENGTH ----------------------------------------

plot.Clim.Length <- function(x, npar=TRUE,print=TRUE) { 
  par(mar=c(5, 11, 4, 1) + 0.1)
  plot(x$LENGTH_mm ~ x$ocean.yr, axes=F, ylim=c(0,max(x$LENGTH_mm,na.rm = TRUE)), xlab="", ylab="",
       type="p",col="grey70", main="",xlim=c(1911,1944), pch=20)
  lines(x$meanL~ x$ocean.yr,pch=20,col="black",lwd=3,lty=3)
  axis(2, ylim=c(0,max(x$LENGTH_mm)),col="black",lwd=2, cex.axis=0.6)
  mtext(2,text="Length(mm)",line=1.9, cex=0.7)
  
  par(new=T)
  plot(x$meanSST ~ x$ocean.yr, axes=F, ylim=c(2,max(x$meanSST)), xlab="", ylab="", 
       type="l", lwd=3, main="",xlim=c(1911,1944), col="coral1")
  axis(2, ylim=c(0,max(x$meanSST)),lwd=2,line=3, cex.axis=0.6, col="coral1")
  #points(x$meanSST ~ x$ocean.yr,pch=17, col="coral1")
  mtext(2,text="Mean SST (C)",line=4.6,col="coral1",cex=0.7)
  
  par(new=T)
  plot(x$avgPDO ~ x$ocean.yr, axes=F, ylim=c(-6,3), xlab="", ylab="", 
       type="l", lwd=3, main="",xlim=c(1911,1944),col="dodgerblue")
  axis(2, ylim=c(0,max(x$avgPDO)),lwd=2,line=5.7,cex.axis=0.6,col="dodgerblue")
  #points(x$avgPDO ~ x$ocean.yr,pch=20)
  mtext(2,text="PDO",line=7.3, col="dodgerblue",cex=0.7)
  
  par(new=T)
  plot(x$ALPI ~ x$ocean.yr, axes=F, ylim=c(-6,20), xlab="", ylab="", 
       type="l", lwd=3, main="",xlim=c(1911,1944),col="limegreen")
  axis(2, ylim=c(0,max(x$ALPI)),lwd=2,line=8.2,cex.axis=0.6,col="limegreen")
  #points(x$avgPDO ~ x$ocean.yr,pch=20)
  mtext(2,text="ALPI",line=9.8, col="limegreen",cex=0.7)
  
  axis(1,pretty(range(x$ocean.yr),10))
  #axis(1, at=seq(1911, 1946, by=1), labels = TRUE)
  mtext("Ocean Entry Year",side=1,col="black",line=2, cex=0.8)
  #mtext(paste(substitute(x), "abundance"),side=3,col="black",line=2, cex=1.8))}
  text(1915, 20, paste(substitute(x)), cex=1.2)}

# Climate & WEIGHT ----------------------------------------

plot.Clim.Wgt <- function(x, npar=TRUE,print=TRUE) { 
  par(mar=c(5, 11, 4, 1) + 0.1)
  plot(x$WEIGHT._g ~ x$ocean.yr, axes=F, ylim=c(0,max(x$WEIGHT._g,na.rm = TRUE)), xlab="", ylab="",
       type="p",col="grey70", main="",xlim=c(1911,1944), pch=20)
  lines(x$meanW~ x$ocean.yr,pch=20,col="black",lwd=3,lty=3)
  axis(2, ylim=c(0,max(x$WEIGHT._g)),col="black",lwd=2, cex.axis=0.6)
  mtext(2,text="Weight (g)",line=1.9, cex=0.7)
  
  par(new=T)
  plot(x$meanSST ~ x$ocean.yr, axes=F, ylim=c(2,max(x$meanSST)), xlab="", ylab="", 
       type="l", lwd=3, main="",xlim=c(1911,1944), col="coral1")
  axis(2, ylim=c(0,max(x$meanSST)),lwd=2,line=3, cex.axis=0.6, col="coral1")
  #points(x$meanSST ~ x$ocean.yr,pch=17, col="coral1")
  mtext(2,text="Mean SST (C)",line=4.6,col="coral1",cex=0.7)
  
  par(new=T)
  plot(x$avgPDO ~ x$ocean.yr, axes=F, ylim=c(-6,3), xlab="", ylab="", 
       type="l", lwd=3, main="",xlim=c(1911,1944),col="dodgerblue")
  axis(2, ylim=c(0,max(x$avgPDO)),lwd=2,line=5.7,cex.axis=0.6,col="dodgerblue")
  #points(x$avgPDO ~ x$ocean.yr,pch=20)
  mtext(2,text="PDO",line=7.3, col="dodgerblue",cex=0.7)
  
  par(new=T)
  plot(x$ALPI ~ x$ocean.yr, axes=F, ylim=c(-6,20), xlab="", ylab="", 
       type="l", lwd=3, main="",xlim=c(1911,1944),col="limegreen")
  axis(2, ylim=c(0,max(x$ALPI)),lwd=2,line=8.2,cex.axis=0.6,col="limegreen")
  #points(x$avgPDO ~ x$ocean.yr,pch=20)
  mtext(2,text="ALPI",line=9.8, col="limegreen",cex=0.7)
  
  axis(1,pretty(range(x$ocean.yr),10))
  #axis(1, at=seq(1911, 1946, by=1), labels = TRUE)
  mtext("Ocean Entry Year",side=1,col="black",line=2, cex=0.8)
  #mtext(paste(substitute(x), "abundance"),side=3,col="black",line=2, cex=1.8))}
  text(1915, 20, paste(substitute(x)), cex=1.2)}

# PLOT REGRESSION of sockeye ages by: 
# Add mtext of slope, pval, etc to graphs (all graphs!)

# Avg climate
plot.avgclim<- function(x, npar=TRUE,print=TRUE) {
par(mfrow = c(2, 2))
plot(Nass1.2$LWcondition ~ Nass1.2$avgPDO)
mod1 <- lm(Nass1.2$LWcondition ~ Nass1.2$avgPDO)
abline(mod1, col="blue",lty=2)
tidy(mod1)
plot(Nass1.2$LWcondition ~ Nass1.2$avgALPI)
mod2 <- lm(Nass1.2$LWcondition ~ Nass1.2$avgALPI)
abline(mod2, col="blue",lty=2)
tidy(mod2)
plot(Nass1.2$LWcondition ~ Nass1.2$avgSSTn)
mod3 <- lm(Nass1.2$LWcondition ~ Nass1.2$avgSSTn)
abline(mod3, col="blue",lty=2)
tidy(mod3)
plot(Nass1.2$LWcondition ~ Nass1.2$avgSSTb)
mod4 <- lm(Nass1.2$LWcondition ~ Nass1.2$avgSSTb)
abline(mod4, col="blue",lty=2)
tidy(mod4)
}


# Calculate the slopes of all the within year conditions

plot.intra <- function(x) { 
  uniq <- unique(unlist(x$ocean.yr))
  slopes <- matrix(NA, nrow=length(uniq), ncol=2)
  
  for (i in 1:length(uniq)){
    data_1 <- subset(x, x== uniq[i])
    mod <- lm(data_1$FultonK_N5 ~ data_1$WEEK)
    b <- mod$coefficients[2]
    slopes[i,1] <- data_1$ocean.yr[1]
    slopes[i,2] <- mod$coefficients[2]}
  print(slopes)
}

# HISTORICAL SOCKEYE FUNCTIONS
# Created by C Freshwater APRIL 26, 2018
# Last revised: ONGOING
# Functions for fullSockeyeModels.R
# -------------------------------------------------


## This model standardizes explanatory variables in the main dataset
standardizeVar <- function(x){
	x$pc2Std <- (x$pc2 - mean(x$pc2))/sd(x$pc2)
	x$tempStd <- (x$temp - mean(x$temp))/sd(x$temp)
	x$pdoStd <- (x$pdo - mean(x$pdo))/sd(x$pdo)
	x$alpiStd <- (x$alpi - mean(x$alpi))/sd(x$alpi)
	x$sockStd <- (x$sockCatch - mean(x$sockCatch))/sd(x$sockCatch)
	x$pinkStd <- (x$pinkCatch - mean(x$pinkCatch))/sd(x$pinkCatch)
	x$totalStd <- (x$totalCatch - mean(x$totalCatch))/sd(x$totalCatch)
	return(x)
}
# ------------------------------------------------------------------------


## This model plots TS of response and std exp variables for each dataset
plotTS <- function(meanDat, dat){
	temp <- data.frame(retYr = seq(from=min(dat$retYr), to=max(dat$retYr), by=1),
		yrFac = as.factor(seq(from=min(dat$retYr), to=max(dat$retYr), by=1)))
	temp$meanFL <- meanDat$meanFL[match(temp$yrFac, meanDat$yrFac)]
	temp$meanP1 <- meanDat$meanP1[match(temp$yrFac, meanDat$yrFac)]
	temp$meanP2 <- meanDat$meanP2[match(temp$yrFac, meanDat$yrFac)]
	figTitle <- paste(unique(dat$watershed), unique(dat$dataSet), sep=" ")

	x <- plot(temp$meanFL ~ temp$retYr, type="l", axes=FALSE, lwd=1.5, xlab="", ylab="")
	axis(1, tick=T, at=pretty(c(seq(from=min(temp$retYr), to=max(temp$retYr), by=5)), n=4))
	axis(2, tick=T, at=pretty(c(seq(from=min(temp$meanFL, na.rm=TRUE), to=max(temp$meanFL, na.rm=TRUE), by=5)), n=4))
	mtext(side=2, line=2.5, 'Mean Fork Length', cex=1.2)
	mtext(side=3, line=1.25, figTitle, cex=1.2)
	par(new=TRUE)
	x <- plot(temp$meanP1 ~ temp$retYr, type="l", col="#1b9e77", axes=F, xlab="", ylab="", lty=2, lwd=1.25,
		ylim=c(min(temp$meanP1, na.rm=T), max(temp$meanP1, na.rm=T)))
	lines(temp$meanP2~ temp$retYr, type="l", col="#7570b3", xlab="", ylab="", lty=2, lwd=1.25,
		ylim=c(min(temp$meanP2, na.rm=T), max(temp$meanP2, na.rm=T)))
	return(x)
}
# ------------------------------------------------------------------------

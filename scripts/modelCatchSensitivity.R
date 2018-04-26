# HISTORICAL SOCKEYE ANALYSIS - SENSITIVITY ANALYSIS
# Created by C Freshwater Apr 22, 2018
# Last revised: ONGOING
# Replicates analyses conducted in fullSockeyeModels.R, but excludes fish data collected
# pre-1925 to test relative influence of including AK/WA catches on model selection and parameter
# -------------------------------------------------

setwd("/Users/cam/github")

library(mgcv); library(dplyr); library(ggplot2); library(reshape2); library(here); 
library(MuMIn); library(corrplot); library(car); library(mgcv.helper)


sstPca <- read.table(here("github/histSockeye/data/sstPCA.txt")) #principal components of SST variation in NE Pacific (170E to 240E, 40N-65N)
sstRaw <- read.table(here("github/histSockeye/data/sstRaw.txt")) # pacific ocean SST
pdo <- read.csv(here("github/histSockeye/data/pdo.csv"), stringsAsFactors=F) 
meanAlpi <- read.csv(here("github/histSockeye/data/alpi.csv"), stringsAsFactors=F) #stops at 2015
sockDat <- read.csv(here("github/histSockeye/data/nassFullSox.csv"), stringsAsFactors=F)
akCatchDat <- read.csv(here("github/histSockeye/data/akCatch.csv"), stringsAsFactors=F)
trimTotCatch <- read.csv(here("github/histSockeye/data/cleanAgTrimCatch.csv"), stringsAsFactors=F)

## ---------------------- Clean ------------------------------
### sst PCA
colnames(sstPca) <- c("id", "retYr", "month", "pc1", "pc2", "pc3", "pc4", "pc5")
months <- c("3","4","5","6")
sstPca <- sstPca[sstPca$month %in% months, c("retYr","month","pc2")] 
meanPca <- sstPca %>%
	group_by(retYr) %>%
	summarise(pc2=mean(pc2))

### sst raw
colnames(sstRaw) <- c("long", "lat", "retYr", "month", "temp")
months <- c("3","4","5","6")
sstRaw <- sstRaw[sstRaw$month %in% months, c("retYr","month","temp")] 
meanSst <- sstRaw %>%
	group_by(retYr) %>%
	summarise(temp=mean(temp))

### pdo
meanPdo <- data.frame(retYr=pdo$YEAR,
					  pdo=apply(pdo[,c("MAR","APR","MAY","JUN")], 1, mean)
					  )

### alpi
names(meanAlpi)[1:2] <- c("retYr", "alpi")

### catch
trimAkCatch <- akCatchDat[akCatchDat$year > 1920, ]

## AK catch
sockAkCatch <- trimAkCatch[trimAkCatch$species=="Sockeye", -c(1,2,4)]
names(sockAkCatch)[c(1,2)] <- c("retYr", "sockCatch")
pinkAkCatch <- trimAkCatch[trimAkCatch$species=="Pink", -c(1,2,4)]
names(pinkAkCatch)[c(1,2)] <- c("retYr", "pinkCatch")
totalAkCatch <- data.frame(retYr = sockAkCatch$retYr, totalCatch = (sockAkCatch$sockCatch + pinkAkCatch$pinkCatch))

## full catch
sockTotCatch <- trimTotCatch[trimCatch$species=="Sockeye", -c(1,2,4)]
names(sockTotCatch)[c(1,2)] <- c("retYr", "sockCatch")
pinkTotCatch <- trimCatch[trimCatch$species=="Pink", -c(1,2,4)]
names(pinkTotCatch)[c(1,2)] <- c("retYr", "pinkCatch")
totalCatch <- data.frame(retYr = sockTotCatch$retYr, totalCatch = (sockTotCatch$sockCatch + pinkTotCatch$pinkCatch))

### merge sox data w/ environmental
#Ak data
# fullDat <- Reduce(function(x, y) merge(x, y, by=c("retYr")), list(sockDat, meanPdo, meanSst, meanPca, meanAlpi, sockAkCatch, pinkAkCatch, totalAkCatch))
#total data
fullDat <- Reduce(function(x, y) merge(x, y, by=c("retYr")), list(sockDat, meanPdo, meanSst, meanPca, meanAlpi, sockAkCatch, pinkAkCatch, totalAkCatch))

nassDat <- subset(fullDat, fullDat$watershed %in% "nass")
nassDatMod <- subset(nassDat, nassDat$dataSet %in% "mod")
nassDat <- subset(nassDat, nassDat$dataSet %in% "hist")
riversDat <- subset(fullDat,fullDat$watershed %in% "rivers")

nassDat$age <- factor(nassDat$age)
nassDat$yrFac <- factor(nassDat$retYr)
nassDatMod$age <- factor(nassDatMod$age)
nassDatMod$yrFac <- factor(nassDatMod$retYr)
riversDat$age <- factor(riversDat$age)
riversDat$yrFac <- factor(riversDat$retYr)

datList <- list(nassDat, nassDatMod, riversDat)

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

datListN <- lapply(datList, function(x) standardizeVar(x))
names(datListN) <- c("nassDat", "nassDatMod", "riversDat")

# correlations among predictor variables
makeCorPlot <- function(df){
	mat <- cor(df[,c(21:27)])
	figTitle <- paste(unique(df$watershed), unique(df$dataSet), sep=" ")
	corrplot.mixed(mat, lower="ellipse", upper="number")
	mtext(side=3, line=1.5, figTitle, cex=1.2)
}

pdf(here("github/histSockeye/outputs/figs/corrPlot.pdf"), height=10, width=10)
par(mfrow=c(2,2), mar=c(0,0,2.75,0)+0.1, oma=c(0,0,0,0))
sapply(datListN, function(x) makeCorPlot(x))
dev.off()

# -----------------------------------------------
# meanDat <- fullDat %>%
# 	group_by(yrFac, age, watershed, dataSet) %>%
# 	summarize(meanFL = mean(fl), pdo = mean(pdo), alpi = mean(alpi), 
# 		rawSst = mean(temp), pcaSst = mean(pc2), pink = mean(pinkCatch), 
# 		sox = mean(sockCatch))

# ## Changes in length through time
# ggplot(meanDat, aes(x = as.numeric(yrFac), y = meanFL)) + 
#     geom_line() + 
#     facet_wrap(~ watershed)


# ------------------------------------------------------
## Full model comparison with different environmental covariates; removed AR1 terms because they don't seem to be doing anything
for(i in seq_along(datListN)){
	dataset <- datListN[[i]]
	dataset$dum <- 1
	
	# models
	null <- gam(fl ~ age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	pdo <- gam(fl ~ s(pdoStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	temp <- gam(fl ~ s(tempStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	pc2 <- gam(fl ~ s(pc2Std, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	alpi <- gam(fl ~ s(alpiStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	pink <- gam(fl ~ s(pinkStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	sock <- gam(fl ~ s(sockStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	total <- gam(fl ~ s(totalStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	pdoPink <- gam(fl ~ s(pdoStd, by=age, k=3) + s(pinkStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	tempPink <- gam(fl ~ s(tempStd, by=age, k=3) + s(pinkStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	pc2Pink <- gam(fl ~ s(pc2Std, by=age, k=3) + s(pinkStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	alpiPink <- gam(fl ~ s(alpiStd, by=age, k=3) + s(pinkStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	pdoSock <- gam(fl ~ s(pdoStd, by=age, k=3) + s(sockStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	tempSock <- gam(fl ~ s(tempStd, by=age, k=3) + s(sockStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	pc2Sock <- gam(fl ~ s(pc2Std, by=age, k=3) + s(sockStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	alpiSock <- gam(fl ~ s(alpiStd, by=age, k=3) + s(sockStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	pdoTotal <- gam(fl ~ s(pdoStd, by=age, k=3) + s(totalStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	tempTotal <- gam(fl ~ s(tempStd, by=age, k=3) + s(totalStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	pc2Total <- gam(fl ~ s(pc2Std, by=age, k=3) + s(totalStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)
	alpiTotal <- gam(fl ~ s(alpiStd, by=age, k=3) + s(totalStd, by=age, k=3) + age + s(yrFac, bs="re", by=dum), method="ML", data=dataset)

	modOutputs <- list(null, pdo, temp, pc2, alpi, pink, sock, total, pdoPink, tempPink, pc2Pink, alpiPink, 
		pdoSock, tempSock, pc2Sock, alpiSock, pdoTotal, tempTotal, pc2Total, alpiTotal)
	names(modOutputs)[1:length(modOutputs)] <- c("null", "pdo", "temp", "pc2", "alpi", "pink", "sock", "total", "pdoPink", "tempPink", "pc2Pink", "alpiPink", 
		"pdoSock", "tempSock", "pc2Sock", "alpiSock", "pdoTotal", "tempTotal", "pc2Total", "alpiTotal")
	assign(paste("fitsAk", dataset$watershed[1], dataset$dataSet[1], sep="_"), modOutputs)
	modRankings <- AICc(null, pdo, temp, pc2, alpi, pink, sock, total, pdoPink, tempPink, pc2Pink, alpiPink, 
		pdoSock, tempSock, pc2Sock, alpiSock, pdoTotal, tempTotal, pc2Total, alpiTotal)
	modRankings <- modRankings[order(modRankings$AICc),]
	print(modRankings)
	assign(paste("rankingAk", dataset$watershed[1], dataset$dataSet[1], sep="_"), modRankings)
}

# save data files
saveRDS(fitsAk_nass_mod, here("github/histSockeye/outputs/data/sensAnalysis/nassModernModFitsAk.rds"))
saveRDS(fitsAk_nass_hist, here("github/histSockeye/outputs/data/sensAnalysis/nassHistoricModFitsAk.rds"))
saveRDS(fitsAk_rivers_hist, here("github/histSockeye/outputs/data/sensAnalysis/riversHistoricModFitsAgg.rds"))

# read in data files
fits_nass_mod <- readRDS(here("github/histSockeye/outputs/data/nassModernModFits.rds"))
fits_nass_hist <- readRDS(here("github/histSockeye/outputs/data/nassHistoricModFits.rds"))
fits_rivers_hist <- readRDS(here("github/histSockeye/outputs/data/riversHistoricModFits.rds"))

# extract top models 
nmPc2Pink <- fitsAgg_nass_mod$pc2Pink
nmPc2PinkAk <- fitsAk_nass_mod$pc2Pink
nhTempPink <- fitsAgg_nass_hist$tempPink
nhTempPinkAk <- fitsAk_nass_hist$tempPink
rhPdoPink <- fitsAgg_rivers_hist$pdoPink
rhPdoPinkAk <- fitsAk_rivers_hist$pdoPink

# NASS SOCKEYE DATA - CLEMENS & GILBERT
# REGRESSION FILE
# R File Created by AMO May15, 2015
# Last updated: Oct 17, 2016
# -------------------------------------------------

## Updated and compiled by Jacob W Dec 18, 2017 ##
## Updated by Cameron F Jan 3, 2018 ##


# Run lines 13 - 155 to clean data; then 218 - 247 to clean data for modelling
# Run 259; and then run the subset of models a) Nass2.2 OR b) Rivers1.2 OR c) Rivers1.3
rm(list = ls()) # clear objects 

setwd("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/NassRivers")
# setwd("/Users/JacobWeil/Desktop/R/Nass_2017") # change working directory
source("Nass2_Source.R")

data<-read.csv("nass_data.csv", stringsAsFactors=F) # in csv, changed WEIGHT._g to WEIGHT_g AND BROOD_YR to BROODYR
rivers<-read.csv("Rivers.Sockeye.csv", stringsAsFactors=F) # in csv, changed age to AGE

nSST<-read.csv("sst_narrow.csv", stringsAsFactors=F) # regional scale SST
bSST<-read.csv("sst_broad.csv", stringsAsFactors=F) # pacific ocean SST
PDO<-read.csv("pdo_historical_data.csv", stringsAsFactors=F) 
ALPI<-read.csv("alpi_historical_data.csv", stringsAsFactors=F)
pink<-read.csv("Pink_Abundance.csv", stringsAsFactors=F, header=TRUE)

head(data)
str(data)

## ---------------------- Clean Up Nass Data ------------------------------
# Omit years with little to no collections & weeks with irregular sampling
nass<-data[!data$YEAR == 1915,]
nass<-nass[!nass$YEAR == 1938,]
nass<-nass[!nass$WEEK == 25,]

# Keep dominant age classes
nass <- subset(nass,nass$AGE %in% c(1.2,1.3,2.2,2.3))
nass$watershed <- "nass"

# Re-calculate condition so units match between watersheds
nass$L3<-nass$LENGTH_mm^3
nass$W3<-nass$WEIGHT_g*1000
nass$K<-(nass$W3/nass$L3)

# ----------------------- Clean up Rivers Data ----------------------------

# Subset dominant age classes
rivers <- subset(rivers,rivers$AGE %in% c(1.2,1.3))
rivers$watershed <-"rivers"

# Condition
rivers$L3<-rivers$LENGTH_mm^3
rivers$W3<-rivers$WEIGHT_g*1000
rivers$K<-rivers$W3/rivers$L3
head(rivers)

# ---------------------- Merge dataframes ---------------------------

sock <- rbind(nass,rivers)
names(sock)[1]<-paste(c("retYr")) #return.yr throws errors at some points because return is a base function
names(pink)[1]<-paste(c("retYr"))

# Calculate Ocean Entry Yr.
sock$seatime <- ifelse(sock$AGE == 1.2 | sock$AGE == 2.2, 2,
                ifelse(sock$AGE == 1.3 | sock$AGE == 2.3, 3,NA))
sock$entry.yr <- sock$retYr - sock$seatime

sock$intyr1 <- sock$entry.yr + 1
sock$intyr2 <- ifelse(sock$AGE == 1.3 | sock$AGE == 2.3, sock$entry.yr + 2, NA)

## ---------------------- Clean Up SST ------------------------------

# Use months March - June from 1911 to 1946 and get the avg over the 4 mnths for each year
bSST <- ddply(bSST, "year",summarize, meanSST = mean(temp))
nSST <- ddply(nSST, "year",summarize, meanSST = mean(temp))

# rename to match up with ocean entry, brood, and return years
nSSTentry <- plyr::rename(nSST, c("year"="entry.yr", "meanSST"="nSSTentry"))
nSSTint1 <- plyr::rename(nSST, c("year"="intyr1", "meanSST"="nSSTint1"))
nSSTint2 <- plyr::rename(nSST, c("year"="intyr2", "meanSST"="nSSTint2"))
nSSTint2<-nSSTint2[3:nrow(nSSTint2),]
nSSTreturn <- plyr::rename(nSST, c("year"="retYr", "meanSST"="nSSTreturn"))

bSSTentry <- plyr::rename(bSST, c("year"="entry.yr", "meanSST"="bSSTentry"))
bSSTint1 <- plyr::rename(bSST, c("year"="intyr1", "meanSST"="bSSTint1"))
bSSTint2 <- plyr::rename(bSST, c("year"="intyr2", "meanSST"="bSSTint2"))
bSSTint2<-bSSTint2[3:nrow(bSSTint2),]
bSSTreturn <- plyr::rename(bSST, c("year"="retYr", "meanSST"="bSSTreturn"))

## ---------------------- Clean Up PDO ------------------------------
# Take months March - June from 1911 - 1946
PDO$annualPDO <- PDO$Mean # describe the mean as annual mean over 12 months
PDO$avgPDO <- rowMeans(PDO[,3:6]) # only take mean from Mar - Jun
PDO <- subset(PDO,select=c("Year","avgPDO"))

# rename to match up with ocean entry, brood, and return years
PDOentry <- plyr::rename(PDO, c("Year"="entry.yr", "avgPDO"="PDOentry"))
PDOint1 <- plyr::rename(PDO, c("Year"="intyr1", "avgPDO"="PDOint1"))
PDOint2 <- plyr::rename(PDO, c("Year"="intyr2", "avgPDO"="PDOint2"))
PDOint2<-na.omit(PDOint2); PDOint2<-PDOint2[3:nrow(PDOint2),]
PDOreturn <- plyr::rename(PDO, c("Year"="retYr", "avgPDO"="PDOreturn"))
PDOreturn<-na.omit(PDOreturn)

## ---------------------- Clean Up ALPI ------------------------------
# Subset years 1911 - 1946
ALPI<- ALPI[12:47,]

# rename to match up with ocean entry, brood, and return years
ALPIentry <- plyr::rename(ALPI, c("year"="entry.yr", "ALPI"="ALPIentry"))
ALPIint1 <- plyr::rename(ALPI, c("year"="intyr1", "ALPI"="ALPIint1"))
ALPIint2 <- plyr::rename(ALPI, c("year"="intyr2", "ALPI"="ALPIint2"))
ALPIint2<-ALPIint2[3:nrow(ALPIint2),]
ALPIreturn <- plyr::rename(ALPI, c("year"="retYr", "ALPI"="ALPIreturn"))

## ---------------------- MERGE SALMON W/ CLIMATE DATA ------------------------------

# Merge w/ narrow SST data
sockeye <-merge(sock, nSSTentry, by="entry.yr")
sockeye <-merge(sockeye, nSSTint1, by="intyr1") 
sockeye <-merge(sockeye, nSSTint2, by="intyr2",all=TRUE) 
sockeye <-merge(sockeye, nSSTreturn, by="retYr") 

# Merge w/ broad SST data
sockeye <-merge(sockeye, bSSTentry, by="entry.yr") 
sockeye <-merge(sockeye, bSSTint1, by="intyr1") 
sockeye <-merge(sockeye, bSSTint2, by="intyr2",all=TRUE) 
sockeye <-merge(sockeye, bSSTreturn, by="retYr") 

# Merge w/ PDO data
sockeye <-merge(sockeye, PDOentry, by="entry.yr") 
sockeye <-merge(sockeye, PDOint1, by="intyr1") 
sockeye <-merge(sockeye, PDOint2, by="intyr2",all=TRUE) 
sockeye <-merge(sockeye, PDOreturn, by="retYr") 

# Merge w/ SST data
sockeye <-merge(sockeye, ALPIentry, by="entry.yr") 
sockeye <-merge(sockeye, ALPIint1, by="intyr1") 
sockeye <-merge(sockeye, ALPIint2, by="intyr2",all=TRUE) 
sockeye <-merge(sockeye, ALPIreturn, by="retYr") 

# Merge w/ Pink Data
sockeye<-merge(sockeye, pink, by ="retYr")


# Calculate avg SST, PDO, and ALPI experienced throughout life

sockeye$avgSSTn <- ifelse(sockeye$AGE == 1.2 | sockeye$AGE == 2.2, rowMeans(sockeye[,c("nSSTentry","nSSTint1","nSSTreturn")]) ,
                          ifelse(sockeye$AGE == 1.3 | sockeye$AGE == 2.3, rowMeans(sockeye[,c("nSSTentry","nSSTint1","nSSTint2","nSSTreturn")]),NA))
sockeye$avgSSTb <- ifelse(sockeye$AGE == 1.2 | sockeye$AGE == 2.2, rowMeans(sockeye[,c("bSSTentry","bSSTint1","bSSTreturn")]) ,
                          ifelse(sockeye$AGE == 1.3 | sockeye$AGE == 2.3, rowMeans(sockeye[,c("bSSTentry","bSSTint1","bSSTint2","bSSTreturn")]),NA))
sockeye$avgPDO <- ifelse(sockeye$AGE == 1.2 | sockeye$AGE == 2.2, rowMeans(sockeye[,c("PDOentry","PDOint1","PDOreturn")]) ,
                         ifelse(sockeye$AGE == 1.3 | sockeye$AGE == 2.3, rowMeans(sockeye[,c("PDOentry","PDOint1","PDOint2","PDOreturn")]),NA))
sockeye$avgALPI <- ifelse(sockeye$AGE == 1.2 | sockeye$AGE == 2.2, rowMeans(sockeye[,c("ALPIentry","ALPIint1","ALPIreturn")]) ,
                          ifelse(sockeye$AGE == 1.3 | sockeye$AGE == 2.3, rowMeans(sockeye[,c("ALPIentry","ALPIint1","ALPIint2","ALPIreturn")]),NA))

sock1 <- subset(sockeye,!is.na(sockeye$WEIGHT_g) & !is.na(sockeye$LENGTH_mm))


# Rename retYr
sock1$retYr <- as.factor(sock1$retYr)

# ------- LW Condition by Resids from Regression ----------------

# plot raw data
par(mfrow = c(1, 1))
mod <- lm(sock1$WEIGHT_g ~ sock1$LENGTH_mm)
summary(mod) # pvalue = <2.2e-16, R2 = 0.68 and b = 1.2
plot(sock1$WEIGHT_g ~ sock1$LENGTH_mm)
# Add residual to salmon
sock1$LWcondition <-resid(mod)

# ----------------------Check for autocorrelation---------------------------

#nSST 
par(mfrow = c(1, 1))
plot(nSST$meanSST ~ nSST$year)
mdl1 <- lm(meanSST ~ year, data=nSST)
summary(mdl1)
par(mfrow=c(2,2))
plot(mdl1)
plot(residuals(mdl1),type="b")
abline(h=0,lty=3)
acf(residuals(mdl1))  # lines do not cross the blue 95% CI lines
acf(nSST$meanSST)

#bSST 
par(mfrow = c(1, 1))
plot(bSST$year,bSST$meanSST)
mdl2 <- lm(meanSST ~ year,data=bSST)
summary(mdl2)
par(mfrow=c(2,2))
plot(mdl2)
par(mfrow=c(1,1))
plot(residuals(mdl2),type="b")
abline(h=0,lty=3)
acf(residuals(mdl2))  # lines do not cross the blue 95% CI lines
acf(bSST$meanSST)

#PDO 
par(mfrow = c(1, 1))
PDO<- na.omit(PDO)
plot(PDO$Year,PDO$avgPDO)
mdl3 <- lm(avgPDO ~ Year,data=PDO)
summary(mdl3)
par(mfrow=c(2,2)); plot(mdl3)
par(mfrow=c(1,1))
plot(residuals(mdl3),type="b")
abline(h=0,lty=3)
acf(residuals(mdl3))  # lines do not cross the blue 95% CI lines
acf(PDO$avgPDO)

#ALPI 
plot(ALPI$year,ALPI$ALPI)
mdl4 <- lm(ALPI ~ year,data=ALPI)
summary(mdl4)
par(mfrow=c(2,2)); plot(mdl4)
par(mfrow=c(1,1))
plot(residuals(mdl4),type="b")
abline(h=0,lty=3)
acf(residuals(mdl4))  # lines do not cross the blue 95% CI lines
acf(ALPI$ALPI)

# Check for covariation among fixed effects
par(mfrow = c(1, 1))
plot(pink_SE~ nSSTreturn , data = sock1)
summary(lm(pink_SE~ nSSTreturn , data = sock1))
sock1$PINKnSSTret_rd <-resid(lm(pink_SE~ nSSTreturn, data = sock1))

plot(pink_SE~ bSSTreturn, data = sock1)
summary(lm(pink_SE~ bSSTreturn , data = sock1))
sock1$PINKbSSTret_rd <-resid(lm(pink_SE~ bSSTreturn, data = sock1))

plot(pink_SE~ ALPIreturn, data = sock1)
summary(lm(pink_SE~ ALPIreturn , data = sock1))
sock1$PINKalpiret_rd <-resid(lm(pink_SE~ ALPIreturn, data = sock1))

plot(pink_SE~ PDOreturn, data = sock1)
summary(lm(pink_SE~ PDOreturn , data = sock1))
sock1$PINKpdoret_rd <-resid(lm(pink_SE~ PDOreturn, data = sock1))

# ----------------  Look at distribution of response variables-------------------- #
# Data
par(mfrow = c(2, 2))
hist(sock1$K, 50,col ="orange", main="Fulton K")
hist(sock1$L3, 50,col ="orange", main = "Length (mm)")
hist(sock1$W3, 50,col ="orange", main = "Weight (g)")
hist(sock1$LWcondition, 50,col ="orange", main="LW Residuals") # skewed

# -----------------------------------------------
# -------------------- Standardize Predictors for modelling ---------------------------#

#change AGE and watershed to factor
sock1$AGE <-as.factor(sock1$AGE)
sock1$AGE <- relevel(sock1$AGE, ref="1.2") # Set dummy ref lev?
sock1$watershed <-as.factor(sock1$watershed)

#standardizing data 
#sock1$avgALPI<-(sock1$avgALPI - mean(sock1$avgALPI))/(sd(sock1$avgALPI))
#sock1$avgSSTn<-(sock1$avgSSTn - mean(sock1$avgSSTn))/(sd(sock1$avgSSTn))
#sock1$avgSSTb<-(sock1$avgSSTb - mean(sock1$avgSSTb))/(sd(sock1$avgSSTb))

#sock1$PDOentry<-(sock1$PDOentry - mean(sock1$PDOentry))/(sd(sock1$PDOentry))
#sock1$ALPIentry<-(sock1$ALPIentry - mean(sock1$ALPIentry))/(sd(sock1$ALPIentry))
#sock1$avgSSTn<-(sock1$avgSSTn - mean(sock1$avgSSTn))/(sd(sock1$avgSSTn))
#sock1$avgSSTb<-(sock1$avgSSTb - mean(sock1$avgSSTb))/(sd(sock1$avgSSTb))

sock1$PDOreturn1<-(sock1$PDOreturn - mean(sock1$PDOreturn))/(sd(sock1$PDOreturn))
sock1$ALPIreturn1<-(sock1$ALPIreturn - mean(sock1$ALPIreturn))/(sd(sock1$ALPIreturn))
sock1$nSSTreturn1<-(sock1$nSSTreturn - mean(sock1$nSSTreturn))/(sd(sock1$nSSTreturn))
sock1$bSSTreturn1<-(sock1$bSSTreturn - mean(sock1$bSSTreturn))/(sd(sock1$bSSTreturn))

sock1$PINKnSSTret_rd <- (sock1$PINKnSSTret_rd - mean(sock1$PINKnSSTret_rd))/(sd(sock1$PINKnSSTret_rd))
sock1$PINKbSSTret_rd <-(sock1$PINKbSSTret_rd - mean(sock1$PINKbSSTret_rd))/(sd(sock1$PINKbSSTret_rd))
sock1$PINKalpiret_rd <- (sock1$PINKalpiret_rd - mean(sock1$PINKalpiret_rd))/(sd(sock1$PINKalpiret_rd))
sock1$PINKpdoret_rd <-(sock1$PINKpdoret_rd - mean(sock1$PINKpdoret_rd))/(sd(sock1$PINKpdoret_rd))


# -------------------- Subset by Watershed and age class ---------------------------#

Nass <- subset(sock1,sock1$watershed %in% "nass")
#Nass1.2 <- subset(Nass,Nass$AGE %in% "1.2")
#Nass1.3 <- subset(Nass,Nass$AGE %in% "1.3")
#Nass2.2 <- subset(Nass,Nass$AGE %in% "2.2")
#Nass2.3 <- subset(Nass,Nass$AGE %in% "2.3")

Rivers <- subset(sock1,sock1$watershed %in% "rivers")
#Rivers1.2 <- subset(Rivers,Rivers$AGE %in% "1.2")
#Rivers1.3 <- subset(Rivers,Rivers$AGE %in% "1.3")

write.csv(sock1, "/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/NassRivers/histSockDat.csv", row.names=FALSE)
write.csv(Nass, "/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/NassRivers/nassDat.csv", row.names=FALSE)
write.csv(Rivers, "/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/NassRivers/riversDat.csv", row.names=FALSE)


#  ---------------- Is there colinearity in the predictors after standardizing? ------------------- #
library(usdm)

# Use VIF <5 as acceptable correlation, as per Jeffrey 2016 and other sources

#PDO and ALPI, VIF = 3.18, FINE
#df <- data.frame(sockeye$PDOreturn, sock1$return); vif(df); v1 <- vifcor(df, th=0.9) # GOOD VIF = 2.8

# PDO and SST, VIF = 2.65 and 1.96, TECHNICALLY NOT FINE, But won't be same model
#df <- data.frame(sock1$PDOentry, sock1$nSSTentry); vif(df); v1 <- vifcor(df, th=0.9) # BAD VIF = 5.7 ???
#df <- data.frame(sock1$PDOentry, sock1$bSSTentry); vif(df); v1 <- vifcor(df, th=0.9) # GOOD VIF = 2.5

# SE Pink VIF fine
#df <- data.frame(sock1$pink_SE, sock1$nSSTreturn); vif(df); v1 <- vifcor(df, th=0.9) # GOOD VIF = 1.328
#df <- data.frame(sock1$pink_SE, sock1$bSSTreturn); vif(df); v1 <- vifcor(df, th=0.9) # GOOD VIF = 1.378
#df <- data.frame(sock1$pink_SE, sock1$ALPIreturn); vif(df); v1 <- vifcor(df, th=0.9) # GOOD VIF = 1.107
#df <- data.frame(sock1$pink_SE, sock1$PDOreturn); vif(df); v1 <- vifcor(df, th=0.9) # GOOD VIF = 1.236

# CEN Pink VIF fine
#df <- data.frame(sock1$pink_CEN, sock1$nSSTreturn); vif(df); v1 <- vifcor(df, th=0.9) # GOOD VIF = 1.233
#df <- data.frame(sock1$pink_CEN, sock1$bSSTreturn); vif(df); v1 <- vifcor(df, th=0.9) # GOOD VIF = 1.404
#df <- data.frame(sock1$pink_CEN, sock1$ALPIreturn); vif(df); v1 <- vifcor(df, th=0.9) # GOOD VIF = 1.109
#df <- data.frame(sock1$pink_CEN, sock1$PDOreturn); vif(df); v1 <- vifcor(df, th=0.9) # GOOD VIF = 1.183



# -----------------------------------------------
# -----------------------------------------------
### ----------------GAMMs NASS-----------------###
# -----------------------------------------------
# -----------------------------------------------


library(mgcv); library(nlme);library(usdm);library(MuMIn); library(itsadug)


### Fulton's K Models ------------------------------- 
### Models labeled m1,m2...mn correspond to Fultonk K as Response Variable

# m1 <- gamm(K ~ s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|as.factor(retYr)))
# m2 <- gamm(K ~ s(PDOreturn, k=3), data = Nass, fx=TRUE,correlation = corAR1(form = ~ 1|retYr))
# m3 <- gamm(K ~ s(nSSTreturn, k=3), data = Nass,fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# m4 <- gamm(K ~ s(bSSTreturn, k=3), data = Nass, fx=TRUE,correlation = corAR1(form = ~ 1|retYr))
# m5 <- gamm(K ~ s(ALPIreturn, k=3), data = Nass,fx=TRUE, correlation = corAR1(form = ~ 1|retYr))

# m6 <- gamm(K ~ s(PDOreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))
# m7 <- gamm(K ~ s(nSSTreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))
# m8 <- gamm(K ~ s(bSSTreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))
# m9 <- gamm(K ~ s(ALPIreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))

# m10 <- gamm(K ~ s(PDOreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# m11 <- gamm(K ~ s(nSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# m12 <- gamm(K ~ s(bSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# m13 <- gamm(K ~ s(ALPIreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))

# model.sel(m1$lme, m2$lme, m3$lme, m4$lme,m5$lme,m6$lme,m7$lme,m8$lme,
#           m9$lme,m10$lme,m11$lme,m12$lme,m13$lme, rank = AIC)

## Rerun models without correlation structure and w/ random effects
ctrl <- lmeControl(maxIter=500) #up max number of iterations to help difficult models converge

m1b <- gamm(K ~ s(pink_SE, k=3), random=list(retYr=~1), data=Nass)
m2b <- gamm(K ~ s(PDOreturn, k=3), random=list(retYr=~1), data=Nass)
m3b <- gamm(K ~ s(nSSTreturn, k=3), random=list(retYr=~1), data=Nass)
m4b <- gamm(K ~ s(bSSTreturn, k=3), random=list(retYr=~1), data=Nass)
m5b <- gamm(K ~ s(ALPIreturn, k=3), random=list(retYr=~1), data=Nass)

m6b <- gamm(K ~ s(PDOreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)
m7b <- gamm(K ~ s(nSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)
m8b <- gamm(K ~ s(bSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)
m9b <- gamm(K ~ s(ALPIreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)

m10b <- gamm(K ~ s(PDOreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass)
m11b <- gamm(K ~ s(nSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass)
m12b <- gamm(K ~ s(bSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass, control=ctrl)
m13b <- gamm(K ~ s(ALPIreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass)

model.sel(m1b$lme, m2b$lme, m3b$lme, m4b$lme, m5b$lme, m6b$lme, m7b$lme, m8b$lme, m9b$lme, m10b$lme, m11b$lme, m12b$lme, m13b$lme,
  rank = AIC)


##Top Model --> 
summary(m6b$gam)

# Effects of bSSTreturn*AGE ("partial effect" bc it's 1 predictor)
par(mfrow = c(3,2))

plot(m6b$gam, select=1, shade=TRUE, scale=0, ylim=c(-0.0005,0.0005), main= "bSST (Return Year)",
     ylab= "Effect on K for Age 1.2", xlab="bSST at return year")
abline(h=0)

plot(m6b$gam, select=2, shade=TRUE, scale=0, ylim=c(-0.0005,0.0005), main= "bSST (Return Year)",
     ylab= "Effect on K for Age 1.3", xlab="bSST at return year")
abline(h=0)

plot(m6b$gam, select=3, shade=TRUE, scale=0, ylim=c(-0.0005,0.0005), main= "bSST (Return Year)",
     ylab= "Effect on K for Age 2.2", xlab="bSST at return year")
abline(h=0)

plot(m6b$gam, select=4, shade=TRUE, scale=0, ylim=c(-0.0005,0.0005), main= "bSST (Return Year)",
     ylab= "Effect on K for Age 2.3", xlab="bSST at return year")
abline(h=0)

# Partial effects of Pink Abundance
plot(m6b$gam, select=5, shade=TRUE, scale=0, ylim=c(-0.0005,0.0005), main= "SE Alaskan Pink Abundance",
     ylab="Effect on K", xlab="Pink Abundance")
abline(h=0)

title("Nass River Condition (K) Top Model Partial Effects", outer=TRUE)


### WEIGHT_g Models ------------------------------- 
### Models labeled Wm1,Wm2...Wmn correspond to Weight in g as Response Variable


# Wm1 <- gamm(WEIGHT_g ~ s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# Wm2 <- gamm(WEIGHT_g ~ s(PDOreturn, k=3), data = Nass, fx=TRUE,correlation = corAR1(form = ~ 1|retYr))
# Wm3 <- gamm(WEIGHT_g ~ s(nSSTreturn, k=3), data = Nass,fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# Wm4 <- gamm(WEIGHT_g ~ s(bSSTreturn, k=3), data = Nass, fx=TRUE,correlation = corAR1(form = ~ 1|retYr))
# Wm5 <- gamm(WEIGHT_g ~ s(ALPIreturn, k=3), data = Nass,fx=TRUE, correlation = corAR1(form = ~ 1|retYr))

# Wm6 <- gamm(WEIGHT_g ~ s(PDOreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))
# Wm7 <- gamm(WEIGHT_g ~ s(nSSTreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))
# Wm8 <- gamm(WEIGHT_g ~ s(bSSTreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))
# Wm9 <- gamm(WEIGHT_g ~ s(ALPIreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))

# Wm10 <- gamm(WEIGHT_g ~ s(PDOreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# Wm11 <- gamm(WEIGHT_g ~ s(nSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# Wm12 <- gamm(WEIGHT_g ~ s(bSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# Wm13 <- gamm(WEIGHT_g ~ s(ALPIreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))

Wm1b <- gamm(WEIGHT_g ~ s(pink_SE, k=3), random=list(retYr=~1), data=Nass)
Wm2b <- gamm(WEIGHT_g ~ s(PDOreturn, k=3), random=list(retYr=~1), data=Nass)
Wm3b <- gamm(WEIGHT_g ~ s(nSSTreturn, k=3), random=list(retYr=~1), data=Nass)
Wm4b <- gamm(WEIGHT_g ~ s(bSSTreturn, k=3), random=list(retYr=~1), data=Nass)
Wm5b <- gamm(WEIGHT_g ~ s(ALPIreturn, k=3), random=list(retYr=~1), data=Nass)

Wm6b <- gamm(WEIGHT_g ~ s(PDOreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)
Wm7b <- gamm(WEIGHT_g ~ s(nSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)
Wm8b <- gamm(WEIGHT_g ~ s(bSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)
Wm9b <- gamm(WEIGHT_g ~ s(ALPIreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)

Wm10b <- gamm(WEIGHT_g ~ s(PDOreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass)
Wm11b <- gamm(WEIGHT_g ~ s(nSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass)
Wm12b <- gamm(WEIGHT_g ~ s(bSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass)
Wm13b <- gamm(WEIGHT_g ~ s(ALPIreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass, control=ctrl)

model.sel(Wm1b$lme, Wm2b$lme, Wm3b$lme, Wm4b$lme, Wm5b$lme, Wm6b$lme, Wm7b$lme, Wm8b$lme,
          Wm9b$lme, Wm10b$lme, Wm11b$lme, Wm12b$lme, Wm13b$lme, rank = AIC)
          
## Top Model --> Wm10
summary(Wm10$gam)

# Effects of bSST ("partial effect" bc it's 1 predictor)
par(mfrow = c(3,2))

plot(Wm10$gam, select=1, shade=TRUE, scale=0, ylim=c(-500,1000), main= "PDO (Return Yr.)",
     ylab= "Effect on Weight (g) for Age 1.2", xlab="PDO (Return)")
abline(h=0)

plot(Wm10$gam, select=2, shade=TRUE, scale=0, ylim=c(-500,1000), main= "PDO (Return Yr.)",
     ylab= "Effect on Weight (g) for Age 1.3", xlab="PDO (Return)")
abline(h=0)

plot(Wm10$gam, select=3, shade=TRUE, scale=0, ylim=c(-500,1000), main= "PDO (Return Yr.)",
     ylab= "Effect on Weight (g) for Age 2.2", xlab="PDO (Return)")
abline(h=0)

plot(Wm10$gam, select=4, shade=TRUE, scale=0, ylim=c(-500,1000), main= "PDO (Return Yr.)",
     ylab= "Effect on Weight (g) for Age 2.3", xlab="PDO (Return)")
abline(h=0)

# Partial effects of Pink Abundance
plot(Wm10$gam, select=5, shade=TRUE, scale=0, ylim=c(-500,1000), main= "SE Alaskan Pink Abundance",
     ylab="Effect on Weight (g)", xlab="Pink Abundance (Return)")
abline(h=0)

title("Nass River Weight (g) Top Model Partial Effects", outer=TRUE)


### Length_mm Models ------------------------------- 
### Models labeled LTm1,LTm2...LTmn correspond to Length in mm as Response Variable


# LTm1 <- gamm(LENGTH_mm ~ s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# LTm2 <- gamm(LENGTH_mm ~ s(PDOreturn, k=3), data = Nass, fx=TRUE,correlation = corAR1(form = ~ 1|retYr))
# LTm3 <- gamm(LENGTH_mm ~ s(nSSTreturn, k=3), data = Nass,fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# LTm4 <- gamm(LENGTH_mm ~ s(bSSTreturn, k=3), data = Nass, fx=TRUE,correlation = corAR1(form = ~ 1|retYr))
# LTm5 <- gamm(LENGTH_mm ~ s(ALPIreturn, k=3), data = Nass,fx=TRUE, correlation = corAR1(form = ~ 1|retYr))

# LTm6 <- gamm(LENGTH_mm ~ s(PDOreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))
# LTm7 <- gamm(LENGTH_mm ~ s(nSSTreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))
# LTm8 <- gamm(LENGTH_mm ~ s(bSSTreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))
# LTm9 <- gamm(LENGTH_mm ~ s(ALPIreturn, by = AGE, k=3), fx=TRUE,data = Nass, correlation = corAR1(form = ~ 1|retYr))

# LTm10 <- gamm(LENGTH_mm ~ s(PDOreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# LTm11 <- gamm(LENGTH_mm ~ s(nSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# LTm12 <- gamm(LENGTH_mm ~ s(bSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))
# LTm13 <- gamm(LENGTH_mm ~ s(ALPIreturn, by = AGE, k=3) + s(pink_SE, k=3), data = Nass, fx=TRUE, correlation = corAR1(form = ~ 1|retYr))

LTm1b <- gamm(LENGTH_mm ~ s(pink_SE, k=3), random=list(retYr=~1), data=Nass)
LTm2b <- gamm(LENGTH_mm ~ s(PDOreturn, k=3), random=list(retYr=~1), data=Nass)
LTm3b <- gamm(LENGTH_mm ~ s(nSSTreturn, k=3), random=list(retYr=~1), data=Nass)
LTm4b <- gamm(LENGTH_mm ~ s(bSSTreturn, k=3), random=list(retYr=~1), data=Nass)
LTm5b <- gamm(LENGTH_mm ~ s(ALPIreturn, k=3), random=list(retYr=~1), data=Nass)

LTm6b <- gamm(LENGTH_mm ~ s(PDOreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)
LTm7b <- gamm(LENGTH_mm ~ s(nSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)
LTm8b <- gamm(LENGTH_mm ~ s(bSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)
LTm9b <- gamm(LENGTH_mm ~ s(ALPIreturn, by = AGE, k=3), random=list(retYr=~1), data=Nass)

LTm10b <- gamm(LENGTH_mm ~ s(PDOreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass)
LTm11b <- gamm(LENGTH_mm ~ s(nSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass)
LTm12b <- gamm(LENGTH_mm ~ s(bSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass)
LTm13b <- gamm(LENGTH_mm ~ s(ALPIreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Nass, control=ctrl)

model.sel(LTm1b$lme, LTm2b$lme, LTm3b$lme, LTm4b$lme, LTm5b$lme, LTm6b$lme, LTm7b$lme, LTm8b$lme,
          LTm9b$lme, LTm10b$lme, LTm11b$lme, LTm12b$lme, LTm13b$lme, rank = AIC)

## Top Model --> LTm13
summary(LTm10b$gam)
summary(LTm10c)
plot_acf(LTm10c)
check_resid(LTm10c)

LTm10c <- gam(LENGTH_mm ~ s(PDOreturn, k=3) + s(pink_SE, k=3) + s(retYr, bs="re"), data=Nass)


# Effects of bSST ("partial effect" bc it's 1 predictor)
par(mfrow = c(3,2))

plot(LTm10b$gam, select=1, shade=TRUE, scale=0, ylim=c(-50,75), main= "Nass ALPI (Return Yr.)",
     ylab= "Effect on Length (mm) for Age 1.2", xlab="ALPI (Return Yr.)")
abline(h=0)

plot(LTm10b$gam, select=2, shade=TRUE, scale=0, ylim=c(-50,75), main= "ALPI (Return Yr.)",
     ylab= "Effect on Length (mm) for Age 1.3", xlab="ALPI (Return Yr.)")
abline(h=0)

plot(LTm10b$gam, select=3, shade=TRUE, scale=0, ylim=c(-50,75), main= "ALPI (Return Yr.)",
     ylab= "Effect on Length (mm) for Age 2.2", xlab="ALPI (Return Yr.)")
abline(h=0)

plot(LTm10b$gam, select=4, shade=TRUE, scale=0, ylim=c(-50,75), main= "ALPI (Return Yr.)",
     ylab= "Effect on Length (mm) for Age 2.3", xlab="ALPI (Return Yr.)")
abline(h=0)

# Partial effects of Pink Abundance
plot(LTm10b$gam, select=5, shade=TRUE, scale=0, ylim=c(-50,75), main= "SE Alaskan Pink Abundance",
     ylab="Effect on Length (mm)", xlab="Pink Abundance")
abline(h=0)

# title("Rivers Inlet Length (mm) Top Model Partial Effects", outer=FALSE)


# -----------------------------------------------
# -----------------------------------------------
### --------------GAMMs RIVERS INLET----------###
# -----------------------------------------------
# -----------------------------------------------


library(mgcv); library(nlme);library(usdm);library(MuMIn)


### Fulton's K Models ------------------------------- 
### Models labeled R_m1,R_m2...R_mn correspond to Fultonk K as Response Variable


R_m1 <- gamm(K ~ s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_m2 <- gamm(K ~ s(PDOreturn, k=3), random=list(retYr=~1), data=Rivers)
R_m3 <- gamm(K ~ s(nSSTreturn, k=3), random=list(retYr=~1), data=Rivers)
R_m4 <- gamm(K ~ s(bSSTreturn, k=3), random=list(retYr=~1), data=Rivers)
R_m5 <- gamm(K ~ s(ALPIreturn, k=3), random=list(retYr=~1), data=Rivers)

R_m6 <- gamm(K ~ s(PDOreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)
R_m7 <- gamm(K ~ s(nSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)
R_m8 <- gamm(K ~ s(bSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)
R_m9 <- gamm(K ~ s(ALPIreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)

R_m10 <- gamm(K ~ s(PDOreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_m11 <- gamm(K ~ s(nSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_m12 <- gamm(K ~ s(bSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_m13 <- gamm(K ~ s(ALPIreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)

model.sel(R_m1$lme, R_m2$lme, R_m3$lme, R_m4$lme, R_m5$lme, R_m6$lme, R_m7$lme, R_m8$lme,
          R_m9$lme, R_m10$lme, R_m11$lme, R_m12$lme, R_m13$lme, rank = AIC)

##Top Model --> R_m10
summary(R_m6$gam)

# Effects of PDOreturn*AGE ("partial effect" bc it's 1 predictor)
par(mfrow = c(3,1))

plot(R_m10$gam, select=1, shade=TRUE, scale=0, ylim=c(-0.0005,0.0005), main= "PDO (Return Year)",
     ylab= "Effect on Sockeye Condition (K) for Age 1.2", xlab="PDO at return year")
abline(h=0)

plot(R_m10$gam, select=2, shade=TRUE, scale=0, ylim=c(-0.0005,0.0005), main= "PDO (Return Year)",
     ylab= "Effect on Sockeye Condition (K) for Age 1.3", xlab="PDO at return year")
abline(h=0)

# Partial effects of Pink Abundance
plot(R_m10$gam, select=3, shade=TRUE, scale=0, ylim=c(-0.0005,0.0005), main= "SE Alaskan Pink Abundance",
     ylab="Effect on Sockeye Condition (K)", xlab="Pink Abundance")
abline(h=0)

title("Rivers Condition (K) Top Model Partial Effects", outer=TRUE)


### WEIGHT_g Models ------------------------------- 
### Models labeled R_Wm1,R_Wm2...R_Wmn correspond to Weight in g as Response Variable


R_Wm1 <- gamm(WEIGHT_g ~ s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_Wm2 <- gamm(WEIGHT_g ~ s(PDOreturn, k=3), random=list(retYr=~1), data=Rivers)
R_Wm3 <- gamm(WEIGHT_g ~ s(nSSTreturn, k=3), random=list(retYr=~1), data=Rivers)
R_Wm4 <- gamm(WEIGHT_g ~ s(bSSTreturn, k=3), random=list(retYr=~1), data=Rivers)
R_Wm5 <- gamm(WEIGHT_g ~ s(ALPIreturn, k=3), random=list(retYr=~1), data=Rivers)

R_Wm6 <- gamm(WEIGHT_g ~ s(PDOreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)
R_Wm7 <- gamm(WEIGHT_g ~ s(nSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)
R_Wm8 <- gamm(WEIGHT_g ~ s(bSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)
R_Wm9 <- gamm(WEIGHT_g ~ s(ALPIreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)

R_Wm10 <- gamm(WEIGHT_g ~ s(PDOreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_Wm11 <- gamm(WEIGHT_g ~ s(nSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_Wm12 <- gamm(WEIGHT_g ~ s(bSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_Wm13 <- gamm(WEIGHT_g ~ s(ALPIreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)

model.sel(R_Wm1$lme, R_Wm2$lme, R_Wm3$lme, R_Wm4$lme, R_Wm5$lme, R_Wm6$lme, R_Wm7$lme, R_Wm8$lme,
          R_Wm9$lme, R_Wm10$lme, R_Wm11$lme, R_Wm12$lme, R_Wm13$lme, rank = AIC)
          
## Top Model --> R_Wm12
summary(R_Wm6$gam)

# Effects of bSST ("partial effect" bc it's 1 predictor)
par(mfrow = c(3,1))

plot(R_Wm12$gam, select=1, shade=TRUE, scale=0, ylim=c(-500,1000), main= "Broad-scale SST",
     ylab= "Effect on Sockeye Weight (g) for Age 1.2", xlab="Broad SST")
abline(h=0)

plot(R_Wm12$gam, select=2, shade=TRUE, scale=0, ylim=c(-500,1000), main= "Broad-scale SST",
     ylab= "Effect on Sockeye Weight (g) for Age 1.3", xlab="Broad SST")
abline(h=0)

# Partial effects of Pink Abundance
plot(R_Wm12$gam, select=3, shade=TRUE, scale=0, ylim=c(-500,1000), main= "SE Alaskan Pink Abundance",
     ylab="Effect on Sockeye Weight (g)", xlab="Pink Abundance")
abline(h=0)

title("Rivers Inlet Weight (g) Top Model Partial Effects", outer=TRUE)


### Length_mm Models ------------------------------- 
### Models labeled R_LTm1,R_LTm2...R_LTmn correspond to Length in mm as Response Variable


R_LTm1 <- gamm(LENGTH_mm ~ s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_LTm2 <- gamm(LENGTH_mm ~ s(PDOreturn, k=3), random=list(retYr=~1), data=Rivers)
R_LTm3 <- gamm(LENGTH_mm ~ s(nSSTreturn, k=3), random=list(retYr=~1), data=Rivers)
R_LTm4 <- gamm(LENGTH_mm ~ s(bSSTreturn, k=3), random=list(retYr=~1), data=Rivers)
R_LTm5 <- gamm(LENGTH_mm ~ s(ALPIreturn, k=3), random=list(retYr=~1), data=Rivers)

R_LTm6 <- gamm(LENGTH_mm ~ s(PDOreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)
R_LTm7 <- gamm(LENGTH_mm ~ s(nSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)
R_LTm8 <- gamm(LENGTH_mm ~ s(bSSTreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)
R_LTm9 <- gamm(LENGTH_mm ~ s(ALPIreturn, by = AGE, k=3), random=list(retYr=~1), data=Rivers)

R_LTm10 <- gamm(LENGTH_mm ~ s(PDOreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_LTm11 <- gamm(LENGTH_mm ~ s(nSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_LTm12 <- gamm(LENGTH_mm ~ s(bSSTreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)
R_LTm13 <- gamm(LENGTH_mm ~ s(ALPIreturn, by = AGE, k=3) + s(pink_SE, k=3), random=list(retYr=~1), data=Rivers)

model.sel(R_LTm1$lme, R_LTm2$lme, R_LTm3$lme, R_LTm4$lme, R_LTm5$lme, R_LTm6$lme, R_LTm7$lme, R_LTm8$lme,
          R_LTm9$lme, R_LTm10$lme, R_LTm11$lme, R_LTm12$lme, R_LTm13$lme, rank = AIC)

## Top Model --> R_LTm12
summary(R_LTm10$gam)
summary(R_LTm2$gam)
check_resid(R_LTm2)

# Effects of bSST ("partial effect" bc it's 1 predictor)
par(mfrow = c(3,3))

plot(R_LTm10$gam, select=1, shade=TRUE, scale=0, ylim=c(-50,75), main= "Broad-scale SST",
     ylab= "Effect on Sockeye Length (mm) for Age 1.2", xlab="Broad SST")
abline(h=0)

plot(R_LTm10$gam, select=1, shade=TRUE, scale=0, ylim=c(-50,75), main= "Broad-scale SST",
     ylab= "Effect on Sockeye Length (mm) for Age 1.3", xlab="Broad SST")
abline(h=0)

# Partial effects of Pink Abundance
plot(R_LTm10$gam, select=2, shade=TRUE, scale=0, ylim=c(-50,75), main= "SE Alaskan Pink Abundance",
     ylab="Effect on Sockeye Length (mm)", xlab="Pink Abundance")
abline(h=0)


title("Rivers Inlet Length (mm) Top Model Partial Effects", outer=FALSE)












# -----------------------------------------------
# -----------------------------------------------
### ------------OLDER FIGURES & WORK----------###
# -----------------------------------------------
# -----------------------------------------------










# REGRESSION GRAPHS by average oceanographic conditions only ---------------

### Condition (L-W residuals) -------------------------- ##

# NASS 1.2's---------------------------------------------------------------------------
#AVERAGE
par(mfrow = c(2, 2))
plot(Nass1.2$LWcondition ~ Nass1.2$avgPDO,col= "grey50")
mod1 <- lm(Nass1.2$LWcondition ~ Nass1.2$avgPDO); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LWcondition ~ Nass1.2$avgALPI,col= "grey50")
mod2 <- lm(Nass1.2$LWcondition ~ Nass1.2$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(4,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LWcondition ~ Nass1.2$avgSSTn,col= "grey50")
mod3 <- lm(Nass1.2$LWcondition ~ Nass1.2$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LWcondition ~ Nass1.2$avgSSTb,col= "grey50")
mod4 <- lm(Nass1.2$LWcondition ~ Nass1.2$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY

par(mfrow = c(2, 2))
plot(Nass1.2$LWcondition ~ Nass1.2$PDOentry,col= "grey50")
mod1 <- lm(Nass1.2$LWcondition ~ Nass1.2$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LWcondition ~ Nass1.2$ALPIentry,col= "grey50")
mod2 <- lm(Nass1.2$LWcondition ~ Nass1.2$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(6,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LWcondition ~ Nass1.2$nSSTentry,col= "grey50")
mod3 <- lm(Nass1.2$LWcondition ~ Nass1.2$nSSTentry); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LWcondition ~ Nass1.2$bSSTentry,col= "grey50")
mod4 <- lm(Nass1.2$LWcondition ~ Nass1.2$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# NASS 1.3's----------------------------------------------------------------

# AVERAGE
par(mfrow = c(2, 2))
plot(Nass1.3$LWcondition ~Nass1.3$avgPDO,col= "grey50")
mod1 <- lm(Nass1.3$LWcondition ~Nass1.3$avgPDO); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", p))

plot(Nass1.3$LWcondition ~Nass1.3$avgALPI,col= "grey50")
mod2 <- lm(Nass1.3$LWcondition ~Nass1.3$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(2,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$LWcondition ~Nass1.3$avgSSTn,col= "grey50")
mod3 <- lm(Nass1.3$LWcondition ~Nass1.3$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$LWcondition ~Nass1.3$avgSSTb,col= "grey50")
mod4 <- lm(Nass1.3$LWcondition ~Nass1.3$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(6.6,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY -------

par(mfrow = c(1, 4), mar = c(4, 0, 0, 0), oma = c(1,5,2,1))
plot(Nass1.3$K ~Nass1.3$PDOentry,col= "grey50", pch=19, xlab="PDO")
mod1 <- lm(Nass1.3$K ~Nass1.3$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,0.0195, paste("R =", round(r, 3), "\nP =", round(p, 3)))

plot(Nass1.3$K ~Nass1.3$ALPIentry,col= "grey50", ylab="", yaxt="n", pch=19, xlab="ALPI")
mod2 <- lm(Nass1.3$K ~Nass1.3$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(4,0.0195, paste("R =", round(r, 3),"\nP =", round(p, 3)))

plot(Nass1.3$K ~Nass1.3$nSSTentry,col= "grey50",ylab="", yaxt="n", pch=19, xlab="Regional SST")
mod3 <- lm(Nass1.3$K ~Nass1.3$nSSTentry); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.5,0.0195, paste("R =", round(r, 3), "\nP =", round(p, 3)))

plot(Nass1.3$K ~Nass1.3$bSSTentry,col= "grey50",ylab="", yaxt="n", pch=19 , xlab="Broad SST")
mod4 <- lm(Nass1.3$K ~Nass1.3$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(6.75,0.0195, paste("R =", round(r, 3), "\nP =", round(p, 3)))

mtext("Sockeye Age Class 1.3 - Condition (K)", side = 2, outer = TRUE, cex = 0.9, line = 2.5, col = "grey20")
mtext("Nass River 1911-1946", side = 3, outer = TRUE, cex = 0.9, line = 1, col = "grey20")


# NASS 2.2's --------------------------------------------------------

par(mfrow = c(2, 2))
plot(Nass2.2$LWcondition ~Nass2.2$avgPDO,col= "grey50")
mod1 <- lm(Nass2.2$LWcondition ~Nass2.2$avgPDO); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$LWcondition ~Nass2.2$avgALPI,col= "grey50")
mod2 <- lm(Nass2.2$LWcondition ~Nass2.2$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(3.5,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$LWcondition ~Nass2.2$avgSSTn,col= "grey50")
mod3 <- lm(Nass2.2$LWcondition ~Nass2.2$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$LWcondition ~Nass2.2$avgSSTb,col= "grey50")
mod4 <- lm(Nass2.2$LWcondition ~Nass2.2$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY-------

par(mfrow = c(1, 4), mar = c(4, 0, 0, 0), oma = c(1,5,2,1))
plot(Nass2.2$K ~Nass2.2$PDOentry,col= "grey50", pch=19, xlab="PDO")
mod1 <- lm(Nass2.2$K ~Nass2.2$PDOentry); abline(mod1, col="blue",lty=1,lwd=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.8,0.0192, paste("R =", round(r, 3), "\nP =", round(p, 3)))

plot(Nass2.2$K ~Nass2.2$ALPIentry,col= "grey50", ylab="", yaxt="n", pch=19, xlab="ALPI")
mod2 <- lm(Nass2.2$K ~Nass2.2$ALPIentry); abline(mod2, col="blue",lty=1,lwd=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(5,0.0192, paste("R =", round(r, 3),"\nP =", round(p, 3)))

plot(Nass2.2$K ~Nass2.2$nSSTentry,col= "grey50",ylab="", yaxt="n", pch=19, xlab="Regional SST")
mod3 <- lm(Nass2.2$K ~Nass2.2$nSSTentry); abline(mod3, col="blue",lty=1,lwd=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,0.0192, paste("R =", round(r, 3), "\nP =", round(p, 3)))

plot(Nass2.2$K ~Nass2.2$bSSTentry,col= "grey50",ylab="", yaxt="n", pch=19 , xlab="Broad SST")
mod4 <- lm(Nass2.2$K ~Nass2.2$bSSTentry); abline(mod4, col="blue",lty=1,lwd=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,0.0192, paste("R =", round(r, 3), "\nP =", round(p, 3)))

mtext("Sockeye Age Class 2.2 \n Condition (K)", side = 2, outer = TRUE, cex = 0.9, line = 2.6, col = "grey20")
mtext("Nass River 1914-1946", side = 3, outer = TRUE, cex = 0.9, line = 1, col = "grey20")

# NASS 2.3's----------------------------------------------------------

# AVERAGE
par(mfrow = c(2, 2))
plot(Nass2.3$LWcondition ~Nass2.3$avgPDO,col= "grey50")
mod1 <- lm(Nass2.3$LWcondition ~Nass2.3$avgPDO); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1,900, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LWcondition ~Nass2.3$avgALPI,col= "grey50")
mod2 <- lm(Nass2.3$LWcondition ~Nass2.3$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(2,900, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LWcondition ~Nass2.3$avgSSTn,col= "grey50")
mod3 <- lm(Nass2.3$LWcondition ~Nass2.3$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9,900, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LWcondition ~Nass2.3$avgSSTb,col= "grey50")
mod4 <- lm(Nass2.3$LWcondition ~Nass2.3$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(6.6,900, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY-------

par(mfrow = c(2, 2))
plot(Nass2.3$LWcondition ~Nass2.3$PDOentry,col= "grey50")
mod1 <- lm(Nass2.3$LWcondition ~Nass2.3$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,900, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LWcondition ~Nass2.3$ALPIentry,col= "grey50")
mod2 <- lm(Nass2.3$LWcondition ~Nass2.3$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(6,900, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LWcondition ~Nass2.3$nSSTentry,col= "grey50")
mod3 <- lm(Nass2.3$LWcondition ~Nass2.3$nSSTentry); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,900, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LWcondition ~Nass2.3$bSSTentry,col= "grey50")
mod4 <- lm(Nass2.3$LWcondition ~Nass2.3$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,900, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# RIVERS INLET---------------------------------------------------------

# 1.2's
# AVERAGE
par(mfrow = c(2, 2))
plot(Rivers1.2$LWcondition ~Rivers1.2$avgPDO,col= "grey50")
mod1 <- lm(Rivers1.2$LWcondition ~Rivers1.2$avgPDO); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$LWcondition ~Rivers1.2$avgALPI,col= "grey50")
mod2 <- lm(Rivers1.2$LWcondition ~Rivers1.2$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(3,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$LWcondition ~Rivers1.2$avgSSTn,col= "grey50")
mod3 <- lm(Rivers1.2$LWcondition ~Rivers1.2$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$LWcondition ~Rivers1.2$avgSSTb,col= "grey50")
mod4 <- lm(Rivers1.2$LWcondition ~Rivers1.2$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# Entry-------

par(mfrow = c(1, 4), mar = c(4, 0, 0, 0), oma = c(1,5,2,1))
plot(Rivers1.2$K ~Rivers1.2$PDOentry,col= "grey50", pch=19, xlab="PDO")
mod1 <- lm(Rivers1.2$K ~Rivers1.2$PDOentry); abline(mod1, col="red",lty=1,lwd=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.8,0.023, paste("R =", round(r, 3), "\nP =", round(p, 3)))

plot(Rivers1.2$K ~Rivers1.2$ALPIentry,col= "grey50", ylab="", yaxt="n", pch=19, xlab="ALPI")
mod2 <- lm(Rivers1.2$K ~Rivers1.2$ALPIentry); abline(mod2, col="red",lty=1, lwd=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(5,0.023, paste("R =", round(r, 3),"\nP =", round(p, 3)))

plot(Rivers1.2$K ~Rivers1.2$nSSTentry,col= "grey50",ylab="", yaxt="n", pch=19, xlab="Regional SST")
mod3 <- lm(Rivers1.2$K ~Rivers1.2$nSSTentry); abline(mod3, col="red",lty=1,lwd=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,0.023, paste("R =", round(r, 3), "\nP =", round(p, 3)))

plot(Rivers1.2$K ~Rivers1.2$bSSTentry,col= "grey50",ylab="", yaxt="n", pch=19 , xlab="Broad SST")
mod4 <- lm(Rivers1.2$K ~Rivers1.2$bSSTentry); abline(mod4, col="red",lty=1, lwd=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,0.023, paste("R =", round(r, 3), "\nP =", round(p, 3)))

mtext("Sockeye Age Class 1.2 \n Condition (K)", side = 2, outer = TRUE, cex = 0.9, line = 2.6, col = "grey20")
mtext("Rivers Inlet 1915-1946", side = 3, outer = TRUE, cex = 0.9, line = 1, col = "grey20")

# RIVERS 1.3---------------------------------------------------------

# AVERAGE
par(mfrow = c(2, 2))
plot(Rivers1.3$LWcondition ~Rivers1.3$avgPDO,col= "grey50"); mod1 <- lm(Rivers1.3$LWcondition ~Rivers1.3$avgPDO)
abline(mod1, col="blue",lty=2); tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$LWcondition ~Rivers1.3$avgALPI,col= "grey50")
mod2 <- lm(Rivers1.3$LWcondition ~Rivers1.3$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(2,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$LWcondition ~Rivers1.3$avgSSTn,col= "grey50")
mod3 <- lm(Rivers1.3$LWcondition ~Rivers1.3$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$LWcondition ~Rivers1.3$avgSSTb,col= "grey50")
mod4 <- lm(Rivers1.3$LWcondition ~Rivers1.3$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(6.6,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY -------

par(mfrow = c(1, 4), mar = c(4, 0, 0, 0), oma = c(1,5,2,1))
plot(Rivers1.3$K ~Rivers1.3$PDOentry,col= "grey50", pch=19, xlab="PDO")
mod1 <- lm(Rivers1.3$K ~Rivers1.3$PDOentry); abline(mod1, col="red",lty=1)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.8,0.025, paste("R =", round(r, 3), "\nP =", round(p, 3)))

plot(Rivers1.3$K ~Rivers1.3$ALPIentry,col= "grey50", ylab="", yaxt="n", pch=19, xlab="ALPI")
mod2 <- lm(Rivers1.3$K ~Rivers1.3$ALPIentry); abline(mod2, col="red",lty=1)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(5,0.025, paste("R =", round(r, 3),"\nP =", round(p, 3)))

plot(Rivers1.3$K ~Rivers1.3$nSSTentry,col= "grey50",ylab="", yaxt="n", pch=19, xlab="Regional SST")
mod3 <- lm(Rivers1.3$K ~Rivers1.3$nSSTentry); abline(mod3, col="red",lty=1)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,0.025, paste("R =", round(r, 3), "\nP =", round(p, 3)))

plot(Rivers1.3$K ~Rivers1.3$bSSTentry,col= "grey50",ylab="", yaxt="n", pch=19 , xlab="Broad SST")
mod4 <- lm(Rivers1.3$K ~Rivers1.3$bSSTentry); abline(mod4, col="red",lty=1)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,0.025, paste("R =", round(r, 3), "\nP =", round(p, 3)))

mtext("Sockeye Age Class 1.3 \n Condition (K)", side = 2, outer = TRUE, cex = 0.9, line = 2.6, col = "grey20")
mtext("Rivers Inlet 1915-1946", side = 3, outer = TRUE, cex = 0.9, line = 1, col = "grey20")

# ----------------------------------------------------------------------
# REGRESSION GRAPHS -------------WEIGHT --------------------------------


# Nass. Riv 1.2's 
# AVERAGE
par(mfrow = c(2, 2))
plot(Nass1.2$WEIGHT_g ~ Nass1.2$avgPDO,col= "grey50"); mod1 <- lm(Nass1.2$WEIGHT_g ~ Nass1.2$avgPDO)
abline(mod1, col="blue",lty=2); tab <- tidy(mod1)
b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$WEIGHT_g ~ Nass1.2$avgALPI,col= "grey50")
mod2 <- lm(Nass1.2$WEIGHT_g ~ Nass1.2$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(4,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$WEIGHT_g ~ Nass1.2$avgSSTn,col= "grey50")
mod3 <- lm(Nass1.2$WEIGHT_g ~ Nass1.2$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$WEIGHT_g ~ Nass1.2$avgSSTb,col= "grey50")
mod4 <- lm(Nass1.2$WEIGHT_g ~ Nass1.2$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY

par(mfrow = c(2, 2))
plot(Nass1.2$WEIGHT_g ~ Nass1.2$PDOentry,col= "grey50")
mod1 <- lm(Nass1.2$WEIGHT_g ~ Nass1.2$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$WEIGHT_g ~ Nass1.2$ALPIentry,col= "grey50")
mod2 <- lm(Nass1.2$WEIGHT_g ~ Nass1.2$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(6,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$WEIGHT_g ~ Nass1.2$nSSTentry,col= "grey50")
mod3 <- lm(Nass1.2$WEIGHT_g ~ Nass1.2$nSSTentry); abline(mod3, col="blue",lty=2); 
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$WEIGHT_g ~ Nass1.2$bSSTentry,col= "grey50")
mod4 <- lm(Nass1.2$WEIGHT_g ~ Nass1.2$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# Nass. Riv 2.2's  ---------------------------------------------------------------
# AVERAGE
par(mfrow = c(2, 2))
plot(Nass2.2$WEIGHT_g ~ Nass2.2$avgPDO,col= "grey50"); mod1 <- lm(Nass2.2$WEIGHT_g ~ Nass2.2$avgPDO)
abline(mod1, col="blue",lty=2); tab <- tidy(mod1)
b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$WEIGHT_g ~ Nass2.2$avgALPI,col= "grey50")
mod2 <- lm(Nass2.2$WEIGHT_g ~ Nass2.2$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(4,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$WEIGHT_g ~ Nass2.2$avgSSTn,col= "grey50")
mod3 <- lm(Nass2.2$WEIGHT_g ~ Nass2.2$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$WEIGHT_g ~ Nass2.2$avgSSTb,col= "grey50")
mod4 <- lm(Nass2.2$WEIGHT_g ~ Nass2.2$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY

par(mfrow = c(2, 2))
plot(Nass2.2$WEIGHT_g ~ Nass2.2$PDOentry,col= "grey50")
mod1 <- lm(Nass2.2$WEIGHT_g ~ Nass2.2$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$WEIGHT_g ~ Nass2.2$ALPIentry,col= "grey50")
mod2 <- lm(Nass2.2$WEIGHT_g ~ Nass2.2$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(6,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$WEIGHT_g ~ Nass2.2$nSSTentry,col= "grey50")
mod3 <- lm(Nass2.2$WEIGHT_g ~ Nass2.2$nSSTentry); abline(mod3, col="blue",lty=2); 
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$WEIGHT_g ~ Nass2.2$bSSTentry,col= "grey50")
mod4 <- lm(Nass2.2$WEIGHT_g ~ Nass2.2$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# Nass. Riv 2.3's  ---------------------------------------------------------------
# AVERAGE
par(mfrow = c(2, 2))
plot(Nass2.3$WEIGHT_g ~ Nass2.3$avgPDO,col= "grey50"); mod1 <- lm(Nass2.3$WEIGHT_g ~ Nass2.3$avgPDO)
abline(mod1, col="blue",lty=2); tab <- tidy(mod1)
b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$WEIGHT_g ~ Nass2.3$avgALPI,col= "grey50")
mod2 <- lm(Nass2.3$WEIGHT_g ~ Nass2.3$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(2.5,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$WEIGHT_g ~ Nass2.3$avgSSTn,col= "grey50")
mod3 <- lm(Nass2.3$WEIGHT_g ~ Nass2.3$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$WEIGHT_g ~ Nass2.3$avgSSTb,col= "grey50")
mod4 <- lm(Nass2.3$WEIGHT_g ~ Nass2.3$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(6.8,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY

par(mfrow = c(2, 2))
plot(Nass2.3$WEIGHT_g ~ Nass2.3$PDOentry,col= "grey50")
mod1 <- lm(Nass2.3$WEIGHT_g ~ Nass2.3$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$WEIGHT_g ~ Nass2.3$ALPIentry,col= "grey50")
mod2 <- lm(Nass2.3$WEIGHT_g ~ Nass2.3$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(6,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$WEIGHT_g ~ Nass2.3$nSSTentry,col= "grey50")
mod3 <- lm(Nass2.3$WEIGHT_g ~ Nass2.3$nSSTentry); abline(mod3, col="blue",lty=2); 
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$WEIGHT_g ~ Nass2.3$bSSTentry,col= "grey50")
mod4 <- lm(Nass2.3$WEIGHT_g ~ Nass2.3$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# Nass. Riv 1.3's  ---------------------------------------------------------------
# AVERAGE
par(mfrow = c(2, 2))
plot(Nass1.3$WEIGHT_g ~ Nass1.3$avgPDO,col= "grey50"); mod1 <- lm(Nass1.3$WEIGHT_g ~ Nass1.3$avgPDO)
abline(mod1, col="blue",lty=2); tab <- tidy(mod1)
b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$WEIGHT_g ~ Nass1.3$avgALPI,col= "grey50")
mod2 <- lm(Nass1.3$WEIGHT_g ~ Nass1.3$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(2.5,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$WEIGHT_g ~ Nass1.3$avgSSTn,col= "grey50")
mod3 <- lm(Nass1.3$WEIGHT_g ~ Nass1.3$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$WEIGHT_g ~ Nass1.3$avgSSTb,col= "grey50")
mod4 <- lm(Nass1.3$WEIGHT_g ~ Nass1.3$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(6.8,1800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY

par(mfrow = c(2, 2))
plot(Nass1.3$WEIGHT_g ~ Nass1.3$PDOentry,col= "grey50")
mod1 <- lm(Nass1.3$WEIGHT_g ~ Nass1.3$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$WEIGHT_g ~ Nass1.3$ALPIentry,col= "grey50")
mod2 <- lm(Nass1.3$WEIGHT_g ~ Nass1.3$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(6,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$WEIGHT_g ~ Nass1.3$nSSTentry,col= "grey50")
mod3 <- lm(Nass1.3$WEIGHT_g ~ Nass1.3$nSSTentry); abline(mod3, col="blue",lty=2); 
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$WEIGHT_g ~ Nass1.3$bSSTentry,col= "grey50")
mod4 <- lm(Nass1.3$WEIGHT_g ~ Nass1.3$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# RIVERS INLET---------------------------------------------------------

# 1.2's
# AVERAGE
par(mfrow = c(2, 2))
plot(Rivers1.2$WEIGHT_g ~Rivers1.2$avgPDO,col= "grey50")
mod1 <- lm(Rivers1.2$WEIGHT_g ~Rivers1.2$avgPDO); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$WEIGHT_g ~Rivers1.2$avgALPI,col= "grey50")
mod2 <- lm(Rivers1.2$WEIGHT_g ~Rivers1.2$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(3,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$WEIGHT_g ~Rivers1.2$avgSSTn,col= "grey50")
mod3 <- lm(Rivers1.2$WEIGHT_g ~Rivers1.2$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$WEIGHT_g ~Rivers1.2$avgSSTb,col= "grey50")
mod4 <- lm(Rivers1.2$WEIGHT_g ~Rivers1.2$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))
#-------

par(mfrow = c(2, 2))
plot(Rivers1.2$WEIGHT_g ~Rivers1.2$PDOentry,col= "grey50")
mod1 <- lm(Rivers1.2$WEIGHT_g ~Rivers1.2$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,4500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$WEIGHT_g ~Rivers1.2$ALPIentry,col= "grey50")
mod2 <- lm(Rivers1.2$WEIGHT_g ~Rivers1.2$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(5,4500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$WEIGHT_g ~Rivers1.2$nSSTentry,col= "grey50")
mod3 <- lm(Rivers1.2$WEIGHT_g ~Rivers1.2$nSSTentry); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(10,4500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$WEIGHT_g ~Rivers1.2$bSSTentry,col= "grey50")
mod4 <- lm(Rivers1.2$WEIGHT_g ~Rivers1.2$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7.25,4500, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# RIVERS 1.3---------------------------------------------------------

# AVERAGE
par(mfrow = c(2, 2))
plot(Rivers1.3$WEIGHT_g ~Rivers1.3$avgPDO,col= "grey50"); mod1 <- lm(Rivers1.3$WEIGHT_g ~Rivers1.3$avgPDO)
abline(mod1, col="blue",lty=2); tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$WEIGHT_g ~Rivers1.3$avgALPI,col= "grey50")
mod2 <- lm(Rivers1.3$WEIGHT_g ~Rivers1.3$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(2,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$WEIGHT_g ~Rivers1.3$avgSSTn,col= "grey50")
mod3 <- lm(Rivers1.3$WEIGHT_g ~Rivers1.3$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$WEIGHT_g ~Rivers1.3$avgSSTb,col= "grey50")
mod4 <- lm(Rivers1.3$WEIGHT_g ~Rivers1.3$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(6.6,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY -------

par(mfrow = c(2, 2))
plot(Rivers1.3$WEIGHT_g ~Rivers1.3$PDOentry,col= "grey50")
mod1 <- lm(Rivers1.3$WEIGHT_g ~Rivers1.3$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$WEIGHT_g ~Rivers1.3$ALPIentry,col= "grey50")
mod2 <- lm(Rivers1.3$WEIGHT_g ~Rivers1.3$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(6,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$WEIGHT_g ~Rivers1.3$nSSTentry,col= "grey50")
mod3 <- lm(Rivers1.3$WEIGHT_g ~Rivers1.3$nSSTentry); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$WEIGHT_g ~Rivers1.3$bSSTentry, col= "grey50")
mod4 <- lm(Rivers1.3$WEIGHT_g ~Rivers1.3$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,1000, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ----------------------------------------------------------------------
# REGRESSION GRAPHS ------------- LENGTH --------------------------------

# Nass. Riv 1.2's 
# AVERAGE
par(mfrow = c(2, 2))
plot(Nass1.2$LENGTH_mm ~ Nass1.2$avgPDO,col= "grey50"); mod1 <- lm(Nass1.2$LENGTH_mm ~ Nass1.2$avgPDO)
abline(mod1, col="blue",lty=2); tab <- tidy(mod1)
b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LENGTH_mm ~ Nass1.2$avgALPI,col= "grey50")
mod2 <- lm(Nass1.2$LENGTH_mm ~ Nass1.2$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(4,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LENGTH_mm ~ Nass1.2$avgSSTn,col= "grey50")
mod3 <- lm(Nass1.2$LENGTH_mm ~ Nass1.2$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.5,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LENGTH_mm ~ Nass1.2$avgSSTb,col= "grey50")
mod4 <- lm(Nass1.2$LENGTH_mm ~ Nass1.2$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY

par(mfrow = c(2, 2))
plot(Nass1.2$LENGTH_mm ~ Nass1.2$PDOentry,col= "grey50")
mod1 <- lm(Nass1.2$LENGTH_mm ~ Nass1.2$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LENGTH_mm ~ Nass1.2$ALPIentry,col= "grey50")
mod2 <- lm(Nass1.2$LENGTH_mm ~ Nass1.2$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(5.5,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LENGTH_mm ~ Nass1.2$nSSTentry,col= "grey50")
mod3 <- lm(Nass1.2$LENGTH_mm ~ Nass1.2$nSSTentry); abline(mod3, col="blue",lty=2); 
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.2$LENGTH_mm ~ Nass1.2$bSSTentry,col= "grey50")
mod4 <- lm(Nass1.2$LENGTH_mm ~ Nass1.2$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# Nass. Riv 2.2's  ---------------------------------------------------------------
# AVERAGE
par(mfrow = c(2, 2))
plot(Nass2.2$LENGTH_mm ~ Nass2.2$avgPDO,col= "grey50"); mod1 <- lm(Nass2.2$LENGTH_mm ~ Nass2.2$avgPDO)
abline(mod1, col="blue",lty=2); tab <- tidy(mod1)
b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$LENGTH_mm ~ Nass2.2$avgALPI,col= "grey50")
mod2 <- lm(Nass2.2$LENGTH_mm ~ Nass2.2$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(4,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$LENGTH_mm ~ Nass2.2$avgSSTn,col= "grey50")
mod3 <- lm(Nass2.2$LENGTH_mm ~ Nass2.2$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$LENGTH_mm ~ Nass2.2$avgSSTb,col= "grey50")
mod4 <- lm(Nass2.2$LENGTH_mm ~ Nass2.2$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY

par(mfrow = c(2, 2))
plot(Nass2.2$LENGTH_mm ~ Nass2.2$PDOentry,col= "grey50")
mod1 <- lm(Nass2.2$LENGTH_mm ~ Nass2.2$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$LENGTH_mm ~ Nass2.2$ALPIentry,col= "grey50")
mod2 <- lm(Nass2.2$LENGTH_mm ~ Nass2.2$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(6,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$LENGTH_mm ~ Nass2.2$nSSTentry,col= "grey50")
mod3 <- lm(Nass2.2$LENGTH_mm ~ Nass2.2$nSSTentry); abline(mod3, col="blue",lty=2); 
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.2$LENGTH_mm ~ Nass2.2$bSSTentry,col= "grey50")
mod4 <- lm(Nass2.2$LENGTH_mm ~ Nass2.2$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,800, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# Nass. Riv 2.3's  ---------------------------------------------------------------
# AVERAGE
par(mfrow = c(2, 2))
plot(Nass2.3$LENGTH_mm ~ Nass2.3$avgPDO,col= "grey50"); mod1 <- lm(Nass2.3$LENGTH_mm ~ Nass2.3$avgPDO)
abline(mod1, col="blue",lty=2); tab <- tidy(mod1)
b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.2,550, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LENGTH_mm ~ Nass2.3$avgALPI,col= "grey50")
mod2 <- lm(Nass2.3$LENGTH_mm ~ Nass2.3$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(2.6,550, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LENGTH_mm ~ Nass2.3$avgSSTn,col= "grey50")
mod3 <- lm(Nass2.3$LENGTH_mm ~ Nass2.3$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.3,550, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LENGTH_mm ~ Nass2.3$avgSSTb,col= "grey50")
mod4 <- lm(Nass2.3$LENGTH_mm ~ Nass2.3$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(6.8,550, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY

par(mfrow = c(2, 2))
plot(Nass2.3$LENGTH_mm ~ Nass2.3$PDOentry,col= "grey50")
mod1 <- lm(Nass2.3$LENGTH_mm ~ Nass2.3$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,550, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LENGTH_mm ~ Nass2.3$ALPIentry,col= "grey50")
mod2 <- lm(Nass2.3$LENGTH_mm ~ Nass2.3$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(5.5,550, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LENGTH_mm ~ Nass2.3$nSSTentry,col= "grey50")
mod3 <- lm(Nass2.3$LENGTH_mm ~ Nass2.3$nSSTentry); abline(mod3, col="blue",lty=2); 
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,550, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass2.3$LENGTH_mm ~ Nass2.3$bSSTentry,col= "grey50")
mod4 <- lm(Nass2.3$LENGTH_mm ~ Nass2.3$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,550, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# Nass. Riv 1.3's  ---------------------------------------------------------------
# AVERAGE
par(mfrow = c(2, 2))
plot(Nass1.3$LENGTH_mm ~ Nass1.3$avgPDO,col= "grey50"); mod1 <- lm(Nass1.3$LENGTH_mm ~ Nass1.3$avgPDO)
abline(mod1, col="blue",lty=2); tab <- tidy(mod1)
b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1,750, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$LENGTH_mm ~ Nass1.3$avgALPI,col= "grey50")
mod2 <- lm(Nass1.3$LENGTH_mm ~ Nass1.3$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(2.5,750, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$LENGTH_mm ~ Nass1.3$avgSSTn,col= "grey50")
mod3 <- lm(Nass1.3$LENGTH_mm ~ Nass1.3$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9,750, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$LENGTH_mm ~ Nass1.3$avgSSTb,col= "grey50")
mod4 <- lm(Nass1.3$LENGTH_mm ~ Nass1.3$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(6.8,750, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY

par(mfrow = c(2, 2))
plot(Nass1.3$LENGTH_mm ~ Nass1.3$PDOentry,col= "grey50")
mod1 <- lm(Nass1.3$LENGTH_mm ~ Nass1.3$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,750, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$LENGTH_mm ~ Nass1.3$ALPIentry,col= "grey50")
mod2 <- lm(Nass1.3$LENGTH_mm ~ Nass1.3$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(5.5,750, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$LENGTH_mm ~ Nass1.3$nSSTentry,col= "grey50")
mod3 <- lm(Nass1.3$LENGTH_mm ~ Nass1.3$nSSTentry); abline(mod3, col="blue",lty=2); 
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,750, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Nass1.3$LENGTH_mm ~ Nass1.3$bSSTentry,col= "grey50")
mod4 <- lm(Nass1.3$LENGTH_mm ~ Nass1.3$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,750, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# RIVERS INLET---------------------------------------------------------

# 1.2's
# AVERAGE
par(mfrow = c(2, 2))
plot(Rivers1.2$LENGTH_mm ~Rivers1.2$avgPDO,col= "grey50")
mod1 <- lm(Rivers1.2$LENGTH_mm ~Rivers1.2$avgPDO); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,650, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$LENGTH_mm ~Rivers1.2$avgALPI,col= "grey50")
mod2 <- lm(Rivers1.2$LENGTH_mm ~Rivers1.2$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(3,650, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$LENGTH_mm ~Rivers1.2$avgSSTn,col= "grey50")
mod3 <- lm(Rivers1.2$LENGTH_mm ~Rivers1.2$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,650, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$LENGTH_mm ~Rivers1.2$avgSSTb,col= "grey50")
mod4 <- lm(Rivers1.2$LENGTH_mm ~Rivers1.2$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,650, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY -------

par(mfrow = c(2, 2))
plot(Rivers1.2$LENGTH_mm ~Rivers1.2$PDOentry,col= "grey50")
mod1 <- lm(Rivers1.2$LENGTH_mm ~Rivers1.2$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.75,650, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$LENGTH_mm ~Rivers1.2$ALPIentry,col= "grey50")
mod2 <- lm(Rivers1.2$LENGTH_mm ~Rivers1.2$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(5.5,650, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$LENGTH_mm ~Rivers1.2$nSSTentry,col= "grey50")
mod3 <- lm(Rivers1.2$LENGTH_mm ~Rivers1.2$nSSTentry); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(10,650, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.2$LENGTH_mm ~Rivers1.2$bSSTentry,col= "grey50")
mod4 <- lm(Rivers1.2$LENGTH_mm ~Rivers1.2$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7.25,650, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# RIVERS 1.3---------------------------------------------------------

# AVERAGE
par(mfrow = c(2, 2))
plot(Rivers1.3$LENGTH_mm ~Rivers1.3$avgPDO,col= "grey50"); mod1 <- lm(Rivers1.3$LENGTH_mm ~Rivers1.3$avgPDO)
abline(mod1, col="blue",lty=2); tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1,700, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$LENGTH_mm ~Rivers1.3$avgALPI,col= "grey50")
mod2 <- lm(Rivers1.3$LENGTH_mm ~Rivers1.3$avgALPI); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(2,700, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$LENGTH_mm ~Rivers1.3$avgSSTn,col= "grey50")
mod3 <- lm(Rivers1.3$LENGTH_mm ~Rivers1.3$avgSSTn); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9,700, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$LENGTH_mm ~Rivers1.3$avgSSTb,col= "grey50")
mod4 <- lm(Rivers1.3$LENGTH_mm ~Rivers1.3$avgSSTb); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(6.6,700, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

# ENTRY -------

par(mfrow = c(2, 2))
plot(Rivers1.3$LENGTH_mm ~Rivers1.3$PDOentry,col= "grey50")
mod1 <- lm(Rivers1.3$LENGTH_mm ~Rivers1.3$PDOentry); abline(mod1, col="blue",lty=2)
tab <- tidy(mod1); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod1)$adj.r.squared
text(1.5,700, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$LENGTH_mm ~Rivers1.3$ALPIentry,col= "grey50")
mod2 <- lm(Rivers1.3$LENGTH_mm ~Rivers1.3$ALPIentry); abline(mod2, col="blue",lty=2)
tab <- tidy(mod2); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod2)$adj.r.squared
text(6,700, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$LENGTH_mm ~Rivers1.3$nSSTentry,col= "grey50")
mod3 <- lm(Rivers1.3$LENGTH_mm ~Rivers1.3$nSSTentry); abline(mod3, col="blue",lty=2)
tab <- tidy(mod3); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod3)$adj.r.squared
text(9.75,700, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))

plot(Rivers1.3$LENGTH_mm ~Rivers1.3$bSSTentry, col= "grey50")
mod4 <- lm(Rivers1.3$LENGTH_mm ~Rivers1.3$bSSTentry); abline(mod4, col="blue",lty=2)
tab <- tidy(mod4); b <- tab[2,2]; p <- tab[2,5]; r <- summary(mod4)$adj.r.squared
text(7,700, paste("R =", round(r, 4), "\nb =", round(b, 4), "\np-val =", round(p, 5)))



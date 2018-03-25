# CONTEMPORARY NASS SOCKEYE DATA CLEAN
# Created by C Freshwater MAR 22, 2018
# Last revised: ONGOING
# Cleans fishwheel data provided by LGL, DFO North Coast
# and Nisqaa
# Note: Due to data sharing agreements these data are not
# publicly available and will NOT be posted to gitHub, hence
# call to local directory
# -------------------------------------------------

library(mgcv); library(dplyr); library(ggplot2); library(reshape2); library(here)

origModDat <- read.csv("C:/Users/FRESHWATERC/Documents/SideProjects/histSockeye/privData/nassModernData.csv", stringsAsFactors=F)
nassDat <- read.csv(here("data/nassDat.csv"))

## ---------------------- Clean ------------------------------
keepAges <- c("42", "52", "53", "63")
origModDat <- origModDat[origModDat$age %in% keepAges,]
origModDat$age <- as.factor(origModDat$age)
levels(origModDat$age) <- c(levels(origModDat$age), "1.2", "1.3", "2.2", "2.3")
origModDat$age[origModDat$age == "42"] <- "1.2"
origModDat$age[origModDat$age == "52"] <- "1.3"
origModDat$age[origModDat$age == "53"] <- "2.2"
origModDat$age[origModDat$age == "63"] <- "2.3"
origModDat$age <- factor(origModDat$age)

dates <- strptime(origModDat$date, "%d-%b-%y")

modDat <- data.frame(retYr = as.numeric(origModDat$year),
                     month = as.numeric(format(strptime(origModDat$date, "%d-%b-%y"), "%m")),
                     day = as.numeric(format(strptime(origModDat$date, "%d-%b-%y"), "%d")),
                     jDay = as.numeric(format(strptime(origModDat$date, "%d-%b-%y"), "%j")),
                     sex = origModDat$sex,
                     age = origModDat$age,
                     fl = as.numeric(origModDat$flCM)*10
                     )
modDat<-modDat[!is.na(modDat$fl),]

# Outlier check
ggplot(modDat, aes(x=fl)) +
  geom_histogram(position="identity") +
  facet_grid(~factor(age))
# looks good


## ---------------------- Exploratory Comparisons With Contemporary ------------------------------
fullDate <- as.Date(paste(nassDat$DAY, nassDat$MONTH, nassDat$retYr, sep="-"), format="%d-%m-%Y")
jDay <- format(fullDate, "%j")
histDat <- data.frame(retYr = as.numeric(nassDat$retYr),
                      month = as.numeric(paste(0, nassDat$MONTH, sep="")),
                      day = as.numeric(nassDat$DAY),
                      jDay = as.numeric(jDay)+6,
                      sex = nassDat$SEX,
                      age = nassDat$AGE,
                      fl = nassDat$LENGTH_mm
)
nassFull <- rbind(modDat, histDat)
nassFull$dataSet <- as.factor(ifelse(nassFull$retYr < 1950, "hist", "mod"))


# Changes in size
ggplot(nassFull, aes(x = as.numeric(retYr), y = fl, colour=dataSet)) + 
  geom_line() + 
  facet_wrap(~ age)


# Changes in age at maturity
temp <- nassFull[,c("retYr", "age")]
temp2 <- prop.table(table(temp), 1)
nassAge <- as.data.frame(melt(temp2))
names(nassAge) <- c("retYr", "age", "ppn")
ggplot(nassAge, aes(x = retYr, y = ppn, fill = as.factor(age)))  +
  geom_area(position = 'stack') +
  xlab("Return Year") +
  ylab("Proportion") +
  scale_fill_discrete(name = "age")


# Return timing
ggplot(nassFull, aes(x=jDay, fill=dataSet)) +
  geom_histogram(position="identity", colour="grey40", alpha=0.2) +
  facet_grid(~factor(age))


## ---------------------- Add Environmental Data ------------------------------

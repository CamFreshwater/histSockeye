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

modDat <- data.frame(retYr = origModDat$year,
                     month = format(strptime(origModDat$date, "%d-%b-%y"), "%m"),
                     day = format(strptime(origModDat$date, "%d-%b-%y"), "%d"),
                     sex = origModDat$sex,
                     age = origModDat$age,
                     fl = as.numeric(origModDat$flCM)*10
                     )


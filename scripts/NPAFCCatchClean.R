# NPAFC CATCH DATA CLEAN
# Created by C Freshwater April 10, 2010
# Last revised: ONGOING
# Cleans catch data from BC, WA, and AK for use in fullSockeyeModels.R
# Original source NPAFC_Catch_Stat_20July2016.xls 
# (downloaded March 31, 2018)
# -------------------------------------------------


setwd("/Users/cam/github")

library(reshape2); library(here); library(tidyr)

catchDatWide <- read.csv(here("github/histSockeye/data/rawCatchData/NPAFC_Catch_Stat_20July2016.csv"), stringsAsFactors=F)

# Convert from wide to long
catchDatLong <- melt(catchDatWide, variable.name = "year", value.name = "catch",
				id.vars = c("Country", "Whole.Country.Province.State", "Reporting.Area", "Species", "Catch.Type", "Data.Type"))
trimCatch <- catchDatLong[,-1]
colnames(trimCatch)[1:5] <- c("region", "area", "species", "fishery", "metric") 
trimCatch$year <- as.factor(substring(trimCatch$year, 2)) #remove introduced x character
trimCatch <- trimCatch[trimCatch$metric == "Number (000's)",] #remove biomass counts
trimCatch <- trimCatch[,-5] #remove metric column



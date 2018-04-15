# NPAFC CATCH DATA CLEAN
# Created by C Freshwater April 10, 2010
# Last revised: ONGOING
# Cleans catch data from BC, WA, and AK for use in fullSockeyeModels.R
# Original source NPAFC_Catch_Stat_20July2016.xls 
# (downloaded March 31, 2018)
# -------------------------------------------------


setwd("/Users/cam/github")

library(reshape2); library(here); library(dplyr); library(ggplot2)


### Import and clean modern catch data
catchDatWide <- read.csv(here("github/histSockeye/data/rawCatchData/NPAFC_Catch_Stat_20July2016.csv"), stringsAsFactors=F)

# Convert from wide to long
catchDatLong <- melt(catchDatWide, variable.name = "year", value.name = "catch",
				id.vars = c("Country", "Whole.Country.Province.State", "Reporting.Area", "Species", "Catch.Type", "Data.Type"))
trimCatch <- catchDatLong[,-1]
colnames(trimCatch)[1:5] <- c("region", "area", "species", "fishery", "metric") 
trimCatch$year <- as.numeric(substring(trimCatch$year, 2)) #remove introduced x character
trimCatch <- trimCatch[trimCatch$metric == "Number (000's)",] #remove biomass counts
trimCatch$catch <- as.numeric(trimCatch$catch) * 1000
trimCatch <- trimCatch[,-5] #remove metric column
trimCatch <- trimCatch[trimCatch$year > 1950,] #drop earlier years that will be redunant with INPFC
trimCatch <- trimCatch[trimCatch$fishery == "Commercial",] #drop all fisheries except commercial to be consistent


aggCatch <- trimCatch %>%
	group_by(region, species, fishery, year) %>%
	summarise(totalCatch = sum(catch, na.rm=TRUE))
aggCatch$data <- "npafc"


### Import and clean historic catch data
oldCatchDat <- read.csv(here("github/histSockeye/data/rawCatchData/INPFC_CatchData_April2018.csv"), stringsAsFactors=F)

oldCatchDat$year <- as.numeric(oldCatchDat$year)
oldCatchDat$catch <- as.numeric(oldCatchDat$catch)

aggOldCatch <- oldCatchDat %>%
	group_by(region, species, fishery, year) %>%
	summarise(totalCatch = sum(catch, na.rm=TRUE))
aggOldCatch$data <- "inpfc"

fullCatch <- rbind(aggCatch, aggOldCatch)
fullCatch[fullCatch$totalCatch == 0,]$totalCatch <- NA

agFullCatch <- fullCatch %>%
	group_by(species, year, data) %>%
	summarise(totalCatch = sum(totalCatch, na.rm=TRUE))
# agFullCatch[agFullCatch$totalCatch == 0,]$totalCatch <- NA

# Catches line up well
ggplot(agFullCatch, aes(x = as.numeric(year), y = totalCatch)) + 
    geom_line(aes(colour = data)) + 
    facet_wrap(~ species)


write.csv(fullCatch, here("github/histSockeye/data/cleanCatch.csv"))
write.csv(agFullCatch, here("github/histSockeye/data/cleanAgCatch.csv"))
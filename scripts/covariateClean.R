# NEW ENVIRONMENTAL DATA CLEAN
# Created by C Freshwater MAR 25, 2018
# Last revised: ONGOING
# Cleans full environmental dataset; based off of NASS_Jan2018_CFedits.R script,
# but revised to include contemporary, as well as historical data
# ------------------------------------------------------------------------------


library(dplyr); library(here)

sst <- read.table(here("github/histSockeye/data/scrnep.txt"))


## ---------------------- Clean ------------------------------

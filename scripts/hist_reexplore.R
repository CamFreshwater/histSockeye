## HISTORICAL SOCKEYE ANALYSIS - REBOOT
## Created by C Freshwater Dec. 4, 2019

# -Stripped down version of fullSockeyeModels
# -Uses data prepped in modernNassClean.R (merges historical and contemporary
#  sockeye data) and NASS_Jan2018_CFedits.R (all else); note these should be 
#  revisited depending on the analyses that are chosen


library(tidyverse); library(ggplot2); library(car); library(reshape2)

source(here::here("scripts/histSockeyeFunc.R"))

#principal components of SST variation in NE Pacific (170E to 240E, 40N-65N)
# sstPca <- read.table(here::here("data/sstPCA.txt"))
# sstRaw <- read.table(here::here("data/sstRaw.txt")) # pacific ocean SST
# pdo <- read.csv(here::here("data/pdo.csv"), stringsAsFactors=F)
# #stops at 2015
# meanAlpi <- read.csv(here::here("data/alpi.csv"), stringsAsFactors=F)
# sockDat <- read.csv(here::here("data/nassFullSox.csv"), stringsAsFactors=F)
# catchDat <- read.csv(here::here("data/akCatch.csv"), stringsAsFactors=F)

full <- readRDS(here::here("data", "combinedData.RDS")) %>% 
  mutate(index = paste(watershed, dataSet, sep = "_"))

#EXPLORATORY -------------------------------------------------------------------
# sockDat %>% 
#   group_by(dataSet, watershed) %>% 
#   tally()

# Plot changes in size among indices
mean_dat <- full %>%
  group_by(retYr, age, index) %>%
  summarize(meanFL = mean(fl), pdo = mean(pdo), alpi = mean(alpi), 
            rawSst = mean(temp), pcaSst = mean(pc2), pink = mean(pinkCatch), 
            sox = mean(sockCatch)) 

lin_size_trend <- ggplot(mean_dat, aes(x = as.numeric(retYr), y = meanFL, 
                                      col = as.factor(age))) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~index, scales = "free_x")

pdf(here::here("outputs", "figs", "linear_size_trend.pdf"))
lin_size_trend
dev.off()

# Plot changes in age structure
age_ppn <- full %>% 
  group_by(index, retYr, age) %>% 
  summarise(n = n()) %>% 
  mutate(ageFreq = n / sum(n))

age_comp <- ggplot(age_ppn) +
  geom_area(aes(x = retYr, y = ageFreq, fill = as.factor(age)), 
            position = 'stack') +
  xlab("Return Year") +
  ylab("Proportion") +
  scale_fill_discrete(name = "age") +
  facet_wrap(~index, scales = "free_x")

pdf(here::here("outputs", "figs", "age_comp.pdf"))
age_comp
dev.off()

## Quick comparison of ocean entry timing by age based on otolith data
# Passed to Skip for his evolutionary stable strategies paper
# April 12, 2020

library(tidyverse); library(plyr)

DataJS <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/JohnstoneStrait/JSdatClean.csv", 
                   stringsAsFactors = FALSE, strip.white = TRUE, 
                   na.strings = c("NA","")) %>% 
  select(FishNumber, Year, Lat, Long, ShipFL, MEC_Radius = MECRadius, 
         Total_Radius = TotalRadius, Age, 
         EntryDate, Stock) %>% 
  mutate(Data = "JS")
Fall08 <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/SockeyeFallMigrants/2007-08OtolithsNew.csv", 
                   stringsAsFactors = FALSE, strip.white = TRUE, 
                   na.strings = c("NA","")) %>% 
  filter(Use == "Y",
         Season == "Fall") %>% 
  mutate(entry.date = JulianDate - Total_Count,
         Data = "Fall08") %>% 
  select(FishNumber, Year, Lat, Long, ShipFL, MEC_Radius, Total_Radius, Age, 
         EntryDate = entry.date, Stock, Data)
SoG1112 <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/SoGPurse/SoGPurseFinalData.csv", 
                    stringsAsFactors = FALSE, strip.white = TRUE, 
                    na.strings = c("NA","")) %>% 
  select(FishNumber = UniversalNumber, Year, Lat, Long, ShipFL, 
         MEC_Radius = MECRadius, Total_Radius = TotalRadius, Age, 
         EntryDate, Stock) %>% 
  mutate(Data = "SoG1112")
Coast0708 <- read.csv("/Users/cam/cam\ stuff/Grad\ School/JuvSockeyeData/Analysis/SockeyeLatGradient/LatGradientOtolithCF.csv", 
                      stringsAsFactors = FALSE, strip.white = TRUE, 
                      na.strings = c("NA","")) %>% 
  select(FishNumber, Year, Lat = lat, Long = long, ShipFL, MEC_Radius, Total_Radius, Age, 
         EntryDate = entry.date, Stock) %>% 
  mutate(Data = "Coast0708")


#correct stock names
all_dat <- rbind(DataJS, Fall08, SoG1112, Coast0708) %>% 
  mutate(Stock = case_when(
    Stock %in% c("CH", "Chilko", "Chilko_south") ~ "Chilko",
    Stock %in% c("L_Adams", "L_Shuswap", "LA", "MiddleShuswap", "Scotch", 
                 "Seymour") ~ "Shuswap",
    Stock == "GC" ~ "Great Central",
    Stock == "SP" ~ "Sproat",
    grepl("Ques", Stock) ~ "Quesnel",
    TRUE ~ Stock
  ))
saveRDS(all_dat, here::here("outputs", "data", "all_otoliths.rds"))

stock_w_age2 <- all_dat %>% 
  filter(Age == 2) %>% 
  pull(Stock) %>% 
  unique()

# Compare entry date by ages and stock
phen_age <- all_dat %>% 
  filter(Stock %in% stock_w_age2,
         !Data == "Fall08") %>% 
  ggplot(.) + 
  geom_boxplot(aes(x = as.factor(Age), y = EntryDate)) +
  facet_wrap(~Stock)
phen_age

pdf(here::here("outputs", "figs", "age_phenology.pdf"))
phen_age
dev.off()

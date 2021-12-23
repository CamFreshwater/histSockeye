# Data clean
# Import Skip's first century of data and clean as necessary
# Created by C Freshwater Dec 23, 2021
# -------------------------------------------------

library(tidyverse)


dat_in <- read.table(here::here("data", "nasscentury.txt"))
colnames(dat_in) <- c("year", "month", "day", "hart_species", "sex", "age", "fl")
dat <- dat_in %>% 
  mutate(
    sex = ifelse(sex == 1, "male", "female"),
    age = as.character(age)
  ) 


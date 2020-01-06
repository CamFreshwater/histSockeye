## Selectivity analysis
# Jan. 5, 2019

library(tidyverse)

dat <- foreign::read.systat(here::here("data", "selectivityData", 
                                       "AREA3sockeye.SYD")) %>% 
  mutate(YEAR.WEEK = paste(YEAR, WEEK, sep = "_")) %>% 
  group_by(YEAR.WEEK) %>%
  #replace missing values with next/prev (opposite doesn't change results)
  fill(GEAR., .direction = "downup")

## Identify years w/ both gear types
both_gear <- dat %>% 
  #remove years with only one gear type
  group_by(YEAR.WEEK) %>%
  filter(!length(unique(GEAR.)) < 2) %>% 
  tally() %>% 
  filter(n < 10)
yrs <- unique(both_gear$YEAR.WEEK)

trim_dat <- dat %>% 
  filter(YEAR.WEEK %in% yrs) %>% 
  select(YEAR, WEEK, AGE., ESTFORK, GEAR., YEAR.WEEK)

ggplot(trim_dat) +
  geom_histogram(aes(x = ESTFORK, fill = GEAR.)) +
  facet_wrap(~YEAR.WEEK)

## Hypural-fork length comparison
# Some years only have HL data; what are consequences of using transformed FL
# data based on these values
# March 24, 2020
# -------------------------------------------------

library(tidyverse)
library(lme4)

fl_dat <- read.csv(here::here("data", "hyp_FL_comparison_bilton.csv")) %>% 
  mutate(year_f = as.factor(YEAR),
         hyp_z = scale(HYP_LEN))

ggplot(fl_dat) +
  geom_point(aes(x = hyp_z, y = FORK_LENGTH)) +
  facet_wrap(~year_f)

mod <- lmer(FORK_LENGTH ~ hyp_z + (1 | year_f), data = fl_dat)
modb <- lmer(FORK_LENGTH ~ hyp_z + (1 + hyp_z | year_f), data = fl_dat)
mod2 <- lm(FORK_LENGTH ~ hyp_z, data = fl_dat)

summary(mod)
summary(mod2)

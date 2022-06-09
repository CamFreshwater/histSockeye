# Covariance Analysis - Online Supplement 2
# June 3, 2022
# Identify correlations among cohorts 
# ------------------------------------------------------------------------------

library(tidyverse)
library(ggcorrplot)


# posterior predictions from brms_fit
post <- readRDS(
  here::here("outputs", "data", "brms_fits", "brms_post_preds.rds")) 


# year means
dat_avg <- read.table(here::here("data", "nasscenturyavgv3FLAT.txt"),
                      header = TRUE) %>%
  pivot_longer(cols = -(iy.)) %>%
  mutate(
    ret_year = iy.,
    year_f = as.factor(iy.),
    sex = ifelse(grepl("m", name), "male", "female") %>% as.factor,
    age = case_when(
      grepl("42", name) ~ "42",
      grepl("52", name) ~ "52",
      grepl("53", name) ~ "53",
      grepl("63", name) ~ "63"
    ) %>%
      as.factor(),
    brood_year = case_when(
      age == "42" ~ ret_year - 4,
      age == "63" ~ ret_year - 6,
      grepl("5", age) ~ ret_year -5
    ),
    entry_year = case_when(
      grepl("2", age) ~ brood_year + 2,
      grepl("3", age) ~ brood_year + 3
    ),
    period = case_when(
      iy. < 1957 ~ "Gilbert-Clemens",
      iy. > 1972 & iy. < 1994 ~ "Monkley Dump",
      iy. > 1993 ~ "Nisga'a",
      TRUE ~ "Bilton"
    ),
    age_f = age,
    period = fct_relevel(as.factor(period),
                         "Gilbert-Clemens",
                         "Bilton",
                         "Monkley Dump",
                         "Nisga'a"),
    value = value / 10
    ) %>% 
  rename(obs_mean = value) %>% 
  distinct()
  
# combine with models est
dat_avg2 <- left_join(dat_avg, 
                      post %>% 
                        filter(year %in% dat_avg$ret_year) %>% 
                        select(ret_year = year,
                               age,
                               sex,
                               median_est = median) %>% 
                        distinct(),
                      by = c("age", "sex", "ret_year"))
levels(dat_avg2$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")


# CORRELATION MATRICES ---------------------------------------------------------

# example with just return year and observed correlations
corr_foo <- function(data, sex_in = "male", year_class, data_class,
                     title) {
  dum <- data %>% 
    filter(sex == sex_in) %>% 
    select({{ year_class }}, age_f, {{ data_class }}) %>% 
    pivot_wider(names_from = age_f, 
                values_from = {{ data_class }}) %>%
    drop_na() %>% 
    select(-{{ year_class }}) %>% 
    cor(.) 

  return(dum)
    
  # ggcorrplot(dum, hc.order = TRUE, type = "lower",
  #            lab = TRUE, title = title) 
}


dat_avg2 %>% 
  filter(sex == sex_in) %>% 
  select({{ year_class }}, age_f, {{ data_class }})
pivot_wider(names_from = age_f, 
            values_from = {{ data_class }}) %>%
  drop_na() %>% 
  select(-{{ year_class }}) %>% 
  cor(.) 


pdf(here::here("outputs", "figs", "cor_plot.pdf"))
corr_foo(dat_avg2, year_class = ret_year, data_class = obs_mean,
         title = "Return Year")
corr_foo(dat_avg2, year_class = brood_year, data_class = obs_mean,
         title = "Brood Year")
corr_foo(dat_avg2, year_class = entry_year, data_class = obs_mean,
         title = "Entry Year")
dev.off()

corr_foo(dat_avg2, year_class = ret_year, data_class = median_est)
corr_foo(dat_avg2, year_class = brood_year, data_class = median_est)
corr_foo(dat_avg2, year_class = entry_year, data_class = median_est)


dum1 <- dat_avg2 %>% 
  filter(sex == "male",
         age_f %in% c("4[2]", "5[3]")) %>%
  select(ret_year, age_f, obs_mean) %>% 
  pivot_wider(names_from = age_f, 
              values_from = obs_mean) %>%
  drop_na() %>% 
  select(-ret_year) %>% 
  cor(.)
dum2 <- dat_avg2 %>% 
  filter(sex == "male",
         age_f %in% c("4[2]", "5[3]")) %>%
  select(entry_year, age_f, obs_mean) %>% 
  pivot_wider(names_from = age_f, 
              values_from = obs_mean) %>%
  drop_na() %>% 
  select(-entry_year) %>%
  cor(.)

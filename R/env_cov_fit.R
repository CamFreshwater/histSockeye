# Covariate analysis
# For context explore patterns in residuals after including PDO and salmon
# abundance (from NPAFC summary statistics) as covariate
# NOTE: post-hoc for presentation, not included in published manuscript
# April 6, 2023
# ------------------------------------------------------------------------------


library(tidyverse)

# import nass data
dat_in <- read.table(here::here("data", "nasscenturyv3.txt"))
colnames(dat_in) <- c("year", "month", "day", "hart_species", "sex", "age", 
                      "fl")

# import pdo
monthly_pdo <- read.csv(here::here("data", "monthly_pdo.csv")) %>% 
  pivot_longer(cols = c(Jan:Dec), names_to = "month") %>% 
  janitor::clean_names() 
pdo <- monthly_pdo %>% 
  group_by(year) %>% 
  summarize(mean_pdo = mean(value)) %>% 
  mutate(mean_pdo_z = scale(mean_pdo) %>% as.numeric) %>% 
  glimpse()

# import salmon abund
salmon <- read.csv(here::here("data", "npafc_abundance.csv")) %>% 
  mutate(
    pink_z = scale(pink) %>% as.numeric(),
    chum_z = scale(chum) %>% as.numeric(),
    sockeye_z = scale(sockeye) %>% as.numeric(),
    total_z = scale(total) %>% as.numeric()
  )

# correlations among salmon species
library(ggcorrplot)
corr <- cor(salmon %>% select(pink_z, chum_z, sockeye_z))



set.seed(123)
dat <- dat_in %>%
  mutate(
    sex = ifelse(sex == 1, "male", "female") %>% as.factor(),
    age = factor(as.character(age)),
    year_f = as.factor(year),
    age_f = age,
    date = as.POSIXct(paste(day, month, year, sep = "-"),
                      format = "%d-%m-%Y"),
    yday = lubridate::yday(date),
    yday_c = yday - mean(yday, na.rm = TRUE)
  ) %>% 
  #drop secondary age classes and years missing capture data
  filter(
    age %in% c("42", "52", "53", "63"),
    !is.na(yday),
    !is.na(fl)
  )  %>%
  # Add a factor representing sample period 
  mutate(
    period = case_when(
      year < 1957 ~ "Gilbert-Clemens",
      year > 1972 & year < 1994 ~ "Monkley Dump",
      year > 1993 ~ "NFWD",
      TRUE ~ "Bilton"
    ),
    period = fct_relevel(as.factor(period), 
                         "Gilbert-Clemens",
                         "Bilton",
                         "Monkley Dump",
                         "NFWD"),
    age_sex = paste(sex, age, sep = "_") %>% as.factor()
  ) %>% 
  droplevels()

mean_dat <- dat %>% 
  group_by(year_f, age_sex) %>% 
  summarize(obs_mean_fl = mean(fl))


# fit hierarchical gam to calculate annual means by age/sex accounting for 
# changes in sampling date
gam_fit <- mgcv::gam(
  fl ~ s(yday_c, by = age_sex, m = 2) + age + sex + s(year_f, bs = "re"),
  data = dat
)

new_dat <- expand.grid(
  year = unique(dat$year),
  age = unique(dat$age),
  sex = unique(dat$sex),
  yday_c = 0
) %>% 
  mutate(year_f = as.factor(year),
         age_sex = paste(sex, age, sep = "_") %>% as.factor()) %>% 
  left_join(., mean_dat, by = c("year_f", "age_sex")) %>% 
  left_join(., pdo %>% select(-mean_pdo), by = "year") %>%
  left_join(., salmon %>% select(year = return_year, ends_with("_z"))) 

preds <- predict(gam_fit, newdata = new_dat)
new_dat$pred_fl <- preds %>% as.numeric()

ggplot(new_dat %>% filter(sex == "female")) +
  geom_point(aes(x = year, y = obs_mean_fl), shape = 21, fill = "red") +
  geom_point(aes(x = year, y = pred_fl), shape = 21, fill = "black") +
  facet_wrap(~age)

# fit second gam to estimate effects of pdo and salmon
new_dat_trim <- new_dat %>% filter(!is.na(pink_z))
gam_fit2 <- mgcv::gam(
  pred_fl ~ s(mean_pdo_z, k = 3) + s(total_z, k = 3) + age_sex,
  data = new_dat_trim
)
new_dat_trim$resids <- resid(gam_fit2)


png(here::here("outputs", "figs", "cov_gam_resids.png"), width = 4, height = 4,
    res = 250, units = "in")
ggplot(new_dat_trim %>% filter(sex == 'female', age == "42")) +
  geom_point(aes(x = year, y = resids), fill = "grey", shape = 21) +
  labs(x = "Year", y = "Residuals") +
  geom_hline(yintercept = 0, lty = 2, colour = "red") +
  ggsidekick::theme_sleek()
dev.off()

pdo_preds <- data.frame(
  mean_pdo_z = seq(-1, 1, length = 100),
  age_sex = "female_42",
  total_z = 0,
  var = "pdo"
)
salmon_preds <- data.frame(
  total_z = seq(-1, 1, length = 100),
  age_sex = "female_42",
  mean_pdo_z = 0,
  var = "salmon abundance"
)
new_dat2 <- rbind(pdo_preds, salmon_preds)
preds2 <- predict(gam_fit2, newdata = new_dat2, se.fit = TRUE)
new_dat2 <- new_dat2 %>% 
  mutate(
    pred_fl = preds2$fit,
    pred_fl_se = preds2$se.fit,
    upr = pred_fl + (1.96 * pred_fl_se),
    lwr = pred_fl - (1.96 * pred_fl_se)
  )
ymin <- min(new_dat2$lwr)
ymax <- max(new_dat2$upr)

p1 <- ggplot(new_dat2 %>% filter(var == "pdo"), aes(x = mean_pdo_z)) +
  geom_line(aes(y = pred_fl)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  ggsidekick::theme_sleek() +
  lims(y = c(ymin, ymax)) +
  labs(x = "Scaled PDO", y = "Fork Length") 
p2 <- ggplot(new_dat2 %>% filter(var == "salmon abundance"), aes(x = total_z)) +
  geom_line(aes(y = pred_fl)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  ggsidekick::theme_sleek() +
  lims(y = c(ymin, ymax)) +
  labs(x = "Scaled Total\nSalmon Abundance", y = "Fork Length") 


png(here::here("outputs", "figs", "cov_gam_fe.png"), width = 8, height = 4,
    res = 250, units = "in")
cowplot::plot_grid(p1, p2, ncol = 2)  
dev.off()
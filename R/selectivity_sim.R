## Dummy simulation evaluate potential selectivity effects
# Dec. 23, 2022


library(tidyverse)

# import nass data
dat_in <- read.table(here::here("data", "nasscenturyv3.txt"))
colnames(dat_in) <- c("year", "month", "day", "hart_species", "sex", "age", 
                      "fl")

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
  # Add a factor representing sample collection to visualize effects 
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
    bin = cut(fl, 
              breaks = c(seq(400, 600, by = 10), max(fl + 5)),
              right = FALSE)
  ) %>% 
  droplevels()
levels(dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")


# fit bare bones model for estimates of age, sex and year effects
fit_lm <- lm(fl ~ year + age + sex, data = dat)

# generate mean estimates for each age-sex-year
new_dat <- expand.grid(
  year = seq(1914, 2019, by = 1),
  age = unique(dat$age),
  sex = unique(dat$sex)
)
pred_mu <- predict(fit_lm, newdata = new_dat)

abund_dat <- data.frame(
  age = c("42", "52", "53", "63"),
  ppn = c(0.240, 0.181, 0.484, 0.0950)
) %>% 
  mutate(
    abund = 10000 * ppn
  )

# selectivity values for 10 mm size bins based on Kendall et al. 2009
sel_foo <- function(l_i, beta, l_50 = 500) {
  1 / (1 + exp(-beta * (l_i - l_50)))
}

fl_seq = c(seq(400, 600, by = 10), max(dat$fl + 5))
sel_seq = sel_foo(fl_seq, beta = 0.1)
plot(sel_seq ~ fl_seq)
sel_plot_dat <- data.frame(
  fl = fl_seq, 
  selectivity = sel_seq
)

png(here::here("outputs", "figs", "selectivity_curve.png"), 
    height = 7, width = 11, units = "in", res = 200)
ggplot(sel_plot_dat %>% filter(fl < 650)) +
  geom_line(aes(x = fl, y = selectivity)) +
  ggsidekick::theme_sleek()
dev.off()


sel_dat <- data.frame(
  fl_bin =  cut(fl_seq, 
                breaks = c(seq(400, 600, by = 10), max(fl_seq)),
                right = FALSE,
                labels = FALSE),
  sel = sel_seq
) 

new_tbl <- new_dat %>% 
  mutate(
    mu_fl = pred_mu %>% as.numeric()
  ) %>% 
  left_join(., abund_dat, by = "age")

# simulate based on predicted mu and estimated sigma for each age and year
set.seed(2023)
nsims <- 250
dat_list <- vector(length = nsims, mode = "list")

for (i in seq_len(nsims)) {
  new_tbl$sim_fl <- purrr::map2(
    new_tbl$abund,
    new_tbl$mu_fl,
    ~ rnorm(n = .x, mean = .y, sd = summary(fit_lm)$sigma)
  )
  
  # assign bins to simulated data
  dum <- unnest(new_tbl, cols = "sim_fl") %>% 
    mutate(
      fl_bin = cut(sim_fl, 
                   breaks = c(seq(400, 600, by = 10), max(sim_fl + 50)),
                   right = FALSE,
                   labels = FALSE)
    ) %>% 
    left_join(., sel_dat, by = "fl_bin") 
  
  # identify number of samples for each size bin
  pool_dum <- dum %>% 
    group_by(year, fl_bin, sel) %>% 
    tally() %>%
    mutate(
      exp_rate = 0.75,
      sample_rate = 0.05,
      # remove selectivity when year > 1974
      sel = ifelse(year > 1972, 1, sel), 
      n_samps = rbinom(1, n, sel * exp_rate * sample_rate)
    ) 
  
  samp_dat <- dum %>% 
    group_by(year, fl_bin) %>% 
    group_nest(.key = "data") %>% 
    left_join(., 
              pool_dum %>% select(year, fl_bin, n_samps), 
              by = c("year", "fl_bin")) %>% 
    mutate(
      # pull individual fish based on sample size
      samp_fish = purrr::map2(data, n_samps, ~ sample_n(.x, .y))
    ) %>% 
    unnest(cols = samp_fish) 
  
  # fit to samples
  fit_lm2 <- lm(sim_fl ~ year + age + sex, data = samp_dat)
  
  dat_list[[i]] <- data.frame(
    iter = i,
    coef = names(coef(fit_lm2)),
    est = coef(fit_lm2) %>% as.numeric()
  )
  
}

true_coef <- data.frame(
  coef = names(coef(fit_lm)),
  est = summary(fit_lm)$coefficients[ , 1],
  se = summary(fit_lm)$coefficients[ , 2] 
) 
coef_dat <- bind_rows(dat_list) 
saveRDS(coef_dat, here::here("outputs", "model_fits", "sim_pars.rds"))


png(here::here("outputs", "figs", "selectivity_boxplot_pars.png"), 
    height = 7, width = 11, units = "in", res = 200)
ggplot() +
  geom_boxplot(data = coef_dat %>% 
                 mutate(dat = "recovered"), 
               aes(x = coef, y = est, fill = dat)) +
  geom_point(data = true_coef %>% 
                    mutate(dat = "true"), 
                  aes(x = coef, y = est, fill = dat), 
                  position = position_jitterdodge(),
                  shape = 21, size = 2) +
  facet_wrap(~coef, scales = "free") +
  ggsidekick::theme_sleek()
dev.off()


## calculate difference in average size
# true difference
new_dat %>% 
  mutate(
    mu_fl = pred_mu %>% as.numeric()
  ) %>% 
  filter(year %in% c("1914", "2019")) %>% 
  pivot_wider(names_from = year, values_from = mu_fl) %>% 
  mutate(
    mean_decline = (`2019` - `1914`),
    mean_rel_decline = (`2019` - `1914`) / `1914`
  )
  
mean_coefs <- coef_dat %>% 
  group_by(coef) %>% 
  summarize(
    mean_est = mean(est)
  )
est_1914 <- mean_coefs$mean_est[1] + (mean_coefs$mean_est[6] * 1914) 
est_2019 <- mean_coefs$mean_est[1] + (mean_coefs$mean_est[6] * 2019) 
(est_1914 - est_2019) / est_1914


# example time series
comb_dat <- rbind(
  dum %>% select(year, age, sex, fl = sim_fl) %>% mutate(dat = "true"),
  samp_dat %>% select(year, age, sex, fl = sim_fl) %>% mutate(dat = "recovered")
)

ggplot() +
  geom_boxplot(data = comb_dat %>%
                 filter(sex == "female",
                        year %in% seq(1920, 2000, by = 20)),
               aes(x = as.factor(year), y = fl, fill = dat)) +
  facet_wrap(~age)

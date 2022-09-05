# test fecundity effects assuming model fit to cm (rather than mm) data


library(tidyverse)
library(dplyr)
library(forcats)
library(sdmTMB)


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
    fl_cm = fl / 10
  ) %>% 
  droplevels()
levels(dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")


# make fake mesh (necessary in current branch)
dat$x <- runif(nrow(dat))
dat$y <- runif(nrow(dat))
dum_mesh <- make_mesh(dat, c("x", "y"), cutoff = 1000)

fit_cm <- sdmTMB(fl_cm ~ s(yday_c, by = age, m = 2) +
                s(year, by = age, m = 2) +
                age + sex,
              dispformula = ~ 0 + period,
              data = dat,
              mesh = dum_mesh,
              spatial = "off",
              spatiotemporal = "off",
              control = sdmTMBcontrol(
                nlminb_loops = 2,
                newton_loops = 2
              ))
sanity(fit_cm)


fec_dat <- expand.grid(
  age = unique(dat$age),
  sex = "female",
  period = "NFWD",
  yday_c = 0,
  year = seq(min(dat$year), max(dat$year), by = 1),
  # dummy spatial variables required 
  x = runif(1),
  y = runif(1)
) %>% 
  mutate(
    age_f = as.factor(age)
  )
levels(fec_dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")

# components for calculating smooths
source(here::here("R", "utils.R"))
formula_no_sm <- remove_s_and_t2(fit_cm$formula[[1]])
X_ij <- model.matrix(formula_no_sm, data = dat)
sm <- parse_smoothers(fit_cm$formula[[1]], data = dat)
sm_start <- sm$b_smooth_start
pred_X_ij_fec <- predict(mgcv::gam(formula_no_sm, data = dat), 
                         fec_dat, type = "lpmatrix")
sm_pred <- parse_smoothers(fit_cm$formula[[1]], data = dat, 
                           newdata = fec_dat)
pred_X1_ij_fec = pred_X_ij_fec
pred_Zs_fec = sm_pred$Zs
pred_Xs_fec = sm_pred$Xs

# cov matrix and parameter estimates concatenated to same order as matrix rownames
cov <- fit_cm$sd_report$cov.fixed
coef_list <- as.list(fit_cm$sd_report, "Estimate")
coef_vec <- c(coef_list$b_j, coef_list$b_disp_k, coef_list$bs,
              coef_list$ln_smooth_sigma)


n_sims <- 1000
cov_sim <- MASS::mvrnorm(n_sims, coef_vec, cov)
b1_j <- cov_sim[ , grep("b_j", colnames(cov))] #fixed effect pars
bs <- cov_sim[ , grep("bs", colnames(cov))] # fixed effect smooths
ssdr <- summary(fit_cm$sd_report)
b_smooth <- ssdr[rownames(ssdr) == "b_smooth", "Estimate"] %>% as.numeric() #random smooths

# simulate as above but assuming sigma in both years equal to NFWD value and
# include calculations; uses same coefficient samples as prediction intervals
# above
sim_list_fec <- vector(mode = "list", length = n_sims)

for (i in seq_len(n_sims)) {
  pred_mu1 <- pred_X1_ij_fec %*% b1_j[i, ]
  
  pred_smooth <- matrix(NA, nrow(pred_X1_ij_fec), length(sm_start))
  
  for (ss in seq_along(sm_start)) {
    beta_ss <- rep(NA, times = ncol(pred_Zs_fec[[ss]]))
    for (j in seq_along(beta_ss)) {
      beta_ss[j] <- b_smooth[sm_start[ss] + j]
    }
    pred_smooth[ , ss] <- pred_Zs_fec[[ss]] %*% beta_ss
  }
  
  pred_smooth1_i <- apply(pred_smooth, 1, sum)
  pred_smooth2_i <- pred_Xs_fec %*% bs[i, ]
  pred_smooth_i <- pred_smooth1_i + pred_smooth2_i
  pred_mu_i <- pred_mu1 + pred_smooth_i
  
  sim_list_fec[[i]] <- fec_dat %>% 
    mutate(
      iter = i,
      pred_mu = as.numeric(pred_mu_i)
    ) %>% 
    mutate(
      # convert to POH length in mm using inverse of equation in main manuscript
      sim_hl =  (0.833 * pred_mu * 10) - 3.508,
      fec_beta = ifelse(age %in% c("42", "53"), 10.67, 9.77),
      fec_int = ifelse(age %in% c("42", "53"), 1811.9, 1244.7),
      sim_fec = (fec_beta * sim_hl - fec_int) / 1000
    ) %>% 
    group_by(age) %>%
    mutate(
      mean_sim_fec = mean(sim_fec)
    ) %>% 
    filter(
      year %in% c(min(dat$year), max(dat$year))
    ) %>% 
    select(year, age, age_f, iter, sim_fec, mean_sim_fec) %>% 
    pivot_wider(names_from = "year", values_from = "sim_fec") %>% 
    mutate(diff = `2019` - `1914`,
           first_diff = diff / `1914`,
           rel_diff = diff / mean_sim_fec) %>% 
    ungroup()
}

sim_dat_fec <- data.table::rbindlist(sim_list_fec) %>% 
  as.data.frame()


diff_fecundity <- sim_dat_fec %>% 
  group_by(age_f) %>% 
  summarize(
    mean = mean(diff),
    low = quantile(diff, probs = 0.025),
    up = quantile(diff, probs = 0.975),
    rel_mean = mean(first_diff),
    rel_low = quantile(first_diff, probs = 0.025),
    rel_up = quantile(first_diff, probs = 0.975),
    .groups = "drop"
  ) 


png(here::here("outputs", "figs", "fecundity_diff_cm.png"), 
    height = 4, width = 5, units = "in", res = 250)
ggplot(diff_fecundity) +
  geom_pointrange(aes(x = age_f, y = rel_mean, ymin = rel_low, 
                      ymax = rel_up)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  ggsidekick::theme_sleek() +
  scale_x_discrete(
    "Age",
    labels = c(expression("4"["2"]), expression("5"["2"]), 
               expression("5"["3"]), expression("6"["3"]))
  ) +
  scale_y_continuous(
    "Percent Decline in Fecundity (2019-1914)",
    breaks = c(-0.3, -0.2, -0.1, 0.0),
    labels = c("-30%", "-20%", "-10%", "0%")
  ) +
  theme(legend.position = "none")
dev.off()


## changes in size
new_dat2 <- expand.grid(
  age = unique(dat$age),
  sex = "female",
  yday_c = 0,
  year = seq(min(dat$year), max(dat$year), length = 100),
  # dummy spatial variables required 
  x = runif(1),
  y = runif(1)
) %>% 
  mutate(
    age_f = as.factor(age)
  )
levels(new_dat2$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")

# smoothed predictions
smooth_preds <- predict(fit_cm, newdata = new_dat2, re_form = NA, se_fit = TRUE)
smooth_preds$low <- smooth_preds$est + (qnorm(0.025) * smooth_preds$est_se)
smooth_preds$up <- smooth_preds$est + (qnorm(0.975) * smooth_preds$est_se)



smooth_list <- split(smooth_preds, smooth_preds$age_f) 
# smooth_period_list <- split(sim_dat, sim_dat$age_f) 
age_names <- c("42", "52", "53", "63")

cm_diff <- purrr::map2(smooth_list, age_names, function (x, y) {
  ts <- x$est
  #calc time series mean for each iteration
  ts_mean = mean(ts)
  data.frame(
    mean_diff = ts[length(ts)] - ts_mean,
    first_diff = ts[length(ts)] - ts[1],
    mean_rel_diff = (ts[length(ts)] - ts_mean) / ts_mean,
    first_rel_diff = (ts[length(ts)] - ts[1]) / ts[1],
    age = y
  )
}) %>% 
  bind_rows


## check sigma estimates
period_sig <- data.frame(
  sig_est = as.list(fit_cm$sd_report, "Estimate")$b_disp_k,
  sig_se = as.list(fit_cm$sd_report, "Std. Error")$b_disp_k,
  period = unique(dat$period)
) %>% 
  mutate(
    low = sig_est + (qnorm(0.025) * sig_se),
    up = sig_est + (qnorm(0.975) * sig_se)
  )

ggplot(period_sig) +
  geom_pointrange(aes(x = period, y = sig_est, shape = period, 
                      ymin = low, ymax = up),
                  fill = "white") +
  scale_shape_manual(values = shape_pal) +
  ggsidekick::theme_sleek() +
  labs(x = "Sampling Period", y = "Estimated Residual SD") +
  theme(legend.position = "none")

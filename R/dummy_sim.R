## simulation exploring correlated categorical and continous predictors


library(tidyverse)


# set model parameters
set.seed(999)
b0 <- 55 # intercept
b1 <- -0.075 # continuous slope
b2 <- 1.4 # factor level 1 coefficient
b3 <- 2 # factor level 1 coefficient
sigma <- 3.0 # residual standard deviation
N <- 10000 # number of data points


# refit 1000 times to estimate potential bias in parameter estimates
# n_sims <- 100
# dat_out <- matrix(NA, nrow = n_sims, ncol = 4)
# 
# for (i in 1:n_sims) {
  dum_dat <- data.frame(
    x1 = runif(N, -50, 50) %>% round(., digits = 0)
  ) %>% 
    mutate(
      # x2 <- rbinom(N, size = 1, prob = 0.4) # binary predictor data
      x2 = ifelse(x1 > 0 & x1 < 25, 1, 0),
      x3 = ifelse(x1 >= 25, 1, 0),
      y = rnorm(N, mean = b0 + b1*x1 + b2*x2 + b3*x3, sd = sigma),
      x_f = case_when(
        x2 == 1 ~ "2",
        x3 == 1 ~ "3",
        TRUE ~ "1"
      )
    ) 
  
  m <- lm(y ~ x1 + x_f, data = dum_dat)
#   dat_out[i, ] <- coef(m)
# }

m2 <- lm(y ~ x_f, data = dum_dat)


apply(dat_out, 2, mean)
apply(dat_out, 2, sd)


# plot of one draw
ggplot(dum_dat, aes(x = x1, y = y, colour = as.factor(x_f))) +
  geom_point()
dum_dat %>% group_by(x_f) %>% summarize(mean(y))


# 
# 
# m2 <- lm(y ~ x1 + x2, data = dum_dat)
# 
# new_dat <- data.frame(
#   x1 = seq(0, 20, by = 0.1),
#   x2f = "1",
#   x2 = 0.5
# )
# new_dat2 <- new_dat
# 
# m_pred <- predict(m, new_dat, se.fit = TRUE)
# new_dat$m_pred <- m_pred$fit
# new_dat$m_up <- m_pred$fit + (qnorm(0.975) * m_pred$se.fit)
# new_dat$m_low <- m_pred$fit + (qnorm(0.025) * m_pred$se.fit)
# m2_pred <- predict(m2, new_dat, se.fit = TRUE)
# new_dat2$m_pred <- m2_pred$fit
# new_dat2$m_up <- m2_pred$fit + (qnorm(0.975) * m2_pred$se.fit)
# new_dat2$m_low <- m2_pred$fit + (qnorm(0.025) * m2_pred$se.fit)
# 
# 
# dum <- rbind(new_dat %>% mutate(model = "uncentered"), 
#              new_dat2 %>% mutate(model = "centered"))
# 
# 
# ggplot() +
#   geom_line(data = dum, aes(x = x1, y = m_pred, colour = model)) +
#   geom_ribbon(data = dum, aes(x = x1, y = m_pred, ymin = m_low, ymax = m_up,
#                               fill = model), alpha = 0.3) +
#   geom_point()
# 
# set.seed(999)
# b0 <- 1.4 # intercept
# b1 <- 0.2 # continuous slope
# b2 <- 1.7 # factor level 1 coefficient
# b1.2 <- 0.5 # interaction between b1 and b2
# sigma <- 2.0 # residual standard deviation
# N <- 50 # number of data points
# 
# x1 <- runif(N, 0, 20) # continuous predictor data
# x2 <- rbinom(N, size = 1, prob = 0.4) # binary predictor data
# 
# # generate response data:
# y <- rnorm(N, mean = b0 +
#              b1 * x1 +
#              b2 * x2 +
#              x1 * x2 * b1.2,
#            sd = sigma)
# ggplot(dat, aes(x = x1, y = y, colour = as.factor(x2))) +
#   geom_point()
# 
# 
# dat <- data.frame(x1, x2, y)
# head(dat)
# 
# m <- lm(y ~ x1 * x2, data = dat)

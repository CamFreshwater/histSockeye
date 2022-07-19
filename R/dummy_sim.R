set.seed(999)
b0 <- 1.4 # intercept
b1 <- 0.4 # continuous slope
b2 <- 1.7 # factor level 1 coefficient
# b1.2 <- 0.5 # interaction between b1 and b2
sigma <- 2.0 # residual standard deviation
N <- 150 # number of data points

x1 <- runif(N, 0, 20) # continuous predictor data
x2 <- rbinom(N, size = 1, prob = 0.4) # binary predictor data

dum_dat <- data.frame(
  x1 = runif(N, 0, 20)
) %>% 
  mutate(
    x2 = ifelse(x1 < 10, 0, 1),
    y = rnorm(N, mean = b0 + b1*x1 + b2*x2, sd = sigma),
    x2f = as.factor(x2)
  ) 

ggplot(dum_dat, aes(x = x1, y = y, colour = as.factor(x2))) +
  geom_point()


m <- lm(y ~ x1 + x2f, data = dum_dat)
m2 <- lm(y ~ x1 + x2, data = dum_dat)

new_dat <- data.frame(
  x1 = seq(0, 20, by = 0.1),
  x2f = "1",
  x2 = 0.5
)
new_dat2 <- new_dat

m_pred <- predict(m, new_dat, se.fit = TRUE)
new_dat$m_pred <- m_pred$fit
new_dat$m_up <- m_pred$fit + (qnorm(0.975) * m_pred$se.fit)
new_dat$m_low <- m_pred$fit + (qnorm(0.025) * m_pred$se.fit)
m2_pred <- predict(m2, new_dat, se.fit = TRUE)
new_dat2$m_pred <- m2_pred$fit
new_dat2$m_up <- m2_pred$fit + (qnorm(0.975) * m2_pred$se.fit)
new_dat2$m_low <- m2_pred$fit + (qnorm(0.025) * m2_pred$se.fit)


dum <- rbind(new_dat %>% mutate(model = "uncentered"), 
             new_dat2 %>% mutate(model = "centered"))


ggplot() +
  geom_line(data = dum, aes(x = x1, y = m_pred, colour = model)) +
  geom_ribbon(data = dum, aes(x = x1, y = m_pred, ymin = m_low, ymax = m_up,
                              fill = model), alpha = 0.3) +
  geom_point()

set.seed(999)
b0 <- 1.4 # intercept
b1 <- 0.2 # continuous slope
b2 <- 1.7 # factor level 1 coefficient
b1.2 <- 0.5 # interaction between b1 and b2
sigma <- 2.0 # residual standard deviation
N <- 50 # number of data points

x1 <- runif(N, 0, 20) # continuous predictor data
x2 <- rbinom(N, size = 1, prob = 0.4) # binary predictor data

# generate response data:
y <- rnorm(N, mean = b0 +
             b1 * x1 +
             b2 * x2 +
             x1 * x2 * b1.2,
           sd = sigma)
ggplot(dat, aes(x = x1, y = y, colour = as.factor(x2))) +
  geom_point()


dat <- data.frame(x1, x2, y)
head(dat)

m <- lm(y ~ x1 * x2, data = dat)

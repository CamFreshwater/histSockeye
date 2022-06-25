# Fit sdmTMB model to century of Nass data 
# Faster alternative to brms_fit while still accounting for variable dispersion
# 1) Import Skip's first century of data and clean as necessary
# 2) Estimate temporal trend using GAM and brms to account for changes in
# variability with sampling regime
# Created by C Freshwater Dec 23, 2021
# Update w/ edited data
# -------------------------------------------------

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
      year > 1993 ~ "Nisga'a",
      TRUE ~ "Bilton"
    ),
    period = fct_relevel(as.factor(period), 
                         "Gilbert-Clemens",
                         "Bilton",
                         "Monkley Dump",
                         "Nisga'a"),
    fl_cm = fl / 10,
    # add random identifier to split datasets for subsetting,
    data_group = sample.int(4, nrow(.), replace = T)
  ) %>% 
  droplevels()
levels(dat$age_f) <- c("4[2]", "5[2]", "5[3]", "6[3]")


# FIT MODEL  -------------------------------------------------------------------

# parallelize based on operating system (should speed up some of the spatial
# processing calculations)
library(parallel)
ncores <- detectCores() - 2
if (Sys.info()['sysname'] == "Windows") {
  library(doParallel)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
} else {
  doMC::registerDoMC(ncores)
}

dat_trim <- dat %>% sample_n(30000)

# make fake mesh
dat$x <- runif(nrow(dat))
dat$y <- runif(nrow(dat))
dum_mesh <- make_mesh(dat, c("x", "y"), cutoff = 1000)

fit <- sdmTMB(fl_cm ~ s(yday_c, m = 2, k = 3) + 
                s(yday_c, by = age, m = 1, k = 3) +
                s(year, m = 2, k = 5) + 
                s(year, by = age, m = 1, k = 5) +
                age + sex + period,
              dispformula = ~ period,
              data = dat_trim,
              mesh = dum_mesh,
              spatial = "off",
              spatiotemporal = "off")
fit2 <- sdmTMB(fl_cm ~  s(yday_c, m = 2, k = 3) + 
                 s(yday_c, by = age, m = 1, k = 3) +
                 s(year, m = 2, k = 5) + 
                 s(year, by = age, m = 1, k = 5) +
                 age + sex + period,
               dispformula = ~ period,
               data = dat,
               mesh = dum_mesh,
               spatial = "off"#, 
               # spatiotemporal = "off"
               )

brm1_old <- readRDS(here::here("outputs", "data", "brms_fits", "ind_ls.rds"))


# fixed effects predictions
new_dat <- expand.grid(
  age = unique(dat_trim$age),
  sex = unique(dat_trim$sex),
  period = unique(dat_trim$period),
  yday_c = 0,
  year = 1969,#mean(dat_trim$year),
  # dummy spatial variables required 
  x = runif(1),
  y = runif(1)
) 

fe_preds <- predict(fit2, newdata = new_dat, re_form = NA, se_fit = TRUE)
fe_preds <- predict(fit2, re_form = NA, se_fit = TRUE)
fe_preds <- predict(fit2, re_form = NA, se_fit = FALSE)



m2 <- sdmTMB(
  density ~ 0 + depth_scaled + depth_scaled2 + as.factor(year),
  data = pcod_2011, time = "year", mesh = pcod_mesh_2011, family = tweedie(link = "log"),
  dispformula = ~ 0 + as.factor(year),
  spatial = "off", spatiotemporal = "off"
)
p <- predict(m2, re_form = NA, se_fit = TRUE)
head(p[,c("est", "est_se")])
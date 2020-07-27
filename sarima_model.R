library(cdcfluview)
library(tsfknn)
library(timeseries)
library(dplyr)
library(ggplot2)
library(caret)
library(tidyverse)
library(forecast)


interpolate_and_clean_missing <- function(y) {
  if(any(is.na(y) | is.infinite(y)) && !is.na(tail(y, 1))) {
    ## drop leading NAs
    if(is.na(y[1]) | is.infinite(y[1])) {
      num_leading_nas <- rle(is.na(y) | is.infinite(y))$lengths[1]
      y <- y[- seq_len(num_leading_nas)]
    }
    
    ## interpolate internal NAs
    while(any(is.na(y) | is.infinite(y))) {
      na_rle <- rle(is.na(y) | is.infinite(y))
      na_run_start_ind <- na_rle$lengths[1] + 1
      na_run_end_ind <- na_run_start_ind + na_rle$lengths[2] - 1
      y[na_run_start_ind:na_run_end_ind] <-
        approx(
          x = c(na_run_start_ind - 1, na_run_end_ind + 1),
          y = y[c(na_run_start_ind - 1, na_run_end_ind + 1)],
          xout = na_run_start_ind:na_run_end_ind,
          method = "linear"
        )$y
    }
  }
  
  return(y)
}

do_difference <- function(y, d = 0, D = 0, frequency = 1) {
  # first differencing
  for(i in seq_len(d)) {
    y <- ts(
      c(NA,
        y[seq(from = 1 + 1, to = length(y))] -
          y[seq(from = 1, to = length(y) - 1)]),
      frequency = frequency)
  }
  
  # seasonal differencing
  if(D > 0 && frequency < 2) {
    stop("It doesn't make sense to do seasonal differencing with a time series frequency of 1.")
  }
  for(i in seq_len(D)) {
    y <- ts(
      c(rep(NA, frequency),
        y[seq(from = frequency + 1, to = length(y))] -
          y[seq(from = 1, to = length(y) - frequency)]),
      frequency = frequency)
  }
  
  return(y)
}

invert_difference <- function(dy, y, d, D, frequency) {
  for(i in seq_len(d)) {
    y_dm1 <- do_difference(y, d = d-i, D = D, frequency = frequency)
    dy_full <- c(y_dm1, dy)
    for(t in seq_len(length(dy))) {
      dy_full[length(y_dm1) + t] <- dy_full[length(y_dm1) + t - 1] + dy_full[length(y_dm1) + t]
    }
    dy <- dy_full[length(y_dm1) + seq_along(dy)]
  }
  
  for(i in seq_len(D)) {
    y_dm1 <- do_difference(y, d = 0, D = D-i, frequency = frequency)
    dy_full <- c(y_dm1, dy)
    for(t in seq_len(length(dy))) {
      dy_full[length(y_dm1) + t] <- dy_full[length(y_dm1) + t - frequency] + dy_full[length(y_dm1) + t]
    }
    dy <- dy_full[length(y_dm1) + seq_along(dy)]
  }
  
  return(ts(dy, frequency = frequency))
}

usflu <- ilinet(region = "national", years = 1997:2018)

national_train_data <- ilinet(region = "national", years = 1997:2019) %>%
  filter(week != 53) %>%
  select(year,week,weighted_ili) %>% 
  mutate(season = ifelse(week %in% (40:52), year, year-1)) %>%
  mutate(log_wili = log(weighted_ili)) %>%
  select(season, week, weighted_ili, log_wili)

national_test_data <- ilinet(region = "national", years = 2019) %>%
  filter(week != 53) %>%
  select(year,week,weighted_ili) %>% 
  mutate(season = ifelse(week %in% (40:52), year, year-1)) %>%
  select(season, week, weighted_ili) 

difference_y <- do_difference(national_train_data$log_wili, d = 0, D = 1, frequency = 52)
interpolated_y <- interpolate_and_clean_missing(national_train_data$log_wili)

model_fit <- auto.arima(interpolated_y,
                        optim.method = "BFGS",
                        d = 0,
                        D = 1,
                        seasonal = TRUE)

nsim <- 100000
h <- 52
trajectory_samples <- matrix(NA_real_, nrow = nsim, ncol = h)

pred <- forecast(model_fit, h = 2)

for (i in 1:nsim) {
  trajectory_samples[i, ] <- simulate(model_fit, nsim = h, bootstrap = TRUE, future = TRUE)
}

for(i in seq_len(nsim)) {
  trajectory_samples[i, ] <-
    invert_difference(
      dy = trajectory_samples[i, ],
      y = national_train_data$log_wili,
      d = 0,
      D = 1,
      frequency = h)
  
  trajectory_samples[i, ] <- exp(trajectory_samples[i,])
}

seasonal_peak <- matrix(NA_real_, nrow = nsim, ncol = 1)
for(i in 1:nrow(trajectory_samples)){
  seasonal_peak[i,] <- max(trajectory_samples[i,])
}

#Sample season peak probability
dnorm(6.21, mean(seasonal_peak[,1], sd(seasonal_peak[,1])))



library(cdcfluview)
library(dplyr)
library(prophet)

##All supporting functions
interpolate_and_clean_missing <- function(y,dates) {
  
  if(any(is.na(y) | is.infinite(y)) && !is.na(tail(y, 1))) {
    ## drop leading NAs
    if(is.na(y[1]) | is.infinite(y[1])) {
      num_leading_nas <- rle(is.na(y) | is.infinite(y))$lengths[1]
      y <- y[- seq_len(num_leading_nas)]
      dates <- dates[-seq_len(num_leading_nas)]
      ds<-dates
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
  
  return(data.frame(y,ds))
}

fit_prophet<- function(
  data,
  test_season
) {
  require(forecast)
  browser()
  test_season_ind <- min(which(data$season == test_season))
  data <- data[seq_len(test_season_ind - 1), , drop = FALSE]
  #view(data)
  
  undifferenced_target_data <- log(data$weighted_ili)
  
  pred_target <- ts(c(rep(NA, 52),
                      undifferenced_target_data[seq(from = 53, to = nrow(data))] -
                        undifferenced_target_data[seq(from = 1, to = nrow(data) - 52)]),
                    frequency = 52)
  temp <- interpolate_and_clean_missing(pred_target,data$week_start)
  
  model = prophet(temp,yearly.seasonality = TRUE,weekly.seasonality = FALSE,daily.seasonality = FALSE ,mcmc.samples = 0,interval.width = 0.95)
  
  future=make_future_dataframe(model,periods = 52,freq = 'weeks',include_history = FALSE)
  
  
  #return(sarima_fit)
  
}

fit_sarima<- function(
  data,
  test_season
) {
  require(forecast)
  
  test_season_ind <- min(which(data$season == test_season))
  data <- data[seq_len(test_season_ind - 1), , drop = FALSE]
  #view(data)
  
  undifferenced_target_data <- log(data$weighted_ili)
  
  pred_target <- ts(c(rep(NA, 52),
                      undifferenced_target_data[seq(from = 53, to = nrow(data))] -
                        undifferenced_target_data[seq(from = 1, to = nrow(data) - 52)]),
                    frequency = 52)
  pred_target <- interpolate_and_clean_missing(pred_target,data$week_start)$y
  
  #print(pred_target)
  sarima_fit <- auto.arima(pred_target,
                           d = 0,
                           D = 0,
                           stationary = TRUE)
  
  return(sarima_fit)
  
}

sample_predictive_trajectories<- function (
  sarima_fit,
  h,
  nsim = 10^5
) {
  browser()
  sim <- matrix(NA, nrow = nsim, ncol = h)
  
  for (i in 1:nsim) {
    try(sim[i, ] <- simulate(sarima_fit, nsim = h))
  }
  
  return(sim)
  
}


create_sample_trajectories <- function (
  n_sim,
  predict_horizon,
  data,
  analysis_time_season,
  analysis_time_season_week,
  sarima_fit
) {
  #browser()
  data <- data %>%
    filter(season < analysis_time_season | 
             (season == analysis_time_season & 
                season_week <= 
                analysis_time_season_week))
  
  
  
  undifferenced_new_target_data <- log(data$weighted_ili)
  new_target_data <- ts(c(rep(NA, 52),
                          undifferenced_new_target_data[seq(from = 53, to = nrow(data))] -
                            undifferenced_new_target_data[seq(from = 1, to = nrow(data) - 52)]),
                        frequency = 52)
  
  new_target_data <- interpolate_and_clean_missing(new_target_data,data$week_start)
  #print(new_target_data)
  
  model = prophet(new_target_data,yearly.seasonality = TRUE,weekly.seasonality = FALSE,daily.seasonality = FALSE ,mcmc.samples = 0,interval.width = 0.95)
  
  future=make_future_dataframe(model,periods = predict_horizon,freq = 'weeks',include_history = FALSE)
  trajectory_samples <- matrix(NA_real_,nrow = 10000,ncol=predict_horizon)
  k=1
  #Creating  samples
  print(predict_horizon)
  for(i in 1:200){
    
    t=prophet::predictive_samples(model,future)
    for(j in 1:50){
      trajectory_samples[k,]<- exp(t$yhat[1:predict_horizon,j])
      #errors_samples[k,]<-(exp(t$yhat[1:52,j]) - usflu[1093:1144,]$weighted_ili)
      k=k+1
    }
    
  }
  
  #view(raw_trajectory_samples)
  raw_trajectory_samples<-trajectory_samples
  for(prediction_horizon in seq_len(predict_horizon)) {
    trajectory_samples[, prediction_horizon] <- 
      raw_trajectory_samples[, prediction_horizon] +
      undifferenced_new_target_data[nrow(data) + prediction_horizon - 52]
  }
  
  trajectory_samples <- round(exp(trajectory_samples),1)
  #View(trajectory_samples)
  return(trajectory_samples)
}

get_observed_seasonal_quantities <- function(
  data,
  season,
  onset_baseline
) {
  first_season_ind <- min(which(data$season == season))
  last_season_ind <- max(which(data$season == season))
  
  obs_inc_in_season_leading_trailing_nas <-
    data[seq(from = first_season_ind, to = last_season_ind),
         'weighted_ili']
  
  ## pad so that we start at season week 1
  if(data$season_week[first_season_ind] != 1) {
    obs_inc_in_season_leading_trailing_nas <- c(
      rep(NA, data$season_week[first_season_ind] - 1),
      obs_inc_in_season_leading_trailing_nas)
  }
  
  observed_peak_inc <- max(
    obs_inc_in_season_leading_trailing_nas,
    na.rm = TRUE)
  
  observed_peak_inc_bin <- get_inc_bin(observed_peak_inc, return_character = TRUE)
  
  ## peak week timing is based on rounded values
  round_to_.1 <- function(inc_val) {
    if(is.na(inc_val)) {
      return(inc_val)
    } else {
      floor_val <- floor(inc_val * 10) / 10
      if(inc_val >= floor_val + 0.05) {
        return(floor_val + 0.1)
      } else {
        return(floor_val)
      }
    }
  }
  
  rounded_observed_peak_inc <- round(observed_peak_inc, 1)
  rounded_obs_inc_in_season <- unlist(round(obs_inc_in_season_leading_trailing_nas,1))
  
  
  observed_peak_week <-
    which(rounded_obs_inc_in_season == as.numeric(rounded_observed_peak_inc))
  
  
  #print(length(rounded_obs_inc_in_season))
  observed_onset_week <- get_onset_week(
    incidence_trajectory = rounded_obs_inc_in_season,
    baseline = onset_baseline,
    onset_length = 3
  )
  
  return(list(observed_onset_week = observed_onset_week,
              observed_peak_week = observed_peak_week,
              observed_peak_inc = observed_peak_inc,
              observed_peak_inc_bin = observed_peak_inc_bin
  ))
}

#function to compute onset week based on a trajectory of incidence values
get_onset_week <- function(incidence_trajectory,
                           baseline,
                           onset_length = 3) {
  
  exceeded_threshold <- sapply(
    seq_len(length(incidence_trajectory) - onset_length),
    function(start_ind) {
      above_baseline <- incidence_trajectory[seq(from = start_ind, length = onset_length)] >= baseline
      length(above_baseline)>0 &&
        all(above_baseline) &&
        !all(is.na(incidence_trajectory))
    }
  )
  
  if(any(exceeded_threshold, na.rm = TRUE)) {
    season_week <- min(which(exceeded_threshold))
    
    return(season_week)
  } else {
    return("none")
  }
}

get_inc_bin <- function(inc,
                        return_character = TRUE) {
  #print(inc)
  inc <- round(inc, 1)
  bin_numeric <- ifelse(inc < 13,
                        floor(inc*10)/10, ## floors to 1st decimal place
                        13)
  if(return_character) {
    return(as.character(bin_numeric))
  } else {
    return(bin_numeric)
  }
}

logspace_sum <- function(logx) {
  dim(logx) <- c(1, length(logx))
  return(logspace_sum_matrix_rows(logx))
}

make_predictions_dataframe <- function(data,
                                       model_name,
                                       incidence_bin_names = as.character(seq(1,13,0.1)),
                                       first_analysis_time_season_week = 1,
                                       last_analysis_time_season_week = 51) {
  
  ## allocate more than enough space up front,
  ## delete extra later
  na_vec <- rep(NA,
                (last_analysis_time_season_week - first_analysis_time_season_week + 1)
  )
  
  onset_week_bins <- c(as.character(1:33), "none")
  peak_week_bins <- as.character(1:33)
  
  predictions_df <- cbind(
    data.frame(
      analysis_time_season = na_vec,
      analysis_time_season_week = na_vec,
      prediction_week_ph_1 = na_vec,
      prediction_week_ph_2 = na_vec,
      prediction_week_ph_3 = na_vec,
      prediction_week_ph_4 = na_vec,
      peak_week_log_score = na_vec,
      onset_log_score = na_vec,
      peak_inc_log_score = na_vec,
      ph_1_inc_log_score = na_vec,
      ph_2_inc_log_score = na_vec,
      ph_3_inc_log_score = na_vec,
      ph_4_inc_log_score = na_vec,
      stringsAsFactors = FALSE), # don't want model to be a factor with 1 level
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(onset_week_bins))) %>%
      `colnames<-`(paste0("onset_bin_", onset_week_bins, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(peak_week_bins))) %>%
      `colnames<-`(paste0("peak_week_bin_", peak_week_bins, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(incidence_bin_names))) %>%
      `colnames<-`(paste0("peak_inc_bin_", incidence_bin_names, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(incidence_bin_names))) %>%
      `colnames<-`(paste0("ph_1_inc_bin_", incidence_bin_names, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(incidence_bin_names))) %>%
      `colnames<-`(paste0("ph_2_inc_bin_", incidence_bin_names, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(incidence_bin_names))) %>%
      `colnames<-`(paste0("ph_3_inc_bin_", incidence_bin_names, "_log_prob")),
    as.data.frame(matrix(NA, nrow = length(na_vec), ncol = length(incidence_bin_names))) %>%
      `colnames<-`(paste0("ph_4_inc_bin_", incidence_bin_names, "_log_prob"))
  )
  
  return(predictions_df)
}

logspace_sum_matrix_rows <- function(logX) {
  return(log(apply(exp(logX), 1, sum)))
}

Prophet_Forecast <- function(predict_result,nsim){
  mean <- predict_result$yhat 
  sigma <- as.numeric(with(predict_result, ((yhat_upper-yhat_lower)/2)/qnorm(0.975)))
  browser()
  samples <- mapply(rlogitnorm,mu = mean,sigma = sigma,n = nsim) * 100
  
  return(samples)
}

get_log_scores_via_trajectory_simulation <- function(
  data,
  observed_seasonal_quantities,
  analysis_time_season,
  first_analysis_time_season_week = 1, # == week 40 of year
  last_analysis_time_season_week = 4, # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead
  n_trajectory_sims,
  model_object
) {
  #browser()
  incidence_bin_names = as.character(seq(1,13,0.1))
  predictions_df <- make_predictions_dataframe(
    data = data,
    model_name = model_name,
    incidence_bin_names = as.character(seq(1,13,0.1)),
    first_analysis_time_season_week = 1,
    last_analysis_time_season_week = 33)
  
  results_save_row <- 1L
  last_analysis_time_season_week_in_data <- max(data$season_week[data$season == analysis_time_season])
  for(analysis_time_season_week in seq(from = first_analysis_time_season_week, 
                                       to = min(last_analysis_time_season_week, 
                                                last_analysis_time_season_week_in_data - 1))) {
    analysis_time_ind <- which(data$season == analysis_time_season &
                                 data$season_week == analysis_time_season_week)
    max_prediction_horizon <- max(4L,
                                  last_analysis_time_season_week + 1 - analysis_time_season_week)
    
    data1 <- data %>%
      filter(season < analysis_time_season | 
               (season == analysis_time_season & 
                  season_week <= 
                  analysis_time_season_week))
    print(tail(data1$week_start))
    
    undifferenced_new_target_data <- log(data1$weighted_ili)
    new_target_data <- ts(c(rep(NA, 52),
                            undifferenced_new_target_data[seq(from = 53, to = nrow(data1))] -
                              undifferenced_new_target_data[seq(from = 1, to = nrow(data1) - 52)]),
                          frequency = 52)
    
    new_target_data <- interpolate_and_clean_missing(new_target_data,data1$week_start)
    
    
    new_data <- data.frame(ds = new_target_data$ds,
                           y = new_target_data$y)
    
    
    model = prophet(new_data,yearly.seasonality = TRUE,weekly.seasonality = FALSE,daily.seasonality = FALSE ,mcmc.samples = 0,interval.width = 0.95)
    
    trajectory_samples <- matrix(NA_real_,nrow = 10000,ncol=max_prediction_horizon)
    future=make_future_dataframe(model,periods = max_prediction_horizon,freq = 'weeks',include_history = FALSE)
    k=1
    
    for(i in 1:200){
      
      t=prophet::predictive_samples(model,future)
      for(j in 1:50){
        trajectory_samples[k,]<- exp(t$yhat[1:max_prediction_horizon,j])
        #errors_samples[k,]<-(exp(t$yhat[1:52,j]) - usflu[1093:1144,]$weighted_ili)
        k=k+1
      }
      
    }
    
   # browser()
    
    trajectory_samples <- round(exp(trajectory_samples),1)
    
    trajectory_samples <- apply(
      trajectory_samples,
      c(1, 2),
      function(inc_val) {
        if(is.na(inc_val)) {
          return(inc_val)
        } else {
          floor_val <- floor(inc_val * 10) / 10
          if(inc_val >= floor_val + 0.05) {
            return(floor_val + 0.1)
          } else {
            return(floor_val)
          }
        }
      }
    )
    #print(trajectory_samples[1])
    sample_inds_with_na <- apply(
      trajectory_samples[,
                         seq_len(last_analysis_time_season_week + 1 - analysis_time_season_week),
                         drop = FALSE],
      1,
      function(x) any(is.na(x)))
    
    ## Predictions for things about the whole season
    if(!all(sample_inds_with_na)) {
      made_predictions <- TRUE
      
      ## subset to sampled trajectories that are usable/do not have NAs
      subset_trajectory_samples <- trajectory_samples[!sample_inds_with_na, ]
      #print(subset_trajectory_samples[1])
      ## Augment trajectory samples with previous observed incidence values
      ## This is where we should be adding in something to account for backfill.
      first_season_obs_ind <- min(which(data$season == analysis_time_season))
      subset_trajectory_samples <- cbind(
        matrix(
          rep(unlist(data[seq(from = first_season_obs_ind, to = analysis_time_ind), "weighted_ili"]), 
              each = n_trajectory_sims),
          nrow = nrow(subset_trajectory_samples)
        ),
        subset_trajectory_samples
      )
      #print(subset_trajectory_samples[1])
      
      ## If first observation for the season was not at season week 1,
      ## augment with leading NAs
      first_season_obs_week <- data$season_week[first_season_obs_ind]
      #print(first_season_obs_week)
      if(first_season_obs_week != 1) {
        subset_trajectory_samples <- cbind(
          matrix(NA, nrow = nrow(subset_trajectory_samples), ncol = first_season_obs_week - 1),
          subset_trajectory_samples
        )
      }
      
      ## values before the first analysis time week are NA so that
      ## onset and peak calculations only look at data within the CDC's definition
      ## of the flu season for purposes of the competition
      #print(subset_trajectory_samples[1])
      ## Convert to binned values
      binned_trajectory_samples <- 
        get_inc_bin(subset_trajectory_samples,
                    return_character = FALSE)
      #On set week
      onset_week_by_sim_ind <-
        apply(binned_trajectory_samples, 1, function(trajectory) {
          get_onset_week(
            incidence_trajectory = trajectory,
            baseline = 2.3,
            onset_length = 3L
          )
        })
      
      ## Get peak incidence for each simulated trajectory
      peak_inc_bin_by_sim_ind <-
        apply(binned_trajectory_samples, 1, function(trajectory) {
          max(trajectory, na.rm = TRUE)
        })
      
      ## get peak week by sim ind
      ## note that some sim inds may have more than 1 peak week...
      peak_weeks_by_sim_ind <- unlist(lapply(
        seq_len(nrow(binned_trajectory_samples)),
        function(sim_ind) {
          bin_val <- peak_inc_bin_by_sim_ind[sim_ind]
          peak_season_weeks <- which(
            binned_trajectory_samples[sim_ind, ] == bin_val)
          return(peak_season_weeks)
        }
      ))
      
      onset_week_bins <- c(as.character(1:33), "none")
      
      onset_bin_log_probs <- log(sapply(
        onset_week_bins,
        function(bin_name) {
          sum(onset_week_by_sim_ind == bin_name)
        })) 
      onset_bin_log_probs <- onset_bin_log_probs - logspace_sum(onset_bin_log_probs)
      predictions_df[results_save_row, paste0("onset_bin_", onset_week_bins, "_log_prob")] <-
        onset_bin_log_probs
      predictions_df[results_save_row, "onset_log_score"] <-
        onset_bin_log_probs[
          as.character(observed_seasonal_quantities$observed_onset_week)]
      peak_week_bins <- as.character(1:33)
      peak_week_bin_log_probs <- log(sapply(
        peak_week_bins,
        function(bin_name) {
          sum(peak_weeks_by_sim_ind == bin_name)
        }))
      peak_week_bin_log_probs <- peak_week_bin_log_probs - logspace_sum(peak_week_bin_log_probs)
      names(peak_week_bin_log_probs) <- as.character(peak_week_bins)
      predictions_df[results_save_row, paste0("peak_week_bin_", peak_week_bins, "_log_prob")] <-
        peak_week_bin_log_probs
      predictions_df[results_save_row, "peak_week_log_score"] <-
        logspace_sum(peak_week_bin_log_probs[
          as.character(observed_seasonal_quantities$observed_peak_week)])
      
      peak_inc_bin_log_probs <- log(sapply(
        incidence_bin_names,
        function(bin_name) {
          sum(peak_inc_bin_by_sim_ind == as.numeric(bin_name))
        }))
      peak_inc_bin_log_probs <- peak_inc_bin_log_probs - logspace_sum(peak_inc_bin_log_probs)
      predictions_df[results_save_row, paste0("peak_inc_bin_", incidence_bin_names, "_log_prob")] <-
        peak_inc_bin_log_probs
      predictions_df[results_save_row, "peak_inc_log_score"] <-
        peak_inc_bin_log_probs[
          as.character(observed_seasonal_quantities$observed_peak_inc_bin)]
      
    }
    for(ph in 1:4) {
      sample_inds_with_na <- is.na(trajectory_samples[, ph])
      
      ## get observed value/bin
      observed_ph_inc <-
        data[analysis_time_ind + ph, 'weighted_ili']
      observed_ph_inc_bin <- get_inc_bin(observed_ph_inc, return_character = TRUE)
      
      if(!all(sample_inds_with_na) && !is.na(observed_ph_inc)) {
        made_predictions <- TRUE
        
        ## get sampled incidence values at prediction horizon that are usable/not NAs
        ph_inc_by_sim_ind <- trajectory_samples[!sample_inds_with_na, ph]
        ph_inc_bin_by_sim_ind <- get_inc_bin(ph_inc_by_sim_ind, return_character = TRUE)
        
        ## get log score
        ph_inc_bin_log_probs <- log(sapply(
          incidence_bin_names,
          function(bin_name) {
            sum(ph_inc_bin_by_sim_ind == bin_name)
          })) -
          log(length(ph_inc_bin_by_sim_ind))
        predictions_df[results_save_row, paste0("ph_", ph, "_inc_bin_", incidence_bin_names, "_log_prob")] <-
          ph_inc_bin_log_probs
        predictions_df[results_save_row, paste0("ph_", ph, "_inc_log_score")] <-
          ph_inc_bin_log_probs[observed_ph_inc_bin]
      }
    } # ph loop
    
    predictions_df[results_save_row, "analysis_time_season"] <- analysis_time_season
    predictions_df[results_save_row, "analysis_time_season_week"] <- analysis_time_season_week
    predictions_df[results_save_row, "prediction_week_ph_1"] <- analysis_time_season_week + 1
    predictions_df[results_save_row, "prediction_week_ph_2"] <- analysis_time_season_week + 2
    predictions_df[results_save_row, "prediction_week_ph_3"] <- analysis_time_season_week + 3
    predictions_df[results_save_row, "prediction_week_ph_4"] <- analysis_time_season_week + 4
    
    results_save_row <- results_save_row + 1
    
  }
  
  ## if there are extra rows in the predictions_df, delete them
  if(results_save_row <= nrow(predictions_df)) {
    predictions_df <- predictions_df[
      -seq(from = results_save_row, to = nrow(predictions_df)),
      ,
      drop = FALSE
      ]
  }
  return(predictions_df)
  
}

###

usflu <- ilinet(region = "national", years = 1997:2018) %>%
  filter(week != 53) %>%
  select(year,week,weighted_ili,week_start) %>% 
  mutate(season = ifelse(week %in% (40:52), year, year-1)) %>%
  mutate(season = ifelse(week %in% (40:53), year, year-1),
         season_week = if_else(week > 39, abs(week - 39), week + 13),
         weighted_ili = round(weighted_ili,1)) %>%
  select(year, week, season, season_week, weighted_ili,week_start)

test_data <- ilinet(region = "national", years = 2018) %>%
  filter(week != 53) %>%
  select(year,week,weighted_ili,week_start) %>% 
  mutate(season = ifelse(week %in% (40:52), year, year-1)) %>%
  mutate(season = ifelse(week %in% (40:53), year, year-1),
         season_week = if_else(week > 39, abs(week - 39), week + 13),
         weighted_ili = round(weighted_ili,1)) %>%
  select(year, week, season, season_week, weighted_ili,week_start) %>%
  filter(season_week <=51)


observed_seasonal_targets <- get_observed_seasonal_quantities(usflu, 2018,2.2)

sarima_fit <- fit_prophet(usflu, 2018) #no need of this for prophet



pred <- get_log_scores_via_trajectory_simulation(usflu,
                                                 observed_seasonal_targets,
                                                 2018,
                                                 1,
                                                 33,
                                                 1,
                                                 sarima_fit)


library(ggplot2)

plot(pred$onset_log_score, type = 'l')
plot(pred$peak_week_log_score, type = 'l')
plot(pred$peak_inc_log_score, type = 'l')
plot(pred$ph_1_inc_log_score, type = 'l')
plot(pred$ph_2_inc_log_score, type = 'l')
plot(pred$ph_3_inc_log_score, type = 'l')
plot(pred$ph_4_inc_log_score, type = 'l')

week_bins <- c(as.character(1:33), "none")
onset_prob <- cbind(analysis_time_season_week = pred$analysis_time_season_week,
                    exp(pred[, paste0("onset_bin_", week_bins, "_log_prob")]))
library(tidyr)

onset_prob <- onset_prob %>%
  gather("onset_week", "probability", 2:35)

onset_prob$onset_week <- ifelse(onset_prob$onset_week == 'onset_bin_none_log_prob',
                                'none', 
                                as.numeric(gsub("[^[:digit:].]", "",  onset_prob$onset_week)))

week_order <- as.character(c(1:33,'none')) 
onset_prob %>%
  ggplot(mapping = aes(x = factor(onset_week, level = week_order), y = probability)) + 
  geom_line(aes(color = as.integer(analysis_time_season_week), group = as.integer(analysis_time_season_week))) +
  geom_vline(xintercept = observed_seasonal_targets$observed_onset_week, color = "Red") +
  labs(title = "Onset Week Probability Distribution") +
  scale_x_discrete(name = "Onset Week",breaks=c("1","10","20",'30','none')) + 
  labs(y = "Predictive Probability", color = "Analysis Time\nSeason Week") + 
  scale_color_gradient(high = "#132B43", low = "#56B1F7")

#Peak Week
peak_week_bins <- c(as.character(1:33))
peak_prob <- cbind(analysis_time_season_week = pred$analysis_time_season_week,
                   exp(pred[, paste0("peak_week_bin_", peak_week_bins, "_log_prob")]))

peak_prob <- peak_prob %>%
  gather("peak_week", "probability", 2:34)

peak_prob$peak_week <- as.numeric(gsub("[^[:digit:].]", "",  peak_prob$peak_week))

peak_prob %>%
  ggplot(mapping = aes(x = factor(peak_week, level = week_order), y = probability)) + 
  geom_line(aes(color = as.integer(analysis_time_season_week), group = as.integer(analysis_time_season_week))) +
  geom_vline(xintercept = observed_seasonal_targets$observed_peak_week, color = "Red") +
  labs(title = "Peak Week Probability Distribution") +
  scale_x_discrete(name = "Peak Week",breaks=c("1","10","20",'30')) + 
  labs(y = "Predictive Probability", color = "Analysis Time\nSeason Week") + 
  scale_color_gradient(high = "#132B43", low = "#56B1F7")

#Peak incidence
peak_inc_bins <- c(as.character(seq(1,13,0.1)))
peak_inc_prob <- cbind(analysis_time_season_week = pred$analysis_time_season_week,
                       exp(pred[, paste0("peak_inc_bin_", peak_inc_bins, "_log_prob")]))

peak_inc_prob <- peak_inc_prob %>%
  gather("peak_inc_bin", "probability", 2:122)

peak_inc_prob$peak_inc_bin <- as.numeric(gsub("[^[:digit:].]", "",  peak_inc_prob$peak_inc_bin))

peak_inc_prob %>%
  ggplot(mapping = aes(x = as.numeric(peak_inc_bin), y = probability)) + 
  geom_line(aes(color = as.integer(analysis_time_season_week), group = as.integer(analysis_time_season_week))) +
  geom_vline(xintercept = as.numeric(observed_seasonal_targets$observed_peak_inc_bin), color = "Red") +
  labs(title = "Peak Incidence Probability Distribution") +
  scale_x_continuous(name = "Peak Incidence", breaks = c(1,5,10,13)) + 
  labs(y = "Predictive Probability", color = "Analysis Time\nSeason Week") + 
  scale_color_gradient(high = "#132B43", low = "#56B1F7")



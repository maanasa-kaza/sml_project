usflu <- ilinet(region = "national", years = 1997:2018) %>%
  filter(week != 53) %>%
  select(year,week,weighted_ili) %>% 
  mutate(season = ifelse(week %in% (40:52), year, year-1)) %>%
  mutate(season = ifelse(week %in% (40:53), year, year-1),
         season_week = if_else(week > 39, abs(week - 39), week + 13),
         weighted_ili = round(weighted_ili,1)) %>%
  select(year, week, season, season_week, weighted_ili)

test_data <- ilinet(region = "national", years = 2018) %>%
  filter(week != 53) %>%
  select(year,week,weighted_ili) %>% 
  mutate(season = ifelse(week %in% (40:52), year, year-1)) %>%
  mutate(season = ifelse(week %in% (40:53), year, year-1),
         season_week = if_else(week > 39, abs(week - 39), week + 13),
         weighted_ili = round(weighted_ili,1)) %>%
  select(year, week, season, season_week, weighted_ili) %>%
  filter(season_week <=51)

observed_seasonal_targets <- get_observed_seasonal_quantities(usflu, 2018,2.2)
sarima_fit <- fit_sarima(usflu, 2018)
pred <- get_log_scores_via_trajectory_simulation(usflu,
                                                 observed_seasonal_targets,
                                                 2018,
                                                 1,
                                                 33,
                                                 10000,
                                             sarima_fit)
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


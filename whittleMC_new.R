# install packages; pacman will install & load any other packages not in your system
install.packages("pacman")
pacman::p_load(tidyverse, fracdiff, forecast, devtools, RCurl)

# source the code script from Github
script <- getURL("https://raw.githubusercontent.com/taylorgrant/TimeSeries/master/whittleFML.R", ssl.verifypeer = FALSE)
# load the functions into memory
eval(parse(text = script))
rm(script)

# function to calculate RMSE 
rmse <- function(error) {
  sqrt(mean(error^2))
  }

## -------- ## 
# function to run the Monte Carlo simulations and pull out coefficients # 
# this uses the fracdiff.sim from package(fracdiff) to simulate long 
# memory series 

est_fn <- function(d, length, ar = NULL, ma = NULL) {
  
  # if no AR or MA parameters # 
  if (is.null(ar) & is.null(ma)) {
    if (d < .5) {
      tmp_x = fracdiff.sim(length, d = d)$series
    } else {
        tmp_x = cumsum(fracdiff.sim(length,d=(d-1))$series)
    }
  
  tmp_est <- suppressWarnings(whittleFML(tmp_x, inits = list(d = 0)))$coefficients[1,1] %>%
    tibble(d_est = .,
           d_act = d, 
           error = d_est - d_act,
           length = length)
  
  # if there is an AR process that is being estimated as well #
  } else if (!is.null(ar) & is.null(ma)) {
    if (d < .5) {
      tmp_x = fracdiff.sim(length, d = d, ar = ar)$series
    } else {
      tmp_x = cumsum(fracdiff.sim(length,d=(d-1), ar = ar)$series)
    }
    
    tmp_est <- suppressWarnings(whittleFML(tmp_x, inits = list(d = 0, AR = 0)))$coefficients[1:2,1] %>% 
      tibble(d_est = .) %>%
      mutate(est = c("d_est", "ar_est")) %>%
      spread(est, d_est) %>%
      mutate(d_act = d,
             ar_act = ar,
             d_error = d_est - d_act,
             ar_error = ar_est - ar_act,
             length = length)
    
    # if there is an MA prcoess being estimated as well #
  } else if (is.null(ar) & !is.null(ma)) {
    if (d < .5) {
      tmp_x = fracdiff.sim(length, d = d, ma = ma)$series
    } else {
      tmp_x = cumsum(fracdiff.sim(length,d=(d-1), ma = ma)$series)
    }
    
    tmp_est <- suppressWarnings(whittleFML(tmp_x, inits = list(d = 0, MA = 0)))$coefficients[1:2,1] %>% 
      tibble(d_est = .) %>%
      mutate(est = c("d_est", "ma_est")) %>%
      spread(est, d_est) %>%
      mutate(d_act = d,
             ma_act = ma,
             d_error = d_est - d_act,
             ma_error = ma_est - ma_act,
             length = length)
  }
}

## --------- MONTE CARLO --------- ##

# monte carlo simulation of straight d  # 
# 1000 simulations of each # 
set.seed(2367)
d_vals <- c(0,0.25,0.45,0.75)
t_length <- c(40, 60, 80, 100)
ar_params <- c(-0.40, 0, 0.40, 0.80)
ma_params <- c(-0.40, 0, 0.40, 0.80)
mc_length <- 1000

# initialize an empty data frame to capture the MC results
df <- NULL
for (d in d_vals) {
  for (t in t_length) {
      for (len in seq_len(mc_length)) {
      tmp <- est_fn(d = d, length = t, ar = NULL, ma = NULL)
      df <- rbind(df, tmp)
      }
    }
}

# calculate the RMSE and the Bias # 
df %>% group_by(d_act, length) %>%
  summarise(bias = mean(d_est),
            rmse = rmse(error), sims = n()) %>%
  mutate(bias = bias - d_act)


# graph the densities of d estimates # 
facet_df <- df %>%
  group_by(d_act, length) %>%
  summarise(d = mean(d_act))

ggplot(df, aes(x = d_est)) + 
  geom_density() + 
  facet_wrap(vars(d_act, length)) + 
  geom_vline(data = facet_df, aes(group = d_act, xintercept=d), color ="red", lwd=0.5, lty="dashed") + 
  theme_minimal() + 
  labs(x = "d estimate", y = "",
       title = "density of d estimate - (0,d,0)",
       caption = "dashed line represents true value of d")

# monte carlo simulation of d with AR parameter  # 
# 1000 simulations of each # 
# initialize an empty data frame to capture the MC results
df_withAR <- NULL
for (d in d_vals) {
  for (t in t_length) {
    for (ars in ar_params) {
    for (len in seq_len(mc_length)) {
      tmp <- est_fn(d = d, length = t, ar = ars, ma = NULL)
      df_withAR <- rbind(df_withAR, tmp)
    }
  }
  }
}

# calculate the RMSE and the Bias of d and AR params # 
df_withAR %>% group_by(d_act, ar_act, length) %>%
  summarise(d_bias = mean(d_est),
    d_rmse = rmse(d_error), 
    ar_bias = mean(ar_est), 
    ar_rmse = rmse(ar_error),
    sims = n()) %>%
  mutate(d_bias = d_bias - d_act,
         ar_bias = ar_bias - ar_act)

# monte carlo simulation of d with MA parameter  # 
# 1000 simulations of each # 
# initialize an empty data frame to capture the MC results
df_withMA <- NULL
for (d in d_vals) {
  for (t in t_length) {
    for (mas in ma_params) {
      for (len in seq_len(mc_length)) {
        tmp <- est_fn(d = d, length = t, ar = NULL, ma = mas)
        df_withMA <- rbind(df_withMA, tmp)
      }
    }
  }
}

# calculate the RMSE and the Bias of d and AR params# 
df_withMA %>% group_by(d_act, ma_act, length) %>%
  summarise(d_bias = mean(d_est),
            d_rmse = rmse(d_error), 
            ma_bias = mean(ma_est), 
            ma_rmse = rmse(ma_error),
            sims = n()) %>%
  mutate(d_bias = d_bias - d_act,
         ma_bias = ma_bias - ma_act)



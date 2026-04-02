# Spatio-Temporal Analysis and Prediction of PM2.5 Concentrations
# in Inner London — ARIMA and STARIMA Modelling Pipeline
# London Air Quality Network (LAQN), 2024-2026

# LOAD LIBRARIES

library(sf)
library(tidyverse)
library(forecast)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyr)
library(geosphere)
library(tseries)
library(corrplot)
library(spdep)
library(starma)


# LOAD AND PREPARE DATA

london_pm25 <- read.csv("LaqnData24-26-central.csv")

# Convert values to numeric (empty cells become NA)
london_pm25 <- london_pm25 %>%
  mutate(Value = as.numeric(Value))

# Parse datetime and extract year-month
london_pm25 <- london_pm25 %>%
  mutate(
    datetime   = dmy_hm(ReadingDateTime, tz = "Europe/London"),
    year_month = floor_date(datetime, "month")
  )

# Compute daily mean PM2.5 per site
london_pm25$Date <- as.Date(london_pm25$ReadingDateTime, format = "%d/%m/%Y %H:%M")

daily_data <- london_pm25 %>%
  group_by(Site, Date, Latitude, Longitude) %>%
  summarise(pm25 = mean(Value, na.rm = TRUE), .groups = "drop")

# Compute monthly mean PM2.5 per site
monthly_pm25 <- london_pm25 %>%
  filter(!is.na(Value)) %>%
  group_by(Site, year_month) %>%
  summarise(
    pm25_mean = mean(Value),
    n_obs     = n(),
    .groups   = "drop"
  ) %>%
  mutate(
    month = month(year_month, label = TRUE),
    year  = year(year_month)
  )


# EXPLORATORY DATA ANALYSIS

# Time Series Plots

# Daily PM2.5 by site
daily_data %>%
  ggplot(aes(x = Date, y = pm25)) +
  geom_line() +
  facet_wrap(~ Site, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Daily PM2.5 by Site",
    y     = expression("PM"[2.5]~(mu*"g/m"^3))
  )

# Monthly average PM2.5 by site
ggplot(monthly_pm25, aes(x = year_month, y = pm25_mean)) +
  geom_line() +
  facet_wrap(~ Site) +
  labs(
    title = "Monthly Average PM2.5 by Monitoring Site",
    x     = "Month",
    y     = expression("PM"[2.5]~(mu*"g/m"^3))
  ) +
  theme_minimal()

# Seasonal Variation 

# Boxplot by month (all sites combined)
ggplot(monthly_pm25, aes(x = month, y = pm25_mean)) +
  geom_boxplot() +
  labs(
    title = "Seasonal Variation in PM2.5 (All Sites)",
    y     = expression("PM"[2.5]~(mu*"g/m"^3))
  ) +
  theme_minimal()

# Summary Statistics

monthly_pm25 %>%
  group_by(Site) %>%
  summarise(
    mean_pm25 = mean(pm25_mean, na.rm = TRUE),
    sd_pm25   = sd(pm25_mean,   na.rm = TRUE),
    min_pm25  = min(pm25_mean,  na.rm = TRUE),
    max_pm25  = max(pm25_mean,  na.rm = TRUE)
  )

# Spatial Structure

# Pearson correlation matrix
wide_data  <- monthly_pm25 %>%
  select(Site, year_month, pm25_mean) %>%
  pivot_wider(names_from = Site, values_from = pm25_mean)

cor_matrix <- cor(wide_data[, -1], use = "complete.obs")
corrplot(cor_matrix, method = "color",
         title = "Correlation Matrix of Monthly PM2.5 by Site")

# Variance by site
monthly_pm25 %>%
  group_by(Site) %>%
  summarise(variance = var(pm25_mean, na.rm = TRUE)) %>%
  arrange(desc(variance))

# Distance Matrix 

coords <- data.frame(
  Site = c("CE2", "WMD", "WM5", "WM6", "LB4", "CW3", "HP1", "CD1", "GR9"),
  lat  = c(51.5076, 51.4922, 51.5120, 51.5139, 51.4641,
           51.5030, 51.4497, 51.5442, 51.4564),
  lon  = c(-0.1331, -0.1471, -0.1216, -0.1528, -0.1146,
           -0.0180, -0.0374, -0.1753,  0.0407)
)

# Geographic distance matrix (km)
dist_matrix <- distm(coords[, c("lon", "lat")], fun = distHaversine) / 1000
rownames(dist_matrix) <- colnames(dist_matrix) <- coords$Site
round(dist_matrix, 2)

# Moran's I (Spatial Autocorrelation)

# Inverse distance weights
inv_dist       <- 1 / dist_matrix
diag(inv_dist) <- 0
lw             <- mat2listw(inv_dist, style = "W")


# first day with complete data across all 9 sites
complete_days <- daily_data %>%
  group_by(Date) %>%
  summarise(n_sites    = n(),
            n_complete = sum(!is.na(pm25))) %>%
  filter(n_sites == 9, n_complete == 9) %>%
  arrange(Date)

# Extract PM2.5 for first complete day, ordered to match coords
one_day <- daily_data %>%
  filter(Date == complete_days$Date[1]) %>%
  arrange(Site)

# Ensure site matches coords
one_day   <- one_day[match(coords$Site, one_day$Site), ]
pm_values <- one_day$pm25

moran.test(pm_values, lw)


# ARIMA MODELLING - ALL SITES

# Helper function: build daily ts, run ADF, return ts object
build_ts <- function(site_code) {
  df <- daily_data %>%
    filter(Site == site_code) %>%
    arrange(Date) %>%
    filter(!is.na(pm25))
  ts(df$pm25, frequency = 7)
}

# Build time series for all sites
ts_ce2 <- build_ts("CE2")
ts_wmd <- build_ts("WMD")
ts_wm5 <- build_ts("WM5")
ts_wm6 <- build_ts("WM6")
ts_lb4 <- build_ts("LB4")
ts_cw3 <- build_ts("CW3")
ts_hp1 <- build_ts("HP1")
ts_cd1 <- build_ts("CD1")
ts_gr9 <- build_ts("GR9")

site_ts_list <- list(
  CE2 = ts_ce2, WMD = ts_wmd, WM5 = ts_wm5,
  WM6 = ts_wm6, LB4 = ts_lb4, CW3 = ts_cw3,
  HP1 = ts_hp1, CD1 = ts_cd1, GR9 = ts_gr9
)

# ADF Tests (all sites)

for (site in names(site_ts_list)) {
  result <- adf.test(site_ts_list[[site]])
  cat(sprintf("%-4s: ADF statistic = %.4f, p-value = %.4f\n",
              site, result$statistic, result$p.value))
}

# ACF/PACF - Original Series (example: WM6)

par(mfrow = c(1, 2))
acf(ts_wm6,  lag.max = 50, main = "ACF: Original (WM6)")
pacf(ts_wm6, lag.max = 50, main = "PACF: Original (WM6)")

# Differencing: First-Order and Seasonal (all sites)

diff_ts <- function(ts_obj) diff(diff(ts_obj, 1), lag = 7)

ts_ce2_diff <- diff_ts(ts_ce2)
ts_wmd_diff <- diff_ts(ts_wmd)
ts_wm5_diff <- diff_ts(ts_wm5)
ts_wm6_diff <- diff_ts(ts_wm6)
ts_lb4_diff <- diff_ts(ts_lb4)
ts_cw3_diff <- diff_ts(ts_cw3)
ts_hp1_diff <- diff_ts(ts_hp1)
ts_cd1_diff <- diff_ts(ts_cd1)
ts_gr9_diff <- diff_ts(ts_gr9)

# ADF on differenced series (d=1, D=1)(all sites)
diff_list <- list(
  CE2 = ts_ce2_diff, WMD = ts_wmd_diff, WM5 = ts_wm5_diff,
  WM6 = ts_wm6_diff, LB4 = ts_lb4_diff, CW3 = ts_cw3_diff,
  HP1 = ts_hp1_diff, CD1 = ts_cd1_diff, GR9 = ts_gr9_diff
)
for (site in names(diff_list)) {
  result <- adf.test(diff_list[[site]])
  cat(sprintf("%-4s: ADF statistic = %.4f, p-value = %.4f\n",
              site, result$statistic, result$p.value))
}

# ACF/PACF after differencing (example: WM6)
par(mfrow = c(1, 2))
acf(ts_wm6_diff,  lag.max = 50, main = "ACF: d=1, D=1 (WM6)")
pacf(ts_wm6_diff, lag.max = 50, main = "PACF: d=1, D=1 (WM6)")

# ARIMA Model Fitting - All Sites 
# Candidate models evaluated:
#   (A) ARIMA(0,1,1)(1,1,1)[7]  
#   (B) ARIMA(1,1,1)(1,1,1)[7]  

fit_arima <- function(ts_obj) {
  fit_a <- Arima(ts_obj, order = c(0, 1, 1), seasonal = c(1, 1, 1))
  fit_b <- Arima(ts_obj, order = c(1, 1, 1), seasonal = c(1, 1, 1))
  list(model_A = fit_a, model_B = fit_b)
}

arima_fits <- lapply(site_ts_list, fit_arima)

# AIC comparison across all sites
cat("\nAIC Comparison: ARIMA(0,1,1)(1,1,1)[7] vs ARIMA(1,1,1)(1,1,1)[7]\n")
aic_table <- data.frame(
  Site    = names(arima_fits),
  AIC_A   = sapply(arima_fits, function(x) AIC(x$model_A)),
  AIC_B   = sapply(arima_fits, function(x) AIC(x$model_B)),
  Selected = sapply(arima_fits, function(x)
    ifelse(AIC(x$model_A) < AIC(x$model_B), "ARIMA(0,1,1)(1,1,1)[7]",
           "ARIMA(1,1,1)(1,1,1)[7]"))
)
print(aic_table)

# ARIMA Diagnostics - Ljung-Box Test (all sites)

cat("\nLjung-Box Test on ARIMA(1,1,1)(1,1,1)[7] Residuals\n")
for (site in names(arima_fits)) {
  lb <- Box.test(residuals(arima_fits[[site]]$model_B),
                 lag = 20, type = "Ljung-Box")
  cat(sprintf("%-4s: statistic = %.4f, p-value = %.4f\n",
              site, lb$statistic, lb$p.value))
}

# Full residual diagnostic plot for selected site (WM6)
checkresiduals(arima_fits$WM6$model_B)

# ARIMA Forecasting (train/test split 80/20, example: WM6)

n_wm6       <- length(ts_wm6)
train_size  <- floor(0.8 * n_wm6)
train_wm6   <- ts_wm6[1:train_size]
test_wm6    <- ts_wm6[(train_size + 1):n_wm6]

fit_wm6_train <- Arima(ts(train_wm6, frequency = 7),
                       order    = c(1, 1, 1),
                       seasonal = c(1, 1, 1))

forecast_wm6 <- forecast(fit_wm6_train, h = length(test_wm6))

autoplot(forecast_wm6) +
  autolayer(ts(test_wm6, frequency = 7), series = "Observed") +
  theme_minimal() +
  labs(
    title = "ARIMA(1,1,1)(1,1,1)[7] Forecast vs Observed — WM6",
    y     = expression("PM"[2.5]~(mu*"g/m"^3))
  )

cat("\ARIMA Forecast Accuracy (WM6, test set)\n")
print(accuracy(forecast_wm6, test_wm6))


# STARIMA MODELLING

# Construct Space-Time Data Matrix

# Wide format: rows = days, columns = sites
wide_daily <- daily_data %>%
  select(Date, Site, pm25) %>%
  pivot_wider(names_from = Site, values_from = pm25) %>%
  arrange(Date)

wide_daily  <- na.omit(wide_daily)
Z           <- as.matrix(wide_daily[, -1])

# Apply first-order and seasonal differencing (d=1, D=1 at lag 7)
Z_diff      <- apply(Z, 2, function(x) diff(diff(x, 1), lag = 7))
Z_diff_clean <- Z_diff[complete.cases(Z_diff), ]

cat("Space-time matrix dimensions (time x sites):", dim(Z_diff_clean), "\n")
cat("NAs in Z_diff_clean:", sum(is.na(Z_diff_clean)), "\n")

# Build Spatial Weight Matrix 

W_mat      <- listw2mat(lw)
W0         <- diag(ncol(W_mat))         # Spatial order 0 = identity
W1         <- W_mat                      # Spatial order 1 = inverse-distance weights
wlist      <- list(W0, W1)

# STACF and STPACF Identification

stacf(Z_diff_clean, wlist, 48)
stpacf(Z_diff_clean, wlist, 48)

# Fit STARIMA Models 

# Candidate 1: MA(1) with spatial order 1 at tlag 1
ma1       <- matrix(c(1, 1), nrow = 1, ncol = 2)
fit_star1 <- starma(Z_diff_clean, wlist, ar = 0, ma = ma1, iterate = 5)

# Candidate 2: AR(1) + MA(1) with spatial order 1
ar2       <- matrix(c(1, 1), nrow = 1, ncol = 2)
ma2       <- matrix(c(1, 1), nrow = 1, ncol = 2)
fit_star2 <- starma(Z_diff_clean, wlist, ar = ar2, ma = ma2, iterate = 5)

# Candidate 3: MA(1) temporal only
ma3       <- matrix(c(1, 0), nrow = 1, ncol = 2)
fit_star3 <- starma(Z_diff_clean, wlist, ar = 0, ma = ma3, iterate = 5)

# Candidate 4: AR(1) + MA(1) temporal only
ar4       <- matrix(c(1, 0), nrow = 1, ncol = 2)
ma4       <- matrix(c(1, 0), nrow = 1, ncol = 2)
fit_star4 <- starma(Z_diff_clean, wlist, ar = ar4, ma = ma4, iterate = 5)

# Candidate 5: MA(1) + seasonal MA(7), temporal only 
# Mirrors ARIMA(0,1,1)(0,1,1)[7] structure in the spatio-temporal domain
ma5       <- matrix(0, nrow = 7, ncol = 2)
ma5[1, ]  <- c(1, 0)   # tlag 1: temporal MA
ma5[7, ]  <- c(1, 0)   # tlag 7: seasonal temporal MA
fit_star5 <- starma(Z_diff_clean, wlist, ar = 0, ma = ma5, iterate = 5)

# BIC comparison
cat("\nSTARIMA Model BIC Comparison\n")
bic_table <- data.frame(
  Model         = c("fit_star1 (MA1 spatial)",
                    "fit_star2 (AR1+MA1 spatial)",
                    "fit_star3 (MA1 temporal)",
                    "fit_star4 (AR1+MA1 temporal)",
                    "fit_star5 (MA1 + seasonal MA7 temporal) *selected*"),
  BIC = round(c(fit_star1$bic, fit_star2$bic, fit_star3$bic,
                fit_star4$bic, fit_star5$bic), 2)
)
print(bic_table)

# Parameter summary for selected model
summary(fit_star3)

# STARIMA Diagnostics

# Residual STACF
stacf(fit_star3$residuals, wlist, 24)

# Residual Moran's I (single time slice)
resid_slice3 <- fit_star3$residuals[nrow(fit_star3$residuals), ]
print(moran.test(resid_slice3, lw))

resid_slice5 <- fit_star5$residuals[nrow(fit_star5$residuals), ]
print(moran.test(resid_slice5, lw))

# Observed vs fitted plot - site WM6
site_idx <- which(colnames(Z_diff_clean) == "WM6")

fitted_df <- data.frame(
  time     = 1:nrow(Z_diff_clean),
  observed = Z_diff_clean[, site_idx],
  fitted   = Z_diff_clean[, site_idx] - fit_star3$residuals[, site_idx]
)

ggplot(fitted_df, aes(x = time)) +
  geom_line(aes(y = observed), colour = "black", alpha = 0.5) +
  geom_line(aes(y = fitted), colour = "steelblue") +
  labs(
    title = "STARIMA Observed vs Fitted — WM6",
    x     = "Time (days)",
    y     = expression("PM"[2.5]~"(differenced)")
  ) +
  theme_minimal()


# STARIMA Forecasting - fit_star3 (STARMA(0,1))

# one parameter, theta10
# theta is a 1x2 matrix: row = tlag 1, columns = slag 0 and slag 1
theta10 <- fit_star3$theta[1, 1]   # tlag 1, slag 0
cat(sprintf("theta10 (tlag=1, slag=0): %.6f\n", theta10))

resids  <- fit_star3$residuals
n       <- nrow(resids)
n_sites <- ncol(resids)
h       <- 30   # forecast horizon (days)

# For MA(1): only the 1-step ahead forecast uses real residuals
# All subsequent steps revert to zero
forecasts <- matrix(0, nrow = h, ncol = n_sites)

for (t in 1:h) {
  idx1 <- n + t - 1
  ma1_contrib <- if (idx1 >= 1 && idx1 <= n) theta10 * resids[idx1, ] else rep(0, n_sites)
  forecasts[t, ] <- ma1_contrib
}

colnames(forecasts) <- colnames(Z_diff_clean)


# Plot 30-day forecast
forecast_df        <- as.data.frame(forecasts)
forecast_df$horizon <- 1:h

forecast_long <- forecast_df %>%
  pivot_longer(-horizon, names_to = "Site", values_to = "forecast")

ggplot(forecast_long, aes(x = horizon, y = forecast, colour = Site)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(
    title = "STARMA(0,1) 30-Day Ahead Forecasts by Site",
    x     = "Forecast Horizon (days)",
    y     = expression("PM"[2.5]~(mu*"g/m"^3))
  ) +
  theme_minimal()

# Forecast Accuracy

# For MA(1), only the 1-step forecast carries real signal
# Both 1-step and 7-step horizons evaluated

test_Z <- tail(Z_diff_clean, 7)

# 1-step accuracy
mae_star_1  <- mean(abs(forecasts[1, ] - test_Z[1, ]), na.rm = TRUE)
rmse_star_1 <- sqrt(mean((forecasts[1, ] - test_Z[1, ])^2, na.rm = TRUE))
cat(sprintf("STARMA(0,1) 1-step MAE:  %.4f\n", mae_star_1))
cat(sprintf("STARMA(0,1) 1-step RMSE: %.4f\n", rmse_star_1))

# 7-step accuracy (steps 2-7 forecast = 0, measures full window)
forecasts_eval <- forecasts[1:7, ]
mae_star_7  <- mean(abs(forecasts_eval - test_Z),     na.rm = TRUE)
rmse_star_7 <- sqrt(mean((forecasts_eval - test_Z)^2, na.rm = TRUE))
cat(sprintf("STARMA(0,1) 7-step MAE:  %.4f\n", mae_star_7))
cat(sprintf("STARMA(0,1) 7-step RMSE: %.4f\n", rmse_star_7))

# Naive seasonal benchmark 
naive_forecast <- matrix(
  rep(Z_diff_clean[nrow(Z_diff_clean) - 7, ], 7),
  nrow = 7, byrow = TRUE
)
mae_naive  <- mean(abs(naive_forecast - test_Z),     na.rm = TRUE)
rmse_naive <- sqrt(mean((naive_forecast - test_Z)^2, na.rm = TRUE))
cat(sprintf("Naive seasonal 7-step MAE:  %.4f\n", mae_naive))
cat(sprintf("Naive seasonal 7-step RMSE: %.4f\n", rmse_naive))



# Additional Model Testing - fit_star5

# STARIMA(0,1)(0,1)[7] Forecasting

# Extract theta parameters by matrix position
# theta is a 7x2 matrix: rows = tlags, columns = spatial orders
theta10 <- fit_star5$theta[1, 1]   # tlag 1, slag 0
theta70 <- fit_star5$theta[7, 1]   # tlag 7, slag 0

cat(sprintf("\ntheta10 (tlag=1, slag=0): %.6f\n", theta10))
cat(sprintf("theta70 (tlag=7, slag=0): %.6f\n", theta70))

resids   <- fit_star5$residuals
n        <- nrow(resids)
n_sites  <- ncol(resids)
h        <- 30   # forecast days

# set to zero beyond observed range
forecasts <- matrix(0, nrow = h, ncol = n_sites)

for (t in 1:h) {
  idx1 <- n + t - 1
  idx7 <- n + t - 7
  
  ma1_contrib <- if (idx1 >= 1 && idx1 <= n) theta10 * resids[idx1, ] else rep(0, n_sites)
  ma7_contrib <- if (idx7 >= 1 && idx7 <= n) theta70 * resids[idx7, ] else rep(0, n_sites)
  
  forecasts[t, ] <- ma1_contrib + ma7_contrib
}

colnames(forecasts) <- colnames(Z_diff_clean)

# Plot 30-day forecast
forecast_df   <- as.data.frame(forecasts)
forecast_df$horizon <- 1:h

forecast_long <- forecast_df %>%
  pivot_longer(-horizon, names_to = "Site", values_to = "forecast")

ggplot(forecast_long, aes(x = horizon, y = forecast, colour = Site)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(
    title = "STARIMA(0,1)(0,1)[7] 30-Day Ahead Forecasts by Site",
    x     = "Forecast Horizon (days)",
    y     = expression("PM"[2.5]~(mu*"g/m"^3))
  ) +
  theme_minimal()

# Forecast Accuracy (7-step horizon where MA signal is active)

test_Z         <- tail(Z_diff_clean, 7)
forecasts_eval <- forecasts[1:7, ]

mae_star  <- mean(abs(forecasts_eval - test_Z),      na.rm = TRUE)
rmse_star <- sqrt(mean((forecasts_eval - test_Z)^2,  na.rm = TRUE))

# Naive seasonal benchmark (repeat last observed week)
naive_forecast <- matrix(
  rep(Z_diff_clean[nrow(Z_diff_clean) - 7, ], 7),
  nrow = 7, byrow = TRUE
)
mae_naive  <- mean(abs(naive_forecast - test_Z),     na.rm = TRUE)
rmse_naive <- sqrt(mean((naive_forecast - test_Z)^2, na.rm = TRUE))

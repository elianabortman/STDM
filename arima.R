#load required libraries
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

#navigate to the csv file containing the data
london_pm25 <- read.csv("LaqnData24-26-central.csv")

#values are numeric and empty cells become NA
london_pm25 <- london_pm25 %>%
  mutate(Value = as.numeric(Value))

#date and time values cleaned
london_pm25 <- london_pm25 %>%
  mutate(datetime = dmy_hm(ReadingDateTime, tz = "Europe/London"),
    year_month = floor_date(datetime, "month"))

#compute daily means pm2.5 value for each site
london_pm25$Date <- as.Date(london_pm25$ReadingDateTime)
daily_data <- london_pm25 %>%
  group_by(Site, Date, Latitude, Longitude) %>%
  summarise(pm25 = mean(Value, na.rm = TRUE),
            .groups = "drop")
#plot daily pm2.5 by site
daily_data %>%
  ggplot(aes(x = Date, y = pm25)) +
  geom_line() +
  facet_wrap(~ Site, scales = "free_y") +
  theme_minimal() +
  labs(title = "Daily PM2.5 by Site")

#compute monthly average pm2.5 value for each site
monthly_pm25 <- london_pm25 %>%
  filter(!is.na(Value)) %>%
  group_by(Site, year_month) %>%
  summarise(pm25_mean = mean(Value),
    n_obs = n(),
    .groups = "drop")

monthly_pm25 <- monthly_pm25 %>%
  mutate(
    month = month(year_month, label = TRUE),
    year  = year(year_month)
  )

#plot monthly average pm2.5 values for each site
ggplot(monthly_pm25, aes(x = year_month, y = pm25_mean)) +
  geom_line() +
  facet_wrap(~ Site) +
  labs(
    x = "Month",
    y = expression("PM"[2.5]~(mu*"g/m"^3)),
    title = "Monthly Average PM2.5 by Monitoring Site"
  ) +
  theme_minimal()


#summary statistics table by site
monthly_pm25 %>%
  group_by(Site) %>%
  summarise(
    mean_pm25 = mean(pm25_mean, na.rm = TRUE),
    sd_pm25   = sd(pm25_mean, na.rm = TRUE),
    min_pm25  = min(pm25_mean, na.rm = TRUE),
    max_pm25  = max(pm25_mean, na.rm = TRUE)
  )

#overall distribution box plots
ggplot(monthly_pm25, aes(x = Site, y = pm25_mean)) +
  geom_boxplot() +
  labs(
    title = "Distribution of Monthly PM2.5 by Site",
    y = expression("PM"[2.5]~(mu*"g/m"^3))
  ) +
  theme_minimal()

#seasonal box plot
ggplot(monthly_pm25, aes(x = month, y = pm25_mean)) +
  geom_boxplot() +
  labs(
    title = "Seasonal Variation in PM2.5 (All Sites)",
    y = expression("PM"[2.5]~(mu*"g/m"^3))
  ) +
  theme_minimal()

#mean seasonal cycle per site
monthly_pm25 %>%
  group_by(Site, month) %>%
  summarise(mean_pm25 = mean(pm25_mean, na.rm = TRUE)) %>%
  ggplot(aes(x = month, y = mean_pm25, group = Site, colour = Site)) +
  geom_line() +
  labs(
    title = "Average Seasonal Cycle by Site",
    y = expression("PM"[2.5]~(mu*"g/m"^3))
  ) +
  theme_minimal()


#correlation matrix
wide_data <- monthly_pm25 %>%
  select(Site, year_month, pm25_mean) %>%
  pivot_wider(names_from = Site, values_from = pm25_mean)

cor_matrix <- cor(wide_data[,-1], use = "complete.obs")

cor_matrix

corrplot(cor_matrix, method = "color")


#variance comparison across sites
monthly_pm25 %>%
  group_by(Site) %>%
  summarise(variance = var(pm25_mean, na.rm = TRUE)) %>%
  arrange(desc(variance))


#spatial gradient: distance matrix
coords <- data.frame(
  Site = c("CE2","WMD","WM5","WM6","LB4","CW3","HP1","CD1","GR9"),
  lat = c(
    51.5075819,
    51.4922482,
    51.5119770,
    51.5139287,
    51.4641135,
    51.5029879,
    51.4496740,
    51.54421,
    51.456357
  ),
  lon = c(
    -0.1330684,
    -0.1471148,
    -0.1216272,
    -0.1527927,
    -0.1145810,
    -0.0180220,
    -0.0374180,
    -0.175269,
    0.040725
  )
)

dist_matrix <- as.matrix(dist(coords[, c("lat","lon")]))
round(dist_matrix, 3)

#distance matrix using geographic distances in km
dist_matrix <- distm(
  coords[, c("lon","lat")],
  fun = distHaversine
)

dist_matrix <- dist_matrix / 1000  #convert to km
round(dist_matrix, 2)

#moran's I for daily regional mean
# extract PM2.5 values for one day
one_day <- daily_data %>%
     filter(Date == min(Date)) %>%
     arrange(Site)
pm_values <- one_day$pm25
moran.test(pm_values, lw)




#daily data, missing days removed, weekly cycle
#Model Identification: ADF test and ACF behavior, ACF/PACF plots
ce2_daily <- daily_data %>%
  filter(Site == "CE2") %>%
  arrange(Date) %>%
  filter(!is.na(pm25))
ts_ce2 <- ts(ce2_daily$pm25, frequency = 7)
adf.test(ts_ce2)
Acf(ts_ce2, lag.max = 400)


wm5_daily <- daily_data %>%
  filter(Site == "WM5") %>%
  arrange(Date) %>%
  filter(!is.na(pm25))
ts_wm5 <- ts(wm5_daily$pm25, frequency = 7)
adf.test(ts_wm5)
Acf(ts_wm5, lag.max = 400)


wm6_daily <- daily_data %>%
  filter(Site == "WM6") %>%
  arrange(Date) %>%
  filter(!is.na(pm25))
ts_wm6 <- ts(wm6_daily$pm25, frequency = 7)
adf.test(ts_wm6)
Acf(ts_wm6, lag.max = 400)


wmd_daily <- daily_data %>%
  filter(Site == "WMD") %>%
  arrange(Date) %>%
  filter(!is.na(pm25))
ts_wmd <- ts(wmd_daily$pm25, frequency = 7)
adf.test(ts_wmd)
Acf(ts_wmd, lag.max = 400)


cw3_daily <- daily_data %>%
  filter(Site == "CW3") %>%
  arrange(Date) %>%
  filter(!is.na(pm25))
ts_cw3 <- ts(cw3_daily$pm25, frequency = 7)
adf.test(ts_cw3)
Acf(ts_cw3, lag.max = 400)


cd1_daily <- daily_data %>%
  filter(Site == "CD1") %>%
  arrange(Date) %>%
  filter(!is.na(pm25))
ts_cd1 <- ts(cd1_daily$pm25, frequency = 7)
adf.test(ts_cd1)
Acf(ts_cd1, lag.max = 400)


lb4_daily <- daily_data %>%
  filter(Site == "LB4") %>%
  arrange(Date) %>%
  filter(!is.na(pm25))
ts_lb4 <- ts(lb4_daily$pm25, frequency = 7)
adf.test(ts_lb4)
Acf(ts_lb4, lag.max = 400)


hp1_daily <- daily_data %>%
  filter(Site == "HP1") %>%
  arrange(Date) %>%
  filter(!is.na(pm25))
ts_hp1 <- ts(hp1_daily$pm25, frequency = 7)
adf.test(ts_hp1)
Acf(ts_hp1, lag.max = 400)


gr9_daily <- daily_data %>%
  filter(Site == "GR9") %>%
  arrange(Date) %>%
  filter(!is.na(pm25))
ts_gr9 <- ts(gr9_daily$pm25, frequency = 7)
adf.test(ts_gr9)
Acf(ts_gr9, lag.max = 400)

#ACF and PACF 
par(mfrow = c(1, 2))
acf(ts_ce2, lag.max = 50, main = "ACF: Original")
pacf(ts_ce2, lag.max = 50, main = "PACF: Original")

par(mfrow = c(1, 2))
acf(ts_wm5, lag.max = 50, main = "ACF: Original")
pacf(ts_wm5, lag.max = 50, main = "PACF: Original")

par(mfrow = c(1, 2))
acf(ts_wm6, lag.max = 50, main = "ACF: Original")
pacf(ts_wm6, lag.max = 50, main = "PACF: Original")

par(mfrow = c(1, 2))
acf(ts_wmd, lag.max = 50, main = "ACF: Original")
pacf(ts_wmd, lag.max = 50, main = "PACF: Original")

par(mfrow = c(1, 2))
acf(ts_cd1, lag.max = 50, main = "ACF: Original")
pacf(ts_cd1, lag.max = 50, main = "PACF: Original")

par(mfrow = c(1, 2))
acf(ts_cw3, lag.max = 50, main = "ACF: Original")
pacf(ts_cw3, lag.max = 50, main = "PACF: Original")

par(mfrow = c(1, 2))
acf(ts_lb4, lag.max = 50, main = "ACF: Original")
pacf(ts_lb4, lag.max = 50, main = "PACF: Original")

par(mfrow = c(1, 2))
acf(ts_hp1, lag.max = 50, main = "ACF: Original")
pacf(ts_hp1, lag.max = 50, main = "PACF: Original")

par(mfrow = c(1, 2))
acf(ts_gr9, lag.max = 50, main = "ACF: Original")
pacf(ts_gr9, lag.max = 50, main = "PACF: Original")


#first difference WM6
ts_wm6_diff1 <- diff(ts_wm6, differences = 1)
adf.test(ts_wm6_diff1)
#seasonal difference WM6
ts_wm6_diff2 <- diff(diff(ts_wm6, 1), lag = 7)
adf.test(ts_wm6_diff2)

acf(ts_wm6_diff2, lag.max = 50, main = "ACF: d=1, D=1")
pacf(ts_wm6_diff2, lag.max = 50, main = "PACF: d=1, D=1")

#first difference CE2
ts_ce2_diff1 <- diff(ts_ce2, differences = 1)
adf.test(ts_ce2_diff1)
#seasonal difference CE2
ts_ce2_diff2 <- diff(diff(ts_ce2, 1), lag = 7)
adf.test(ts_ce2_diff2)

acf(ts_ce2_diff2, lag.max = 50, main = "ACF: d=1, D=1")
pacf(ts_ce2_diff2, lag.max = 50, main = "PACF: d=1, D=1")

#first difference WM5
ts_wm5_diff1 <- diff(ts_wm5, differences = 1)
adf.test(ts_wm5_diff1)
#seasonal difference WM5
ts_wm5_diff2 <- diff(diff(ts_wm5, 1), lag = 7)
adf.test(ts_wm5_diff2)

acf(ts_wm5_diff2, lag.max = 50, main = "ACF: d=1, D=1")
pacf(ts_wm5_diff2, lag.max = 50, main = "PACF: d=1, D=1")

#first difference WMD
ts_wmd_diff1 <- diff(ts_wmd, differences = 1)
adf.test(ts_wmd_diff1)
#seasonal difference WMD
ts_wmd_diff2 <- diff(diff(ts_wmd, 1), lag = 7)
adf.test(ts_wmd_diff2)

acf(ts_wmd_diff2, lag.max = 50, main = "ACF: d=1, D=1")
pacf(ts_wmd_diff2, lag.max = 50, main = "PACF: d=1, D=1")

#first difference HP1
ts_hp1_diff1 <- diff(ts_hp1, differences = 1)
adf.test(ts_hp1_diff1)
#seasonal difference HP1
ts_hp1_diff2 <- diff(diff(ts_hp1, 1), lag = 7)
adf.test(ts_hp1_diff2)

acf(ts_hp1_diff2, lag.max = 50, main = "ACF: d=1, D=1")
pacf(ts_hp1_diff2, lag.max = 50, main = "PACF: d=1, D=1")

#first difference LB4
ts_lb4_diff1 <- diff(ts_lb4, differences = 1)
adf.test(ts_lb4_diff1)
#seasonal difference LB4
ts_lb4_diff2 <- diff(diff(ts_lb4, 1), lag = 7)
adf.test(ts_lb4_diff2)

acf(ts_lb4_diff2, lag.max = 50, main = "ACF: d=1, D=1")
pacf(ts_lb4_diff2, lag.max = 50, main = "PACF: d=1, D=1")

#first difference GR9
ts_gr9_diff1 <- diff(ts_gr9, differences = 1)
adf.test(ts_gr9_diff1)
#seasonal difference GR9
ts_gr9_diff2 <- diff(diff(ts_gr9, 1), lag = 7)
adf.test(ts_gr9_diff2)

acf(ts_gr9_diff2, lag.max = 50, main = "ACF: d=1, D=1")
pacf(ts_gr9_diff2, lag.max = 50, main = "PACF: d=1, D=1")

#first difference CD1
ts_cd1_diff1 <- diff(ts_cd1, differences = 1)
adf.test(ts_cd1_diff1)
#seasonal difference CD1
ts_cd1_diff2 <- diff(diff(ts_cd1, 1), lag = 7)
adf.test(ts_cd1_diff2)

acf(ts_cd1_diff2, lag.max = 50, main = "ACF: d=1, D=1")
pacf(ts_cd1_diff2, lag.max = 50, main = "PACF: d=1, D=1")

#first difference CW3
ts_cw3_diff1 <- diff(ts_cw3, differences = 1)
adf.test(ts_cw3_diff1)
#seasonal difference CW3
ts_cw3_diff2 <- diff(diff(ts_cw3, 1), lag = 7)
adf.test(ts_cw3_diff2)

acf(ts_cw3_diff2, lag.max = 50, main = "ACF: d=1, D=1")
pacf(ts_cw3_diff2, lag.max = 50, main = "PACF: d=1, D=1")




# WM6 ARIMA modelling pipeline
# Model 1: simple AR
fit1 <- Arima(ts_wm6, order = c(1,0,0))

# Model 2: AR + seasonal AR
fit2 <- Arima(ts_wm6, order = c(1,0,0),
              seasonal = c(1,0,0))

# Model 3: ARMA
fit3 <- Arima(ts_wm6, order = c(1,0,1),
              seasonal = c(1,0,0))

# Model 4: include seasonal MA
fit4 <- Arima(ts_wm6, order = c(1,0,1),
              seasonal = c(1,0,1))

fit_auto <- auto.arima(ts_wm6,
                       seasonal = TRUE,
                       stepwise = FALSE,
                       approximation = FALSE)

fit_seasonal <- Arima(ts_wm6,
                      order = c(1,0,1),
                      seasonal = c(1,1,1))

fit <- Arima(ts_wm6,
             order = c(0,1,1),
             seasonal = c(1,1,1))

fits <- Arima(ts_wm6,
              order = c(1,1,1),
              seasonal = c(1,1,1))

AIC(fit1, fit2, fit3, fit4)
AIC(fit, fits)

#Ljung-Box test
checkresiduals(fit4)
checkresiduals(fit_auto)
checkresiduals(fit_seasonal)
checkresiduals(fit)
checkresiduals(fits)

# CW3 ARIMA fitting
cw3fit_auto <- auto.arima(ts_cw3,
                       seasonal = TRUE,
                       stepwise = FALSE,
                       approximation = FALSE)

cw3fit_seasonal <- Arima(ts_cw3,
                      order = c(1,0,1),
                      seasonal = c(1,1,1))

cw3fit <- Arima(ts_cw3,
             order = c(0,1,1),
             seasonal = c(1,1,1))

cw3fits <- Arima(ts_cw3,
              order = c(1,1,1),
              seasonal = c(1,1,1))

AIC(cw3fit4)
AIC(cw3fit, cw3fits)
AIC(cw3fit_auto)
AIC(cw3fit_seasonal)

#Ljung-Box test
checkresiduals(cw3fit_auto)












#train/test split 80/20
n <- length(ts_cw3)
train_size <- floor(0.8 * n)

train <- ts_cw3[1:train_size]
test  <- ts_cw3[(train_size+1):n]

#forecasting
forecast_values <- forecast(cw3fit_auto,
                            h = length(test))
autoplot(forecast_values) +
  autolayer(as.ts(test), series = "Observed") +
  theme_minimal()

#RMSE, MAE, MAPE
accuracy(forecast_values, test)



------------------------------------------------------------------------------
# STARIMA pipeline





W_mat <- listw2mat(lw)
W0 <- diag(ncol(W_mat))
W1 <- W_mat
wlist <- list(W0, W1)

cat("Number of sites:", ncol(Z_diff_clean), "\n")
cat("wlist length:", length(wlist), "\n")
cat("W0 dimensions:", dim(W0), "\n")
cat("W1 dimensions:", dim(W1), "\n")

# Candidate 1: MA(1) with spatial order 1 — no AR terms
# ar matrix: 0 rows needed, so use ar=0 (integer, means no AR)
# ma matrix: 1 row (tlag 1), 2 cols (slag 0 and slag 1), both = 1
ma1 <- matrix(c(1, 1), nrow = 1, ncol = 2)
fit_star1 <- starma(Z_diff_clean, wlist, ar = 0, ma = ma1, iterate = 5)
summary(fit_star1)

# Candidate 2: AR(1) with spatial order 1, MA(1) with spatial order 1
ar2 <- matrix(c(1, 1), nrow = 1, ncol = 2)
ma2 <- matrix(c(1, 1), nrow = 1, ncol = 2)
fit_star2 <- starma(Z_diff_clean, wlist, ar = ar2, ma = ma2, iterate = 5)
summary(fit_star2)

# Candidate 3: MA(1) temporal only, no spatial MA
ma3 <- matrix(c(1, 0), nrow = 1, ncol = 2)
fit_star3 <- starma(Z_diff_clean, wlist, ar = 0, ma = ma3, iterate = 5)
summary(fit_star3)

# Candidate 4: AR(1) temporal only, MA(1) temporal only — no spatial
ar4 <- matrix(c(1, 0), nrow = 1, ncol = 2)
ma4 <- matrix(c(1, 0), nrow = 1, ncol = 2)
fit_star4 <- starma(Z_diff_clean, wlist, ar = ar4, ma = ma4, iterate = 5)
summary(fit_star4)

# Compare BIC
cat("BIC fit_star1 (MA1 spatial):", fit_star1$bic, "\n")
cat("BIC fit_star2 (AR1+MA1 spatial):", fit_star2$bic, "\n")
cat("BIC fit_star3 (MA1 temporal only):", fit_star3$bic, "\n")
cat("BIC fit_star4 (AR1+MA1 temporal only):", fit_star4$bic, "\n")


# Residual STACF — want all lags within confidence bands
stacf(fit_star3$residuals, wlist, 24)

# Residual Moran's I on a single time slice
resid_slice <- fit_star3$residuals[nrow(fit_star1$residuals), ]
moran.test(resid_slice, lw)


# Candidate 5: MA(1) at tlag 1 + MA(1) at tlag 7, temporal only
# Directly mirrors ARIMA(0,1,1)(0,1,1)[7] in STARIMA form
ma5 <- matrix(0, nrow = 7, ncol = 2)
ma5[1, ] <- c(1, 0)   # tlag 1: temporal MA only
ma5[7, ] <- c(1, 0)   # tlag 7: seasonal MA temporal only

fit_star5 <- starma(Z_diff_clean, wlist, ar = 0, ma = ma5, iterate = 5)
summary(fit_star5)
cat("BIC fit_star5:", fit_star5$bic, "\n")

# Check residuals
stacf(fit_star5$residuals, wlist, 24)

resid_slice5 <- fit_star5$residuals[nrow(fit_star5$residuals), ]
moran.test(resid_slice5, lw)

# forecasting
theta10 <- fit_star5$theta[1]
theta70 <- fit_star5$theta[2]

# Get last residuals for MA terms
resids <- fit_star5$residuals
n <- nrow(resids)

# h-step ahead forecast for each site
h <- 30  # forecast horizon in days

# Create extended residuals matrix padded with zeros for future steps
resids_extended <- rbind(resids, matrix(0, nrow = h, ncol = 9))

# Forecast loop
forecasts <- matrix(NA, nrow = h, ncol = 9)

for (t in 1:h) {
  ma1_contrib <- theta10 * resids_extended[n + t - 1, ]
  ma7_contrib <- theta70 * resids_extended[n + t - 7, ]
  forecasts[t, ] <- ma1_contrib + ma7_contrib
}

# Check forecasts look reasonable
head(forecasts)
colnames(forecasts) <- colnames(Z_diff_clean)

# Plot forecasts for all sites
forecast_df <- as.data.frame(forecasts)
forecast_df$horizon <- 1:h

forecast_long <- forecast_df %>%
  pivot_longer(-horizon, names_to = "Site", values_to = "forecast")

ggplot(forecast_long, aes(x = horizon, y = forecast, colour = Site)) +
  geom_line() +
  labs(
    title = "STARIMA 30-Day Ahead Forecasts by Site",
    x = "Forecast Horizon (days)",
    y = expression("PM"[2.5]~(mu*"g/m"^3))
  ) +
  theme_minimal()





# Use last 30 rows of Z_diff_clean as test set
# and everything before as training
n_total <- nrow(Z_diff_clean)
test_Z   <- Z_diff_clean[(n_total - h + 1):n_total, ]
train_Z  <- Z_diff_clean[1:(n_total - h), ]

# Re-fit on training data only
resids_train <- fit_star5$residuals[1:(nrow(fit_star5$residuals) - h), ]
n_train <- nrow(resids_train)

resids_train_ext <- rbind(resids_train, matrix(0, nrow = h, ncol = 9))

forecasts_eval <- matrix(NA, nrow = h, ncol = 9)
for (t in 1:h) {
  ma1_contrib <- theta10 * resids_train_ext[n_train + t - 1, ]
  ma7_contrib <- theta70 * resids_train_ext[n_train + t - 7, ]
  forecasts_eval[t, ] <- ma1_contrib + ma7_contrib
}

# MAE and RMSE vs actual test values
mae_star  <- mean(abs(forecasts_eval - test_Z), na.rm = TRUE)
rmse_star <- sqrt(mean((forecasts_eval - test_Z)^2, na.rm = TRUE))

cat("STARIMA MAE:", round(mae_star, 4), "\n")
cat("STARIMA RMSE:", round(rmse_star, 4), "\n")





# Diagnose first — check what residuals look like
cat("Dimensions of residuals:", dim(resids), "\n")
cat("Any NAs in residuals:", sum(is.na(resids)), "\n")
cat("Last few rows of residuals:\n")
tail(resids)

# Check theta values extracted correctly
cat("theta10:", theta10, "\n")
cat("theta70:", theta70, "\n")


# First inspect the full theta object
print(fit_star5$theta)
print(fit_star5$theta_sd)

# Check the full summary to see parameter names and order
summary(fit_star5)

# Extract by name rather than position
theta10 <- fit_star5$theta["theta10"]
theta70 <- fit_star5$theta["theta70"]

cat("theta10:", theta10, "\n")
cat("theta70:", theta70, "\n")

# If names don't work, print all contents and extract positionally
cat("All theta values:\n")
print(fit_star5$theta)
cat("Length of theta:", length(fit_star5$theta), "\n")
cat("Names of theta:", names(fit_star5$theta), "\n")





# theta is a 7x1 matrix — extract by row position
# tlag 1 = row 1, tlag 7 = row 7, column 1 = slag 0
theta10 <- fit_star5$theta[1, 1]   # tlag 1, slag 0
theta70 <- fit_star5$theta[7, 1]   # tlag 7, slag 0

cat("theta10:", theta10, "\n")
cat("theta70:", theta70, "\n")

# Now rebuild forecasts
resids <- fit_star5$residuals
n <- nrow(resids)
n_sites <- ncol(resids)
h <- 30

forecasts <- matrix(0, nrow = h, ncol = n_sites)

for (t in 1:h) {
  
  # MA(1) at tlag 1
  idx1 <- n + t - 1
  if (idx1 <= n && idx1 >= 1) {
    ma1_contrib <- theta10 * resids[idx1, ]
  } else {
    ma1_contrib <- rep(0, n_sites)
  }
  
  # Seasonal MA at tlag 7
  idx7 <- n + t - 7
  if (idx7 <= n && idx7 >= 1) {
    ma7_contrib <- theta70 * resids[idx7, ]
  } else {
    ma7_contrib <- rep(0, n_sites)
  }
  
  forecasts[t, ] <- ma1_contrib + ma7_contrib
}

colnames(forecasts) <- colnames(Z_diff_clean)

# Verify non-NA
cat("Any NAs in forecasts:", sum(is.na(forecasts)), "\n")
cat("Head of forecasts:\n")
head(forecasts, 10)

# Plot
forecast_df <- as.data.frame(forecasts)
forecast_df$horizon <- 1:h

forecast_long <- forecast_df %>%
  pivot_longer(-horizon, names_to = "Site", values_to = "forecast")

ggplot(forecast_long, aes(x = horizon, y = forecast, colour = Site)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(
    title = "STARIMA 30-Day Ahead Forecasts by Site",
    x = "Forecast Horizon (days)",
    y = expression("PM"[2.5]~(mu*"g/m"^3))
  ) +
  theme_minimal()

# Evaluate accuracy on first 7 steps where MA signal exists
test_Z <- tail(Z_diff_clean, 7)
forecasts_eval <- forecasts[1:7, ]

mae_star  <- mean(abs(forecasts_eval - test_Z), na.rm = TRUE)
rmse_star <- sqrt(mean((forecasts_eval - test_Z)^2, na.rm = TRUE))

cat("STARIMA 7-step MAE:", round(mae_star, 4), "\n")
cat("STARIMA 7-step RMSE:", round(rmse_star, 4), "\n")


# 1. Get your ARIMA accuracy on same differenced scale and same 7-step horizon
#    for direct comparison
# Check what scale your ARIMA accuracy was computed on
# If on original scale, recompute ARIMA on differenced scale for fair comparison

# Get ARIMA residuals on differenced series for fair comparison
# Re-run ARIMA on differenced data
fit_arima_diff <- Arima(ts(Z_diff_clean[,1], frequency = 7),
                        order = c(0,0,1),
                        seasonal = c(0,0,1))

# Or extract accuracy from your existing ARIMA forecast
# and confirm which scale it used
accuracy(forecast_values, test)

# For a fairer longer-horizon comparison, consider also evaluating
# a naive seasonal benchmark
naive_forecast <- matrix(rep(Z_diff_clean[nrow(Z_diff_clean) - 7, ], 7),
                         nrow = 7, byrow = TRUE)
mae_naive  <- mean(abs(naive_forecast - test_Z), na.rm = TRUE)
rmse_naive <- sqrt(mean((naive_forecast - test_Z)^2, na.rm = TRUE))

cat("Naive seasonal MAE:", round(mae_naive, 4), "\n")
cat("Naive seasonal RMSE:", round(rmse_naive, 4), "\n")
# 2. Build final comparison table
results_table <- data.frame(
  Model = c("ARIMA(0,1,1)(1,1,1)[7]", "STARMA(0,1)(0,1)[7]", "Naive Seasonal"),
  BIC = c("—", "20267.32", "—"),
  MAE_7step = c("—", "5.6021", "5.1651"),
  RMSE_7step = c("—", "6.4882", "6.3856"),
  Residual_MoranI_p = c("< 0.001", "0.445", "—")
)
print(results_table)

# 3. Plot observed vs fitted for a single site to show model fit visually
site_idx <- which(colnames(Z_diff_clean) == "WM6")

fitted_df <- data.frame(
  time = 1:nrow(Z_diff_clean),
  observed = Z_diff_clean[, site_idx],
  fitted = Z_diff_clean[, site_idx] - fit_star5$residuals[, site_idx]
)

ggplot(fitted_df, aes(x = time)) +
  geom_line(aes(y = observed), colour = "black", alpha = 0.5) +
  geom_line(aes(y = fitted), colour = "steelblue") +
  labs(
    title = "STARIMA Observed vs Fitted — WM6",
    x = "Time (days)",
    y = expression("PM"[2.5]~"(differenced)")
  ) +
  theme_minimal()


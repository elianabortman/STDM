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

#navigate to the csv file containing the data
london_pm25 <- read.csv("STDM/Coursework/LaqnData24-26-central.csv")

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
# Extract PM2.5 values for one day
one_day <- daily_data %>%
     filter(Date == min(Date)) %>%
     arrange(Site)
pm_values <- one_day$pm25
moran.test(pm_values, lw)

#monthly data stationarity: ADF test per site
ce2_data <- monthly_pm25 %>%
  filter(Site == "CE2") %>%
  arrange(year_month)
wm5_data <- monthly_pm25 %>%
  filter(Site == "WM5") %>%
  arrange(year_month)
wm6_data <- monthly_pm25 %>%
  filter(Site == "WM6") %>%
  arrange(year_month)
wmd_data <- monthly_pm25 %>%
  filter(Site == "WMD") %>%
  arrange(year_month)
cw3_data <- monthly_pm25 %>%
  filter(Site == "CW3") %>%
  arrange(year_month)
cd1_data <- monthly_pm25 %>%
  filter(Site == "CD1") %>%
  arrange(year_month)
lb4_data <- monthly_pm25 %>%
  filter(Site == "LB4") %>%
  arrange(year_month)
hp1_data <- monthly_pm25 %>%
  filter(Site == "HP1") %>%
  arrange(year_month)
gr9_data <- monthly_pm25 %>%
  filter(Site == "GR9") %>%
  arrange(year_month)

ts_ce2 <- ts(ce2_data$pm25_mean,
             frequency = 12)
ts_wm5 <- ts(wm5_data$pm25_mean,
             start = c(2024, 1),
             frequency = 12)
ts_wm6 <- ts(wm6_data$pm25_mean,
             start = c(2024, 1),
             frequency = 12)
ts_wmd <- ts(wmd_data$pm25_mean,
             start = c(2024, 1),
             frequency = 12)
ts_cw3 <- ts(cw3_data$pm25_mean,
             start = c(2024, 1),
             frequency = 12)
ts_cd1 <- ts(cd1_data$pm25_mean,
             start = c(2024, 1),
             frequency = 12)
ts_lb4 <- ts(lb4_data$pm25_mean,
             start = c(2024, 1),
             frequency = 12)
ts_hp1 <- ts(hp1_data$pm25_mean,
             start = c(2024, 1),
             frequency = 12)
ts_gr9 <- ts(gr9_data$pm25_mean,
             start = c(2024, 1),
             frequency = 12)

adf.test(ts_ce2)
adf.test(ts_wm5)
adf.test(ts_wm6)
adf.test(ts_wmd)
adf.test(ts_cw3)
adf.test(ts_cd1)
adf.test(ts_lb4)
adf.test(ts_hp1)
adf.test(ts_gr9)

ts_ce2_diff <- diff(ts_ce2, lag = 12)
ts_wm5_diff <- diff(ts_wm5, lag = 12)
ts_wm6_diff <- diff(ts_wm6, lag = 12)
ts_wmd_diff <- diff(ts_wmd, lag = 12)
ts_cw3_diff <- diff(ts_cw3, lag = 12)
ts_cd1_diff <- diff(ts_cd1, lag = 12)
ts_lb4_diff <- diff(ts_lb4, lag = 12)
ts_hp1_diff <- diff(ts_hp1, lag = 12)
ts_gr9_diff <- diff(ts_gr9, lag = 12)

adf.test(ts_ce2_diff)
adf.test(ts_wm5_diff)
adf.test(ts_wm6_diff)
adf.test(ts_wmd_diff)
adf.test(ts_cw3_diff)
adf.test(ts_cd1_diff)
adf.test(ts_lb4_diff)
adf.test(ts_hp1_diff)
adf.test(ts_gr9_diff)

#first difference d=1 and seasonal difference D=1, lag=12
ts_ce2_diff2 <- diff(ts_ce2, differences = 1)
ts_ce2_diff2 <- diff(ts_ce2_diff2, lag = 12)
adf.test(ts_ce2_diff2)

ts_wm5_diff2 <- diff(ts_wm5, differences = 1)
ts_wm5_diff2 <- diff(ts_wm5_diff2, lag = 12)
adf.test(ts_wm5_diff2)

ts_wm6_diff2 <- diff(ts_wm6, differences = 1)
ts_wm6_diff2 <- diff(ts_wm6_diff2, lag = 12)
adf.test(ts_wm6_diff2)

ts_wmd_diff2 <- diff(ts_wmd, differences = 1)
ts_wmd_diff2 <- diff(ts_wmd_diff2, lag = 12)
adf.test(ts_wmd_diff2)

ts_cw3_diff2 <- diff(ts_cw3, differences = 1)
ts_cw3_diff2 <- diff(ts_cw3_diff2, lag = 12)
adf.test(ts_cw3_diff2)

ts_cd1_diff2 <- diff(ts_cd1, differences = 1)
ts_cd1_diff2 <- diff(ts_cd1_diff2, lag = 12)
adf.test(ts_cd1_diff2)

ts_lb4_diff2 <- diff(ts_lb4, differences = 1)
ts_lb4_diff2 <- diff(ts_lb4_diff2, lag = 12)
adf.test(ts_lb4_diff2)

ts_hp1_diff2 <- diff(ts_hp1, differences = 1)
ts_hp1_diff2 <- diff(ts_hp1_diff2, lag = 12)
adf.test(ts_hp1_diff2)

ts_gr9_diff2 <- diff(ts_gr9, differences = 1)
ts_gr9_diff2 <- diff(ts_gr9_diff2, lag = 12)
adf.test(ts_gr9_diff2)

#ACF and PACF
par(mfrow = c(1, 2))
acf(ts_ce2, main = "ACF: Original Series")
pacf(ts_ce2, main = "PACF: Original Series")
par(mfrow = c(1, 1))

#first difference
ts_ce2_diff1 <- diff(ts_ce2, differences = 1)
par(mfrow = c(1, 2))
acf(ts_ce2_diff1, main = "ACF: First Difference")
pacf(ts_ce2_diff1, main = "PACF: First Difference")
#seasonal difference
ts_ce2_seasdiff <- diff(ts_ce2, lag = 12)
par(mfrow = c(1, 2))
acf(ts_ce2_seasdiff, main = "ACF: Seasonal Difference")
pacf(ts_ce2_seasdiff, main = "PACF: Seasonal Difference")
#first and seasonal differencing
ts_ce2_diff2 <- diff(diff(ts_ce2, differences = 1), lag = 12)
par(mfrow = c(1, 2))
acf(ts_ce2_diff2, main = "ACF: d=1, D=1")
pacf(ts_ce2_diff2, main = "PACF: d=1, D=1")





#daily data, missing days removed, weekly cycle
ce2_daily <- daily_data %>%
  filter(Site == "CE2") %>%
  arrange(Date) %>%
  filter(!is.na(pm25))

ts_ce2 <- ts(ce2_daily$pm25, frequency = 7)

adf.test(ts_ce2)

Acf(ts_ce2, lag.max = 400)

#first difference
ts_ce2_diff1 <- diff(ts_ce2, differences = 1)
adf.test(ts_ce2_diff1)
#seasonal difference
ts_ce2_diff2 <- diff(diff(ts_ce2, 1), lag = 7)
adf.test(ts_ce2_diff2)
#ACF and PACF 
par(mfrow = c(1, 2))
acf(ts_ce2, lag.max = 50, main = "ACF: Original")
pacf(ts_ce2, lag.max = 50, main = "PACF: Original")

par(mfrow = c(1, 2))
acf(ts_ce2_diff2, lag.max = 50, main = "ACF: d=1, D=1")
pacf(ts_ce2_diff2, lag.max = 50, main = "PACF: d=1, D=1")

par(mfrow = c(1,1))




#train/test split 80/20
n <- length(ts_ce2)
train_size <- floor(0.8 * n)

train <- ts_ce2[1:train_size]
test  <- ts_ce2[(train_size+1):n]

#model 1: auto arima
fit_auto <- auto.arima(train,
                       seasonal = TRUE,
                       stepwise = FALSE,
                       approximation = FALSE)

summary(fit_auto)
#model 2
fit_manual <- Arima(train,
                    order = c(0,1,1),
                    seasonal = c(0,1,1))

summary(fit_manual)
#compare AIC
AIC(fit_auto, fit_manual)

checkresiduals(fit_auto)
#forecasting
forecast_values <- forecast(fit_auto,
                            h = length(test))
autoplot(forecast_values) +
  autolayer(test, series = "Observed") +
  theme_minimal()
#RMSE, MAE, MAPE
accuracy(forecast_values, test)







#model 1: auto.arima
auto.arima(ts_ce2,
           seasonal = TRUE,
           stepwise = FALSE,
           approximation = FALSE)

#model 2: Seasonal Arima
Arima(ts_ce2, order=c(0,1,1), seasonal=c(0,1,1))

fit1 <- Arima(ts_ce2, order=c(0,1,1), seasonal=c(0,1,0))
fit2 <- Arima(ts_ce2, order=c(0,1,1), seasonal=c(1,1,0))
fit3 <- Arima(ts_ce2, order=c(1,1,1), seasonal=c(0,1,0))

AIC(fit1, fit3)
checkresiduals(fit3)


-----------------------------------------------------------------------------
#autocorrelation
site_data <- monthly_pm25 %>%
  filter(Site == "WM5") %>%
  arrange(year_month)

ts_pm25 <- ts(
  site_data$pm25_mean,
  start = c(year(site_data$year_month[1]),
            month(site_data$year_month[1])),
  frequency = 12
)
#plot original time series
autoplot(ts_pm25) +
  labs(
    title = "Monthly PM2.5 (Original)",
    y = expression("PM"[2.5]~(mu*"g/m"^3)),
    x = "Year"
  ) +
  theme_minimal()

#seasonal differencing
ts_pm25_seasdiff <- diff(ts_pm25, lag = 12)
autoplot(ts_pm25_seasdiff) +
  labs(
    title = "Seasonally Differenced PM2.5 (lag = 12)",
    y = "Differenced PM2.5",
    x = "Year"
  ) +
  theme_minimal()

#ACF and PACF of original series
par(mfrow = c(1, 2))
acf(ts_pm25, main = "ACF: Original Series")
pacf(ts_pm25, main = "PACF: Original Series")
par(mfrow = c(1, 1))
#ACF and PACF after seasonal differencing
par(mfrow = c(1, 2))
acf(ts_pm25_seasdiff, main = "ACF: Seasonal Difference")
pacf(ts_pm25_seasdiff, main = "PACF: Seasonal Difference")
par(mfrow = c(1, 1))

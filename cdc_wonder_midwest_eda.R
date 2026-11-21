# cdc_wonder_midwest_eda.R
#
# Data from Midwest 2018-2024, by county

library(here)
library(readxl)
library(sf)
library(tigris)
library(ggplot2)
library(tidyverse)
library(spdep)
# pak::pak("connordonegan/geostan")
library(geostan)

data_path <- yaml::read_yaml(here("auth.yaml"))

mid_counties <- counties(
  state = c(
    "MI", "IL", "IN", "IA", "KS", "MN", "MO", "NE", "ND", "OH", "SD", "WI"
  ), 
  cb = T,
  year = 2024
)

mid_states <- states(cb = T, year = 2024) |> 
  filter(NAME %in% c("Michigan", "Illinois", "Indiana", "Iowa", "Kansas",
                     "Minnesota", "Missouri", "Nebrasha", "North Dakota",
                     "Ohio", "South Dakota", "Wisconsin"))

# 2018-2024
cdc_data <- read_xlsx(
  here(data_path, "overdose_deaths_by_county_midwest_2018-2024.xlsx"),
  range = "B1:L1056"
)

# Post-COVID
# Issues with population missing in this data pull
# cdc_data <- read_xlsx(
#   here(data_path, "overdose_deaths_by_county_midwest_03.2020-2024.xlsx"),
#   range = "B1:L1056"
# )

# Pre-COVID
# Issues with population missing in this data pull
# cdc_data <- read_xlsx(
#   here(data_path, "overdose_deaths_by_county_midwest_2015-02.2020.xlsx"),
#   range = "B1:J1056"
# )

# Process data
cdc_data <- cdc_data |> 
  # CDC has counties formatted as "County, ST" - remove ", ST"
  mutate(
    County = County |> str_remove_all(", .."),
    # Convert "Suppressed" to NA
    across(
      .cols = Deaths:`Crude Rate Upper 95% Confidence Interval`, 
      .fns = ~replace_when(.x, .x == "Suppressed" ~ NA)
    ),
    # Change some data columns from chr to num
    across(
      .cols = Deaths:`Crude Rate Upper 95% Confidence Interval`, 
      .fns = as.numeric
    )
  )

# Merge cdc data to county
map_data <- left_join(
  mid_counties, cdc_data,
  join_by(NAMELSAD == County, STATE_NAME == State)
)
  

#### Generate Maps ####

# Map - Crude Rate of Overdose-Related Deaths between 2018-2024, by County
map_data |> 
  ggplot() + 
  geom_sf(aes(fill = `Crude Rate`), color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    na.value = "gray", name = "Crude Rate\nof Overdose-\nRelated Death"
  ) +
  geom_sf(
    data = mid_states, 
    fill = NA, color = "white", linewidth = 0.5, 
    inherit.aes = F
  ) +
  theme_void() +
  labs(title = "Overdose-Related Deaths by County between 2018-2024")

#### Use geostan to fill in missing counties ####

map_data <- map_data |> 
  # filtering out 4 counties with NA populations until I can find a better
  # solution for these (1 in SD, 3 in IN)
  filter(!is.na(Population))

# Adjacency matrix
C <- shape2mat(map_data, style = "B")

# Quick spatial diagnostics
sp_diag(map_data$`Crude Rate`, map_data, w = C, name = "Crude Rate")

# Prepare list of data for the conditional autoregressive (CAR) model
cars <- prep_car_data(C)

# Fit model
fit <- stan_car(Deaths ~ offset(log(Population)),
                censor_point = 9,
                data = map_data,
                car_parts = cars,
                family = poisson(),
                iter = 2000, # number of MCMC samples
                quiet = TRUE) # to silence printing

# Quick model diagnostics
sp_diag(fit, map_data)

# Summary of model
print(fit)

# Get predictions
pred <- fitted(fit)

# Predicted map
cbind(map_data, pred) |> 
  mutate(mean_per_100k = mean*100000) |> 
  ggplot() +
  geom_sf(aes(fill = mean_per_100k), color = "white", size = 0.2) +
  geom_sf(
    data = mid_states, 
    fill = NA, color = "white", linewidth = 0.5, 
    inherit.aes = F
  ) +
  theme_void() +
  scale_fill_viridis_c(
    na.value = "gray", 
    name = "Predicted \nCrude Rate",
    limits = c(0, 120)
  ) +
  labs(title = "Overdose-Related Deaths by County between 2018-2024")

# Observed map, again
map_data |> 
  ggplot() + 
  geom_sf(aes(fill = `Crude Rate`), color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    na.value = "gray", 
    name = "Observed \nCrude Rate",
    limits = c(0, 120)
  ) +
  geom_sf(
    data = mid_states, 
    fill = NA, color = "white", linewidth = 0.5, 
    inherit.aes = F
  ) +
  theme_void() +
  labs(title = "Overdose-Related Deaths by County between 2018-2024")

#### Clustering ####

# Also have to filter out the 4 counties missing population again here...
mid_counties <- mid_counties |> 
  filter(GEOID %in% map_data$GEOID)

neighbors <- poly2nb(mid_counties)

# Calculate weights, standardized by row
weights <- nb2listw(neighbors, style = "W")

pred_map_data <- cbind(map_data, pred)

# Calculate global Moran's I
# Global Moran's I gives a clustering value for the state overall
g_moran_i <- moran.test(pred_map_data$mean, weights)
g_moran_i

moran.plot(x = pred_map_data$mean, 
           listw = weights, 
           labels = pred_map_data$NAME, 
           ylab = "Spatially-lagged mean crude rate")

# Local Moran's I gives a clustering value for each county
l_moran_i <- localmoran(pred_map_data$mean, weights)

# LISA Analysis - following: 
# https://mkram01.github.io/EPI563-SpatialEPI/spatial-structure-and-clustering-i-morans-i-and-lisa.html
#
# Another resource:
# https://rpubs.com/Jqedwards/808832

pred_map_data$lag <- lag.listw(weights, var = pred_map_data$mean)
pred_map_data$lm_pv <- l_moran_i[,5]

# Create a new dataset that includes standardized values, and then creates a new variable
# 'lm_quad' which takes on the above categorical values.
pred_map_data <- pred_map_data |> 
  mutate(raw_std = as.numeric(scale(mean)), # scale means standardize to mean 0, 1 SD
         lag_std = as.numeric(scale(lag)),
         lm_quad = factor(case_when(  # All of this is assigning labels based on values
           raw_std >= 0 & lag_std >= 0 & lm_pv < 0.05 ~ 'High-High',
           raw_std <= 0 & lag_std <= 0 & lm_pv < 0.05 ~ 'Low-Low',
           raw_std <= 0 & lag_std >= 0 & lm_pv < 0.05 ~ 'Low-High',
           raw_std >= 0 & lag_std <= 0 & lm_pv < 0.05 ~ 'High-Low',
           lm_pv >= 0.05 ~ 'Non-significant'),
           levels = c('High-High','Low-Low','Low-High','High-Low','Non-significant')))

pred_map_data |> 
  ggplot() +
  geom_sf(aes(fill = lm_quad), color = "white", size = 0.2) +
  geom_sf(
    data = mid_states, 
    fill = NA, color = "white", linewidth = 0.5, 
    inherit.aes = F
  ) +
  theme_void()

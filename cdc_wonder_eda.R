# cdc_wonder_eda.R
#
# A script exploring data on overdose-related deaths by county in Michigan.
# Overdose-related deaths were identified using ICD-10 codes X40–X45 
# (Accidental), X60–X65 (Suicide), X85 (Homicide), and Y10–Y15 (Unknown). 
# These data span 2018-2024.


#### Setup ####

library(here)
library(readxl)
library(sf)
library(tigris)
library(ggplot2)
library(tidyverse)
library(spdep)
# pak::pak("connordonegan/geostan")
library(geostan)

i_am("cdc_wonder_eda.R")

# Path to data in my Dropbox
data_path <- yaml::read_yaml(here("auth.yaml"))


#### Read and process data ####

# Get county outlines (cb = T uses lower res file)
mi_counties <- counties(state = "MI", cb = T, year = 2024)

# CDC Wonder Data

# 2018-2024
cdc_data <- read_xlsx(
  here(data_path, "overdose_deaths_by_county_2018-2024.xlsx"),
  range = "B1:H84"
)

# 1999-2008
# cdc_data <- read_xlsx(
#   here(data_path, "overdose_deaths_by_county_1999-2008.xlsx"),
#   range = "B1:F84"
# )

# 2009-2017
# cdc_data <- read_xlsx(
#   here(data_path, "overdose_deaths_by_county_2009-2017.xlsx"),
#   range = "B1:F84"
# )

# 2024 Substance Use Vulnerability Index (SUVI)
# PROBLEM: They split up "Out-Wayne" and the City of Detroit
# TEMP FIX: Using City of Detroit SUVI for all of Wayne county
# mi_county_suvi <- read_xlsx(
#   here(data_path, "2024 MI County Substance Use Vulnerability Index Results.xlsx"),
#   sheet = "County Ranks",
#   range = "A3:M88"
# ) |> 
#   mutate(
#     # Change county names with "Saint" to "St." for later merging
#     County = County |> str_replace("Saint", "St.") |> 
#       # TEMP FIX: Using City of Detroit SUVI for all of Wayne county
#       replace_when(County == "City of Detroit" ~ "Wayne")
#     ) |> 
#   filter(County %in% mi_counties$NAME) |> 
#   select(County, `MI-SUVI Score (2024)`)

# Process data
cdc_data <- cdc_data |> 
  # CDC has counties formatted as "County, MI" - remove ", MI"
  mutate(
    County = County |> str_remove_all(", MI"),
    # Convert "Suppressed" to NA
    across(
      .cols = Deaths:`Crude Rate Upper 95% Confidence Interval`, 
      .fns = ~replace_when(.x, .x %in% c("Suppressed", "Unreliable") ~ NA)
    ),
    # Change some data columns from chr to num
    across(
      .cols = Deaths:`Crude Rate Upper 95% Confidence Interval`, 
      .fns = as.numeric
    )
  )

# Processing data where counties are "X County, MI" but cdc is just "X"
cdc_data <- cdc_data |> 
  # CDC has counties formatted as "County, MI" - remove ", MI"
  mutate(
    County = County |> str_remove_all(" County, MI"),
    # Convert "Suppressed" to NA
    across(
      .cols = Deaths:`Crude Rate`, 
      .fns = ~replace_when(.x, .x %in% c("Suppressed", "Unreliable") ~ NA)
    ),
    # Change some data columns from chr to num
    across(
      .cols = Deaths:`Crude Rate`, 
      .fns = as.numeric
    )
  )

# Merge cdc data to county
# map_data <- left_join(
#   mi_counties, cdc_data,
#   join_by(NAMELSAD == County)
# ) |> 
#   left_join(
#     mi_county_suvi, 
#     join_by(NAME == County)
#   )

map_data <- left_join(
  mi_counties, cdc_data,
  join_by(NAMELSAD == County)
)

#### Generate maps ####

# Map - Crude Rate of Overdose-Related Deaths between 2018-2024, by County
map_data |> 
  ggplot() +
  geom_sf(aes(fill = `Crude Rate`), color = "white", size = 0.2) +
  theme_void() +
  scale_fill_viridis_c(na.value = "gray") +
  labs(title = "Overdose-Related Deaths by County between 2018-2024")

# Map - Counts of Overdose-Related Deaths between 2018-2024, by County
# map_data |> 
#   ggplot() +
#   geom_sf(aes(fill = Deaths), color = "white", size = 0.2) +
#   theme_minimal() +
#   scale_fill_viridis_c(name = "Counts") +
#   labs(title = "Overdose-Related Deaths by County between 2018-2024")


#### Use geostan to fill in missing counties ####

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
  theme_void() +
  scale_fill_viridis_c(
    na.value = "gray", 
    name = "Predicted \nCrude Rate",
    limits = c(8, 51)
  ) +
  labs(title = "Overdose-Related Deaths by County between 2018-2024")

# Observed map, again
map_data |> 
  ggplot() + 
  geom_sf(aes(fill = `Crude Rate`), color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    na.value = "gray", 
    name = "Observed \nCrude Rate",
    limits = c(8, 51)
  ) +
  theme_void() +
  labs(title = "Overdose-Related Deaths by County between 2018-2024")


#### Explore Clustering ####

# Make neighbors graph - which counties share a border?
neighbors <- poly2nb(mi_counties)

# Calculate weights, standardized by row
weights <- nb2listw(neighbors, style = "W")

# Join prediction data to map data
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
  theme_minimal()

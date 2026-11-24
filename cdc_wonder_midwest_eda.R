# cdc_wonder_midwest_eda.R
#
# A script exploring CDC Wonder data on overdose-related deaths by county in the 
# U.S. Midwest. Overdose-related deaths were identified using ICD-10 codes 
# X40–X45 (Accidental), X60–X65 (Suicide), X85 (Homicide), and 
# Y10–Y15 (Unknown). The datasets explored span 2018-2024 (initial exploratory
# analysis), 2015-2019 (pre-covid period), and 2020-2024 (post-covid period).

#### Setup ####

library(here)
library(readxl)
library(sf)
library(tigris)
library(ggplot2)
library(tidyverse)
library(tidycensus)
library(spdep)
# pak::pak("connordonegan/geostan")
library(geostan)

# Auth
auth <- yaml::read_yaml(here("auth.yaml"))
data_path <- auth$data_path
# NOTE: The code below only needs to be run once. You may need to restart R 
#       before using your API key for the first time.
# census_api_key(key = auth$api, install = TRUE)


#### Read and process data ####

# 2018-2024 CDC Wonder Data
cdc_data_18_24 <- read_xlsx(
  here(data_path, "overdose_deaths_by_county_midwest_2018-2024.xlsx"),
  range = "B1:L1056"
)

# 2015-2019 CDC Wonder Data
cdc_data_15_19 <- read_xlsx(
  here(data_path, "overdose_deaths_by_county_midwest_2015-2019.xlsx"),
  range = "B1:J1056"
)

# 2020-2024 CDC Wonder Data
cdc_data_20_24 <- read_xlsx(
  here(data_path, "overdose_deaths_by_county_midwest_2020-2024.xlsx"),
  range = "B1:L1056"
)

# Process data
# Helper function that cleans the CDC Wonder data
clean_cdc <- function(cdc_data) {
  cdc_data <- cdc_data |> 
    mutate(
      County = County |> 
        # CDC has counties formatted as "XYZ County, ST" - remove " County, ST"
        str_remove_all(" County, .."),
      County = County |> 
        # Three county names are spelled differently in CDC data vs. tigris 
        str_replace_all("Lagrange", "LaGrange") |> 
        str_replace_all("La Porte", "LaPorte") |> 
        str_replace_all("De Kalb", "DeKalb") |> 
        # One county in SD underwent a name and GEOID change
        replace_when(
          County == "Shannon" & State == "South Dakota" ~ 
            "Oglala Lakota"
        ),
      # One county in SD underwent a name and GEOID change
      `County Code` = `County Code` |> replace_when(
        County == "Oglala Lakota" ~ 46102
      ) |> as.character(),
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
  
  return(cdc_data)
}

cdc_data_18_24 <- clean_cdc(cdc_data_18_24)
# NOTE: 15_19 data has many crude rates suppressed as "Unreliable"
#       Overriding this by manually re-calculating crude rate
cdc_data_15_19 <- clean_cdc(cdc_data_15_19) |> 
  mutate(`Crude Rate` = round(Deaths/Population*100000, 1))
cdc_data_20_24 <- clean_cdc(cdc_data_20_24)

# Midwest states, by year
state_names <- cdc_data_18_24$State |> unique()
mid_states24 <- states(cb = T, year = 2024) |> 
  filter(NAME %in% state_names)
mid_states19 <- states(cb = T, year = 2019) |> 
  filter(NAME %in% state_names)

# Midwest counties, by year
mid_counties24 <- counties(state = state_names, cb = T, year = 2024)
mid_counties19 <- counties(state = state_names, cb = T, year = 2019)

# Merge cdc data to county
map_data_18_24 <- left_join(
  mid_counties24, cdc_data_18_24,
  join_by(GEOID == `County Code`)
)

map_data_15_19 <- left_join(
  mid_counties19, cdc_data_15_19,
  join_by(GEOID == `County Code`)
)

map_data_20_24 <- left_join(
  mid_counties24, cdc_data_20_24,
  join_by(GEOID == `County Code`)
)


#### Generate Maps ####

# Helper function that generates a figure of crude OD rates by county
crude_fig <- function(map_data, states, years_label) {
  map_data |> 
    ggplot() + 
    geom_sf(aes(fill = `Crude Rate`), color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(
      na.value = "gray", name = "Crude Rate\nof Overdose-\nRelated Death"
    ) +
    geom_sf(
      data = states,
      fill = NA, color = "white", linewidth = 0.5,
      inherit.aes = F
    ) +
    theme_void() +
    labs(title = paste("Overdose-Related Deaths by County between", years_label))
}

# Map - Crude Rate of Overdose-Related Deaths, by County
crude_fig(map_data_18_24, mid_states24, "2018-2024")
crude_fig(map_data_15_19, mid_states19, "2015-2019")
crude_fig(map_data_20_24, mid_states24, "2020-2024")


#### Use geostan to fill in missing counties ####

# Adjacency matrix
# C <- shape2mat(map_data_15_19, style = "B")
C_19 <- shape2mat(mid_counties19, style = "B")
C_24 <- shape2mat(mid_counties24, style = "B")

# Quick spatial diagnostics
sp_diag(map_data_18_24$`Crude Rate`, map_data_18_24, w = C_24, name = "Crude Rate")
sp_diag(map_data_15_19$`Crude Rate`, map_data_15_19, w = C_19, name = "Crude Rate")
sp_diag(map_data_20_24$`Crude Rate`, map_data_20_24, w = C_24, name = "Crude Rate")

# Prepare list of data for the conditional autoregressive (CAR) model
cars_19 <- prep_car_data(C_19)
cars_24 <- prep_car_data(C_24)

# Fit model
fit_18_24 <- stan_car(Deaths ~ offset(log(Population)),
                      censor_point = 9,
                      data = map_data_18_24,
                      car_parts = cars_24,
                      family = poisson(),
                      iter = 2000, # number of MCMC samples
                      quiet = TRUE) # to silence printing

fit_15_19 <- stan_car(Deaths ~ offset(log(Population)),
                      censor_point = 9,
                      data = map_data_15_19,
                      car_parts = cars_19,
                      family = poisson(),
                      iter = 2000, # number of MCMC samples
                      quiet = TRUE) # to silence printing

fit_20_24 <- stan_car(Deaths ~ offset(log(Population)),
                      censor_point = 9,
                      data = map_data_20_24,
                      car_parts = cars_24,
                      family = poisson(),
                      iter = 2000, # number of MCMC samples
                      quiet = TRUE) # to silence printing

# Quick model diagnostics
sp_diag(fit_18_24, map_data_18_24)
sp_diag(fit_15_19, map_data_15_19)
sp_diag(fit_20_24, map_data_20_24)

# Summary of model
print(fit_18_24)
print(fit_15_19)
print(fit_20_24)

# Get predictions
pred_map_data_18_24 <- cbind(map_data_18_24, fitted(fit_18_24))
pred_map_data_15_19 <- cbind(map_data_15_19, fitted(fit_15_19))
pred_map_data_20_24 <- cbind(map_data_20_24, fitted(fit_20_24))

# Predicted map
# Helper function that generates a figure of predicted crude OD rates by county
pred_fig <- function(pred_map_data, states, years_label) {
  pred_map_data |> 
    mutate(mean_per_100k = mean*100000) |> 
    ggplot() +
    geom_sf(aes(fill = mean_per_100k), color = "white", size = 0.2) +
    geom_sf(
      data = states, 
      fill = NA, color = "white", linewidth = 0.5, 
      inherit.aes = F
    ) +
    theme_void() +
    scale_fill_viridis_c(
      na.value = "gray", 
      name = "Predicted \nCrude Rate"
      # limits = c(0, 120)
    ) +
    labs(title = paste("Overdose-Related Deaths by County between", years_label))
}

pred_fig(pred_map_data_18_24, mid_states24, "2018-2024")
pred_fig(pred_map_data_15_19, mid_states19, "2015-2019")
pred_fig(pred_map_data_20_24, mid_states24, "2020-2024")


#### Clustering ####

# Generate neighbors list
neighbors19 <- poly2nb(mid_counties19)
neighbors24 <- poly2nb(mid_counties24)

# Calculate weights, standardized by row
weights19 <- nb2listw(neighbors19, style = "W")
weights24 <- nb2listw(neighbors24, style = "W")

# Calculate global Moran's I
# Global Moran's I gives a clustering value for the state overall
g_moran_i_18_24 <- moran.test(pred_map_data_18_24$mean, weights24)
g_moran_i_18_24
g_moran_i_15_19 <- moran.test(pred_map_data_15_19$mean, weights19)
g_moran_i_15_19
g_moran_i_20_24 <- moran.test(pred_map_data_20_24$mean, weights24)
g_moran_i_20_24

moran.plot(x = pred_map_data_18_24$mean, 
           listw = weights24, 
           labels = pred_map_data_18_24$NAME, 
           ylab = "Spatially-lagged mean crude rate")
moran.plot(x = pred_map_data_15_19$mean, 
           listw = weights19, 
           labels = pred_map_data_15_19$NAME, 
           ylab = "Spatially-lagged mean crude rate")
moran.plot(x = pred_map_data_20_24$mean, 
           listw = weights24, 
           labels = pred_map_data_20_24$NAME, 
           ylab = "Spatially-lagged mean crude rate")


# Local Moran's I gives a clustering value for each county
l_moran_i_18_24 <- localmoran(pred_map_data_18_24$mean, weights24)
l_moran_i_15_19 <- localmoran(pred_map_data_15_19$mean, weights19)
l_moran_i_20_24 <- localmoran(pred_map_data_20_24$mean, weights24)

# LISA Analysis - following: 
# https://mkram01.github.io/EPI563-SpatialEPI/spatial-structure-and-clustering-i-morans-i-and-lisa.html
#
# Another resource:
# https://rpubs.com/Jqedwards/808832

# A helper function to prep data for LISA analysis
prep_lisa <- function(pred_map_data, weights, l_moran_i) {
  pred_map_data |> 
    mutate(
      lag = lag.listw(weights, var = pred_map_data$mean),
      lm_pv = l_moran_i[,5],
      raw_std = as.numeric(scale(mean)), # scale means standardize to mean 0, 1 SD
      lag_std = as.numeric(scale(lag)),
      lm_quad = factor(case_when(  # All of this is assigning labels based on values
        raw_std >= 0 & lag_std >= 0 & lm_pv < 0.05 ~ 'High-High',
        raw_std <= 0 & lag_std <= 0 & lm_pv < 0.05 ~ 'Low-Low',
        raw_std <= 0 & lag_std >= 0 & lm_pv < 0.05 ~ 'Low-High',
        raw_std >= 0 & lag_std <= 0 & lm_pv < 0.05 ~ 'High-Low',
        lm_pv >= 0.05 ~ 'Non-significant'),
        levels = c('High-High','Low-Low','High-Low','Low-High','Non-significant'))
    )
}

# A helper function to plot the LISA results
lisa_fig <- function(lisa_data, states, years_label) {
  lisa_data |> 
    ggplot() +
    geom_sf(aes(fill = lm_quad), color = "white", size = 0.2, show.legend = T) +
    geom_sf(
      data = states, 
      fill = NA, color = "white", linewidth = 0.5, 
      inherit.aes = F
    ) +
    scale_fill_discrete(
      palette = c("#A61D08", "#0815A6", "#D67D78", "#787FD6", "gray75"), 
      drop = F,
      name = "LISA \nClusters"
    ) +
    theme_void() +
    labs(title = paste("LISA Analysis of Overdose Deaths between", years_label))
}

# LISA analysis figures
prep_lisa(pred_map_data_18_24, weights24, l_moran_i_18_24) |> 
  lisa_fig(mid_states24, "2018-2024")

prep_lisa(pred_map_data_15_19, weights24, l_moran_i_15_19) |> 
  lisa_fig(mid_states19, "2015-2019")

prep_lisa(pred_map_data_20_24, weights24, l_moran_i_20_24) |> 
  lisa_fig(mid_states24, "2020-2024")

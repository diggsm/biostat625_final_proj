# cdc_wonder_eda.R
#
# A script exploring data on overdose-related deaths by county in Michigan.
# Overdose-related deaths were identified using ICD-10 codes X40–X45 
# (Accidental), X60–X65 (Suicide), X85 (Homicide), and Y10–Y15 (Unknown). 
# These data span 2018-2024.


#### Setup ####

# Install INLA package, which is not hosted on CRAN
# install.packages(
#   "INLA",
#   repos = "https://inla.r-inla-download.org/R/stable", 
#   dep = TRUE
# )

# May also need to install package fmesher
# install.packages("fmesher")

library(here)
library(readxl)
library(sf)
library(tigris)
library(ggplot2)
library(tidyverse)
library(INLA)
library(spdep)

i_am("cdc_wonder_eda.R")

# Path to data in my Dropbox
data_path <- yaml::read_yaml(here("auth.yaml"))


#### Read and process data ####

# Get county outlines (cb = T uses lower res file)
mi_counties <- counties(state = "MI", cb = T)

# CDC Wonder Data
cdc_data <- read_xlsx(
  here(data_path, "overdose_deaths_by_county_2018-2024.xlsx"),
  range = "B1:H84"
  )

# 2024 Substance Use Vulnerability Index (SUVI)
# PROBLEM: They split up "Out-Wayne" and the City of Detroit
# TEMP FIX: Using City of Detroit SUVI for all of Wayne county
mi_county_suvi <- read_xlsx(
  here(data_path, "2024 MI County Substance Use Vulnerability Index Results.xlsx"),
  sheet = "County Ranks",
  range = "A3:M88"
) |> 
  mutate(
    # Change county names with "Saint" to "St." for later merging
    County = County |> str_replace("Saint", "St.") |> 
      # TEMP FIX: Using City of Detroit SUVI for all of Wayne county
      replace_when(County == "City of Detroit" ~ "Wayne")
    ) |> 
  filter(County %in% mi_counties$NAME) |> 
  select(County, `MI-SUVI Score (2024)`)

# Process data
cdc_data <- cdc_data |> 
  # CDC has counties formatted as "County, MI" - remove ", MI"
  mutate(
    County = County |> str_remove_all(", MI"),
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
  mi_counties, cdc_data,
  join_by(NAMELSAD == County)
) |> 
  left_join(
    mi_county_suvi, 
    join_by(NAME == County)
  )


#### Generate maps ####

# Map - Crude Rate of Overdose-Related Deaths between 2018-2024, by County
map_data |> 
  ggplot() +
  geom_sf(aes(fill = `Crude Rate`), color = "white", size = 0.2) +
  theme_minimal() +
  scale_fill_viridis_c() +
  labs(title = "Overdose-Related Deaths by County between 2018-2024")

# Map - Counts of Overdose-Related Deaths between 2018-2024, by County
map_data |> 
  ggplot() +
  geom_sf(aes(fill = Deaths), color = "white", size = 0.2) +
  theme_minimal() +
  scale_fill_viridis_c(name = "Counts") +
  labs(title = "Overdose-Related Deaths by County between 2018-2024")


#### Bayesian Hierarchical Model? ####

# Map Data - adding id variable
map_data <- map_data |> 
  mutate(
    # id column
    id = row_number()
  ) |> 
  rename(suvi24 = `MI-SUVI Score (2024)`)

# Make neighbors graph - which counties share a border?
mi_counties <- mi_counties |> 
  arrange(GEOID)

neighbors <- poly2nb(mi_counties)
nb2INLA("county_neighbors.adj", neighbors)


# Fit Poisson model with BYM spatial random effect and population offset
# Deaths explained by substance use vulnerability index
# NOTE: Could consider adding other demographic variables as predictors
bhm_deaths_suvi <- inla(
  Deaths ~ suvi24 + offset(log(Population)) + 
    f(id, model = "bym2", graph = "county_neighbors.adj"), 
  family = "poisson", 
  data = map_data
)

# Process data
bhm_map_data <- cbind(
  map_data,
  bhm_deaths_suvi$summary.fitted.values
) |> 
  mutate(
    bhm_crude_rate = mean/Population*100000
  )

# Plot results
bhm_map_data |> 
  ggplot() +
  geom_sf(aes(fill = bhm_crude_rate), color = "white", size = 0.2) +
  theme_minimal() +
  scale_fill_viridis_c(name = "Predicted \nCrude Death \nRate") +
  labs(title = "BHM Predicted Overdose-Related Crude Death Rate \nby County between 2018-2024")

# app.R
#
# A Shiny app to explore the results of our analysis.

library(shiny)
library(shinyWidgets)
library(bslib)
library(here)
library(tidyverse)
library(tigris)

# Auth
auth <- yaml::read_yaml(here("auth.yaml"))
data_path <- auth$data_path

#### Read in data ####

results_18_24 <- readRDS(here(data_path, "results_18_24.rds"))
results_15_19 <- readRDS(here(data_path, "results_15_19.rds"))
results_20_24 <- readRDS(here(data_path, "results_20_24.rds"))

results <- bind_rows(
  list(
    "2018-2024" = results_18_24, 
    "2015-2019" = results_15_19, 
    "2020-2024" = results_20_24
  ),
  .id = "years"
)

results <- results |> 
  mutate(
    mean_per_100k = mean*100000,
    metric = `Crude.Rate`
  )

mid_states24 <- states(cb = T, year = 2024) |> 
  filter(NAME %in% unique(results$State))
mid_states19 <- states(cb = T, year = 2019) |> 
  filter(NAME %in% unique(results$State))


#### UI ####

ui <- page_fluid(
  page_navbar(
    title = "A Spatial Analysis of Overdose Deaths in the U.S. Midwest",
    id = "page",
    theme = bs_theme(bootswatch = "flatly", version = 5),
    
    ##### Sidebar #####
    
    sidebar = sidebar(
      open = "always",
      width = "275px",
      title = "Inputs",
      
      conditionalPanel(
        "input.page === 'main_panel'",
        # Years
        radioButtons(
          "years", 
          label = "Time Frame",
          choices = c("2018-2024", "2015-2019", "2020-2024"), 
          selected = "2018-2024"
        ),
        # Input: Disease
        pickerInput(
          "state",
          "States",
          choices = sort(unique(results$State)),
          options = pickerOptions(
            liveSearch = TRUE,
            actionsBox = TRUE,
            deselectAllText = "Unselect All",
            selectAllText = "Select All"
          ),
          multiple = TRUE,
          selected = unique(results$State)
        ),
        # Input: Metric
        radioButtons(
          "metric", 
          label = "Metric",
          choices = c(
            "Observed Data", 
            "Predicted Data", 
            "LISA Analysis"
          ), 
          selected = "Observed Data"
        )
      )
    ),
    
    ##### Map Figure #####
    
    nav_panel(
      "Home",
      value = "main_panel",
      # Figure
      card(plotOutput("map_fig"))
    )
    
    ##### Other Panels #####
    # Could add a section to display global moran's I results?
  )
)


#### Server ####

server <- function(input, output, session) {
  
  ##### Reactives ####
  
  # States selected
  states_selected <- reactive({input$state})
  mid_states_filtered <- reactive({
    mid_states() |> 
      filter(NAME %in% input$state)
  })
  
  # Years selected
  years_selected <- reactive({input$years})
  results_filtered <- reactive({
    results |> 
      filter(years == input$years)
  })
  
  # Version of midwest states
  mid_states <- reactive({
    if (input$years %in% c("2018-2024", "2020-2024")) {
      mid_states24
    } else if (input$years == "2015-2019") {
      mid_states19
    }
  })

  # Type of metric to display
  metric_data <- reactive({
    if (input$metric == "Observed Data") {
      results_filtered()$Crude.Rate
    } else if (input$metric == "Predicted Data") {
      results_filtered()$mean_per_100k
    } else if (input$metric == "LISA Analysis") {
      results_filtered()$lm_quad
    }
  })
  
  # Max rates for color scales
  max_rate <- reactive({
    if (input$metric == "Observed Data") {
      max(results$Crude.Rate, na.rm = T)
    } else if (input$metric == "Predicted Data") {
      max(results$mean_per_100k, na.rm = T)
    } 
  })
  
  # Map Figure
  output$map_fig <- renderPlot({
    p <- results_filtered() |>
      mutate(metric = metric_data()) |> 
      filter(State %in% states_selected()) |> 
      ggplot() +
      geom_sf(aes(fill = metric), color = "white", size = 0.2, show.legend = T) +
      geom_sf(
        data = mid_states_filtered(),
        fill = NA, color = "white", linewidth = 0.7,
        inherit.aes = F
      ) +
      theme_void() +
      labs(title = paste("Overdose-Related Deaths by County between", years_selected()))
    
    if (input$metric == "Observed Data") {
      p <- p +
        scale_fill_viridis_c(
          na.value = "gray",
          name = "Observed Crude \nOverdose Rates",
          limits = c(0, max_rate())
        )
    }
    
    if (input$metric == "Predicted Data") {
      p <- p +
        scale_fill_viridis_c(
          na.value = "gray",
          name = "Predicted Crude \nOverdose Rates",
          limits = c(0, max_rate())
        )
    }
    
    if (input$metric == "LISA Analysis") {
      p <- p +
        scale_fill_discrete(
          palette = c("#A61D08", "#0815A6", "#D67D78", "#787FD6", "gray75"), 
          drop = F,
          name = "LISA \nClusters"
        )
    }
    p
  })
  
}

#### Run App ####

shinyApp(ui = ui, server = server)

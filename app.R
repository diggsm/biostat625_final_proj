# app.R
#
# A Shiny app to explore the results of our analysis.

library(shiny)
library(shinyWidgets)
library(bslib)
library(here)
library(tidyverse)
library(tigris)


#### Read in data ####

results <- readRDS(here("midwest_predictions.rds"))

results <- results |> 
  mutate(
    mean_per_100k = mean*100000,
    metric = `Crude.Rate`
  )

mid_states <- states(cb = T, year = 2024) |> 
  filter(NAME %in% c("Michigan", "Illinois", "Indiana", "Iowa", "Kansas",
                     "Minnesota", "Missouri", "Nebrasha", "North Dakota",
                     "Ohio", "South Dakota", "Wisconsin"))


#### UI ####

ui <- page_fluid(
  page_navbar(
    title = "Our Final Project",
    id = "page",
    theme = bs_theme(bootswatch = "flatly", version = 5),
    
    ##### Sidebar #####
    
    sidebar = sidebar(
      open = "always",
      width = "275px",
      title = "Inputs",
      
      conditionalPanel(
        "input.page === 'main_panel'",
        # Input: Disease
        pickerInput(
          "state",
          "States",
          choices = sort(unique(results$STATE_NAME)),
          options = pickerOptions(
            liveSearch = TRUE,
            actionsBox = TRUE,
            deselectAllText = "Unselect All",
            selectAllText = "Select All"
          ),
          multiple = TRUE,
          selected = unique(results$STATE_NAME)
        ),
        # Input: Metric
        radioButtons(
          "metric", 
          label = "Metric",
          choices = c("Observed", "Predicted", "LISA"), 
          selected = "Observed"
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
    mid_states |> 
      filter(NAME %in% input$state)
  })
  
  
  # Type of metric to display
  metric_data <- reactive({
    if (input$metric == "Observed") {
      results$Crude.Rate
    } else if (input$metric == "Predicted") {
      results$mean_per_100k
    } else if (input$metric == "LISA") {
      results$lm_quad
    }
  })
  
  
  # Map Figure
  output$map_fig <- renderPlot({
    p <- results |> 
      mutate(metric = metric_data()) |> 
      filter(STATE_NAME %in% states_selected()) |> 
      ggplot() +
      geom_sf(aes(fill = metric), color = "white", size = 0.2) +
      geom_sf(
        data = mid_states_filtered(), 
        fill = NA, color = "white", linewidth = 0.7, 
        inherit.aes = F
      ) +
      theme_void() +
      labs(title = "Overdose-Related Deaths by County between 2018-2024")
    
    if (input$metric != "LISA") {
      p <- p +
        scale_fill_viridis_c(
          na.value = "gray",
          name = "Predicted \nCrude Rate"
        )
    }
    p
  })
  
}

#### Run App ####

shinyApp(ui = ui, server = server)

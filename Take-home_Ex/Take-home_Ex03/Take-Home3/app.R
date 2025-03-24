
library(shiny)

choices_ordered <- weather_sf %>%
  mutate(Date = as.Date(paste("01", format(Date, "%b-%Y")), format = "%d %b-%Y")) %>%
  pull(Date) %>%
  unique() %>%
  sort(decreasing = TRUE)

# Convert to Month-Year format for display
choices_ordered <- format(choices_ordered, "%b-%Y")

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Weather Data Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("month_year", "Select Month-Year:",
                  choices = choices_ordered,
                  selected = choices_ordered[1],  # Default is the most recent
                  selectize = TRUE)
    ),
    
    mainPanel(
      plotOutput("weatherPlot")
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive expression to load raw weather data (assuming `weather_sf` is loaded already)
  weather_data <- reactive({
    weather_sf %>%
      filter(Year == 2021, Month == 4)
  })
  
  # Reactive expression to aggregate data based on user input
  aggregated_weather <- reactive({
    req(input$time_resolution)  # Ensure input is available
    
    if (input$time_resolution == "Monthly") {
      weather_data() %>%
        group_by(Station, Year, Month) %>%
        summarise(MonthlyRainfall = sum(DailyRainfall, na.rm = TRUE),
                  MonthlyTemp = mean(MeanTemperature, na.rm = TRUE)) %>%
        ungroup()
    } else if (input$time_resolution == "Yearly") {
      weather_data() %>%
        group_by(Station, Year) %>%
        summarise(YearlyRainfall = sum(DailyRainfall, na.rm = TRUE),
                  YearlyTemp = mean(MeanTemperature, na.rm = TRUE)) %>%
        ungroup()
    } else {
      # Daily resolution (no aggregation)
      weather_data()
    }
  })
  
  # Output: Plot or other UI element based on aggregated data
  output$weatherPlot <- renderPlot({
    data_to_plot <- aggregated_weather()
    
    ggplot(data_to_plot, aes(x = Station, y = MonthlyRainfall)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = paste("Weather Data -", input$time_resolution, "Resolution"))
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)

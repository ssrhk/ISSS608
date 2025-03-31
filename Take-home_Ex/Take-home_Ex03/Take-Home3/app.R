pacman::p_load(shiny,readr,dplyr,tidyverse,SmartEDA,sf,plotly,tmap,terra,gstat,automap)

weather_filtered <- read_rds("rds/weather_filtered.rds")
mpsz2019 <- read_rds("rds/mpsz2019.rds")

ui <- fluidPage(
  
  titlePanel("Spatial Interpolation of Weather Data"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("variable", 
                  "Select Variable:", 
                  choices = c("Total Rainfall (mm)" = "MonthlyRainfall",
                              "Mean Temperature (°C)" = "MonthlyMeanTemp",
                              "Mean Wind Speed (km/h)" = "MonthlyMeanWindSpeed"),
                  selected = "MonthlyMeanTemp"),
      
      selectInput("month_year", 
                  "Select Month-Year:", 
                  choices = format(seq(as.Date("2021-01-01"), as.Date("2024-04-01"), by = "month"), "%b-%Y"),
                  selected = "Jan-2021"),
      
      # Tabset for Auto and Manual
      tabsetPanel(
        id = "tabset_variogram",
        tabPanel("Automatic Variogram",
                 h4("Auto-Fitted Variogram"),
                 p("Auto-fitted variogram automatically determines the optimal variogram model and parameters (psill, range, nugget) based on the data, so users don't need to manually input these values.")
        ),
        tabPanel("Manual Adjustment", 
                 h4("Manual Variogram Adjustment"),
                 sliderInput("psill", "Psill:", min = 0, max = 5, value = 0.5),
                 selectInput("model", "Model Type:",
                             choices = c("Spherical" = "Sph", 
                                         "Exponential" = "Exp", 
                                         "Gaussian" = "Gau"), 
                             selected = "Sph"),
                 sliderInput("range", "Range:", min = 100, max = 10000, value = 5000),
                 sliderInput("nugget", "Nugget:", min = 0, max = 1, value = 0.1)
        )
      ),
      actionButton("show_result", "Show Result")
    ),
    
    mainPanel(
      fluidRow(
        column(6, tmapOutput("mean_temp_plot")),
        column(6, tmapOutput("variance_plot"))
      ),
       
      wellPanel(
        h4("Introduction to Spatial Interpolation"),
        p("Spatial interpolation is the process of predicting values for unmeasured locations based on known values from nearby locations."),
        p("Kriging is a method of spatial interpolation that uses statistical models to predict the value of a variable at unmeasured points, taking into account the spatial correlation between data points."),
        h4("How to Interpret the Plots"),
        p("The first plot represents the predicted values (e.g., temperature, rainfall, or wind speed) for the selected variable at different locations."),
        p("The second plot shows the Kriging variance, which represents the uncertainty of the predictions. Higher variance indicates less confidence in the prediction."),
        p("These plots allow you to visualize the distribution of weather data across Singapore for a specific month and year, and understand the spatial patterns and uncertainties.")
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Create the grid data
  grid <- terra::rast(mpsz2019, nrows = 400, ncols = 700)  
  
  # Create xy from the grid
  xy <- terra::xyFromCell(grid, 1:ncell(grid))
  
  # Create coop spatial points data frame 
  coop <- st_as_sf(as.data.frame(xy), coords = c("x", "y"), crs = st_crs(mpsz2019))
  coop <- st_filter(coop, mpsz2019)
  
  # Reactive dataset based on selected Month-Year
  weather_month <- eventReactive(input$show_result, {
    selected_date <- as.Date(paste("01", input$month_year), format = "%d %b-%Y")
    selected_year <- format(selected_date, "%Y")
    selected_month <- format(selected_date, "%m")
    
    weather_filtered %>%
      filter(Year == as.numeric(selected_year), Month == as.numeric(selected_month)) %>%
      group_by(Station) %>%
      summarise(
        MonthlyRainfall = sum(DailyRainfall, na.rm = TRUE),
        MonthlyMeanTemp = mean(MeanTemperature, na.rm = TRUE),
        MonthlyMeanWindSpeed = mean(MeanWindSpeed, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      filter(Station %in% c("Admiralty", "Ang Mo Kio", "Changi", "Choa Chu Kang (South)", 
                            "East Coast Parkway", "Jurong (West)", "Jurong Island", "Newton", 
                            "Pasir Panjang", "Pulau Ubin", "Seletar", "Sentosa Island", 
                            "Tai Seng", "Tuas South"))
  })
  
  # Kriging 
  kriging_results <- eventReactive(input$show_result, {
    req(weather_month())
    
    model_type <- switch(input$model, "Sph" = "Sph", "Exp" = "Exp", "Gau" = "Gau")
    
    # Check which tab is active
    if (input$tabset_variogram == "Automatic Variogram") {
      # Automatically compute variogram model
      auto_vgm <- variogram(as.formula(paste(input$variable, "~ 1")), data = weather_month())
      v_model <- fit.variogram(auto_vgm, vgm(model = model_type))
    } else {
      # Use manually specified values for variogram parameters
      v_model <- vgm(psill = input$psill, model = model_type, range = input$range, nugget = input$nugget)
    }
    
    krige_model <- gstat(formula = as.formula(paste(input$variable, "~ 1")), 
                         model = v_model, 
                         data = weather_month())
    
    predictions <- predict(krige_model, coop)
    predictions$x <- st_coordinates(predictions)[,1]
    predictions$y <- st_coordinates(predictions)[,2]
    predictions$pred <- predictions$var1.pred
    predictions$variance <- predictions$var1.var
    
    kpred <- terra::rasterize(predictions, grid, field = "pred")
    kpred_var <- terra::rasterize(predictions, grid, field = "variance")
    
    list(pred_raster = kpred, var_raster = kpred_var, selected_variable = input$variable, selected_month_year = input$month_year)
  })
  
  
  
  # Render Kriging Plot
  output$mean_temp_plot <- renderTmap({
    req(kriging_results())
    
    variable_title <- switch(kriging_results()$selected_variable,
                             "MonthlyMeanTemp" = "Mean Temperature (°C)",
                             "MonthlyRainfall" = "Total Rainfall (mm)",
                             "MonthlyMeanWindSpeed" = "Mean Wind Speed (km/h)")
    
    tm_shape(kriging_results()$pred_raster) + 
      tm_raster(col_alpha = 0.6, palette = "YlOrRd", title = paste(variable_title)) +
      tm_layout(main.title = paste("Distribution of", variable_title, "for", kriging_results()$selected_month_year), frame = TRUE, 
                legend.position = c("left", "top"), legend.frame = FALSE,
                asp = 1) +
      tm_compass(type = "8star", size = 2) +
      tm_grid(alpha = 0.2)
  })
  
  # Render Kriging Variance Plot
  output$variance_plot <- renderTmap({
    req(kriging_results())
    
    variable_title <- switch(kriging_results()$selected_variable,
                             "MonthlyMeanTemp" = "Mean Temperature (°C)",
                             "MonthlyRainfall" = "Total Rainfall (mm)",
                             "MonthlyMeanWindSpeed" = "Mean Wind Speed (km/h)")
    
    tm_shape(kriging_results()$var_raster) + 
      tm_raster(col_alpha = 0.6, palette = "YlGnBu", title = paste("Kriging Variance of", variable_title)) +
      tm_layout(main.title = paste("Kriging Variance of", variable_title, "for", kriging_results()$selected_month_year), frame = TRUE, 
                legend.position = c("left", "top"), legend.frame = FALSE,
                asp = 1) +
      tm_compass(type = "8star", size = 2) +
      tm_grid(alpha = 0.2)
  })
}


shinyApp(ui, server)

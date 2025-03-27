pacman::p_load(shiny, readr, dplyr, tidyverse, sf, sfdep, tmap, terra, gstat, automap)

# Define UI
ui <- fluidPage(
  titlePanel("Local Measure of Spatial Autocorrelation"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("variable", "Select Variable:",
                  choices = c("Total Rainfall (mm)" = "MonthlyRainfall",
                              "Mean Temperature (°C)" = "MonthlyMeanTemp",
                              "Mean Wind Speed (km/h)" = "MonthlyMeanWindSpeed"),
                  selected = "MonthlyMeanTemp"),
      
      selectInput("month_year", "Select Month-Year:",
                  choices = format(seq(as.Date("2021-01-01"), as.Date("2024-04-01"), by = "month"), "%b-%Y"),
                  selected = "Jan-2021"),
      
      selectInput("stat", "Select Statistic:",
                  choices = c("Local Moran I" = "ii", "P-value" = "p_ii", "Std Deviation" = "z_ii", "Variance" = "var_ii", "Expectation" = "eii"),
                  selected = "ii"),
      
      actionButton("show_result", "Show Result"),
      
      # Text box explaining the maps
      br(),
      wellPanel(
        h5("Interpretation of Maps:"),
        p("The first map displays the selected statistic (e.g., Local Moran's I, P-value, Std Deviation, etc.) for the chosen variable (e.g., Total Rainfall, Mean Temperature, Mean Wind Speed). Each color in the map represents the value of the statistic for a particular geographic region, where darker or lighter colors indicate higher or lower values of the statistic, respectively. This allows you to observe the spatial distribution and variability of the statistic across the different regions of Singapore."),
        p("The second map is the LISA (Local Indicators of Spatial Association) map, which shows significant clusters of similar values for the selected variable. LISA highlights areas where high values are clustered together, and similarly, low values are grouped in specific locations. The LISA map helps to identify whether the spatial patterns are random or whether they exhibit some form of spatial autocorrelation (i.e., whether nearby areas tend to have similar or dissimilar values).")
      )
    ),
    
    mainPanel(
      tmapOutput("stat_plot"),
      tmapOutput("lisa_map")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Load required datasets
  weather <- read_rds("rds/weather.rds")
  mpsz2019 <- read_rds("rds/mpsz.rds")
  
  filtered_data <- eventReactive(input$show_result, {
    # Capture all the inputs inside the eventReactive to trigger only on button click
    month_year_selected <- input$month_year
    variable_selected <- input$variable
    stat_selected <- input$stat 
    
    weather_month <- weather %>%
      filter(format(as.Date(paste(Year, Month, "01", sep = "-")), "%b-%Y") == month_year_selected)
    
    keepstations <- c("Admiralty", "Ang Mo Kio", "Changi", "Choa Chu Kang (South)", "East Coast Parkway", 
                      "Jurong (West)", "Jurong Island", "Newton", "Pasir Panjang", "Pulau Ubin", 
                      "Seletar", "Sentosa Island", "Tai Seng", "Tuas South")
    
    weather_month <- weather_month %>% filter(Station %in% keepstations)
    
    # Aggregate data
    weather_aggregated <- weather_month %>%
      group_by(Station) %>%
      summarise(
        MonthlyRainfall = sum(DailyRainfall, na.rm = TRUE),
        MonthlyMeanTemp = mean(MeanTemperature, na.rm = TRUE),
        MonthlyMeanWindSpeed = mean(MeanWindSpeed, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Merge with spatial data
    weather_sf <- left_join(weather_aggregated, mpsz2019, by = "Station") %>%
      st_as_sf()
    
    if (!(variable_selected %in% colnames(weather_sf))) {
      stop(paste("Error: Variable", variable_selected, "not found in dataset."))
    }
    
    formula <- as.formula(paste(variable_selected, "~ 1"))
    
    variogram_model <- variogram(formula, weather_sf)
    fit_model <- tryCatch({
      fit.variogram(variogram_model, model = vgm("Exp"), fit.method = 6)
    }, error = function(e) {
      NULL
    })
    
    if (is.null(fit_model)) {
      stop("Error: Variogram fitting failed.")
    }
    
    # Perform Kriging interpolation
    kriging_result <- krige(formula, weather_sf, mpsz2019, model = fit_model)
    
    # Convert kriging result to sf object
    kriged_sf <- st_as_sf(kriging_result)
    
    # Compute spatial weights using sfdep
    mpsz_centroids <- st_centroid(mpsz2019)
    knn_result <- st_knn(mpsz_centroids, k = 3)
    knn_weights <- st_weights(knn_result, style = "W")
    
    # Compute Local Moran’s I
    local_moran_res <- local_moran(as.numeric(kriged_sf$var1.pred), nb = knn_result, wt = knn_weights, nsim = 99)
    local_moran_df <- as.data.frame(local_moran_res)
    
    kriged_sf <- kriged_sf %>%
      mutate(
        ii = local_moran_df$ii,
        p_ii = local_moran_df$p_ii,
        z_ii = local_moran_df$z_ii,
        var_ii = local_moran_df$var_ii,
        eii = local_moran_df$eii
      )
    
    return(kriged_sf)
  })
  
  output$stat_plot <- renderTmap({
    req(filtered_data())  
    data <- filtered_data()
    
    # Use input values to build titles inside the eventReactive context
    variable_title <- switch(input$variable,
                             "MonthlyMeanTemp" = "Mean Temperature (°C)",
                             "MonthlyRainfall" = "Total Rainfall (mm)",
                             "MonthlyMeanWindSpeed" = "Mean Wind Speed (km/h)")
    
    stat_title <- switch(input$stat,
                         "ii" = "Local Moran I",
                         "p_ii" = "P-value",
                         "z_ii" = "Std Deviation",
                         "var_ii" = "Variance",
                         "eii" = "Expectation")
    
    # Create the map with fixed color theme (green)
    tm_shape(data) +
      tm_fill(input$stat, palette = "brewer.blues", legend.show = TRUE) +  
      tm_borders() +
      tm_title(paste(stat_title, "of", variable_title,",", input$month_year))  
  })
  
  output$lisa_map <- renderTmap({
    req(filtered_data())  
    data <- filtered_data()
    
    # Define variable_title inside the lisa_map output function
    variable_title <- switch(input$variable,
                             "MonthlyMeanTemp" = "Mean Temperature (°C)",
                             "MonthlyRainfall" = "Total Rainfall (mm)",
                             "MonthlyMeanWindSpeed" = "Mean Wind Speed (km/h)")
    
    if (!"p_ii" %in% colnames(data)) {
      stop("Error: p_ii column not found in dataset.")
    }
    
    lisa_sig <- data %>% filter(p_ii < 0.05)
    
    # LISA map with fixed green color theme
    tm_shape(data) +
      tm_polygons() +
      tm_borders(fill_alpha = 0.5) +
      tm_shape(lisa_sig) +
      tm_fill("ii", palette = "brewer.rd_bu", legend.show = TRUE)+ 
      tm_borders(fill_alpha = 0.4) +
      tm_title(paste("LISA Map of", variable_title))
  })
}

shinyApp(ui = ui, server = server)



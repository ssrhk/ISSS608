library(shiny)
library(tmap)
library(sf)
library(dplyr)
library(gstat)
library(spdep)

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
                  selected = "Apr-2021"),
      
      selectInput("stat", "Select Statistic:",
                  choices = c("ii" = "ii", "p_ii" = "p_ii", "z_ii" = "z_ii", "var_ii" = "var_ii", "eii" = "eii"),
                  selected = "ii")
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
  
  # Filter data based on user input
  filtered_data <- reactive({
    month_year_selected <- input$month_year
    variable_selected <- input$variable
    
    weather_month <- weather %>%
      filter(format(as.Date(paste(Year, Month, "01", sep = "-")), "%b-%Y") == month_year_selected)
    
    # Filter for relevant stations
    keepstations <- c("Admiralty", "Ang Mo Kio", "Changi", "Choa Chu Kang (South)", "East Coast Parkway", 
                      "Jurong (West)", "Jurong Island", "Newton", "Pasir Panjang", "Pulau Ubin", 
                      "Seletar", "Sentosa Island", "Tai Seng", "Tuas South")
    
    weather_month <- weather_month %>%
      filter(Station %in% keepstations)
    
    # Convert to spatial data
    weather_sf <- st_as_sf(weather_month, coords = c("Longitude", "Latitude"), crs = 4326)
    
    # Build the formula dynamically based on selected variable
    formula <- as.formula(paste(variable_selected, "~ 1"))
    
    # Fit the variogram model
    variogram_model <- variogram(formula, weather_sf)
    fit_model <- fit.variogram(variogram_model, model = vgm("Exp"))
    
    # Perform Kriging interpolation
    kriging_result <- krige(formula, weather_sf, mpsz2019, model = fit_model)
    
    # Convert kriging result to sf object
    kriged_sf <- st_as_sf(kriging_result)
    
    # Create centroids for the polygons in mpsz2019
    mpsz_centroids <- st_centroid(mpsz2019)
    
    # Apply KNN for spatial weights
    knn_result <- st_knn(mpsz_centroids, k = 3)
    knn_weights <- st_weights(knn_result, style = "W")
    
    # Calculate Local Moran's I using kriged predicted values (var1.pred)
    local_moran_res <- local_moran(kriged_sf$var1.pred, knn_result, knn_weights, nsim = 99)
    local_moran_df <- as.data.frame(local_moran_res)
    
    # Add results to kriged_sf
    kriged_sf <- kriged_sf %>%
      mutate(
        ii = local_moran_df$ii,
        p_ii_sim = local_moran_df$p_ii_sim,
        eii = local_moran_df$eii,
        var_ii = local_moran_df$var_ii,
        z_ii = local_moran_df$z_ii,
        p_ii = local_moran_df$p_ii,
        mean = local_moran_df$mean
      )
    
    kriged_sf
  })
  
  # Render Stat Plot
  output$stat_plot <- renderTmap({
    data <- filtered_data()
    tm_shape(data) +
      tm_fill(input$stat, palette = "Blues", legend.show = TRUE) +  # Use dynamic statistic based on user input
      tm_borders() +
      tm_title(paste(input$stat))
  })
  
  # Render LISA Map
  output$lisa_map <- renderTmap({
    data <- filtered_data()
    # Filter significant Local Moran's I results (p-value < 0.05)
    lisa_sig <- data %>%
      filter(p_ii_sim < 0.05)
    
    tm_shape(data) +
      tm_polygons() +
      tm_borders(fill_alpha = 0.5) +
      tm_shape(lisa_sig) +
      tm_fill("mean", palette = "RdBu", legend.show = TRUE) +  # Use "mean" for LISA significance
      tm_borders(fill_alpha = 0.4) +
      tm_title("LISA Map")
  })
}

# Run the application
shinyApp(ui = ui, server = server)

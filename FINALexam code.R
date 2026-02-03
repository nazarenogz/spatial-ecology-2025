library(rgbif) #Downloading raw GBIF data.
library(sf) #To transform data frames into spatial objects, and other useful functions.
library(spatstat) #To calculate KDE, Lcross and ppp. 
library(viridis) #To use color-blind friendly color palettes in the plots. 
library(rnaturalearth) #To load the Italy map.
library(ggplot2) #To visualize the occurrence and density in Italy. 
library(patchwork) #To snap separate plots into a single grid easily. 

# 1. Spatial Setup
#We load the italian map and immediately transform it to a metric system, better for calculus. 
italy <- ne_countries(country = "Italy", scale = "medium", returnclass = "sf") |>
  st_transform(32632)
#We transform the map into the observation window used later for the ppp
italy_poly <- as.owin(italy)
#Let's visualize the italy map without our data plugged in
ggplot() +
    geom_sf(data = italy, fill = "#f8f9fa", color = "grey80", linewidth = 0.2) + 
theme_minimal() + theme(panel.grid = element_blank())

# 2. Data Loading Function
#We create the function to download and clean the species data. 
load_species_sf <- function(taxonKey, exclude_pattern = NULL) {
  #We trarnsform the downloaded data into a data table with $data, and filter for data that has coordinates and is in Italy.
  data <- occ_search(taxonKey = taxonKey, country = "IT", hasCoordinate = TRUE, limit = 10000)$data
  #We eliminate useless columns and keep just these 3. 
  data <- data[, c("decimalLongitude", "decimalLatitude", "scientificName")]
  #We filter to eliminate NA rows. 
  data <- data[!is.na(data$decimalLongitude) & !is.na(data$decimalLatitude), ]
  #The grep1 function with ! keeps the rows without the specified pattern, so we can filter for domestic animals. 
  #ignore.case is true to ensure that the pattern ignores upper or lower cases. 
  if (!is.null(exclude_pattern)) {
    data <- data[!grepl(exclude_pattern, data$scientificName, ignore.case = TRUE), ]
  }
  #We filter for duplicated data and just keep the first row to avoid problems with KDE. 
  data <- data[!duplicated(data[, c("decimalLongitude","decimalLatitude")]), ]
  #We turn the data table into an sf object, and again transfor it to metric units. 
  sf_points <- st_as_sf(data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) |>
    st_transform(32632)
  #Finally we just keep the points inside the Italy border with a logical vector. 
  #Sparse is false to get back a logical vector. 
  return(sf_points[st_intersects(sf_points, italy, sparse = FALSE), ])
}

# 3. Data Acquisition
#We apply the function and download Boar and Wolf data, excluding dogs and pigs. 
wolf_sf <- load_species_sf(5219173, "familiaris")
boar_sf <- load_species_sf(7705930, "domestic|familiaris")

# 4. Density Calculation (Sigma: 20km)
#We create our ppp objects necessary for the KDE. We extract X and Y coordinates in meters from the spatial object, 1 being the longitude and 2 the latitude.
#Then, we also define the observation window created before, our italy polygon. We do this process for wolf and boar.
wolf_ppp <- ppp(st_coordinates(wolf_sf)[,1], st_coordinates(wolf_sf)[,2], window = italy_poly)
boar_ppp <- ppp(st_coordinates(boar_sf)[,1], st_coordinates(boar_sf)[,2], window = italy_poly)
#With our ppp objects ready, we compute our density calculus. We use a common sigma of 20km for both species, with a 512x512 grid.
wolf_dens <- density(wolf_ppp, sigma = 20000, dimyx = 512)
boar_dens <- density(boar_ppp, sigma = 20000, dimyx = 512)

# 5. Log-Normalization
#We create the Log-Normalization function for our densities, to visualize and compare them better in the plots.
apply_log_norm <- function(dens_obj) {
  #We first create and add a small offset to add to every pixel, to avoid log(0) problem, so every pixel has a value.
  #na.rm ignores any NA values. 
  offset <- max(dens_obj$v, na.rm = TRUE) / 1000
  dens_obj$v <- log(dens_obj$v + offset)
  #We apply a Min-Max scaling for the normalization. Now every value is contained between 1.0 and 0.0 for both density scales. 
  dens_obj$v <- (dens_obj$v - min(dens_obj$v, na.rm=T)) / 
    (max(dens_obj$v, na.rm=T) - min(dens_obj$v, na.rm=T))
  return(dens_obj)
}
#Apply the Log-Normalization to our densities. 
wolf_dens_log <- apply_log_norm(wolf_dens)
boar_dens_log <- apply_log_norm(boar_dens)

# 6. Plotting Functions
#The first plotting function is used to plot the occurrence points in Italy.
plot_occ <- function(sf_points, species_label, color_p) {
  ggplot() +
    #We first draw the Italy map to put the points into.
    geom_sf(data = italy, fill = "#f8f9fa", color = "grey80", linewidth = 0.2) +
    #We draw the points, taking the data from our sf object created before.
    geom_sf(data = sf_points, color = color_p, size = 0.3, alpha = 0.4) +
    labs(title = paste("Occurrences:", species_label)) + 
    #We remove the default grey backgorund and grid lines of R, making the map look cleaner. 
    theme_minimal() + theme(panel.grid = element_blank())
}
#The second plotting function is used to plot the density estimations in Italy with a "heatmap". 
plot_dens <- function(dens_obj, species_label, palette) {
  #The density function creates an image, but we need a data table for ggplot2. We use the as.data.frame function to obtain it. 
  df <- as.data.frame(dens_obj)
  #We put names onto the columns. 
  colnames(df) <- c("x", "y", "value")
  
  ggplot() +
    #We first draw the shape of italy, filled with a light grey, acting as a background, so areas with zero density still look like part of the country.
    geom_sf(data = italy, fill = "#eeeeee", color = NA) +
    #We then draw our grid of pixels, the aes function tells R that the color of the pixel should be determined by the density number (0 to 1). Alpha is at 85% so you can still see the map.
    geom_raster(data = df, aes(x=x, y=y, fill=value), alpha = 0.85) +
    #We set our viridis color palette
    scale_fill_viridis(option = palette, name = "Log Density") +
    #We draw the italy border again, but this time with fill NA, just putting the outline on top of the heatmap so the borders look sharp. 
    geom_sf(data = italy, fill = NA, color = "white", linewidth = 0.1) +
    labs(title = paste("KDE:", species_label), subtitle = "Sigma: 20km") +
    #We remove the default grey backgorund and grid lines of R, making the map look cleaner. 
    theme_minimal() theme(panel.grid = element_blank())
}

# 7. Final Layout
#We create all 4 plots
p1 <- plot_occ(wolf_sf, "Wolf", "#35b779")
p2 <- plot_dens(wolf_dens_log, "Wolf", "viridis")
p3 <- plot_occ(boar_sf, "Boar", "#440154")
p4 <- plot_dens(boar_dens_log, "Boar", "magma")
#We join the created plots in one single image using patchwork. 
(p1 + p2) / (p3 + p4)

# 8. Statistical Analysis
#We compute our Spearman Analysis. We convert our density matrix into a single column vector, so each grid cell becomes a single observation. 
#We use complete.obs to ignore NA values. 
#We calculate the correlation coefficient between the two vectors with the Spearman method. 
spearman_rho <- cor(as.vector(wolf_dens_log$v), as.vector(boar_dens_log$v), 
                    method = "spearman", use = "complete.obs")
#We round our result to 2 decimal places, and then see the result. 
print(paste("Spearman Correlation:", round(spearman_rho, 2)))

# 9. Density Difference Map
#We create the Density Difference Map, showing areas where the Boar or the Wolf dominate. 
#We again turn the wolf density matrix into a data table, giving it specific names. 
diff_df <- as.data.frame(wolf_dens_log)
colnames(diff_df) <- c("x", "y", "wolf_val")
#Since both the Wolf and Boar density maps were created using the exact same grid, we can simply "paste" the Boar density values as a new column in our table. They line up perfectly pixel-for-pixel.
diff_df$boar_val <- as.data.frame(boar_dens_log)$value
#Since both variables are normalized, here is where we calculate the density. With this order, positive values will mean higher olf density, and negative will mean higher boar density.
diff_df$diff <- diff_df$wolf_val - diff_df$boar_val
#We use the same logic used for the normal density maps, just changing the data used (The differences). 
ggplot() +
  geom_sf(data = italy, fill = "grey95", color = NA) +
  geom_raster(data = diff_df, aes(x=x, y=y, fill=diff), alpha = 0.9) +
  scale_fill_viridis_c(option = "mako", name = "Difference") +
  geom_sf(data = italy, fill = NA, color = "white", linewidth = 0.1) +
  labs(title = "Spatial Dominance: Wolf vs Boar", 
       subtitle = "Brighter = Wolf dominance | Darker = Boar dominance") +
theme_minimal() + theme(panel.grid = element_blank())








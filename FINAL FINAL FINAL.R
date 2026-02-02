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
  if (!is.null(exclude_pattern)) {
    data <- data[!grepl(exclude_pattern, data$scientificName, ignore.case = TRUE), ]
  }
  #We filter for duplicated data and just keep the first row to avoid problems with KDE. 
  data <- data[!duplicated(data[, c("decimalLongitude","decimalLatitude")]), ]
  #We turn the data table into an sf object, and again transfor it to metric units. 
  sf_points <- st_as_sf(data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) |>
    st_transform(32632)
  #Finally we just keep the points inside the Italy border with a logical vector. 
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
apply_log_norm <- function(dens_obj) {
  offset <- max(dens_obj$v, na.rm = TRUE) / 1000
  dens_obj$v <- log(dens_obj$v + offset)
  dens_obj$v <- (dens_obj$v - min(dens_obj$v, na.rm=T)) / 
    (max(dens_obj$v, na.rm=T) - min(dens_obj$v, na.rm=T))
  return(dens_obj)
}

wolf_dens_log <- apply_log_norm(wolf_dens)
boar_dens_log <- apply_log_norm(boar_dens)

# 6. Plotting Functions
plot_occ <- function(sf_points, species_label, color_p) {
  ggplot() +
    geom_sf(data = italy, fill = "#f8f9fa", color = "grey80", linewidth = 0.2) +
    geom_sf(data = sf_points, color = color_p, size = 0.3, alpha = 0.4) +
    labs(title = paste("Occurrences:", species_label)) + 
    theme_minimal() + theme(panel.grid = element_blank())
}

plot_dens <- function(dens_obj, species_label, palette) {
  df <- as.data.frame(dens_obj)
  colnames(df) <- c("x", "y", "value")
  
  ggplot() +
    geom_sf(data = italy, fill = "#eeeeee", color = NA) +
    geom_raster(data = df, aes(x=x, y=y, fill=value), alpha = 0.85) +
    scale_fill_viridis(option = palette, name = "Log Density") +
    geom_sf(data = italy, fill = NA, color = "white", linewidth = 0.1) +
    labs(title = paste("KDE:", species_label), subtitle = "Sigma: 20km") +
    theme_minimal() + theme(panel.grid = element_blank())
}

# 7. Final Layout
p1 <- plot_occ(wolf_sf, "Wolf", "#35b779")
p2 <- plot_dens(wolf_dens_log, "Wolf", "viridis")
p3 <- plot_occ(boar_sf, "Boar", "#440154")
p4 <- plot_dens(boar_dens_log, "Boar", "magma")

(p1 + p2) / (p3 + p4)

# 8. Statistical Analysis (L-cross and Spearman)
multi_ppp <- superimpose(wolf = wolf_ppp, boar = boar_ppp)
ck <- Lcross(multi_ppp, "wolf", "boar", correction="border")
plot(ck, . - r ~ r, main="Spatial Interaction: Wolf vs Boar")

spearman_rho <- cor(as.vector(wolf_dens_log$v), as.vector(boar_dens_log$v), 
                    method = "spearman", use = "complete.obs")
print(paste("Spearman Correlation:", round(spearman_rho, 4)))

# 9. Density Difference Map
diff_df <- as.data.frame(wolf_dens_log)
colnames(diff_df) <- c("x", "y", "wolf_val")
diff_df$boar_val <- as.data.frame(boar_dens_log)$value
diff_df$diff <- diff_df$wolf_val - diff_df$boar_val

ggplot() +
  geom_sf(data = italy, fill = "grey95", color = NA) +
  geom_raster(data = diff_df, aes(x=x, y=y, fill=diff), alpha = 0.9) +
  scale_fill_viridis_c(option = "mako", name = "Difference") +
  geom_sf(data = italy, fill = NA, color = "white", linewidth = 0.1) +
  labs(title = "Spatial Dominance: Wolf vs Boar", 
       subtitle = "Brighter = Wolf dominance | Darker = Boar dominance") +

  theme_minimal() + theme(panel.grid = element_blank())


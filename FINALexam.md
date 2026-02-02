# Predator-Prey Spatial Interaction Analysis: Wolf (*Canis lupus*) vs. Wild Boar (*Sus scrofa*) in Italy

This project explores the spatial relationship between wolves (predator) and wild boars (prey) across the Italian territory using GBIF occurrence data and Spatial Point Pattern Analysis.

![WolfBoar](Lupo-e-Cinghiale.png)

## Research Question
Do high density areas for wolves coincide with high density areas for Boar? If so, At what distance ($r$) is the spatial attraction between the two species most significant?

# Data and Methodology
Data was downloaded using GBIF through the rgbif package. Domestic dogs and pigs were cleaned from the dataset to ensure data integrity. The Kernel Density estimation was performed using a sigma equal to 20km for both species, because it represents a realistic movement range for these large mammals. For the coordinate system, all data was projected to EPSG:32632 (UTM 32N) to allow for accurate distance measurements in meters. A log-transformation was applied to the density surfaces to handle the high variance in occurrence intensity and highlight subtle spatial trends. The analysis was performed entirely in R.

## Packages used
Here are the packages used in the project. 
- `rgbif` allows R to access directly to the GBIF servers to download occurrence recorods.
- `sf` treats geographic data (points, polygons) like a data frame, making it easy to crop, project and transform coordinates.
- `spatstat` was used for the Point Porcess objects (ppp) and calculating the KDE.
- `rnaturalearth` provided italy's borders used as the window of the analysis.
- `viridis` provided color scales designed to be read by everyone, including color blind people.
- `ggplot2` was used to build the maps and charts.
- `patchwork` was used to combine the obtained plots into a single image. 

## Study Area
Before loading the data, we must define the study area. The `rnaturalearth` package was used to get Italy's borders and project them into UTM Zone 32N (EPSG:32632). This ensures that all distance-based calculations (like the 20km sigma) are measured accurately in meters rather than degrees. This is better due to degrees not being a consistent unit of measurement, and could influence the analysis later on. 

Then the observation window was defined using the spatial data acquired. This will work the container for the `ppp`.

```R
#We load the italian map and immediately transform it to a metric system, better for calculus. 
italy <- ne_countries(country = "Italy", scale = "medium", returnclass = "sf") |>
  st_transform(32632)
#We transform the map into the observation window used later for the ppp
italy_poly <- as.owin(italy)
```
## Data Acquisition 
We retrieve occurrence data from GBIF. To ensure data integrity, we filter out domestic variants (e.g., "familiaris" or "domestic") and remove duplicate coordinates that would otherwise cause artificial density "spikes."

```R
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
#We apply the function and download Boar and Wolf data, excluding dogs and pigs. 
wolf_sf <- load_species_sf(5219173, "familiaris")
boar_sf <- load_species_sf(7705930, "domestic|familiaris")
```

## Kernel Density Estimation and Normalization
To analyze the spatial relationship between species, we must move from discrete "points" to a continuous "intensity surface." This process allows us to identify hotspots and areas of low activity across the entire Italian landscape.

We convert the coordinate data into a Point Pattern Object. This is required by the `spatstat` package that links the occurrence points to our defined geographic window (Italy).
```R
#We create our ppp objects necessary for the KDE. We extract X and Y coordinates in meters from the spatial object, 1 being the longitude and 2 the latitude.
#Then, we also define the observation window created before, our italy polygon. We do this process for wolf and boar.
wolf_ppp <- ppp(st_coordinates(wolf_sf)[,1], st_coordinates(wolf_sf)[,2], window = italy_poly)
boar_ppp <- ppp(st_coordinates(boar_sf)[,1], st_coordinates(boar_sf)[,2], window = italy_poly)
```
We then apply a KDE with a Sigma of 20km. This acts as a smoothing radius, also reflecting a realistic ecological scale for large mammals. `dimyx` is set to 51 to create a high resolution grid for the final maps. 
```R
#With our ppp objects ready, we compute our density calculus. We use a common sigma of 20km for both species, with a 512x512 grid.
wolf_dens <- density(wolf_ppp, sigma = 20000, dimyx = 512)
boar_dens <- density(boar_ppp, sigma = 20000, dimyx = 512)
```
Species occurrence data is often highly skewed, with a few areas having massive numbers of sightings while most have very few. To account for this and compare both datasets fairly, we created a custom function `apply_log_norm`, that performs a logarithmic transformation. This makes subtle patterns in lower-density areas more visible alongside high-density hotspots. This function also normalizes the data, scaling the values between 0 and 1. 

The normalization performed is a Min-Max normalization, performed using the following formula:
![Formula](normalization.png)
With this method, every number in the resulting matrix is now ranging from 0.0 to 1.0. 

```R
#We create the Log-Normalization function for our densities, to visualize and compare them better in the plots.
apply_log_norm <- function(dens_obj) {
  #We first create and add a small offset to add to every pixel, to avoid log(0) problem, so every pixel has a value.
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
```

## Plotting functions
Functions are used to guarantee that both Wolf and Boar maps have the exact same criteria used. It's also efficient, we can generate all four maps with lesser lines of code. 

### Occurence plots
The `plot_occ` function is used to plot the individual occurrence points, so raw GBIF data, over Italy. We set `size = 0.3` and `alpha = 0.4`. Using a low alpha (transparency) is crucial; it prevents "overplotting" where points stack on top of each other, allowing us to see where sightings are most densely clustered.
```R
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
```
### Density plots
The `plot_dens` function converts our mathematical KDE results into a visual heatmap. We convert the `spatstat` density object into a coordinate grid. `geom_raster` is used to draw the continuous color surface, covering the entire study area. High density areas are recognized thanks to the `scale_fill_viridis`. There are two `geom_sf` calls, the first is to draw the background and the second to draw a thin white border over the density colors, to clearly define the coastline of italy. 
```R
plot_dens <- function(dens_obj, species_label, palette) {
  #The density function creates an image, but we need a data table for ggplot2. 
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
```
## Final layout
Finally, we use the `patchwork` package, so we can see the plots side by side, ready to be compared instantly. The plots were placed into a single 2x2 grid. 
```R
p1 <- plot_occ(wolf_sf, "Wolf", "#35b779")
p2 <- plot_dens(wolf_dens_log, "Wolf", "viridis")
p3 <- plot_occ(boar_sf, "Boar", "#440154")
p4 <- plot_dens(boar_dens_log, "Boar", "magma")
(p1 + p2) / (p3 + p4)
```

## Statistical Analysis
Maps are useful for visual estimations of our data, but we need statistics to confirm our observations. We use two different methods: one for the raw points and one for the density surfaces.

### The Spearman Rank Correlation
We perform a pixel-by-pixel correlation between the two KDE surfaces. Spearman was chosen due to it's non-parametric nature. It looks at the rank of the density rather than the raw values, making it much better at handling the "clumpy" nature of wildlife data and any remaining outliers. A value closer to +1 indicates that as Boar density increases, Wolf density increases predictably. Inversly, values closer to -1 indicate the contrary. 0 means no correlation at all. 
```R
#We compute our Spearman Analysis. We convert our density matrix into a single column vector, so each grid cell becomes a single observation. 
#We use complete.obs to ignore NA values. 
#We calculate the correlation coefficient between the two vectors with the Spearman method. 
spearman_rho <- cor(as.vector(wolf_dens_log$v), as.vector(boar_dens_log$v), 
                    method = "spearman", use = "complete.obs")
#We round our result to 4 decimal places, and then see the result. 
print(paste("Spearman Correlation:", round(spearman_rho, 4)))
```

## The Density Difference Map
This maps tries to analyze in which parts of Italy is the wolf more established relative to its prey, and vice versa. We treat the two normalized surfaces as layers in a "spatial calculator." By subtracting the Boar values from the Wolf values. Positive values (Brighter colors) area ssociated with areas where the wolf's relative intensity is higher than the boar's. Negative values (darker colors) are a reas where the Boar's intensity is higher than the Wolf's.
```R
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
```

# Results + Discussion
## The Plot Grid
These were all the plots produced by the functions, placed on a single 2x2 grid. 

![Occurence and Density Plots](comparison_grid.png)

### Interpretation
Both species show some spatial clustering rathere than uniform coverage in the territory, in the Western Alps for wolves and boars and also the eastern-most tip of Friuli-Venezia Giulia for Boars. This is probably not the real ecological distribution of these species and there is sampling bias, more records where people are and where there are more monitoring projects. 

The Wolf shows high-density areas in the appennines and parts of the alps, with weaker but continuous density along the appennine. Wolves usually prefer low human density and good connectivity along its territory, which makes sense in this distribution. The Appennines acts as a long corridor along Itally, serving as a dispersal route for the species. 

The Po Plain and the south, on the contrary, are low density. This could be explained with more human presence and activity, limiting the distribution. The Isalnds are completely empty. The KDE suggests range expansion with strong reliance on mountainous corridors.

We could conclude that the Wolf's distribution looks like a recovering apex predator, that recently came back to the territory and whose distribution is more constrained by human pressure rather than climatic factors. 

The Boar shows a more broader and continuous coverage in the North and along the Appennine corridor. This makes sense due to the species being more of a generalist, with high reproductive rate and tollerating more human presence. Boars, differently from wolves, can thrive in fragmented habitats and mixed agro-forest landscapes. 

The Po Plain still shows some moderate density, unlike wolves. Boars could be exploiting crops and edge habitats in their favor. 

## Density Difference Map + Spearman 

![Density diff](difference_map.png)

### Interpretation

The analysis yielded a Spearman’s rank correlation coefficient ($\rho$) of 0.62. This indicates a statistically significant, strong positive correlation between the two species' density surfaces. It is also evidence of a numerical response, the mechanism where a predator population increases in density in response to an increase in prey density. 

The Density Difference Map reveals the fine-scale spatial balance between predator and prey. Brighter regions indicate areas of high relative Wolf establishment, particularly along the mountainous backbone of Italy. Conversely, darker regions highlight 'prey reservoirs'—high-density Boar areas that may represent future expansion fronts for the Wolf population. This spatial heterogeneity, combined with a positive Spearman correlation, illustrates a dynamic recolonization process where the predator follows the prey's footprint.

# Conclusion

Ultimately, the spatial analysis confirms that high density areas for both species coincide. The Spearman’s of 0.62 confirms that as Boar density increases, Wolf density follows a predictable upward trend. The spatial attraction between the two species is most significant at a distance of 20km to 60km. As shown by the L-cross function, this range represents the peak deviation from randomness, suggesting that the predator-prey link is strongest at the scale of regional pack territories and dispersal distances. The data proves that the Wolf’s recolonization of Italy is not a random expansion, but a statistically measurable "tracking" of its primary prey's footprint.















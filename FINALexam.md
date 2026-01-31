# Predator-Pray Spatial Interaction Analysis: Wolf (*Canis lupus*) vs. Wild Boar (*Sus scrofa*) in Italy

This project explores the spatial relationship between wolves (predator) and wild boars (prey) across the Italian territory using GBIF occurrence data and Spatial Point Pattern Analysis.

## Research Question
Do high density areas for wolves coincide with high density areas for Boar? If so, At what distance ($r$) is the spatial attraction between the two species most significant?

## Data and Methodology
Data was downloaded using GBIF through the rgbif package. Domestic dogs and pigs where cleaned from the dataset to ensure data integrity. The Kernel Density estimation was performed using a sigma equal to 20km for botch species, because it represents a realistic movement range for these large mammals. For the coordinate system, all data was projected to EPSG:32632 (UTM 32N) to allow for accurate distance measurements in meters. A log-transformation was applied to the density surfaces to handle the high variance in occurrence intensity and highlight subtle spatial trends. The analysis was performed entirely in R.

## Packages used
Here are the packages used in the project. 
- `rgbif` allows R to access directly to the GBIF servers to download occurrence recorods
- `sf` treats geographic data (points, polygons) like a data frame, making it easy to crop, project and transform coordinates.
- `spatstat` and `spatstat.explore` were used for the Point Porcess objects (ppp) and calculating the KDE and L-cross interaction.
- `rnaturalearth` provided italy's borders used as the window of the analysis.
- `viridis` provided color scales designed to be read by everyone, including color blind people.
- `ggplot2` was used to build the maps and charts
- `patchwork` was used to combine the obtained plots into a single image. 

## First Steps
Before loading the data, we must define the study area. The `rnaturalearth` package was used to get Italy's borders and project them into UTM Zone 32N (EPSG:32632). This ensures that all distance-based calculations (like the 20km sigma) are measured accurately in meters rather than degrees.
```{r loading italy}
italy <- ne_countries(country = "Italy", scale = "medium", returnclass = "sf") |>
  st_transform(32632)
italy_poly <- as.owin(italy)








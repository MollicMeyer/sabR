# 🧪 sabR: Example Workflow

This example demonstrates how to use the `sabR` package to:

- Load raster stacks of soil properties
- Compute zonal statistics from polygons
- Extract raster values at points
- Convert results into `SoilProfileCollection` objects

```r
# 📦 Load Required Libraries
library(sabR)
library(sf)
library(ggplot2)
library(terra)

# 📁 Define Inputs
raster_dir <- "T:/IA-Kitchen/results/maps/Batch1"
shapefile_path <- system.file("extdata", "SABRplots.shp", package = "sabR")
plots <- sf::st_read(shapefile_path)

# 🧮 Depth-weighted Zonal Statistics
results <- sabr.zonalstats(
  raster_path = raster_dir,
  tdepth = 0,
  bdepth = 20,
  props = c("sand", "clay", "TOC", "pH"),
  shapes = plots,
  plots = c("NW_plot", "SW_plot", "NE_plot", "SE_plot"),
  stats = c("mean", "sd", "min", "max"),
  output_dir = tempdir(),
  PSDnormalize = TRUE,
  MakePlot = FALSE,
  wtd.mean = TRUE
)

# 📊 Plot Weighted Means ± SD
ggplot(results$Weighted, aes(x = plot, y = weighted_mean, color = SoilProperty)) +
  geom_errorbar(
    aes(ymin = weighted_mean - weighted_sd, ymax = weighted_mean + weighted_sd),
    width = 0.2, position = position_dodge(0.6)
  ) +
  geom_point(position = position_dodge(0.6), size = 4) +
  facet_wrap(~SoilProperty, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    x = "Plot", y = "Weighted Mean ± SD",
    title = "Depth-weighted Zonal Stats by Plot and Soil Property",
    color = "Soil Property"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 🌐 Fetch SABR Raster Stack
sabr <- fetch_sabR(
  raster_path = raster_dir,
  props = c("sand", "clay", "TOC", "pH"),
  depths = c("0-5", "5-15", "15-30", "30-60", "60-100", "100-150", "150-200"),
  PSDnormalize = TRUE,
  aoi = plots,
  resample_method = "cubic",
  resample_res = 12
)

plot(sabr$stack)

# 🧱 Raster Stack → SoilProfileCollection
sabr_spc <- sabRas_to_spc(
  stack = sabr$stack,
  props = c("sand", "clay", "TOC", "pH"),
  depths = c("0-5", "5-15", "15-30", "30-60", "60-100", "100-150", "150-200"),
  new_depths = c(0, 20, 40, 60, 100, 200),
  source = sabr$product
)

# 📐 Polygons → SoilProfileCollection
zones_spc <- sabRzs_to_spc(
  stack = sabr$stack,
  aoi = plots,
  zones = c("NW_plot", "SW_plot", "NE_plot", "SE_plot"),
  id_column = "Name",
  props = c("sand", "clay", "TOC", "pH"),
  depths = c("0-5", "5-15", "15-30", "30-60", "60-100", "100-150", "150-200"),
  new_depths = c(0, 20, 40, 60, 100, 200),
  stats = c("mean", "sd"),
  source = "SABR"
)

# 📍 Random Points → SoilProfileCollection
set.seed(42)
aoi_sf <- if (!inherits(plots, "sf")) sf::st_as_sf(plots) else plots
rand_pts <- sf::st_sample(aoi_sf, size = 25, type = "random") %>% sf::st_as_sf()
rand_pts$Name <- paste0("pt_", seq_len(nrow(rand_pts)))

pts_spc <- sabRpts_to_spc(
  stack = sabr$stack,
  pts = rand_pts,
  props = c("sand", "clay", "TOC", "pH"),
  depths = c("0-5", "5-15", "15-30", "30-60", "60-100", "100-150", "150-200"),
  new_depths = c(0, 20, 40, 60, 100, 150),
  source = "SABR"
)
```

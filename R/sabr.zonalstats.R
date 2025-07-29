#' Zonal Statistics for SABR Raster Stacks
#'
#' Computes polygon-wise summary statistics for SABR-formatted raster data, including optional depth weighting.
#'
#' @param raster_path Character. Path to the directory containing raster subfolders.
#' @param tdepth Numeric. Top depth of interest in cm.
#' @param bdepth Numeric. Bottom depth of interest in cm.
#' @param props Character vector of soil property names to include.
#' @param shapes `sf` or `SpatVector` polygon object with plot geometry. Must include a `Name` field.
#' @param plots Character vector of plot IDs matching `shapes$Name`.
#' @param stats Character vector of statistics to compute (e.g., `"mean"`, `"sd"`, `"min"`, `"max"`).
#' @param output_dir Output folder for CSV results.
#' @param PSDnormalize Logical. If `TRUE`, use normalized PSD layers for sand/silt/clay.
#' @param MakePlot Logical. If `TRUE`, plots raster + overlay for each polygon.
#' @param wtd.mean Logical. If `TRUE`, computes depth-weighted summary across horizons.
#'
#' @return A list with:
#' \itemize{
#'   \item{Unweighted: Data frame of raw zonal statistics for all raster layers.}
#'   \item{Weighted: (if wtd.mean = TRUE) Depth-weighted summary statistics.}
#' }
#'
#' @export
sabr.zonalstats <- function(
  raster_path = maps_dir,
  tdepth = 0,
  bdepth = 20,
  props = c("sand", "clay"),
  shapes = plots,
  plots = c("NW_plot", "SW_plot"),
  stats = c("mean", "min", "max", "sd"),
  output_dir = output_dir,
  PSDnormalize = TRUE,
  MakePlot = TRUE,
  wtd.mean = TRUE
) {
  require(terra)
  require(sf)
  require(dplyr)
  require(tidyr)

  # Allowable property names
  allowed_attrs <- c(
    "BD",
    "buf_pH",
    "Ca",
    "CCE",
    "CEC",
    "cf",
    "clay",
    "Cstock",
    "EC",
    "K",
    "Mg",
    "Na",
    "P_Bray",
    "P_Olsen",
    "pH",
    "sand",
    "silt",
    "TOC",
    "Total_N"
  )
  props <- props[props %in% allowed_attrs]

  # Depth intervals lookup
  depth_intervals <- list(
    `0_5` = c(0, 5),
    `5_15` = c(5, 15),
    `15_30` = c(15, 30),
    `30_60` = c(30, 60),
    `60_100` = c(60, 100),
    `100_150` = c(100, 150),
    `150_200` = c(150, 200)
  )

  statdf <- data.frame()

  for (plot in plots) {
    poly <- shapes[shapes$Name == plot, ]
    poly <- terra::vect(st_geometry(poly))

    for (prop in props) {
      base_dir <- if (PSDnormalize && prop %in% c("sand", "silt", "clay")) {
        file.path(raster_path, "PSD_normalized")
      } else {
        file.path(raster_path, prop)
      }

      raster_files <- list.files(
        base_dir,
        pattern = paste0(prop, ".*\\.tif$"),
        full.names = TRUE,
        recursive = TRUE
      )

      for (di in names(depth_intervals)) {
        d_range <- depth_intervals[[di]]
        if (d_range[1] >= bdepth || d_range[2] <= tdepth) {
          next
        }

        file_match <- grep(di, raster_files, value = TRUE, fixed = TRUE)
        if (length(file_match) == 0) {
          next
        }

        stack <- rast(file_match[1])
        names(stack) <- tools::file_path_sans_ext(basename(file_match[1]))

        for (stat in stats) {
          zs <- terra::zonal(stack, poly, stat)[[1]]
          if (MakePlot) {
            plot(
              stack,
              col = terrain.colors(100),
              main = paste(plot, prop, di, stat, round(zs, 2))
            )
            plot(st_geometry(shapes), add = TRUE, border = "black")
            plot(poly, add = TRUE, border = "red", lwd = 2)
          }

          statdf <- rbind(
            statdf,
            data.frame(
              raster = names(stack),
              soilproperty = prop,
              raster_tdepth = d_range[1],
              raster_bdepth = d_range[2],
              tdepth = tdepth,
              bdepth = bdepth,
              plot = plot,
              statistic = stat,
              zonal_stats = zs
            )
          )
        }
      }
    }
  }

  result <- tidyr::spread(statdf, key = statistic, value = zonal_stats)
  result <- result %>%
    mutate(
      depth_range = raster_bdepth - raster_tdepth,
      overlap = pmin(bdepth, raster_bdepth) - pmax(tdepth, raster_tdepth),
      overlap = ifelse(overlap < 0, 0, overlap),
      weight = overlap / depth_range
    ) %>%
    arrange(soilproperty, raster_tdepth)

  # Write unweighted output
  f_out <- paste0(
    output_dir,
    "/ZonalStats_",
    paste0(props, collapse = "_"),
    "_",
    tdepth,
    "_",
    bdepth,
    ".csv"
  )
  write.csv(result, f_out, row.names = FALSE)

  if (wtd.mean) {
    weighted <- result %>%
      group_by(soilproperty, plot) %>%
      summarise(
        weighted_mean = sum(mean * weight, na.rm = TRUE) / sum(weight),
        max = max(mean, na.rm = TRUE),
        min = min(mean, na.rm = TRUE),
        sd = sd(mean, na.rm = TRUE),
        weighted_sd = sum(sd * weight, na.rm = TRUE) / sum(weight),
        .groups = "drop"
      )
    f_wtd <- paste0(
      output_dir,
      "/ZonalStatsWeighted_",
      paste0(props, collapse = "_"),
      "_",
      tdepth,
      "_",
      bdepth,
      ".csv"
    )
    write.csv(weighted, f_wtd, row.names = FALSE)
  } else {
    weighted <- NULL
  }

  return(list(Unweighted = result, Weighted = weighted))
}

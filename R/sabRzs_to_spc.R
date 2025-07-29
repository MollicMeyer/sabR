#' Convert Zonal Statistics from Polygons to SoilProfileCollection
#'
#' Computes zonal summary statistics for polygons and converts them to a list of SoilProfileCollections.
#'
#' @param stack SpatRaster; raster stack of SABR properties.
#' @param aoi Character path to shapefile, or an `sf`, `SpatVector`, or `SpatRaster` object.
#' @param zones Character vector; values in the `id_column` to subset from the AOI.
#' @param id_column Character; column name in AOI used for profile ID (e.g., "Name").
#' @param props Character vector; soil properties to include (e.g., c("sand", "clay")).
#' @param depths Character vector; depth intervals to include (e.g., c("0-5", "5-15")).
#' @param new_depths Numeric vector; target depth bins to aggregate into (e.g., c(0, 20, 40, 60)).
#' @param stats Character vector; statistics to compute via zonal stats (e.g., c("mean", "sd")).
#' @param source Character; prefix string for profile ID (e.g., "SABR").
#'
#' @return A named list of SoilProfileCollections (e.g., $mean, $sd).
#' @export
sabRzs_to_spc <- function(
  stack = sabr$stack,
  aoi = plots,
  zones = c("NW_plot", "SE_plot"),
  id_column = "Name",
  props = c("sand", "clay"),
  depths = c("0-5", "5-15"),
  new_depths = c(0, 15, 30),
  stats = c("mean", "sd"),
  source = "SABR"
) {
  require(terra)
  require(aqp)
  require(dplyr)
  require(tidyr)
  require(stringr)

  depth_lookup <- list(
    "0-5" = c(0, 5),
    "5-15" = c(5, 15),
    "15-30" = c(15, 30),
    "30-60" = c(30, 60),
    "60-100" = c(60, 100),
    "100-150" = c(100, 150),
    "150-200" = c(150, 200)
  )

  valid_depths <- names(depth_lookup)
  if (!all(depths %in% valid_depths)) {
    stop("All depths must be one of: ", paste(valid_depths, collapse = ", "))
  }

  wanted <- unlist(lapply(props, function(p) paste0(p, "_", depths)))
  stack <- stack[[wanted]]

  # Load or convert AOI
  if (is.character(aoi)) {
    aoi <- terra::vect(aoi)
  } else if (inherits(aoi, "sf")) {
    aoi <- terra::vect(aoi)
  } else if (inherits(aoi, "SpatRaster")) {
    aoi <- terra::as.polygons(aoi)
  } else if (!inherits(aoi, "SpatVector")) {
    stop("`aoi` must be a file path, sf object, SpatVector, or SpatRaster.")
  }

  # CRS align
  if (terra::crs(stack) != terra::crs(aoi)) {
    aoi <- terra::project(aoi, terra::crs(stack))
  }

  aoi_df <- as.data.frame(aoi)

  if (!id_column %in% names(aoi_df)) {
    stop("id_column not found in AOI.")
  }
  if (!all(zones %in% aoi_df[[id_column]])) {
    stop("Some requested zones not found in AOI.")
  }

  aoi_sub <- aoi[aoi_df[[id_column]] %in% zones, ]

  result_list <- list()

  for (stat in stats) {
    zonal_result <- terra::zonal(stack, aoi_sub, fun = stat)

    if (!"zone" %in% names(zonal_result)) {
      zonal_result$zone <- seq_len(nrow(zonal_result))
    }

    zone_map <- data.frame(
      zone = seq_len(nrow(aoi_sub)),
      id_val = as.character(aoi_sub[[id_column]])
    )

    zonal_result <- dplyr::left_join(zonal_result, zone_map, by = "zone")
    zonal_result$peiid <- paste0(source, "_", zonal_result$id_val)

    zonal_result <- zonal_result %>%
      select(-zone, -id_val)

    long_df <- zonal_result %>%
      pivot_longer(-peiid, names_to = "layer", values_to = "value") %>%
      separate(
        layer,
        into = c("variable", "depth"),
        sep = "_(?=\\d)",
        remove = FALSE
      ) %>%
      mutate(
        hzdept = sapply(depth, function(d) depth_lookup[[d]][1]),
        hzdepb = sapply(depth, function(d) depth_lookup[[d]][2])
      )

    wide_df <- long_df %>%
      select(peiid, hzdept, hzdepb, variable, value) %>%
      pivot_wider(names_from = variable, values_from = value)

    depths(wide_df) <- peiid ~ hzdept + hzdepb
    spc <- wide_df

    diced <- dice(
      spc,
      fm = as.formula(paste0(
        "1:",
        max(new_depths),
        " ~ ",
        paste(props, collapse = " + ")
      )),
      SPC = FALSE
    )

    #Remove slices outside new depth bins
    diced <- diced %>%
      filter(hzdept >= min(new_depths) & hzdept < max(new_depths))

    diced$depth_bin <- cut(
      diced$hzdept,
      breaks = new_depths,
      right = FALSE,
      labels = FALSE
    )

    bin_table <- data.frame(
      top = new_depths[-length(new_depths)],
      bottom = new_depths[-1],
      depth_bin = seq_along(new_depths[-1])
    )

    agg <- diced %>%
      left_join(bin_table, by = "depth_bin") %>%
      group_by(peiid, top, bottom) %>%
      summarise(
        across(all_of(props), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )

    depths(agg) <- peiid ~ top + bottom
    site(agg) <- data.frame(peiid = unique(agg$peiid))

    result_list[[stat]] <- agg
  }

  return(result_list)
}

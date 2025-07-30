#' Convert Zonal Statistics from Polygons to SoilProfileCollection
#'
#' Computes zonal summary statistics from raster stack over polygons and converts them to an `aqp::SoilProfileCollection`.
#'
#' @param stack SpatRaster; raster stack (e.g., from `fetch_sabR()`).
#' @param aoi Character path, sf object, SpatVector, or SpatRaster; area of interest (polygon zones or raster extent).
#' @param zones Character vector; values from `id_column` to extract (e.g., `"Zone1"`, `"Zone2"`).
#' @param id_column Character; column name in `aoi` used to identify zones (e.g., `"Name"`).
#' @param props Character vector; soil properties to include.
#' @param depths Character vector; raster depth intervals to include. Must be a subset of:
#' `"0-5"`, `"5-15"`, `"15-30"`, `"30-60"`, `"60-100"`, `"100-150"`, `"150-200"`.
#' @param new_depths Numeric vector; standard depths to aggregate to.
#' @param stats Character vector; statistics to compute (`"mean"`, `"sd"`, etc.).
#' @param source Character; prefix for profile IDs (e.g., `"SABR"`).
#'
#' @return A named list of `SoilProfileCollection`s, one for each statistic requested.
#'
#' @details
#' Raster layers are aggregated within each zone using specified statistics, then converted to profiles.
#' Profiles are sliced into 1 cm intervals and aggregated to `new_depths`.
#'
#' @import terra
#' @importFrom dplyr select mutate group_by summarise filter across left_join
#' @importFrom tidyr pivot_longer pivot_wider separate
#' @importFrom aqp depths site dice
#' @export
#'
sabRzs_to_spc <- function(
  stack = sabr$stack,
  aoi = NULL,
  zones = NULL,
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

  # Validate depths
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

  # Subset raster
  wanted <- unlist(lapply(props, function(p) paste0(p, "_", depths)))
  stack <- stack[[wanted]]

  # Load and coerce AOI
  if (is.character(aoi)) {
    if (!file.exists(aoi)) {
      stop("Shapefile not found at provided path.")
    }
    aoi <- terra::vect(aoi)
  } else if (inherits(aoi, "sf")) {
    aoi <- terra::vect(aoi)
  } else if (inherits(aoi, "SpatRaster")) {
    aoi <- terra::as.polygons(aoi)
  } else if (!inherits(aoi, "SpatVector")) {
    stop("aoi must be a path, sf, SpatVector, or SpatRaster.")
  }

  # Reproject AOI to match raster
  aoi <- terra::project(aoi, terra::crs(stack))
  aoi_df <- as.data.frame(aoi)

  if (!id_column %in% names(aoi_df)) {
    stop("id_column not found in AOI.")
  }

  if (!is.null(zones)) {
    if (!all(zones %in% aoi_df[[id_column]])) {
      stop("Some requested zones not found in AOI.")
    }
    aoi_sub <- aoi[aoi_df[[id_column]] %in% zones, ]
    zone_ids <- as.character(aoi_df[[id_column]][
      aoi_df[[id_column]] %in% zones
    ])
  } else {
    aoi_sub <- aoi
    zone_ids <- as.character(aoi_df[[id_column]])
  }

  result_list <- list()

  for (stat in stats) {
    zonal_result <- terra::zonal(stack, aoi_sub, fun = stat)

    # Validate and assign zone index
    if (!"zone" %in% names(zonal_result)) {
      zonal_result$zone <- seq_len(nrow(zonal_result))
    }

    # Ensure zone_ids is long enough
    if (length(zone_ids) < max(zonal_result$zone)) {
      stop(
        "Zone ID mapping mismatch. Check that aoi_sub matches zonal result order."
      )
    }

    # Construct peiid
    zonal_result$peiid <- paste0(source, "_", zone_ids[zonal_result$zone])

    # Long format and prep for SPC
    long_df <- zonal_result %>%
      select(-zone) %>%
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

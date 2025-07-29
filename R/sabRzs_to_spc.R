#' Convert Zonal Statistics from Polygons to SoilProfileCollection
#'
#' Computes zonal summary statistics for polygons and converts them to SoilProfileCollection.
#'
#' @param stack SpatRaster; the raster stack (e.g., stack = sabr$stack)
#' @param aoi SpatVector or sf object; polygons to extract zonal statistics from (e.g., aoi = plots)
#' @param zones Character vector; values in the `id_column` to extract (e.g., zones = c("Zone1", "Zone2"))
#' @param id_column Character; name of column in `aoi` used as profile ID (e.g., id_column = "Name")
#' @param props Character vector; soil properties to include (e.g., props = c("sand", "clay"))
#' @param depths Character vector; depth intervals (e.g., depths = c("0-5", "5-15"))
#' @param new_depths Numeric vector; new depth bins to aggregate to (e.g., new_depths = c(0, 15, 30))
#' @param stats Character vector; statistics to compute (e.g., stats = c("mean", "sd"))
#' @param source Character; string to prefix each profile ID (e.g., source = "SABR")
#'
#' @return A list of SoilProfileCollections, one per statistic
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

  # Filter raster stack to desired layers
  wanted <- unlist(lapply(props, function(p) paste0(p, "_", depths)))
  stack <- stack[[wanted]]

  # Ensure AOI is a SpatVector and filter to requested zones
  aoi <- terra::vect(aoi)
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

    # Join proper zone ID values
    zone_map <- data.frame(
      zone = zonal_result$zone,
      id_val = as.character(aoi_sub[[id_column]])
    )

    zonal_result <- zonal_result %>%
      left_join(zone_map, by = "zone") %>%
      mutate(peiid = paste0(source, "_", id_val)) %>%
      select(-zone, -id_val)

    # Long format + depth parsing
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

    # Wide format for SPC
    wide_df <- long_df %>%
      select(peiid, hzdept, hzdepb, variable, value) %>%
      pivot_wider(names_from = variable, values_from = value)

    depths(wide_df) <- peiid ~ hzdept + hzdepb
    spc <- wide_df

    # Dice to 1cm
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

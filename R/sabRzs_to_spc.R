#' Convert Zonal Statistics from Polygons to SoilProfileCollection
#'
#' Computes zonal summary statistics from a raster stack using input polygons (AOI)
#' and converts the results into SoilProfileCollection objects.
#'
#' @param stack SpatRaster; the raster stack (e.g., output from `fetch_sabR()$stack`)
#' @param aoi sf or SpatVector; polygon shapefile or object containing zones
#' @param zones Character vector; subset of zone IDs (e.g., c("NW_plot", "SE_plot")) to use from the `id_column`
#' @param id_column Character; name of column in `aoi` used as unique zone ID (e.g., "Name")
#' @param props Character vector; soil properties to extract (e.g., c("sand", "clay"))
#' @param depths Character vector; depth intervals to use (must match available raster names)
#' @param new_depths Numeric vector; new target depth bins (e.g., c(0, 15, 30))
#' @param stats Character vector; statistics to compute per zone (e.g., c("mean", "sd"))
#' @param source Character; prefix to append to profile ID (e.g., "SABR")
#'
#' @return A named list of SoilProfileCollection objects (e.g., $mean, $sd, etc.)
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

  # Depth label -> numeric top/bottom lookup
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

  # Filter raster stack to desired properties/depths
  wanted <- unlist(lapply(props, function(p) paste0(p, "_", depths)))
  stack <- stack[[wanted]]

  # Convert AOI to terra SpatVector and subset selected zones
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
    # Run zonal stats
    zonal_result <- terra::zonal(stack, aoi_sub, fun = stat)

    # Get IDs from zones (zone = row number)
    zone_ids <- as.data.frame(aoi_sub)[[id_column]]
    zonal_result$peiid <- paste0(source, "_", zone_ids[zonal_result$zone])

    # Reshape
    long_df <- zonal_result %>%
      pivot_longer(-c(zone, peiid), names_to = "layer", values_to = "value") %>%
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

    # Wide format
    wide_df <- long_df %>%
      select(peiid, hzdept, hzdepb, variable, value) %>%
      pivot_wider(names_from = variable, values_from = value)

    # Construct SPC
    depths(wide_df) <- peiid ~ hzdept + hzdepb
    spc <- wide_df

    # Dice to 1cm slices
    diced <- dice(
      spc,
      fm = as.formula(paste(
        "1:",
        max(new_depths),
        "~",
        paste(props, collapse = " + ")
      )),
      SPC = FALSE
    )

    # Bin and aggregate
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

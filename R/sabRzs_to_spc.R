#' Convert Zonal Statistics from Polygons to SoilProfileCollection
#'
#' Computes zonal summary statistics for polygons and converts them to SoilProfileCollection.
#'
#' @param stack SpatRaster; the raster stack (e.g., output from `fetch_sabR()$stack`)
#' @param zones SpatVector or sf object; polygons with ID field
#' @param id_column Character; name of column in `zones` used as unique profile ID (e.g., id_column = "Name")
#' @param props Character vector; soil properties to include (e.g., props = c("sand", "clay"))
#' @param depths Character vector; depth intervals (e.g., depths = c("0-5", "5-15"))
#' @param new_depths Numeric vector; new depth bins to aggregate to (e.g., c(0, 15, 30))
#' @param stats Character vector; statistics to compute (e.g., stats = c("mean", "sd"))
#' @param source Character; identifier to append to profile ID (e.g., source = "SABR")
#'
#' @return A list of SoilProfileCollections, one per statistic (e.g., $mean, $sd)
#' @export
sabRzs_to_spc <- function(
  stack = sabr$stack,
  zones = plots,
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

  lyrnames <- names(stack)
  wanted <- unlist(lapply(props, function(p) paste0(p, "_", depths)))
  stack <- stack[[wanted]]

  zones <- terra::vect(zones)

  result_list <- list()

  for (stat in stats) {
    # Perform zonal stats
    zonal_result <- terra::zonal(stack, zones, fun = stat)

    # Extract zone ID column from zones
    zone_ids <- as.data.frame(zones)[[id_column]]

    # Attach correct ID to result (terra uses a numeric 'zone' index)
    zonal_result$peiid <- paste0(source, "_", zone_ids[zonal_result$zone])
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
        hzdepb = sapply(depth, function(d) depth_lookup[[d]][2]),
        peiid = paste0(source, "_", peiid)
      )

    wide_df <- long_df %>%
      select(peiid, hzdept, hzdepb, variable, value) %>%
      pivot_wider(names_from = variable, values_from = value)

    depths(wide_df) <- peiid ~ hzdept + hzdepb
    spc <- wide_df

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

#' Extract Raster Values from Points to SoilProfileCollection
#'
#' Converts raster data at specified point locations into a SoilProfileCollection.
#'
#' @param stack SpatRaster; the raster stack to extract values from (e.g., output from `fetch_sabR()$stack`)
#' @param pts SpatVector or sf object; spatial points for extraction
#' @param props Character vector; soil properties to include (e.g., props = c("sand", "clay"))
#' @param depths Character vector; depth intervals to use (e.g., depths = c("0-5", "5-15"))
#' @param new_depths Numeric vector; new standard horizons (e.g., new_depths = c(0, 20, 40, 60))
#' @param source Character; used to append to unique ID (e.g., source = "SABR")
#'
#' @return A SoilProfileCollection with values from raster layers at point locations
#' @export
sabRpts_to_spc <- function(
  stack,
  pts,
  props = c("sand", "clay", "TOC", "pH"),
  depths = c("0-5", "5-15", "15-30"),
  new_depths = c(0, 15, 30),
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

  # Extract raster values at point locations
  vals <- terra::extract(stack, pts)

  # Create peiid from row number
  vals$peiid <- paste0(source, "_pt_", seq_len(nrow(vals)))

  # Convert to long and join with depth lookup
  long_df <- vals %>%
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

  # Convert to SPC and aggregate
  depths(wide_df) <- peiid ~ hzdept + hzdepb
  spc <- wide_df
  spc <- dice(
    spc,
    fm = as.formula(paste(
      "1:",
      max(new_depths),
      "~",
      paste(props, collapse = " + ")
    )),
    SPC = TRUE
  )

  binned <- dice(
    spc,
    fm = as.formula(paste(
      "1:",
      max(new_depths),
      "~",
      paste(props, collapse = " + ")
    )),
    SPC = FALSE
  )
  binned$depth_bin <- cut(
    binned$hzdept,
    breaks = new_depths,
    right = FALSE,
    labels = FALSE
  )

  bin_table <- data.frame(
    top = new_depths[-length(new_depths)],
    bottom = new_depths[-1],
    depth_bin = seq_along(new_depths[-1])
  )

  agg <- binned %>%
    left_join(bin_table, by = "depth_bin") %>%
    group_by(peiid, top, bottom) %>%
    summarise(across(all_of(props), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  depths(agg) <- peiid ~ top + bottom
  site(agg) <- data.frame(peiid = unique(agg$peiid))

  return(agg)
}

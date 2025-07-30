#' Extract Raster Values from Points and Convert to SoilProfileCollection
#'
#' Extracts raster data at specified point locations and converts them into an `aqp::SoilProfileCollection`.
#'
#' @param stack SpatRaster; raster stack to extract values from.
#' @param pts Character path, sf object, SpatVector, or SpatRaster; point locations for extraction.
#' @param props Character vector; soil properties to extract (e.g., `"sand"`, `"clay"`).
#' @param depths Character vector; raster depth intervals to include. Must be a subset of:
#' `"0-5"`, `"5-15"`, `"15-30"`, `"30-60"`, `"60-100"`, `"100-150"`, `"150-200"`.
#' @param new_depths Numeric vector; standard depths (in cm) to aggregate to, e.g., `c(0, 15, 30, 60, 100)`.
#' @param source Character; prefix for profile IDs (e.g., `"SABR"`).
#'
#' @return A `SoilProfileCollection` with one profile per point location.
#'
#' @details
#' Raster layers are matched by name to the specified properties and depths.
#' Extracted values are reshaped, converted to 1 cm slices, and aggregated by mean within `new_depths`.
#'
#' @import terra
#' @importFrom dplyr select mutate group_by summarise filter across left_join
#' @importFrom tidyr pivot_longer pivot_wider separate
#' @importFrom aqp depths site dice
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

  # Subset stack to desired layers
  wanted <- unlist(lapply(props, function(p) paste0(p, "_", depths)))
  stack <- stack[[wanted]]

  # Coerce pts to SpatVector
  if (is.character(pts)) {
    if (!file.exists(pts)) {
      stop("File does not exist at the provided path.")
    }
    pts <- terra::vect(pts)
  } else if (inherits(pts, "sf")) {
    pts <- terra::vect(pts)
  } else if (inherits(pts, "SpatRaster")) {
    pts <- terra::as.points(pts)
  } else if (!inherits(pts, "SpatVector")) {
    stop("pts must be a path, sf, SpatVector, or SpatRaster.")
  }

  # Reproject to raster CRS
  pts <- terra::project(pts, terra::crs(stack))

  # Extract values
  vals <- terra::extract(stack, pts)

  # Remove `.row` column if present
  vals <- vals[, !names(vals) %in% ".row", drop = FALSE]

  # Add ID
  vals$peiid <- paste0(source, "_pt_", seq_len(nrow(vals)))

  # Reshape and split variable names safely
  long_df <- vals %>%
    select(-ID) %>%
    pivot_longer(cols = -peiid, names_to = "layer", values_to = "value") %>%
    filter(!is.na(value))

  # Extract variable and depth using regex that respects internal underscores
  layer_info <- str_match(long_df$layer, "^(.*)_((\\d+-\\d+))$")

  # Assign variable and depth
  long_df$variable <- layer_info[, 2]
  long_df$depth <- layer_info[, 3]

  # Filter out any layers not in depth_lookup (if needed)
  long_df <- long_df %>%
    filter(depth %in% names(depth_lookup)) %>%
    mutate(
      hzdept = vapply(depth, function(d) depth_lookup[[d]][1], numeric(1)),
      hzdepb = vapply(depth, function(d) depth_lookup[[d]][2], numeric(1))
    )

  # Wide format
  wide_df <- long_df %>%
    select(peiid, hzdept, hzdepb, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value)

  # Create SPC
  depths(wide_df) <- peiid ~ hzdept + hzdepb
  spc <- wide_df

  # Slice to 1cm
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

  # Filter out NA-sliced layers
  diced <- diced %>%
    filter(hzdept >= min(new_depths) & hzdept < max(new_depths))

  # Bin to new depths
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

  # Aggregate
  agg <- diced %>%
    left_join(bin_table, by = "depth_bin") %>%
    group_by(peiid, top, bottom) %>%
    summarise(across(all_of(props), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  # Final SPC
  depths(agg) <- peiid ~ top + bottom
  site(agg) <- data.frame(peiid = unique(agg$peiid))

  return(agg)
}

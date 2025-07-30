#' Extract Raster Values from Points and Convert to SoilProfileCollection
#'
#' Extracts raster data at specified point locations and converts them into an `aqp::SoilProfileCollection`,
#' with optional assignment of group IDs from point attributes or spatial join with polygons.
#'
#' @param stack SpatRaster; raster stack to extract values from.
#' @param pts Character path, sf object, SpatVector, or SpatRaster; point locations for extraction.
#' @param props Character vector; soil properties to extract (e.g., `"sand"`, `"clay"`).
#' @param depths Character vector; raster depth intervals to include. Must be a subset of:
#' `"0-5"`, `"5-15"`, `"15-30"`, `"30-60"`, `"60-100"`, `"100-150"`, `"150-200"`.
#' @param new_depths Numeric vector; standard depths (in cm) to aggregate to.
#' @param source Character; prefix for profile IDs (e.g., `"SABR"`).
#' @param group_column Optional character; name of a column in `pts` to assign groups.
#' @param group_poly Optional polygon layer (sf, SpatVector, or path); spatial polygons for assigning groups via spatial join.
#' @param id_column Character; used with `group_poly` to assign group ID from that column.
#'
#' @return A `SoilProfileCollection` with site-level group metadata for later grouped depth function plots.
#'
#' @import terra
#' @importFrom dplyr select mutate group_by summarise filter across left_join
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_match
#' @importFrom aqp depths site dice
#' @export
sabRpts_to_spc <- function(
  stack,
  pts,
  props = c("sand", "clay", "TOC", "pH"),
  depths = c("0-5", "5-15", "15-30"),
  new_depths = c(0, 15, 30),
  source = "SABR",
  group_column = NULL,
  group_poly = NULL,
  id_column = "Name"
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

  # Coerce input
  if (is.character(pts)) {
    pts <- terra::vect(pts)
  } else if (inherits(pts, "sf")) {
    pts <- terra::vect(pts)
  } else if (inherits(pts, "SpatRaster")) {
    pts <- terra::as.points(pts)
  } else if (!inherits(pts, "SpatVector")) {
    stop("pts must be a path, sf, SpatVector, or SpatRaster.")
  }

  pts <- terra::project(pts, terra::crs(stack))

  # Optional: Join group from polygon
  if (!is.null(group_poly)) {
    if (is.character(group_poly)) {
      group_poly <- terra::vect(group_poly)
    } else if (inherits(group_poly, "sf")) {
      group_poly <- terra::vect(group_poly)
    } else if (!inherits(group_poly, "SpatVector")) {
      stop("group_poly must be a path, sf, or SpatVector.")
    }
    group_poly <- terra::project(group_poly, terra::crs(stack))
    join_result <- terra::extract(group_poly, pts, bind = TRUE)
    pts$group_id <- join_result[[id_column]]
  } else if (!is.null(group_column)) {
    if (!(group_column %in% names(pts))) {
      stop("group_column not found in pts.")
    }
    pts$group_id <- pts[[group_column]]
  } else {
    pts$group_id <- NA
  }

  # Extract raster values
  vals <- terra::extract(stack, pts)
  vals$peiid <- paste0(source, "_pt_", seq_len(nrow(vals)))
  vals$group_id <- pts$group_id

  # Reshape and parse variable/depth names
  long_df <- vals %>%
    select(-ID) %>%
    pivot_longer(
      cols = c(-peiid, -group_id),
      names_to = "layer",
      values_to = "value"
    ) %>%
    filter(!is.na(value))

  layer_info <- str_match(long_df$layer, "^(.*)_((\\d+-\\d+))$")
  long_df$variable <- layer_info[, 2]
  long_df$depth <- layer_info[, 3]

  long_df <- long_df %>%
    filter(depth %in% names(depth_lookup)) %>%
    mutate(
      hzdept = vapply(depth, function(d) depth_lookup[[d]][1], numeric(1)),
      hzdepb = vapply(depth, function(d) depth_lookup[[d]][2], numeric(1))
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
    summarise(across(all_of(props), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  depths(agg) <- peiid ~ top + bottom
  site(agg) <- data.frame(
    peiid = unique(agg$peiid),
    group_id = unique(vals$group_id)
  )

  return(agg)
}

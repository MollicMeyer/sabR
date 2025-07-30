#' Convert Raster Stack to SoilProfileCollection with Depth Aggregation
#'
#' Converts a stack of soil property rasters into an `aqp::SoilProfileCollection` object,
#' then aggregates values to user-defined standard depth intervals.
#'
#' @param stack SpatRaster; the raster stack where layer names follow the pattern `property_depth`, e.g., `sand_0-5`.
#' @param props Character vector; soil properties to include (e.g., `"sand"`, `"clay"`).
#' @param depths Character vector; raster depth intervals to include (e.g., `"0-5"`, `"5-15"`). Must be a subset of:
#' `"0-5"`, `"5-15"`, `"15-30"`, `"30-60"`, `"60-100"`, `"100-150"`, `"150-200"`.
#' @param new_depths Numeric vector; standard depths (in cm) to aggregate to, e.g., `c(0, 15, 30, 60, 100)`.
#' @param source Character; optional prefix for profile IDs (e.g., `"SABR"`).
#'
#' @return A `SoilProfileCollection` aggregated to the specified `new_depths`.
#'
#' @details
#' Each raster pixel becomes a profile. Spatial coordinates are preserved in the site data.
#' Raster values are diced into 1 cm slices and aggregated to the specified depth bins using mean values.
#'
#' @import terra
#' @importFrom dplyr select mutate group_by summarise filter across left_join
#' @importFrom tidyr pivot_longer pivot_wider separate
#' @importFrom aqp depths site dice
#' @export

sabRas_to_spc <- function(
  stack = rstack,
  props = c("sand", "clay", "silt", "TOC", "BD"),
  depths = c("0-5", "5-15", "15-30", "30-60", "60-100", "100-150", "150-200"),
  new_depths = c(0, 5, 15, 30, 60, 100, 150, 200),
  source = "SABR"
) {
  require(terra)
  require(dplyr)
  require(tidyr)
  require(aqp)
  require(stringr)

  # Allowed depth labels and attributes
  allowed_depths <- c(
    "0-5",
    "5-15",
    "15-30",
    "30-60",
    "60-100",
    "100-150",
    "150-200"
  )
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

  # Validate
  if (any(!props %in% allowed_attrs)) {
    stop(
      "Invalid soil property name(s): ",
      paste(props[!props %in% allowed_attrs], collapse = ", ")
    )
  }

  if (any(!depths %in% allowed_depths)) {
    stop(
      "All `depths` must be one of: ",
      paste(allowed_depths, collapse = ", ")
    )
  }

  # Lookup table for depths
  depth_dict <- list(
    "0-5" = c(0, 5),
    "5-15" = c(5, 15),
    "15-30" = c(15, 30),
    "30-60" = c(30, 60),
    "60-100" = c(60, 100),
    "100-150" = c(100, 150),
    "150-200" = c(150, 200)
  )

  # Convert raster to data.frame
  pt_df <- as.data.frame(stack, xy = TRUE, cells = TRUE, na.rm = TRUE)
  pt_df$peiid <- paste0("cell_", pt_df$cell, "_", source)

  # Site-level metadata
  site_data <- pt_df %>%
    select(peiid, x, y)

  # Convert to long format and safely split variable/depth
  long_df <- pt_df %>%
    select(-x, -y, -cell) %>%
    pivot_longer(cols = -peiid, names_to = "layer", values_to = "value")

  # Extract variable and depth label from the layer name
  layer_info <- str_match(long_df$layer, "^(.*)_((\\d+-\\d+))$")

  # Assign variable and depth_label
  long_df$variable <- layer_info[, 2]
  long_df$depth_label <- layer_info[, 3]

  # Filter and attach horizon depths
  long_df <- long_df %>%
    filter(variable %in% props, depth_label %in% depths) %>%
    mutate(
      hzdept = vapply(depth_label, function(d) depth_dict[[d]][1], numeric(1)),
      hzdepb = vapply(depth_label, function(d) depth_dict[[d]][2], numeric(1))
    )

  # Wide format
  hz_data <- long_df %>%
    select(peiid, hzdept, hzdepb, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value)

  # Convert to SPC
  depths(hz_data) <- peiid ~ hzdept + hzdepb
  site(hz_data) <- site_data

  # Dice to 1-cm increments up to max(new_depths)
  max_depth <- max(new_depths)
  vars <- intersect(props, horizonNames(hz_data))
  fm <- as.formula(paste0(
    "1:",
    max_depth,
    " ~ ",
    paste(vars, collapse = " + ")
  ))
  diced <- dice(hz_data, fm = fm, fill = TRUE, strict = FALSE, SPC = FALSE)

  # Depth bin assignment
  bin_labels <- paste0(new_depths[-length(new_depths)], "-", new_depths[-1])
  diced$depth_bin <- cut(
    diced$hzdept,
    breaks = new_depths,
    right = FALSE,
    labels = bin_labels
  )

  bin_df <- data.frame(
    depth_bin = bin_labels,
    top = new_depths[-length(new_depths)],
    bottom = new_depths[-1]
  )

  id_col <- idname(hz_data)

  # Aggregate by profile + bin
  agg_df <- diced %>%
    group_by(!!sym(id_col), depth_bin) %>%
    summarise(
      across(all_of(vars), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    left_join(bin_df, by = "depth_bin") %>%
    select(!!sym(id_col), top, bottom, all_of(vars)) %>%
    filter(!is.na(top) & !is.na(bottom)) %>%
    filter(rowSums(is.na(across(all_of(vars)))) < length(vars))

  # Final SPC
  agg_df <- as.data.frame(agg_df)
  depths(agg_df) <- as.formula(paste(id_col, "~ top + bottom"))
  site(agg_df) <- site(hz_data)

  return(agg_df)
}

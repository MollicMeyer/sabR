#' Fetch SABR Soil Property Rasters
#'
#' Loads a local directory of SABR soil raster data into a unified SpatRaster stack.
#' Supports optional normalization, clipping, and resampling.
#'
#' @param raster_path Path to directory containing raster folders (e.g., "T:/IA-Kitchen/results/maps").
#' @param props Character vector of soil properties to load (e.g., "sand", "clay").
#' @param PSDnormalize Logical; whether to use normalized PSD rasters (sand, silt, clay sum to 100%).
#' @param clip_shapefile Optional path to shapefile to clip rasters to.
#' @param resample_res Resolution to resample to (numeric); set to NULL to skip resampling.
#' @param resample_method Method used for resampling (e.g., "bilinear", "cubic").
#' @param output_dir Optional path to save combined raster stack (GeoTIFF).
#'
#' @return A list with elements:
#' \describe{
#'   \item{stack}{A `SpatRaster` of selected rasters.}
#'   \item{product}{The string `"SABR"` for consistency.}
#'   \item{file_paths}{Paths of included raster files.}
#' }
#'
#' @export
fetch_sabR <- function(
  raster_path = "T:/IA-Kitchen/results/maps/Batch1",
  props = c("sand", "clay", "silt", "TOC", "BD"),
  PSDnormalize = TRUE,
  clip_shapefile = NULL,
  resample_res = 12,
  resample_method = "bilinear",
  output_dir = NULL
) {
  require(terra)

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

  if (any(!props %in% allowed_attrs)) {
    warning(paste(
      "Warning: The following soil property map(s) don't exist:",
      paste("(", props[!props %in% allowed_attrs], ")", collapse = ", "),
      "Please check your spelling."
    ))
    props <- props[props %in% allowed_attrs]
  }

  depth_intervals <- list(
    `0_5` = "0-5",
    `5_15` = "5-15",
    `15_30` = "15-30",
    `30_60` = "30-60",
    `60_100` = "60-100",
    `100_150` = "100-150",
    `150_200` = "150-200"
  )

  all_layers <- list()
  file_paths <- list()

  for (prop in props) {
    base_dir <- if (PSDnormalize && prop %in% c("sand", "silt", "clay")) {
      file.path(raster_path, "PSD_normalized")
    } else {
      file.path(raster_path, prop)
    }

    raster_files <- list.files(
      base_dir,
      pattern = "\\.tif$",
      recursive = TRUE,
      full.names = TRUE
    )
    raster_files <- grep(prop, raster_files, value = TRUE, fixed = TRUE)

    for (depth_code in names(depth_intervals)) {
      depth_str <- depth_intervals[[depth_code]]
      depth_file <- grep(depth_code, raster_files, value = TRUE, fixed = TRUE)

      if (length(depth_file) > 0) {
        r <- rast(depth_file[1])
        names(r) <- paste0(prop, "_", depth_str)
        all_layers[[length(all_layers) + 1]] <- r
        file_paths[[length(file_paths) + 1]] <- depth_file[1]
      }
    }
  }

  if (length(all_layers) == 0) {
    stop("No raster layers were found for the requested properties.")
  }

  stack <- do.call(base::c, all_layers)

  if (!is.null(clip_shapefile)) {
    if (!file.exists(clip_shapefile)) {
      stop("Shapefile not found at provided path.")
    }
    shape <- terra::vect(clip_shapefile)
    shape <- terra::project(shape, terra::crs(stack))
    stack <- terra::mask(terra::crop(stack, shape), shape)
  }

  if (!is.null(resample_res)) {
    ext <- terra::ext(stack)
    crs_val <- terra::crs(stack)
    template <- terra::rast(
      extent = ext,
      resolution = resample_res,
      crs = crs_val
    )
    stack <- terra::resample(stack, template, method = resample_method)
  }

  attr(stack, "product") <- "SABR"

  if (!is.null(output_dir)) {
    writeRaster(
      stack,
      filename = file.path(output_dir, "sabR_stack.tif"),
      overwrite = TRUE
    )
  }

  return(list(
    stack = stack,
    product = "SABR",
    file_paths = unlist(file_paths)
  ))
}

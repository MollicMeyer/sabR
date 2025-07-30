#' Convert SoilProfileCollection to APSIM-X Soil Profiles
#'
#' This function converts a SoilProfileCollection object to a list of APSIM-X compatible soil profiles.
#' It supports estimation of hydraulic parameters using Saxton-Rawls or ROSETTA pedotransfer functions.
#'
#' @param spc A SoilProfileCollection object
#' @param KSAT.method Method to estimate saturated hydraulic conductivity: `"SaxtonRawls"` or `"ROSETTA"` in mm/day
#' @param BD.method Method to estimate bulk density: `"SaxtonRawls"`, `"map"` (use observed), or `"hybrid"` (choose observed if lower than SR estimate)
#' @param SAT.method Method to estimate saturation: `"SaxtonRawls"` or `"ROSETTA"`
#' @param carbon_conv Conversion factor for TOC to OM. Choose `"Van Bemmelen 1.724"` (default) or `"Pribyl 2.0"`
#' @param depthconv Unit conversion factor for depth. Default = 10 (e.g., cm to mm)
#' @param soil.bottom Bottom of soil profile in mm. Default = 180
#'
#' @return A named list of APSIM-X soil profiles with metadata (peiid, x, y, group_id if available)
#' @export
#'
#' @details
#' - **ROSETTA** refers to the pedotransfer model by Zhang & Schaap (2017) using sand, silt, clay, and optionally bulk density.
#' - **SaxtonRawls** refers to the model from Saxton & Rawls (2006) using sand, clay, and organic matter.
#' - **BD.method = "map"** uses existing BD values in the SPC; `"hybrid"` uses whichever is lower: observed BD or SaxtonRawls.
#' - Missing or incomplete data will cause the function to error.
#'
#' @references
#' - Saxton, K. E., & Rawls, W. J. (2006). Soil water characteristic estimates by texture and organic matter for hydrologic solutions. *Soil Sci. Soc. Am. J.*, 70(5), 1569–1578.
#' - Zhang, Y., & Schaap, M. G. (2017). Weighted recalibration of the Rosetta pedotransfer model with improved estimates of hydraulic parameter distributions and summary statistics. *J. Hydrology*, 547, 39–53.

spc_to_asp <- function(
  spc,
  KSAT.method = c("SaxtonRawls", "ROSETTA"),
  BD.method = c("SaxtonRawls", "map", "hybrid"),
  SAT.method = c("SaxtonRawls", "ROSETTA"),
  carbon_conv = "Van Bemmelen 1.724",
  depthconv = 10,
  soil.bottom = 180
) {
  KSAT.method <- match.arg(KSAT.method)
  BD.method <- match.arg(BD.method)
  SAT.method <- match.arg(SAT.method)

  stopifnot(inherits(spc, "SoilProfileCollection"))

  asp <- vector("list", length(spc))
  names(asp) <- profile_id(spc)

  for (i in seq_along(spc)) {
    prof <- spc[i, ]
    h <- horizons(prof)

    required_cols <- c("top", "bottom", "sand", "clay", "TOC")
    if (KSAT.method == "ROSETTA" || SAT.method == "ROSETTA") {
      required_cols <- unique(c(required_cols, "silt"))
    }
    if (BD.method != "SaxtonRawls") {
      required_cols <- unique(c(required_cols, "BD"))
    }

    missing_cols <- required_cols[!required_cols %in% names(h)]
    if (length(missing_cols) > 0) {
      stop(paste(
        "Missing required columns:",
        paste(missing_cols, collapse = ", ")
      ))
    }

    if (anyNA(h[, required_cols])) {
      stop(paste("Missing values detected in profile", profile_id(prof)))
    }

    h <- h %>%
      mutate(
        Thickness = (bottom - top) * depthconv,
        ParticleSizeSand = sand,
        ParticleSizeSilt = silt,
        ParticleSizeClay = clay,
        Carbon = TOC,
        PH = pH,
        AirDry = NA,
        LL15 = NA,
        DUL = NA,
        SAT = NA,
        KS = NA
      )

    for (j in seq_len(nrow(h))) {
      s <- as.numeric(h$ParticleSizeSand[j])
      c <- as.numeric(h$ParticleSizeClay[j])
      si <- as.numeric(h$ParticleSizeSilt[j])
      carbon <- as.numeric(h$Carbon[j])
      om <- carbon * carbon_conv

      h$LL15[j] <- SaxtonRawlsLL15(s, c, om)
      h$DUL[j] <- SaxtonRawlsDUL(s, c, om)
      h$AirDry[j] <- h$LL15[j] * ifelse(h$top[j] == 0, 0.5, 1)

      bd_rawls <- SaxtonRawlsBD(s, c, om)
      h$BD[j] <- switch(
        BD.method,
        "SaxtonRawls" = bd_rawls,
        "map" = h$BD[j],
        "hybrid" = ifelse(
          !is.na(h$BD[j]) && h$BD[j] < bd_rawls,
          h$BD[j],
          bd_rawls
        )
      )

      h$SAT[j] <- switch(
        SAT.method,
        "SaxtonRawls" = SaxtonRawlsSAT(s, c, om),
        "ROSETTA" = {
          ro <- try(
            ROSETTA(
              data.frame(sand = s, silt = si, clay = c, db = h$BD[j]),
              vars = c("sand", "silt", "clay", "db")
            ),
            silent = TRUE
          )
          if (
            !inherits(ro, "try-error") &&
              !is.null(ro$theta_s[1]) &&
              is.finite(ro$theta_s[1])
          ) {
            ro$theta_s[1]
          } else {
            warning(paste(
              "SAT invalid at profile",
              profile_id(prof),
              "layer",
              j
            ))
            NA
          }
        }
      )

      h$KS[j] <- switch(
        KSAT.method,
        "SaxtonRawls" = SaxtonRawlsKSAT(s, c, om),
        "ROSETTA" = {
          ro <- try(
            ROSETTA(
              data.frame(sand = s, silt = si, clay = c, db = h$BD[j]),
              vars = c("sand", "silt", "clay", "db")
            ),
            silent = TRUE
          )
          if (
            !inherits(ro, "try-error") &&
              !is.null(ro$ksat[1]) &&
              is.finite(ro$ksat[1])
          ) {
            10^ro$ksat[1] * 10 # mm/day
          } else {
            warning(paste(
              "KSAT invalid at profile",
              profile_id(prof),
              "layer",
              j
            ))
            NA
          }
        }
      )
    }

    asp[[i]] <- apsimx::apsimx_soil_profile(
      nlayers = nrow(h),
      Thickness = h$Thickness,
      BD = h$BD,
      AirDry = h$AirDry,
      LL15 = h$LL15,
      DUL = h$DUL,
      SAT = h$SAT,
      KS = h$KS,
      Carbon = h$Carbon,
      PH = h$PH,
      ParticleSizeClay = h$ParticleSizeClay,
      ParticleSizeSand = h$ParticleSizeSand,
      ParticleSizeSilt = h$ParticleSizeSilt,
      soil.bottom = soil.bottom
    )

    # Attach metadata
    site_data <- site(prof)
    asp[[i]]$metadata <- list(
      peiid = site_data$peiid[1],
      group_id = if ("group_id" %in% names(site_data)) {
        site_data$group_id[1]
      } else {
        NA
      },
      x = if ("x" %in% names(site_data)) site_data$x[1] else NA,
      y = if ("y" %in% names(site_data)) site_data$y[1] else NA
    )
  }

  return(asp)
}

#' Full Saxton & Rawls Soil Hydraulic Property Calculator
#'
#' Estimates BD, LL15, DUL, SAT, SWCON, and KSAT from sand, clay, and organic matter.
#'
#' @param pSand Sand content (%)
#' @param pClay Clay content (%)
#' @param pOM Organic matter content (%)
#'
#' @return Named list of BD, LL15, DUL, SAT, SWCON, and KSAT
#' @export
SaxtonRawls <- function(pSand = 50, pClay = 25, pOM = 4) {
  pSand <- pSand / 100
  pClay <- pClay / 100
  pOM <- pOM / 100

  theta_1500t <- -0.024 *
    pSand +
    0.487 * pClay +
    0.006 * pOM +
    0.005 * pSand * pOM -
    0.013 * pClay * pOM +
    0.068 * pSand * pClay +
    0.031
  LL15 <- theta_1500t + (0.14 * theta_1500t - 0.02)
  LL15 <- round(pmax(0.01, pmin(0.99, LL15)), 3)

  theta_33t <- -0.251 *
    pSand +
    0.195 * pClay +
    0.011 * pOM +
    0.006 * pSand * pOM -
    0.027 * pClay * pOM +
    0.452 * pSand * pClay +
    0.299
  DUL <- theta_33t + (1.283 * theta_33t^2 - 0.374 * theta_33t - 0.015)
  DUL <- round(pmax(0.01, pmin(0.99, DUL)), 3)

  theta_sat33t <- 0.278 *
    pSand +
    0.034 * pClay +
    0.022 * pOM -
    0.018 * pSand * pOM -
    0.027 * pClay * pOM -
    0.584 * pSand * pClay +
    0.078
  theta_sat33 <- theta_sat33t + (0.636 * theta_sat33t - 0.107)

  SAT <- DUL + theta_sat33 - 0.097 * pSand + 0.043
  SAT <- round(pmax(0.01, pmin(0.99, SAT)), 3)

  BD <- round(pmax(1.0, pmin(2.1, (1 - SAT) * 2.65)), 3)

  lambda <- (log(DUL) - log(LL15)) / (log(1500) - log(33))
  ksat <- 1930 * ((SAT - DUL)^(3 - lambda))
  SWCON <- round(0.15 + pmin(ksat, 75) / 100, 3)

  res <- list(
    BD = BD,
    LL15 = LL15 * 100,
    DUL = DUL * 100,
    SAT = SAT * 100,
    SWCON = SWCON * 100,
    KSAT = ksat * 24
  )
  return(res)
}

#' Estimate Bulk Density (BD) using Saxton & Rawls (2006)
#' @inheritParams SaxtonRawls
#' @return Bulk density (g/cm³)
#' @export
SaxtonRawlsBD <- function(pSand = 50, pClay = 25, pOM = 4) {
  SAT <- SaxtonRawlsSAT(pSand, pClay, pOM)
  BD <- round(pmax(1.0, pmin(2.1, (1 - SAT / 100) * 2.65)), 3)
  return(BD)
}

#' Estimate LL15 using Saxton & Rawls (2006)
#' @inheritParams SaxtonRawls
#' @return Lower limit (LL15) in %
#' @export
SaxtonRawlsLL15 <- function(pSand = 50, pClay = 25, pOM = 4) {
  res <- SaxtonRawls(pSand, pClay, pOM)
  return(res$LL15)
}

#' Estimate DUL using Saxton & Rawls (2006)
#' @inheritParams SaxtonRawls
#' @return Drained upper limit (DUL) in %
#' @export
SaxtonRawlsDUL <- function(pSand = 50, pClay = 25, pOM = 4) {
  res <- SaxtonRawls(pSand, pClay, pOM)
  return(res$DUL)
}

#' Estimate SAT using Saxton & Rawls (2006)
#' @inheritParams SaxtonRawls
#' @return Saturation (SAT) in %
#' @export
SaxtonRawlsSAT <- function(pSand = 50, pClay = 25, pOM = 4) {
  res <- SaxtonRawls(pSand, pClay, pOM)
  return(res$SAT)
}

#' Estimate SWCON using Saxton & Rawls (2006)
#' @inheritParams SaxtonRawls
#' @return SWCON (drainage coefficient) in %
#' @export
SaxtonRawlsSWCON <- function(pSand = 50, pClay = 25, pOM = 4) {
  res <- SaxtonRawls(pSand, pClay, pOM)
  return(res$SWCON)
}

#' Estimate Saturated Hydraulic Conductivity (KSAT) using Saxton & Rawls (2006)
#' @inheritParams SaxtonRawls
#' @return Saturated hydraulic conductivity (mm/day)
#' @export
SaxtonRawlsKSAT <- function(pSand = 50, pClay = 25, pOM = 4) {
  res <- SaxtonRawls(pSand, pClay, pOM)
  return(res$KSAT)
}

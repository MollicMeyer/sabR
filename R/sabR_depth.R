#' Plot Grouped Soil Depth Functions from a SoilProfileCollection
#'
#' Creates grouped depth function plots for selected soil properties using slab summaries.
#' Allows plotting of means ± SD, medians ± IQR, or medians ± percentiles by group.
#'
#' @param spc A `SoilProfileCollection` object.
#' @param variables Character vector of soil property names to plot (must match horizon names).
#' @param slab_structure Numeric vector of depth breaks to use for aggregation (e.g., `c(0, 20, 40, 60)`).
#' @param group_id Optional character name of the grouping column in `site(spc)` (e.g., `"plot"` or `"source"`). If `NULL`, all profiles are treated as one group.
#' @param stat Character; statistical summary type. One of:
#'   - `"mean_sd"`: mean ± standard deviation,
#'   - `"med_sd"`: median ± standard deviation,
#'   - `"med_qr"`: median ± interquartile range (25th–75th),
#'   - `"med_pr"`: median ± percentiles (set via `intervals`).
#' @param intervals Numeric vector of length 2; used when `stat = "med_pr"` or `"med_qr"` to set quantile boundaries (e.g., `c(0.1, 0.9)` or `c(0.25, 0.75)`).
#'
#' @return A `lattice` plot object (multi-panel depth function plot).
#'
#' @details
#' - Grouping variable (`group_id`) must be a factor. If not, it will be coerced.
#' - Limited to 8 groups for color palette compatibility.
#' - Internally uses `aqp::slab()` and `lattice::xyplot()`.
#'
#' @examples
#' \dontrun{
#' sabR_depth(
#'   spc = pts_spc_groupcol,
#'   variables = c("sand", "clay", "TOC", "pH"),
#'   slab_structure = c(0, 20, 40, 60),
#'   group_id = "group",
#'   stat = "med_pr",
#'   intervals = c(0.1, 0.9)
#' )
#' }
#'
#' @importFrom lattice xyplot
#' @importFrom aqp slab site horizonDepths horizonNames prepanel.depth_function panel.depth_function
#' @importFrom dplyr group_by summarise across ungroup
#' @importFrom tidyr pivot_longer
#' @export

sabR_depth <- function(
  spc,
  variables = c("sand", "clay", "TOC", "pH"),
  slab_structure = c(0, 20, 40, 60),
  group_id = NULL,
  stat = "mean_sd",
  intervals = c(5, 95)
) {
  if (!inherits(spc, "SoilProfileCollection")) {
    stop("`spc` must be a SoilProfileCollection object.")
  }

  # Validate inputs
  stat_options <- c("mean_sd", "med_sd", "med_qr", "med_pr")
  if (!stat %in% stat_options) {
    stop("Invalid `stat`. Choose from: ", paste(stat_options, collapse = ", "))
  }

  if (!is.null(group_id)) {
    if (!group_id %in% names(site(spc))) {
      stop(paste0("`", group_id, "` not found in site(spc)."))
    }
    site(spc)$source <- as.factor(site(spc)[[group_id]])
  } else {
    site(spc)$source <- factor("All Profiles")
  }

  n_groups <- length(unique(site(spc)$source))
  if (n_groups > 8) {
    stop("Too many groups: this function supports up to 8 groups for plotting.")
  }

  # Define statistical summaries
  slab_fun <- switch(
    stat,
    "mean_sd" = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      c(mean = m, lower = m - s, upper = m + s)
    },
    "med_sd" = function(x) {
      m <- median(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      c(median = m, lower = m - s, upper = m + s)
    },
    "med_qr" = function(x) {
      q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
      m <- median(x, na.rm = TRUE)
      c(median = m, lower = q[1], upper = q[2])
    },
    "med_pr" = function(x) {
      probs <- intervals / 100
      p <- quantile(x, probs = probs, na.rm = TRUE)
      m <- median(x, na.rm = TRUE)
      c(median = m, lower = p[1], upper = p[2])
    }
  )

  # Run slab
  slab_df <- slab(
    spc,
    fm = as.formula(paste0(
      "site(source) ~ ",
      paste(variables, collapse = " + ")
    )),
    slab.structure = slab_structure,
    slab.fun = slab_fun
  )

  slab_df$source <- factor(slab_df$source)

  # Set color palette
  okabe_ito <- c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#999999"
  )
  palette_colors <- okabe_ito[seq_len(n_groups)]

  # Force slab bottom to 0 for top slice
  slab_df$bottom[slab_df$top == 0] <- 0

  # Plot label logic
  label_map <- list(
    "mean_sd" = "Mean ± SD",
    "med_sd" = "Median ± SD",
    "med_qr" = "Median ± IQR",
    "med_pr" = paste0("Median ± P", intervals[1], "-", intervals[2])
  )
  y_label <- label_map[[stat]]

  # Generate plot
  xyplot(
    bottom ~ mean | variable,
    data = slab_df,
    lower = slab_df$lower,
    upper = slab_df$upper,
    groups = source,
    sync.colors = TRUE,
    alpha = 0.5,
    ylab = "Depth (cm)",
    xlab = y_label,
    ylim = c(max(slab_structure), 0),
    layout = c(length(variables), 1),
    scales = list(
      x = list(tick.number = 4, alternating = 2, relation = 'free'),
      y = list(tick.number = 6, alternating = 3, relation = 'free')
    ),
    par.settings = list(
      superpose.line = list(lwd = 3, col = palette_colors, lty = 1)
    ),
    panel = panel.depth_function,
    prepanel = prepanel.depth_function,
    strip = strip.custom(bg = grey(0.9)),
    auto.key = list(columns = 2, lines = TRUE, points = FALSE)
  )
}

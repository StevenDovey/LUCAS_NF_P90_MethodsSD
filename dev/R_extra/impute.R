#07.06.26 19:45 NZST
# Gap interpolation for missing diameter, used only where a stem genuinely
# existed but a measurement was missed. Interpolates between two bracketing
# measurements when both are available. Every fill is tagged with its source
# and bounded by the growth cut-offs (no negative growth, capped rate).

default_impute_config <- function() {
  list(dbh_growth_max_cm_yr = 5,   # cap on back/forward extrapolation
       allow_negative_growth = FALSE)
}

# Linear interpolation of one stem's DBH series over cycles (years), in place.
# `years` aligns cycles to time so the interpolation respects spacing.
interpolate_bracketed <- function(value, years, source, icfg) {
  n <- length(value)
  if (n < 3) return(list(value = value, source = source))
  for (i in 2:(n - 1)) {
    if (is.na(value[i])) {
      lo <- max(which(!is.na(value[1:(i - 1)])), -Inf)
      hi_idx <- which(!is.na(value)); hi_idx <- hi_idx[hi_idx > i]
      if (is.finite(lo) && length(hi_idx)) {
        hi <- hi_idx[1]
        frac <- (years[i] - years[lo]) / (years[hi] - years[lo])
        est  <- value[lo] + frac * (value[hi] - value[lo])
        if (!icfg$allow_negative_growth) est <- max(est, value[lo])
        value[i]  <- est
        source[i] <- "imputed_interp"
      }
    }
  }
  list(value = value, source = source)
}

# Apply bracketed interpolation per stem to a long table (sorted by cycle).
impute_dbh_gaps <- function(stems, icfg = default_impute_config()) {
  parts <- split(seq_len(nrow(stems)), paste(stems$plot, stems$tag))
  for (idx in parts) {
    g <- stems[idx, , drop = FALSE]
    o <- order(g$cycle)
    res <- interpolate_bracketed(g$dbh[o], g$year[o], g$dbh_source[o], icfg)
    stems$dbh[idx][o]        <- res$value
    stems$dbh_source[idx][o] <- res$source
  }
  stems
}

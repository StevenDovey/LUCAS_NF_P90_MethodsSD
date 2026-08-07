#08.06.26 #1
# Cross-cycle and reference checks, extending the per-row checks in validate.R. Adds a per-stem trajectory pass across the rounds, a height-to-diameter consistency check, a species reference check, and a plot-level stem density check. All are warnings, since each can have a legitimate cause and needs a person to review it, so none are quarantined. Input: long stem table (plot, tag, cycle, dbh, height, status, decay_class, species_code, year); `known_species` is the set of codes in the species table. Output: a flags data.frame on the validate.R schema.

deep_validate <- function(stems, known_species, cfg = default_validation_config()) {
  mk <- function(plot, tag, cycle, rule, detail) {
    if (length(plot) == 0) return(NULL)
    data.frame(plot = plot, tag = tag, cycle = cycle, rule = rule,
               severity = "warning", detail = detail, stringsAsFactors = FALSE)
  }
  out <- list()
  alive <- function(x) !is.na(x) & x == "A"

  # --- per-stem trajectory (light match-based pivot, no reshape) -------------
  key <- paste(stems$plot, stems$tag)
  uk  <- unique(key)
  pull <- function(col, k) { i <- stems$cycle == k; stems[[col]][i][match(uk, key[i])] }
  w  <- data.frame(plot = stems$plot[match(uk, key)], tag = stems$tag[match(uk, key)],
                   stringsAsFactors = FALSE)
  d1 <- pull("dbh",1); d2 <- pull("dbh",2); d3 <- pull("dbh",3)
  s1 <- pull("status",1); s2 <- pull("status",2); s3 <- pull("status",3)
  y1 <- pull("year",1); y2 <- pull("year",2); y3 <- pull("year",3)
  c1 <- pull("decay_class",1); c2 <- pull("decay_class",2); c3 <- pull("decay_class",3)

  pair <- function(da, db, sa, sb, ya, yb, iv) {
    dt <- yb - ya; both <- alive(sa) & alive(sb) & !is.na(da) & !is.na(db) & !is.na(dt) & dt > 0
    rate <- (db - da) / dt
    grow <- both & rate > cfg$dbh_growth_max_cm_yr
    shr  <- both & (db - da) < -cfg$dbh_shrink_tol_cm
    list(
      mk(w$plot[grow], w$tag[grow], NA, "dbh_growth_implausible",
         sprintf("diameter +%.2f cm/yr over %s", rate[grow], iv)),
      mk(w$plot[shr], w$tag[shr], NA, "dbh_shrink_while_alive",
         sprintf("diameter %.2f cm over %s while alive", (db - da)[shr], iv)))
  }
  out <- c(out, pair(d1,d2,s1,s2,y1,y2,"1-2"), pair(d2,d3,s2,s3,y2,y3,"2-3"))

  mono <- !is.na(d1) & !is.na(d2) & !is.na(d3) &
    sign(d2 - d1) * sign(d3 - d2) < 0 &
    abs(d2 - d1) > cfg$dbh_shrink_tol_cm & abs(d3 - d2) > cfg$dbh_shrink_tol_cm
  out <- c(out, list(mk(w$plot[mono], w$tag[mono], NA, "dbh_non_monotonic",
                        "diameter rises then falls (or the reverse) across rounds")))

  thr <- cfg$large_first_appearance_cm
  fb2 <- is.na(d1) & !is.na(d2) & d2 >= thr & alive(s2)
  fb3 <- is.na(d1) & is.na(d2) & !is.na(d3) & d3 >= thr & alive(s3)
  fb  <- fb2 | fb3
  out <- c(out, list(mk(w$plot[fb], w$tag[fb], NA, "large_first_appearance",
    sprintf("first recorded at %.0f cm; likely present and missed earlier", ifelse(fb2, d2, d3)[fb]))))

  dg <- (!alive(s1) & !alive(s2) & !is.na(d1) & !is.na(d2) & (d2 - d1) > cfg$dbh_shrink_tol_cm) |
        (!alive(s2) & !alive(s3) & !is.na(d2) & !is.na(d3) & (d3 - d2) > cfg$dbh_shrink_tol_cm)
  out <- c(out, list(mk(w$plot[dg], w$tag[dg], NA, "dead_stem_growing",
                        "standing dead diameter increases across rounds")))

  dd <- (!is.na(c1) & !is.na(c2) & c2 < c1) | (!is.na(c2) & !is.na(c3) & c3 < c2)
  out <- c(out, list(mk(w$plot[dd], w$tag[dd], NA, "decay_class_decreasing",
                        "decay class decreases over time (wood does not un-decay)")))

  # --- height to diameter outlier (long) ------------------------------------
  hd <- stems[!is.na(stems$height) & stems$height > 1.35 & !is.na(stems$dbh) & stems$dbh > 0, ]
  if (nrow(hd) > 50) {
    m <- lm(log(height - 1.35) ~ log(dbh), data = hd)
    r <- residuals(m); o <- abs(r) > cfg$height_resid_sd * sd(r)
    out <- c(out, list(mk(hd$plot[o], hd$tag[o], hd$cycle[o], "height_diameter_outlier",
      "height far from the height-diameter curve; may be a broken or snapped stem")))
  }

  # --- species reference completeness ---------------------------------------
  miss <- !is.na(stems$species_code) & nzchar(stems$species_code) & !(stems$species_code %in% known_species)
  out <- c(out, list(mk(stems$plot[miss], stems$tag[miss], stems$cycle[miss],
    "species_not_in_reference", "species code absent from the species table; default density used")))

  # --- plot-level live stem density outlier ---------------------------------
  liv <- stems$status == "A" & !is.na(stems$dbh) & stems$dbh >= cfg$dbh_min_live_cm
  cnt <- aggregate(liv, by = list(plot = stems$plot, cycle = stems$cycle), FUN = sum)
  z <- abs(cnt$x - median(cnt$x)) / mad(cnt$x)
  od <- is.finite(z) & z > cfg$density_mad_z
  out <- c(out, list(mk(cnt$plot[od], NA, cnt$cycle[od], "plot_density_outlier",
    sprintf("live stem count %d is a plot-level outlier", cnt$x[od]))))

  do.call(rbind, out)
}

# The rule names this stage can raise, for a coverage report that shows zeros.
deep_rule_names <- function() c(
  "dbh_growth_implausible", "dbh_shrink_while_alive", "dbh_non_monotonic",
  "large_first_appearance", "dead_stem_growing", "decay_class_decreasing",
  "height_diameter_outlier", "species_not_in_reference", "plot_density_outlier")

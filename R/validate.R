#08.06.26 #3
# Data-cleaning stage: flag records against a rule set and quarantine hard errors before the carbon calculation.
# Input: long stem table, one row per visited plot/tag/cycle.
# Output: list of
#   flags      - data.frame(plot, tag, cycle, rule, severity, detail)
#   quarantine - rows that hit an error (excluded from the carbon calc)
#   clean      - the rest (warnings and infos kept, still flagged)
# An error quarantines the row; a warning or info is kept. Thresholds live in `cfg` so they are tunable without editing code.

default_validation_config <- function() {
  list(
    dbh_max_cm            = 400,
    dbh_shrink_tol_cm     = 0.5,
    dbh_growth_max_cm_yr  = 5,
    height_max_m          = 60,
    height_resid_sd       = 4,
    mat_min_c             = -5,
    mat_max_c             = 20,
    large_first_appearance_cm = 10,
    density_mad_z         = 5,
    dbh_min_live_cm       = 2.5,
    status_allowed        = c("A", "X", "F", "P", "S", "")
  )
}

# Internal helper: append a flag row.
.add_flag <- function(flags, df, mask, rule, severity, detail) {
  mask[is.na(mask)] <- FALSE
  if (!any(mask)) return(flags)
  rbind(flags, data.frame(
    plot = df$plot[mask], tag = df$tag[mask], cycle = df$cycle[mask],
    rule = rule, severity = severity, detail = detail,
    stringsAsFactors = FALSE))
}

validate_stems <- function(stems, cfg = default_validation_config()) {
  flags <- data.frame(plot = character(), tag = character(), cycle = integer(),
                      rule = character(), severity = character(),
                      detail = character(), stringsAsFactors = FALSE)

  # --- Identity / structural ------------------------------------------------
  # A duplicate (plot,tag,cycle) is flagged on every copy, but only the most complete row (present diameter, then present species code) is kept; the rest are quarantined. The flag records which happened for each row.
  key <- paste(stems$plot, stems$tag, stems$cycle)
  is_dup <- duplicated(key) | duplicated(key, fromLast = TRUE)
  completeness <- (!is.na(stems$dbh)) + (!is.na(stems$species_code) & nzchar(stems$species_code))
  keep_dup <- !logical(nrow(stems))
  for (idx in split(which(is_dup), key[is_dup])) {
    keep_dup[idx] <- FALSE
    keep_dup[idx[which.max(completeness[idx])]] <- TRUE
  }
  dup_detail <- ifelse(keep_dup, "duplicate (plot,tag,cycle); the more complete row was kept",
                                  "duplicate (plot,tag,cycle); dropped, a more complete duplicate was kept")
  flags <- .add_flag(flags, stems, is_dup, "duplicate_stem_cycle", "error", dup_detail[is_dup])

  # A blank or unrecognised AliveState maps to no status (NA) rather than being dropped before validation, so it is caught here and quarantined like any other structural error, with the raw value kept for the record.
  has_raw_state <- "alive_state_raw" %in% names(stems)
  bad_status <- !(stems$status %in% cfg$status_allowed)
  status_detail <- if (has_raw_state)
    sprintf("status not in allowed set (raw AliveState: %s)",
            ifelse(is.na(stems$alive_state_raw) | !nzchar(stems$alive_state_raw),
                   "blank", stems$alive_state_raw))
    else rep("status not in allowed set", nrow(stems))
  flags <- .add_flag(flags, stems, bad_status, "status_code_invalid", "error", status_detail[bad_status])

  # AliveState "Unknown" is treated as dead (X) for the carbon calculation, the same as a confirmed death; flagged for review since Unknown is not itself evidence of death.
  if (has_raw_state) {
    flags <- .add_flag(flags, stems,
                       !is.na(stems$alive_state_raw) & stems$alive_state_raw == "Unknown",
                       "alive_state_unknown_treated_dead", "warning",
                       "AliveState recorded as Unknown; treated as dead (X) for the carbon calculation")
  }

  # --- DBH plausibility -----------------------------------------------------
  live <- stems$status == "A"
  stem_key <- paste(stems$plot, stems$tag)
  has_measured_dbh <- stem_key %in%
    names(which(tapply(!is.na(stems$dbh) & stems$dbh > 0, stem_key, any)))
  flags <- .add_flag(flags, stems, live & !is.na(stems$dbh) & stems$dbh <= 0,
                     "dbh_nonpositive_live", "error", "live stem with recorded DBH <= 0")
  flags <- .add_flag(flags, stems, live & is.na(stems$dbh) & !has_measured_dbh,
                     "dbh_missing_unfillable", "warning",
                     "live stem with no measured DBH at any cycle; carbon zero-filled")
  flags <- .add_flag(flags, stems, !is.na(stems$dbh) & stems$dbh > cfg$dbh_max_cm,
                     "dbh_above_max", "error",
                     sprintf("DBH > %g cm", cfg$dbh_max_cm))

  # --- Height plausibility --------------------------------------------------
  flags <- .add_flag(flags, stems,
                     live & !is.na(stems$height) & !is.na(stems$dbh) &
                       stems$height <= 1.35,
                     "height_below_breast", "error", "live DBH stem with height <= 1.35 m")
  flags <- .add_flag(flags, stems, !is.na(stems$height) & stems$height > cfg$height_max_m,
                     "height_above_max", "warning",
                     sprintf("height > %g m", cfg$height_max_m))

  # --- Dead-wood / CWD ------------------------------------------------------
  flags <- .add_flag(flags, stems,
                     !is.na(stems$decay_class) & !(stems$decay_class %in% 0:4),
                     "decay_class_invalid", "error", "decay class outside 0-4")
  if (all(c("led", "sed") %in% names(stems))) {
    flags <- .add_flag(flags, stems,
                       !is.na(stems$led) & !is.na(stems$sed) & stems$sed > stems$led,
                       "cwd_taper_inverted", "warning",
                       "CWD small-end > large-end (flagged, not swapped)")
  }

  # --- Spatial --------------------------------------------------------------
  if ("mat" %in% names(stems)) {
    flags <- .add_flag(flags, stems,
                       !is.na(stems$mat) & (stems$mat < cfg$mat_min_c | stems$mat > cfg$mat_max_c),
                       "mat_out_of_range", "warning", "MAT outside plausible NZ range")
  }

  # --- Cross-cycle life history (resurrection vs missed sample) -------------
  flags <- rbind(flags, .life_history_flags(stems))

  # --- Partition ------------------------------------------------------------
  # Any other error quarantines every row sharing its key, as before. A duplicate quarantines only the copies that were not kept.
  other_err_keys <- with(flags[flags$severity == "error" & flags$rule != "duplicate_stem_cycle", , drop = FALSE],
                         paste(plot, tag, cycle))
  is_err <- (key %in% other_err_keys) | (is_dup & !keep_dup)
  list(flags = flags,
       clean = stems[!is_err, , drop = FALSE],
       quarantine = stems[is_err, , drop = FALSE])
}

# Resurrection = dead then alive (warning). Missed sample = alive after a gap (info). Vectorised: sort by stem and cycle, then compare each row with the previous row of the same stem. No per-stem split or growing rbind.
.life_history_flags <- function(stems) {
  dead_codes <- c("X", "F", "P", "S")
  o   <- order(stems$plot, stems$tag, stems$cycle)
  st  <- stems[o, , drop = FALSE]
  n   <- nrow(st)
  key <- paste(st$plot, st$tag)
  same <- c(FALSE, key[-1] == key[-n])
  prev_status <- c(NA_character_, st$status[-n])
  prev_cycle  <- c(NA_integer_,   st$cycle[-n])
  res <- same & st$status == "A" & prev_status %in% dead_codes
  gap <- same & st$status == "A" & prev_status == "A" & (st$cycle - prev_cycle) > 1
  res[is.na(res)] <- FALSE; gap[is.na(gap)] <- FALSE

  empty <- data.frame(plot = character(), tag = character(), cycle = integer(),
                      rule = character(), severity = character(),
                      detail = character(), stringsAsFactors = FALSE)
  out <- rbind(
    if (any(res)) data.frame(plot = st$plot[res], tag = st$tag[res], cycle = st$cycle[res],
      rule = "resurrection_dead_to_alive", severity = "warning",
      detail = "dead at earlier cycle, alive later; back-propagated to alive (field sampling correction)",
      stringsAsFactors = FALSE) else empty,
    if (any(gap)) data.frame(plot = st$plot[gap], tag = st$tag[gap], cycle = NA_integer_,
      rule = "missed_sample_gap", severity = "info",
      detail = "alive after a cycle gap; interpolate", stringsAsFactors = FALSE) else empty)
  out
}

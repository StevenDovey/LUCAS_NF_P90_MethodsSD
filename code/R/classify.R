#08.06.26 #1
# Per-interval stem classification: growth, mortality, ingrowth, net change. Functions take aligned start/end vectors for one interval. "Effective DBH" is the measured field diameter for that cycle, or NA where there was no visit. NA counts as below threshold, since a stem absent at the start cannot have been above the size cut-off; this is the correct ingrowth gate.

# Net change in carbon over an interval (per stem).
carbon_change <- function(carbon_start, carbon_end) {
  alive_start <- !is.na(carbon_start) & carbon_start > 0
  alive_end   <- !is.na(carbon_end)   & carbon_end   > 0
  ifelse(alive_start & alive_end,  carbon_end - carbon_start,
  ifelse(alive_start & !alive_end, -carbon_start,
         NA_real_))
}

# Growth: a survivor already above the size threshold at both ends. Stems below threshold at the start are recruitment (see ingrowth), not growth, so the two partition the change with no overlap.
#   dbh_start_eff - effective start DBH (measured, else NA)
#   threshold_cm  - live size cut-off (default 2.5)
carbon_growth <- function(carbon_start, carbon_end,
                          dbh_start_eff = NULL, threshold_cm = 2.5) {
  above_start <- if (is.null(dbh_start_eff)) TRUE
                 else (!is.na(dbh_start_eff) & dbh_start_eff >= threshold_cm)
  ifelse(!is.na(carbon_start) & carbon_start > 0 &
         !is.na(carbon_end)   & carbon_end   > 0 & above_start,
         carbon_end - carbon_start, NA_real_)
}

# Net change per stem: end minus start, with an absent stem treated as zero. This is the true stock difference and does not depend on how recruits are split between growth and ingrowth.
carbon_net_change <- function(carbon_start, carbon_end) {
  s <- ifelse(is.na(carbon_start), 0, carbon_start)
  e <- ifelse(is.na(carbon_end),   0, carbon_end)
  d <- e - s
  ifelse(s == 0 & e == 0, NA_real_, d)
}

# Mortality: alive at start, gone/dead at end.
carbon_mortality <- function(carbon_start, carbon_end) {
  ifelse(!is.na(carbon_start) & carbon_start > 0 &
         (is.na(carbon_end) | carbon_end == 0),
         -carbon_start, NA_real_)
}

# Ingrowth: crossed the live size threshold during the interval.
#   visited_start    - TRUE if the start cycle was a real field visit
#   dbh_start_eff    - effective start DBH (measured, else NA)
#   dbh_end          - end-cycle DBH
#   threshold_cm     - live size cut-off (default 2.5)
carbon_ingrowth <- function(carbon_start, carbon_end, dbh_start_eff, dbh_end,
                            visited_start, threshold_cm = 2.5) {
  below_start <- is.na(dbh_start_eff) | dbh_start_eff < threshold_cm
  gate <- visited_start &
          !is.na(carbon_end) & carbon_end > 0 &
          !is.na(dbh_end) & dbh_end >= threshold_cm &
          below_start
  ifelse(gate, carbon_end - ifelse(is.na(carbon_start), 0, carbon_start),
         NA_real_)
}

# Convenience: all AGB fluxes for one interval as a data.frame.
classify_interval <- function(carbon_start, carbon_end,
                              dbh_start_eff, dbh_end, visited_start,
                              threshold_cm = 2.5) {
  data.frame(
    change    = carbon_net_change(carbon_start, carbon_end),
    growth    = carbon_growth(carbon_start, carbon_end, dbh_start_eff, threshold_cm),
    mortality = carbon_mortality(carbon_start, carbon_end),
    ingrowth  = carbon_ingrowth(carbon_start, carbon_end, dbh_start_eff,
                                dbh_end, visited_start, threshold_cm)
  )
}

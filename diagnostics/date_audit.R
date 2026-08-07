#08.06.26 #3
# Cross-checks the round date in plot_summary.csv, and the CWD file's own date, against the stems file (the baseline), per plot per round. Flags a missing date on either side; reports the day gap when both are present rather than judging a tolerance itself.
# Run from the repository root:  Rscript diagnostics/date_audit.R
root <- getwd()
rd <- function(name) read.csv(file.path(root, "mfedata", name), stringsAsFactors = FALSE,
                              check.names = FALSE, na.strings = c("NA", ""))

cycle_ranges <- list(
  `1` = as.Date(c("2002-01-01", "2007-12-31")),
  `2` = as.Date(c("2009-01-01", "2014-07-31")),
  `3` = as.Date(c("2014-08-01", "2024-12-31")))
assign_cycle <- function(obs_date) {
  out <- rep(NA_integer_, length(obs_date))
  for (cyc in names(cycle_ranges)) {
    r <- cycle_ranges[[cyc]]
    out[!is.na(obs_date) & obs_date >= r[1] & obs_date <= r[2]] <- as.integer(cyc)
  }
  out
}

# --- stems: the baseline, one representative date per plot per round -------
stm <- rd("MFESensitivity_Stems_20260209_reduced_dbl_removed3.csv")
stm_date  <- as.Date(stm$PlotObsStartDate, "%d/%m/%Y")
stm_cycle <- assign_cycle(stm_date)
stm_med   <- aggregate(stm_date, by = list(plot = stm$ParentPlotName, cycle = stm_cycle),
                       FUN = function(x) median(x))
names(stm_med)[3] <- "stems_date"
stm_med <- stm_med[!is.na(stm_med$cycle), ]

# --- plot_summary: one date per plot per round, already labelled by round. YYMMDD digit order (year = 2000 + d %/% 10000, month = (d %% 10000) %/% 100, day = d %% 100) confirmed against every non-empty value in the real file. --
ps <- rd("plot_summary.csv")
decode <- function(d) {
  d <- as.numeric(d)
  out <- as.Date(rep(NA_character_, length(d)))
  ok <- !is.na(d)
  out[ok] <- as.Date(sprintf("%04d-%02d-%02d", 2000 + d[ok] %/% 10000, (d[ok] %% 10000) %/% 100, d[ok] %% 100))
  out
}
ps_long <- do.call(rbind, lapply(1:3, function(k) data.frame(
  plot = ps$PLOT, cycle = k,
  ps_date = decode(ps[[sprintf("Date %s meas", c("1st","2nd","3rd")[k])]]),
  stringsAsFactors = FALSE)))

# --- CWD: same aggregation as stems -----------------------------------------
cw <- rd("MFESensitivity_CWD_20260209_reduced_dblV2.csv")
cw_date  <- as.Date(cw$PlotObsStartDate, "%d/%m/%Y")
cw_cycle <- assign_cycle(cw_date)
cw_med   <- aggregate(cw_date, by = list(plot = cw$ParentPlotName, cycle = cw_cycle),
                      FUN = function(x) median(x))
names(cw_med)[3] <- "cwd_date"
cw_med <- cw_med[!is.na(cw_med$cycle), ]

# --- every plot x round in play, from all three sources ---------------------
all_pc <- unique(rbind(stm_med[c("plot","cycle")], ps_long[c("plot","cycle")], cw_med[c("plot","cycle")]))
all_pc <- merge(all_pc, stm_med, by = c("plot","cycle"), all.x = TRUE)
all_pc <- merge(all_pc, ps_long, by = c("plot","cycle"), all.x = TRUE)
all_pc <- merge(all_pc, cw_med,  by = c("plot","cycle"), all.x = TRUE)

# base = stems (the baseline); other = plot_summary or CWD. missing_in_stems and both_missing are major; missing_in_other is minor, since the stems date is available to infill it; a date_gap row always reports the day count, no tolerance judged here.
flag_pair <- function(base, other, other_name) {
  both_na    <- is.na(base) & is.na(other)
  miss_base  <- is.na(base) & !is.na(other)
  miss_other <- !is.na(base) & is.na(other)
  gap        <- as.numeric(other - base)
  data.frame(
    plot = all_pc$plot, cycle = all_pc$cycle, comparison = other_name,
    stems_date = base, other_date = other,
    days_apart = ifelse(both_na | miss_base | miss_other, NA, gap),
    flag = ifelse(both_na, "both_missing",
            ifelse(miss_base, "missing_in_stems",
            ifelse(miss_other, "missing_in_other", "date_gap"))),
    severity = ifelse(both_na | miss_base, "major",
                ifelse(miss_other, "minor", "info")),
    stringsAsFactors = FALSE)
}

out <- rbind(
  flag_pair(all_pc$stems_date, all_pc$ps_date,  "plot_summary"),
  flag_pair(all_pc$stems_date, all_pc$cwd_date, "cwd"))

dc <- file.path(root, "output", "diagnostics")
dir.create(dc, recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(dc, "date_diagnostic.csv"), row.names = FALSE, na = "")
cat(sprintf("date_diagnostic.csv written : %d rows, %d major, %d minor\n",
            nrow(out), sum(out$severity == "major"), sum(out$severity == "minor")))

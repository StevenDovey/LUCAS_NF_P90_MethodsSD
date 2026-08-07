#08.06.26 #6
# Copies the raw stems file exactly, plus a derived cycle column and the diagnostic flag columns, so the file can be reviewed, the underlying data fixed, the added columns deleted, and the result put back as the original. Flags: the validation rules (error/warning/info), the mortality double-fire flag, and the date-check flag for that plot/round. A cross-round rule (cycle = NA in the validation logs, e.g. a diameter trajectory check) is attached to every cycle-row of that stem, not just one. The output is sorted by plot, tag, cycle so a duplicate lands on adjacent rows; the original row order is kept in the sort column, so sorting by that column restores it before the file replaces the original. Run after run.R, diagnostics/mortality_audit.R, and diagnostics/date_audit.R. Written to output/diagnostics/mfedata/, alongside the CWD and plot_summary working copies.
# Run from the repository root:  Rscript diagnostics/stems_annotated.R
root <- getwd()
rd <- function(name) read.csv(file.path(root, "mfedata", name), stringsAsFactors = FALSE,
                              check.names = FALSE, na.strings = c("NA", ""))
dc      <- file.path(root, "output", "diagnostics")
dc_mfe  <- file.path(dc, "mfedata")

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

stm <- rd("MFESensitivity_Stems_20260209_reduced_dbl_removed3.csv")
plot <- stm$ParentPlotName
tag  <- stm$ItemID
stm$cycle <- assign_cycle(as.Date(stm$PlotObsStartDate, "%d/%m/%Y"))
key      <- paste(plot, tag, stm$cycle)
stem_key <- paste(plot, tag)

# --- validation rules, by severity, cycle-specific and cross-round combined -
collapse <- function(rule, k) vapply(split(rule, k), function(x) paste(unique(x), collapse = ";"), character(1))
join_flags <- function(sv) {
  logs <- read.csv(file.path(dc, sprintf("validation_%ss.csv", sv)), stringsAsFactors = FALSE, check.names = FALSE)
  cyc_specific <- logs[!is.na(logs$cycle), ]
  cross_round  <- logs[is.na(logs$cycle), ]
  by_cycle <- collapse(cyc_specific$rule, paste(cyc_specific$plot, cyc_specific$tag, cyc_specific$cycle))
  by_stem  <- collapse(cross_round$rule,  paste(cross_round$plot, cross_round$tag))
  a <- unname(by_cycle[key]);      a[is.na(a)] <- ""
  b <- unname(by_stem[stem_key]);  b[is.na(b)] <- ""
  ifelse(nzchar(a) & nzchar(b), paste(a, b, sep = ";"), ifelse(nzchar(a), a, b))
}
stm$flags_error   <- join_flags("error")
stm$flags_warning <- join_flags("warning")
stm$flags_info    <- join_flags("info")

# --- mortality double-fire, by stem (plot,tag), no cycle --------------------
mort <- read.csv(file.path(dc, "mortality_diagnostic.csv"), stringsAsFactors = FALSE, check.names = FALSE)
dbl  <- mort[mort$section == "mortality_double_fire_detail", ]
stm$flags_mortality_double_fire <- ifelse(stem_key %in% paste(dbl$key, dbl$value), "mortality_double_fire", "")

# --- date check, by plot/cycle, major and minor only -------------------------
dat <- read.csv(file.path(dc, "date_diagnostic.csv"), stringsAsFactors = FALSE, check.names = FALSE)
dat <- dat[dat$severity %in% c("major", "minor"), ]
by_pc <- collapse(paste0(dat$comparison, ":", dat$flag), paste(dat$plot, dat$cycle))
stm$flags_date <- unname(by_pc[paste(plot, stm$cycle)])
stm$flags_date[is.na(stm$flags_date)] <- ""

# sort holds the original row position, added last so it is the final column;
# sort by it ascending to restore the original order. The rows themselves are
# written sorted by plot, tag, cycle so a duplicate_stem_cycle pair lands on
# adjacent rows.
stm$sort <- seq_len(nrow(stm))
stm <- stm[order(plot, tag, stm$cycle), ]

out_name <- "MFESensitivity_Stems_20260209_reduced_dbl_removed3.csv"
dir.create(dc_mfe, recursive = TRUE, showWarnings = FALSE)
write.csv(stm, file.path(dc_mfe, out_name), row.names = FALSE, na = "")
has_flag <- nzchar(stm$flags_error) | nzchar(stm$flags_warning) | nzchar(stm$flags_info) |
            nzchar(stm$flags_mortality_double_fire) | nzchar(stm$flags_date)
cat(sprintf("mfedata/%s written : %d rows, %d with at least one flag\n", out_name, nrow(stm), sum(has_flag)))

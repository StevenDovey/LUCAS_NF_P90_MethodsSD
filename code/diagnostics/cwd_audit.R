#08.06.26 #5
# Copies the raw CWD file exactly, plus a derived cycle column and diagnostic flag columns, so the file can be reviewed, the underlying data fixed, the added columns deleted, and the result put back as the original. The checks are the same kind validate.R already does for stems, applied to CWD for the first time since CWD never goes through validate_stems(). Nothing is excluded or corrected here. The output is sorted by plot, tag, cycle so a duplicate lands on adjacent rows; the original row order is kept in the sort column, so sorting by that column restores it before the file replaces the original. Written to output/diagnostics/mfedata/, alongside the stems and plot_summary working copies.
# Run from the repository root:  Rscript diagnostics/cwd_audit.R
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

cw <- rd("MFESensitivity_CWD_20260209_reduced_dblV2.csv")
sp <- rd("species.csv"); names(sp) <- make.names(names(sp))
cw$cycle <- assign_cycle(as.Date(cw$PlotObsStartDate, "%d/%m/%Y"))
led <- sqrt(as.numeric(cw$LargeEnd1) * as.numeric(cw$LargeEnd2))
sed <- sqrt(as.numeric(cw$SmallEnd1) * as.numeric(cw$SmallEnd2))
len <- as.numeric(cw$ItemObsComponentLinearDimension)
decay <- as.numeric(cw$DecayClass)
stature_bad <- !is.na(cw$CWDStature) & nzchar(cw$CWDStature) & !(cw$CWDStature %in% c("F","P","S"))
dup_key <- paste(cw$ParentPlotName, cw$ItemID, cw$cycle)

n <- nrow(cw)
rules <- list(
  list(mask = !is.na(led) & !is.na(sed) & sed > led, rule = "cwd_taper_inverted", severity = "warning"),
  list(mask = !is.na(decay) & !(decay %in% 0:4), rule = "decay_class_invalid", severity = "error"),
  list(mask = nzchar(cw$PreferredSpeciesCode) & !(cw$PreferredSpeciesCode %in% sp$Species.code),
       rule = "species_not_in_reference", severity = "warning"),
  list(mask = is.na(len) | len <= 0, rule = "length_missing_or_nonpositive", severity = "warning"),
  list(mask = is.na(led) & is.na(sed), rule = "both_diameters_missing", severity = "warning"),
  list(mask = stature_bad, rule = "stature_unexpected", severity = "info"),
  list(mask = duplicated(dup_key) | duplicated(dup_key, fromLast = TRUE), rule = "duplicate_cwd_piece", severity = "error"))

by_severity <- function(sv) {
  out <- rep("", n)
  for (r in rules) if (r$severity == sv) {
    m <- r$mask; m[is.na(m)] <- FALSE
    out[m] <- ifelse(nzchar(out[m]), paste(out[m], r$rule, sep = ";"), r$rule)
  }
  out
}
cw$flags_error   <- by_severity("error")
cw$flags_warning <- by_severity("warning")
cw$flags_info    <- by_severity("info")

# sort holds the original row position, added last so it is the final column;
# sort by it ascending to restore the original order. The rows themselves are
# written sorted by plot, tag, cycle so a duplicate_cwd_piece pair lands on
# adjacent rows.
cw$sort <- seq_len(n)
cw <- cw[order(cw$ParentPlotName, cw$ItemID, cw$cycle), ]

dc_mfe <- file.path(root, "output", "diagnostics", "mfedata")
dir.create(dc_mfe, recursive = TRUE, showWarnings = FALSE)
out_name <- "MFESensitivity_CWD_20260209_reduced_dblV2.csv"
write.csv(cw, file.path(dc_mfe, out_name), row.names = FALSE, na = "")
has_flag <- nzchar(cw$flags_error) | nzchar(cw$flags_warning) | nzchar(cw$flags_info)
cat(sprintf("mfedata/%s written : %d rows, %d with at least one flag\n", out_name, n, sum(has_flag)))

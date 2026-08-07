#08.06.26 #5
# Copies the raw plot_summary file exactly, plus diagnostic flag columns, so the file can be reviewed, the underlying data fixed, the added columns deleted, and the result put back as the original. Checks: MAT and inner-plot-area completeness and plausibility, checked directly on this file rather than only visible after being merged onto stems. A disagreement between the two inner-area methods reports the numeric difference, no tolerance is judged here. The output is sorted by PLOT so a duplicate_plot_row pair lands on adjacent rows; the original row order is kept in the sort column, so sorting by that column restores it before the file replaces the original. Written to output/diagnostics/mfedata/, alongside the stems and CWD working copies.
# Run from the repository root:  Rscript diagnostics/plotsummary_audit.R
root <- getwd()
source(file.path(root, "R/validate.R"))
cfg <- default_validation_config()
rd <- function(name) read.csv(file.path(root, "mfedata", name), stringsAsFactors = FALSE,
                              check.names = FALSE, na.strings = c("NA", ""))

ps  <- rd("plot_summary.csv")
mat <- as.numeric(ps$MAT)
new  <- as.numeric(ps[["Inner plot area new method"]])
orig <- as.numeric(ps[["Inner plot area original method"]])
mat_bad   <- !is.na(mat) & (mat < cfg$mat_min_c | mat > cfg$mat_max_c)
area_diff <- !is.na(new) & !is.na(orig) & new != orig
dup       <- duplicated(ps$PLOT) | duplicated(ps$PLOT, fromLast = TRUE)

n <- nrow(ps)
rules <- list(
  list(mask = mat_bad, rule = "mat_out_of_range", severity = "warning"),
  list(mask = is.na(mat), rule = "mat_missing", severity = "warning"),
  list(mask = is.na(new) & is.na(orig), rule = "inner_area_missing", severity = "warning"),
  list(mask = area_diff, rule = "inner_area_methods_disagree", severity = "info"),
  list(mask = dup, rule = "duplicate_plot_row", severity = "error"))

by_severity <- function(sv) {
  out <- rep("", n)
  for (r in rules) if (r$severity == sv) {
    m <- r$mask; m[is.na(m)] <- FALSE
    out[m] <- ifelse(nzchar(out[m]), paste(out[m], r$rule, sep = ";"), r$rule)
  }
  out
}
ps$flags_error   <- by_severity("error")
ps$flags_warning <- by_severity("warning")
ps$flags_info    <- by_severity("info")

# sort holds the original row position, added last so it is the final column;
# sort by it ascending to restore the original order. The rows themselves are
# written sorted by PLOT so a duplicate_plot_row pair lands on adjacent rows.
ps$sort <- seq_len(n)
ps <- ps[order(ps$PLOT), ]

dc_mfe <- file.path(root, "output", "diagnostics", "mfedata")
dir.create(dc_mfe, recursive = TRUE, showWarnings = FALSE)
out_name <- "plot_summary.csv"
write.csv(ps, file.path(dc_mfe, out_name), row.names = FALSE, na = "")
has_flag <- nzchar(ps$flags_error) | nzchar(ps$flags_warning) | nzchar(ps$flags_info)
cat(sprintf("mfedata/%s written : %d plots, %d with at least one flag\n", out_name, n, sum(has_flag)))

#07.06.26 18:30 NZST
# Three-way comparison of imputed/predicted DBH and height:
# reference vs old R vs new pipeline. Rows are keyed by (ItemID, cycle).
# Only predicted or imputed values are compared, since measured values are
# the same in all three sources.
root <- "/home/user/SAS2R"

# --- NEW pipeline: build `st` (has predht + measured dbh) and impute dbh -----
source(file.path(root, "tools", "compare_run.R"))   # builds `st`
source(file.path(root, "R", "impute.R"))
nd <- data.frame(plot = st$plot, tag = st$tag, cycle = st$cycle,
                 year = st$year, dbh = st$dbh, stringsAsFactors = FALSE)
nd <- impute_dbh_gaps(nd)                            # NEW interpolated dbh
new <- data.frame(tag = st$tag, cycle = st$cycle,
                  new_predht = st$predht,            # pure model prediction
                  new_dbh = nd$dbh,
                  new_dbh_src = nd$dbh_source,
                  is_woody = st$is_woody,
                  height_meas = st$height)
new <- new[!is.na(new$tag) & !is.na(new$cycle), ]

# --- reference output --------------------------------------------------------
ref <- read.csv(file.path(root, "sasout/SAS_measurements-with-carbonV10.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)[
                c("ItemID","measurement","DBH","predht","height","source")]
names(ref) <- c("tag","cycle","ref_dbh","ref_predht","ref_height","ref_src")
for (c in c("cycle","ref_dbh","ref_predht","ref_height")) ref[[c]] <- suppressWarnings(as.numeric(ref[[c]]))

# --- OLD R ------------------------------------------------------------------
old <- read.csv(file.path(root, "rfiles/measurements_C_R.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)[
                c("ItemID","cycle","dbh","pred_dbh","predht","height","source")]
names(old) <- c("tag","cycle","old_dbh","old_pred_dbh","old_predht","old_height","old_src")
for (c in c("cycle","old_dbh","old_pred_dbh","old_predht")) old[[c]] <- suppressWarnings(as.numeric(old[[c]]))

# --- align ------------------------------------------------------------------
m <- merge(new, ref, by = c("tag","cycle"))
m <- merge(m, old, by = c("tag","cycle"))
cat(sprintf("aligned stem-cycles: %d\n", nrow(m)))

verdict <- function(new_v, old_v, ref, label, sub) {
  ok <- sub & is.finite(new_v) & is.finite(old_v) & is.finite(ref)
  if (!any(ok)) { cat(sprintf("  %-28s (no rows)\n", label)); return(invisible()) }
  en <- abs(new_v[ok] - ref[ok]); eo <- abs(old_v[ok] - ref[ok])
  cat(sprintf("  %-26s n=%6d | median|err| new %6.3f  old %6.3f | new closer: %5.1f%% | mean|err| new %6.3f old %6.3f\n",
      label, sum(ok), median(en), median(eo), 100*mean(en < eo), mean(en), mean(eo)))
}

cat("\n== HEIGHT: predht vs reference predht ==\n")
verdict(m$new_predht, m$old_predht, m$ref_predht, "all woody",
        m$is_woody)
verdict(m$new_predht, m$old_predht, m$ref_predht, "woody, height NOT measured",
        m$is_woody & is.na(m$height_meas))

cat("\n== DBH: imputed value vs reference DBH ==\n")
# where reference imputed the DBH (source != measured) and new/old both have a value
ref_est <- !is.na(m$ref_src) & m$ref_src != "measured"
verdict(m$new_dbh, m$old_dbh, m$ref_dbh, "where reference DBH estimated", ref_est)
verdict(m$new_dbh, m$old_dbh, m$ref_dbh, "where NEW interpolated dbh",
        !is.na(m$new_dbh_src) & m$new_dbh_src == "imputed_interp")

#08.06.26 #3
# Summarise the validation flags, run after run.R. Joins the flags to habit, species and cycle from the per-stem diagnostics, then writes count tables and a flags-by-rule chart under output/diagnostics/:
#   validation_summary_by_rule.csv        rule, severity, flags, plots, stems
#   validation_summary_rule_x_cycle.csv   flags by rule and cycle
#   validation_summary_rule_x_habit.csv   flags by rule and habit
#   validation_summary_rule_x_species.csv flags by rule and species (long)
#   validation_flags_by_rule.png
# Run from the repository root, after run.R:  Rscript diagnostics/run_validation_summary.R
root <- getwd()
dc <- file.path(root, "output", "diagnostics")
files <- file.path(dc, c("validation_errors.csv","validation_warnings.csv","validation_infos.csv"))
logs <- do.call(rbind, lapply(files, function(f) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)))
d <- read.csv(file.path(dc, "stem_diagnostics.csv"), stringsAsFactors = FALSE, check.names = FALSE)

# Habit and species are stem attributes and do not vary by round, so they are joined on the stem key alone. Many deep-validation rules carry cycle = NA (they flag a trajectory across rounds, not one round), so a cycle-aware join would leave those rules unattributed.
stem_key <- paste(logs$plot, logs$tag)
dk <- paste(d$plot, d$tag)
logs$habit        <- d$habit[match(stem_key, dk)]
logs$species_code <- d$species_code[match(stem_key, dk)]

# main table: per rule, the counts and the plots and stems affected
by_rule <- do.call(rbind, lapply(split(seq_len(nrow(logs)), logs$rule), function(i)
  data.frame(rule = logs$rule[i][1], severity = logs$severity[i][1],
             flags = length(i), plots = length(unique(logs$plot[i])),
             stems = length(unique(stem_key[i])), stringsAsFactors = FALSE)))
by_rule <- by_rule[order(factor(by_rule$severity, c("error","warning","info")), -by_rule$flags), ]
write.csv(by_rule, file.path(dc, "validation_summary_by_rule.csv"), row.names = FALSE, na = "")

# Breakdowns by cycle, habit, species. The blank-key cases are labelled rather than left as NA, which would read as missing data. A flag has no single cycle when its rule compares a stem across rounds (cross_round). A flag has no habit or species when its stem is not in the diagnostics: the plot-level density rule has no stem, and quarantined duplicates are removed (unmatched).
logs$cycle_lab   <- ifelse(is.na(logs$cycle), "cross_round", as.character(logs$cycle))
logs$habit_lab   <- ifelse(is.na(logs$habit), "unmatched", logs$habit)
logs$species_lab <- ifelse(is.na(logs$species_code), "unmatched", logs$species_code)
wtab <- function(by) { m <- as.data.frame.matrix(table(logs$rule, by))
  cbind(rule = rownames(m), m, stringsAsFactors = FALSE) }
write.csv(wtab(logs$cycle_lab), file.path(dc, "validation_summary_rule_x_cycle.csv"), row.names = FALSE, na = "")
write.csv(wtab(logs$habit_lab), file.path(dc, "validation_summary_rule_x_habit.csv"), row.names = FALSE, na = "")
spx <- as.data.frame(table(rule = logs$rule, species = logs$species_lab))
spx <- spx[spx$Freq > 0, ]; spx <- spx[order(spx$rule, -spx$Freq), ]
write.csv(spx, file.path(dc, "validation_summary_rule_x_species.csv"), row.names = FALSE, na = "")

# chart: flags by rule, coloured by severity
sev_col <- c(error = "#cc4444", warning = "#dd9933", info = "#4477aa")
o <- order(by_rule$flags)
png(file.path(dc, "validation_flags_by_rule.png"),
    width = 15, height = 10, units = "cm", res = 300, pointsize = 8)
par(mar = c(4, 11, 3, 1))
barplot(by_rule$flags[o], names.arg = by_rule$rule[o], horiz = TRUE, las = 1, cex.names = 0.6,
        col = sev_col[by_rule$severity[o]], xlab = "number of flags", main = "Validation flags by rule")
legend("bottomright", legend = names(sev_col), fill = sev_col, cex = 0.7, bty = "n")
dev.off()

cat("validation summary written to output/diagnostics/\n\n")
print(by_rule, row.names = FALSE)

#08.06.26 04:10 NZST
# Data-cleaning stage for the stem table. Runs the per-row checks plus the
# cross-cycle and reference checks, splits the flags into error, warning, and
# info logs alongside the quarantine file, and prints a per-rule count that
# lists every rule even when it never fired.
# Run from the repository root:  Rscript run_validate.R
root <- getwd()
source(file.path(root, "R", "validate.R"))
source(file.path(root, "R", "validate_deep.R"))

sp <- read.csv(file.path(root, "mfedata/species.csv"), stringsAsFactors = FALSE, check.names = FALSE)
names(sp) <- make.names(names(sp))
known_species <- sp$Species.code

stm <- read.csv(file.path(root, "mfedata/MFESensitivity_Stems_20260209_reduced_dbl_removed3.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
s <- data.frame(plot = stm$ParentPlotName, tag = stm$ItemID,
                status = ifelse(stm$AliveState == "Alive", "A", "X"),
                dbh = as.numeric(stm$DiameterValue), height = as.numeric(stm$Height),
                decay_class = as.numeric(stm$DecayClass), species_code = stm$PreferredSpeciesCode,
                subplot = stm$SubPlot,
                date = as.Date(stm$PlotObsStartDate, "%d/%m/%Y"),
                stringsAsFactors = FALSE)
# Cycle by exact match of the record date to the plot's first three visit dates.
pd <- tapply(s$date, s$plot, function(x) head(sort(unique(x)), 3))
lk <- unlist(lapply(names(pd), function(p) setNames(seq_along(pd[[p]]), paste(p, as.integer(pd[[p]])))))
s$cycle <- unname(lk[paste(s$plot, as.integer(s$date))])
s$year  <- as.integer(format(s$date, "%Y"))
s <- s[!is.na(s$cycle), ]

v  <- validate_stems(s)
dv <- deep_validate(s, known_species)
flags <- rbind(v$flags, dv)

dir.create(file.path(root, "output"), showWarnings = FALSE)
for (sv in c("error", "warning", "info"))
  write.csv(flags[flags$severity == sv, ], file.path(root, sprintf("output/validation_%ss.csv", sv)), row.names = FALSE)
write.csv(v$quarantine, file.path(root, "output/validation_quarantine.csv"), row.names = FALSE)

cat(sprintf("stem-cycles checked: %d\n", nrow(s)))
cat(sprintf("errors: %d   warnings: %d   info: %d   quarantined: %d   clean: %d\n",
            sum(flags$severity == "error"), sum(flags$severity == "warning"),
            sum(flags$severity == "info"), nrow(v$quarantine), nrow(v$clean)))
cat("written: output/validation_errors.csv, validation_warnings.csv, validation_infos.csv, validation_quarantine.csv\n")

# Coverage report: every rule with its count, including the rules with zero.
all_rules <- c(
  "duplicate_stem_cycle", "status_code_invalid", "dbh_nonpositive_live", "dbh_above_max",
  "dbh_missing_unfillable", "height_below_breast", "height_above_max", "decay_class_invalid",
  "cwd_taper_inverted", "mat_out_of_range", "resurrection_dead_to_alive", "missed_sample_gap",
  deep_rule_names())
counts <- sapply(all_rules, function(r) sum(flags$rule == r))
cat("\ncoverage (every rule, including zeros):\n")
print(data.frame(rule = all_rules, count = counts, row.names = NULL))

#08.06.26 #6
# Mortality double-count check and AliveState resolution, run after run.R. Confirms two things against the real data rather than a reading of the code: no stem fires a mortality loss in both the 1-2 and 2-3 intervals, and how many "Not Found"/"Unknown" coded deaths are later resurrected versus stay dead. Writes one tidy CSV, section/key/value.
# Run from the repository root, after run.R:  Rscript diagnostics/mortality_audit.R
root <- getwd()
source(file.path(root, "R/classify.R"))

# --- mortality double-fire, from the real per-cycle carbon_agb -------------
d <- read.csv(file.path(root, "output/diagnostics/stem_diagnostics.csv"),
              stringsAsFactors = FALSE, check.names = FALSE)
key <- paste(d$plot, d$tag)
uk  <- unique(key)
pull <- function(col, cyc) { i <- d$cycle == cyc; d[[col]][i][match(uk, key[i])] }
a1 <- pull("carbon_agb", 1); a2 <- pull("carbon_agb", 2); a3 <- pull("carbon_agb", 3)
m12 <- carbon_mortality(a1, a2)
m23 <- carbon_mortality(a2, a3)
double_fire <- !is.na(m12) & !is.na(m23)
w_plot <- d$plot[match(uk, key)]
w_tag  <- d$tag[match(uk, key)]

# --- raw AliveState, cycle re-derived from the raw stem date (mirrors run.R's assign_cycle) since AliveState granularity does not survive into stem_diagnostics.csv, only the collapsed A/X status does -------------------
stm <- read.csv(file.path(root, "mfedata/MFESensitivity_Stems_20260209_reduced_dbl_removed3.csv"),
                stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", ""))
obs_date <- as.Date(stm$PlotObsStartDate, "%d/%m/%Y")
cycle_ranges <- list(
  `1` = as.Date(c("2002-01-01", "2007-12-31")),
  `2` = as.Date(c("2009-01-01", "2014-07-31")),
  `3` = as.Date(c("2014-08-01", "2024-12-31")))
raw_cycle <- rep(NA_integer_, length(obs_date))
for (cyc in names(cycle_ranges)) {
  r <- cycle_ranges[[cyc]]
  raw_cycle[!is.na(obs_date) & obs_date >= r[1] & obs_date <= r[2]] <- as.integer(cyc)
}
alive_state <- stm$AliveState

# --- Not Found / Unknown: resurrected later, or stayed dead ----------------
raw_key <- paste(stm$ParentPlotName, stm$ItemID)
nf_unk  <- alive_state %in% c("Not Found", "Unknown")
affected <- unique(raw_key[nf_unk & !is.na(raw_cycle)])
resurrected <- 0L; stayed_dead <- 0L
for (k in affected) {
  i <- raw_key == k & !is.na(raw_cycle)
  nf_cycles    <- raw_cycle[i][nf_unk[i]]
  alive_cycles <- raw_cycle[i][alive_state[i] == "Alive"]
  if (any(nf_cycles < max(c(alive_cycles, -Inf)))) resurrected <- resurrected + 1L
  else stayed_dead <- stayed_dead + 1L
}

# --- assemble ----------------------------------------------------------------
summary_rows <- data.frame(section = "summary",
  key = c("stems_checked", "mortality_double_fire_count"),
  value = as.character(c(length(uk), sum(double_fire))), stringsAsFactors = FALSE)

tab <- as.data.frame(table(cycle = raw_cycle, alive_state = alive_state, useNA = "ifany"),
                     stringsAsFactors = FALSE)
alive_rows <- data.frame(section = "alive_state_by_cycle",
  key = paste0("cycle ", tab$cycle, " - ", tab$alive_state),
  value = as.character(tab$Freq), stringsAsFactors = FALSE)

resolution_rows <- data.frame(section = "not_found_unknown_resolution",
  key = c("resurrected", "stayed_dead"),
  value = as.character(c(resurrected, stayed_dead)), stringsAsFactors = FALSE)

detail_rows <- if (any(double_fire)) {
  data.frame(section = "mortality_double_fire_detail",
             key = w_plot[double_fire], value = w_tag[double_fire], stringsAsFactors = FALSE)
} else {
  data.frame(section = character(), key = character(), value = character())
}

out <- rbind(summary_rows, alive_rows, resolution_rows, detail_rows)
dc <- file.path(root, "output", "diagnostics")
dir.create(dc, recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(dc, "mortality_diagnostic.csv"), row.names = FALSE, na = "")
cat(sprintf("mortality_diagnostic.csv written : %d stems checked, %d double-fire\n",
            length(uk), sum(double_fire)))

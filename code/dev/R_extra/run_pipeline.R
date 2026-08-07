#07.06.26 19:45 NZST
# Orchestrator that wires the pipeline stages with explicit paths, no setwd, no
# global state:
#   stems_long -> validate -> impute -> per-stem carbon -> per-ha -> plot
# The whole flow is one call: run_pipeline(stems_long).

# Directory of this file's R/ folder.
pipeline_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) dirname(normalizePath(f)) else "R"
}

source_pipeline <- function(dir = "R") {
  for (f in c("coefficients.R", "allometry.R", "classify.R", "validate.R",
              "expand.R", "summarise.R", "impute.R", "ingest.R"))
    sys.source(file.path(dir, f), envir = globalenv())
}

# Per-stem carbon for a long table with canonical columns.
compute_stem_carbon <- function(stems, co = carbon_coefficients()) {
  woody  <- c("canopy tree", "subcanopy tree", "shrub", "unknown")
  ferns  <- c("tree fern", "palm", "cabbage tree")
  live   <- stems$status %in% c("A", "U", " ")

  agb <- rep(NA_real_, nrow(stems))
  iw  <- live & stems$habit %in% woody & !is.na(stems$dbh) & !is.na(stems$height)
  agb[iw] <- live_woody_agb(stems$dbh[iw], stems$height[iw],
                            stems$wood_density[iw], stems$phylum[iw], co)
  iff <- live & stems$habit %in% ferns & !is.na(stems$dbh) & !is.na(stems$height)
  agb[iff] <- fern_aboveground_carbon(stems$dbh[iff], stems$height[iff], co)

  stems$carbon_agb <- agb
  stems$carbon_bgb <- belowground_carbon(agb, stems$phylum, stems$habit, co)
  stems
}

# Full pipeline from a tidy long table.
run_pipeline <- function(stems_long, co = carbon_coefficients(),
                         vcfg = default_validation_config(),
                         icfg = default_impute_config()) {
  v        <- validate_stems(stems_long, vcfg)
  imputed  <- impute_dbh_gaps(v$clean, icfg)
  withc    <- compute_stem_carbon(imputed, co)
  list(flags = v$flags, quarantine = v$quarantine, stems = withc)
}

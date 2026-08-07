#07.06.26 19:45 NZST
# Reads raw stem and CWD tables into the canonical long schema, one row per
# plot/tag/cycle visited:
#   plot, tag, cycle, status, species_code, phylum, habit, subplot,
#   dbh, height, decay_class, wood_density, mat, led, sed, piece_length, year
# Base R only, no dependencies, for an auditable annual run. `column_map`
# renames raw source columns to the canonical names above.

# Map raw status text to canonical codes.
normalise_status <- function(x) {
  x <- trimws(as.character(x))
  out <- rep("", length(x))
  out[x %in% c("Alive")]                          <- "A"
  out[x %in% c("Unknown", "Not Found", "Dead")]   <- "X"
  out
}

# Map free-text phylum to C (conifer) / B (fern/bryo) / P (broadleaf default).
normalise_phylum <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(x %in% c("pinophyta", "conifer", "conifers"), "C",
  ifelse(grepl("bryo|pterido", x), "B", "P"))
}

# Read a CSV and rename columns per `column_map` (named vector: canonical=raw).
# A mapped source column that is absent stops the read on the native error.
read_canonical <- function(path, column_map) {
  raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  out <- data.frame(row.names = seq_len(nrow(raw)))
  for (canon in names(column_map)) out[[canon]] <- raw[, column_map[[canon]]]
  out
}

canonical_columns <- c("plot", "tag", "cycle", "status", "species_code",
                       "phylum", "habit", "subplot", "dbh", "height",
                       "decay_class", "wood_density", "mat")

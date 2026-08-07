# Known-answer tests for the carbon pipeline. Each check compares a function
# output against a value worked out by hand, so any drift in the maths shows up.
# Base R only, no test framework needed.
# Run from the repository root:  Rscript dev/tests/test_pipeline.R
for (f in c("coefficients.R", "allometry.R", "classify.R", "validate.R", "expand.R"))
  source(file.path("R", f))
source("dev/R_extra/impute.R")   # impute_dbh_gaps lives in dev now

co <- carbon_coefficients()
pass <- 0L; fail <- 0L
check <- function(label, cond) {
  if (isTRUE(cond)) { pass <<- pass + 1L; cat(sprintf("  ok  %s\n", label)) }
  else              { fail <<- fail + 1L; cat(sprintf("FAIL  %s\n", label)) }
}
near <- function(a, b, tol = 1e-8) is.finite(a) & is.finite(b) & abs(a - b) < tol

cat("== allometry ==\n")
# Hand-computed live stem volume: 0.0000483 * (30^2 * 20)^0.978
v_expect <- 0.0000483 * (30^2 * 20)^0.978
check("live_stem_volume", near(live_stem_volume(30, 20, co), v_expect))
check("branch_carbon",    near(branch_carbon(30, co), 0.0175 * 30^2.20))
check("foliage_carbon",   near(foliage_carbon(30, co), 0.0171 * 30^1.75))

# Live woody AGB = density*vol*cf + branch + foliage  (broadleaf cf=0.48)
agb_expect <- 500 * v_expect * 0.48 + 0.0175 * 30^2.20 + 0.0171 * 30^1.75
check("live_woody_agb broadleaf", near(live_woody_agb(30, 20, 500, "P", co), agb_expect))
# Conifer uses cf 0.51
agb_con <- 500 * v_expect * 0.51 + 0.0175 * 30^2.20 + 0.0171 * 30^1.75
check("live_woody_agb conifer", near(live_woody_agb(30, 20, 500, "C", co), agb_con))

cat("== below-ground (the renamed coefficient) ==\n")
check("BGB broadleaf tree -> 0.234", near(belowground_carbon(100, "P", "canopy tree", co), 23.4))
check("BGB conifer tree   -> 0.245", near(belowground_carbon(100, "C", "canopy tree", co), 24.5))
check("BGB tree fern      -> 0.194", near(belowground_carbon(100, "P", "tree fern", co), 19.4))
check("BGB cabbage tree   -> 0.437", near(belowground_carbon(100, "P", "cabbage tree", co), 43.7))

cat("== decay modifier ==\n")
check("decay class 0 -> 1",    near(decay_modifier(0, co), 1))
check("decay class 3 -> 0.47", near(decay_modifier(3, co), 0.47))
check("decay class 4 -> 0",    near(decay_modifier(4, co), 0))

cat("== classification ==\n")
check("growth",    near(carbon_growth(10, 15), 5))
check("mortality", near(carbon_mortality(10, NA), -10))
check("change loss",near(carbon_change(10, 0), -10))
# Ingrowth: visited start, below threshold start, above at end
check("ingrowth recruit (no start carbon)",
      near(carbon_ingrowth(NA, 8, NA, 3.0, TRUE), 8))
check("ingrowth crossing 2.5",
      near(carbon_ingrowth(2, 8, 2.0, 3.0, TRUE), 6))
check("no ingrowth when start visit missing",
      is.na(carbon_ingrowth(NA, 8, NA, 3.0, FALSE)))
check("no ingrowth when already large at start",
      is.na(carbon_ingrowth(2, 8, 3.0, 4.0, TRUE)))

cat("== validation: resurrection vs missed sample ==\n")
stems <- rbind(
  data.frame(plot="p1", tag="dead2alive", cycle=1, status="X", dbh=12),
  data.frame(plot="p1", tag="dead2alive", cycle=2, status="A", dbh=13),
  data.frame(plot="p1", tag="missed",     cycle=1, status="A", dbh=20),
  data.frame(plot="p1", tag="missed",     cycle=3, status="A", dbh=24),
  data.frame(plot="p1", tag="ok",         cycle=1, status="A", dbh=5),
  data.frame(plot="p1", tag="ok",         cycle=2, status="A", dbh=7))
stems$height <- 10; stems$decay_class <- NA; stems$mat <- 11
stems$species_code <- "z"; stems$phylum <- "P"; stems$habit <- "canopy tree"
stems$subplot <- "IN"; stems$wood_density <- 500
v <- validate_stems(stems)
check("resurrection flagged as warning",
      any(v$flags$rule == "resurrection_dead_to_alive" & v$flags$severity == "warning"))
check("missed sample flagged as info (not error)",
      any(v$flags$rule == "missed_sample_gap" & v$flags$severity == "info"))
check("resurrection NOT quarantined (kept for review)",
      !any(v$quarantine$tag == "dead2alive"))
check("missed sample NOT quarantined",
      !any(v$quarantine$tag == "missed"))

cat("== gap interpolation (missed middle cycle) ==\n")
g <- data.frame(plot="p1", tag="m", cycle=c(1,2,3), year=c(2002,2010,2019),
                status="A", dbh=c(20, NA, 30),
                dbh_source=c("measured", NA, "measured"), stringsAsFactors = FALSE)
gi <- impute_dbh_gaps(g)
# linear in time: 2010 is (2010-2002)/(2019-2002)=8/17 of the way
expect_mid <- 20 + (8/17) * (30 - 20)
check("interpolated middle DBH", near(gi$dbh[2], expect_mid, 1e-6))
check("interpolated value tagged", gi$dbh_source[2] == "imputed_interp")
check("measured values untouched", gi$dbh_source[1] == "measured")

cat("== per-hectare expansion ==\n")
# live 30 cm broadleaf in inner subplot -> inner area; per-ha = kg/area/1000
area <- sampling_area_ha("A", 30, "IN", co$plot_area_inner_default_ha, co)
check("inner area chosen for 30 cm live", near(area, 0.033975))
check("circular area for 70 cm live",
      near(sampling_area_ha("A", 70, "IN", 0.033975, co), co$plot_area_circular_ha))
check("per-ha conversion", near(carbon_per_ha(1000, 0.1257), 1000/0.1257/1000))

cat("== taper (Newton-Raphson to 10 cm) ==\n")
h10 <- taper_height_at_10cm(0.5, 20, co)
check("taper h_10 in (0,1)", is.finite(h10) && h10 > 0 && h10 < 1)
check("taper fraction monotonic", taper_volume_fraction_above(0.8, co) >
                                  taper_volume_fraction_above(0.2, co))

cat(sprintf("\n%d passed, %d failed\n", pass, fail))
if (fail > 0) quit(status = 1)

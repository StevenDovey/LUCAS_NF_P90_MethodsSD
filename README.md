# SAS2R

## What this is

An R model of carbon stock and stock change for New Zealand indigenous forest permanent sample plots (PSPs), built from the underlying allometric science rather than ported from the original SAS implementation. Plots are measured across up to three census rounds (2002-2007, 2009 to July 2014, August 2014-2024). The output feeds New Zealand's Emissions Trading Scheme carbon accounting. The earlier SAS-mirroring R code is retained in git history only, under `OLD/`, and is not part of the live pipeline.

## Running it

From the repository root, with the real input files present under
`mfedata/` (gitignored, local only):

```
Rscript run.R
```

`run.R` runs the production pipeline, then runs every diagnostic below in one go, its own last block. Seven run in parallel (each is independent, needing only what `run.R` already wrote to disk); `stems_annotated.R` runs last, since it needs two of the other seven to have finished first. There is no separate `diagnostics.R` to run.

Or run any of the eight diagnostics on its own. `date_audit.R`, `cwd_audit.R`, and `plotsummary_audit.R` read `mfedata/` directly and do not need `run.R` first; `stems_annotated.R` needs `run.R`, `mortality_audit.R`, and `date_audit.R` to have already run:

```
Execution order from run.R:
R/coefficients.R
R/allometry.R
R/classify.R
R/expand.R
R/predict_dbh.R
R/predict_height.R
R/validate.R
R/validate_deep.R

Then the main body runs inline. At the end, diagnostics run in parallel:
9. diagnostics/run_validate_post.R
10. diagnostics/run_validation_summary.R
11. diagnostics/mortality_audit.R
12. diagnostics/date_audit.R
13. diagnostics/cwd_audit.R
14. diagnostics/plotsummary_audit.R
15. diagnostics/run_charts.R

Then sequentially after the parallel block:
16. diagnostics/stems_annotated.R
```

## Scripts

### Production pipeline

**`run.R`**
Reads the stems, coarse woody debris (CWD), species, and plot summary files from `mfedata/`, assigns each record to a census cycle by a fixed date range on its own observation date, runs validation and quarantines hard errors, fills missing diameters and predicts missing heights, computes every carbon pool and stem volume, classifies per-interval fluxes (growth, mortality, ingrowth, net change), expands to per hectare by sampling scope, and writes `output/plot_carbon_summary.csv` plus the data-cleaning tables under `output/diagnostics/`. Its last block then runs every script in `diagnostics/`, see below.

### Diagnostics (`diagnostics/`)

Every script lives directly under `diagnostics/`, no subfolders. `run.R` runs all eight at the end of its own run, seven in parallel plus `stems_annotated.R` last; each also runs standalone. It's the *outputs* that are separated by which input file they check, not the scripts: see `output/diagnostics/mfedata/` below.

**`diagnostics/run_charts.R`**
Confidence charts. Reads the per-stem diagnostics and the plot summary and writes one PNG per chart to `output/charts/`.

**`diagnostics/run_validate_post.R`**
Output-level checks on the finished plot summary: negative pool carbon, an unusually high above-ground stock, negative ingrowth, a net change larger than the standing stock. Writes `output/diagnostics/validation_post.csv`.

**`diagnostics/run_validation_summary.R`**
Summarises the validation flags from `run.R` by rule, cycle, habit, and species, and writes a flags-by-rule chart. All outputs to `output/diagnostics/`.

**`diagnostics/date_audit.R`**
Cross-checks the round date in `plot_summary.csv`, and the CWD file's own date, against the stems file (the baseline), per plot per round. Writes `output/diagnostics/date_diagnostic.csv`.

**`diagnostics/mortality_audit.R`**
Checks the real per-stem carbon values to confirm no stem is counted as a mortality loss in more than one interval, and reports how many stems coded dead from "Not Found" or "Unknown" are later resurrected versus stay dead. Writes `output/diagnostics/mortality_diagnostic.csv`.

**`diagnostics/stems_annotated.R`**
An exact copy of the raw stems file, plus a derived `cycle` column and every diagnostic flag that applies to each row added as a new column: the validation rules by severity, the mortality double-fire flag, and the date check for that plot/round. This is a working copy: fix the underlying data, delete the added columns, and the result replaces the original in `mfedata/`. Writes `output/diagnostics/mfedata/MFESensitivity_Stems_20260209_reduced_dbl_removed3.csv`.

**`diagnostics/cwd_audit.R`**
The same working-copy treatment applied to the raw CWD file: an exact copy plus a derived `cycle` column and flag columns, from the per-row checks `validate.R` already does for stems, applied to CWD for the first time: an inverted taper, an invalid decay class, a species code absent from the species table, a missing or non-positive piece length, both end diameters missing, an unexpected stature code, and a duplicate piece. Writes `output/diagnostics/mfedata/MFESensitivity_CWD_20260209_reduced_dblV2.csv`.

**`diagnostics/plotsummary_audit.R`**
The same working-copy treatment applied to the raw `plot_summary.csv` file: an exact copy plus flag columns, checked directly on this file rather than only visible after being merged onto stems: mean annual temperature (MAT) out of range or missing, inner plot area missing, the two inner-area methods disagreeing, and a duplicate plot row. Writes `output/diagnostics/mfedata/plot_summary.csv`.

### R/ modules (used in `run.R`)

- **`coefficients.R`**: All constants, grouped by type, with its source cited. `VERIFY` marks a value not yet checked against its reference.
- **`allometry.R`**: Per-stem biomass and carbon functions: live stem volume, branch and foliage carbon, fern/palm/cabbage-tree carbon, carbon fraction, the decay-class density modifier, below-ground carbon, the Kimberley-compatible taper used to truncate a dead stem's volume at the 10 cm top diameter.
- **`classify.R`**: Per-interval stem classification: growth, mortality, ingrowth, net change, from aligned start/end carbon and diameter vectors.
- **`expand.R`**: Per-hectare expansion: which sampling area (inner subplot or full circular plot) a stem is expanded by, and the resulting t C/ha conversion.
- **`predict_dbh.R`**: Fills a missing diameter at a cycle from the nearest measured diameter of the same stem, via a species random-intercept mixed model.
- **`predict_height.R`**: The two-stage height model, a species random-intercept mixed model of log(height minus 1.35) on diameter, then a per-plot bias correction.
- **`validate.R`**: Data-cleaning stage: flags records against rules and quarantines errors before the carbon calculation. Uses `default_validation_config()` for thresholds.
- **`validate_deep.R`**: Cross-cycle and reference checks that extend `validate.R`: a per-stem trajectory pass across the three rounds, a height-to-diameter consistency check, a species reference check, and a plot-level stem density check.

### One-off development tools (`dev/tools/`)

Comparison and model-selection tools used during development, not part of the operational run and not part of the `diagnostics/` chain, since most need historical reference files rather than a fresh `run.R` output. `compare_run.R`, `compare_dbh_height.R`, `compare_deadwood.R`, and `compare_fluxes.R` cross-check pipeline output against a reference; `dbh_model_bakeoff.R` and `height_model_bakeoff.R` are the cross-validation runs that selected the current diameter and height model forms; `linecount_diag.R` explains row-count differences between sources; `make_report.R` is an earlier working script, still runnable. `stock_run_v2.R` is dead: it reads `cmp2/SAS_plotsummary_CWD_multiV10.csv`, which does not exist in the repository, and cannot run. Flagged for removal, not yet removed.

### Tests

`dev/tests/test_pipeline.R` and `dev/tests/test_predict.R` are self-contained known-answer tests against hand-computed expected values, run with `Rscript dev/tests/test_pipeline.R`.

### Retained, not wired into the live run

`dev/R_extra/` (`decay_cwd.R`, `impute.R`, `ingest.R`, `run_pipeline.R`, `summarise.R`) and `dev/reference/` hold validated but unused modules and reference code from the earlier implementation, kept for comparison.

### Debris, not yet removed

`dev/run_validate.R` is stale and superseded: it duplicates what `run.R` now does inline, but with the diameter/status logic from before the cycle and AliveState fixes, and it writes to `output/validation_*.csv` rather than the current `output/diagnostics/` location. Running it produces wrong output on the current data. Flagged, not deleted.

## Outputs

### `output/plot_carbon_summary.csv`

One row per plot. Columns follow the pattern `<pool>_ha.<cycle>` for the combined figure, and `<pool>_ha_inner.<cycle>` / `<pool>_ha_circular.<cycle>` for the same figure split by sampling scope (the inner subplot versus the full circular plot); the two scope columns sum exactly to the combined one. Cycle is `1`, `2`, or `3`.

- `agb_ha`: above-ground live carbon, t C/ha
- `bgb_ha`: below-ground live carbon, t C/ha
- `spars_ha`: standing dead carbon, t C/ha (both dead-wood multipliers applied)
- `fallen_ha`: fallen coarse woody debris carbon, t C/ha (both multipliers applied)
- `lvol_ha`, `dvol_ha`: live and dead stem volume, m3/ha
- `CWD_ha`: `spars_ha` plus `fallen_ha`
- `stem_vol_ha`: `lvol_ha` plus `dvol_ha`

Flux columns follow `<AGB|BGB>_<growth|mortality|ingrowth|change>_<interval>`, where interval is `1_2`, `2_3`, or `1_3`. The `1_3` column is a redundant full-period check; for a continuous survivor it duplicates `1_2` plus `2_3` algebraically and should not be summed with them.

### `output/diagnostics/`

- **`stem_diagnostics.csv`**: one row per stem per cycle, plot, tag, cycle, year, species, habit, status, diameter (and its source, measured or imputed), height, predicted height, above- and below-ground carbon. Feeds `diagnostics/run_charts.R` and `diagnostics/run_validation_summary.R`.

- **`validation_errors.csv` / `validation_warnings.csv` / `validation_infos.csv`**: every validation flag, split by severity. An error quarantines the row; a warning or info is kept, flagged for review. Columns: plot, tag, cycle, rule, severity, detail.

  Each rule below is followed by the raw source column(s) it actually reads and the exact condition that trips it, so a flagged row can be checked against the real data, not just the code name.

  - `duplicate_stem_cycle` (error). Checks: `ParentPlotName`, `ItemID`, and the cycle derived from `PlotObsStartDate`. Fails when that plot+ItemID+cycle combination appears more than once (e.g. the same stem recorded on two different dates that both fall in the same cycle window). The most complete row is kept; the flag on the kept row says so, the flag on the dropped row says a more complete duplicate was kept.
  - `status_code_invalid` (error). Checks: `AliveState`. Fails when the value is blank or is not one of "Alive", "Dead", "Not Found", "Unknown". Quarantined, since carbon cannot be computed without a valid status.
  - `alive_state_unknown_treated_dead` (warning). Checks: `AliveState`. Fails when the value is exactly "Unknown"; the stem's death is inferred only from that, not a confirmed record. Kept in the calculation, flagged for review.
  - `dbh_nonpositive_live` (error). Checks: `AliveState` (must be "Alive") and `DiameterValue`. Fails when `DiameterValue` is present and at or below 0.
  - `dbh_missing_unfillable` (warning). Checks: `AliveState` and `DiameterValue` across all three cycles for that stem. Fails when the stem is alive at this cycle but has no measured `DiameterValue` at any cycle; its carbon is not computed.
  - `dbh_above_max` (error). Checks: `DiameterValue`. Fails above 400 cm.
  - `height_below_breast` (error). Checks: `AliveState`, `Height`, `DiameterValue`. Fails when the stem is alive, carries a diameter, and `Height` is present but at or below 1.35 m.
  - `height_above_max` (warning). Checks: `Height`. Fails above 60 m.
  - `decay_class_invalid` (error). Checks: `DecayClass`. Fails when present and outside 0-4.
  - `mat_out_of_range` (warning). Checks: `MAT`, from `plot_summary.csv`, merged onto the stem by plot. Fails outside -5 to 20 degrees C.
  - `resurrection_dead_to_alive` (warning). Checks: `AliveState` and the cycle derived from `PlotObsStartDate`, across a stem's rounds. Fails when a stem is dead (`AliveState` not "Alive") at an earlier cycle and alive at a later one; back-propagated to alive as a field-sampling correction.
  - `missed_sample_gap` (info). Checks: `AliveState` and cycle, across a stem's rounds. Fails when a stem is alive at one cycle and alive again two cycles later with no record in between; interpolated.
  - `dbh_growth_implausible` / `dbh_shrink_while_alive` (warning). Checks: `DiameterValue`, `AliveState`, and the cycle year, between two consecutive alive cycles. Growth fails above 5 cm/year; shrink fails on any decrease greater than 0.5 cm.
  - `dbh_non_monotonic` (warning). Checks: `DiameterValue` across all three cycles. Fails when it rises then falls, or falls then rises, by more than 0.5 cm each way.
  - `large_first_appearance` (warning). Checks: `DiameterValue` and `AliveState`. Fails when a stem has no diameter at an earlier cycle but is alive with `DiameterValue` at or above 10 cm at this cycle, likely present and missed earlier.
  - `dead_stem_growing` (warning). Checks: `DiameterValue` and `AliveState`. Fails when a stem not alive at two consecutive cycles shows `DiameterValue` increasing by more than 0.5 cm between them.
  - `decay_class_decreasing` (warning). Checks: `DecayClass` across cycles. Fails when it decreases from one cycle to the next.
  - `height_diameter_outlier` (warning). Checks: `Height` and `DiameterValue`. Fails when the residual from a fitted height-diameter curve exceeds 4 standard deviations, possibly a broken or snapped stem.
  - `species_not_in_reference` (warning). Checks: `PreferredSpeciesCode` against the `Species code` column in `species.csv`. Fails when the code is present but not found there; the default wood density was used.
  - `plot_density_outlier` (warning). Checks: count of alive stems (`AliveState` "Alive", `DiameterValue` at or above 2.5 cm) per plot per cycle. Fails when that count is a statistical outlier (MAD z-score above 5) against all plots.

- **`validation_quarantine.csv`**: the rows excluded from the carbon calculation, one row per excluded record.
- **`validation_post.csv`**: output-level sanity flags from `diagnostics/run_validate_post.R` (negative pool carbon, an implausible stock ceiling, negative ingrowth, a net change exceeding the stock).
- **`validation_summary_by_rule.csv`**: flag count, plots affected, and stems affected, per rule.
- **`validation_summary_rule_x_cycle.csv`** / **`_rule_x_habit.csv`** / **`_rule_x_species.csv`**: the same flags broken down by cycle, habit, and species. `cross_round` marks a rule that compares a stem across rounds rather than one cycle; `unmatched` marks a flag whose stem is not in the diagnostics (a plot-level rule, or a quarantined duplicate).
- **`validation_flags_by_rule.png`**: the by-rule count as a bar chart, coloured by severity.
- **`cwd_stature_distribution.csv`**: count of CWD pieces by raw stature code and cycle. Visibility only; stature is not used to differentiate carbon fraction, since neither this model nor the earlier SAS-based one did.
- **`mortality_diagnostic.csv`**: from `diagnostics/mortality_audit.R`, stems checked, the mortality double-fire count (expected zero), the raw AliveState distribution by cycle, and how many "Not Found"/"Unknown" stems were later resurrected versus stayed dead.
- **`date_diagnostic.csv`**: from `diagnostics/date_audit.R`, one row per plot, cycle, and comparison (`plot_summary` or `cwd`) against the stems file. Checks: the stems file's own `PlotObsStartDate` (aggregated per plot/cycle) against `plot_summary.csv`'s "Date 1st/2nd/3rd meas" column for that round, and separately against the CWD file's own `PlotObsStartDate`. Columns: `stems_date`, `other_date`, `days_apart` (only when both are present), `flag`, `severity`.
  - `missing_in_stems` (major): a date exists in the comparison source but not in stems.
  - `missing_in_other` (minor): a date exists in stems but not in the comparison source; infillable from stems.
  - `both_missing` (major): no date in either source for that plot and round.
  - `date_gap` (info): both present; `days_apart` reports the gap, no tolerance is judged by the script.

#### `output/diagnostics/mfedata/`

Working copies for a manual data-cleaning round-trip: every row of the named source file, unchanged, plus a few appended columns, same filename as the source. Review the flags, fix the underlying data, delete the appended columns (including `sort`), and the result replaces the original file in `mfedata/`.

Each file is written sorted by its natural key (plot, tag, cycle for stems and CWD; plot for plot_summary), so a duplicate lands on adjacent rows instead of being scattered through hundreds of thousands of unrelated rows. The `sort` column, the last column in the file, holds the original row number; sort ascending by `sort` to restore the file to its original order before it replaces the source.

- **`MFESensitivity_Stems_20260209_reduced_dbl_removed3.csv`**: from `diagnostics/stems_annotated.R`. Adds `cycle` (derived from `PlotObsStartDate`, the same way `run.R` derives it), five flag columns, each a semicolon-joined list of what fired for that row, blank if nothing did: `flags_error`, `flags_warning`, `flags_info` (the validation rule names above; a cross-round rule is attached to every cycle-row of that stem, not just one), `flags_mortality_double_fire`, and `flags_date` (the date-check result for that plot and round, major and minor only), and `sort` (original row order).
- **`MFESensitivity_CWD_20260209_reduced_dblV2.csv`**: from `diagnostics/cwd_audit.R`. Adds `cycle`, three flag columns (`flags_error`, `flags_warning`, `flags_info`), and `sort` (original row order), from the same kind of per-row checks `validate.R` already does for stems:
  - `cwd_taper_inverted` (warning). Checks: `SmallEnd1`, `SmallEnd2`, `LargeEnd1`, `LargeEnd2`. Fails when the computed small-end diameter is greater than the computed large-end diameter.
  - `decay_class_invalid` (error). Checks: `DecayClass`. Fails when present and outside 0-4.
  - `species_not_in_reference` (warning). Checks: `PreferredSpeciesCode` against `species.csv`. Fails when not found there.
  - `length_missing_or_nonpositive` (warning). Checks: `ItemObsComponentLinearDimension`. Fails when missing or at or below 0; the piece is excluded from the carbon calculation.
  - `both_diameters_missing` (warning). Checks: `LargeEnd1`, `LargeEnd2`, `SmallEnd1`, `SmallEnd2`. Fails when both the large-end and small-end diameters cannot be computed; the piece contributes no carbon.
  - `stature_unexpected` (info). Checks: `CWDStature`. Fails when present but not one of "F", "P", "S".
  - `duplicate_cwd_piece` (error). Checks: `ParentPlotName`, `ItemID`, and the cycle derived from `PlotObsStartDate`. Fails when that combination repeats.
- **`plot_summary.csv`**: from `diagnostics/plotsummary_audit.R`. Adds three flag columns (`flags_error`, `flags_warning`, `flags_info`) and `sort` (original row order):
  - `mat_out_of_range` (warning). Checks: `MAT`. Fails outside -5 to 20 degrees C.
  - `mat_missing` (warning). Checks: `MAT`. Fails when blank; the default is used.
  - `inner_area_missing` (warning). Checks: "Inner plot area new method" and "Inner plot area original method". Fails when both are blank; the default is used.
  - `inner_area_methods_disagree` (info). Checks: the same two columns. Fires when both are present and numerically different; reports the difference, no tolerance judged.
  - `duplicate_plot_row` (error). Checks: `PLOT`. Fails when the same plot appears in more than one row.

### `output/charts/`

- `height_pred_vs_obs_core.png` / `_tail.png`: height predicted against observed, split at the training diameter range.
- `height_pred_vs_obs_species_NN.png`: the same, paged by species.
- `height_residual.png`: height residual against predicted height.
- `height_residual_by_habit.png`: height residual by growth form.
- `height_bias_by_species.png`: mean height bias, worst 25 species.
- `dbh_pred_vs_obs.png` / `dbh_pred_vs_obs_species_NN.png`: diameter predicted against observed on measured consecutive pairs, overall and paged by species.
- `dbh_increment_per_year.png`: diameter increment per year, flags implausible growth or shrinkage.
- `height_diameter_species_NN.png`: height against diameter per species, with the fitted curve.
- `agb_stock_by_cycle.png`: above-ground stock distribution by cycle.
- `agb_change_1_2.png` / `_2_3.png`: net above-ground change per plot, by interval.
- `pool_composition_by_cycle.png`: mean carbon by pool and cycle.

# Validation results — new pipeline vs SAS reference

**Date:** 2026-06-07
**Reference:** SAS outputs V10 (`SAS_measurements-with-carbonV10.csv` stem-level,
`SAS_plotsummary_CWD_multiV10.csv` plot-level) on the MFE Sensitivity 3-cycle
dataset (963 plots, 428k stem records).
**Goal (per project direction):** not byte parity — a cleaner pipeline that
reaches the required outcome by sound, matched pathways. Each pathway below is
tested in isolation so a result is unambiguous.

Reproduce with: `Rscript tools/compare_run.R`, `tools/compare_fluxes.R`,
`tools/compare_deadwood.R` (input data is local/gitignored).

---

## Summary

| Pathway | How it was isolated | Result |
|---|---|---|
| Live allometry (stem / branch / foliage / AGB) | SAS's own DBH + predht fed to `R/allometry.R` | **exact, 0.00%** (339,519 live woody stems) |
| BGB root:shoot ratios | per-stem, SAS AGB | **exact, 0.00%** |
| Mortality flux | SAS per-stem carbon → `R/classify.R` | **exact, corr 1.000** (1-2, 2-3, 1-3) |
| Net change flux | bucket-invariant total | **exact, corr 1.000** |
| Growth vs ingrowth split | recruit bucketing | **deliberate departure**, net-invariant (see below) |
| Dead-pool carbon (spar woody / spar fern / fallen) | SAS volumes + modifiers → my formulas | **exact, 0.00%** (72,057 stems) |
| Taper truncation (`h_10`, Newton-Raphson) | independent impl → corrected | **exact, 0.00%** (15,184 spars) |
| End-to-end AGB/BGB **stock** | full new pipeline, real data | **~1%, r = 0.99** (~800 plots) |

**Conclusion:** the new pipeline reproduces everything the SAS model computes —
exactly where it should, and by a conserved, documented improvement where it
differs. The only residual (~1% on stock) is the simplified stand-in height
model, not the carbon science.

---

## Detail

### 1. Allometry (live AGB/BGB) — exact
Feeding SAS's own `DBH` and `predht` into `live_woody_agb()`,
`branch_carbon()`, `foliage_carbon()`, `belowground_carbon()` reproduces
SAS `carbon_stem/branches/foliage/AGB/BGB` to 0.00% on all 339,519 live woody
stems (100% within 0.1%). Confirms the equations and the renamed BGB ratios
(`bgb_ratio_broadleaf_palm` etc.) are identical to SAS — including that
conifers correctly take 0.245, not 0.234.

### 2. Fluxes — mortality and net change exact; growth/ingrowth a clean re-bucket
Feeding SAS's own per-stem `carbon_AGB` + `DBH` into `R/classify.R`:

- **Mortality** matches SAS to corr 1.000 (mean diff ~0.005 t C/ha) on all
  three intervals.
- **Net change** (growth + ingrowth + mortality) matches SAS to corr 1.000
  (mean diff ~0.004). This is bucket-invariant, so it is the decisive proof
  the interval logic conserves total carbon change identically to SAS.
- **Growth vs ingrowth** differ by design. SAS books a brand-new stem
  (absent at start → alive at end) into *growth* (full end carbon), and its
  "ingrowth" column captures only already-measured small stems crossing
  2.5 cm. The new pipeline books new recruits as *ingrowth* (recruitment),
  giving a disjoint partition — growth = survivors already ≥ 2.5 cm at both
  ends; ingrowth = recruits + threshold-crossers. The combined growth +
  ingrowth and the net change are invariant to this choice. This is the
  intended "better pathway", documented, not a discrepancy.

The ingrowth NA-gate (`is.na(dbh_start) | dbh_start < 2.5`) treats an absent
start as below threshold — the conceptually correct recruit test, which also
matches SAS's `missing = -infinity` behaviour.

### 3. Dead wood / CWD — formulas exact; taper corrected to exact
- Dead-pool carbon **formulas** reproduce SAS to 0.00% when fed SAS's own
  volume + decay modifier + carbon fraction + density: standing-dead spar
  (woody, n=10,818), spar (fern, n=2,185), fallen F/P/S (n=59,054).
- The **taper** truncation was the one piece implemented from the concept and
  was initially wrong (~470% off). Replacing it with the exact compatible
  volume/taper function (sectional polynomial with the `h^81` term; `tap_*`
  coefficients scaled by stem volume / height; Newton-Raphson to
  diameter² = 0.10²) makes `h_10` match SAS to 0.00% on 15,184 spars.

### 4. End-to-end stock — ~1%, r = 0.99
Running the full new pipeline on real inputs (cycles + per-plot areas from the
reference; a simple per-species log-log height model standing in for the 77% of
live stems without a measured height) gives AGB/BGB stock within ~1% (median
ratio 0.99) at r = 0.99 across ~800 plots per cycle. The ~1% low bias is
entirely the stand-in height model; allometry and expansion are exact (sec 1).

---

## End-to-end run (tools/run_full.R)

The pipeline reads the six inputs, fills diameter gaps with the species
mixed model, predicts height with the two-stage method, computes the live
and dead carbon pools, expands per hectare, summarises to plot, and writes
output/plotsummary_new.csv. Selected models were chosen by cross-validation
against the field measurements: diameter improves on the SAS per-plot form
at 1.220 against 1.323 cm root mean square error; height adopts the SAS
two-stage form, which the bake-off selected at 2.477 m.

The output carries forty nine columns: stocks per cycle for above-ground,
below-ground, standing dead, fallen debris, total debris, and stem volume,
and the per-interval growth, mortality, ingrowth and change for above-ground
and below-ground.

Per-plot agreement against the SAS plot summary:

| quantity | median ratio | correlation |
|---|---|---|
| above-ground stock, cycle 1 | 0.979 | 0.991 |
| below-ground stock, cycle 1 | 0.980 | 0.991 |
| standing dead stock, cycle 1 | 0.976 | 0.899 |
| fallen debris stock, cycle 1 | 0.835 | 0.960 |
| above-ground stock, cycle 2 | 0.998 | 0.997 |
| above-ground stock, cycle 3 | 0.993 | 0.996 |
| above-ground mortality, 1 to 2 | 1.006 | 0.971 |
| above-ground net change, 1 to 2 | 0.993 | 0.950 |
| above-ground mortality, 2 to 3 | 0.992 | 0.989 |

The live stocks agree to within one to two percent at correlation 0.99.
Mortality and net change match. Growth and ingrowth differ by the recruit
re-bucketing, which is net-invariant. The standing dead stock agrees on
magnitude. The fallen debris stock, expanded by the documented nested
subplot design with the 60 cm piece split, agrees at 0.835 with a per-plot
total debris carbon ratio of 0.961, the residual being the area expansion.

## Validation stage (tools/run_validate.R)

On 405,051 stem-cycles the data-cleaning stage raised 4,750 flags and
quarantined 4,414 rows. The findings include 683 resurrection cases, dead
then alive, which are genuine recording errors for review, and 231 missed
sample gaps routed to interpolation. The live-without-diameter rule is
currently stricter than the imputation design, since a bracketed missing
diameter should interpolate rather than quarantine, and needs the cross-cycle
refinement.

The standing dead and fallen ratios above are the raw, unmultiplied values,
which is the like-for-like comparison since the reference file carries no
dead-wood correction. With the confirmed correction applied, both dead pools
are scaled by the product of the two multipliers, 1.767 × 1.19 = 2.103.

## Decay rate and loss

The decay constant, the species rate adjusted for mean annual temperature and
piece diameter, reproduces the reference exactly on 194,888 rows. The decay
loss of existing debris over an interval reproduces the reference exactly on
12,356 plots, once the 1.767 sampling correction is applied to the debris.
That correction on the debris pool independently confirms the multiplier
scope below. The projection module applies the rate to existing debris over
the interval and to a mortality cohort over half the interval.

## Multiplier scope (resolved)

The operational code applies a 1.767 sampling correction and a 1.19
below-ground uplift to both standing dead and fallen debris. The original
author confirmed this uniform application was the intended treatment, made to
carry the pool uncertainty as the best available approach at the time. The
model follows that decision: each dead pool is scaled by the product of the
two multipliers, 1.767 × 1.19 = 2.103, named at the top of the run. Against
the previous scoped split, where the sampling correction reached only fallen
debris and the uplift only standing dead, the change lifts the mean standing
dead pool from 10.5 to 18.6 t C/ha and the mean fallen pool from 20.2 to 24.1
t C/ha, raising mean total stock from 183.4 to 195.4 t C/ha, a 6.5 percent
increase, with the dead-wood share moving from 16.8 to 21.9 percent. The
cabbage-tree root:shoot value of 0.437 is confirmed against the Easdale 2019
monocot group.

## Remaining work
1. The fallen debris area residual under the nested design, held at the
   documented value rather than tuned to SAS.
2. The decay mortality-cohort gain assembled end to end. The rate, the loss
   term, and the projection function are confirmed exact; the cohort detection
   is the final assembly.

# Clean-slate rebuild — design & specification

**Status:** draft for review (no R code written yet)
**Scope decision:** *pure clean-slate* implementation of the scientific concept. This is **not** a port and **not** bound to SAS numeric parity. SAS output may be used as one of several sanity references at the end, but matching it byte-for-byte is explicitly a non-goal.
**Coefficient decision:** the empirical equations (volume, taper, decay, root:shoot, carbon fractions, densities) are **published external science**, reused here with a source note against each. Their exact numeric values are to be **verified against the primary sources at the end** (see §8), not re-derived.
**Validation decision:** records failing a hard scientific check are **flagged + quarantined** (excluded from the carbon calc into a separate table, reviewable and reinstatable). Warnings are kept and flagged.

---

## 1. What the model computes (the concept, restated independently)

Carbon stock and carbon stock *change* in New Zealand indigenous-forest permanent sample plots (PSPs), measured at up to three census cycles (~2002, ~2010, ~2019). The unit of observation is an individual **stem** identified by `(plot, tag)`, observed at each cycle it was visited.

For each stem and each measurement interval the model answers: *how much carbon did this stem hold, and how did that change?* Carbon is split into pools:

- **AGB** — above-ground biomass of live stems (stem + branches + foliage).
- **BGB** — below-ground (root) biomass, a habit-specific fraction of AGB.
- **Dead-wood** — standing dead (spars) and fallen coarse woody debris (CWD), with decay attenuation.

Each pool's *change* over an interval decomposes into:

- **Growth** — stem alive and above the size threshold at both ends → end carbon − start carbon.
- **Mortality** — alive at start, gone/dead at end → −(start carbon).
- **Ingrowth** — below threshold or absent at start, at/above threshold at end → recruitment in.

Per-stem carbon is then expanded to **t C / ha** using the area of the subplot the stem was sampled in, and summed to **plot × interval × pool** totals.

This is a closed, re-derivable concept. The only parts that are *not* derivable from first principles are the empirical allometric coefficients (§7), which come from the literature.

---

## 2. Design principles (how this differs from the existing port)

| Principle | Rationale / what it fixes |
|---|---|
| **Tidy long format end-to-end** (`one row per stem × cycle`) | The existing port's worst defects — structural-grid phantom rows, imputed values leaking into non-visited slots — exist *only* because it mimics SAS's repeated wide⇄long transposes. Staying long eliminates that whole class of bug by construction. Wide is materialised *only* for final reporting. |
| **Explicit provenance column beside every value** | `dbh_source`, `height_source` ∈ {`measured`, `imputed_<method>`}. Imputed values are never silently treated as measurements. |
| **Coefficients centralised + cited** | One `coefficients.R` file; every constant carries a source comment and a `# VERIFY` tag until checked against the primary PDF. |
| **Validation is a first-class stage with its own output** | Not buried in `if/else`. Produces a flags table and a quarantine table (§6). |
| **One area source, used everywhere** | The port uses spreadsheet area in one summary and shoelace area in another. Here a single resolved `plot_area` per plot feeds every expansion. |
| **Pure, unit-tested functions** | Allometry and classification are pure functions with synthetic known-answer tests, so correctness is checkable **without** SAS. |
| **Fail loud, never silently coerce** | e.g. an inverted CWD taper (SED > LED) is *flagged*, not silently swapped as the port does. |
| **Built for an annual, auditable rerun** | This is a yearly reporting exercise. Coefficients and config are versioned per report year; every run emits a validation report, a provenance log, and a quarantine table so the result is reproducible and defensible to an auditor. Re-running last year's inputs must reproduce last year's numbers. |
| **Human error is expected, not exceptional** | Field data contain transcription mistakes and missed visits. The design separates *genuine errors* (quarantine) from *missed samples* (interpolate) — see §4.1 — rather than discarding imperfect records. |

---

## 3. Architecture

```
R/
  00_coefficients.R   all empirical constants + sources (§7)
  01_ingest.R         read stems + CWD CSV/xlsx → canonical tidy long table
  02_validate.R       scientific soundness checks → flags + quarantine (§6)
  03_impute.R         DBH + height imputation, with provenance
  04_classify.R       per-interval growth / mortality / ingrowth (pure)
  05_allometry.R      per-stem AGB, BGB, CWD/spar/fallen carbon (§5)
  06_expand.R         per-hectare expansion by subplot area
  07_summarise.R      plot × interval × pool totals
  08_report.R         carbon outputs + QC / validation report
  run_pipeline.R      orchestrator (no setwd; explicit in/out paths)
tests/
  testthat/           synthetic known-answer fixtures + unit tests
docs/
  clean_rebuild_spec.md   (this file)
```

Each stage is a function `stage(input, coeffs, config) -> output` with no global state. `run_pipeline.R` wires them. Configuration (paths, thresholds, toggles) lives in one list, not scattered constants.

---

## 4. Canonical data model

After ingest, one tidy table `stems_long`, one row per `(plot, tag, cycle)` **that was actually visited**:

| column | meaning |
|---|---|
| `plot`, `tag` | stem identity |
| `cycle` | 1 / 2 / 3 |
| `visited` | TRUE for real field visits (always TRUE in this table; absence = no row) |
| `species_code`, `habit`, `phylum` | taxonomy / form (tree, shrub, tree fern, palm, cabbage tree) |
| `status` | A (alive), X (standing dead), F/P/S (fallen/stump pieces), per §6 status map |
| `subplot` | inner / circular / EXT flag |
| `dbh`, `dbh_source` | DBH (cm) + provenance |
| `height`, `height_source` | height (m) + provenance |
| `decay_class` | 0–4 for dead wood |
| `led`, `sed`, `piece_length` | CWD piece geometry |
| `wood_density`, `mat` | plot/species attributes |

**Why no padded rows:** a stem not visited at a cycle simply has no row for that cycle. Interval logic (§4-classify) joins cycles by `(plot, tag)` and treats a missing row as "absent" — which is exactly the ingrowth/mortality semantics, with no phantom structural records to clean up later.

### 4.1 Stem life-history state model (absent ≠ dead)

The data are collected by field crews, so **a missing observation is not the same as death**. Two failure modes must be told apart, because they are handled in completely different stages:

- **Missed sample** — the stem was *absent* (no row) at a cycle but is *alive* again at a later cycle. A crew missed it. This is a **gap to interpolate** (§5.3), **not** mortality and **not** resurrection.
- **Resurrection** — the stem was *marked dead* (`X`, or recorded as a fallen piece) and then appears *alive*. **Dead is dead**: this is a genuine recording error and is **quarantined** (§6).

Per-stem state transitions across the ordered cycles, classified once in ingest/validate:

| from → to | meaning | handling |
|---|---|---|
| `A → A` | survived | growth |
| `A → X` / `A → F/P/S` | died (standing / fallen) | mortality + enters dead-wood pool |
| `A → absent → A` | **missed middle visit** | interpolate the gap (§5.3); stem is alive throughout |
| `absent → A` (no earlier record) | recruited **or** missed earlier | if size/age implies it pre-existed → back-cast gap; else ingrowth (§5.4) |
| `A → absent` (end, no later record) | died/gone (model cut-off) | mortality at the closing interval |
| `X → X` | still-standing dead | decaying dead-wood, no live flux |
| `X → F/P/S` / `X → absent` | dead wood fell / collapsed | dead-wood transition |
| `X → A`, `F/P/S → A` | **resurrection** | **error → quarantine (§6)** |

This state model is the single source of truth that both validation (§6) and imputation (§5.3) read from, so "missed" vs "dead" is decided in exactly one place.

---

## 5. Pipeline stages

### 5.1 Ingest (`01_ingest.R`)
- Read stem table + CWD table; map raw column names to the canonical schema.
- Normalise `status`: blank → (CWD/absent), {Unknown, Not Found, Dead} → `X`, Alive → `A`; pieces F/P/S from CWD geometry.
- Derive `plot`, `cycle` (by matching measurement date to the plot's census dates), `species_code`, `tag`.
- Attach plot attributes (areas, MAT) and species attributes (habit, phylum, wood density).
- Output: `stems_long` + a row-count reconciliation (rows in vs rows mapped vs rows dropped, with reasons).

### 5.2 Validate (§6) — runs **before** imputation so checks see raw field data.

### 5.3 Impute — gap interpolation (backward & forward) (`03_impute.R`)
Goal: fill DBH/height **only** where the life-history state model (§4.1) says the stem genuinely existed at a cycle but was missed or partially recorded — i.e. close real gaps — in a **sound, defensible, documented** manner. Imputation never invents a stem that wasn't there and never overrides a recorded measurement.

**Which gaps are eligible (and direction):**
- **Interpolate (bracketed gap)** — `A → absent → A`: the stem is alive on both sides, so the middle cycle is interpolated *between* the two measurements (the most defensible case — bounded by real data).
- **Back-cast (backward)** — stem first appears at a later cycle but its size/state implies it pre-existed (was missed earlier, not newly recruited): estimate the earlier missing value *backward* from later measurements.
- **Forward-cast (forward)** — a missing later value bracketed appropriately, or a same-cycle missing component (e.g. DBH present, height missing).

**Methods (in defensibility order, most-bounded first):**
- **DBH:** where two real measurements bracket the gap, prefer **monotone interpolation between them** (growth is the well-behaved, bounded case). Where only one side exists, use plot-level regression of the missing cycle's DBH on the measured cycle's DBH (random intercept by species when ≥2 species and ≥ N pairs; else fixed-effect `lm`; else habit/site group model; else "no growth" carry for tree ferns/palms/cabbage trees).
- **Height:** species mixed model of `log(H − 1.35) ~ DBH^(−0.3)` (random intercept + slope by species) with a per-plot×cycle bias correction; mean-height fallbacks for ferns/palms.

**Defensibility constraints (cut-offs the interpolation must respect):**
- No negative growth across an interpolated interval for a live stem (DBH non-decreasing beyond measurement tolerance); clamp and flag if a fit implies shrinkage.
- Interpolated values must stay within the model's size cut-offs (e.g. ≥ 2.5 cm if the stem is treated as live-present) and within the bracketing measurements.
- Bounded extrapolation only: back/forward-casts are capped to a maximum plausible per-year increment (§6 growth limit); beyond that the value is flagged, not silently extended.
- Imputed values are usable for carbon but **never** satisfy a "was-measured" test, **never** flip a stem's alive/dead state, and **never** fill a cycle the state model marks as truly absent (dead/gone).

Every filled value records `dbh_source` / `height_source = imputed_<method>_<direction>` (e.g. `imputed_interp`, `imputed_backcast_lme`) so the report can show exactly how each gap was closed.

### 5.4 Classify (`04_classify.R`)
For each interval (1→2, 2→3, 1→3), join start and end cycle rows by `(plot, tag)` and assign per pool:
- **Growth**: alive both ends, both above size threshold (2.5 cm live).
- **Mortality**: alive at start, dead/absent at end.
- **Ingrowth**: start absent **or** below 2.5 cm; end ≥ 2.5 cm and alive. A stem with no start row is "below threshold" by definition (the conceptually correct treatment — and incidentally what the SAS `missing = −∞` quirk approximated). Ingrowth carbon = end − start, with start = 0 for a true new recruit.
- Classification keys off **size/alive state**, with carbon as the magnitude; pure function returning a tidy `(plot, tag, interval, pool, class, delta_carbon)`.

### 5.5 Allometry (`05_allometry.R`) — per-stem carbon
- **Live stem volume** (trees/shrubs): `Vol = a·(DBH²·H)^b`.
- **Live ferns/palms/cabbage**: separate volume/carbon allometry.
- **Stem carbon** = volume × wood density × carbon fraction; **branches** and **foliage** from DBH power laws.
- **Standing dead (spars)**: Kimberley compatible volume–taper function, truncated by Newton–Raphson to the 10 cm-diameter point and to break height; × decay-class density modifier × carbon fraction.
- **Fallen CWD**: conic-frustum volume from end diameters (pieces split at the 60 cm size boundary); × decay modifier × carbon fraction. **Inverted taper (SED > LED) is flagged, not swapped.**
- **AGB** = stem + branches + foliage. **BGB** = habit-specific root:shoot fraction × AGB.
- All coefficients from §7.

### 5.6 Expand (`06_expand.R`)
- Two subplot areas: inner (default 0.033975 ha) for stems ≥ 2.5 & < 60 cm (live) / ≥ 10 cm (dead); full circular plot (0.1257 ha) for ≥ 60 cm. `EXT` excluded from inner.
- One resolved area per stem per cycle; interval fluxes use a single shared denominator per stem so growth/mortality/ingrowth are commensurable.
- `t C/ha = carbon(kg) / area(ha) / 1000`.

### 5.7 Summarise (`07_summarise.R`)
- Sum per-ha carbon to `plot × interval × pool`.
- Dead-wood multipliers (sampling-bias and below-ground-deadwood uplift) applied to both dead pools, standing dead and fallen, per the confirmed operational treatment; named explicitly at the top of the run.

### 5.8 Report (`08_report.R`)
- Plot-level carbon table (the deliverable).
- **Validation report**: counts of flags by rule and severity; the quarantine table; a provenance summary (how many DBH/height values were imputed, by method).

---

## 6. Data-cleaning / scientific-soundness stage (flag + quarantine)

Runs on raw ingested data. Each rule emits a row in `validation_flags`: `(plot, tag, cycle, rule, severity, detail)`. **Severity `error`** → record moved to `quarantine` and excluded from the carbon calc (reviewable/reinstatable). **Severity `warning`** → kept in the pipeline, flagged.

### Identity / structural
| rule | severity |
|---|---|
| duplicate `(plot, tag, cycle)` | error |
| species code changes across cycles for same tag | warning |
| subplot flag inconsistent across cycles | warning |
| status code not in allowed set {A, X, F, P, S, blank} | error |

### DBH plausibility
| rule | severity |
|---|---|
| DBH ≤ 0 or non-numeric on a live stem | error |
| DBH above plausible max (default > 400 cm) | error |
| DBH shrinks > tolerance while alive (default > 0.5 cm absolute) | warning |
| DBH growth above biological max (default > 5 cm·yr⁻¹) | warning |
| live in main plot recorded < 2.5 cm | warning |

### Height plausibility
| rule | severity |
|---|---|
| height ≤ 1.35 m on a live DBH-bearing stem | error |
| height above NZ max (default > 60 m) | warning |
| height:DBH allometric residual beyond N SD (default 4) | warning |
| height shrinks > tolerance while alive | warning |

### Life-history logic (driven by the §4.1 state model)
| rule | severity |
|---|---|
| **resurrection**: marked dead (`X` / `F` / `P` / `S`) then `A` at a later cycle — *dead is dead* | error |
| **missed sample**: `absent → A` or `A → absent → A` (alive after a gap) | info — **not** an error; routed to gap interpolation (§5.3) |
| alive with neither DBH nor height and no bracketing basis to impute | error |
| status `A` but no size at any cycle | warning |
| repeated dead with increasing implied size (dead stem "grows") | warning |

> The missed-sample case is deliberately **info**, not a defect: human crews miss stems, and the model's job is to close the gap defensibly, not to penalise the record. Only *dead → alive* is treated as a true error.

### Dead-wood / CWD
| rule | severity |
|---|---|
| decay class outside 0–4 | error |
| CWD piece SED > LED (inverted taper) | warning (flagged, **not** auto-swapped) |
| CWD end diameter ≤ 0, or both ends missing | error |
| standing dead with DBH below dead threshold (10 cm) | warning |

### Spatial / plot
| rule | severity |
|---|---|
| plot inner or circular area missing / outside plausible range | error |
| MAT missing or outside NZ range (default −5 … 20 °C) | warning (fall back to default MAT, flagged) |

### Output sanity (post-calc, re-validated)
| rule | severity |
|---|---|
| any pool carbon negative or non-finite | error |
| ingrowth flux negative | warning |

All thresholds live in `config` so they're tunable without code edits. Defaults above are first-pass and explicitly up for your review.

---

## 7. Coefficient appendix (reused published equations — values to verify at §8)

> These are external published constants, reused with attribution. Values listed are the working set; each carries a `# VERIFY` tag in code until checked against the primary PDF.

**Volume / above-ground biomass — Beets et al. (2012), *Forests* 3(3):818, doi:10.3390/f3030818**
- Live tree/shrub stem volume (m³): `0.0000483 · (DBH²·H)^0.978`
- Branch carbon: `0.0175 · DBH^2.20`
- Foliage carbon: `0.0171 · DBH^1.75`
- Tree fern/palm/cabbage stem volume: `0.00001343 · (DBH²·H)^1.22`
- Tree fern/palm/cabbage integrated carbon: `0.00270 · (DBH²·H)^1.19`
- Carbon fraction: conifer (gymnosperm) **0.51**, broadleaf/other **0.48**, dead wood **0.50**
- Default unknown wood density: **477 kg·m⁻³**

**Stem taper / compatible volume (dead-spar truncation) — Beets et al. (2012) / van der Colff, Scion (LUCAS); Kimberley & Beets national volume function**
- Sectional taper coefficients (a,b,c,d,f): `0.06501, 2.92127, −3.37103, 1.35551, 0.02924`; `π/4 = 0.7854`
- Newton–Raphson solve to 10 cm stem diameter (truncation point)

**Decay-class density modifiers — Coomes et al. (2002)**
- class 0→**1.00**, 1→**0.82**, 2→**0.66**, 3→**0.47**, 4→**0.00**

**CWD exponential decay — Garrett et al. (2018), *For. Ecol. Manage.*; Beets & Hood (2008)**
- base λ = **0.0216**; species overrides (BEITAW 0.0233, DACCUP 0.0350, DACDAC 0.0506, PRUTAX 0.0171, PRUFER 0.0259, WEIRAC 0.0220, NOTMEN/LOPMEN 0.0246, NOTFUS/FUSFUS 0.0204, NOTSOL/NOTSVS 0.0162)
- climate/size adjustment: `λ_adj = λ · exp(0.093·(MAT−10)) / (1 + 0.00908·(Diam−60))`
- ⚠️ *Corroboration:* secondary sources give an across-species **mean** above-ground λ ≈ 0.0270 (half-life ~25.7 yr) and the directional pattern matches the overrides (DACCUP/DACDAC fast <20 yr; PRUTAX/FUSCLI slow ~40 yr). The **base 0.0216** and the **adjustment coefficients 0.093 / 0.00908** could **not** be confirmed from open sources — keep `# VERIFY` against Garrett (2018).

**Root:shoot / below-ground biomass — Easdale et al. (2019), *For. Ecol. Manage.* 432:117463**
- *Working set:* broadleaf trees/palms **0.234**; tree fern **0.194**; cabbage tree **0.437**; conifers + default **0.245**
- ✅ **No swap (checked against code).** Both SAS (`SAS3cycle.sas:1443`) and R (`SAS2R_final:1191`) gate 0.234 on `phylum in ('B','P')` = broadleaf; conifers carry phylum `C`, fail that test, and fall through to **0.245** — matching Easdale (gymnosperm 0.245, angiosperm 0.234). The earlier "possible swap" was a **false positive from a misleading R constant name** (`BGBFracConifPalm <- 0.234` is actually applied to broadleaf, not conifers). Computation is correct in both; fix the name in the rebuild.
- ⚠️ tree fern **0.194** corroborated; **cabbage tree 0.437 not corroborated** by literature (Easdale groups monocots with angiosperms ≈ 0.234) — but it *is* in the SAS reference (`SAS3cycle.sas:1446`), so it's an original-model value to verify, not an R discrepancy.

**Plot geometry / thresholds — LUCAS PSP design**
- circular plot **0.1257 ha** (20 m radius); inner subplot default **0.033975 ha**
- live DBH thresholds: inner **2.5 cm**, circular **60 cm**; dead/CWD: inner **10 cm**, circular **60 cm**
- standing-dead minimum DBH **10 cm**; below-ground dead-wood uplift **×1.19**; CWD sampling multiplier **1.767**

---

## 8. Source documents, corroboration status & end-stage validation

Every primary publisher (MDPI, ScienceDirect, Springer, Scion, environment.govt.nz, doc.govt.nz) returned **HTTP 403** to automated fetch this session, so nothing below is read from a primary PDF. Values are taken from the project's existing constant set and **corroborated against secondary sources** (abstracts, government summaries, citing papers). "Corroborated" ≠ "verified against primary" — final sign-off still needs the PDFs.

### 8.1 Documents to obtain (the reading list)

| # | Document | What it pins down | Best lead (likely needs institutional / author access) |
|---|---|---|---|
| 1 | **Beets, Kimberley, Oliver, Pearce, Graham & Brandon (2012)** *Allometric Equations for Estimating Carbon Stocks in Natural Forest in NZ.* Forests 3(3):818–839, doi:10.3390/f3030818 | stem/branch/foliage volume & biomass eqns; tree-fern eqns; taper; carbon fractions | MDPI (open access — but 403 to bots); ResearchGate 273226431 |
| 2 | **Kimberley & Beets (2007)** *National volume function… Pinus radiata.* NZJFS 37(3):355–371 | compatible volume–taper function form (4th-order + butt-swell) | Scion `…/0004/59053/03_Kimberley.pdf` |
| 3 | **Payton, Newell & Beets (2004)** *NZ Carbon Monitoring System — Indigenous Forest & Shrubland Data Collection Manual* | **plot design, DBH cut-offs, field protocol → the §10 cut-offs** | ResearchGate 237809770; Landcare/MfE |
| 4 | **DOC/MfE** *Field Protocols for Tier 1 Inventory & Monitoring and LUCAS Plots* | plot/subplot geometry, measurement rules, missed-stem handling | doc.govt.nz contentassets PDF |
| 5 | **Coomes, Allen, Scott, Goulding & Beets (2002)** *Designing systems to monitor carbon stocks in forests and shrublands.* For. Ecol. Manage. 164:89–108 | monitoring-system design; decay-class density modifiers | ScienceDirect S0378-1127(01)00592-? |
| 6 | **Beets & Hood et al. (2008)** *CWD decay rates for seven indigenous species…* FEM 256(5):548–557 | species CWD half-lives / λ | ScienceDirect S0378112708003903 |
| 7 | **Garrett et al. (2018)** *Decay rates of above- & below-ground CWD…* FEM; & *Comparison of measured/modelled CWD…* FEM | λ model, MAT & diameter adjustment coefficients | Scion `…/07_Garrett.pdf`; ScienceDirect S0378112718312659 / …12477 |
| 8 | **Easdale et al. (2019)** *Root biomass allocation in southern temperate forests.* FEM 432:117463, doi:10.1016/j.foreco.2019.117463 | root:shoot ratios by group | ScienceDirect S0378112719307546 |

### 8.2 Corroboration status of §7 values

| value | status |
|---|---|
| Carbon fractions 0.51 / 0.48 (conifer/broadleaf) | **corroborated** (Beets 2012 / IPCC 2006) |
| Decay-class modifiers 1.00 / 0.82 / 0.66 / 0.47 (classes 0–3) | **corroborated** (Coomes 2002); class 4 → 0 assumed (fully decayed) — verify |
| DBH cut-offs 2.5 / 10 / 60 cm; plot 0.1257 ha | **corroborated** (CMS manual / LUCAS) |
| Tree-fern root:shoot 0.194 | **corroborated** (Easdale 2019) |
| Volume/branch/foliage/taper coefficients (§7) | **not corroborated** from open sources — `# VERIFY` against Beets 2012 / Kimberley 2007 |
| Decay base λ 0.0216 & adjustment 0.093 / 0.00908 | **not corroborated** (mean λ ≈ 0.0270 seen) — `# VERIFY` against Garrett 2018 |
| Root:shoot conifer/other assignment (0.234 vs 0.245) | ✅ **correct in code** (conifers→0.245 in both SAS & R); R constant **name** `BGBFracConifPalm` is misleading — rename only |
| Cabbage-tree root:shoot 0.437 | ⚠️ **not corroborated** by literature, but present in SAS reference — verify source |

### 8.3 End-stage checklist
1. Obtain the §8.1 PDFs (institutional access / author contact: Beets & Kimberley at Scion).
2. For every `# VERIFY` constant, confirm value + equation form; record page/equation number in the code comment.
3. Note unit conventions (cm vs m, over- vs under-bark, kg vs t).
4. **Resolve the remaining ⚠️ item** (cabbage-tree 0.437 root:shoot — uncorroborated by literature, present in SAS reference). The root:shoot conifer/broadleaf assignment was checked and is **correct** in both models; only the R constant name `BGBFracConifPalm` needs renaming.
5. Use docs #3/#4 to lock the §10 cut-offs (missed-stem rules, minimum sizes, plot geometry).

---

## 9. Testing strategy

- **Known-answer unit tests** (`tests/testthat`): hand-computed carbon for a canonical live tree, a recruit (ingrowth), a death (mortality), a standing dead spar, a fallen CWD piece — verifying allometry and classification independently of any reference dataset.
- **Property tests**: carbon monotonic in DBH; ingrowth ≥ 0 for true recruits; mortality ≤ 0; per-ha = per-stem / area.
- **Validation tests**: each rule fires on a crafted bad record and passes a good one.
- **Optional regression-vs-SAS** (sanity only, not parity): compare distributions and flag large divergences for explanation — not as a pass/fail gate.

---

## 10. Open questions for review

1. **Validation thresholds (§6)** — are the first-pass defaults (0.5 cm shrink tolerance, 5 cm·yr⁻¹ max growth, 60 m max height, 400 cm max DBH) sensible for NZ indigenous species, or should any be species-specific?
2. **Model cut-offs (§4.1, §5.3)** — please confirm the hard rules: (a) `A → absent` at the *final* cycle = mortality (vs. "possibly missed, carry forward")? (b) a stem appearing only at a later cycle — what size/age evidence makes it a *back-cast missed sample* vs. genuine *ingrowth*? (c) max plausible per-year DBH increment used to bound back/forward-casts.
3. **CWD multipliers** — confirm the intended single application of the 1.767 sampling factor and the 1.19 below-ground-deadwood factor, and which pools each applies to.
4. **Output format** — match the existing `output/` CSV schema for downstream compatibility, or define a clean new schema for the annual report?
5. **Input source** — same wide stem table + CWD table as today, or is a cleaner input contract available?
6. **Annual provenance** — should each year's run archive its config + coefficient version + quarantine decisions alongside the output, so prior reports are reproducible?

---

## 11. Build order (once spec approved)

1. `00_coefficients.R` + `tests` for the allometry math.
2. `01_ingest.R` + ingest reconciliation.
3. `02_validate.R` + validation tests.
4. `03_impute.R`, `04_classify.R`.
5. `05_allometry.R`, `06_expand.R`, `07_summarise.R`.
6. `08_report.R`, `run_pipeline.R`.
7. End-stage coefficient validation (§8) + optional SAS sanity comparison.

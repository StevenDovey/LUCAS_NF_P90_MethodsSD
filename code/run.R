#08.06.26 #5
# End-to-end carbon run. Reads the inputs from mfedata/, fills diameter gaps, predicts height, computes every carbon pool and stem volume, classifies the per-interval fluxes, expands per hectare, summarises to plot, and writes output/plot_carbon_summary.csv. The data-cleaning and validation tables go to output/diagnostics/.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
root <- getwd()

# Inputs (mfedata/): stems CSV, CWD CSV, species.csv, plot_summary.csv.
f_cwd     <- "MFESensitivity_CWD_20260209_reduced_dblV2.csv"
f_stems   <- "MFESensitivity_Stems_20260209_reduced_dbl_removed3.csv"
f_species <- "species.csv"
f_plots   <- "plot_summary.csv"


for (f in c("coefficients.R","allometry.R","classify.R","expand.R","predict_dbh.R",
            "predict_height.R","validate.R","validate_deep.R"))
  source(file.path(root, "R", f))
co <- carbon_coefficients()
# Read blanks as NA, not as empty strings, so a blank text cell is missing rather than a value that silently passes a comparison.
rd <- function(name) read.csv(file.path(root, "mfedata", name), stringsAsFactors = FALSE,
                              check.names = FALSE, na.strings = c("NA", ""))

# Dead-wood corrections. Both multipliers apply to both dead pools (standing dead and fallen): the line-intersect sampling correction and the below-ground dead-root uplift each scale both pools, the approach taken to carry the pool uncertainty. Set either to 1 to disable.
deadwood_factor <- co$cwd_sampling_multiplier * co$cwd_belowground_uplift  # 1.767 * 1.19

sp  <- rd(f_species); names(sp) <- make.names(names(sp))
ps  <- rd(f_plots)

# Per-plot inner subplot area and mean annual temperature, from the spreadsheet.
inner_new <- as.numeric(ps[["Inner plot area new method"]])
inner_old <- as.numeric(ps[["Inner plot area original method"]])
plot_meta <- data.frame(plot = ps$PLOT,
  inner_plot_area = ifelse(is.finite(inner_new), inner_new,
                    ifelse(is.finite(inner_old), inner_old, co$plot_area_inner_default_ha)),
  mat = ifelse(is.finite(as.numeric(ps$MAT)), as.numeric(ps$MAT), co$mat_default_c),
  stringsAsFactors = FALSE)

spA <- function(code, col) sp[[col]][match(code, sp$Species.code)]
# Cycle by a fixed date range on the record's own observation date, the documented national survey rounds. 2008 is a genuine gap between round 1 and round 2, and a date matching none of the three ranges gets no cycle, which the current data never triggers.
cycle_ranges <- list(
  `1` = as.Date(c("2002-01-01", "2007-12-31")),
  `2` = as.Date(c("2009-01-01", "2014-07-31")),
  `3` = as.Date(c("2014-08-01", "2024-12-31")))
assign_cycle <- function(obs_date) {
  out <- rep(NA_integer_, length(obs_date))
  for (cyc in names(cycle_ranges)) {
    r <- cycle_ranges[[cyc]]
    out[!is.na(obs_date) & obs_date >= r[1] & obs_date <= r[2]] <- as.integer(cyc)
  }
  out
}
woody <- c("canopy tree","subcanopy tree","shrub","unknown")
ferns <- c("tree fern","palm","cabbage tree")

# --- live stems -------------------------------------------------------------
# Status: alive to A, the dead states (including Unknown) to X, anything else to NA. The raw AliveState is kept so validation can flag a blank or unrecognised value and an Unknown-derived death, rather than the row disappearing before validation ever sees it.
stm <- rd(f_stems)
s <- data.frame(plot = stm$ParentPlotName, tag = stm$ItemID,
                status = ifelse(stm$AliveState == "Alive", "A",
                         ifelse(stm$AliveState %in% c("Dead","Not Found","Unknown"), "X", NA_character_)),
                alive_state_raw = stm$AliveState,
                species_code = stm$PreferredSpeciesCode, dbh = as.numeric(stm$DiameterValue),
                height = as.numeric(stm$Height), lean = as.numeric(stm$LeanAngle),
                decay = as.numeric(stm$DecayClass), subplot = stm$SubPlot,
                obs_date = as.Date(stm$PlotObsStartDate, "%d/%m/%Y"),
                stringsAsFactors = FALSE)
s$year <- as.integer(format(s$obs_date, "%Y"))
s$habit   <- spA(s$species_code, "Habit")
.ph <- spA(s$species_code, "Phylum")
s$phylum  <- ifelse(!is.na(.ph) & .ph == "C", "C", "P")   # unknown defaults to broadleaf
s$density <- as.numeric(spA(s$species_code, "Whole.stem.density"))
s$density[is.na(s$density)] <- co$wood_density_unknown_kg_m3   # documented unknown-species density
s <- merge(s, plot_meta, by = "plot", all.x = TRUE)
s$cycle <- assign_cycle(s$obs_date)
s <- s[!is.na(s$cycle), ]

# Validation: run the per-row and deep checks on the raw status, quarantine the error-level records before the carbon, and write the separated logs. The deep checks are warnings, so they inform review but do not exclude.
s$decay_class <- s$decay
.v  <- validate_stems(s)
s   <- .v$clean
.dv <- deep_validate(s, sp$Species.code)
.flags <- rbind(.v$flags, .dv)
dc <- file.path(root, "output", "diagnostics")
dir.create(dc, recursive = TRUE, showWarnings = FALSE)
for (.sv in c("error","warning","info"))
  write.csv(.flags[.flags$severity == .sv, ], file.path(dc, sprintf("validation_%ss.csv", .sv)), row.names = FALSE, na = "")
write.csv(.v$quarantine, file.path(dc, "validation_quarantine.csv"), row.names = FALSE, na = "")

# A stem alive at a later cycle is set alive at the earlier cycles, the field sampling correction for a stem missed or mis-recorded as dead.
.skey <- paste(s$plot, s$tag)
.max_alive <- ave(ifelse(s$status == "A", s$cycle, NA_integer_), .skey,
                  FUN = function(x) if (all(is.na(x))) NA_integer_ else max(x, na.rm = TRUE))
s$status[s$status != "A" & !is.na(.max_alive) & s$cycle < .max_alive] <- "A"
s$dbh_source <- ifelse(is.na(s$dbh), NA_character_, "measured")
s <- fill_dbh_gaps(s)
s <- predict_height_two_stage(s)

fern_ht <- tapply(s$height[s$habit %in% ferns], s$species_code[s$habit %in% ferns], function(x) mean(x, na.rm = TRUE))
s$hu <- s$height
wf <- s$habit %in% woody & is.na(s$hu); s$hu[wf] <- s$predht[wf]
ff <- s$habit %in% ferns & is.na(s$hu); s$hu[ff] <- fern_ht[s$species_code[ff]]
s$dmod <- unname(co$decay_class_density_modifier[as.character(s$decay)])
oh <- pmax(s$height, s$predht, na.rm = TRUE)

s$carbon_agb <- NA_real_; s$carbon_bgb <- NA_real_; s$carbon_spars <- NA_real_
s$live_vol <- NA_real_; s$dead_vol <- NA_real_
iw <- s$status == "A" & s$habit %in% woody & !is.na(s$dbh) & s$dbh > 0 & !is.na(s$hu)
s$carbon_agb[iw] <- live_woody_agb(s$dbh[iw], s$hu[iw], s$density[iw], s$phylum[iw], co)
s$live_vol[iw]   <- live_stem_volume(s$dbh[iw], s$hu[iw], co)
lf <- s$status == "A" & s$habit %in% ferns & !is.na(s$dbh) & s$dbh > 0 & !is.na(s$hu)
s$carbon_agb[lf] <- fern_aboveground_carbon(s$dbh[lf], s$hu[lf], co)
s$live_vol[lf]   <- fern_stem_volume(s$dbh[lf], s$hu[lf], co)
s$carbon_bgb <- belowground_carbon(s$carbon_agb, s$phylum, s$habit, co)
sw <- s$status == "X" & s$habit %in% woody & !is.na(s$dbh) & s$dbh >= 10 & is.finite(oh)
s$dead_vol[sw]     <- stem_volume_to_10cm(s$dbh[sw], oh[sw], co)
s$carbon_spars[sw] <- s$density[sw] * s$dmod[sw] * s$dead_vol[sw] * co$carbon_fraction_deadwood
sx <- s$status == "X" & s$habit %in% ferns & !is.na(s$dbh) & s$dbh >= 10 & is.finite(oh)
s$carbon_spars[sx] <- fern_aboveground_carbon(s$dbh[sx], oh[sx], co) * s$dmod[sx]
s$dead_vol[sx]     <- fern_stem_volume(s$dbh[sx], oh[sx], co)

# --- coarse woody debris (fallen), nested-subplot split ---------------------
cw  <- rd(f_cwd)
c2 <- data.frame(plot = cw$ParentPlotName, species_code = cw$PreferredSpeciesCode,
                 decay = as.numeric(cw$DecayClass), len = as.numeric(cw$ItemObsComponentLinearDimension),
                 le1 = as.numeric(cw$LargeEnd1), le2 = as.numeric(cw$LargeEnd2),
                 se1 = as.numeric(cw$SmallEnd1), se2 = as.numeric(cw$SmallEnd2), subplot = cw$SubPlot,
                 stature = cw$CWDStature,
                 obs_date = as.Date(cw$PlotObsStartDate, "%d/%m/%Y"),
                 stringsAsFactors = FALSE)
c2$led <- sqrt(c2$le1 * c2$le2); c2$sed <- sqrt(c2$se1 * c2$se2)
c2$density <- as.numeric(spA(c2$species_code, "Whole.stem.density"))
c2$density[is.na(c2$density)] <- co$wood_density_unknown_kg_m3
c2 <- merge(c2, plot_meta, by = "plot", all.x = TRUE)
c2$cycle <- assign_cycle(c2$obs_date)
c2 <- c2[!is.na(c2$cycle) & is.finite(c2$len), ]
# Stature is read but not yet used to differentiate carbon fraction; the distribution is written for review before any fraction is sourced and applied to it.
write.csv(as.data.frame(table(stature = c2$stature, cycle = c2$cycle, useNA = "ifany")),
          file.path(dc, "cwd_stature_distribution.csv"), row.names = FALSE, na = "")
c2$dmod <- unname(co$decay_class_density_modifier[as.character(c2$decay)])
fru <- function(d1, d2, L) (pi * L / 3) * ((d1/200)^2 + (d2/200)^2 + (d1/200)*(d2/200))
cyl <- function(d, L) pi * L * (d/200)^2
circ <- co$plot_area_circular_ha
ext  <- !is.na(c2$subplot) & c2$subplot == "EXT"
big_both <- !is.na(c2$led) & c2$led >= 60 & (is.na(c2$sed) | c2$sed >= 60)
straddle <- !is.na(c2$led) & !is.na(c2$sed) & c2$led >= 60 & c2$sed < 60
small    <- !big_both & !straddle
lenA <- ifelse(straddle, (c2$led - 60) / (c2$led - c2$sed) * c2$len, 0)
volA <- ifelse(straddle, fru(c2$led, 60, lenA),
        ifelse(big_both, ifelse(is.na(c2$sed), cyl(c2$led, c2$len), fru(c2$led, c2$sed, c2$len)), 0))
volB <- ifelse(straddle, fru(60, c2$sed, c2$len - lenA),
        ifelse(small, cwd_piece_volume(c2$led, c2$sed, c2$len), 0))
cf <- co$carbon_fraction_deadwood
cA <- c2$density * c2$dmod * volA * cf
cB <- c2$density * c2$dmod * volB * cf

# --- per-hectare stocks per cycle -------------------------------------------
s$area  <- sampling_area_ha(s$status, s$dbh, s$subplot, s$inner_plot_area, co)
s$agb_ha   <- carbon_per_ha(s$carbon_agb, s$area)
s$bgb_ha   <- carbon_per_ha(s$carbon_bgb, s$area)
s$spars_ha <- carbon_per_ha(s$carbon_spars, s$area) * deadwood_factor
s$lvol_ha  <- ifelse(is.na(s$area) | s$area <= 0, NA, s$live_vol / s$area)
s$dvol_ha  <- ifelse(is.na(s$area) | s$area <= 0, NA, s$dead_vol / s$area)
haA <- ifelse(volA > 0, cA / circ / 1000, 0)
haB <- ifelse(volB > 0 & !ext & c2$inner_plot_area > 0, cB / c2$inner_plot_area / 1000, 0)
c2$fallen_ha <- (haA + haB) * deadwood_factor

# Inner vs circular sampling scope, as additional columns alongside the combined figures. A stem's scope follows the area it was actually expanded by: the full circular plot for large stems, the inner subplot for the rest. Each scope column is the additive component of the combined total contributed by stems of that scope, and scope_inner + scope_circular reconciles exactly to the combined column.
s$scope <- ifelse(is.na(s$area), NA_character_,
            ifelse(abs(s$area - circ) < 1e-9, "circular", "inner"))
pools <- c("agb_ha","bgb_ha","spars_ha","lvol_ha","dvol_ha")
for (p in pools) {
  s[[paste0(p,"_inner")]]    <- ifelse(!is.na(s$scope) & s$scope == "inner",    s[[p]], NA_real_)
  s[[paste0(p,"_circular")]] <- ifelse(!is.na(s$scope) & s$scope == "circular", s[[p]], NA_real_)
}
c2$fallen_ha_circular <- haA * deadwood_factor
c2$fallen_ha_inner    <- haB * deadwood_factor

# Same inner/circular split, as a raw plot total in t C rather than a per-hectare density, so the area figure itself can be checked independently. carbon_agb/carbon_bgb/carbon_spars are in kg (expand.R), so /1000 alone gives t C, with no area divisor.
s$agb_tc   <- s$carbon_agb / 1000
s$bgb_tc   <- s$carbon_bgb / 1000
s$spars_tc <- s$carbon_spars / 1000 * deadwood_factor
tc_pools <- c("agb_tc","bgb_tc","spars_tc")
for (p in tc_pools) {
  s[[paste0(p,"_inner")]]    <- ifelse(!is.na(s$scope) & s$scope == "inner",    s[[p]], NA_real_)
  s[[paste0(p,"_circular")]] <- ifelse(!is.na(s$scope) & s$scope == "circular", s[[p]], NA_real_)
}
c2$fallen_tc_circular <- ifelse(volA > 0, cA / 1000, 0) * deadwood_factor
c2$fallen_tc_inner    <- ifelse(volB > 0 & !ext, cB / 1000, 0) * deadwood_factor

ag <- function(df, val) { a <- aggregate(df[[val]], by = df[c("plot","cycle")], FUN = function(x) sum(x, na.rm = TRUE)); names(a)[3] <- val; a }
stem_cols <- c(pools, paste0(rep(pools, each = 2), c("_inner","_circular")),
               paste0(rep(tc_pools, each = 2), c("_inner","_circular")))
stockL <- Reduce(function(a,b) merge(a,b,by=c("plot","cycle"),all=TRUE),
                 lapply(stem_cols, function(v) ag(s, v)))
stockL <- merge(stockL, ag(c2, "fallen_ha"), by = c("plot","cycle"), all = TRUE)
stockL <- merge(stockL, ag(c2, "fallen_ha_inner"), by = c("plot","cycle"), all = TRUE)
stockL <- merge(stockL, ag(c2, "fallen_ha_circular"), by = c("plot","cycle"), all = TRUE)
stockL <- merge(stockL, ag(c2, "fallen_tc_inner"), by = c("plot","cycle"), all = TRUE)
stockL <- merge(stockL, ag(c2, "fallen_tc_circular"), by = c("plot","cycle"), all = TRUE)
Wst <- reshape(stockL, idvar = "plot", timevar = "cycle", direction = "wide")

# --- per-interval fluxes ----------------------------------------------------
s$key <- paste(s$plot, s$tag)
mwide <- function(m) { x <- s[s$cycle == m, c("key","plot","carbon_agb","carbon_bgb","dbh","status","area")]
  names(x)[3:7] <- paste0(c("agb","bgb","dbh","status","area"), m); x }
W <- merge(merge(mwide(1), mwide(2)[ -2], by = "key", all = TRUE), mwide(3)[ -2], by = "key", all = TRUE)
W$plot <- s$plot[match(W$key, s$key)]
isC <- function(a) !is.na(a) & abs(a - co$plot_area_circular_ha) < 1e-9
inner <- plot_meta$inner_plot_area[match(W$plot, plot_meta$plot)]
W$denom <- ifelse(isC(W$area1) | isC(W$area2) | isC(W$area3), co$plot_area_circular_ha, inner)
th <- co$dbh_min_live_inner_cm
flux <- function(a_s, a_e, d_s, d_e) data.frame(
  growth = carbon_growth(a_s, a_e, d_s, th), mortality = carbon_mortality(a_s, a_e),
  ingrowth = carbon_ingrowth(a_s, a_e, d_s, d_e, TRUE, th), change = carbon_net_change(a_s, a_e))
ph <- function(x) ifelse(is.na(W$denom) | W$denom <= 0, NA, x / W$denom / 1000)
fa <- list(`1_2` = flux(W$agb1,W$agb2,W$dbh1,W$dbh2), `2_3` = flux(W$agb2,W$agb3,W$dbh2,W$dbh3), `1_3` = flux(W$agb1,W$agb3,W$dbh1,W$dbh3))
fb <- list(`1_2` = flux(W$bgb1,W$bgb2,W$dbh1,W$dbh2), `2_3` = flux(W$bgb2,W$bgb3,W$dbh2,W$dbh3), `1_3` = flux(W$bgb1,W$bgb3,W$dbh1,W$dbh3))
Fd <- data.frame(plot = W$plot)
for (nm in c("growth","mortality","ingrowth","change")) for (iv in c("1_2","2_3","1_3")) {
  Fd[[paste0("AGB_",nm,"_",iv)]] <- ph(fa[[iv]][[nm]]); Fd[[paste0("BGB_",nm,"_",iv)]] <- ph(fb[[iv]][[nm]]) }
Wfl <- aggregate(. ~ plot, data = Fd, FUN = function(x) sum(x, na.rm = TRUE), na.action = na.pass)

out <- merge(Wst, Wfl, by = "plot", all = TRUE)
for (k in 1:3) {
  out[[paste0("CWD_ha.",k)]]      <- rowSums(cbind(out[[paste0("spars_ha.",k)]], out[[paste0("fallen_ha.",k)]]), na.rm = TRUE)
  out[[paste0("stem_vol_ha.",k)]] <- rowSums(cbind(out[[paste0("lvol_ha.",k)]],  out[[paste0("dvol_ha.",k)]]),  na.rm = TRUE)
  for (sc in c("inner","circular")) {
    out[[paste0("CWD_ha_",sc,".",k)]]      <- rowSums(cbind(out[[paste0("spars_ha_",sc,".",k)]], out[[paste0("fallen_ha_",sc,".",k)]]), na.rm = TRUE)
    out[[paste0("stem_vol_ha_",sc,".",k)]] <- rowSums(cbind(out[[paste0("lvol_ha_",sc,".",k)]],  out[[paste0("dvol_ha_",sc,".",k)]]),  na.rm = TRUE)
    out[[paste0("CWD_tc_",sc,".",k)]]      <- rowSums(cbind(out[[paste0("spars_tc_",sc,".",k)]], out[[paste0("fallen_tc_",sc,".",k)]]), na.rm = TRUE)
  }
}
dir.create(file.path(root, "output"), showWarnings = FALSE)
write.csv(out, file.path(root, "output/plot_carbon_summary.csv"), row.names = FALSE, na = "")
# Per-stem diagnostics for the chart and validation-summary scripts.
write.csv(s[, c("plot","tag","cycle","year","species_code","habit","status",
                "dbh","dbh_source","height","predht","carbon_agb","carbon_bgb")],
          file.path(dc, "stem_diagnostics.csv"), row.names = FALSE, na = "")
cat(sprintf("written output/plot_carbon_summary.csv : %d plots, %d columns\n", nrow(out), ncol(out)))


# # # # # # # # # # # # # # # # #
# Run optional diagnostics. run in parallel. .
root <- getwd()
library(parallel)
parallel_scripts <- c("run_validate_post.R", "run_validation_summary.R", "mortality_audit.R",
                      "date_audit.R", "cwd_audit.R", "plotsummary_audit.R", "run_charts.R")
n_workers <- max(1, min(length(parallel_scripts), detectCores() - 1))
cl <- makeCluster(n_workers)
clusterExport(cl, "root")
results <- parLapply(cl, parallel_scripts, function(f) {
  setwd(root)
  capture.output(source(file.path(root, "diagnostics", f)))
})
stopCluster(cl)
for (i in seq_along(parallel_scripts))
  cat(sprintf("--- %s ---\n%s\n", parallel_scripts[i], paste(results[[i]], collapse = "\n")))
cat("--- stems_annotated.R ---\n")
source(file.path(root, "diagnostics", "stems_annotated.R"))

#07.06.26 21:45 NZST
# End-to-end AGB/BGB stock with the selected diameter gap fill and two-stage
# height model. Compares per-plot stock against the old R output and the
# reference plot summary.
root <- "/home/user/SAS2R"
for (f in c("coefficients.R","allometry.R","expand.R","predict_dbh.R","predict_height.R"))
  source(file.path(root, "R", f))
co <- carbon_coefficients()

stm <- read.csv(file.path(root, "mfedata/MFESensitivity_Stems_20260209_reduced_dbl_removed3.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
sp  <- read.csv(file.path(root, "mfedata/species.csv"), stringsAsFactors = FALSE, check.names = FALSE)
names(sp) <- make.names(names(sp))
ref <- read.csv(file.path(root, "dev/compare/plotsummary_3Cycles_LUCAS_plotsV2_R.csv"), check.names = FALSE)
ref <- read.csv(file.path(root, "cmp2/SAS_plotsummary_CWD_multiV10.csv"), check.names = FALSE)
ref <- ref[!is.na(ref$plot) & ref$plot != "", ]

s <- data.frame(plot = stm$ParentPlotName, tag = stm$ItemID,
                date = stm$PlotObsStartDate, status = ifelse(stm$AliveState == "Alive", "A", "X"),
                species_code = stm$PreferredSpeciesCode, dbh = as.numeric(stm$DiameterValue),
                height = as.numeric(stm$Height), lean = as.numeric(stm$LeanAngle),
                subplot = stm$SubPlot, stringsAsFactors = FALSE)
s$habit   <- sp$Habit[match(s$species_code, sp$Species.code)]
s$phylum  <- ifelse(sp$Phylum[match(s$species_code, sp$Species.code)] == "C", "C", "P")
s$density <- as.numeric(sp$Whole.stem.density[match(s$species_code, sp$Species.code)])
s$density[is.na(s$density)] <- co$wood_density_unknown_kg_m3
s$year <- as.integer(format(as.Date(s$date, "%d/%m/%Y"), "%Y"))
ry <- ref[, c("plot","plotyear1","plotyear2","plotyear3","inner_plot_area")]
s <- merge(s, ry, by = "plot", all.x = TRUE)
s$cycle <- with(s, ifelse(year == plotyear1, 1L, ifelse(year == plotyear2, 2L,
                   ifelse(year == plotyear3, 3L, NA_integer_))))
s <- s[!is.na(s$cycle), ]
s$dbh_source <- ifelse(is.na(s$dbh), NA_character_, "measured")

s <- fill_dbh_gaps(s)
s <- predict_height_two_stage(s)

ferns <- c("tree fern","palm","cabbage tree")
fern_mean <- tapply(s$height[s$habit %in% ferns], s$species_code[s$habit %in% ferns],
                    function(x) mean(x, na.rm = TRUE))
s$height_use <- s$height
wf <- s$habit %in% c("canopy tree","subcanopy tree","shrub","unknown") & is.na(s$height_use)
s$height_use[wf] <- s$predht[wf]
ff <- s$habit %in% ferns & is.na(s$height_use)
s$height_use[ff] <- fern_mean[s$species_code[ff]]

live <- s$status == "A" & !is.na(s$dbh) & s$dbh > 0 & !is.na(s$height_use)
s$carbon_agb <- NA_real_
iw <- live & s$habit %in% c("canopy tree","subcanopy tree","shrub","unknown")
s$carbon_agb[iw] <- live_woody_agb(s$dbh[iw], s$height_use[iw], s$density[iw], s$phylum[iw], co)
fl <- live & s$habit %in% ferns
s$carbon_agb[fl] <- fern_aboveground_carbon(s$dbh[fl], s$height_use[fl], co)
s$carbon_bgb <- belowground_carbon(s$carbon_agb, s$phylum, s$habit, co)

s$area  <- sampling_area_ha(s$status, s$dbh, s$subplot, s$inner_plot_area, co)
s$agb_ha <- carbon_per_ha(s$carbon_agb, s$area)
s$bgb_ha <- carbon_per_ha(s$carbon_bgb, s$area)

agg <- aggregate(cbind(agb_ha, bgb_ha) ~ plot + cycle, data = s, FUN = function(x) sum(x, na.rm = TRUE))
w <- reshape(agg, idvar = "plot", timevar = "cycle", direction = "wide")

report <- function(mine, theirs, label) {
  ok <- is.finite(mine) & is.finite(theirs) & theirs > 0
  cat(sprintf("%-22s n=%4d  median ratio %.3f  corr %.3f\n",
      label, sum(ok), median(mine[ok]/theirs[ok]), cor(mine[ok], theirs[ok])))
}
cR <- merge(w, ref[, c("plot","carbon_AGB_ha1","carbon_AGB_ha2","carbon_AGB_ha3","carbon_BGB_ha1")], by = "plot")
cS <- merge(w, ref[, c("plot","carbon_AGB_ha1","carbon_AGB_ha2","carbon_AGB_ha3")], by = "plot")
cat("vs R reference:\n")
report(cR$agb_ha.1, cR$carbon_AGB_ha1, "  AGB cycle 1")
report(cR$agb_ha.2, cR$carbon_AGB_ha2, "  AGB cycle 2")
report(cR$agb_ha.3, cR$carbon_AGB_ha3, "  AGB cycle 3")
report(cR$bgb_ha.1, cR$carbon_BGB_ha1, "  BGB cycle 1")
cat("vs reference plot summary:\n")
report(cS$agb_ha.1, cS$carbon_AGB_ha1, "  AGB cycle 1")
report(cS$agb_ha.2, cS$carbon_AGB_ha2, "  AGB cycle 2")
report(cS$agb_ha.3, cS$carbon_AGB_ha3, "  AGB cycle 3")

# Which is closer to the reference: old R or the new model, on identical plots.
names(ref)[names(ref) %in% c("carbon_AGB_ha1","carbon_AGB_ha2","carbon_AGB_ha3")] <-
  c("ref1","ref2","ref3")
names(ref)[names(ref) %in% c("carbon_AGB_ha1","carbon_AGB_ha2","carbon_AGB_ha3")] <-
  c("old1","old2","old3")
T <- merge(merge(w, ref[, c("plot","old1","old2","old3")], by = "plot"),
           ref[, c("plot","ref1","ref2","ref3")], by = "plot")
closer <- function(new, old, ref, label) {
  ok <- is.finite(new) & is.finite(old) & is.finite(ref) & ref > 0
  en <- abs(new[ok] - ref[ok]); eo <- abs(old[ok] - ref[ok])
  cat(sprintf("%-12s n=%4d | median|new-reference| %6.3f  median|old-reference| %6.3f | new closer %4.1f%% | RMSE new %6.3f old %6.3f\n",
      label, sum(ok), median(en), median(eo), 100*mean(en < eo),
      sqrt(mean(en^2)), sqrt(mean(eo^2))))
}
cat("\nCloser to reference (AGB t C/ha), new model vs old R, identical plots:\n")
closer(T$agb_ha.1, T$old1, T$ref1, "cycle 1")
closer(T$agb_ha.2, T$old2, T$ref2, "cycle 2")
closer(T$agb_ha.3, T$old3, T$ref3, "cycle 3")

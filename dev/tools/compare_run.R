# Rough first-light run on the real MFE data, checked against the reference
# plot-summary output.
#   Faithful:     allometry (R/allometry.R), per-plot inner area and cycle years
#                 (taken from the reference file), real DBH/species/status.
# Approximated: height for the roughly 78% of live stems with no measured
#                 height (a simple per-species log-log model plus a plot by
#                 cycle bias correction, not the full mixed model); DBH gaps
#                 not imputed; CWD, dead-wood and fluxes not computed.
# So systematic differences are expected; this checks order of magnitude and
# correlation of AGB/BGB stock, not an exact match.

root <- "/home/user/SAS2R"
source(file.path(root, "R/coefficients.R"))
source(file.path(root, "R/allometry.R"))
source(file.path(root, "R/expand.R"))
co <- carbon_coefficients()

# --- Load -------------------------------------------------------------------
S  <- file.path(root, "mfedata/MFESensitivity_Stems_20260209_reduced_dbl_removed3.csv")
d  <- read.csv(S, stringsAsFactors = FALSE, check.names = FALSE)
sp <- read.csv(file.path(root, "mfedata/species.csv"), stringsAsFactors = FALSE, check.names = FALSE)
names(sp) <- make.names(names(sp))
ref <- read.csv(file.path(root, "dev/compare/plotsummary_3Cycles_LUCAS_plotsV2_R.csv"), check.names = FALSE)

# --- Canonical stem frame ---------------------------------------------------
st <- data.frame(
  plot   = d$ParentPlotName,
  tag    = d$ItemID,
  date   = d$PlotObsStartDate,
  status = ifelse(d$AliveState == "Alive", "A",
            ifelse(d$AliveState %in% c("Dead","Not Found","Unknown",""), "X", "X")),
  species_code = d$PreferredSpeciesCode,
  dbh    = suppressWarnings(as.numeric(d$DiameterValue)),
  height = suppressWarnings(as.numeric(d$Height)),
  lean   = suppressWarnings(as.numeric(d$LeanAngle)),
  subplot = d$SubPlot,
  stringsAsFactors = FALSE)

# species attributes
st$habit   <- sp$Habit[match(st$species_code, sp$Species.code)]
st$phylum  <- ifelse(sp$Phylum[match(st$species_code, sp$Species.code)] == "C", "C", "P")
st$density <- suppressWarnings(as.numeric(sp$Whole.stem.density[match(st$species_code, sp$Species.code)]))
st$density[is.na(st$density)] <- co$wood_density_unknown_kg_m3

# observation year, and per-plot cycle years from the reference
st$year <- as.integer(format(as.Date(st$date, "%d/%m/%Y"), "%Y"))
ry <- ref[, c("plot","plotyear1","plotyear2","plotyear3","inner_plot_area")]
st <- merge(st, ry, by = "plot", all.x = TRUE)
st$cycle <- with(st, ifelse(year == plotyear1, 1L,
                     ifelse(year == plotyear2, 2L,
                     ifelse(year == plotyear3, 3L, NA_integer_))))
st <- st[!is.na(st$cycle), ]

# --- Height model (approximate) ---------------------------------------------
woody <- c("canopy tree","subcanopy tree","shrub","unknown")
st$is_woody <- st$habit %in% woody
tr <- st[st$status=="A" & st$is_woody & !is.na(st$height) & st$height>1.35 &
           is.na(st$lean) & !is.na(st$dbh) & st$dbh>0, ]
tr$lny <- log(tr$height - 1.35); tr$x3 <- tr$dbh^(-0.3)
glob <- lm(lny ~ x3, data = tr)
# per-species slopes/intercepts where enough data
fit_species <- function(g) if (nrow(g) >= 10) coef(lm(lny ~ x3, data=g)) else NULL
sp_models <- lapply(split(tr, tr$species_code), fit_species)
pred_lny <- function(species, x3) {
  m <- sp_models[[as.character(species[1])]]
  if (is.null(m)) m <- coef(glob)
  m[1] + m[2]*x3
}
st$x3 <- st$dbh^(-0.3)
st$predht <- NA_real_
iw <- st$is_woody & !is.na(st$dbh) & st$dbh>0
# predict by species group
for (s in unique(st$species_code[iw])) {
  sel <- iw & st$species_code==s
  m <- sp_models[[as.character(s)]]; if (is.null(m)) m <- coef(glob)
  st$predht[sel] <- 1.35 + exp(pmin(m[1] + m[2]*st$x3[sel], log(80)))
}
# plot x cycle bias correction toward measured heights
st$h135 <- st$height - 1.35
bias_tab <- aggregate(cbind(meas=h135, pred=predht-1.35) ~ plot+cycle,
                      data=st[st$status=="A" & st$is_woody & !is.na(st$height),],
                      FUN=function(x) mean(x,na.rm=TRUE))
bias_tab$bias <- bias_tab$meas / bias_tab$pred
st <- merge(st, bias_tab[,c("plot","cycle","bias")], by=c("plot","cycle"), all.x=TRUE)
st$bias[is.na(st$bias) | !is.finite(st$bias)] <- 1
# final height: measured where present, else bias-corrected prediction
st$height_use <- ifelse(!is.na(st$height) & st$height>1.35, st$height,
                        1.35 + (st$predht-1.35)*st$bias)

# --- Per-stem carbon (live only; AGB + BGB) ---------------------------------
ferns <- c("tree fern","palm","cabbage tree")
live  <- st$status=="A" & !is.na(st$dbh) & st$dbh>0 & !is.na(st$height_use)
st$carbon_agb <- NA_real_
iw2 <- live & st$is_woody
st$carbon_agb[iw2] <- live_woody_agb(st$dbh[iw2], st$height_use[iw2], st$density[iw2], st$phylum[iw2], co)
iff <- live & st$habit %in% ferns
st$carbon_agb[iff] <- fern_aboveground_carbon(st$dbh[iff], st$height_use[iff], co)
st$carbon_bgb <- belowground_carbon(st$carbon_agb, st$phylum, st$habit, co)

# --- Per-hectare expansion --------------------------------------------------
st$area <- sampling_area_ha(st$status, st$dbh, st$subplot, st$inner_plot_area, co)
st$agb_ha <- carbon_per_ha(st$carbon_agb, st$area)
st$bgb_ha <- carbon_per_ha(st$carbon_bgb, st$area)

# --- Plot x cycle stock -----------------------------------------------------
agg <- aggregate(cbind(agb_ha, bgb_ha) ~ plot+cycle, data=st, FUN=function(x) sum(x,na.rm=TRUE))
wide <- reshape(agg, idvar="plot", timevar="cycle", direction="wide")

cmp <- merge(ref[,c("plot","carbon_AGB_ha1","carbon_AGB_ha2","carbon_AGB_ha3",
                    "carbon_BGB_ha1","carbon_BGB_ha2","carbon_BGB_ha3")],
             wide, by="plot", all.x=TRUE)

report <- function(label, mine, theirs) {
  ok <- is.finite(mine) & is.finite(theirs) & theirs>0
  cat(sprintf("\n%s  (n=%d plots)\n", label, sum(ok)))
  cat(sprintf("  mine   mean %7.2f   theirs mean %7.2f\n", mean(mine[ok]), mean(theirs[ok])))
  cat(sprintf("  ratio mine/theirs  median %.3f   mean %.3f\n",
      median((mine[ok])/(theirs[ok])), mean((mine[ok])/(theirs[ok]))))
  cat(sprintf("  correlation %.3f\n", cor(mine[ok], theirs[ok])))
}
report("AGB stock cycle1", cmp$agb_ha.1, cmp$carbon_AGB_ha1)
report("AGB stock cycle2", cmp$agb_ha.2, cmp$carbon_AGB_ha2)
report("AGB stock cycle3", cmp$agb_ha.3, cmp$carbon_AGB_ha3)
report("BGB stock cycle1", cmp$bgb_ha.1, cmp$carbon_BGB_ha1)

cat(sprintf("\nstems used (live, carbon computed): %d\n", sum(!is.na(st$carbon_agb))))
cat(sprintf("heights: measured %.1f%%, modelled %.1f%%\n",
    100*mean(!is.na(st$height[live])), 100*mean(is.na(st$height[live]))))

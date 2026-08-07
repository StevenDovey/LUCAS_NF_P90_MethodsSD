# Validates the flux pathway (growth, mortality, ingrowth, change) against the
# reference plot-summary output. Takes the reference's own per-stem carbon and
# diameter by measurement, pivots to wide per stem, applies the interval logic
# from R/classify.R and the shared area denominator from R/expand.R, sums to
# plot, and compares to the reference plot summary. This isolates the interval
# logic and the ingrowth NA gate from the allometry (already exact) and the
# height model.

root <- "/home/user/SAS2R"
source(file.path(root, "R/coefficients.R"))
source(file.path(root, "R/classify.R"))
source(file.path(root, "R/expand.R"))
co <- carbon_coefficients()

M <- file.path(root, "sasout/SAS_measurements-with-carbonV10.csv")
keep <- c("plot","tag","ItemID","measurement","carbon_AGB","carbon_BGB",
          "DBH","status","subplot","Inner_plot_area")
d <- read.csv(M, stringsAsFactors = FALSE, check.names = FALSE)[keep]
for (c in c("measurement","carbon_AGB","carbon_BGB","DBH","Inner_plot_area"))
  d[[c]] <- suppressWarnings(as.numeric(d[[c]]))
d$key <- paste(d$plot, d$ItemID)

# --- pivot to wide per stem (merge per measurement) -------------------------
d <- d[!is.na(d$measurement) & d$measurement %in% 1:3, ]
# collapse any duplicate (key, measurement) by taking the first
d <- d[!duplicated(paste(d$key, d$measurement)), ]
mwide <- function(m) {
  s <- d[d$measurement == m, c("key","plot","carbon_AGB","carbon_BGB","DBH","status","subplot","Inner_plot_area")]
  names(s) <- c("key","plot", paste0(c("agb","bgb","dbh","status","subplot","inner"), m))
  s
}
w <- mwide(1)
w <- merge(w, mwide(2)[setdiff(names(mwide(2)),"plot")], by = "key", all = TRUE)
w <- merge(w, mwide(3)[setdiff(names(mwide(3)),"plot")], by = "key", all = TRUE)
# plot + stable subplot/inner across cycles
w$plot   <- d$plot[match(w$key, d$key)]
w$subplot <- ifelse(!is.na(w$subplot1), w$subplot1,
              ifelse(!is.na(w$subplot2), w$subplot2, w$subplot3))
w$inner  <- ifelse(!is.na(w$inner1), w$inner1, ifelse(!is.na(w$inner2), w$inner2, w$inner3))
for (m in 1:3) w[[paste0("status",m)]][is.na(w[[paste0("status",m)]])] <- ""

# --- shared area denominator (circular if any cycle circular, else inner) ---
area1 <- sampling_area_ha(ifelse(is.na(w$status1),"",w$status1), w$dbh1, w$subplot, w$inner, co)
area2 <- sampling_area_ha(ifelse(is.na(w$status2),"",w$status2), w$dbh2, w$subplot, w$inner, co)
area3 <- sampling_area_ha(ifelse(is.na(w$status3),"",w$status3), w$dbh3, w$subplot, w$inner, co)
is_circ <- function(a) !is.na(a) & abs(a - co$plot_area_circular_ha) < 1e-9
w$denom <- ifelse(is_circ(area1) | is_circ(area2) | is_circ(area3),
                  co$plot_area_circular_ha,
                  ifelse(!is.na(w$inner), w$inner, NA))

# --- my interval logic (AGB) ; visited_start = plot measured at start (TRUE) -
mk <- function(s, e, dse, de) {
  data.frame(
    growth    = carbon_growth(s, e, dse, threshold_cm = co$dbh_min_live_inner_cm),
    mortality = carbon_mortality(s, e),
    ingrowth  = carbon_ingrowth(s, e, dse, de, visited_start = TRUE,
                                threshold_cm = co$dbh_min_live_inner_cm),
    change    = carbon_net_change(s, e))
}
f12 <- mk(w$agb1, w$agb2, w$dbh1, w$dbh2)
f23 <- mk(w$agb2, w$agb3, w$dbh2, w$dbh3)
f13 <- mk(w$agb1, w$agb3, w$dbh1, w$dbh3)

perha <- function(x) ifelse(is.na(w$denom) | w$denom <= 0, NA, x / w$denom / 1000)
mine <- data.frame(plot = w$plot,
  g12 = perha(f12$growth), m12 = perha(f12$mortality), i12 = perha(f12$ingrowth), c12 = perha(f12$change),
  g23 = perha(f23$growth), m23 = perha(f23$mortality), i23 = perha(f23$ingrowth), c23 = perha(f23$change),
  g13 = perha(f13$growth), m13 = perha(f13$mortality), i13 = perha(f13$ingrowth), c13 = perha(f13$change))
agg <- aggregate(. ~ plot, data = mine, FUN = function(x) sum(x, na.rm = TRUE), na.action = na.pass)

# --- reference plot summary -------------------------------------------------------
ref <- read.csv(file.path(root, "cmp2/SAS_plotsummary_CWD_multiV10.csv"), check.names = FALSE)
ref <- ref[!is.na(ref$plot) & ref$plot != "", ]
# reference "growth+ingrowth" combines recruits the way mine splits them; compare both.
ref$gi12 <- ref$carbon_AGB_ha_growth1_2 + ref$carbon_AGB_ha_ingrowth1_2
ref$gi23 <- ref$carbon_AGB_ha_growth2_3 + ref$carbon_AGB_ha_ingrowth2_3
ref$gi13 <- ref$carbon_AGB_ha_growth1_3 + ref$carbon_AGB_ha_ingrowth1_3
agg$gi12 <- agg$g12 + agg$i12; agg$gi23 <- agg$g23 + agg$i23; agg$gi13 <- agg$g13 + agg$i13
ref$chg12 <- ref$carbon_AGB_ha_change1_2; ref$chg23 <- ref$carbon_AGB_ha_change2_3; ref$chg13 <- ref$carbon_AGB_ha_change1_3

pairs <- list(
  c("m12","carbon_AGB_ha_mortality1_2"), c("m23","carbon_AGB_ha_mortality2_3"), c("m13","carbon_AGB_ha_mortality1_3"),
  c("gi12","gi12"), c("gi23","gi23"), c("gi13","gi13"),
  c("c12","chg12"), c("c23","chg23"), c("c13","chg13"))
need <- unique(c(sapply(pairs, `[`, 2)))
cmp <- merge(agg, ref[, c("plot", intersect(need, names(ref)))], by = "plot")

cat(sprintf("plots compared: %d\n\n", nrow(cmp)))
cat("Disjoint-partition check vs reference (mortality; growth+ingrowth combined; net change)\n")
cat(sprintf("%-22s %8s %8s %8s %7s\n","quantity (AGB t C/ha)","mine","reference","meandiff","corr"))
lab <- c(m12="mortality 1-2", m23="mortality 2-3", m13="mortality 1-3",
         gi12="growth+ingrowth 1-2", gi23="growth+ingrowth 2-3", gi13="growth+ingrowth 1-3",
         c12="net change 1-2", c23="net change 2-3", c13="net change 1-3")
for (p in pairs) {
  a <- suppressWarnings(as.numeric(cmp[[p[1]]])); b <- suppressWarnings(as.numeric(cmp[[p[2]]]))
  ok <- is.finite(a) & is.finite(b)
  cat(sprintf("%-22s %8.3f %8.3f %8.4f %7.3f\n", lab[[p[1]]],
      mean(a[ok]), mean(b[ok]), mean(a[ok]-b[ok]), suppressWarnings(cor(a[ok],b[ok]))))
}

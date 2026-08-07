#08.06.26 #2
# Confidence charts, run after run.R. Reads the per-stem diagnostics and the plot summary and writes one PNG per chart to output/charts/. Per-stem charts are also split into per-species pages. Cycle is shown by colour.
# Run from the repository root, after run.R:  Rscript diagnostics/run_charts.R
library(nlme)
root <- getwd()
d  <- read.csv(file.path(root, "output/diagnostics/stem_diagnostics.csv"), stringsAsFactors = FALSE, check.names = FALSE)
ps <- read.csv(file.path(root, "output/plot_carbon_summary.csv"), check.names = FALSE)
dir.create(file.path(root, "output/charts"), recursive = TRUE, showWarnings = FALSE)

cf   <- function(name) file.path(root, "output/charts", name)
full <- function(f) png(f, width = 15, height = 10, units = "cm", res = 300, pointsize = 8)
half <- function(f) png(f, width = 7.5, height = 5, units = "cm", res = 300, pointsize = 6)
cyc_col   <- c(`1` = "#1b9e77", `2` = "#d95f02", `3` = "#7570b3")
habit_pal <- c(`canopy tree` = "#1b9e77", `subcanopy tree` = "#d95f02", shrub = "#7570b3",
               `tree fern` = "#e7298a", palm = "#66a61e", `cabbage tree` = "#e6ab02", unknown = "#999999")

scatter11 <- function(obs, pred, col, xlab, ylab, main) {
  ok <- is.finite(obs) & is.finite(pred)
  plot(obs[ok], pred[ok], pch = 16, col = col[ok], cex = 0.4, xlab = xlab, ylab = ylab, main = main)
  abline(0, 1, col = "red", lwd = 1.5)
  bias <- mean(pred[ok] - obs[ok]); rmse <- sqrt(mean((pred[ok] - obs[ok])^2))
  mtext(sprintf("n = %d   bias = %.2f   RMSE = %.2f", sum(ok), bias, rmse), side = 3, cex = 0.6)
}

# Per-species pages, coloured by cycle, paged. The reference line is supplied by the caller: the identity line for predicted-against-observed panels, a fitted curve where the two axes are different quantities (height against diameter).
ref_identity <- function(x, y) abline(0, 1, col = "red", lwd = 1)
ref_height_diameter <- function(x, y) {
  b <- coef(lm(log(y - 1.35) ~ log(x)))
  xs <- seq(min(x), max(x), length.out = 100)
  lines(xs, 1.35 + exp(b[1] + b[2] * log(xs)), col = "red", lwd = 1.5)
}
facet_species <- function(df, xcol, ycol, prefix, xlab, ylab, refline = ref_identity,
                          per_page = 16L, min_n = 30L) {
  tab <- sort(table(df$species_code), decreasing = TRUE)
  sps <- names(tab)[tab >= min_n & nzchar(names(tab))]
  if (!length(sps)) return(invisible())
  nc <- 4L; nr <- ceiling(per_page / nc); np <- ceiling(length(sps) / per_page)
  for (pg in seq_len(np)) {
    these <- sps[((pg - 1L) * per_page + 1L):min(pg * per_page, length(sps))]
    full(cf(sprintf("%s_%02d.png", prefix, pg)))
    par(mfrow = c(nr, nc), mar = c(3, 3, 2, 0.6), mgp = c(1.6, 0.5, 0), cex.main = 0.8)
    for (s in these) {
      g <- df[df$species_code == s, ]; x <- g[[xcol]]; y <- g[[ycol]]; ok <- is.finite(x) & is.finite(y)
      plot(x[ok], y[ok], pch = 16, col = cyc_col[as.character(g$cycle[ok])], cex = 0.5,
           xlab = xlab, ylab = ylab, main = sprintf("%s (n=%d)", s, sum(ok)))
      refline(x[ok], y[ok])
    }
    dev.off()
  }
}

# --- height predicted vs observed -------------------------------------------
hm <- d[!is.na(d$height) & d$height > 1.35 & !is.na(d$predht) & !is.na(d$dbh) & d$dbh > 0, ]
woody_h <- c("canopy tree","subcanopy tree","shrub","unknown")
trng <- quantile(hm$dbh[hm$habit %in% woody_h], c(0.025, 0.975), na.rm = TRUE)
hcol <- habit_pal[hm$habit]; hcol[is.na(hcol)] <- "#999999"
core <- hm$dbh >= trng[1] & hm$dbh <= trng[2]
full(cf("height_pred_vs_obs_core.png"))
scatter11(hm$height[core], hm$predht[core], paste0(hcol[core], "99"),
          "observed height (m)", "predicted height (m)",
          sprintf("Height predicted vs observed, within training diameter (%.0f-%.0f cm)", trng[1], trng[2]))
legend("bottomright", legend = names(habit_pal), col = habit_pal, pch = 16, cex = 0.6, bty = "n")
dev.off()
full(cf("height_pred_vs_obs_tail.png"))
scatter11(hm$height[!core], hm$predht[!core], paste0(hcol[!core], "99"),
          "observed height (m)", "predicted height (m)", "Height predicted vs observed, diameter tail (sparse data)")
legend("bottomright", legend = names(habit_pal), col = habit_pal, pch = 16, cex = 0.6, bty = "n")
dev.off()
facet_species(hm, "height", "predht", "height_pred_vs_obs_species", "observed height (m)", "predicted height (m)")

# height residual vs predicted (bias structure)
full(cf("height_residual.png"))
plot(hm$predht, hm$predht - hm$height, pch = 16, col = "#33669944", cex = 0.4,
     xlab = "predicted height (m)", ylab = "residual: predicted minus observed (m)",
     main = "Height residual vs predicted")
abline(h = 0, col = "red", lwd = 1.5)
dev.off()

# height residual by habit (compare growth forms)
hm$resid <- hm$predht - hm$height
full(cf("height_residual_by_habit.png"))
par(mar = c(7, 4, 3, 1))
boxplot(resid ~ habit, data = hm, las = 2, cex.axis = 0.7, ylab = "residual: predicted minus observed (m)",
        xlab = "", main = "Height residual by habit", col = "#99bbaa")
abline(h = 0, col = "red", lwd = 1.5)
dev.off()

# mean bias by species (worst species rank)
bsp <- tapply(hm$predht - hm$height, hm$species_code, mean)
nsp <- tapply(hm$predht - hm$height, hm$species_code, length)
bsp <- bsp[nsp >= 30]; bsp <- bsp[order(abs(bsp), decreasing = TRUE)][1:min(25, length(bsp))]
half(cf("height_bias_by_species.png"))
par(mar = c(3, 7, 3, 1))
barplot(rev(bsp), horiz = TRUE, las = 1, cex.names = 0.5, xlab = "mean height bias (m)",
        main = "Height bias by species (worst 25)", col = "#cc7755")
abline(v = 0, col = "red")
dev.off()

# --- diameter predicted vs observed (model on measured consecutive pairs) ----
m <- d[d$status == "A" & !is.na(d$dbh) & d$dbh > 0 & d$dbh_source == "measured" & nzchar(d$species_code), ]
pair <- function(a, b) {
  ea <- m[m$cycle == a, ]; eb <- m[m$cycle == b, ]
  j <- merge(ea[c("plot","tag","dbh","species_code","year")], eb[c("plot","tag","dbh","year")], by = c("plot","tag"), suffixes = c("_f","_t"))
  data.frame(from = j$dbh_f, to = j$dbh_t, species_code = j$species_code, cycle = b,
             dt = j$year_t - j$year_f, stringsAsFactors = FALSE)
}
P <- rbind(pair(1, 2), pair(2, 3))
P <- P[is.finite(P$from) & is.finite(P$to) & nzchar(P$species_code), ]
fit <- lme(to ~ from, random = ~1 | species_code, data = P, control = lmeControl(opt = "optim"))
re <- ranef(fit)[[1]]; names(re) <- rownames(ranef(fit)); b <- fixef(fit)
r <- re[P$species_code]; r[is.na(r)] <- 0
P$pred <- b[1] + b[2] * P$from + r
full(cf("dbh_pred_vs_obs.png"))
scatter11(P$to, P$pred, paste0(cyc_col[as.character(P$cycle)], "99"),
          "observed diameter (cm)", "predicted diameter (cm)", "Diameter predicted vs observed")
legend("bottomright", legend = paste("to cycle", names(cyc_col)), col = cyc_col, pch = 16, cex = 0.6, bty = "n")
dev.off()
facet_species(P, "to", "pred", "dbh_pred_vs_obs_species", "observed diameter (cm)", "predicted diameter (cm)")

# diameter increment per year (flags implausible growth and shrinkage)
inc <- (P$to - P$from) / P$dt; inc <- inc[is.finite(inc)]
half(cf("dbh_increment_per_year.png"))
hist(inc[inc > -3 & inc < 8], breaks = 60, xlab = "diameter increment (cm per year)",
     main = "Diameter increment per year", col = "#7799bb")
abline(v = 0, col = "red", lwd = 1.5); abline(v = 5, col = "red", lty = 2)
mtext(sprintf("below 0: %d   above 5: %d", sum(inc < 0), sum(inc > 5)), side = 3, cex = 0.6)
dev.off()

# --- height to diameter relationship, per species only ----------------------
hd <- d[!is.na(d$height) & d$height > 1.35 & !is.na(d$dbh) & d$dbh > 0, ]
facet_species(hd, "dbh", "height", "height_diameter_species", "diameter (cm)", "height (m)",
              refline = ref_height_diameter)

# --- plot-level results (separate by cycle) ---------------------------------
half(cf("agb_stock_by_cycle.png"))
boxplot(list(`cycle 1` = ps$agb_ha.1, `cycle 2` = ps$agb_ha.2, `cycle 3` = ps$agb_ha.3),
        ylab = "above-ground carbon (t C/ha)", main = "Plot stock by cycle", col = "#88bb88")
dev.off()
for (iv in c("1_2", "2_3")) {
  half(cf(sprintf("agb_change_%s.png", iv)))
  x <- ps[[paste0("AGB_change_", iv)]]; x <- x[is.finite(x)]
  hist(x, breaks = 40, xlab = "net above-ground change (t C/ha)",
       main = sprintf("Net change per plot, interval %s", sub("_", " to ", iv)), col = "#cc9966")
  abline(v = 0, col = "red", lwd = 1.5)
  dev.off()
}
pm <- rbind(
  `above-ground`  = c(mean(ps$agb_ha.1, na.rm=TRUE), mean(ps$agb_ha.2, na.rm=TRUE), mean(ps$agb_ha.3, na.rm=TRUE)),
  `below-ground`  = c(mean(ps$bgb_ha.1, na.rm=TRUE), mean(ps$bgb_ha.2, na.rm=TRUE), mean(ps$bgb_ha.3, na.rm=TRUE)),
  `standing dead` = c(mean(ps$spars_ha.1, na.rm=TRUE), mean(ps$spars_ha.2, na.rm=TRUE), mean(ps$spars_ha.3, na.rm=TRUE)),
  fallen          = c(mean(ps$fallen_ha.1, na.rm=TRUE), mean(ps$fallen_ha.2, na.rm=TRUE), mean(ps$fallen_ha.3, na.rm=TRUE)))
colnames(pm) <- c("cycle 1", "cycle 2", "cycle 3")
half(cf("pool_composition_by_cycle.png"))
barplot(t(pm), beside = TRUE, legend.text = colnames(pm), ylab = "mean carbon (t C/ha)",
        main = "Carbon pools by cycle", col = c("#88bb88","#bb8866","#aaaaaa"),
        args.legend = list(cex = 0.6, bty = "n"))
dev.off()

cat(sprintf("charts written to output/charts/ (%d files)\n", length(list.files(cf(".")))))

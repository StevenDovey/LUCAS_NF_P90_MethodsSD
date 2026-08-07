# Validates the dead-wood and CWD pathway against the reference output, in two
# layers:
#   (1) FORMULA: feed the reference's own volume, decay_modifier, carbon_fraction
#       and density into the dead-pool carbon formulas here; should match exactly.
#   (2) TAPER: compare the independent Newton-Raphson truncated dead-stem volume
#       against the reference's dead_stem_vol; this part was built from the
#       scientific concept, so some divergence is expected here.

root <- "/home/user/SAS2R"
source(file.path(root, "R/coefficients.R")); source(file.path(root, "R/allometry.R"))
co <- carbon_coefficients()

M <- file.path(root, "sasout/SAS_measurements-with-carbonV10.csv")
keep <- c("status","habit","Phylum","DBH","height","predht","original_height",
          "DecayClass","decay_modifier","stem_vol","dead_stem_vol",
          "Whole_stem_density","carbon_fraction","carbon_spars","carbon_fallen",
          "h_10","stem_vol_to_10cm_dia")
d <- read.csv(M, stringsAsFactors=FALSE, check.names=FALSE)[keep]
for (c in setdiff(keep, c("status","habit","Phylum")))
  d[[c]] <- suppressWarnings(as.numeric(d[[c]]))

pct <- function(my, ref, lab){ ok<-is.finite(my)&is.finite(ref)&ref!=0
  if(!any(ok)){cat(sprintf("  %-28s (no rows)\n",lab));return(invisible())}
  e<-(my[ok]-ref[ok])/ref[ok]*100
  cat(sprintf("  %-28s n=%6d  median %+.3f%%  |max| %7.2f%%  within0.1%%: %5.1f%%\n",
      lab, sum(ok), median(e), max(abs(e)), 100*mean(abs(e)<0.1))) }

woody <- c("canopy tree","subcanopy tree","shrub","unknown")
ferns <- c("tree fern","palm","cabbage tree")

cat("== (1) FORMULA check (using the reference's own volumes/modifiers) ==\n")
# Standing dead woody spar: density * decay_modifier * dead_stem_vol * carbon_fraction
sp <- d$status=="X" & d$habit %in% woody
my_spar <- d$Whole_stem_density * d$decay_modifier * d$dead_stem_vol * d$carbon_fraction
pct(my_spar[sp], d$carbon_spars[sp], "spar carbon (woody)")

# Standing dead fern: fern integrated carbon * decay_modifier
fr <- d$status=="X" & d$habit %in% ferns
hfern <- ifelse(!is.na(d$height), d$height, d$predht)
my_fern <- fern_aboveground_carbon(d$DBH, hfern, co) * d$decay_modifier
pct(my_fern[fr], d$carbon_spars[fr], "spar carbon (fern)")

# Fallen pieces: density * decay_modifier * stem_vol * carbon_fraction
fl <- d$status %in% c("F","P","S")
my_fallen <- d$Whole_stem_density * d$decay_modifier * d$stem_vol * d$carbon_fraction
pct(my_fallen[fl], d$carbon_fallen[fl], "fallen carbon (F/P/S)")

cat("\n== (2) TAPER check (my Newton-Raphson vs reference) ==\n")
keep2 <- c("original_stem_vol")
d2 <- read.csv(M, stringsAsFactors=FALSE, check.names=FALSE)[keep2]
d$original_stem_vol <- suppressWarnings(as.numeric(d2$original_stem_vol))
tp <- d$status=="X" & d$habit %in% woody & is.finite(d$original_stem_vol) &
      is.finite(d$original_height) & d$original_height>1.35 & is.finite(d$h_10)
my_h10 <- taper_height_at_10cm(d$original_stem_vol[tp], d$original_height[tp], co)
pct(my_h10, d$h_10[tp], "h_10 (relative height)")

#08.06.26 #1
# Per-stem biomass and carbon functions. All are vectorised and side-effect free, so they unit test against hand-computed values (see tests/). Inputs in cm, m and kg. `co` is the list from carbon_coefficients().

# Live stem volume (m3) for woody trees/shrubs.
live_stem_volume <- function(dbh_cm, height_m, co) {
  co$stem_volume_coef * (dbh_cm^2 * height_m)^co$stem_volume_exp
}

# Live stem volume (m3) for tree ferns / palms / cabbage trees.
fern_stem_volume <- function(dbh_cm, height_m, co) {
  co$fern_volume_coef * (dbh_cm^2 * height_m)^co$fern_volume_exp
}

branch_carbon  <- function(dbh_cm, co) co$branch_carbon_coef  * dbh_cm^co$branch_carbon_exp
foliage_carbon <- function(dbh_cm, co) co$foliage_carbon_coef * dbh_cm^co$foliage_carbon_exp

# Integrated above-ground carbon for tree ferns / palms / cabbage trees (kg).
fern_aboveground_carbon <- function(dbh_cm, height_m, co) {
  co$fern_carbon_coef * (dbh_cm^2 * height_m)^co$fern_carbon_exp
}

# Carbon fraction by live/dead and conifer/broadleaf.
carbon_fraction <- function(status, phylum, co) {
  is_dead <- status %in% c("X", "F", "P", "S")
  ifelse(is_dead, co$carbon_fraction_deadwood,
         ifelse(phylum == "C", co$carbon_fraction_conifer,
                co$carbon_fraction_broadleaf))
}

# Decay-class density modifier (0..4); unknown class -> NA.
decay_modifier <- function(decay_class, co) {
  unname(co$decay_class_density_modifier[as.character(decay_class)])
}

# Below-ground (root) carbon = group root:shoot ratio * AGB. Broadleaf trees/palms get the angiosperm ratio; conifers fall to default.
belowground_carbon <- function(carbon_agb, phylum, habit, co) {
  ratio <- ifelse(
    phylum %in% c("B", "P") & habit %in% c("canopy tree", "subcanopy tree", "palm"),
    co$bgb_ratio_broadleaf_palm,
    ifelse(habit == "tree fern",    co$bgb_ratio_tree_fern,
    ifelse(habit == "cabbage tree", co$bgb_ratio_cabbage_tree,
           co$bgb_ratio_default)))
  ratio * carbon_agb
}

# Above-ground carbon for a live woody stem (kg).
live_woody_agb <- function(dbh_cm, height_m, wood_density, phylum, co) {
  vol  <- live_stem_volume(dbh_cm, height_m, co)
  cf   <- ifelse(phylum == "C", co$carbon_fraction_conifer, co$carbon_fraction_broadleaf)
  stem <- wood_density * vol * cf
  stem + branch_carbon(dbh_cm, co) + foliage_carbon(dbh_cm, co)
}

# Kimberley compatible volume and sectional taper. Fraction of stem volume above relative height h (Beets 2012 form).
taper_volume_fraction_above <- function(h, co) {
  k <- co$taper_coef
  k["a"] * h^2 + k["b"] * h^3 + k["c"] * h^4 + k["d"] * h^5 + k["f"] * h^81
}

# Relative height at which stem diameter reaches 10 cm, solved by Newton-Raphson from total stem volume and height. The tap_* terms scale the sectional polynomial by volume over height; the solve target is diameter^2 = 0.10^2. Ported from the taper in functions.R.
taper_height_at_10cm <- function(stem_vol, height_m, co) {
  k <- co$taper_coef; pq <- co$taper_pi_quarter
  tap_a <- k["a"] * stem_vol * 2  / pq / height_m
  tap_b <- k["b"] * stem_vol * 3  / pq / height_m
  tap_c <- k["c"] * stem_vol * 4  / pq / height_m
  tap_d <- k["d"] * stem_vol * 5  / pq / height_m
  tap_f <- k["f"] * stem_vol * 81 / pq / height_m
  h <- rep(0.5, length(stem_vol))
  for (i in seq_len(co$taper_newton_iterations)) {
    f  <- tap_a*h + tap_b*h^2 + tap_c*h^3 + tap_d*h^4 + tap_f*h^80 - 0.01
    fp <- tap_a + 2*tap_b*h + 3*tap_c*h^2 + 4*tap_d*h^3 + 80*tap_f*h^79
    h  <- h - f / fp
  }
  unname(h)
}

# Stem volume truncated at the 10 cm top diameter (m3), woody live or dead.
stem_volume_to_10cm <- function(dbh_cm, height_m, co) {
  total <- live_stem_volume(dbh_cm, height_m, co)
  h10   <- taper_height_at_10cm(total, height_m, co)
  total * (1 - taper_volume_fraction_above(h10, co))
}

# Coarse woody debris piece volume (m3): conic frustum if both end diameters are present, otherwise a cylinder from the single end. Diameters in cm.
cwd_piece_volume <- function(led_cm, sed_cm, length_m) {
  led <- led_cm / 200; sed <- sed_cm / 200
  both <- !is.na(led_cm) & !is.na(sed_cm)
  one  <- ifelse(is.na(led_cm), sed, led)
  ifelse(both, (pi * length_m / 3) * (led^2 + sed^2 + led * sed),
               pi * length_m * one^2)
}

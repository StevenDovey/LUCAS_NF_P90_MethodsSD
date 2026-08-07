#08.06.26 01:30 NZST
# Projects coarse woody debris carbon forward through decay (Garrett et al. 2018).
# Each piece decays at a species rate adjusted for temperature and diameter.
# Existing debris decays over the full interval; debris from stems that died
# mid-interval decays for half the interval. Input carbon is already
# sampling-corrected.

# Per-piece decay constant (per year), adjusted for temperature and diameter.
# diam_cm in cm, mat_c in degrees C.
adjusted_decay_lambda <- function(species_code, diam_cm, mat_c, co) {
  base <- co$decay_lambda_by_species[species_code]
  base[is.na(base)] <- co$decay_lambda_base
  unname(base) * exp(co$decay_mat_coef * (mat_c - co$decay_mat_ref_c)) /
    (1 + co$decay_diameter_coef * (diam_cm - co$decay_diameter_ref_cm))
}

# Fraction of debris still remaining after delta_time_yr years of decay.
decay_fraction <- function(delta_time_yr, lambda) exp(-delta_time_yr * lambda)

# Carbon added to debris when a woody stem dies: stem volume up to the 10 cm
# top, times wood density and deadwood carbon fraction.
mortality_debris_carbon_woody <- function(dbh_cm, height_m, wood_density, co) {
  wood_density * stem_volume_to_10cm(dbh_cm, height_m, co) * co$carbon_fraction_deadwood
}

# Carbon added to debris when a tree fern, palm or cabbage tree dies: its
# above-ground carbon (no separate fraction applies).
mortality_debris_carbon_fern <- function(dbh_cm, height_m, co) {
  fern_aboveground_carbon(dbh_cm, height_m, co)
}

# Carry debris carbon from one cycle to the next: old debris decayed over the
# full interval, plus new debris from deaths decayed over half the interval.
project_cwd <- function(carbon_existing, carbon_mortality_gain, delta_time_yr, lambda) {
  carbon_existing * decay_fraction(delta_time_yr, lambda) +
    carbon_mortality_gain * decay_fraction(delta_time_yr / 2, lambda)
}

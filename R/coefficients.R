#08.06.26 #1
# Empirical constants for the carbon pipeline, grouped by what they describe. Values come from published sources, cited inline. A VERIFY tag marks a value not yet checked against the primary reference (see docs/clean_rebuild_spec.md).

carbon_coefficients <- function() {
  list(

    # --- Plot geometry (LUCAS PSP design) -----------------------------------
    plot_area_circular_ha      = 0.1257,     # 20 m radius full plot
    plot_area_inner_default_ha = 0.033975,   # inner subplot (small stems)

    # DBH size cut-offs (cm)
    dbh_min_live_inner_cm = 2.5,   # smallest live stem counted (inner subplot)
    dbh_min_dead_inner_cm = 10,    # smallest dead/CWD counted (inner subplot)
    dbh_circular_cm       = 60,    # >= this measured on full circular plot

    # --- Live stem volume, m3  (Beets et al. 2012, Forests 3(3):818) --------
    # VERIFY: coefficient + exponent against Beets 2012
    stem_volume_coef = 0.0000483,  # VERIFY
    stem_volume_exp  = 0.978,      # VERIFY

    # --- Branch & foliage carbon from DBH (Beets et al. 2012) ---------------
    branch_carbon_coef  = 0.0175,  # VERIFY
    branch_carbon_exp   = 2.20,    # VERIFY
    foliage_carbon_coef = 0.0171,  # VERIFY
    foliage_carbon_exp  = 1.75,    # VERIFY

    # --- Tree fern / palm / cabbage tree (Beets et al. 2012) ----------------
    fern_volume_coef = 0.00001343, # VERIFY  stem volume
    fern_volume_exp  = 1.22,       # VERIFY
    fern_carbon_coef = 0.00270,    # VERIFY  integrated above-ground carbon
    fern_carbon_exp  = 1.19,       # VERIFY

    # --- Carbon fractions (Beets 2012 / IPCC 2006) --------------------------
    carbon_fraction_conifer   = 0.51,  # corroborated
    carbon_fraction_broadleaf = 0.48,  # corroborated
    carbon_fraction_deadwood  = 0.50,

    # --- Wood density fallback ----------------------------------------------
    wood_density_unknown_kg_m3 = 477,  # VERIFY

    # --- Standing-dead / CWD decay-class density modifiers (Coomes 2002) ----
    # index by decay class 0..4 ; class 4 = fully decayed (assumed 0)
    decay_class_density_modifier = c(`0` = 1.00, `1` = 0.82, `2` = 0.66,
                                     `3` = 0.47, `4` = 0.00),  # 0-3 corroborated

    # --- Below-ground (root:shoot) fractions of AGB (Easdale et al. 2019) ----
    # NOTE: applied to BROADLEAF trees/palms (phylum B/P), NOT conifers. Conifers (phylum C) take the default 0.245 - matches Easdale gymnosperm.
    bgb_ratio_broadleaf_palm = 0.234,  # corroborated (angiosperm)
    bgb_ratio_tree_fern      = 0.194,  # corroborated
    bgb_ratio_cabbage_tree   = 0.437,  # Easdale 2019 monocots (confirmed)
    bgb_ratio_default        = 0.245,  # corroborated (gymnosperm/conifer)

    # --- Kimberley compatible volume / sectional taper (Beets 2012) ---------
    # fraction of stem volume above relative height h ; used to truncate dead-spar volume at the 10 cm top-diameter point via Newton-Raphson.
    taper_coef = c(a = 0.06501, b = 2.92127, c = -3.37103,
                   d = 1.35551, f = 0.02924),  # VERIFY
    taper_pi_quarter           = 0.7854,
    taper_top_diameter_cm      = 10,
    taper_newton_iterations    = 5L,

    # --- CWD exponential decay (Garrett et al. 2018) ------------------------
    decay_lambda_base = 0.0216,  # VERIFY (secondary mean ~0.0270)
    decay_lambda_by_species = c(
      BEITAW = 0.0233, DACCUP = 0.0350, DACDAC = 0.0506, PRUTAX = 0.0171,
      PRUFER = 0.0259, WEIRAC = 0.0220, NOTMEN = 0.0246, LOPMEN = 0.0246,
      NOTFUS = 0.0204, FUSFUS = 0.0204, NOTSOL = 0.0162, NOTSVS = 0.0162
    ),  # VERIFY
    decay_mat_coef        = 0.093,    # VERIFY
    decay_diameter_coef   = 0.00908,  # VERIFY
    decay_diameter_ref_cm = 60,
    decay_mat_ref_c       = 10,
    mat_default_c         = 10.49,

    # --- Dead-wood scaling applied once at plot summary ----------------------
    cwd_belowground_uplift  = 1.19,   # below-ground dead wood add-on
    cwd_sampling_multiplier = 1.767   # sampling-bias correction
  )
}

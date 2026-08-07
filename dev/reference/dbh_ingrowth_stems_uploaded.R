# DBH ingrowth and back-prediction (sourced from SAS2R_final_*.R).
# Wide `stems` has dbh1/dbh2/dbh3. This file fills missing DBHs using within-plot
# regression (see `dbh_fit_pred_plot` in functions.R): e.g. predict cycle-2 DBH from
# cycle-3 when both are woody live stems.
#
# constants.R: grid_row / *_grid = structural cycle slot (not a raw visit row).
# input_visit_keys (built in SAS2R before sourcing this file) = one row per
# measured stem/CWD visit (plot, ItemID, database_table, cycle). Revert any
# cycle column to its pre-ingrowth state when that stem has no visit in that
# cycle, the same outcome as a blank slot for model fills, but this also clears
# slots that were non-blank only from wide-table or back-fill steps (not an
# actual measurement).

if (!exists("exclude_dbh_neg_growth_from_training", inherits = TRUE)) {
  exclude_dbh_neg_growth_from_training <- FALSE
}
if (!exists("DbhNegGrowthDropCm", inherits = TRUE)) {
  DbhNegGrowthDropCm <- 0.5
}

if (!exists("input_visit_keys", inherits = TRUE) || nrow(input_visit_keys) == 0L) {
  stems <- stems %>%
    dplyr::mutate(
      stem_visit1 = TRUE,
      stem_visit2 = TRUE,
      stem_visit3 = TRUE
    )
} else {
  STEM_VISIT_K <- input_visit_keys %>%
    dplyr::mutate(cycle = as.integer(.data$cycle)) %>%
    dplyr::group_by(.data$plot, .data$ItemID, .data$database_table) %>%
    dplyr::summarise(
      stem_visit1 = any(.data$cycle == 1L, na.rm = TRUE),
      stem_visit2 = any(.data$cycle == 2L, na.rm = TRUE),
      stem_visit3 = any(.data$cycle == 3L, na.rm = TRUE),
      .groups = "drop"
    )
  stems <- stems %>%
    dplyr::left_join(STEM_VISIT_K, by = c("plot", "ItemID", "database_table")) %>%
    dplyr::mutate(
      stem_visit1 = dplyr::coalesce(.data$stem_visit1, FALSE),
      stem_visit2 = dplyr::coalesce(.data$stem_visit2, FALSE),
      stem_visit3 = dplyr::coalesce(.data$stem_visit3, FALSE)
    )
}

stems <- stems %>%
  mutate(
    .src1_pre_dbh = source1,
    .src2_pre_dbh = source2,
    .src3_pre_dbh = source3,
    .dbh1_pre_dbh = dbh1,
    .dbh2_pre_dbh = dbh2,
    .dbh3_pre_dbh = dbh3,
    .species_code1_pre_dbh = species_code1,
    .species_code2_pre_dbh = species_code2,
    .species_code3_pre_dbh = species_code3,
    .habit1_pre_dbh = habit1,
    .habit2_pre_dbh = habit2,
    .habit3_pre_dbh = habit3,
    .subplot1_pre_dbh = subplot1,
    .subplot2_pre_dbh = subplot2,
    .subplot3_pre_dbh = subplot3,
    .status1_pre_dbh = status1,
    .status2_pre_dbh = status2,
    .status3_pre_dbh = status3,
    .height1_pre_dbh = height1,
    .decayclass1_pre_dbh = decayclass1
  )

# A→A consecutive DBH shrink (training exclusion only; predictions still use all rows).
stems <- stems %>%
  mutate(
    dbh_fit_excl_12 = exclude_dbh_neg_growth_from_training &
      status1 == "A" & status2 == "A" &
      !is.na(dbh1) & !is.na(dbh2) & dbh2 < dbh1 - DbhNegGrowthDropCm,
    dbh_fit_excl_23 = exclude_dbh_neg_growth_from_training &
      status2 == "A" & status3 == "A" &
      !is.na(dbh2) & !is.na(dbh3) & dbh3 < dbh2 - DbhNegGrowthDropCm,
    dbh_fit_excl_13 = exclude_dbh_neg_growth_from_training &
      status1 == "A" & status3 == "A" &
      !is.na(dbh1) & !is.na(dbh3) & dbh3 < dbh1 - DbhNegGrowthDropCm,
    dbh_neg_growth_training_exclude = dbh_fit_excl_12 | dbh_fit_excl_23 | dbh_fit_excl_13
  )

#########################################################################
############# Estimating DBH2 from DBH 3 (first: later cycle informs mid) #
#########################################################################
# Order of blocks matters: mid-cycle gaps filled before end-cycle-only paths.

stems_ingrowth2 <- stems %>%
  filter(
    habit3 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status3 == "A"
  )

preddbh2 <- stems_ingrowth2 %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh2", "dbh3", "species_code3", "pred_dbh2", "dbh_fit_excl_23")
    dplyr::transmute(
      fit$df,
      model_type2 = fit$model_type,
      .dbh_n = fit$diag$n_complete,
      .dbh_n_sp = fit$diag$n_species,
      .dbh_eq = fit$diag$equation,
      .dbh_r2 = fit$diag$r2,
      ItemID,
      database_table,
      pred_dbh2
    )
  }) %>%
  dplyr::ungroup()

dbh_diag_c2_c3 <- preddbh2 %>%
  dplyr::distinct(plot, model_type2, .dbh_n, .dbh_n_sp, .dbh_eq, .dbh_r2) %>%
  dplyr::rename(
    model_type = model_type2,
    n = .dbh_n,
    n_species = .dbh_n_sp,
    equation = .dbh_eq,
    r2 = .dbh_r2
  ) %>%
  dplyr::mutate(
    ingrowth_stage = "cycle2←cycle3",
    resp = "dbh2",
    pred = "dbh3",
    species_term = "species_code3"
  )

preddbh2 <- preddbh2 %>%
  dplyr::select(-.dbh_n, -.dbh_n_sp, -.dbh_eq, -.dbh_r2)

preddbh2 <- preddbh2 %>%
  arrange(plot, ItemID, database_table)

stems <- stems %>%
  arrange(plot, ItemID, database_table)

stems <- stems %>%
  left_join(preddbh2, by = c("plot", "ItemID", "database_table"))

stems <- stems %>%
  mutate(
    dbh2_was_missing = is.na(dbh2) & !is.na(pred_dbh2)
  ) %>%
  mutate(
    dbh2 = if_else(dbh2_was_missing, dbh_use_non_negative(pred_dbh2), dbh2),
    species_code2 = if_else(dbh2_was_missing & (is.na(species_code2) | species_code2 == ""),
                            species_code3, species_code2),
    habit2 = if_else(dbh2_was_missing & (is.na(habit2) | habit2 == ""),
                     habit3, habit2),
    subplot2 = if_else(dbh2_was_missing & (is.na(subplot2) | subplot2 == ""),
                       subplot3, subplot2),
    status2 = if_else(dbh2_was_missing & (is.na(status2) | status2 == ""),
                      status3, status2),
    database_table = if_else(dbh2_was_missing & (is.na(database_table) | database_table == ""),
                             "stems", database_table),
    source2 = if_else(dbh2_was_missing, "m_2fr3", source2)
  ) %>%
  select(-dbh2_was_missing)

stems <- stems %>%
  mutate(
    is_estimated = habit3 %in% c("tree fern", "palm", "cabbage tree") &
      status3 == "A" &
      is.na(dbh2)
  ) %>%
  mutate(
    dbh2 = if_else(is_estimated, dbh_use_non_negative(dbh3), dbh2),
    species_code2 = if_else(is_estimated & (is.na(species_code2) | species_code2 == ""),
                            species_code3, species_code2),
    habit2 = if_else(is_estimated & (is.na(habit2) | habit2 == ""),
                     habit3, habit2),
    subplot2 = if_else(is_estimated & (is.na(subplot2) | subplot2 == ""),
                       subplot3, subplot2),
    status2 = if_else(is_estimated & (is.na(status2) | status2 == ""),
                      status3, status2),
    database_table = if_else(is_estimated & (is.na(database_table) | database_table == ""),
                             "stems", database_table),
    source2 = if_else(is_estimated, "f_3to2", source2)
  ) %>%
  select(-is_estimated)

#########################################################################
############# Forward: DBH2 from DBH1 (SAS %estimate_DBH2 parity) #########
#########################################################################

stems_ingrowth2_fwd <- stems %>%
  filter(
    habit1 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status2 == "A"
  )

preddbh2_c1 <- stems_ingrowth2_fwd %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh2", "dbh1", "species_code1", "pred_dbh2_c1", "dbh_fit_excl_12")
    dplyr::transmute(
      fit$df,
      model_type2c1 = fit$model_type,
      .dbh_n = fit$diag$n_complete,
      .dbh_n_sp = fit$diag$n_species,
      .dbh_eq = fit$diag$equation,
      .dbh_r2 = fit$diag$r2,
      ItemID,
      database_table,
      pred_dbh2_c1
    )
  }) %>%
  dplyr::ungroup()

dbh_diag_c2_c1 <- preddbh2_c1 %>%
  dplyr::distinct(plot, model_type2c1, .dbh_n, .dbh_n_sp, .dbh_eq, .dbh_r2) %>%
  dplyr::rename(
    model_type = model_type2c1,
    n = .dbh_n,
    n_species = .dbh_n_sp,
    equation = .dbh_eq,
    r2 = .dbh_r2
  ) %>%
  dplyr::mutate(
    ingrowth_stage = "cycle2←cycle1",
    resp = "dbh2",
    pred = "dbh1",
    species_term = "species_code1"
  )

preddbh2_c1 <- preddbh2_c1 %>%
  dplyr::select(plot, ItemID, database_table, pred_dbh2_c1, model_type2c1)

preddbh2_c1 <- preddbh2_c1 %>%
  arrange(plot, ItemID, database_table)

stems <- stems %>%
  arrange(plot, ItemID, database_table) %>%
  left_join(preddbh2_c1, by = c("plot", "ItemID", "database_table"))

stems <- stems %>%
  mutate(
    dbh2_was_fwd = is.na(dbh2) & is.finite(pred_dbh2_c1)
  ) %>%
  mutate(
    dbh2 = if_else(dbh2_was_fwd, dbh_use_non_negative(pred_dbh2_c1), dbh2),
    species_code2 = if_else(
      dbh2_was_fwd & (is.na(species_code2) | species_code2 == ""),
      species_code1,
      species_code2
    ),
    habit2 = if_else(
      dbh2_was_fwd & (is.na(habit2) | habit2 == ""),
      habit1,
      habit2
    ),
    subplot2 = if_else(
      dbh2_was_fwd & (is.na(subplot2) | subplot2 == ""),
      subplot1,
      subplot2
    ),
    status2 = if_else(
      dbh2_was_fwd & (is.na(status2) | status2 == ""),
      status1,
      status2
    ),
    database_table = if_else(
      dbh2_was_fwd & (is.na(database_table) | database_table == ""),
      "stems",
      database_table
    ),
    source2 = if_else(dbh2_was_fwd, "m_2fr1", source2)
  ) %>%
  select(-dbh2_was_fwd)

stems <- stems %>%
  mutate(
    is_estimated = habit2 %in% c("tree fern", "palm", "cabbage tree") &
      status2 == "A" &
      is.na(dbh2)
  ) %>%
  mutate(
    dbh2 = if_else(is_estimated, dbh_use_non_negative(dbh1), dbh2),
    species_code2 = if_else(
      is_estimated & (is.na(species_code2) | species_code2 == ""),
      species_code1,
      species_code2
    ),
    habit2 = if_else(
      is_estimated & (is.na(habit2) | habit2 == ""),
      habit1,
      habit2
    ),
    subplot2 = if_else(
      is_estimated & (is.na(subplot2) | subplot2 == ""),
      subplot1,
      subplot2
    ),
    status2 = if_else(
      is_estimated & (is.na(status2) | status2 == ""),
      status1,
      status2
    ),
    database_table = if_else(
      is_estimated & (is.na(database_table) | database_table == ""),
      "stems",
      database_table
    ),
    source2 = if_else(is_estimated, "f_1to2", source2)
  ) %>%
  select(-is_estimated)

#########################################################################
############# Forward: DBH3 from DBH2 (three-cycle extension) #############
#########################################################################

stems_ingrowth3_fwd <- stems %>%
  filter(
    habit2 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status3 == "A"
  )

preddbh3_c2 <- stems_ingrowth3_fwd %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh3", "dbh2", "species_code2", "pred_dbh3_c2", "dbh_fit_excl_23")
    dplyr::transmute(
      fit$df,
      model_type3c2 = fit$model_type,
      .dbh_n = fit$diag$n_complete,
      .dbh_n_sp = fit$diag$n_species,
      .dbh_eq = fit$diag$equation,
      .dbh_r2 = fit$diag$r2,
      ItemID,
      database_table,
      pred_dbh3_c2
    )
  }) %>%
  dplyr::ungroup()

dbh_diag_c3_c2 <- preddbh3_c2 %>%
  dplyr::distinct(plot, model_type3c2, .dbh_n, .dbh_n_sp, .dbh_eq, .dbh_r2) %>%
  dplyr::rename(
    model_type = model_type3c2,
    n = .dbh_n,
    n_species = .dbh_n_sp,
    equation = .dbh_eq,
    r2 = .dbh_r2
  ) %>%
  dplyr::mutate(
    ingrowth_stage = "cycle3←cycle2",
    resp = "dbh3",
    pred = "dbh2",
    species_term = "species_code2"
  )

preddbh3_c2 <- preddbh3_c2 %>%
  dplyr::select(plot, ItemID, database_table, pred_dbh3_c2, model_type3c2)

preddbh3_c2 <- preddbh3_c2 %>%
  arrange(plot, ItemID, database_table)

stems <- stems %>%
  arrange(plot, ItemID, database_table) %>%
  left_join(preddbh3_c2, by = c("plot", "ItemID", "database_table"))

stems <- stems %>%
  mutate(
    dbh3_was_fwd = is.na(dbh3) & is.finite(pred_dbh3_c2)
  ) %>%
  mutate(
    dbh3 = if_else(dbh3_was_fwd, dbh_use_non_negative(pred_dbh3_c2), dbh3),
    species_code3 = if_else(
      dbh3_was_fwd & (is.na(species_code3) | species_code3 == ""),
      species_code2,
      species_code3
    ),
    habit3 = if_else(
      dbh3_was_fwd & (is.na(habit3) | habit3 == ""),
      habit2,
      habit3
    ),
    subplot3 = if_else(
      dbh3_was_fwd & (is.na(subplot3) | subplot3 == ""),
      subplot2,
      subplot3
    ),
    status3 = if_else(
      dbh3_was_fwd & (is.na(status3) | status3 == ""),
      status2,
      status3
    ),
    database_table = if_else(
      dbh3_was_fwd & (is.na(database_table) | database_table == ""),
      "stems",
      database_table
    ),
    source3 = if_else(dbh3_was_fwd, "m_3fr2", source3)
  ) %>%
  select(-dbh3_was_fwd)

stems <- stems %>%
  mutate(
    is_estimated = habit3 %in% c("tree fern", "palm", "cabbage tree") &
      status3 == "A" &
      is.na(dbh3)
  ) %>%
  mutate(
    .dbh_src = dplyr::coalesce(dbh2, dbh1),
    dbh3 = if_else(is_estimated, dbh_use_non_negative(.dbh_src), dbh3),
    species_code3 = if_else(
      is_estimated & (is.na(species_code3) | species_code3 == ""),
      dplyr::coalesce(species_code2, species_code1),
      species_code3
    ),
    habit3 = if_else(
      is_estimated & (is.na(habit3) | habit3 == ""),
      dplyr::coalesce(habit2, habit1),
      habit3
    ),
    subplot3 = if_else(
      is_estimated & (is.na(subplot3) | subplot3 == ""),
      dplyr::coalesce(subplot2, subplot1),
      subplot3
    ),
    status3 = if_else(
      is_estimated & (is.na(status3) | status3 == ""),
      dplyr::coalesce(status2, status1),
      status3
    ),
    database_table = if_else(
      is_estimated & (is.na(database_table) | database_table == ""),
      "stems",
      database_table
    ),
    source3 = dplyr::if_else(
      is_estimated,
      dplyr::if_else(!is.na(dbh2), "f_2to3", "f_1to3"),
      source3
    )
  ) %>%
  select(-is_estimated, -.dbh_src)

#########################################################################
############# ingrowth DBH for Cycle 1 from cycle 2 (join preds only) ####
#########################################################################

stems_ingrowth1 <- stems %>%
  filter(
    habit2 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status2 == "A"
  )

preddbh1 <- stems_ingrowth1 %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh1", "dbh2", "species_code2", "pred_dbh1", "dbh_fit_excl_12")
    dplyr::transmute(
      fit$df,
      model_type1 = fit$model_type,
      .dbh_n = fit$diag$n_complete,
      .dbh_n_sp = fit$diag$n_species,
      .dbh_eq = fit$diag$equation,
      .dbh_r2 = fit$diag$r2,
      ItemID,
      database_table,
      pred_dbh1
    )
  }) %>%
  dplyr::ungroup()

dbh_diag_c1_c2 <- preddbh1 %>%
  dplyr::distinct(plot, model_type1, .dbh_n, .dbh_n_sp, .dbh_eq, .dbh_r2) %>%
  dplyr::rename(
    model_type = model_type1,
    n = .dbh_n,
    n_species = .dbh_n_sp,
    equation = .dbh_eq,
    r2 = .dbh_r2
  ) %>%
  dplyr::mutate(
    ingrowth_stage = "cycle1←cycle2",
    resp = "dbh1",
    pred = "dbh2",
    species_term = "species_code2"
  )

preddbh1 <- preddbh1 %>%
  dplyr::select(-.dbh_n, -.dbh_n_sp, -.dbh_eq, -.dbh_r2)

preddbh1 <- preddbh1 %>%
  arrange(plot, ItemID, database_table)

stems <- stems %>%
  arrange(plot, ItemID, database_table)

stems <- stems %>%
  left_join(preddbh1, by = c("plot", "ItemID", "database_table"))

stems <- stems %>%
  mutate(
    is_estimated = habit2 %in% c("tree fern", "palm", "cabbage tree") &
      status2 == "A" &
      is.na(dbh1)
  ) %>%
  mutate(
    dbh1 = if_else(is_estimated, dbh_use_non_negative(dbh2), dbh1),
    species_code1 = if_else(is_estimated & (is.na(species_code1) | species_code1 == ""),
                            species_code2, species_code1),
    habit1 = if_else(is_estimated & (is.na(habit1) | habit1 == ""),
                     habit2, habit1),
    subplot1 = if_else(is_estimated & (is.na(subplot1) | subplot1 == ""),
                       subplot2, subplot1),
    status1 = if_else(is_estimated & (is.na(status1) | status1 == ""),
                      status2, status1),
    database_table = if_else(is_estimated & (is.na(database_table) | database_table == ""),
                             "stems", database_table),
    source1 = if_else(is_estimated, "f_2to1", source1)
  ) %>%
  select(-is_estimated)

###########################################################################################
#################### Populate missing Spars (Decay <3) that occurred in Cycle 2 ###########
###########################################################################################

stems <- stems %>%
  mutate(flag = !is.na(status2) & status2 == "X" & !is.na(decayclass2) & decayclass2 < 3 & is.na(dbh1)) %>%
  mutate(
    dbh1 = if_else(flag, dbh_use_non_negative(dbh2), dbh1),
    height1 = if_else(flag, height2, height1),
    species_code1 = if_else(flag, species_code2, species_code1),
    habit1 = if_else(flag, habit2, habit1),
    subplot1 = if_else(flag, subplot2, subplot1),
    status1 = if_else(flag, status2, status1),
    decayclass1 = if_else(flag, decayclass2, decayclass1),
    database_table = if_else(flag, "stems", database_table),
    source1 = if_else(flag, "sp_21", source1)
  ) %>%
  select(-flag)

#########################################################################
############# ingrowth estimation DBH for Cycle 1 from cycle3 #######################
#########################################################################

stems_ingrowth1_3 <- stems %>%
  filter(
    habit3 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status3 == "A"
  )

preddbh1_3 <- stems_ingrowth1_3 %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh1", "dbh3", "species_code3", "pred_dbh1a", "dbh_neg_growth_training_exclude")
    dplyr::transmute(
      fit$df,
      model_type1a = fit$model_type,
      .dbh_n = fit$diag$n_complete,
      .dbh_n_sp = fit$diag$n_species,
      .dbh_eq = fit$diag$equation,
      .dbh_r2 = fit$diag$r2,
      ItemID,
      database_table,
      pred_dbh1a
    )
  }) %>%
  dplyr::ungroup()

dbh_diag_c1_c3 <- preddbh1_3 %>%
  dplyr::distinct(plot, model_type1a, .dbh_n, .dbh_n_sp, .dbh_eq, .dbh_r2) %>%
  dplyr::rename(
    model_type = model_type1a,
    n = .dbh_n,
    n_species = .dbh_n_sp,
    equation = .dbh_eq,
    r2 = .dbh_r2
  ) %>%
  dplyr::mutate(
    ingrowth_stage = "cycle1←cycle3",
    resp = "dbh1",
    pred = "dbh3",
    species_term = "species_code3"
  )

preddbh1_3 <- preddbh1_3 %>%
  dplyr::select(-.dbh_n, -.dbh_n_sp, -.dbh_eq, -.dbh_r2)

preddbh1_3 <- preddbh1_3 %>%
  arrange(plot, ItemID, database_table)

stems <- stems %>%
  arrange(plot, ItemID, database_table)

stems <- stems %>%
  left_join(preddbh1_3, by = c("plot", "ItemID", "database_table"))

stems <- stems %>%
  mutate(
    .p1 = pred_dbh1,
    .p3 = pred_dbh1a,
    .has1 = is.finite(.p1),
    .has3 = is.finite(.p3),
    .pred_combined = dplyr::case_when(
      .has1 & .has3 & .p1 > 0 & .p3 <= 0 ~ .p1,
      .has1 & .has3 & .p3 > 0 & .p1 <= 0 ~ .p3,
      .has1 & .has3 & .p1 > 0 & .p3 > 0 ~ pmax(.p1, .p3),
      .has1 & .has3 ~ 0,
      .has1 & .p1 > 0 ~ .p1,
      .has1 ~ 0,
      .has3 & .p3 > 0 ~ .p3,
      .has3 ~ 0,
      TRUE ~ NA_real_
    ),
    dbh1_was_missing = is.na(dbh1) & (.has1 | .has3),
    .used_c2 = dbh1_was_missing & .pred_combined > 0 & (
      (.has1 & .p1 > 0 & (!.has3 | .p3 <= 0)) |
        (.has1 & .has3 & .p1 > 0 & .p3 > 0 & .p1 >= .p3)
    ),
    .used_c3 = dbh1_was_missing & .pred_combined > 0 & !.used_c2,
    .used_dbh0 = dbh1_was_missing & !is.na(.pred_combined) & .pred_combined == 0
  ) %>%
  mutate(
    dbh1 = if_else(dbh1_was_missing, dbh_use_non_negative(.pred_combined), dbh1),
    species_code1 = dplyr::case_when(
      .used_c2 & (is.na(species_code1) | species_code1 == "") ~ species_code2,
      .used_c3 & (is.na(species_code1) | species_code1 == "") ~ species_code3,
      .used_dbh0 & (is.na(species_code1) | species_code1 == "") ~
        dplyr::coalesce(species_code3, species_code2),
      TRUE ~ species_code1
    ),
    habit1 = dplyr::case_when(
      .used_c2 & (is.na(habit1) | habit1 == "") ~ habit2,
      .used_c3 & (is.na(habit1) | habit1 == "") ~ habit3,
      .used_dbh0 & (is.na(habit1) | habit1 == "") ~ dplyr::coalesce(habit3, habit2),
      TRUE ~ habit1
    ),
    subplot1 = dplyr::case_when(
      .used_c2 & (is.na(subplot1) | subplot1 == "") ~ subplot2,
      .used_c3 & (is.na(subplot1) | subplot1 == "") ~ subplot3,
      .used_dbh0 & (is.na(subplot1) | subplot1 == "") ~ dplyr::coalesce(subplot3, subplot2),
      TRUE ~ subplot1
    ),
    status1 = dplyr::case_when(
      .used_c2 & (is.na(status1) | status1 == "") ~ status2,
      .used_c3 & (is.na(status1) | status1 == "") ~ status3,
      .used_dbh0 & (is.na(status1) | status1 == "") ~ dplyr::coalesce(status3, status2),
      TRUE ~ status1
    ),
    database_table = if_else(
      dbh1_was_missing & (is.na(database_table) | database_table == ""),
      "stems",
      database_table
    ),
    source1 = dplyr::case_when(
      .used_dbh0 ~ "m_1_zero",
      .used_c2 ~ "m_1fr2",
      .used_c3 ~ "m_1fr3",
      TRUE ~ source1
    )
  ) %>%
  select(
    -dbh1_was_missing, -.p1, -.p3, -.has1, -.has3, -.pred_combined,
    -.used_c2, -.used_c3, -.used_dbh0
  )

stems <- stems %>%
  mutate(
    is_estimated = habit3 %in% c("tree fern", "palm", "cabbage tree") &
      status3 == "A" &
      is.na(dbh1)
  ) %>%
  mutate(
    dbh1 = if_else(is_estimated, dbh_use_non_negative(dbh3), dbh1),
    species_code1 = if_else(is_estimated & (is.na(species_code1) | species_code1 == ""),
                            species_code3, species_code1),
    habit1 = if_else(is_estimated & (is.na(habit1) | habit1 == ""),
                     habit3, habit1),
    subplot1 = if_else(is_estimated & (is.na(subplot1) | subplot1 == ""),
                       subplot3, subplot1),
    status1 = if_else(is_estimated & (is.na(status1) | status1 == ""),
                      status3, status1),
    database_table = if_else(is_estimated & (is.na(database_table) | database_table == ""),
                             "stems", database_table),
    source1 = if_else(is_estimated, "f_3to1", source1)
  ) %>%
  select(-is_estimated)

#########################################################################
############# Forward fallback: DBH3 from DBH1 (after C1 ingrowth) #######
#########################################################################

stems_ingrowth3_c1 <- stems %>%
  filter(
    habit1 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status3 == "A"
  )

preddbh3_c1 <- stems_ingrowth3_c1 %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh3", "dbh1", "species_code1", "pred_dbh3_c1", "dbh_neg_growth_training_exclude")
    dplyr::transmute(
      fit$df,
      model_type3c1 = fit$model_type,
      .dbh_n = fit$diag$n_complete,
      .dbh_n_sp = fit$diag$n_species,
      .dbh_eq = fit$diag$equation,
      .dbh_r2 = fit$diag$r2,
      ItemID,
      database_table,
      pred_dbh3_c1
    )
  }) %>%
  dplyr::ungroup()

dbh_diag_c3_c1 <- preddbh3_c1 %>%
  dplyr::distinct(plot, model_type3c1, .dbh_n, .dbh_n_sp, .dbh_eq, .dbh_r2) %>%
  dplyr::rename(
    model_type = model_type3c1,
    n = .dbh_n,
    n_species = .dbh_n_sp,
    equation = .dbh_eq,
    r2 = .dbh_r2
  ) %>%
  dplyr::mutate(
    ingrowth_stage = "cycle3←cycle1",
    resp = "dbh3",
    pred = "dbh1",
    species_term = "species_code1"
  )

preddbh3_c1 <- preddbh3_c1 %>%
  dplyr::select(plot, ItemID, database_table, pred_dbh3_c1, model_type3c1)

preddbh3_c1 <- preddbh3_c1 %>%
  arrange(plot, ItemID, database_table)

stems <- stems %>%
  arrange(plot, ItemID, database_table) %>%
  left_join(preddbh3_c1, by = c("plot", "ItemID", "database_table"))

stems <- stems %>%
  mutate(
    dbh3_was_fwd1 = is.na(dbh3) & is.finite(pred_dbh3_c1)
  ) %>%
  mutate(
    dbh3 = if_else(dbh3_was_fwd1, dbh_use_non_negative(pred_dbh3_c1), dbh3),
    species_code3 = if_else(
      dbh3_was_fwd1 & (is.na(species_code3) | species_code3 == ""),
      species_code1,
      species_code3
    ),
    habit3 = if_else(
      dbh3_was_fwd1 & (is.na(habit3) | habit3 == ""),
      habit1,
      habit3
    ),
    subplot3 = if_else(
      dbh3_was_fwd1 & (is.na(subplot3) | subplot3 == ""),
      subplot1,
      subplot3
    ),
    status3 = if_else(
      dbh3_was_fwd1 & (is.na(status3) | status3 == ""),
      status1,
      status3
    ),
    database_table = if_else(
      dbh3_was_fwd1 & (is.na(database_table) | database_table == ""),
      "stems",
      database_table
    ),
    source3 = if_else(dbh3_was_fwd1, "m_3fr1", source3)
  ) %>%
  select(-dbh3_was_fwd1)

# Effective DBH for SAS-parity ingrowth tests.
# Thresholds are computed HERE, inside the same mutate as the stem_visit revert,
# so that non-visit stems get NA thresholds (matching SAS behaviour where dbh is
# simply missing for non-measurement events).
# SAS %estimate_ingrowth_DBH only predicts dbh1 from dbh2, not from dbh3.
# SAS %estimate_DBH2 only predicts dbh2 from dbh1 (forward), not from dbh3.
# NOTE: threshold uses stem_visit flags before they are dropped by select() below.
stems <- stems %>%
  mutate(
    # No raw visit for this stem in cycle k: strip ingrowth (revert to pre-file state).
    dbh1 = if_else(!stem_visit1, .dbh1_pre_dbh, dbh1),
    species_code1 = if_else(!stem_visit1, .species_code1_pre_dbh, species_code1),
    habit1 = if_else(!stem_visit1, .habit1_pre_dbh, habit1),
    subplot1 = if_else(!stem_visit1, .subplot1_pre_dbh, subplot1),
    status1 = if_else(!stem_visit1, .status1_pre_dbh, status1),
    source1 = if_else(!stem_visit1, .src1_pre_dbh, source1),
    height1 = if_else(!stem_visit1, .height1_pre_dbh, height1),
    decayclass1 = if_else(!stem_visit1, .decayclass1_pre_dbh, decayclass1),
    pred_dbh1 = if_else(!stem_visit1, NA_real_, pred_dbh1),
    pred_dbh1a = if_else(!stem_visit1, NA_real_, pred_dbh1a),
    model_type1 = if_else(!stem_visit1, NA_character_, model_type1),
    model_type1a = if_else(!stem_visit1, NA_character_, model_type1a),
    dbh2 = if_else(!stem_visit2, .dbh2_pre_dbh, dbh2),
    species_code2 = if_else(!stem_visit2, .species_code2_pre_dbh, species_code2),
    habit2 = if_else(!stem_visit2, .habit2_pre_dbh, habit2),
    subplot2 = if_else(!stem_visit2, .subplot2_pre_dbh, subplot2),
    status2 = if_else(!stem_visit2, .status2_pre_dbh, status2),
    source2 = if_else(!stem_visit2, .src2_pre_dbh, source2),
    pred_dbh2 = if_else(!stem_visit2, NA_real_, pred_dbh2),
    pred_dbh2_c1 = if_else(!stem_visit2, NA_real_, pred_dbh2_c1),
    model_type2 = if_else(!stem_visit2, NA_character_, model_type2),
    model_type2c1 = if_else(!stem_visit2, NA_character_, model_type2c1),
    dbh3 = if_else(!stem_visit3, .dbh3_pre_dbh, dbh3),
    species_code3 = if_else(!stem_visit3, .species_code3_pre_dbh, species_code3),
    habit3 = if_else(!stem_visit3, .habit3_pre_dbh, habit3),
    subplot3 = if_else(!stem_visit3, .subplot3_pre_dbh, subplot3),
    status3 = if_else(!stem_visit3, .status3_pre_dbh, status3),
    source3 = if_else(!stem_visit3, .src3_pre_dbh, source3),
    pred_dbh3_c2 = if_else(!stem_visit3, NA_real_, pred_dbh3_c2),
    pred_dbh3_c1 = if_else(!stem_visit3, NA_real_, pred_dbh3_c1),
    model_type3c2 = if_else(!stem_visit3, NA_character_, model_type3c2),
    model_type3c1 = if_else(!stem_visit3, NA_character_, model_type3c1),
    # SAS parity: when a measured DBH was present before ingrowth, keep original source label.
    source1 = if_else(stem_visit1 & !is.na(.dbh1_pre_dbh), .src1_pre_dbh, source1),
    source2 = if_else(stem_visit2 & !is.na(.dbh2_pre_dbh), .src2_pre_dbh, source2),
    source3 = if_else(stem_visit3 & !is.na(.dbh3_pre_dbh), .src3_pre_dbh, source3),
    # Ingrowth thresholds: built AFTER revert so non-visit stems get NA, not a
    # predicted value. NA threshold matches SAS missing-DBH behaviour (SAS treats
    # missing numeric as < any threshold, so NA here is handled in SAS2R_final
    # by treating is.na(threshold) as below-threshold, see the ingrowth gate there).
    # Use original dbh1 (not filled-in) + pred_dbh1 from dbh2 only (matching SAS
    # %estimate_ingrowth_DBH which predicts dbh1 from dbh2, not dbh3).
    dbh1_ingrowth_threshold = if_else(
      stem_visit1,
      dplyr::coalesce(.dbh1_pre_dbh, pred_dbh1),
      NA_real_
    ),
    # Use original dbh2 + pred_dbh2_c1 forward from dbh1 (matching SAS %estimate_DBH2).
    dbh2_ingrowth_threshold = if_else(
      stem_visit2,
      dplyr::coalesce(.dbh2_pre_dbh, pred_dbh2_c1),
      NA_real_
    ),
    # dbh3 threshold: original + both forward predictions (no SAS macro restriction here).
    dbh3_ingrowth_threshold = if_else(
      stem_visit3,
      dplyr::coalesce(.dbh3_pre_dbh, pred_dbh3_c2, pred_dbh3_c1),
      NA_real_
    )
  ) %>%
  select(
    -stem_visit1, -stem_visit2, -stem_visit3,
    -.src1_pre_dbh, -.src2_pre_dbh, -.src3_pre_dbh,
    -.dbh1_pre_dbh, -.dbh2_pre_dbh, -.dbh3_pre_dbh,
    -.species_code1_pre_dbh, -.species_code2_pre_dbh, -.species_code3_pre_dbh,
    -.habit1_pre_dbh, -.habit2_pre_dbh, -.habit3_pre_dbh,
    -.subplot1_pre_dbh, -.subplot2_pre_dbh, -.subplot3_pre_dbh,
    -.status1_pre_dbh, -.status2_pre_dbh, -.status3_pre_dbh,
    -.height1_pre_dbh, -.decayclass1_pre_dbh,
    -dbh_fit_excl_12, -dbh_fit_excl_23, -dbh_fit_excl_13
  )

ingrowth_dbh_threshold <- stems %>%
  dplyr::select(
    plot,
    ItemID,
    database_table,
    dbh1_ingrowth_threshold,
    dbh2_ingrowth_threshold,
    dbh3_ingrowth_threshold
  )

stems <- stems %>%
  dplyr::select(-dplyr::any_of(c(
    "dbh1_ingrowth_threshold",
    "dbh2_ingrowth_threshold",
    "dbh3_ingrowth_threshold"
  )))

dbh_ingrowth_model_summary <- dplyr::bind_rows(
  dbh_diag_c2_c3,
  dbh_diag_c2_c1,
  dbh_diag_c1_c2,
  dbh_diag_c1_c3,
  dbh_diag_c3_c2,
  dbh_diag_c3_c1
) %>%
  dplyr::relocate(plot, ingrowth_stage, resp, pred, species_term, model_type, n, n_species, equation, r2)

dbh_ingrowth_model_wide <- dbh_ingrowth_model_summary %>%
  tidyr::pivot_wider(
    id_cols = plot,
    names_from = ingrowth_stage,
    values_from = c(model_type, n, n_species, equation, r2, resp, pred, species_term),
    names_glue = "{ingrowth_stage}_{.value}",
    names_sort = TRUE
  )

readr::write_csv(dbh_ingrowth_model_wide, "dbh_ingrowth_models_by_plot_wide.csv", na = "")
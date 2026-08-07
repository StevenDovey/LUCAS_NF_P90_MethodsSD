# --- DBH imputation (aligned with diagnose_dataprep.R production logic; nlme::lme per plot) ---
dbh_fit_pred_plot <- function(df, resp, pred, species, out_name) {
  sp_chr <- trimws(as.character(df[[species]]))
  ok <- !is.na(df[[resp]]) & !is.na(df[[pred]]) & !is.na(df[[species]]) & nzchar(sp_chr)
  fit_df <- df[ok, , drop = FALSE]
  n_complete <- sum(ok)
  n_species <- if (n_complete == 0L) 0L else length(unique(as.character(fit_df[[species]])))

  pred_vec <- rep(NA_real_, nrow(df))
  model_type <- "none"
  eq_str <- NA_character_
  r2_val <- NA_real_

  fmt_eq_lm <- function(fit_lm_obj, resp_nm, pred_nm) {
    cf <- stats::coef(fit_lm_obj)
    b0 <- unname(cf["(Intercept)"])
    b1 <- unname(cf[[pred_nm]])
    if (length(b1) != 1L || any(is.na(c(b0, b1)))) {
      return(NA_character_)
    }
    paste0(
      resp_nm, " = ", formatC(b1, digits = 5, format = "f", flag = "#"),
      "*", pred_nm, " + ",
      formatC(b0, digits = 5, format = "f", flag = "#")
    )
  }

  diag_out <- function() {
    list(
      n_complete = n_complete,
      n_species = n_species,
      equation = eq_str,
      r2 = r2_val
    )
  }

  if (n_complete < 2L) {
    out <- df
    out[[out_name]] <- pred_vec
    return(list(model_type = model_type, df = out, diag = diag_out()))
  }

  if (n_species <= 1L) {
    fit_lm <- tryCatch(
      stats::lm(stats::reformulate(pred, response = resp), data = fit_df),
      error = function(e) NULL
    )
    if (!is.null(fit_lm)) {
      model_type <- "fixed"
      eq_str <- fmt_eq_lm(fit_lm, resp, pred)
      r2_val <- unname(summary(fit_lm)$r.squared)
      pred_vec <- tryCatch(
        as.numeric(stats::predict(fit_lm, newdata = df)),
        error = function(e) pred_vec
      )
    }
  } else {
    fdf <- fit_df
    fdf[[species]] <- factor(as.character(fdf[[species]]))
    fit_lme <- tryCatch(
      nlme::lme(
        stats::as.formula(paste(resp, "~", pred)),
        random = stats::as.formula(paste0("~ 1 | ", species)),
        data = fdf,
        method = "REML"
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_lme)) {
      model_type <- "mixed"
      fe <- nlme::fixef(fit_lme)
      b0 <- unname(fe[["(Intercept)"]])
      b1 <- unname(fe[[pred]])
      if (length(b1) == 1L && !any(is.na(c(b0, b1)))) {
        eq_str <- paste0(
          resp, " = ", formatC(b1, digits = 5, format = "f", flag = "#"),
          "*", pred, " + ",
          formatC(b0, digits = 5, format = "f", flag = "#"),
          " + u(", species, ")"
        )
      }
      nd <- df
      nd[[species]] <- factor(as.character(nd[[species]]), levels = levels(fdf[[species]]))
      pred_vec <- tryCatch(
        as.numeric(stats::predict(fit_lme, newdata = nd, level = 1)),
        error = function(e) rep(NA_real_, nrow(nd))
      )
      na1 <- is.na(pred_vec)
      if (any(na1)) {
        pred0 <- tryCatch(
          as.numeric(stats::predict(fit_lme, newdata = nd, level = 0)),
          error = function(e) NULL
        )
        if (!is.null(pred0)) pred_vec[na1] <- pred0[na1]
      }
      fitv <- tryCatch(
        as.numeric(stats::predict(fit_lme, newdata = fit_df, level = 1)),
        error = function(e) NULL
      )
      if (!is.null(fitv) && length(fitv) == nrow(fit_df)) {
        r2_val <- stats::cor(fit_df[[resp]], fitv, use = "complete.obs")^2
        if (!is.finite(r2_val)) r2_val <- NA_real_
      }
    }
    if (identical(model_type, "none")) {
      fit_lm <- tryCatch(
        stats::lm(stats::reformulate(pred, response = resp), data = fit_df),
        error = function(e) NULL
      )
      if (!is.null(fit_lm)) {
        model_type <- "fixed"
        eq_str <- fmt_eq_lm(fit_lm, resp, pred)
        r2_val <- unname(summary(fit_lm)$r.squared)
        pred_vec <- tryCatch(
          as.numeric(stats::predict(fit_lm, newdata = df)),
          error = function(e) pred_vec
        )
      }
    }
  }

  out <- df
  out[[out_name]] <- pred_vec
  list(model_type = model_type, df = out, diag = diag_out())
}

stems <- stems %>%
  mutate(
    .src1_pre_dbh = source1,
    .src2_pre_dbh = source2,
    .src3_pre_dbh = source3,
    .dbh1_pre_dbh = dbh1,
    .dbh2_pre_dbh = dbh2,
    .dbh3_pre_dbh = dbh3
  )

#########################################################################
############# Estimating DBH2 from DBH 3 (later cycle informs mid) #####
#########################################################################

stems_ingrowth2 <- stems %>%
  filter(
    habit3 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status3 == "A"
  )

preddbh2 <- stems_ingrowth2 %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh2", "dbh3", "species_code3", "pred_dbh2")
    dplyr::select(fit$df, ItemID, database_table, pred_dbh2)
  }) %>%
  dplyr::ungroup()

preddbh2 <- preddbh2 %>% arrange(plot, ItemID, database_table)
stems <- stems %>% arrange(plot, ItemID, database_table)
stems <- stems %>% left_join(preddbh2, by = c("plot", "ItemID", "database_table"))

stems <- stems %>%
  mutate(dbh2_was_missing = is.na(dbh2) & !is.na(pred_dbh2)) %>%
  mutate(
    dbh2 = if_else(dbh2_was_missing, pmax(0, pred_dbh2), dbh2),
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
      status3 == "A" & is.na(dbh2)
  ) %>%
  mutate(
    dbh2 = if_else(is_estimated, pmax(0, dbh3), dbh2),
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
############# Forward: DBH2 from DBH1 (SAS %estimate_DBH2 parity) ########
#########################################################################

stems_ingrowth2_fwd <- stems %>%
  filter(
    habit1 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status2 == "A"
  )

preddbh2_c1 <- stems_ingrowth2_fwd %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh2", "dbh1", "species_code1", "pred_dbh2_c1")
    dplyr::select(fit$df, ItemID, database_table, pred_dbh2_c1)
  }) %>%
  dplyr::ungroup()

preddbh2_c1 <- preddbh2_c1 %>% arrange(plot, ItemID, database_table)
stems <- stems %>% arrange(plot, ItemID, database_table) %>%
  left_join(preddbh2_c1, by = c("plot", "ItemID", "database_table"))

stems <- stems %>%
  mutate(dbh2_was_fwd = is.na(dbh2) & is.finite(pred_dbh2_c1)) %>%
  mutate(
    dbh2 = if_else(dbh2_was_fwd, pmax(0, pred_dbh2_c1), dbh2),
    species_code2 = if_else(
      dbh2_was_fwd & (is.na(species_code2) | species_code2 == ""),
      species_code1, species_code2
    ),
    habit2 = if_else(dbh2_was_fwd & (is.na(habit2) | habit2 == ""), habit1, habit2),
    subplot2 = if_else(dbh2_was_fwd & (is.na(subplot2) | subplot2 == ""), subplot1, subplot2),
    status2 = if_else(dbh2_was_fwd & (is.na(status2) | status2 == ""), status1, status2),
    database_table = if_else(
      dbh2_was_fwd & (is.na(database_table) | database_table == ""),
      "stems", database_table
    ),
    source2 = if_else(dbh2_was_fwd, "m_2fr1", source2)
  ) %>%
  select(-dbh2_was_fwd)

stems <- stems %>%
  mutate(
    is_estimated = habit2 %in% c("tree fern", "palm", "cabbage tree") &
      status2 == "A" & is.na(dbh2)
  ) %>%
  mutate(
    dbh2 = if_else(is_estimated, pmax(0, dbh1), dbh2),
    species_code2 = if_else(
      is_estimated & (is.na(species_code2) | species_code2 == ""),
      species_code1, species_code2
    ),
    habit2 = if_else(is_estimated & (is.na(habit2) | habit2 == ""), habit1, habit2),
    subplot2 = if_else(is_estimated & (is.na(subplot2) | subplot2 == ""), subplot1, subplot2),
    status2 = if_else(is_estimated & (is.na(status2) | status2 == ""), status1, status2),
    database_table = if_else(
      is_estimated & (is.na(database_table) | database_table == ""),
      "stems", database_table
    ),
    source2 = if_else(is_estimated, "f_1to2", source2)
  ) %>%
  select(-is_estimated)

#########################################################################
############# Forward: DBH3 from DBH2 (three-cycle extension) ##########
#########################################################################

stems_ingrowth3_fwd <- stems %>%
  filter(
    habit2 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status3 == "A"
  )

preddbh3_c2 <- stems_ingrowth3_fwd %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh3", "dbh2", "species_code2", "pred_dbh3_c2")
    dplyr::select(fit$df, ItemID, database_table, pred_dbh3_c2)
  }) %>%
  dplyr::ungroup()

preddbh3_c2 <- preddbh3_c2 %>% arrange(plot, ItemID, database_table)
stems <- stems %>% arrange(plot, ItemID, database_table) %>%
  left_join(preddbh3_c2, by = c("plot", "ItemID", "database_table"))

stems <- stems %>%
  mutate(dbh3_was_fwd = is.na(dbh3) & is.finite(pred_dbh3_c2)) %>%
  mutate(
    dbh3 = if_else(dbh3_was_fwd, pmax(0, pred_dbh3_c2), dbh3),
    species_code3 = if_else(
      dbh3_was_fwd & (is.na(species_code3) | species_code3 == ""),
      species_code2, species_code3
    ),
    habit3 = if_else(dbh3_was_fwd & (is.na(habit3) | habit3 == ""), habit2, habit3),
    subplot3 = if_else(dbh3_was_fwd & (is.na(subplot3) | subplot3 == ""), subplot2, subplot3),
    status3 = if_else(dbh3_was_fwd & (is.na(status3) | status3 == ""), status2, status3),
    database_table = if_else(
      dbh3_was_fwd & (is.na(database_table) | database_table == ""),
      "stems", database_table
    ),
    source3 = if_else(dbh3_was_fwd, "m_3fr2", source3)
  ) %>%
  select(-dbh3_was_fwd)

stems <- stems %>%
  mutate(
    is_estimated = habit3 %in% c("tree fern", "palm", "cabbage tree") &
      status3 == "A" & is.na(dbh3)
  ) %>%
  mutate(
    .dbh_src = dplyr::coalesce(dbh2, dbh1),
    dbh3 = if_else(is_estimated, pmax(0, .dbh_src), dbh3),
    species_code3 = if_else(
      is_estimated & (is.na(species_code3) | species_code3 == ""),
      dplyr::coalesce(species_code2, species_code1), species_code3
    ),
    habit3 = if_else(
      is_estimated & (is.na(habit3) | habit3 == ""),
      dplyr::coalesce(habit2, habit1), habit3
    ),
    subplot3 = if_else(
      is_estimated & (is.na(subplot3) | subplot3 == ""),
      dplyr::coalesce(subplot2, subplot1), subplot3
    ),
    status3 = if_else(
      is_estimated & (is.na(status3) | status3 == ""),
      dplyr::coalesce(status2, status1), status3
    ),
    database_table = if_else(
      is_estimated & (is.na(database_table) | database_table == ""),
      "stems", database_table
    ),
    source3 = dplyr::if_else(
      is_estimated,
      dplyr::if_else(!is.na(dbh2), "f_2to3", "f_1to3"),
      source3
    )
  ) %>%
  select(-is_estimated, -.dbh_src)

#########################################################################
############# ingrowth DBH for Cycle 1 from cycle 2 (join preds only) ###
#########################################################################

stems_ingrowth1 <- stems %>%
  filter(
    habit2 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status2 == "A"
  )

preddbh1 <- stems_ingrowth1 %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh1", "dbh2", "species_code2", "pred_dbh1")
    dplyr::select(fit$df, ItemID, database_table, pred_dbh1)
  }) %>%
  dplyr::ungroup()

preddbh1 <- preddbh1 %>% arrange(plot, ItemID, database_table)
stems <- stems %>% arrange(plot, ItemID, database_table)
stems <- stems %>% left_join(preddbh1, by = c("plot", "ItemID", "database_table"))

stems <- stems %>%
  mutate(
    is_estimated = habit2 %in% c("tree fern", "palm", "cabbage tree") &
      status2 == "A" & is.na(dbh1)
  ) %>%
  mutate(
    dbh1 = if_else(is_estimated, pmax(0, dbh2), dbh1),
    species_code1 = if_else(
      is_estimated & (is.na(species_code1) | species_code1 == ""),
      species_code2, species_code1
    ),
    habit1 = if_else(is_estimated & (is.na(habit1) | habit1 == ""), habit2, habit1),
    subplot1 = if_else(is_estimated & (is.na(subplot1) | subplot1 == ""), subplot2, subplot1),
    status1 = if_else(is_estimated & (is.na(status1) | status1 == ""), status2, status1),
    database_table = if_else(
      is_estimated & (is.na(database_table) | database_table == ""),
      "stems", database_table
    ),
    source1 = if_else(is_estimated, "f_2to1", source1)
  ) %>%
  select(-is_estimated)

###########################################################################################
#################### Populate missing Spars (Decay <3) that occurred in Cycle 2 ###########
###########################################################################################

stems <- stems %>%
  mutate(flag = !is.na(status2) & status2 == "X" & !is.na(decayclass2) & decayclass2 < 3 & is.na(dbh1)) %>%
  mutate(
    dbh1 = if_else(flag, pmax(0, dbh2), dbh1),
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
############# ingrowth estimation DBH for Cycle 1 from cycle3 ############
#########################################################################

stems_ingrowth1_3 <- stems %>%
  filter(
    habit3 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status3 == "A"
  )

preddbh1_3 <- stems_ingrowth1_3 %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh1", "dbh3", "species_code3", "pred_dbh1a")
    dplyr::select(fit$df, ItemID, database_table, pred_dbh1a)
  }) %>%
  dplyr::ungroup()

preddbh1_3 <- preddbh1_3 %>% arrange(plot, ItemID, database_table)
stems <- stems %>% arrange(plot, ItemID, database_table)
stems <- stems %>% left_join(preddbh1_3, by = c("plot", "ItemID", "database_table"))

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
    dbh1 = if_else(dbh1_was_missing, pmax(0, .pred_combined), dbh1),
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
      status3 == "A" & is.na(dbh1)
  ) %>%
  mutate(
    dbh1 = if_else(is_estimated, pmax(0, dbh3), dbh1),
    species_code1 = if_else(
      is_estimated & (is.na(species_code1) | species_code1 == ""),
      species_code3, species_code1
    ),
    habit1 = if_else(is_estimated & (is.na(habit1) | habit1 == ""), habit3, habit1),
    subplot1 = if_else(is_estimated & (is.na(subplot1) | subplot1 == ""), subplot3, subplot1),
    status1 = if_else(is_estimated & (is.na(status1) | status1 == ""), status3, status1),
    database_table = if_else(
      is_estimated & (is.na(database_table) | database_table == ""),
      "stems", database_table
    ),
    source1 = if_else(is_estimated, "f_3to1", source1)
  ) %>%
  select(-is_estimated)

#########################################################################
############# Forward fallback: DBH3 from DBH1 (after C1 ingrowth) ######
#########################################################################

stems_ingrowth3_c1 <- stems %>%
  filter(
    habit1 %in% c("canopy tree", "subcanopy tree", "shrub", "unknown"),
    status3 == "A"
  )

preddbh3_c1 <- stems_ingrowth3_c1 %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(function(.x, .y) {
    fit <- dbh_fit_pred_plot(.x, "dbh3", "dbh1", "species_code1", "pred_dbh3_c1")
    dplyr::select(fit$df, ItemID, database_table, pred_dbh3_c1)
  }) %>%
  dplyr::ungroup()

preddbh3_c1 <- preddbh3_c1 %>% arrange(plot, ItemID, database_table)
stems <- stems %>% arrange(plot, ItemID, database_table) %>%
  left_join(preddbh3_c1, by = c("plot", "ItemID", "database_table"))

stems <- stems %>%
  mutate(dbh3_was_fwd1 = is.na(dbh3) & is.finite(pred_dbh3_c1)) %>%
  mutate(
    dbh3 = if_else(dbh3_was_fwd1, pmax(0, pred_dbh3_c1), dbh3),
    species_code3 = if_else(
      dbh3_was_fwd1 & (is.na(species_code3) | species_code3 == ""),
      species_code1, species_code3
    ),
    habit3 = if_else(dbh3_was_fwd1 & (is.na(habit3) | habit3 == ""), habit1, habit3),
    subplot3 = if_else(dbh3_was_fwd1 & (is.na(subplot3) | subplot3 == ""), subplot1, subplot3),
    status3 = if_else(dbh3_was_fwd1 & (is.na(status3) | status3 == ""), status1, status3),
    database_table = if_else(
      dbh3_was_fwd1 & (is.na(database_table) | database_table == ""),
      "stems", database_table
    ),
    source3 = if_else(dbh3_was_fwd1, "m_3fr1", source3)
  ) %>%
  select(-dbh3_was_fwd1)

stems <- stems %>%
  mutate(
    source1 = if_else(!is.na(.dbh1_pre_dbh), .src1_pre_dbh, source1),
    source2 = if_else(!is.na(.dbh2_pre_dbh), .src2_pre_dbh, source2),
    source3 = if_else(!is.na(.dbh3_pre_dbh), .src3_pre_dbh, source3)
  ) %>%
  select(
    -.src1_pre_dbh, -.src2_pre_dbh, -.src3_pre_dbh,
    -.dbh1_pre_dbh, -.dbh2_pre_dbh, -.dbh3_pre_dbh
  )

##################################### Ingrowth DBH1 from C3 finished ##############################

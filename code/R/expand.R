#08.06.26 #1
# Per-hectare expansion. Each stem is scaled by the area of the subplot it was sampled in, chosen by status and size. Carbon in kg divided by area in ha and by 1000 gives t C/ha.

# Sampling area (ha) for one stem at one cycle.
sampling_area_ha <- function(status, dbh_cm, subplot, inner_area_ha, co) {
  is_ext <- !is.na(subplot) & subplot == "EXT"
  circ   <- co$plot_area_circular_ha
  out    <- rep(NA_real_, length(status))

  live <- status == "A" & !is.na(dbh_cm)
  out[live & dbh_cm >= co$dbh_circular_cm]                              <- circ
  out[live & dbh_cm >= co$dbh_min_live_inner_cm &
        dbh_cm < co$dbh_circular_cm & !is_ext]                         <- inner_area_ha[live & dbh_cm >= co$dbh_min_live_inner_cm & dbh_cm < co$dbh_circular_cm & !is_ext]

  dead <- status %in% c("X", "F", "P", "S") & !is.na(dbh_cm)
  out[dead & dbh_cm >= co$dbh_circular_cm]                              <- circ
  out[dead & dbh_cm >= co$dbh_min_dead_inner_cm &
        dbh_cm < co$dbh_circular_cm & !is_ext]                         <- inner_area_ha[dead & dbh_cm >= co$dbh_min_dead_inner_cm & dbh_cm < co$dbh_circular_cm & !is_ext]
  out
}

# Convert stem carbon (kg) to t C/ha.
carbon_per_ha <- function(carbon_kg, area_ha) {
  ifelse(is.na(area_ha) | area_ha <= 0, NA_real_, carbon_kg / area_ha / 1000)
}

# Shared interval denominator: circular if any cycle used circular, else inner.
interval_area_ha <- function(areas, co) {
  if (any(!is.na(areas) & abs(areas - co$plot_area_circular_ha) < 1e-9))
    return(co$plot_area_circular_ha)
  inner <- areas[!is.na(areas)]
  if (length(inner)) inner[1] else NA_real_
}

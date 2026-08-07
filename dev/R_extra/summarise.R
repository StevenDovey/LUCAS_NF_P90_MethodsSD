#07.06.26 19:45 NZST
# Plot-level aggregation of per-hectare carbon. Sums per-ha carbon to plot by
# interval by pool, then applies the dead-wood scaling factors to the CWD,
# standing-dead and fallen pools.

# Sum the named per-ha columns within each plot.
summarise_to_plot <- function(stems_ha, value_cols) {
  value_cols <- intersect(value_cols, names(stems_ha))
  agg <- aggregate(stems_ha[value_cols],
                   by = list(plot = stems_ha$plot),
                   FUN = function(x) sum(x, na.rm = TRUE))
  agg
}

# Apply dead-wood scaling once to the dead-wood pools.
apply_deadwood_scaling <- function(plot_summary, deadwood_cols, co,
                                   include_sampling = FALSE) {
  factor <- co$cwd_belowground_uplift *
            (if (include_sampling) co$cwd_sampling_multiplier else 1)
  for (cl in intersect(deadwood_cols, names(plot_summary)))
    plot_summary[[cl]] <- plot_summary[[cl]] * factor
  plot_summary
}

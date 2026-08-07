#08.06.26 #3
# Output checks on the plot summary, run after run.R. Catches results that should not occur: negative pool carbon, an implausibly high stock, negative ingrowth, or a net change larger than the standing stock. Writes output/diagnostics/validation_post.csv.
# Run from the repository root, after run.R:  Rscript diagnostics/run_validate_post.R
root <- getwd()
out <- read.csv(file.path(root, "output/plot_carbon_summary.csv"), check.names = FALSE)
agb_ceiling_t_ha <- 700   # plausible upper bound on above-ground stock

flags <- list()
add <- function(mask, rule, severity, detail) {
  mask[is.na(mask)] <- FALSE
  if (!any(mask)) return(NULL)
  data.frame(plot = out$plot[mask], rule = rule, severity = severity, detail = detail,
             stringsAsFactors = FALSE)
}
num <- function(c) suppressWarnings(as.numeric(out[[c]]))
has <- function(c) c %in% names(out)

# Negative stocks (above-ground, below-ground, dead-wood) must not occur.
for (k in 1:3) for (p in c("agb_ha","bgb_ha","spars_ha","fallen_ha","CWD_ha")) {
  cc <- paste0(p, ".", k)
  if (has(cc)) flags <- c(flags, list(add(num(cc) < 0, "negative_pool_carbon", "error",
                                          sprintf("%s negative", cc))))
}
# Implausibly high above-ground stock.
for (k in 1:3) { cc <- paste0("agb_ha.", k)
  if (has(cc)) flags <- c(flags, list(add(num(cc) > agb_ceiling_t_ha, "agb_above_ceiling", "warning",
                                          sprintf("%s exceeds %g t C/ha", cc, agb_ceiling_t_ha)))) }
# Ingrowth should not be negative.
for (iv in c("1_2","2_3","1_3")) { cc <- paste0("AGB_ingrowth_", iv)
  if (has(cc)) flags <- c(flags, list(add(num(cc) < 0, "ingrowth_negative", "warning",
                                          sprintf("%s negative", cc)))) }
# Net change larger in magnitude than the standing stock points at an input error.
for (iv in c("1_2","2_3")) {
  cc <- paste0("AGB_change_", iv); st <- paste0("agb_ha.", substr(iv,1,1))
  if (has(cc) && has(st))
    flags <- c(flags, list(add(abs(num(cc)) > num(st) & num(st) > 0, "change_exceeds_stock", "warning",
                               sprintf("%s magnitude exceeds the standing stock", cc))))
}

res <- do.call(rbind, flags)
if (is.null(res)) res <- data.frame(plot = character(), rule = character(),
                                    severity = character(), detail = character())
dc <- file.path(root, "output", "diagnostics")
dir.create(dc, recursive = TRUE, showWarnings = FALSE)
write.csv(res, file.path(dc, "validation_post.csv"), row.names = FALSE, na = "")
cat(sprintf("plots checked: %d   post-run flags: %d (errors %d, warnings %d)\n",
            nrow(out), nrow(res), sum(res$severity == "error"), sum(res$severity == "warning")))
cat("\nflags by rule:\n"); print(as.data.frame(table(rule = res$rule)))

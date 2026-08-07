#08.06.26 #1
# Diameter gap fill, form selected by cross-validation. A missing diameter is filled from the nearest measured diameter of the same stem, via a species random-intercept mixed model fitted on measured consecutive-cycle pairs. The forward model predicts a later diameter from an earlier one, the backward model the reverse. Unseen species fall to the population mean. Filled values are tagged in `dbh_source`. Input: long table, one row per plot/tag/cycle, columns plot, tag, cycle, dbh, species_code, dbh_source.
library(nlme)

# Fit one direction: pair each cycle with the next, keep measured pairs, fit later diameter on earlier (or the reverse when earlier_is_from is FALSE).
.fit_pair_model <- function(d, earlier_is_from, lc) {
  cyc <- sort(unique(d$cycle))
  pr <- do.call(rbind, lapply(cyc[-length(cyc)], function(a) {
    ea <- d[d$cycle == a, ]; eb <- d[d$cycle == a + 1, ]
    m <- merge(ea[c("tag","dbh","species_code")], eb[c("tag","dbh")], by = "tag",
               suffixes = c("_e","_l"))
    m[which(!is.na(m$dbh_e) & !is.na(m$dbh_l) & m$dbh_e > 0 & m$dbh_l > 0 & nzchar(m$species_code)), ]
  }))
  if (earlier_is_from) { from <- pr$dbh_e; to <- pr$dbh_l } else { from <- pr$dbh_l; to <- pr$dbh_e }
  fit <- lme(to ~ from, random = ~1 | species_code,
             data = data.frame(to = to, from = from, species_code = pr$species_code), control = lc)
  list(b = fixef(fit), re = setNames(ranef(fit)[[1]], rownames(ranef(fit))))
}

fill_dbh_gaps <- function(stems, lc = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)) {
  fwd <- .fit_pair_model(stems, earlier_is_from = TRUE,  lc)
  bwd <- .fit_pair_model(stems, earlier_is_from = FALSE, lc)
  pred1 <- function(mod, sp, from) {
    r <- mod$re[sp]; r[is.na(r)] <- 0
    mod$b[1] + mod$b[2] * from + r
  }
  parts <- split(seq_len(nrow(stems)), paste(stems$plot, stems$tag))
  for (idx in parts) {
    g <- stems[idx, , drop = FALSE]
    o <- order(g$cycle); g <- g[o, ]; ix <- idx[o]
    if (all(!is.na(g$dbh)) || all(is.na(g$dbh))) next
    meas <- which(!is.na(g$dbh))
    for (k in which(is.na(g$dbh))) {
      j <- meas[which.min(abs(meas - k))]
      mod <- if (j < k) fwd else bwd
      stems$dbh[ix[k]]        <- pred1(mod, g$species_code[j], g$dbh[j])
      stems$dbh_source[ix[k]] <- "imputed_model"
    }
  }
  stems
}

#08.06.26 #1
# Two-stage height model, form selected by cross-validation.
#   Stage 1: species random-intercept mixed model of log(height - 1.35) on the diameter term diameter^(-0.3), fitted across all measured stems.
#   Stage 2: per-plot multiplicative bias correction from that plot's measured stems. Unseen species fall to the population mean, unseen plots to a correction of one.
# Returns the input table with a `predht` column for woody stems that carry a diameter. Inputs: status, habit, dbh, height, lean, species_code, plot.
library(nlme)

predict_height_two_stage <- function(stems, lc = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)) {
  woody <- c("canopy tree", "subcanopy tree", "shrub", "unknown")
  is_woody <- stems$habit %in% woody
  stems$x3 <- stems$dbh^(-0.3)

  train <- stems[which(is_woody & stems$status == "A" & !is.na(stems$height) &
                 stems$height > 1.35 & is.na(stems$lean) & stems$dbh > 0 &
                 nzchar(stems$species_code)), ]
  train$lh <- log(train$height - 1.35)

  m  <- lme(lh ~ x3, random = ~1 | species_code, data = train, control = lc)
  b  <- fixef(m)
  re <- ranef(m)[[1]]; names(re) <- rownames(ranef(m))

  eta_pred <- function(sp, x3) {
    r <- re[sp]; r[is.na(r)] <- 0
    b[1] + b[2] * x3 + r
  }
  train$p135 <- exp(eta_pred(train$species_code, train$x3))
  bias <- tapply(train$height - 1.35, train$plot, sum) /
          tapply(train$p135,          train$plot, sum)

  fillable <- is_woody & !is.na(stems$dbh) & stems$dbh > 0
  p135 <- exp(eta_pred(stems$species_code[fillable], stems$x3[fillable]))
  bp <- bias[stems$plot[fillable]]; bp[is.na(bp)] <- 1
  stems$predht <- NA_real_
  stems$predht[fillable] <- 1.35 + p135 * bp
  stems
}

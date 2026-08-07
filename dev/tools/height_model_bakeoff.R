#07.06.26 21:00 NZST
# Cross-validated bake-off of candidate height models. Task: predict measured
# height from measured diameter, truth is the field height of live,
# non-leaning woody stems. Folds split by stem, and every candidate is scored
# on the same folds for mean bias, RMSE and MAE. The reference file's own
# predicted height is also scored in-sample as a benchmark.
library(nlme)
root <- "/home/user/SAS2R"
lc <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)

M <- "sasout/SAS_measurements-with-carbonV10.csv"
keep <- c("ItemID","plot","DBH","height","LeanAngle","predht","species_code","habit","status","source")
d <- read.csv(file.path(root, M), stringsAsFactors = FALSE, check.names = FALSE)[keep]
for (c in c("DBH","height","LeanAngle","predht")) d[[c]] <- as.numeric(d[[c]])
woody <- c("canopy tree","subcanopy tree","shrub","unknown")
d <- d[d$status == "A" & d$habit %in% woody & d$source == "measured" &
       d$DBH > 0 & !is.na(d$height) & d$height > 1.35 & is.na(d$LeanAngle) &
       nzchar(d$species_code), ]
H <- data.frame(stem = d$ItemID, plot = d$plot, sp = d$species_code,
                dbh = d$DBH, H = d$height, ref = d$predht)
H <- H[!is.na(H$plot) & !is.na(H$dbh), ]
H$x3 <- H$dbh^(-0.3); H$h135 <- H$H - 1.35; H$lh <- log(H$h135)
cat(sprintf("measured non-leaning woody heights: %d  stems: %d  plots: %d\n",
            nrow(H), length(unique(H$stem)), length(unique(H$plot))))

cand <- list(
  linear      = function(tr){ m <- lm(H ~ dbh, tr); function(te) predict(m, te) },
  loglog      = function(tr){ m <- lm(log(H) ~ log(dbh), tr); function(te) exp(predict(m, te)) },
  lucas_fixed = function(tr){ m <- lm(lh ~ x3, tr); function(te) 1.35 + exp(predict(m, te)) },
  power_nls   = function(tr){ m <- nls(H ~ 1.35 + a * dbh^b, tr, start = list(a = 2, b = 0.5),
                                       algorithm = "port", lower = c(0.01, 0.01), upper = c(60, 3))
                              function(te) predict(m, te) },
  ri_species  = function(tr){ m <- lme(lh ~ x3, random = ~1 | sp, data = tr, control = lc)
                              b <- fixef(m); re <- ranef(m)[[1]]; names(re) <- rownames(ranef(m))
                              function(te){ r <- re[te$sp]; r[is.na(r)] <- 0; 1.35 + exp(b[1] + b[2]*te$x3 + r) } },
  rs_species  = function(tr){ m <- lme(lh ~ x3, random = ~1 + x3 | sp, data = tr, control = lc)
                              b <- fixef(m); re <- ranef(m); ri <- re[[1]]; rs <- re[[2]]
                              names(ri) <- rownames(re); names(rs) <- rownames(re)
                              function(te){ a <- ri[te$sp]; s <- rs[te$sp]; a[is.na(a)] <- 0; s[is.na(s)] <- 0
                                1.35 + exp(b[1] + b[2]*te$x3 + a + s*te$x3) } },
  ref_2stage  = function(tr){ m <- lme(lh ~ x3, random = ~1 | sp, data = tr, control = lc)
                              b <- fixef(m); re <- ranef(m)[[1]]; names(re) <- rownames(ranef(m))
                              eta <- b[1] + b[2]*tr$x3 + ifelse(is.na(re[tr$sp]), 0, re[tr$sp])
                              p135 <- exp(eta)
                              bias <- tapply(tr$h135, tr$plot, sum) / tapply(p135, tr$plot, sum)
                              function(te){ r <- re[te$sp]; r[is.na(r)] <- 0
                                bp <- bias[te$plot]; bp[is.na(bp)] <- 1
                                1.35 + exp(b[1] + b[2]*te$x3 + r) * bp } }
)

set.seed(1)
stems <- unique(H$stem)
fold <- setNames(sample(rep(1:5, length.out = length(stems))), stems)
H$fold <- fold[as.character(H$stem)]

res <- list()
for (nm in names(cand)) {
  pe <- numeric(0); ye <- numeric(0)
  for (k in 1:5) {
    tr <- H[H$fold != k, ]; te <- H[H$fold == k, ]
    pred <- cand[[nm]](tr)(te)
    pe <- c(pe, pred); ye <- c(ye, te$H)
  }
  e <- pe - ye
  res[[nm]] <- c(bias = mean(e), RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)))
}
es <- H$ref - H$H
res[["reference_predht_insample"]] <- c(bias = mean(es), RMSE = sqrt(mean(es^2)), MAE = mean(abs(es)))

tab <- do.call(rbind, res)
cat("\nHeight bake-off (m), 5-fold CV, truth = measured height\n")
print(round(tab[order(abs(tab[,"bias"])), ], 4))

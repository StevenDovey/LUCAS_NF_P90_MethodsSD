#07.06.26 20:30 NZST
# Cross-validated bake-off of candidate diameter gap-fill models. Task: predict
# a stem's measured diameter at a cycle from its measured diameter at the
# adjacent cycle, truth is the field measurement. Folds split by stem, and
# every candidate is scored on the same folds for mean bias, RMSE and MAE.
library(nlme)
root <- "/home/user/SAS2R"
lc <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)

M <- "sasout/SAS_measurements-with-carbonV10.csv"
keep <- c("ItemID","plot","measurement","DBH","species_code","year","source","status")
d <- read.csv(file.path(root, M), stringsAsFactors = FALSE, check.names = FALSE)[keep]
d$DBH <- as.numeric(d$DBH); d$measurement <- as.numeric(d$measurement); d$year <- as.numeric(d$year)
d <- d[d$source == "measured" & d$status == "A" & d$DBH > 0 & d$measurement %in% 1:3, ]

# consecutive-cycle pairs (from -> to)
mk_pairs <- function(a, b) {
  ea <- d[d$measurement == a, ]; eb <- d[d$measurement == b, ]
  m <- merge(ea[c("ItemID","DBH","year","species_code","plot")],
             eb[c("ItemID","DBH","year")], by = "ItemID", suffixes = c("_from","_to"))
  data.frame(stem = m$ItemID, plot = m$plot, sp = m$species_code,
             x = m$DBH_from, y = m$DBH_to, dt = m$year_to - m$year_from)
}
P <- rbind(mk_pairs(1,2), mk_pairs(2,3))
P <- P[P$dt > 0, ]
n0 <- nrow(P)
keep_rows <- !is.na(P$x) & !is.na(P$y) & !is.na(P$dt) & !is.na(P$plot) & !is.na(P$sp) & nzchar(P$sp)
P <- P[keep_rows, ]
cat(sprintf("pairs: %d (dropped %d incomplete)  stems: %d  plots: %d\n",
            nrow(P), n0 - nrow(P), length(unique(P$stem)), length(unique(P$plot))))

# candidates: each returns a fitted predictor given a training frame
cand <- list(
  carry        = function(tr) function(te) te$x,
  incr_global  = function(tr){ i <- median((tr$y - tr$x) / tr$dt); function(te) te$x + i * te$dt },
  incr_species = function(tr){ i <- tapply((tr$y - tr$x) / tr$dt, tr$sp, median); g <- median((tr$y - tr$x) / tr$dt)
                               function(te){ v <- i[te$sp]; v[is.na(v)] <- g; te$x + v * te$dt } },
  lm_global    = function(tr){ m <- lm(y ~ x, tr); function(te) predict(m, te) },
  loglinear    = function(tr){ m <- lm(log(y) ~ log(x), tr); function(te) exp(predict(m, te)) },
  lm_per_plot  = function(tr){ gco <- coef(lm(y ~ x, tr))
                               pc <- lapply(split(tr, tr$plot), function(g) if (nrow(g) >= 3) coef(lm(y ~ x, g)) else gco)
                               cf <- do.call(rbind, pc)
                               function(te){ idx <- match(te$plot, rownames(cf))
                                 ic <- cf[idx, 1]; sl <- cf[idx, 2]
                                 ic[is.na(idx)] <- gco[1]; sl[is.na(idx)] <- gco[2]; ic + sl * te$x } },
  ri_species   = function(tr){ m <- lme(y ~ x, random = ~1 | sp, data = tr, control = lc); b <- fixef(m)
                               re <- ranef(m)[[1]]; names(re) <- rownames(ranef(m))
                               function(te){ r <- re[te$sp]; r[is.na(r)] <- 0; b[1] + b[2] * te$x + r } },
  ri_plot      = function(tr){ m <- lme(y ~ x, random = ~1 | plot, data = tr, control = lc); b <- fixef(m)
                               re <- ranef(m)[[1]]; names(re) <- rownames(ranef(m))
                               function(te){ r <- re[te$plot]; r[is.na(r)] <- 0; b[1] + b[2] * te$x + r } }
)

set.seed(1)
stems <- unique(P$stem)
fold <- setNames(sample(rep(1:5, length.out = length(stems))), stems)
P$fold <- fold[as.character(P$stem)]

res <- list()
for (nm in names(cand)) {
  pe <- numeric(0); ye <- numeric(0)
  for (k in 1:5) {
    tr <- P[P$fold != k, ]; te <- P[P$fold == k, ]
    pred <- cand[[nm]](tr)(te)
    pe <- c(pe, pred); ye <- c(ye, te$y)
  }
  e <- pe - ye
  res[[nm]] <- c(bias = mean(e), RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)))
}
tab <- do.call(rbind, res)
cat("\nDBH bake-off (cm), 5-fold CV, truth = measured diameter\n")
print(round(tab[order(abs(tab[,"bias"])), ], 4))

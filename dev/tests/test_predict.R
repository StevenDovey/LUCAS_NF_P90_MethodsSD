#08.06.26 00:20 NZST
# Property tests for the diameter and height models. Builds one synthetic stand
# with a known shape, then checks predictions behave as expected (finite, above
# breast height, rising with diameter, close to the true values).
# Run from the repository root:  Rscript dev/tests/test_predict.R
root <- if (file.exists("R/coefficients.R")) "R" else file.path("..", "R")
for (f in c("coefficients.R","allometry.R","predict_dbh.R","predict_height.R")) source(file.path(root, f))
co <- carbon_coefficients()
pass <- 0L; fail <- 0L
check <- function(label, cond) { if (isTRUE(cond)) { pass <<- pass + 1L; cat(sprintf("  ok  %s\n", label)) }
  else { fail <<- fail + 1L; cat(sprintf("FAIL  %s\n", label)) } }

set.seed(7)
nsp <- 4; nstem <- 400
species <- paste0("SP", sample(seq_len(nsp), nstem, replace = TRUE))
plot    <- paste0("P", sample(1:8, nstem, replace = TRUE))
dbh1    <- round(runif(nstem, 5, 60), 1)
grow    <- runif(nstem, 0.5, 3)
dbh2    <- dbh1 + grow
# true height from a species allometry with noise
a <- c(SP1 = 2.2, SP2 = 2.6, SP3 = 1.9, SP4 = 2.4)
h1 <- 1.35 + exp(a[species] + (-2.6) * dbh1^(-0.3)) + rnorm(nstem, 0, 0.3)

cat("== height model ==\n")
sd <- data.frame(plot = plot, tag = seq_len(nstem), status = "A",
                 habit = "canopy tree", species_code = species, dbh = dbh1,
                 height = h1, lean = NA_real_, stringsAsFactors = FALSE)
hd <- predict_height_two_stage(sd)
ok <- is.finite(hd$predht)
check("predht finite for woody with dbh", all(ok))
check("predht above breast height", all(hd$predht[ok] > 1.35))
g1 <- ok & hd$species_code == "SP1"
check("predht rises with dbh within a species",
      cor(hd$dbh[g1], hd$predht[g1], method = "spearman") > 0.8)
check("predht near measured (median abs err < 1.5 m)",
      median(abs(hd$predht - hd$height)) < 1.5)

cat("== diameter gap fill ==\n")
long <- rbind(
  data.frame(plot = plot, tag = seq_len(nstem), cycle = 1L, dbh = dbh1, species_code = species),
  data.frame(plot = plot, tag = seq_len(nstem), cycle = 2L, dbh = dbh2, species_code = species))
long$dbh_source <- "measured"
hole <- long$cycle == 2 & long$tag %% 5 == 0
truth <- long$dbh[hole]
long$dbh[hole] <- NA; long$dbh_source[hole] <- NA
filled <- fill_dbh_gaps(long)
got <- filled$dbh[hole]
check("all gaps filled", all(!is.na(got)))
check("filled tagged imputed_model", all(filled$dbh_source[hole] == "imputed_model"))
check("filled diameter near truth (median abs err < 1.5 cm)", median(abs(got - truth)) < 1.5)
check("filled diameter exceeds earlier (growth, >= 90%)",
      mean(got >= long$dbh[long$cycle == 1][match(filled$tag[hole], long$tag[long$cycle == 1])]) > 0.9)

cat(sprintf("\n%d passed, %d failed\n", pass, fail))
if (fail > 0) quit(status = 1)

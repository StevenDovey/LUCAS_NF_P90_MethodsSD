#07.06.26 22:15 NZST
# Explains why row counts differ between the reference output, the old R
# output, and the raw input. Counts rows, distinct stems, rows per cycle, and
# the keys present in one source but not the other.
root <- "/home/user/SAS2R"

refm <- read.csv(file.path(root, "sasout/SAS_measurements-with-carbonV10.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)[c("ItemID","measurement","status","source","DBH")]
oldm <- read.csv(file.path(root, "rfiles/measurements_C_R.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)[c("ItemID","cycle","status","source","dbh")]
stm  <- read.csv(file.path(root, "mfedata/MFESensitivity_Stems_20260209_reduced_dbl_removed3.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)[c("ItemID","AliveState","DiameterValue")]

cat(sprintf("raw stems input rows : %d  distinct ItemID: %d\n", nrow(stm), length(unique(stm$ItemID))))
cat(sprintf("reference measurements rows: %d  distinct ItemID: %d\n", nrow(refm), length(unique(refm$ItemID))))
cat(sprintf("old R measurements   : %d  distinct ItemID: %d\n", nrow(oldm), length(unique(oldm$ItemID))))

cat("\nreference rows per measurement:\n"); print(table(refm$measurement, useNA = "ifany"))
cat("\nold R rows per cycle:\n");     print(table(oldm$cycle, useNA = "ifany"))

cat("\nreference source breakdown:\n");   print(table(refm$source, useNA = "ifany"))
cat("\nold R source breakdown:\n"); print(table(oldm$source, useNA = "ifany"))

# rows that carry no diameter (candidate padded / structural rows)
cat(sprintf("\nreference rows with NA DBH : %d\n", sum(is.na(refm$DBH))))
cat(sprintf("old R rows with NA dbh: %d\n", sum(is.na(oldm$dbh))))

# key overlap by (ItemID, cycle)
ks <- paste(refm$ItemID, refm$measurement)
ko <- paste(oldm$ItemID, oldm$cycle)
cat(sprintf("\n(ItemID,cycle) keys: reference %d  old R %d\n", length(unique(ks)), length(unique(ko))))
cat(sprintf("in old R but not reference: %d\n", length(setdiff(unique(ko), unique(ks)))))
cat(sprintf("in reference but not old R: %d\n", length(setdiff(unique(ks), unique(ko)))))

# of the old-R-only keys, how many carry a real diameter
extra <- setdiff(unique(ko), unique(ks))
oldm$k <- ko
ex <- oldm[oldm$k %in% extra, ]
cat(sprintf("old-R-only keys with NA dbh: %d of %d\n", sum(is.na(ex$dbh)), nrow(ex)))
cat("old-R-only status breakdown:\n"); print(table(ex$status, useNA = "ifany"))

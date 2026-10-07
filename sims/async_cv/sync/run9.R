## Rscript run9.R <gen> <link> <n> <from> <to>
source("cov9.R"); a <- commandArgs(TRUE); n <- as.integer(a[3]); reps <- as.integer(a[4]):as.integer(a[5])
res <- lapply(reps, function(r) tryCatch(one_rep9(n, a[2], r, a[1]), error = function(e) conditionMessage(e)))
saveRDS(res, sprintf("results9/%s_%s_%d_%d.rds", a[1], a[2], n, reps[1]))

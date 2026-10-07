## Rscript run7.R <link> <n> <from> <to>
source("cov7.R"); a <- commandArgs(TRUE); n <- as.integer(a[2]); reps <- as.integer(a[3]):as.integer(a[4])
res <- lapply(reps, function(r) tryCatch(one_rep7(n, a[1], r), error = function(e) conditionMessage(e)))
saveRDS(res, sprintf("results7/%s_%d_%d.rds", a[1], n, reps[1]))

## Rscript run6.R <link> <n> <from> <to>
source("cov6.R"); a <- commandArgs(TRUE); n <- as.integer(a[2]); reps <- as.integer(a[3]):as.integer(a[4])
res <- lapply(reps, function(r) tryCatch({ .last_lo <<- NA_real_; o <- one_rep_link(n, a[1], r); cbind(o, h_lo = .last_lo) },
  error = function(e) conditionMessage(e)))
saveRDS(res, sprintf("results6/%s_%d_%d.rds", a[1], n, reps[1]))

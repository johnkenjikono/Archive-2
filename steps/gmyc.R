# Single-threshold GMYC (splits::gmyc) on a rooted Newick tree.
# Usage: Rscript gmyc.R tree.nwk out_dir [outgroup]
# Writes out_dir/gmyc_clusters.csv (genome,cluster) and out_dir/gmyc_summary.csv.
suppressPackageStartupMessages({ library(ape); library(splits) })

args <- commandArgs(trailingOnly = TRUE)
tree_path <- args[1]; out_dir <- args[2]; outgroup <- if (length(args) > 2) args[3] else NA

tr <- read.tree(tree_path)
if (!is.na(outgroup)) {
  if (!outgroup %in% tr$tip.label) stop("outgroup '", outgroup, "' is not a tip of ", tree_path)
  tr <- drop.tip(tr, outgroup)
}
if (Ntip(tr) < 3) stop("GMYC needs at least 3 ingroup tips, got ", Ntip(tr))

set.seed(1)  # multi2di resolves polytomies at random
tr <- multi2di(tr)
# chronos and GMYC fail on zero-length branches (polytomies, collapsed clones): give them a tiny length.
min_len <- max(tr$edge.length) * 1e-6
tr$edge.length[tr$edge.length < min_len] <- min_len
invisible(capture.output(tr <- chronos(tr, lambda = 1, quiet = TRUE)))
class(tr) <- "phylo"

res <- gmyc(tr, method = "single", quiet = TRUE)
spec <- spec.list(res)
write.csv(data.frame(genome = spec$sample_name, cluster = spec$GMYC_spec),
          file.path(out_dir, "gmyc_clusters.csv"), row.names = FALSE, quote = FALSE)

best <- which.max(res$likelihood)
in_ci <- res$likelihood > res$likelihood[best] - 2  # 2 log-likelihood units, as summary.gmyc
lr <- 2 * (res$likelihood[best] - res$likelihood[1])
write.csv(data.frame(
  tips = Ntip(tr),
  entities = res$entity[best], entities_ci_low = min(res$entity[in_ci]), entities_ci_high = max(res$entity[in_ci]),
  clusters = res$cluster[best], clusters_ci_low = min(res$cluster[in_ci]), clusters_ci_high = max(res$cluster[in_ci]),
  null_loglik = res$likelihood[1], gmyc_loglik = res$likelihood[best],
  likelihood_ratio = lr, lr_test_p = 1 - pchisq(lr, 2), threshold_time = res$threshold.time[best]
), file.path(out_dir, "gmyc_summary.csv"), row.names = FALSE)

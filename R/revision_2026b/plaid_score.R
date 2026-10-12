# PLAID-style gene-set score, the fast scorer of spatial-scribe (src/spatialscribe/analysis/programs.py,
# score_programs): log1p(CP10k) expression, each gene z-scored across all scored cells (population sd,
# sd 0 -> 1), score = mean z over the set's genes. No ranking, no permutation.
# Sparse throughout: mean_g (x_g - mu_g) / sd_g = t(X) %*% W - colSums(mu * W), with W[g, k] = 1 / (sd_g |set_k|),
# so the genes x cells matrix is never densified.
# counts: genes x cells (dgCMatrix, raw counts) of ALL cells to be compared - z-scores are relative to that pool.
plaid_score <- function(counts, sets) {
  stopifnot(inherits(counts, "dgCMatrix"), !is.null(rownames(counts)), !is.null(colnames(counts)))
  sets <- lapply(sets, function(g) intersect(unique(g), rownames(counts)))
  stopifnot("a gene set has no gene in the matrix" = all(lengths(sets) > 0))
  genes <- sort(unique(unlist(sets)))
  X <- counts[genes, , drop = FALSE] %*% Matrix::Diagonal(x = 1e4 / Matrix::colSums(counts))
  X@x <- log1p(X@x)
  dimnames(X) <- list(genes, colnames(counts))          # %*% Diagonal() drops the cell names
  n <- ncol(X)
  mu <- Matrix::rowSums(X) / n
  sdv <- sqrt(pmax(Matrix::rowSums(X^2) / n - mu^2, 0))
  sdv[sdv == 0] <- 1
  W <- sapply(sets, function(g) ifelse(genes %in% g, 1 / (sdv * length(g)), 0))
  W <- matrix(W, nrow = length(genes), dimnames = list(genes, names(sets)))
  S <- as.matrix(Matrix::t(X) %*% W) - matrix(colSums(mu * W), n, ncol(W), byrow = TRUE)
  dimnames(S) <- list(colnames(counts), names(sets))
  S
}

# Self-check against the direct dense computation on a toy matrix: Rscript plaid_score.R --selftest
if (sys.nframe() == 0 && identical(commandArgs(trailingOnly = TRUE), "--selftest")) {
  set.seed(1)
  m <- Matrix::rsparsematrix(50, 300, 0.3, rand.x = function(k) rpois(k, 5) + 1)
  m <- as(m, "CsparseMatrix"); dimnames(m) <- list(paste0("g", 1:50), paste0("c", 1:300))
  sets <- list(A = paste0("g", 1:10), B = paste0("g", c(5:20, 99)))
  fast <- plaid_score(m, sets)
  ln <- log1p(t(t(as.matrix(m)) / colSums(as.matrix(m)) * 1e4))
  z <- (ln - rowMeans(ln)) / sqrt(rowMeans((ln - rowMeans(ln))^2))
  slow <- sapply(sets, function(g) colMeans(z[intersect(g, rownames(z)), , drop = FALSE]))
  stopifnot(max(abs(fast - slow)) < 1e-10, all(rownames(fast) == colnames(m)))
  cat("plaid_score selftest OK\n")
}

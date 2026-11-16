#!/usr/bin/env Rscript
# =============================================================================
# rnaseq_pipeline.R
# End-to-end bulk RNA-Seq differential expression module (simulated data)
#
# Target: R >= 4.4
# Dependencies (all CRAN, all available as Windows binaries -> no Rtools):
#   ggplot2, ggrepel, pheatmap, matrixStats
#   + MASS (ships with R as a "recommended" package)
#
# Design goals
#   * Reproducible: one seed, explicit parameters, sessionInfo() saved.
#   * Modular: every step is a function with a single responsibility.
#   * Verifiable: data are simulated with KNOWN ground truth, so the
#     statistical procedure can be benchmarked (precision / recall / FDR).
#
# IMPORTANT: this is a transparent, didactic re-implementation of the core
# ideas behind DESeq2 (Love et al., 2014). For real data use DESeq2/edgeR/
# limma-voom. See the commented real-data block at the end of the file.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(matrixStats)
  library(MASS)        # negative.binomial() family
})

# ----------------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------------
CONFIG <- list(
  seed        = 42L,
  n_genes     = 6000L,
  n_per_group = 6L,
  frac_de     = 0.10,   # fraction of truly differentially expressed genes
  fdr_cutoff  = 0.05,
  lfc_cutoff  = 1,      # |log2FC| threshold for "biologically relevant"
  n_top_var   = 500L,   # genes used for PCA
  n_labels    = 12L,    # genes labelled in the volcano plot
  out_dir     = "results"
)

# ----------------------------------------------------------------------------
# 1. Data simulation (Negative Binomial with a mean-dispersion trend)
# ----------------------------------------------------------------------------
# Why Negative Binomial? RNA-Seq counts are overdispersed relative to Poisson
# (Var = mu + phi * mu^2). Biological replicates vary more than pure sampling
# noise, so a Poisson model would give anti-conservative p-values.
simulate_rnaseq <- function(n_genes, n_per_group, frac_de, seed) {
  set.seed(seed)
  n_samples <- 2L * n_per_group
  genes     <- sprintf("GENE%05d", seq_len(n_genes))
  condition <- factor(rep(c("Control", "Treated"), each = n_per_group),
                      levels = c("Control", "Treated"))
  samples   <- paste0(ifelse(condition == "Control", "Ctrl", "Trt"), "_",
                      sprintf("%02d", rep(seq_len(n_per_group), 2)))
  
  # Log-normal baseline expression: few highly expressed genes, many low ones
  base_mean <- rlnorm(n_genes, meanlog = 5, sdlog = 1.6)
  
  # True dispersion follows the typical decreasing trend phi ~ a0 + a1/mu,
  # plus gene-wise scatter around it.
  true_disp <- (0.05 + 1 / base_mean) * rlnorm(n_genes, 0, 0.3)
  
  # Sample-specific sequencing depth (library size) differences
  size_factors <- rlnorm(n_samples, 0, 0.25)
  
  # Ground-truth effects: 10% of genes, up/down with |log2FC| >= 0.75
  n_de   <- round(frac_de * n_genes)
  de_idx <- sample.int(n_genes, n_de)
  true_lfc <- numeric(n_genes)
  true_lfc[de_idx] <- sample(c(-1, 1), n_de, replace = TRUE) *
    (0.75 + rexp(n_de, rate = 1.5))
  
  grp <- as.integer(condition) - 1L                   # 0 = Control, 1 = Treated
  mu  <- outer(base_mean, size_factors) * 2^outer(true_lfc, grp)
  size <- matrix(1 / true_disp, nrow = n_genes, ncol = n_samples)
  counts <- matrix(rnbinom(length(mu), mu = mu, size = size),
                   nrow = n_genes, dimnames = list(genes, samples))
  
  # Synthetic gene sets (NOT real pathways): planted signal (so enrichment can be validated) + null sets
  up_genes   <- genes[de_idx[true_lfc[de_idx] > 0]]
  down_genes <- genes[de_idx[true_lfc[de_idx] < 0]]
  gene_sets <- list(
    SYNTHETIC_UP_SET = c(sample(up_genes, 40),
                            sample(setdiff(genes, up_genes), 60)),
    SYNTHETIC_DOWN_SET = c(sample(down_genes, 40),
                            sample(setdiff(genes, down_genes), 60))
  )
  for (i in 1:8) gene_sets[[sprintf("RANDOM_SET_%d", i)]] <- sample(genes, 100)
  
  list(counts = counts,
       coldata = data.frame(condition = condition, row.names = samples),
       truth = data.frame(gene = genes, true_lfc = true_lfc,
                          is_de = true_lfc != 0, row.names = genes),
       gene_sets = gene_sets)
}

# ----------------------------------------------------------------------------
# 2. Normalisation: median-of-ratios size factors (Anders & Huber, 2010)
# ----------------------------------------------------------------------------
# Total-count scaling is biased when a few strongly regulated genes absorb a
# large share of reads. The median-of-ratios assumes most genes are NOT DE, so
# the median ratio to a pseudo-reference (geometric mean) is robust to them.
estimate_size_factors <- function(counts) {
  log_geo_means <- rowMeans(log(counts))
  ok <- is.finite(log_geo_means)           # genes with no zeros in any sample
  stopifnot(sum(ok) > 100)
  apply(log(counts[ok, , drop = FALSE]) - log_geo_means[ok], 2,
        function(x) exp(median(x)))
}

# ----------------------------------------------------------------------------
# 3. Dispersion estimation with shrinkage towards the mean-dispersion trend
# ----------------------------------------------------------------------------
# With n = 6 per group, gene-wise dispersion estimates are very noisy. We borrow
# information across genes: fit the trend phi ~ a0 + a1/mu and shrink each
# gene-wise estimate towards it (empirical-Bayes idea, simplified: fixed prior
# degrees of freedom, weighted average on the log scale).
estimate_dispersions <- function(counts, size_factors, condition,
                                 prior_df = 4) {
  norm <- sweep(counts, 2, size_factors, "/")
  mu   <- rowMeans(norm)
  df   <- ncol(counts) - nlevels(condition)
  
  # Within-group residuals -> pooled variance (variance NOT driven by the effect)
  resid <- norm
  for (lv in levels(condition)) {
    idx <- which(condition == lv)
    resid[, idx] <- norm[, idx] - rowMeans(norm[, idx, drop = FALSE])
  }
  s2 <- rowSums(resid^2) / df
  
  # Method-of-moments: Var = mu/s + phi*mu^2  =>  phi = (s2 - mu * mean(1/s)) / mu^2
  phi_raw <- pmax((s2 - mu * mean(1 / size_factors)) / mu^2, 1e-8)
  
  keep <- is.finite(phi_raw) & mu > 1
  fit  <- lm(phi_raw[keep] ~ I(1 / mu[keep]))
  a0   <- max(unname(coef(fit)[1]), 1e-4)       # asymptotic dispersion
  a1   <- max(unname(coef(fit)[2]), 0)
  phi_trend <- a0 + a1 / pmax(mu, 1e-3)
  
  w <- df / (df + prior_df)                      # weight on gene-wise estimate
  phi_final <- exp(w * log(phi_raw) + (1 - w) * log(phi_trend))
  phi_final[!is.finite(phi_final)] <- phi_trend[!is.finite(phi_final)]
  
  list(mu = mu, raw = phi_raw, trend = phi_trend, final = phi_final,
       a0 = a0, a1 = a1)
}

# ----------------------------------------------------------------------------
# 4. Variance-stabilising transformation (for EDA ONLY, never for testing)
# ----------------------------------------------------------------------------
# For NB data Var grows with the mean, so highly expressed genes dominate
# distances and PCA. A plain log2(x+1) over-inflates the variance of low counts.
# This is the NB VST of Anders & Huber (2010) in its simplified form with a
# single asymptotic dispersion a0 (DESeq2's parametric VST also includes the
# a1 term of the trend). It behaves like log2 for large counts and like a
# square root for small ones, giving roughly constant variance across the range.
vst_nb <- function(norm_counts, asymp_disp) {
  a <- asymp_disp
  log2((1 + 2 * a * norm_counts + 2 * sqrt(a * norm_counts * (1 + a * norm_counts))) /
         (4 * a))
}

# ----------------------------------------------------------------------------
# 5. Differential expression: per-gene NB GLM + Wald test + BH correction
# ----------------------------------------------------------------------------
# Model:   log(mu_ij) = log(s_j) + b0_i + b1_i * treated_j
# Offset log(s_j) accounts for depth; b1 is the natural-log fold change.
# Testing is done on RAW counts (never on VST/log values): the NB likelihood
# needs the count scale to be valid.
run_de <- function(counts, coldata, size_factors, dispersions,
                   min_mean_norm = 1) {
  X <- model.matrix(~ condition, data = coldata)
  p <- ncol(X)
  off <- log(size_factors)
  n_genes <- nrow(counts)
  
  out <- matrix(NA_real_, n_genes, 3,
                dimnames = list(rownames(counts), c("beta", "se", "stat")))
  
  for (i in seq_len(n_genes)) {
    if (dispersions$mu[i] < min_mean_norm) next      # independent filtering
    fit <- tryCatch(
      suppressWarnings(
        glm.fit(X, counts[i, ], offset = off,
                family = MASS::negative.binomial(theta = 1 / dispersions$final[i]))),
      error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$converged)) next
    cov_unscaled <- chol2inv(fit$qr$qr[seq_len(p), seq_len(p), drop = FALSE])
    se <- sqrt(diag(cov_unscaled))[2]
    out[i, ] <- c(fit$coefficients[2], se, fit$coefficients[2] / se)
  }
  
  res <- data.frame(
    gene       = rownames(counts),
    baseMean   = dispersions$mu,
    dispersion = dispersions$final,
    log2FC     = out[, "beta"] / log(2),
    lfcSE      = out[, "se"] / log(2),
    stat       = out[, "stat"],
    row.names  = rownames(counts)
  )
  res$pvalue <- 2 * pnorm(-abs(res$stat))
  # Benjamini-Hochberg across TESTED genes only (NA are excluded automatically)
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res[order(res$padj), ]
}

classify_genes <- function(res, fdr, lfc) {
  res$significance <- "Not significant"
  ok <- !is.na(res$padj) & res$padj < fdr
  res$significance[ok & res$log2FC >=  lfc] <- "Upregulated"
  res$significance[ok & res$log2FC <= -lfc] <- "Downregulated"
  res$significance <- factor(res$significance,
                             levels = c("Upregulated", "Downregulated", "Not significant"))
  res
}

# ----------------------------------------------------------------------------
# 6. Functional enrichment: over-representation analysis (hypergeometric)
# ----------------------------------------------------------------------------
# P(X >= k) where k = DE genes in the set, drawn from a universe of N tested
# genes containing K set members, with n DE genes in total.
# The universe must be the genes that COULD have been detected (tested genes),
# not the whole genome - a classic source of inflated enrichment.
run_ora <- function(sig_genes, universe, gene_sets) {
  n <- length(sig_genes); N <- length(universe)
  res <- do.call(rbind, lapply(names(gene_sets), function(nm) {
    set <- intersect(gene_sets[[nm]], universe)
    k <- length(intersect(set, sig_genes)); K <- length(set)
    data.frame(term = nm, set_size = K, overlap = k,
               fold_enrichment = (k / n) / (K / N),
               pvalue = phyper(k - 1, K, N - K, n, lower.tail = FALSE))
  }))
  res$padj <- p.adjust(res$pvalue, "BH")
  res[order(res$padj), ]
}

# ----------------------------------------------------------------------------
# 7. Plotting
# ----------------------------------------------------------------------------
theme_portfolio <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top",
          panel.grid.major = element_line(colour = "grey92"))
}
PAL_COND <- c(Control = "#4C72B0", Treated = "#C44E52")
PAL_SIG  <- c(Upregulated = "#C44E52", Downregulated = "#4C72B0",
              `Not significant` = "grey70")

plot_sample_distances <- function(vst, coldata, file) {
  d <- as.matrix(dist(t(vst)))
  pheatmap(d, annotation_col = coldata,
           annotation_colors = list(condition = PAL_COND),
           color = colorRampPalette(c("#08306B", "white"))(100),
           main = "Sample-to-sample Euclidean distances (VST)",
           filename = file, width = 7, height = 6)
}

# PCA on centred (not scaled) VST values of the most variable genes.
# Not scaling keeps the VST's variance-stabilising property intact; scaling would
# give every gene (incl. noisy ones) equal weight and amplify noise.
# Loadings (rotation) tell which genes drive each PC; % variance tells how much
# of the structure each PC captures - PC1 aligned with condition = strong effect.
plot_pca <- function(vst, coldata, n_top, file) {
  top <- order(rowVars(vst), decreasing = TRUE)[seq_len(n_top)]
  pca <- prcomp(t(vst[top, ]), center = TRUE, scale. = FALSE)
  pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  df  <- data.frame(pca$x[, 1:2], coldata, sample = rownames(coldata))
  
  p <- ggplot(df, aes(PC1, PC2, colour = condition, label = sample)) +
    stat_ellipse(level = 0.9, linetype = 2, linewidth = 0.4, show.legend = FALSE) +
    geom_point(size = 4) +
    geom_text_repel(size = 3, show.legend = FALSE) +
    scale_colour_manual(values = PAL_COND) +
    labs(title = "PCA of variance-stabilised expression",
         subtitle = sprintf("Top %d most variable genes", n_top),
         x = sprintf("PC1 (%.1f%% variance)", pve[1]),
         y = sprintf("PC2 (%.1f%% variance)", pve[2]), colour = NULL) +
    theme_portfolio()
  ggsave(file, p, width = 7, height = 5.5, dpi = 300)
  invisible(list(plot = p, pca = pca, pve = pve))
}

plot_volcano <- function(res, fdr, lfc, n_labels, file) {
  d <- res[!is.na(res$padj), ]
  d$neglog10 <- -log10(pmax(d$pvalue, 1e-300))
  lab <- head(d[d$significance != "Not significant", ][
    order(d$padj[d$significance != "Not significant"]), ], n_labels)
  
  p <- ggplot(d, aes(log2FC, neglog10, colour = significance)) +
    geom_point(alpha = 0.6, size = 1.3) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = 2, colour = "grey40") +
    geom_hline(yintercept = -log10(max(d$pvalue[d$padj < fdr], na.rm = TRUE)),
               linetype = 2, colour = "grey40") +
    geom_text_repel(data = lab, aes(label = gene), size = 3,
                    max.overlaps = Inf, box.padding = 0.4,
                    show.legend = FALSE, colour = "black") +
    scale_colour_manual(values = PAL_SIG) +
    labs(title = "Treated vs Control",
         subtitle = sprintf("FDR < %.2f and |log2FC| >= %g  |  Up: %d, Down: %d",
                            fdr, lfc, sum(d$significance == "Upregulated"),
                            sum(d$significance == "Downregulated")),
         x = expression(log[2]~fold~change), y = expression(-log[10]~p),
         colour = NULL) +
    theme_portfolio()
  ggsave(file, p, width = 7.5, height = 6, dpi = 300)
  invisible(p)
}

plot_ora <- function(ora, file) {
  ora$term <- factor(ora$term, levels = rev(ora$term))
  p <- ggplot(ora, aes(fold_enrichment, term, size = overlap, colour = -log10(padj))) +
    geom_point() +
    scale_colour_gradient(low = "grey70", high = "#C44E52") +
    labs(title = "Over-representation analysis", x = "Fold enrichment",
         y = NULL, colour = expression(-log[10]~FDR), size = "Overlap") +
    theme_portfolio()
  ggsave(file, p, width = 7.5, height = 4.5, dpi = 300)
  invisible(p)
}

# ----------------------------------------------------------------------------
# 8. Benchmark against ground truth (possible only because data are simulated)
# ----------------------------------------------------------------------------
benchmark <- function(res, truth, fdr) {
  m <- merge(res, truth, by = "gene")
  called <- !is.na(m$padj) & m$padj < fdr
  tp <- sum(called & m$is_de); fp <- sum(called & !m$is_de)
  fn <- sum(!called & m$is_de)
  data.frame(called = sum(called), TP = tp, FP = fp, FN = fn,
             precision = tp / max(tp + fp, 1), recall = tp / (tp + fn),
             observed_FDR = fp / max(tp + fp, 1), nominal_FDR = fdr,
             lfc_correlation = cor(m$log2FC[called & m$is_de],
                                   m$true_lfc[called & m$is_de]))
}

# ----------------------------------------------------------------------------
# 9. Main
# ----------------------------------------------------------------------------
main <- function(cfg = CONFIG) {
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
  path <- function(f) file.path(cfg$out_dir, f)
  
  message("[1/6] Simulating counts")
  sim <- simulate_rnaseq(cfg$n_genes, cfg$n_per_group, cfg$frac_de, cfg$seed)
  
  message("[2/6] Normalising and estimating dispersions")
  sf   <- estimate_size_factors(sim$counts)
  disp <- estimate_dispersions(sim$counts, sf, sim$coldata$condition)
  norm <- sweep(sim$counts, 2, sf, "/")
  vst  <- vst_nb(norm, disp$a0)
  
  message("[3/6] QC / EDA")
  plot_sample_distances(vst, sim$coldata, path("01_sample_distances.png"))
  plot_pca(vst, sim$coldata, cfg$n_top_var, path("02_pca.png"))
  
  message("[4/6] Differential expression")
  res <- run_de(sim$counts, sim$coldata, sf, disp)
  res <- classify_genes(res, cfg$fdr_cutoff, cfg$lfc_cutoff)
  write.csv(res, path("de_results.csv"), row.names = FALSE)
  
  message("[5/6] Volcano plot and benchmark")
  plot_volcano(res, cfg$fdr_cutoff, cfg$lfc_cutoff, cfg$n_labels,
               path("03_volcano.png"))
  bm <- benchmark(res, sim$truth, cfg$fdr_cutoff)
  print(bm, digits = 3)
  write.csv(bm, path("benchmark.csv"), row.names = FALSE)
  
  message("[6/6] Enrichment")
  tested   <- res$gene[!is.na(res$padj)]
  sig      <- res$gene[res$significance != "Not significant"]
  ora      <- run_ora(sig, tested, sim$gene_sets)
  write.csv(ora, path("ora_results.csv"), row.names = FALSE)
  plot_ora(ora, path("04_ora.png"))
  print(head(ora, 5), digits = 3)
  
  writeLines(capture.output(sessionInfo()), path("sessionInfo.txt"))
  invisible(list(res = res, ora = ora, benchmark = bm))
}

# Runs when called as `Rscript rnaseq_pipeline.R`. With RStudio's Source button
# sys.nframe() is > 0, so call main() manually.
if (sys.nframe() == 0L) main()

# =============================================================================
# REAL-DATA BLOCK (not run): same workflow with the Bioconductor standard stack.
# Use this on the `airway` dataset or your own counts, and annotate with
# org.Hs.eg.db. Kept commented so the script runs with CRAN-only dependencies.
# =============================================================================
# suppressPackageStartupMessages({
#   library(DESeq2); library(airway); library(AnnotationDbi)
#   library(org.Hs.eg.db); library(clusterProfiler)
# })
#
# data(airway)
# dds <- DESeqDataSet(airway, design = ~ cell + dex)
# dds$dex <- relevel(dds$dex, ref = "untrt")
# dds <- dds[rowSums(counts(dds) >= 10) >= 4, ]      # low-count prefilter
# dds <- DESeq(dds)
# resultsNames(dds)                                  # check the coefficient name
#
# # Unshrunken results: keep them, they hold the Wald statistic used for GSEA.
# res_raw <- results(dds, name = "dex_trt_vs_untrt", alpha = 0.05)
# # Shrunken LFCs (apeglm) for plots/ranking; NOTE: they have no valid `stat`.
# res_shr <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm", res = res_raw)
# vsd <- vst(dds, blind = TRUE)                      # for PCA / sample distances
#
# # Annotation: Ensembl -> SYMBOL / ENTREZ
# df <- as.data.frame(res_shr)
# df$stat   <- res_raw$stat
# df$symbol <- mapIds(org.Hs.eg.db, rownames(df), "SYMBOL",   "ENSEMBL", multiVals = "first")
# df$entrez <- mapIds(org.Hs.eg.db, rownames(df), "ENTREZID", "ENSEMBL", multiVals = "first")
#
# # Over-representation, GO Biological Process, universe = tested genes
# tested <- !is.na(df$padj) & !is.na(df$entrez)
# sig    <- tested & df$padj < 0.05 & abs(df$log2FoldChange) >= 1
# ego <- enrichGO(gene = unique(df$entrez[sig]), universe = unique(df$entrez[tested]),
#                 OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
#                 qvalueCutoff = 0.05, readable = TRUE)
# dotplot(ego, showCategory = 15)
#
# # Threshold-free alternative: GSEA on the ranked Wald statistic
# rk <- df[!is.na(df$entrez) & !is.na(df$stat), ]
# rk <- rk[!duplicated(rk$entrez), ]
# ranks <- sort(setNames(rk$stat, rk$entrez), decreasing = TRUE)
# gse <- gseGO(geneList = ranks, OrgDb = org.Hs.eg.db, ont = "BP",
#              keyType = "ENTREZID", pAdjustMethod = "BH", seed = TRUE)

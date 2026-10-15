# ============================================================
# REPRODUCIBLE FIGURE & TABLE GENERATION SCRIPT
# Manuscript: "Human Kidney Precision Cut Slice Cultures Are a
#              Flexible Platform to Study Mechanisms of Glomerular Disease"
#
# INSTRUCTIONS:
#   1. Place this script in the same directory as
#      "GSE345813_GeoMx_Analysis_Complete_Workspace.RData.gz"
#   2. Source this script: source("reproduce_figures_tables.R")
#   3. All manuscript PDF figures and .xlsx supplemental tables
#      will be generated in the same directory.
#
# Required R packages: ggplot2, pheatmap, patchwork, reshape2, limma,
#                      fgsea, openxlsx, ggrepel, msigdbr, GSVA
#
# MANUSCRIPT MAP (Slice-culture-manuscript_10092026.docx):
# Figures 3, 4A/B, 5A/B: culture media and temporal analyses.
# Figure 6A: APOL1 forest plot; 6B: IFN-gamma ssGSEA across treatments.
# Figure 7A/B: IFN/control and baricitinib/DMSO fGSEA.
# Supplemental Figures 2/3: normalization and PCA; 5A: APOL1 expression;
# 6: IFN vs IFN+baricitinib fGSEA; 7: compartment pathway comparison.
# Supplemental Tables 3-10 follow current manuscript numbering.
# Separate source assets: Figures 1, 2, 6C, 7C; Supplemental Figures 1,
# 4, 5B/C; Supplemental Tables 1/2 and Movies 1/2.
# These images, qPCR data, clinical data and primers are not generated from
# the GeoMx workspace. See Manuscript_Output_Checklist.csv after execution.
# Optional environment variables: SLICE_WORKSPACE, SLICE_OUTPUT_DIR.
# ============================================================
cat("=== Reproducing Manuscript Figures & Supplemental Tables ===\n")
cat("Loading workspace...\n")
script_dir <- NULL
if (is.null(script_dir) && requireNamespace("rstudioapi", quietly = TRUE)) {
  tryCatch({
    if (rstudioapi::isAvailable()) {
      path <- rstudioapi::getActiveDocumentContext()$path
      if (nzchar(path)) script_dir <- dirname(path)
    }
  }, error = function(e) NULL)
}
if (is.null(script_dir)) {
  tryCatch({
    frame <- sys.frame(1)
    if (exists("ofile", envir = frame, inherits = FALSE)) {
      script_dir <- dirname(frame$ofile)
    }
  }, error = function(e) NULL)
}
if (is.null(script_dir)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  }
}
if (is.null(script_dir)) {
  script_dir <- getwd()
}
setwd(script_dir)
cat(paste0("Working directory set to: ", getwd(), "\n"))
workspace_file <- Sys.getenv("SLICE_WORKSPACE", file.path(script_dir, "GSE345813_GeoMx_Analysis_Complete_Workspace.RData.gz"))
if (!file.exists(workspace_file)) {
  stop(paste0("Workspace not found: ", workspace_file,
              "\nDownload GSE345813_GeoMx_Analysis_Complete_Workspace.RData.gz from GEO GSE345813 and place it beside this script, or set SLICE_WORKSPACE."))
}
# R load() accepts gzip-compressed files directly; no decompression is needed.
load(workspace_file)
output_dir <- Sys.getenv("SLICE_OUTPUT_DIR", file.path(script_dir, "manuscript_10092026_outputs"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)
set.seed(20261009)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(reshape2)
library(limma)
library(fgsea)
library(ggrepel)
library(openxlsx)
cat("Workspace loaded. Defining helper functions...\n\n")
# ============================================================
# HELPER FUNCTION DEFINITIONS
# ============================================================
plot_pca <- function(mat, meta_df, title = "", color_by = "Batch") {
  pca <- prcomp(t(mat), scale. = TRUE)
  var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
  df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                   Batch = meta_df$Batch, Donor = meta_df$Donor_ID,
                   Treatment = meta_df$Treatment, Timepoint = meta_df$Timepoint,
                   Medium = meta_df$Medium)
  ggplot(df, aes(x = PC1, y = PC2, color = .data[[color_by]])) +
    geom_point(size = 2, alpha = 0.7) +
    labs(x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"), title = title) +
    theme_bw() + theme(plot.title = element_text(size = 10))
}
plot_media_volcano <- function(de_results, title_label, fdr_thresh = 0.05, lfc_thresh = 0.5) {
  de_results$Significant <- "NS"
  de_results$Significant[de_results$adj.P.Val < fdr_thresh & de_results$logFC > lfc_thresh] <- "Up"
  de_results$Significant[de_results$adj.P.Val < fdr_thresh & de_results$logFC < -lfc_thresh] <- "Down"
  n_up <- sum(de_results$Significant == "Up", na.rm = TRUE)
  n_dn <- sum(de_results$Significant == "Down", na.rm = TRUE)
  top_genes <- de_results[de_results$Significant != "NS", ]
  top_labels <- head(top_genes[order(top_genes$adj.P.Val), ], 15)
  ggplot(de_results, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("NS" = "grey70", "Up" = "red3", "Down" = "steelblue")) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey50") +
    geom_text_repel(data = top_labels, aes(label = Gene), size = 2.3, max.overlaps = 15,
                    color = "black", segment.color = "grey40", segment.size = 0.3) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3,
             label = paste0(n_up, " up\n", n_dn, " down"), size = 3, fontface = "bold") +
    ggtitle(title_label) + xlab("log2 Fold Change") + ylab("-log10(P-value)") +
    theme_bw() + theme(plot.title = element_text(size = 9, face = "bold"), legend.position = "none")
}
run_media_3group_DE <- function(expr_mat, meta_df, compartment_label) {
  design <- model.matrix(~ 0 + Condition + Donor_ID, data = meta_df)
  colnames(design) <- make.names(colnames(design))
  fit <- lmFit(expr_mat, design)
  contrast_mat <- makeContrasts(
    A_vs_T0 = ConditionA - ConditionT0,
    BO2_vs_T0 = ConditionB.O2 - ConditionT0,
    BO2_vs_A = ConditionB.O2 - ConditionA,
    levels = design
  )
  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2)
  results_list <- list()
  contrasts <- c("A_vs_T0", "BO2_vs_T0", "BO2_vs_A")
  for (j in seq_along(contrasts)) {
    res <- topTable(fit2, coef = contrasts[j], number = Inf, sort.by = "none")
    res$Gene <- rownames(res)
    res$Compartment <- compartment_label
    results_list[[contrasts[j]]] <- res
  }
  return(results_list)
}
plot_media_delta_barplot <- function(stats_df, ssgsea_mat, meta_df, compartment_label) {
  conditions <- levels(meta_df$Condition)
  mean_scores <- sapply(conditions, function(cond) {
    idx <- which(meta_df$Condition == cond)
    rowMeans(ssgsea_mat[, idx, drop = FALSE])
  })
  colnames(mean_scores) <- conditions
  delta_df <- data.frame(
    Pathway = gsub("HALLMARK_", "", rownames(mean_scores)),
    Delta_A = mean_scores[, "A"] - mean_scores[, "T0"],
    Delta_BO2 = mean_scores[, "B-O2"] - mean_scores[, "T0"],
    stringsAsFactors = FALSE
  )
  delta_df$sig_A <- stats_df$FDR_A_vs_T0[match(delta_df$Pathway, stats_df$Pathway)]
  delta_df$sig_BO2 <- stats_df$FDR_BO2_vs_T0[match(delta_df$Pathway, stats_df$Pathway)]
  make_sig_label <- function(fdr) {
    ifelse(is.na(fdr), "", ifelse(fdr < 0.001, "***", ifelse(fdr < 0.01, "**", ifelse(fdr < 0.05, "*", ""))))
  }
  delta_df$label_A <- make_sig_label(delta_df$sig_A)
  delta_df$label_BO2 <- make_sig_label(delta_df$sig_BO2)
  delta_long <- melt(delta_df[, c("Pathway", "Delta_A", "Delta_BO2")],
                      id.vars = "Pathway", variable.name = "Medium", value.name = "Delta_Score")
  delta_long$Medium <- ifelse(delta_long$Medium == "Delta_A", "A vs T0", "B-O2 vs T0")
  pathway_order <- delta_df$Pathway[order(delta_df$Delta_BO2)]
  delta_long$Pathway <- factor(delta_long$Pathway, levels = pathway_order)
  sig_annot <- rbind(
    data.frame(Pathway = delta_df$Pathway, Medium = "A vs T0",
               Delta_Score = delta_df$Delta_A, label = delta_df$label_A, stringsAsFactors = FALSE),
    data.frame(Pathway = delta_df$Pathway, Medium = "B-O2 vs T0",
               Delta_Score = delta_df$Delta_BO2, label = delta_df$label_BO2, stringsAsFactors = FALSE)
  )
  sig_annot <- sig_annot[sig_annot$label != "", ]
  sig_annot$Pathway <- factor(sig_annot$Pathway, levels = pathway_order)
  sig_annot$x_pos <- ifelse(sig_annot$Delta_Score > 0,
                             sig_annot$Delta_Score + 0.01, sig_annot$Delta_Score - 0.01)
  ggplot(delta_long, aes(x = Pathway, y = Delta_Score, fill = Medium)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
    coord_flip() + geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("A vs T0" = "steelblue", "B-O2 vs T0" = "darkorange")) +
    geom_text(data = sig_annot, aes(x = Pathway, y = x_pos, label = label, group = Medium),
              position = position_dodge(width = 0.8), inherit.aes = FALSE, size = 3, color = "black") +
    ggtitle(paste0(compartment_label, ": Media Effects on Pathways (vs T0, 3 Donors)\n",
                   "* FDR<0.05, **FDR<0.01,** *FDR<0.001")) +
    xlab("") + ylab("Delta ssGSEA Score (vs T0)") +
    theme_bw() + theme(plot.title = element_text(size = 10, face = "bold"),
                       axis.text.y = element_text(size = 6.5), legend.position = "bottom")
}
plot_pseudotime_heatmap_ROIs <- function(expr_mat, meta_df, temporal_obj, compartment_label, filename) {
  genes <- temporal_obj$gene_order
  patterns <- temporal_obj$pattern
  timepoints <- c("T0", "T16", "T48", "T72")
  meta_df$Timepoint <- factor(meta_df$Timepoint, levels = timepoints)
  col_order <- order(meta_df$Timepoint)
  hm_mat <- expr_mat[genes, col_order]
  hm_meta <- meta_df[col_order, ]
  hm_z <- t(scale(t(hm_mat)))
  hm_z[hm_z > 2.5] <- 2.5
  hm_z[hm_z < -2.5] <- -2.5
  annotation_col <- data.frame(Timepoint = hm_meta$Timepoint, row.names = colnames(hm_z))
  annotation_row <- data.frame(Pattern = patterns, row.names = genes)
  time_colors <- c("T0" = "#E41A1C", "T16" = "#377EB8", "T48" = "#4DAF4A", "T72" = "#984EA3")
  pattern_palette <- c("Early_Peak_T0" = "#E41A1C", "Monotonic_Down" = "#FF7F00",
                       "Mid_Peak_T16" = "#377EB8", "Transient_Dip" = "#A6CEE3",
                       "Mid_Peak_T48" = "#4DAF4A", "Peak_T0" = "#FDBF6F",
                       "Peak_T16" = "#80B1D3", "Peak_T48" = "#B2DF8A",
                       "Peak_T72" = "#CAB2D6", "Late_Peak_T72" = "#984EA3",
                       "Monotonic_Up" = "#6A3D9A")
  unique_patterns <- unique(patterns)
  pattern_colors <- pattern_palette[unique_patterns]
  pattern_colors <- pattern_colors[!is.na(pattern_colors)]
  cluster_rle <- rle(as.character(patterns))
  row_gaps <- cumsum(cluster_rle$lengths[-length(cluster_rle$lengths)])
  tp_rle <- rle(as.character(hm_meta$Timepoint))
  col_gaps <- cumsum(tp_rle$lengths[-length(tp_rle$lengths)])
  pdf(filename, width = 10, height = 16)
  pheatmap(hm_z, cluster_rows = FALSE, cluster_cols = FALSE,
           show_colnames = FALSE, show_rownames = TRUE, fontsize_row = 4.5,
           annotation_col = annotation_col, annotation_row = annotation_row,
           gaps_row = as.numeric(row_gaps), gaps_col = as.numeric(col_gaps),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           breaks = seq(-2.5, 2.5, length.out = 101),
           main = paste0(compartment_label, ": Gene Expression Dynamics Over 72 hrs in Culture\n",
                         "(MTS028, B-O2 | Individual ROIs | Row Z-score | n=", length(genes), " genes)"),
           annotation_colors = list(Timepoint = time_colors, Pattern = pattern_colors))
  dev.off()
  cat("  ✓ Saved:", filename, "\n")
}
plot_fgsea_bar <- function(fgsea_results, title_label, padj_thresh = 0.05) {
  fgsea_sig <- fgsea_results[fgsea_results$padj < padj_thresh, ]
  if (nrow(fgsea_sig) == 0) {
    fgsea_sig <- head(fgsea_results, 10)
    title_label <- paste0(title_label, "\n(None FDR<0.05; showing top 10 by p-value)")
  }
  fgsea_sig <- fgsea_sig[order(fgsea_sig$NES), ]
  fgsea_sig$pathway_clean <- factor(fgsea_sig$pathway_clean, levels = fgsea_sig$pathway_clean)
  fgsea_sig$padj_label <- ifelse(fgsea_sig$padj < 1e-10,
                                  formatC(fgsea_sig$padj, format = "e", digits = 1),
                                  formatC(fgsea_sig$padj, format = "e", digits = 2))
  ggplot(fgsea_sig, aes(x = pathway_clean, y = NES, fill = NES > 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "steelblue"), guide = "none") +
    coord_flip() + geom_hline(yintercept = 0) +
    geom_text(aes(label = padj_label, y = ifelse(NES > 0, NES + 0.08, NES - 0.08),
                  hjust = ifelse(NES > 0, 0, 1)), size = 2.5, color = "black") +
    ggtitle(title_label) + xlab("") + ylab("Normalized Enrichment Score (NES)") +
    theme_bw() + theme(plot.title = element_text(size = 10, face = "bold"),
                       axis.text.y = element_text(size = 7))
}
plot_forest <- function(gene, logFC_mat, se_mat, meta_result, compartment = "Glomerulus") {
  lfc <- logFC_mat[gene, ]; se <- se_mat[gene, ]; donors <- names(lfc)
  meta_row <- meta_result[meta_result$Gene == gene, ]
  df_plot <- data.frame(Donor = c(donors, "META"), logFC = c(lfc, meta_row$meta_logFC),
                        SE = c(se, meta_row$meta_SE), stringsAsFactors = FALSE)
  df_plot$lower <- df_plot$logFC - 1.96 * df_plot$SE
  df_plot$upper <- df_plot$logFC + 1.96 * df_plot$SE
  df_plot$is_meta <- c(rep(FALSE, length(donors)), TRUE)
  df_plot$Donor <- factor(df_plot$Donor, levels = rev(df_plot$Donor))
  ggplot(df_plot, aes(x = logFC, y = Donor)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    geom_point(aes(size = is_meta, color = is_meta)) +
    scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red3")) +
    guides(size = "none", color = "none") +
    labs(title = paste0(gene, " (", compartment, ")"), x = "log2 FC (IFN vs Control)", y = "") +
    theme_bw() + theme(plot.title = element_text(size = 11, face = "bold"))
}
cat("Helper functions defined.\n\n")
# ============================================================
# RECOMPUTE: Objects that may not have been saved in workspace
# ============================================================
if (!exists("de_media_3g_glom")) {
  cat("Recomputing media DE results...\n")
  media_gene_full_glom_idx <- which(
    meta_glom$Donor_ID %in% c("MTS023", "MTS024", "MTS028") &
    meta_glom$Treatment == "None" &
    ((meta_glom$Timepoint == "T0") |
     (meta_glom$Timepoint == "T16" & meta_glom$Medium %in% c("A", "B-O2")))
  )
  meta_mg_glom <- meta_glom[media_gene_full_glom_idx, ]
  expr_mg_glom <- qnorm_glom_genes[, media_gene_full_glom_idx]
  meta_mg_glom$Condition <- factor(ifelse(meta_mg_glom$Timepoint == "T0", "T0", meta_mg_glom$Medium),
                                    levels = c("T0", "A", "B-O2"))
  meta_mg_glom$Donor_ID <- factor(meta_mg_glom$Donor_ID)
  de_media_3g_glom <- run_media_3group_DE(expr_mg_glom, meta_mg_glom, "Glomerulus")
  media_gene_full_ti_idx <- which(
    meta_ti$Donor_ID %in% c("MTS023", "MTS024", "MTS028") &
    meta_ti$Treatment == "None" &
    ((meta_ti$Timepoint == "T0") |
     (meta_ti$Timepoint == "T16" & meta_ti$Medium %in% c("A", "B-O2")))
  )
  meta_mg_ti <- meta_ti[media_gene_full_ti_idx, ]
  expr_mg_ti <- qnorm_ti_genes[, media_gene_full_ti_idx]
  meta_mg_ti$Condition <- factor(ifelse(meta_mg_ti$Timepoint == "T0", "T0", meta_mg_ti$Medium),
                                  levels = c("T0", "A", "B-O2"))
  meta_mg_ti$Donor_ID <- factor(meta_mg_ti$Donor_ID)
  de_media_3g_ti <- run_media_3group_DE(expr_mg_ti, meta_mg_ti, "Tubulointerstitium")
  cat("  ✓ Media DE results recomputed\n")
}
if (!exists("media_stats_glom")) {
  cat("Recomputing media ssGSEA stats...\n")
  media_full_glom_idx <- which(
    meta_glom$Donor_ID %in% c("MTS023", "MTS024", "MTS028") &
    ((meta_glom$Timepoint == "T0") |
     (meta_glom$Timepoint == "T16" & meta_glom$Medium %in% c("A", "B-O2") & meta_glom$Treatment == "None"))
  )
  meta_media_full_glom <- meta_glom[media_full_glom_idx, ]
  ssgsea_media_full_glom <- ssgsea_glom[, media_full_glom_idx]
  meta_media_full_glom$Condition <- factor(ifelse(meta_media_full_glom$Timepoint == "T0", "T0",
                                                   meta_media_full_glom$Medium), levels = c("T0", "A", "B-O2"))
  media_full_ti_idx <- which(
    meta_ti$Donor_ID %in% c("MTS023", "MTS024", "MTS028") &
    ((meta_ti$Timepoint == "T0") |
     (meta_ti$Timepoint == "T16" & meta_ti$Medium %in% c("A", "B-O2") & meta_ti$Treatment == "None"))
  )
  meta_media_full_ti <- meta_ti[media_full_ti_idx, ]
  ssgsea_media_full_ti <- ssgsea_ti[, media_full_ti_idx]
  meta_media_full_ti$Condition <- factor(ifelse(meta_media_full_ti$Timepoint == "T0", "T0",
                                                 meta_media_full_ti$Medium), levels = c("T0", "A", "B-O2"))
  run_media_pairwise_tests <- function(ssgsea_mat, meta_df, compartment_label) {
    pathways <- rownames(ssgsea_mat)
    results <- data.frame(Pathway = gsub("HALLMARK_", "", pathways),
                          Delta_A_vs_T0 = NA_real_, Delta_BO2_vs_T0 = NA_real_, Delta_BO2_vs_A = NA_real_,
                          pval_A_vs_T0 = NA_real_, pval_BO2_vs_T0 = NA_real_, pval_BO2_vs_A = NA_real_,
                          stringsAsFactors = FALSE)
    for (i in seq_along(pathways)) {
      pw <- pathways[i]
      df_tmp <- data.frame(Score = ssgsea_mat[pw, ], Condition = meta_df$Condition, Donor = meta_df$Donor_ID)
      fit <- lm(Score ~ Condition + Donor, data = df_tmp)
      coefs <- summary(fit)$coefficients
      if ("ConditionA" %in% rownames(coefs)) {
        results$Delta_A_vs_T0[i] <- coefs["ConditionA", "Estimate"]
        results$pval_A_vs_T0[i] <- coefs["ConditionA", "Pr(>|t|)"]
      }
      if ("ConditionB-O2" %in% rownames(coefs)) {
        results$Delta_BO2_vs_T0[i] <- coefs["ConditionB-O2", "Estimate"]
        results$pval_BO2_vs_T0[i] <- coefs["ConditionB-O2", "Pr(>|t|)"]
      }
      df_tmp$Condition2 <- relevel(df_tmp$Condition, ref = "A")
      fit2 <- lm(Score ~ Condition2 + Donor, data = df_tmp)
      coefs2 <- summary(fit2)$coefficients
      if ("Condition2B-O2" %in% rownames(coefs2)) {
        results$Delta_BO2_vs_A[i] <- coefs2["Condition2B-O2", "Estimate"]
        results$pval_BO2_vs_A[i] <- coefs2["Condition2B-O2", "Pr(>|t|)"]
      }
    }
    results$FDR_A_vs_T0 <- p.adjust(results$pval_A_vs_T0, method = "BH")
    results$FDR_BO2_vs_T0 <- p.adjust(results$pval_BO2_vs_T0, method = "BH")
    results$FDR_BO2_vs_A <- p.adjust(results$pval_BO2_vs_A, method = "BH")
    results$Compartment <- compartment_label
    return(results)
  }
  media_stats_glom <- run_media_pairwise_tests(ssgsea_media_full_glom, meta_media_full_glom, "Glomerulus")
  media_stats_ti <- run_media_pairwise_tests(ssgsea_media_full_ti, meta_media_full_ti, "Tubulointerstitium")
  cat("  ✓ Media ssGSEA stats recomputed\n")
}
if (!exists("pathway_comp_df")) {
  cat("Recomputing pathway comparison df...\n")
  glom_ifn_idx <- which(meta_glom$Donor_ID == "MTS029" & meta_glom$Timepoint == "T16" & meta_glom$Treatment == "IFN")
  glom_dmso_idx <- which(meta_glom$Donor_ID == "MTS029" & meta_glom$Timepoint == "T16" & meta_glom$Treatment == "DMSO")
  ti_ifn_idx <- which(meta_ti$Donor_ID == "MTS029" & meta_ti$Timepoint == "T16" & meta_ti$Treatment == "IFN")
  ti_dmso_idx <- which(meta_ti$Donor_ID == "MTS029" & meta_ti$Timepoint == "T16" & meta_ti$Treatment == "DMSO")
  glom_delta <- rowMeans(ssgsea_glom[, glom_ifn_idx]) - rowMeans(ssgsea_glom[, glom_dmso_idx])
  ti_delta <- rowMeans(ssgsea_ti[, ti_ifn_idx]) - rowMeans(ssgsea_ti[, ti_dmso_idx])
  pathway_comp_df <- data.frame(Pathway = gsub("HALLMARK_", "", names(glom_delta)),
                                 Delta_Glom = glom_delta, Delta_TI = ti_delta[names(glom_delta)],
                                 stringsAsFactors = FALSE)
  cat("  ✓ Pathway comparison df recomputed\n")
}

# Recalculate compartment interaction statistics from the final analysis-log model.
gi <- which(meta_glom$Donor_ID == "MTS029" & meta_glom$Timepoint == "T16" & meta_glom$Treatment %in% c("DMSO", "IFN"))
ti <- which(meta_ti$Donor_ID == "MTS029" & meta_ti$Timepoint == "T16" & meta_ti$Treatment %in% c("DMSO", "IFN"))
interaction_p <- sapply(paste0("HALLMARK_", pathway_comp_df$Pathway), function(pw) {
  scores <- data.frame(Score = c(ssgsea_glom[pw, gi], ssgsea_ti[pw, ti]),
                       Treatment = factor(c(as.character(meta_glom$Treatment[gi]), as.character(meta_ti$Treatment[ti])), levels = c("DMSO", "IFN")),
                       Compartment = factor(c(rep("Glomerulus", length(gi)), rep("Tubulointerstitium", length(ti)))))
  anova(lm(Score ~ Treatment * Compartment, data = scores))["Treatment:Compartment", "Pr(>F)"]
})
pathway_comp_df$Interaction_FDR <- p.adjust(interaction_p, "BH")
pathway_comp_df$Larger_in <- ifelse(abs(pathway_comp_df$Delta_Glom) > abs(pathway_comp_df$Delta_TI), "Glomerulus", "Tubulointerstitium")
pathway_comp_df$sig_label <- ifelse(pathway_comp_df$Interaction_FDR < .001, "***", ifelse(pathway_comp_df$Interaction_FDR < .01, "**", ifelse(pathway_comp_df$Interaction_FDR < .05, "*", "")))

cat("\nGenerating figures...\n\n")
# ============================================================
# SUPPLEMENTAL FIGURE 2: Quantile normalization boxplots
# ============================================================
cat("Generating Supplemental Figure 2...\n")
pdf("SuppFig2_QNorm_Boxplots_Glomerulus.pdf", width = 14, height = 8)
par(mfrow = c(2, 1), mar = c(2, 4, 2, 1))
batch_colors_glom <- ifelse(meta_glom$Batch == "1", "steelblue", "darkorange")
boxplot(log2_glom, las = 2, cex.axis = 0.3, outline = FALSE,
        col = batch_colors_glom, border = batch_colors_glom, xaxt = "n",
        main = "Glomerulus - BEFORE Quantile Normalization (log2)", ylab = "log2(counts + 1)")
legend("topright", legend = c("Batch 1", "Batch 2"), fill = c("steelblue", "darkorange"), cex = 0.8)
boxplot(qnorm_glom, las = 2, cex.axis = 0.3, outline = FALSE,
        col = batch_colors_glom, border = batch_colors_glom, xaxt = "n",
        main = "Glomerulus - AFTER Quantile Normalization", ylab = "Normalized expression")
legend("topright", legend = c("Batch 1", "Batch 2"), fill = c("steelblue", "darkorange"), cex = 0.8)
dev.off()
pdf("SuppFig2_QNorm_Boxplots_Tubulointerstitium.pdf", width = 14, height = 8)
par(mfrow = c(2, 1), mar = c(2, 4, 2, 1))
batch_colors_ti <- ifelse(meta_ti$Batch == "1", "steelblue", "darkorange")
boxplot(log2_ti, las = 2, cex.axis = 0.3, outline = FALSE,
        col = batch_colors_ti, border = batch_colors_ti, xaxt = "n",
        main = "Tubulointerstitium - BEFORE Quantile Normalization (log2)", ylab = "log2(counts + 1)")
legend("topright", legend = c("Batch 1", "Batch 2"), fill = c("steelblue", "darkorange"), cex = 0.8)
boxplot(qnorm_ti, las = 2, cex.axis = 0.3, outline = FALSE,
        col = batch_colors_ti, border = batch_colors_ti, xaxt = "n",
        main = "Tubulointerstitium - AFTER Quantile Normalization", ylab = "Normalized expression")
legend("topright", legend = c("Batch 1", "Batch 2"), fill = c("steelblue", "darkorange"), cex = 0.8)
dev.off()
cat("  ✓ SuppFig2_QNorm_Boxplots_Glomerulus.pdf\n  ✓ SuppFig2_QNorm_Boxplots_Tubulointerstitium.pdf\n")
# ============================================================
# SUPPLEMENTAL FIGURE 3: PCA before/after quantile normalization
# ============================================================
cat("Generating Supplemental Figure 3...\n")
pdf("SuppFig3_PCA_Glomerulus_BeforeAfter_QNorm.pdf", width = 12, height = 10)
p1 <- plot_pca(log2_glom, meta_glom, "BEFORE QNorm - Batch", "Batch")
p2 <- plot_pca(log2_glom, meta_glom, "BEFORE QNorm - Donor", "Donor")
p3 <- plot_pca(log2_glom, meta_glom, "BEFORE QNorm - Treatment", "Treatment")
p4 <- plot_pca(qnorm_glom, meta_glom, "AFTER QNorm - Batch", "Batch")
p5 <- plot_pca(qnorm_glom, meta_glom, "AFTER QNorm - Donor", "Donor")
p6 <- plot_pca(qnorm_glom, meta_glom, "AFTER QNorm - Treatment", "Treatment")
print((p1 | p2 | p3) / (p4 | p5 | p6))
dev.off()
pdf("SuppFig3_PCA_Tubulointerstitium_BeforeAfter_QNorm.pdf", width = 12, height = 10)
p1 <- plot_pca(log2_ti, meta_ti, "BEFORE QNorm - Batch", "Batch")
p2 <- plot_pca(log2_ti, meta_ti, "BEFORE QNorm - Donor", "Donor")
p3 <- plot_pca(log2_ti, meta_ti, "BEFORE QNorm - Treatment", "Treatment")
p4 <- plot_pca(qnorm_ti, meta_ti, "AFTER QNorm - Batch", "Batch")
p5 <- plot_pca(qnorm_ti, meta_ti, "AFTER QNorm - Donor", "Donor")
p6 <- plot_pca(qnorm_ti, meta_ti, "AFTER QNorm - Treatment", "Treatment")
print((p1 | p2 | p3) / (p4 | p5 | p6))
dev.off()
cat("  ✓ SuppFig3_PCA_Glomerulus_BeforeAfter_QNorm.pdf\n  ✓ SuppFig3_PCA_Tubulointerstitium_BeforeAfter_QNorm.pdf\n")
# ============================================================
# FIGURE 3: Media DEG volcano plots (all comparisons)
# ============================================================
cat("Generating Figure 3...\n")
pdf("Fig3_Media_Effects_Volcano_DEGs_AllComparisons.pdf", width = 14, height = 12)
print(
  (plot_media_volcano(de_media_3g_glom$A_vs_T0, "Glomerulus: A vs T0") |
   plot_media_volcano(de_media_3g_glom$BO2_vs_T0, "Glomerulus: B-O2 vs T0") |
   plot_media_volcano(de_media_3g_glom$BO2_vs_A, "Glomerulus: B-O2 vs A")) /
  (plot_media_volcano(de_media_3g_ti$A_vs_T0, "TI: A vs T0") |
   plot_media_volcano(de_media_3g_ti$BO2_vs_T0, "TI: B-O2 vs T0") |
   plot_media_volcano(de_media_3g_ti$BO2_vs_A, "TI: B-O2 vs A"))
)
dev.off()
cat("  ✓ Fig3_Media_Effects_Volcano_DEGs_AllComparisons.pdf\n")
# ============================================================
# FIGURE 4: Media effects ssGSEA barplots (T0 vs A vs B-O2)
# ============================================================
cat("Generating Figure 4...\n")
pdf("Fig4A_Media_T0vsAvsBO2_Pathways_Glomerulus.pdf", width = 10, height = 11)
print(plot_media_delta_barplot(media_stats_glom, ssgsea_media_full_glom, meta_media_full_glom, "Glomerulus"))
dev.off()
pdf("Fig4B_Media_T0vsAvsBO2_Pathways_Tubulointerstitium.pdf", width = 10, height = 11)
print(plot_media_delta_barplot(media_stats_ti, ssgsea_media_full_ti, meta_media_full_ti, "Tubulointerstitium"))
dev.off()
cat("  ✓ Fig4A_Media_T0vsAvsBO2_Pathways_Glomerulus.pdf\n  ✓ Fig4B_Media_T0vsAvsBO2_Pathways_Tubulointerstitium.pdf\n")
# ============================================================
# FIGURE 5: Pseudotime gene heatmaps (72hr time course)
# ============================================================
cat("Generating Figure 5...\n")
plot_pseudotime_heatmap_ROIs(expr_tc_glom, meta_tc_glom, temporal_glom, "Glomerulus",
                              "Fig5A_Pseudotime_GeneHeatmap_IndividualROIs_Glom.pdf")
plot_pseudotime_heatmap_ROIs(expr_tc_ti, meta_tc_ti, temporal_ti, "Tubulointerstitium",
                              "Fig5B_Pseudotime_GeneHeatmap_IndividualROIs_TI.pdf")
# ============================================================
# FIGURE 7A: fGSEA — IFN vs Control (5-donor meta-analysis)
# ============================================================
cat("Generating Figure 7A...\n")
pdf("Fig7A_MetaAnalysis_fGSEA_IFN_vs_Control_5Donors.pdf", width = 9, height = 8)
print(plot_fgsea_bar(fgsea_meta, "fGSEA: IFN vs Control (5-donor meta-analysis, Glomerulus, T16)"))
dev.off()
cat("  ✓ Fig7A_MetaAnalysis_fGSEA_IFN_vs_Control_5Donors.pdf\n")
# ============================================================
# FIGURE 7B: fGSEA — Baricitinib vs DMSO (2-donor meta)
# ============================================================
cat("Generating Figure 7B...\n")
pdf("Fig7B_fGSEA_Bari_vs_DMSO_2Donors_Glom.pdf", width = 9, height = 8)
print(plot_fgsea_bar(fgsea_bari, "fGSEA: Baricitinib vs DMSO (2-donor meta: MTS008 + MTS029, Glom)"))
dev.off()
cat("  ✓ Fig7B_fGSEA_Bari_vs_DMSO_2Donors_Glom.pdf\n")
# ============================================================
# FIGURE 6A: APOL1 forest plot
# ============================================================
cat("Generating Figure 6A...\n")
stopifnot("APOL1" %in% meta_results$Gene)
pdf("Fig6A_APOL1_ForestPlot_Glomerulus.pdf", width = 6, height = 5)
print(plot_forest("APOL1", logFC_mat, se_mat, meta_results))
dev.off()
# ============================================================
# FIGURE 6B: IFN-gamma ssGSEA scores across treatments (Glom)
# ============================================================
cat("Generating Figure 6B...\n")
ifn_pathways_to_plot <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
treat_idx_glom <- which(meta_glom$Donor_ID %in% c("MTS008", "MTS029") &
                          meta_glom$Timepoint == "T16" &
                          meta_glom$Treatment %in% c("DMSO", "IFN", "Baricitinib", "IFN-Baricitinib"))
treat_meta <- meta_glom[treat_idx_glom, ]
treat_ssgsea <- ssgsea_glom[, treat_idx_glom]
treat_plot_list <- list()
for (pw in ifn_pathways_to_plot[ifn_pathways_to_plot %in% rownames(treat_ssgsea)]) {
  treat_plot_list[[pw]] <- data.frame(Score = treat_ssgsea[pw, ], Treatment = treat_meta$Treatment,
                                       Donor = treat_meta$Donor_ID, Pathway = gsub("HALLMARK_", "", pw))
}
treat_plot_df <- do.call(rbind, treat_plot_list)
treat_plot_df$Treatment <- factor(treat_plot_df$Treatment,
                                   levels = c("DMSO", "IFN", "IFN-Baricitinib", "Baricitinib"))
pdf("Fig6B_Baricitinib_ssGSEA_Pathways_Glom.pdf", width = 6, height = 8)
print(ggplot(treat_plot_df, aes(x = Treatment, y = Score, fill = Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  facet_grid(Donor ~ Pathway, scales = "free_y") +
  scale_fill_manual(values = c("DMSO" = "khaki", "IFN" = "red3",
                                "IFN-Baricitinib" = "turquoise", "Baricitinib" = "purple")) +
  ggtitle("Pathway Scores Across Treatment Conditions (Glomerulus, T16)") +
  ylab("ssGSEA Score") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 8, face = "bold"), legend.position = "none"))
dev.off()
cat("  ✓ Fig6B_Baricitinib_ssGSEA_Pathways_Glom.pdf\n")
# ============================================================
# SUPPLEMENTAL FIGURE 6: fGSEA — IFN vs IFN-Baricitinib
# ============================================================
cat("Generating Supplemental Figure 6...\n")
pdf("SuppFig6_fGSEA_IFN_vs_IFNBari_2Donors_Glom.pdf", width = 9, height = 8)
print(plot_fgsea_bar(fgsea_rescue, "fGSEA: IFN vs IFN-Baricitinib (2-donor meta: MTS008 + MTS029, Glom)"))
dev.off()
cat("  ✓ SuppFig6_fGSEA_IFN_vs_IFNBari_2Donors_Glom.pdf\n")
# ============================================================
# SUPPLEMENTAL FIGURE 7: Glom vs TI pathway delta comparison
# ============================================================
cat("Generating Supplemental Figure 7...\n")
pathway_comp_long <- melt(pathway_comp_df[, c("Pathway", "Delta_Glom", "Delta_TI")],
                           id.vars = "Pathway", variable.name = "Compartment", value.name = "Delta_Score")
pathway_comp_long$Compartment <- ifelse(pathway_comp_long$Compartment == "Delta_Glom",
                                         "Glomerulus", "Tubulointerstitium")
pathway_order <- pathway_comp_df$Pathway[order(pathway_comp_df$Delta_Glom)]
pathway_comp_long$Pathway <- factor(pathway_comp_long$Pathway, levels = pathway_order)
pdf("SuppFig7_GlomVsTI_Pathway_Delta_Comparison.pdf", width = 10, height = 11)
print(ggplot(pathway_comp_long, aes(x = Pathway, y = Delta_Score, fill = Compartment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
  coord_flip() + geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Glomerulus" = "steelblue", "Tubulointerstitium" = "darkorange")) +
  ggtitle("IFN Effect on Pathways: Glom vs TI (MTS029, T16)") +
  xlab("") + ylab("Delta ssGSEA Score (IFN - DMSO)") +
  theme_bw() + theme(axis.text.y = element_text(size = 6.5), legend.position = "bottom"))
dev.off()
cat("  ✓ SuppFig7_GlomVsTI_Pathway_Delta_Comparison.pdf\n")
# ============================================================
#          SUPPLEMENTAL TABLES (Excel .xlsx output)
#
# ============================================================
# ============================================================
cat("\n=== Generating Supplemental Tables (.xlsx) ===\n\n")
# ============================================================
# SUPPLEMENTAL TABLE 3: GeoMx ROI Summary Table
# ============================================================
cat("Generating Supplemental Table 3...\n")
wb3 <- createWorkbook()
meta_glom_tmp <- meta_glom
meta_glom_tmp$Compartment <- "Glomerulus"
meta_ti_tmp <- meta_ti
meta_ti_tmp$Compartment <- "Tubulointerstitium"
common_cols <- intersect(colnames(meta_glom_tmp), colnames(meta_ti_tmp))
meta_all <- rbind(meta_glom_tmp[, common_cols], meta_ti_tmp[, common_cols])
summary_df <- as.data.frame(
  table(Donor_ID = meta_all$Donor_ID,
        Compartment = meta_all$Compartment,
        Timepoint = meta_all$Timepoint,
        Treatment = meta_all$Treatment,
        Medium = meta_all$Medium),
  stringsAsFactors = FALSE
)
colnames(summary_df)[ncol(summary_df)] <- "n_ROIs"
summary_df <- summary_df[summary_df$n_ROIs > 0, ]
summary_df <- summary_df[order(summary_df$Donor_ID, summary_df$Compartment,
                                summary_df$Timepoint, summary_df$Treatment), ]
rownames(summary_df) <- NULL
addWorksheet(wb3, "ROI_Summary")
writeData(wb3, "ROI_Summary", summary_df)
setColWidths(wb3, "ROI_Summary", cols = 1:ncol(summary_df), widths = "auto")
summary_simple <- as.data.frame(
  table(Donor_ID = meta_all$Donor_ID, Compartment = meta_all$Compartment),
  stringsAsFactors = FALSE
)
colnames(summary_simple)[3] <- "Total_ROIs"
summary_simple <- summary_simple[summary_simple$Total_ROIs > 0, ]
addWorksheet(wb3, "ROI_Totals_Per_Donor")
writeData(wb3, "ROI_Totals_Per_Donor", summary_simple)
setColWidths(wb3, "ROI_Totals_Per_Donor", cols = 1:ncol(summary_simple), widths = "auto")
saveWorkbook(wb3, "Suppl_Table_03_GeoMx_ROI_Summary.xlsx", overwrite = TRUE)
cat("  ✓ Suppl_Table_03_GeoMx_ROI_Summary.xlsx\n")
# ============================================================
# SUPPLEMENTAL TABLE 4: Media Effects DEG Results (Multi-sheet)
# ============================================================
cat("Generating Supplemental Table 4...\n")
wb5 <- createWorkbook()
for (contrast_name in names(de_media_3g_glom)) {
  sheet_name <- paste0("Glom_", contrast_name)
  addWorksheet(wb5, sheet_name)
  de_out <- de_media_3g_glom[[contrast_name]]
  de_out <- de_out[order(de_out$adj.P.Val), ]
  writeData(wb5, sheet_name, de_out)
  setColWidths(wb5, sheet_name, cols = 1:ncol(de_out), widths = "auto")
}
for (contrast_name in names(de_media_3g_ti)) {
  sheet_name <- paste0("TI_", contrast_name)
  addWorksheet(wb5, sheet_name)
  de_out <- de_media_3g_ti[[contrast_name]]
  de_out <- de_out[order(de_out$adj.P.Val), ]
  writeData(wb5, sheet_name, de_out)
  setColWidths(wb5, sheet_name, cols = 1:ncol(de_out), widths = "auto")
}
saveWorkbook(wb5, "Suppl_Table_04_Media_Effects_DEG_Results.xlsx", overwrite = TRUE)
cat("  ✓ Suppl_Table_04_Media_Effects_DEG_Results.xlsx\n")
# ============================================================
# SUPPLEMENTAL TABLE 5: Media Effects ssGSEA Pathway Statistics
# ============================================================
cat("Generating Supplemental Table 5...\n")
wb6 <- createWorkbook()
addWorksheet(wb6, "Glomerulus")
glom_stats_out <- media_stats_glom[order(media_stats_glom$FDR_BO2_vs_T0), ]
writeData(wb6, "Glomerulus", glom_stats_out)
setColWidths(wb6, "Glomerulus", cols = 1:ncol(glom_stats_out), widths = "auto")
addWorksheet(wb6, "Tubulointerstitium")
ti_stats_out <- media_stats_ti[order(media_stats_ti$FDR_BO2_vs_T0), ]
writeData(wb6, "Tubulointerstitium", ti_stats_out)
setColWidths(wb6, "Tubulointerstitium", cols = 1:ncol(ti_stats_out), widths = "auto")
saveWorkbook(wb6, "Suppl_Table_05_Media_Effects_ssGSEA_Pathway_Stats.xlsx", overwrite = TRUE)
cat("  ✓ Suppl_Table_05_Media_Effects_ssGSEA_Pathway_Stats.xlsx\n")
# ============================================================
# SUPPLEMENTAL TABLE 6: Temporal Gene Clusters (72hr Time Course)
# ============================================================
cat("Generating Supplemental Table 6...\n")
wb7 <- createWorkbook()
addWorksheet(wb7, "Glomerulus")
temporal_glom_df <- data.frame(
  Gene = temporal_glom$gene_order,
  Temporal_Pattern = temporal_glom$pattern,
  Compartment = "Glomerulus",
  stringsAsFactors = FALSE
)
if (exists("expr_tc_glom") && exists("meta_tc_glom")) {
  for (tp in c("T0", "T16", "T48", "T72")) {
    tp_idx <- which(meta_tc_glom$Timepoint == tp)
    if (length(tp_idx) > 0) {
      temporal_glom_df[[paste0("Mean_Expr_", tp)]] <- rowMeans(expr_tc_glom[temporal_glom$gene_order, tp_idx, drop = FALSE])
    }
  }
}
writeData(wb7, "Glomerulus", temporal_glom_df)
setColWidths(wb7, "Glomerulus", cols = 1:ncol(temporal_glom_df), widths = "auto")
addWorksheet(wb7, "Tubulointerstitium")
temporal_ti_df <- data.frame(
  Gene = temporal_ti$gene_order,
  Temporal_Pattern = temporal_ti$pattern,
  Compartment = "Tubulointerstitium",
  stringsAsFactors = FALSE
)
if (exists("expr_tc_ti") && exists("meta_tc_ti")) {
  for (tp in c("T0", "T16", "T48", "T72")) {
    tp_idx <- which(meta_tc_ti$Timepoint == tp)
    if (length(tp_idx) > 0) {
      temporal_ti_df[[paste0("Mean_Expr_", tp)]] <- rowMeans(expr_tc_ti[temporal_ti$gene_order, tp_idx, drop = FALSE])
    }
  }
}
writeData(wb7, "Tubulointerstitium", temporal_ti_df)
setColWidths(wb7, "Tubulointerstitium", cols = 1:ncol(temporal_ti_df), widths = "auto")
addWorksheet(wb7, "Pattern_Summary")
pattern_summary <- rbind(
  data.frame(Compartment = "Glomerulus",
             as.data.frame(table(Pattern = temporal_glom$pattern), stringsAsFactors = FALSE)),
  data.frame(Compartment = "Tubulointerstitium",
             as.data.frame(table(Pattern = temporal_ti$pattern), stringsAsFactors = FALSE))
)
colnames(pattern_summary)[3] <- "Gene_Count"
writeData(wb7, "Pattern_Summary", pattern_summary)
setColWidths(wb7, "Pattern_Summary", cols = 1:ncol(pattern_summary), widths = "auto")
saveWorkbook(wb7, "Suppl_Table_06_Temporal_Gene_Clusters_72hr.xlsx", overwrite = TRUE)
cat("  ✓ Suppl_Table_06_Temporal_Gene_Clusters_72hr.xlsx\n")
# ============================================================
# SUPPLEMENTAL TABLE 7: IFN vs Control Meta-Analysis DEGs (5 Donors)
# ============================================================
cat("Generating Supplemental Table 7...\n")
wb8 <- createWorkbook()
addWorksheet(wb8, "IFN_vs_Control_MetaAnalysis")
meta_results_out <- meta_results[order(meta_results$meta_pval), ]
writeData(wb8, "IFN_vs_Control_MetaAnalysis", meta_results_out)
setColWidths(wb8, "IFN_vs_Control_MetaAnalysis", cols = 1:ncol(meta_results_out), widths = "auto")
addWorksheet(wb8, "Per_Donor_logFC")
logFC_df <- as.data.frame(logFC_mat)
logFC_df <- cbind(Gene = rownames(logFC_df), logFC_df)
writeData(wb8, "Per_Donor_logFC", logFC_df)
setColWidths(wb8, "Per_Donor_logFC", cols = 1:ncol(logFC_df), widths = "auto")
addWorksheet(wb8, "Per_Donor_SE")
se_df <- as.data.frame(se_mat)
se_df <- cbind(Gene = rownames(se_df), se_df)
writeData(wb8, "Per_Donor_SE", se_df)
setColWidths(wb8, "Per_Donor_SE", cols = 1:ncol(se_df), widths = "auto")
saveWorkbook(wb8, "Suppl_Table_07_IFN_vs_Control_MetaAnalysis_DEGs.xlsx", overwrite = TRUE)
cat("  ✓ Suppl_Table_07_IFN_vs_Control_MetaAnalysis_DEGs.xlsx\n")
# ============================================================
# SUPPLEMENTAL TABLE 9: fGSEA Results (Multi-sheet)
# ============================================================
cat("Generating Supplemental Table 9...\n")
wb9 <- createWorkbook()
clean_fgsea_for_export <- function(fgsea_df) {
  out <- fgsea_df
  if ("leadingEdge" %in% colnames(out)) {
    out$leadingEdge <- sapply(out$leadingEdge, function(x) paste(x, collapse = ", "))
  }
  out <- out[order(out$padj), ]
  return(out)
}
addWorksheet(wb9, "IFN_vs_Control_5Donor")
fgsea_meta_out <- clean_fgsea_for_export(fgsea_meta)
writeData(wb9, "IFN_vs_Control_5Donor", fgsea_meta_out)
setColWidths(wb9, "IFN_vs_Control_5Donor", cols = 1:ncol(fgsea_meta_out), widths = "auto")
addWorksheet(wb9, "Bari_vs_DMSO_2Donor")
fgsea_bari_out <- clean_fgsea_for_export(fgsea_bari)
writeData(wb9, "Bari_vs_DMSO_2Donor", fgsea_bari_out)
setColWidths(wb9, "Bari_vs_DMSO_2Donor", cols = 1:ncol(fgsea_bari_out), widths = "auto")
addWorksheet(wb9, "IFN_vs_IFNBari_2Donor")
fgsea_rescue_out <- clean_fgsea_for_export(fgsea_rescue)
writeData(wb9, "IFN_vs_IFNBari_2Donor", fgsea_rescue_out)
setColWidths(wb9, "IFN_vs_IFNBari_2Donor", cols = 1:ncol(fgsea_rescue_out), widths = "auto")
saveWorkbook(wb9, "Suppl_Table_09_fGSEA_Results_All_Comparisons.xlsx", overwrite = TRUE)
cat("  ✓ Suppl_Table_09_fGSEA_Results_All_Comparisons.xlsx\n")
# ============================================================
# SUPPLEMENTAL TABLE 10: Glom vs TI Pathway Comparison (IFN Effect)
# ============================================================
cat("Generating Supplemental Table 10...\n")
wb10 <- createWorkbook()
addWorksheet(wb10, "Glom_vs_TI_Pathway_Delta")
pathway_comp_out <- pathway_comp_df
pathway_comp_out$Abs_Difference <- abs(pathway_comp_out$Delta_Glom - pathway_comp_out$Delta_TI)
pathway_comp_out <- pathway_comp_out[order(-pathway_comp_out$Abs_Difference), ]
writeData(wb10, "Glom_vs_TI_Pathway_Delta", pathway_comp_out)
setColWidths(wb10, "Glom_vs_TI_Pathway_Delta", cols = 1:ncol(pathway_comp_out), widths = "auto")
saveWorkbook(wb10, "Suppl_Table_10_Glom_vs_TI_IFN_Pathway_Comparison.xlsx", overwrite = TRUE)
cat("  ✓ Suppl_Table_10_Glom_vs_TI_IFN_Pathway_Comparison.xlsx\n")
# ============================================================

# Supplemental Figure 5A: spatial APOL1 expression in both baricitinib donors.
stopifnot("APOL1" %in% rownames(qnorm_glom_genes))
apol1_df <- data.frame(Expression = as.numeric(qnorm_glom_genes["APOL1", treat_idx_glom]),
                       Treatment = treat_meta$Treatment, Donor = treat_meta$Donor_ID)
apol1_df$Treatment <- factor(apol1_df$Treatment,
                             levels = c("DMSO", "IFN", "IFN-Baricitinib", "Baricitinib"))
pdf("SuppFig5A_APOL1_Expression_Glomerulus.pdf", width = 8, height = 5)
print(ggplot(apol1_df, aes(Treatment, Expression, fill = Treatment)) +
        geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 1) +
        facet_wrap(~ Donor) + theme_bw() +
        labs(y = "Quantile-normalized log2 expression", x = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none"))
dev.off()

# Supplemental Table 8: baricitinib DEGs (analysis log, final export section).
stopifnot(exists("meta_bari"), exists("logFC_bari"), exists("se_bari"))
wb_bari <- createWorkbook()
bari_out <- meta_bari[order(meta_bari$meta_pval), ]
bari_out$n_positive <- rowSums(logFC_bari[bari_out$Gene, , drop = FALSE] > 0, na.rm = TRUE)
bari_out$n_negative <- rowSums(logFC_bari[bari_out$Gene, , drop = FALSE] < 0, na.rm = TRUE)
bari_out$concordance <- pmax(bari_out$n_positive, bari_out$n_negative) / ncol(logFC_bari)
for (sheet in c("Bari_vs_DMSO_MetaAnalysis", "Per_Donor_logFC", "Per_Donor_SE")) addWorksheet(wb_bari, sheet)
writeData(wb_bari, "Bari_vs_DMSO_MetaAnalysis", bari_out)
writeData(wb_bari, "Per_Donor_logFC", data.frame(Gene = rownames(logFC_bari), logFC_bari, check.names = FALSE))
writeData(wb_bari, "Per_Donor_SE", data.frame(Gene = rownames(se_bari), se_bari, check.names = FALSE))
saveWorkbook(wb_bari, "Suppl_Table_08_Baricitinib_vs_DMSO_MetaAnalysis_DEGs.xlsx", overwrite = TRUE)

# DONE
# ============================================================

external_outputs <- data.frame(
  Item = c("Figure 1", "Figure 2", "Figure 6C", "Figure 7C", "Supplemental Figure 1", "Supplemental Figure 4", "Supplemental Figure 5B", "Supplemental Figure 5C", "Supplemental Table 1", "Supplemental Table 2", "Supplemental Movie 1", "Supplemental Movie 2"),
  Source = c("Workflow and lightsheet images", "LDH assay and PAS histology", "APOL1 RNA-ISH", "Baricitinib JAK-STAT-OxPhos.pdf", "PAS histology", "NPHS1 and SLC6A19 immunohistochemistry", "APOL1 qPCR", "APOL1 immunohistochemistry", "Donor clinical information workbook", "Primer sequences workbook", "Lightsheet z-stack movie", "Lightsheet flyby movie"),
  Status = "Separate source asset required; not generated from GeoMx workspace")
write.csv(external_outputs, "Manuscript_Output_Checklist.csv", row.names = FALSE)

cat("\n=== GeoMx-derived manuscript outputs generated; see checklist for separate source assets ===\n")
cat(paste0("\nOutput directory: ", getwd(), "\n"))
cat("\n--- PDF Figures ---\n")
print(sort(list.files(pattern = "\\.pdf$")))
cat("\n--- Excel Supplemental Tables ---\n")
print(sort(list.files(pattern = "\\.xlsx$")))
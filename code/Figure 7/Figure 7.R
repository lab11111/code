# ==============================================================================
# Figure 7 — Microglial trajectory, SFMBT2-associated state, and CellChat analysis
#
# Figure panels:
#   7A  AD microglial trajectories: HTS vs MTS vs LTS
#   7B  Candidate-gene trajectories: Control vs AD
#   7C  SFMBT2+ vs SFMBT2- CytoTRACE2 distribution
#   7D  SFMBT2 expression vs TWAS activity correlation
#   7E  Global CellChat interaction networks
#   7F  Differential signaling roles
#   7G  CellChat signaling-role scatter
#   7H  CellChat interaction-count heatmap
#   7I  Incoming communication to SFMBT2+ / SFMBT2-
#   7J  Outgoing communication from SFMBT2+ / SFMBT2-
# Publication-ready analysis script
#
# Required input:
#   scRNA_harmony : Seurat object containing at least:
#     - celltype_global
#     - diagnosis
#     - IC_class
#     - Scoring
#     - preKNN_CytoTRACE2_Score
#     - orig.ident
#     - RNA assay
# ===============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Biobase)
  library(monocle)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(gghalves)
  library(cowplot)
  library(scales)
  library(svglite)
  library(qs)
  library(CellChat)
  library(ComplexHeatmap)
  library(pheatmap)
})

set.seed(1234)

output_dir <- "Figure7"
trajectory_dir <- file.path(output_dir, "trajectory")
cellchat_dir <- file.path(output_dir, "cellchat")
dir.create(trajectory_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cellchat_dir, showWarnings = FALSE, recursive = TRUE)

diagnosis_colors <- c(Control = "#80BCC8", AD = "#D88F91")
sfmbt2_colors <- c("SFMBT2+" = "#D88F91", "SFMBT2-" = "#80BCC8")
ic_colors <- c(HTS = "#C95A71", MTS = "#BFBFBF", LTS = "#2AA7A1")

# ==============================================================================
# PREPARATION. Monocle trajectory construction in microglia
# ==============================================================================

micro_obj <- subset(scRNA_harmony, subset = celltype_global == "Micro")

micro_obj <- FindVariableFeatures(
  micro_obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

count_matrix <- GetAssayData(micro_obj, assay = "RNA", layer = "counts")
count_matrix <- as(count_matrix, "sparseMatrix")

phenotype_data <- new("AnnotatedDataFrame", data = micro_obj@meta.data)
feature_data_df <- data.frame(
  gene_short_name = rownames(count_matrix),
  row.names = rownames(count_matrix),
  stringsAsFactors = FALSE
)
feature_data <- new("AnnotatedDataFrame", data = feature_data_df)

cds <- newCellDataSet(
  count_matrix,
  phenoData = phenotype_data,
  featureData = feature_data,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- rownames(
  subset(fData(cds), num_cells_expressed >= 10)
)

Idents(micro_obj) <- "IC_class"
trajectory_de <- FindMarkers(
  micro_obj,
  ident.1 = "HTS",
  ident.2 = "LTS",
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0,
  verbose = FALSE
)

trajectory_de_sig <- trajectory_de %>%
  rownames_to_column("gene") %>%
  filter(p_val < 0.05) %>%
  arrange(p_val)

write.table(
  trajectory_de_sig,
  file = file.path(trajectory_dir, "trajectory_ordering_DEGs.tsv"),
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

ordering_genes <- head(trajectory_de_sig$gene, 500)
ordering_genes <- intersect(ordering_genes, expressed_genes)

cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)

if ("State" %in% colnames(pData(cds)) && 2 %in% unique(pData(cds)$State)) {
  cds <- orderCells(cds, root_state = 2)
}

# ==============================================================================
# PREPARATION. Candidate-gene expression data for trajectory plots
# ==============================================================================

genes_use <- c("LCORL", "ZFP36L2", "SFMBT2")

get_gene_rows <- function(cds, genes) {
  expression_matrix <- exprs(cds)
  feature_data <- fData(cds)

  gene_map <- data.frame(
    gene = genes,
    row_id = NA_character_,
    stringsAsFactors = FALSE
  )

  for (gene in genes) {
    if (gene %in% rownames(expression_matrix)) {
      gene_map$row_id[gene_map$gene == gene] <- gene
    } else if (
      "gene_short_name" %in% colnames(feature_data) &&
      gene %in% feature_data$gene_short_name
    ) {
      matched_row <- rownames(feature_data)[feature_data$gene_short_name == gene][1]
      gene_map$row_id[gene_map$gene == gene] <- matched_row
    }
  }

  if (any(is.na(gene_map$row_id))) {
    stop(
      "Genes not found in cds: ",
      paste(gene_map$gene[is.na(gene_map$row_id)], collapse = ", ")
    )
  }

  gene_map
}

gene_map <- get_gene_rows(cds, genes_use)
trajectory_pd <- pData(cds)
trajectory_mat <- exprs(cds)

stopifnot(
  "Pseudotime" %in% colnames(trajectory_pd),
  "IC_class" %in% colnames(trajectory_pd),
  "diagnosis" %in% colnames(trajectory_pd)
)

expr_df <- as.data.frame(
  t(as.matrix(trajectory_mat[gene_map$row_id, , drop = FALSE]))
)
colnames(expr_df) <- gene_map$gene
expr_df$cell_id <- rownames(expr_df)

trajectory_plot_df <- trajectory_pd %>%
  mutate(cell_id = rownames(trajectory_pd)) %>%
  left_join(expr_df, by = "cell_id") %>%
  pivot_longer(
    cols = all_of(genes_use),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  mutate(
    expression = as.numeric(expression),
    expression_plot = log1p(expression),
    gene = factor(gene, levels = genes_use)
  )

theme_trajectory <- function() {
  theme_bw(base_size = 10.5) +
    theme(
      panel.grid.major = element_line(color = "#E9E9E9", linewidth = 0.35),
      panel.grid.minor = element_line(color = "#F3F3F3", linewidth = 0.25),
      strip.background = element_rect(fill = "#D9D9D9", color = "#D9D9D9"),
      strip.text = element_text(size = 9.5),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.background = element_blank()
    )
}

pseudotime_range <- range(trajectory_plot_df$Pseudotime, na.rm = TRUE)
pseudotime_breaks <- pretty(pseudotime_range, n = 5)

# ==============================================================================
# 7A. Candidate-gene trajectories in AD microglia: HTS vs MTS vs LTS
# ==============================================================================

trajectory_ic_df <- trajectory_plot_df %>%
  filter(
    diagnosis == "AD",
    IC_class %in% c("HTS", "MTS", "LTS")
  ) %>%
  mutate(group = factor(IC_class, levels = c("HTS", "MTS", "LTS")))

p_trajectory_ic <- ggplot(
  trajectory_ic_df,
  aes(x = Pseudotime, y = expression_plot, color = group, fill = group)
) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    span = 0.75,
    linewidth = 1.15,
    alpha = 0.22
  ) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  scale_color_manual(values = ic_colors) +
  scale_fill_manual(values = ic_colors) +
  scale_x_continuous(breaks = pseudotime_breaks, limits = pseudotime_range) +
  labs(x = "Pseudotime", y = "log1p expression") +
  theme_trajectory() +
  theme(legend.position = "right")

ggsave(
  filename = file.path(trajectory_dir, "Fig_7A_AD_ICclass_candidate_gene_trajectory.svg"),
  plot = p_trajectory_ic,
  width = 12,
  height = 3.8,
  bg = "white"
)

# ==============================================================================
# 7B. Candidate-gene trajectories by diagnosis: Control vs AD
# ==============================================================================

trajectory_diagnosis_df <- trajectory_plot_df %>%
  filter(diagnosis %in% c("Control", "AD")) %>%
  mutate(group = factor(diagnosis, levels = c("Control", "AD")))

p_trajectory_diagnosis <- ggplot(
  trajectory_diagnosis_df,
  aes(x = Pseudotime, y = expression_plot, color = group, fill = group)
) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    span = 0.75,
    linewidth = 1.15,
    alpha = 0.22
  ) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  scale_color_manual(values = diagnosis_colors) +
  scale_fill_manual(values = diagnosis_colors) +
  scale_x_continuous(breaks = pseudotime_breaks, limits = pseudotime_range) +
  labs(x = "Pseudotime", y = "log1p expression") +
  theme_trajectory() +
  theme(legend.position = "right")

ggsave(
  filename = file.path(trajectory_dir, "Fig_7B_diagnosis_candidate_gene_trajectory.svg"),
  plot = p_trajectory_diagnosis,
  width = 12,
  height = 3.8,
  bg = "white"
)

# ==============================================================================
# 7C. SFMBT2-associated CytoTRACE2 state
# ==============================================================================

micro_sfmbt2 <- subset(scRNA_harmony, subset = celltype_global == "Micro")
sfmbt2_expression <- FetchData(micro_sfmbt2, vars = "SFMBT2")

micro_sfmbt2$SFMBT2_status <- ifelse(
  sfmbt2_expression$SFMBT2 > 0,
  "SFMBT2+",
  "SFMBT2-"
)

cytotrace_df <- micro_sfmbt2@meta.data %>%
  transmute(
    SFMBT2_status = factor(SFMBT2_status, levels = c("SFMBT2+", "SFMBT2-")),
    CytoTRACE2_score = as.numeric(preKNN_CytoTRACE2_Score)
  ) %>%
  filter(!is.na(SFMBT2_status), !is.na(CytoTRACE2_score)) %>%
  mutate(CytoTRACE2_percentile = percent_rank(CytoTRACE2_score))

cytotrace_summary <- cytotrace_df %>%
  group_by(SFMBT2_status) %>%
  summarise(
    q1 = quantile(CytoTRACE2_percentile, 0.25, na.rm = TRUE),
    q3 = quantile(CytoTRACE2_percentile, 0.75, na.rm = TRUE),
    median = median(CytoTRACE2_percentile, na.rm = TRUE),
    .groups = "drop"
  )

# Wilcoxon rank-sum test: SFMBT2+ vs SFMBT2-
cytotrace_wilcox <- wilcox.test(
  CytoTRACE2_percentile ~ SFMBT2_status,
  data = cytotrace_df,
  exact = FALSE
)

cytotrace_stats <- data.frame(
  comparison = "SFMBT2+ vs SFMBT2-",
  n_SFMBT2_positive = sum(
    cytotrace_df$SFMBT2_status == "SFMBT2+"
  ),
  n_SFMBT2_negative = sum(
    cytotrace_df$SFMBT2_status == "SFMBT2-"
  ),
  median_SFMBT2_positive = median(
    cytotrace_df$CytoTRACE2_percentile[
      cytotrace_df$SFMBT2_status == "SFMBT2+"
    ],
    na.rm = TRUE
  ),
  median_SFMBT2_negative = median(
    cytotrace_df$CytoTRACE2_percentile[
      cytotrace_df$SFMBT2_status == "SFMBT2-"
    ],
    na.rm = TRUE
  ),
  p_value = cytotrace_wilcox$p.value,
  stringsAsFactors = FALSE
)

cytotrace_stats$significance <- dplyr::case_when(
  cytotrace_stats$p_value < 0.0001 ~ "****",
  cytotrace_stats$p_value < 0.001  ~ "***",
  cytotrace_stats$p_value < 0.01   ~ "**",
  cytotrace_stats$p_value < 0.05   ~ "*",
  TRUE ~ "ns"
)

write.csv(
  cytotrace_stats,
  file = file.path(
    output_dir,
    "Fig_7C_SFMBT2_CytoTRACE2_Wilcoxon.csv"
  ),
  row.names = FALSE
)

cytotrace_p_label <- ifelse(
  cytotrace_stats$p_value < 0.001,
  "Wilcoxon p < 0.001",
  paste0(
    "Wilcoxon p = ",
    formatC(
      cytotrace_stats$p_value,
      format = "f",
      digits = 3
    )
  )
)

rect_width <- 0.02
rect_shift <- 0.10
cytotrace_summary <- cytotrace_summary %>%
  mutate(
    x = as.numeric(SFMBT2_status),
    xmin = x - rect_shift - rect_width,
    xmax = x - rect_shift + rect_width,
    xmid = x - rect_shift
  )

p_cytotrace <- ggplot(
  cytotrace_df,
  aes(x = SFMBT2_status, y = CytoTRACE2_percentile, fill = SFMBT2_status)
) +
  gghalves::geom_half_violin(
    side = "r",
    trim = FALSE,
    color = NA,
    alpha = 0.85,
    width = 0.75
  ) +
  geom_rect(
    data = cytotrace_summary,
    aes(xmin = xmin, xmax = xmax, ymin = q1, ymax = q3, fill = SFMBT2_status),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.0
  ) +
  geom_linerange(
    data = cytotrace_summary,
    aes(x = xmid, ymin = 0, ymax = q1),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  geom_linerange(
    data = cytotrace_summary,
    aes(x = xmid, ymin = q3, ymax = 1),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  geom_segment(
    data = cytotrace_summary,
    aes(x = xmin, xend = xmax, y = median, yend = median),
    inherit.aes = FALSE,
    linewidth = 1.0
  ) +
  scale_fill_manual(values = sfmbt2_colors) +
  scale_y_continuous(
    limits = c(0, 1.12),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 1.075,
    label = cytotrace_p_label,
    size = 3.8
  ) +
  labs(
    x = NULL,
    y = "CytoTRACE2 score percentile"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = file.path(output_dir, "Fig_7C_SFMBT2_CytoTRACE2.svg"),
  plot = p_cytotrace,
  width = 4.2,
  height = 3.8,
  bg = "white"
)

# ==============================================================================
# 7D. Correlation between SFMBT2 expression and TWAS activity score
# ==============================================================================

cor_df <- FetchData(micro_sfmbt2, vars = c("SFMBT2", "Scoring")) %>%
  transmute(
    SFMBT2 = as.numeric(SFMBT2),
    Scoring = as.numeric(Scoring)
  ) %>%
  drop_na()

pearson_test <- cor.test(cor_df$SFMBT2, cor_df$Scoring, method = "pearson")
cor_annotation <- sprintf(
  "Pearson r = %.2f, 95%% CI [%.2f, %.2f], p = %.2e, n = %d",
  unname(pearson_test$estimate),
  pearson_test$conf.int[1],
  pearson_test$conf.int[2],
  pearson_test$p.value,
  nrow(cor_df)
)

zero_margin <- theme(plot.margin = margin(0, 0, 0, 0))

p_cor_top <- ggplot(cor_df, aes(x = SFMBT2)) +
  geom_histogram(bins = 40, fill = "#2CA25F", color = "black", linewidth = 0.25) +
  labs(x = NULL, y = "Count") +
  theme_classic(base_size = 9) +
  zero_margin +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), axis.line.x = element_blank()
  )

p_cor_scatter <- ggplot(cor_df, aes(x = SFMBT2, y = Scoring)) +
  geom_point(size = 0.5, alpha = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = FALSE, linewidth = 0.8) +
  labs(x = "SFMBT2 expression", y = "TWAS activity score") +
  theme_classic(base_size = 10) +
  zero_margin

p_cor_right <- ggplot(cor_df, aes(x = Scoring)) +
  geom_histogram(bins = 40, fill = "#FDAE6B", color = "black", linewidth = 0.25) +
  coord_flip() +
  labs(x = NULL, y = "Count") +
  theme_classic(base_size = 9) +
  zero_margin +
  theme(
    axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), axis.line.y = element_blank()
  )

p_cor_blank <- ggplot() + theme_void() + zero_margin
p_cor_bottom <- cowplot::plot_grid(
  p_cor_scatter, p_cor_right,
  nrow = 1, rel_widths = c(1, 0.18), align = "h", axis = "tb"
)
p_cor_toprow <- cowplot::plot_grid(
  p_cor_top, p_cor_blank,
  nrow = 1, rel_widths = c(1, 0.18), align = "h"
)
p_cor_main <- cowplot::plot_grid(
  p_cor_toprow, p_cor_bottom,
  ncol = 1, rel_heights = c(0.18, 1), align = "v"
)
p_correlation <- ggdraw(p_cor_main) +
  draw_label(cor_annotation, x = 0.01, y = 0.99, hjust = 0, vjust = 1, size = 9)

ggsave(
  filename = file.path(output_dir, "Fig_7D_SFMBT2_Scoring_correlation.pdf"),
  plot = p_correlation,
  width = 4.2,
  height = 3.8,
  units = "in"
)

# ==============================================================================
# PREPARATION. SFMBT2+ / SFMBT2- grouping and CellChat analysis in AD
# ==============================================================================

if (!"RNA" %in% Assays(scRNA_harmony)) {
  stop("RNA assay is required for SFMBT2 grouping and CellChat.")
}

DefaultAssay(scRNA_harmony) <- "RNA"

sfmbt2_values <- FetchData(
  scRNA_harmony,
  vars = "SFMBT2"
)

scRNA_harmony$SFMBT2_status <- as.character(
  scRNA_harmony$celltype_global
)

micro_cells <- WhichCells(
  scRNA_harmony,
  expression = celltype_global == "Micro"
)

scRNA_harmony$SFMBT2_status[micro_cells] <- ifelse(
  sfmbt2_values[
    micro_cells,
    "SFMBT2"
  ] > 0,
  "SFMBT2+",
  "SFMBT2-"
)

scRNA_harmony$SFMBT2_status <- factor(
  scRNA_harmony$SFMBT2_status
)

write.csv(
  data.frame(
    cell = colnames(scRNA_harmony),
    celltype_global = scRNA_harmony$celltype_global,
    diagnosis = scRNA_harmony$diagnosis,
    SFMBT2_status = scRNA_harmony$SFMBT2_status,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    output_dir,
    "SFMBT2_status_assignment.csv"
  ),
  row.names = FALSE
)

cellchat_obj <- scRNA_harmony
cellchat_ad <- subset(
  cellchat_obj,
  subset = diagnosis == "AD"
)
cellchat_ad$samples <- cellchat_ad$orig.ident

cellchat <- createCellChat(
  object = cellchat_ad,
  group.by = "SFMBT2_status",
  assay = "RNA"
)

CellChatDB.use <- subsetDB(CellChatDB.human)
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_network <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)

qsave(cellchat, file.path(cellchat_dir, "cellchat_SFMBT2.qs"))
write.csv(
  cellchat_network,
  file = file.path(cellchat_dir, "cellchat_communication_table.csv"),
  row.names = FALSE
)

# ==============================================================================
# 7E. Global CellChat interaction networks
# ==============================================================================

group_size <- as.numeric(table(cellchat@idents))

svglite(file.path(cellchat_dir, "Fig_7E_cellchat_interaction_count.svg"), width = 4, height = 4)
par(mar = c(1, 1, 1, 1))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  vertex.label.cex = 0.001,
  vertex.label.color = NA,
  title.name = NULL
)
dev.off()

svglite(file.path(cellchat_dir, "Fig_7E_cellchat_interaction_weight.svg"), width = 4, height = 4)
par(mar = c(1, 1, 1, 1))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  vertex.label.cex = 0.001,
  vertex.label.color = NA,
  title.name = NULL
)
dev.off()

# ==============================================================================
# 7F. Differential outgoing and incoming signaling roles
# ==============================================================================

positive_group <- "SFMBT2+"
negative_group <- "SFMBT2-"
centrality_list <- cellchat@netP$centr
pathway_names <- names(centrality_list)

extract_pathway_difference <- function(pathway, metric) {
  values <- centrality_list[[pathway]][[metric]]
  if (is.null(values)) return(NULL)
  if (!(positive_group %in% names(values)) || !(negative_group %in% names(values))) return(NULL)

  positive_value <- as.numeric(values[positive_group])
  negative_value <- as.numeric(values[negative_group])
  if (!is.finite(positive_value) || !is.finite(negative_value)) return(NULL)

  data.frame(
    pathway = pathway,
    SFMBT2_positive = positive_value,
    SFMBT2_negative = negative_value,
    difference = positive_value - negative_value,
    stringsAsFactors = FALSE
  )
}

outgoing_df <- do.call(rbind, lapply(pathway_names, extract_pathway_difference, metric = "outdeg"))
incoming_df <- do.call(rbind, lapply(pathway_names, extract_pathway_difference, metric = "indeg"))

# Effect-size ranking used for the displayed top pathways.
top_outgoing <- outgoing_df %>%
  filter(difference > 0) %>%
  arrange(desc(difference)) %>%
  slice_head(n = 18)

top_incoming <- incoming_df %>%
  filter(difference > 0) %>%
  arrange(desc(difference)) %>%
  slice_head(n = 18)

# ------------------------------------------------------------------------------
# 7F statistical significance for SFMBT2+ vs SFMBT2- pathway communication
# ------------------------------------------------------------------------------
# CellChat centrality is an aggregate quantity, so a conventional two-group
# test cannot be applied directly to one centrality value per group. Instead,
# pathway-specific communication probabilities are compared across the
# SFMBT2+ and SFMBT2- side of the network by permutation.

run_cellchat_pathway_permutation <- function(
    cellchat_obj,
    pathway_names,
    direction = c("outgoing", "incoming"),
    n_perm = 1000,
    seed = 1234
) {

  direction <- match.arg(direction)
  set.seed(seed)

  prob_array <- cellchat_obj@netP$prob
  dim_names <- dimnames(prob_array)

  if (is.null(dim_names)) {
    return(data.frame())
  }

  group_levels <- dim_names[[2]]

  if (!(positive_group %in% group_levels) || !(negative_group %in% group_levels)) {
    return(data.frame())
  }

  other_groups <- setdiff(
    group_levels,
    c(positive_group, negative_group)
  )

  if (length(other_groups) == 0) {
    return(data.frame())
  }

  results <- lapply(
    pathway_names,
    function(pathway) {

      pathway_idx <- match(pathway, dim_names[[3]])
      if (is.na(pathway_idx)) return(NULL)

      if (direction == "outgoing") {
        pos_values <- as.numeric(prob_array[positive_group, other_groups, pathway_idx])
        neg_values <- as.numeric(prob_array[negative_group, other_groups, pathway_idx])
      } else {
        pos_values <- as.numeric(prob_array[other_groups, positive_group, pathway_idx])
        neg_values <- as.numeric(prob_array[other_groups, negative_group, pathway_idx])
      }

      pos_values <- pos_values[is.finite(pos_values)]
      neg_values <- neg_values[is.finite(neg_values)]

      if (length(pos_values) == 0 || length(neg_values) == 0) return(NULL)

      observed_difference <- mean(pos_values) - mean(neg_values)
      combined_values <- c(pos_values, neg_values)
      n_pos <- length(pos_values)

      perm_diffs <- replicate(
        n_perm,
        {
          shuffled <- sample(combined_values, length(combined_values), replace = FALSE)
          mean(shuffled[seq_len(n_pos)]) - mean(shuffled[(n_pos + 1):length(shuffled)])
        }
      )

      p_value <- (1 + sum(abs(perm_diffs) >= abs(observed_difference), na.rm = TRUE)) / (n_perm + 1)

      data.frame(
        pathway = pathway,
        direction = direction,
        SFMBT2_positive_mean = mean(pos_values),
        SFMBT2_negative_mean = mean(neg_values),
        difference = observed_difference,
        p_value = p_value,
        stringsAsFactors = FALSE
      )
    }
  )

  bind_rows(results) %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      significance = case_when(
        p_adj_BH < 0.0001 ~ "****",
        p_adj_BH < 0.001  ~ "***",
        p_adj_BH < 0.01   ~ "**",
        p_adj_BH < 0.05   ~ "*",
        TRUE ~ "ns"
      )
    )
}

outgoing_stats_7F <- run_cellchat_pathway_permutation(
  cellchat, pathway_names, direction = "outgoing", n_perm = 1000, seed = 1234
)

incoming_stats_7F <- run_cellchat_pathway_permutation(
  cellchat, pathway_names, direction = "incoming", n_perm = 1000, seed = 1234
)

write.csv(
  outgoing_stats_7F,
  file.path(cellchat_dir, "Fig_7F_outgoing_SFMBT2plus_vs_minus_significance.csv"),
  row.names = FALSE
)

write.csv(
  incoming_stats_7F,
  file.path(cellchat_dir, "Fig_7F_incoming_SFMBT2plus_vs_minus_significance.csv"),
  row.names = FALSE
)

write.csv(
  top_outgoing,
  file.path(cellchat_dir, "Fig_7F_SFMBT2plus_enriched_outgoing_pathways.csv"),
  row.names = FALSE
)
write.csv(
  top_incoming,
  file.path(cellchat_dir, "Fig_7F_SFMBT2plus_enriched_incoming_pathways.csv"),
  row.names = FALSE
)

selected_pathways <- c(
  "SPP1", "PTPRM", "GAS", "TGFb", "BMP", "CSF", "APP", "ADGRL",
  "VISTA", "PSAP", "CD45", "CNTN", "Netrin", "TENASCIN", "COLLAGEN", "LAMININ"
)

available_pathways <- intersect(selected_pathways, names(cellchat@netP$prob))

if (length(available_pathways) > 0) {
  ht_out <- netAnalysis_signalingRole_heatmap(
    cellchat,
    pattern = "outgoing",
    signaling = available_pathways
  )

  ht_in <- netAnalysis_signalingRole_heatmap(
    cellchat,
    pattern = "incoming",
    signaling = available_pathways
  )

  pdf(file.path(cellchat_dir, "Fig_7F_signaling_role_heatmaps.pdf"), width = 10, height = 7)
  draw(ht_out + ht_in)
  dev.off()
}

# ==============================================================================
# 7G. CellChat signaling-role scatter plot
# ==============================================================================

p_cellchat_scatter <- netAnalysis_signalingRole_scatter(cellchat)
ggsave(
  filename = file.path(cellchat_dir, "Fig_7G_cellchat_signaling_role_scatter.svg"),
  plot = p_cellchat_scatter,
  width = 12,
  height = 6
)

# ==============================================================================
# 7H. CellChat interaction-count heatmap
# ==============================================================================

pdf(file.path(cellchat_dir, "Fig_7H_cellchat_interaction_count_heatmap.pdf"), width = 8, height = 7)
pheatmap::pheatmap(
  cellchat@net$count,
  border_color = "black",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontsize = 10,
  display_numbers = TRUE,
  number_color = "black",
  number_format = "%.0f"
)
dev.off()


# ==============================================================================
# 7I. Incoming communication to SFMBT2+ and SFMBT2- microglia
#     Other cell populations -> SFMBT2+ / SFMBT2-
# ==============================================================================

target_groups_7I <- c(
  "SFMBT2+",
  "SFMBT2-"
)

all_identity_levels <- levels(cellchat@idents)

other_identity_levels <- setdiff(
  all_identity_levels,
  target_groups_7I
)

# Display SFMBT2+ before SFMBT2-.
cellchat@idents <- factor(
  cellchat@idents,
  levels = c(
    target_groups_7I,
    other_identity_levels
  )
)

source_groups_7I <- setdiff(
  levels(cellchat@idents),
  target_groups_7I
)

incoming_to_sfmbt2 <- subsetCommunication(
  cellchat,
  sources.use = source_groups_7I,
  targets.use = target_groups_7I
)

incoming_to_sfmbt2_sig <- incoming_to_sfmbt2 %>%
  filter(pval < 0.05) %>%
  arrange(desc(prob))

incoming_top50 <- head(
  incoming_to_sfmbt2_sig,
  50
)

pair_use_7I <- data.frame(
  interaction_name = unique(
    incoming_top50$interaction_name
  )
)

write.csv(
  incoming_top50,
  file.path(
    cellchat_dir,
    "Fig_7I_top50_incoming_to_SFMBT2.csv"
  ),
  row.names = FALSE
)

svg(
  file.path(
    cellchat_dir,
    "Fig_7I_incoming_to_SFMBT2_bubble.svg"
  ),
  width = 10,
  height = 8
)

netVisual_bubble(
  cellchat,
  sources.use = source_groups_7I,
  targets.use = target_groups_7I,
  pairLR.use = pair_use_7I,
  remove.isolate = FALSE
)

dev.off()


# ==============================================================================
# 7J. Outgoing communication from SFMBT2+ and SFMBT2- microglia
#     SFMBT2+ / SFMBT2- -> Other cell populations
# ==============================================================================

source_groups_7J <- c(
  "SFMBT2+",
  "SFMBT2-"
)

target_groups_7J <- setdiff(
  levels(cellchat@idents),
  source_groups_7J
)

outgoing_from_sfmbt2 <- subsetCommunication(
  cellchat,
  sources.use = source_groups_7J,
  targets.use = target_groups_7J
)

outgoing_from_sfmbt2_sig <- outgoing_from_sfmbt2 %>%
  filter(pval < 0.05) %>%
  arrange(desc(prob))

outgoing_top50 <- head(
  outgoing_from_sfmbt2_sig,
  50
)

pair_use_7J <- data.frame(
  interaction_name = unique(
    outgoing_top50$interaction_name
  )
)

write.csv(
  outgoing_top50,
  file.path(
    cellchat_dir,
    "Fig_7J_top50_outgoing_from_SFMBT2.csv"
  ),
  row.names = FALSE
)

svg(
  file.path(
    cellchat_dir,
    "Fig_7J_outgoing_from_SFMBT2_bubble.svg"
  ),
  width = 10,
  height = 8
)

netVisual_bubble(
  cellchat,
  sources.use = source_groups_7J,
  targets.use = target_groups_7J,
  pairLR.use = pair_use_7J,
  remove.isolate = FALSE
)

dev.off()


# ------------------------------------------------------------------------------
# Reproducibility information
# ------------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

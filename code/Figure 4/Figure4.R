# ============================================================================== 
# Figure 4 | Cell-cell communication, pathway activity, and metabolism analysis
# Groups: HTS, MTS, and LTS
# ============================================================================== 

# ---- 0. Environment -----------------------------------------------------------
# No local directory is included in this publication-ready script.
# Required objects/files should be loaded or placed in the current R session.

suppressPackageStartupMessages({
  library(qs)
  library(CellChat)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(GSVA)
  library(GSEABase)
  library(limma)
  library(scMetabolism)
})

# Group labels used throughout the analysis.
HTS <- "HTS"
MTS <- "MTS"
LTS <- "LTS"
GROUP_COL <- "IC_class"

# Set to TRUE only when a precomputed CellChat object generated with
# HTS/MTS/LTS is available in the current working directory.
USE_PRECOMPUTED_CELLCHAT <- FALSE

# The Seurat object should already be loaded as `scRNA_harmony`.
if (!exists("scRNA_harmony")) {
  stop("Please load the Seurat object as `scRNA_harmony` before running this script.")
}

# NOTE:
# This script assumes that scRNA_harmony@meta.data$IC_class already contains
# the labels HTS, MTS, and LTS.


# ============================================================================== 
# 1. CellChat preprocessing
# ============================================================================== 

metadata <- scRNA_harmony@meta.data
is_empty <- is.na(metadata[[GROUP_COL]]) | metadata[[GROUP_COL]] == ""
metadata[[GROUP_COL]][is_empty] <- metadata$celltype_global[is_empty]
scRNA_harmony@meta.data <- metadata

required_groups <- c(HTS, MTS, LTS)
missing_groups <- setdiff(required_groups, unique(metadata[[GROUP_COL]]))
if (length(missing_groups) > 0) {
  warning(
    "The following expected labels are absent from ", GROUP_COL, ": ",
    paste(missing_groups, collapse = ", ")
  )
}

scRNA_AD <- subset(scRNA_harmony, subset = diagnosis == "AD")
scRNA_AD$samples <- scRNA_AD$orig.ident

# Microglial subset reused in Figures 4G and 4H.
micro_obj <- subset(scRNA_harmony, subset = celltype_global == "Micro")

if (USE_PRECOMPUTED_CELLCHAT) {
  cellchat <- qread("cellchat.qs")
} else {
  cellchat <- createCellChat(
    object = scRNA_AD,
    group.by = GROUP_COL,
    assay = "RNA"
  )

  CellChatDB.use <- subsetDB(CellChatDB.human)
  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  message(
    "Signaling matrix: ",
    nrow(cellchat@data.signaling), " genes x ",
    ncol(cellchat@data.signaling), " cells"
  )

  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- smoothData(cellchat, adj = PPI.human)
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    raw.use = FALSE
  )
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat <- aggregateNet(cellchat)

  qsave(cellchat, "cellchat.qs")
}

# Export all inferred communications.
df_net <- subsetCommunication(cellchat)
save(df_net, file = "df_net.RData")
write.csv(
  df_net,
  "df_net.csv",
  row.names = FALSE
)


# ============================================================================== 
# Figure 4A | Global CellChat interaction network
# ============================================================================== 

group_size <- as.numeric(table(cellchat@idents))

grDevices::svg(
  "Figure4A_CellChat_global_network.svg",
  width = 12,
  height = 6
)
par(mfrow = c(1, 2), xpd = TRUE)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights/strength"
)

dev.off()


# ============================================================================== 
# Figure 4B | Differential signaling roles between HTS and LTS
# ============================================================================== 

HTS_NAME <- HTS
LTS_NAME <- LTS

centrality <- cellchat@netP$centr
pathways <- names(centrality)

extract_pathway_difference <- function(pathway, metric) {
  values <- centrality[[pathway]][[metric]]

  if (is.null(values)) {
    return(NULL)
  }

  if (!(HTS_NAME %in% names(values)) || !(LTS_NAME %in% names(values))) {
    return(NULL)
  }

  hts_value <- as.numeric(values[HTS_NAME])
  lts_value <- as.numeric(values[LTS_NAME])

  if (!is.finite(hts_value) || !is.finite(lts_value)) {
    return(NULL)
  }

  data.frame(
    pathway = pathway,
    HTS = hts_value,
    LTS = lts_value,
    difference = hts_value - lts_value,
    stringsAsFactors = FALSE
  )
}

outgoing_df <- do.call(
  rbind,
  lapply(pathways, extract_pathway_difference, metric = "outdeg")
)

incoming_df <- do.call(
  rbind,
  lapply(pathways, extract_pathway_difference, metric = "indeg")
)

# Retain pathways with stronger signaling in HTS than in LTS.
outgoing_positive <- outgoing_df[outgoing_df$difference > 0, , drop = FALSE]
incoming_positive <- incoming_df[incoming_df$difference > 0, , drop = FALSE]

# The original analysis retained the top 18 pathways.
top_outgoing <- head(
  outgoing_positive[order(outgoing_positive$difference, decreasing = TRUE), ],
  18
)

top_incoming <- head(
  incoming_positive[order(incoming_positive$difference, decreasing = TRUE), ],
  18
)

write.csv(
  top_outgoing,
  "HTS_gt_LTS_top18_outgoing_pathways.csv",
  row.names = FALSE
)

write.csv(
  top_incoming,
  "HTS_gt_LTS_top18_incoming_pathways.csv",
  row.names = FALSE
)

selected_pathways <- c(
  "SPP1", "TGFβ", "CSF", "GAS", "COMPLEMENT", "PTPRM", "CD45",
  "ADGRL", "GALECTIN", "APP", "TENASCIN", "CNTN", "LAMININ",
  "COLLAGEN", "PSAP"
)

heatmap_outgoing <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing",
  signaling = selected_pathways
)

heatmap_incoming <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "incoming",
  signaling = selected_pathways
)

grDevices::svg(
  "Figure4B_signaling_role_heatmaps.svg",
  width = 12,
  height = 7
)
ComplexHeatmap::draw(heatmap_outgoing + heatmap_incoming)
dev.off()


# ============================================================================== 
# Figure 4C | Dominant signaling senders and receivers
# ============================================================================== 

p4c <- netAnalysis_signalingRole_scatter(cellchat)

grDevices::svg(
  "Figure4C_signaling_role_scatter.svg",
  width = 7,
  height = 6
)
print(p4c)
dev.off()


# ============================================================================== 
# Figure 4D-E | Ligand-receptor pair fold changes: HTS versus LTS
# ============================================================================== 

communication_df <- subsetCommunication(cellchat)
communication_df$source <- as.character(communication_df$source)
communication_df$target <- as.character(communication_df$target)

# Construct a unique ligand-receptor pair identifier.
if (all(c("ligand", "receptor") %in% colnames(communication_df))) {
  communication_df$LR_pair <- paste(
    communication_df$ligand,
    communication_df$receptor,
    sep = "_"
  )
} else if ("interaction_name_2" %in% colnames(communication_df)) {
  communication_df$LR_pair <- communication_df$interaction_name_2
} else if ("interaction_name" %in% colnames(communication_df)) {
  communication_df$LR_pair <- communication_df$interaction_name
} else {
  communication_df$LR_pair <- paste(
    communication_df$source,
    communication_df$target,
    communication_df$pathway_name,
    sep = "_"
  )
}

all_celltypes <- sort(
  unique(c(communication_df$source, communication_df$target))
)
other_celltypes <- setdiff(all_celltypes, c(HTS, MTS, LTS))

epsilon <- 1e-9

make_pair_fold_long <- function(direction = c("incoming", "outgoing")) {
  direction <- match.arg(direction)

  if (direction == "incoming") {
    # Other cell types -> HTS/LTS; HTS/LTS act as receivers.
    directional_df <- communication_df %>%
      filter(
        source %in% other_celltypes,
        target %in% c(HTS, LTS)
      ) %>%
      mutate(
        celltype = source,
        group = target
      )
  } else {
    # HTS/LTS -> other cell types; HTS/LTS act as senders.
    directional_df <- communication_df %>%
      filter(
        source %in% c(HTS, LTS),
        target %in% other_celltypes
      ) %>%
      mutate(
        celltype = target,
        group = source
      )
  }

  pair_counts <- directional_df %>%
    group_by(celltype, group) %>%
    summarise(
      n_pairs = n_distinct(LR_pair),
      .groups = "drop"
    )

  complete_grid <- expand.grid(
    celltype = other_celltypes,
    group = c(LTS, HTS),
    stringsAsFactors = FALSE
  )

  pair_counts <- complete_grid %>%
    left_join(pair_counts, by = c("celltype", "group"))
  pair_counts$n_pairs[is.na(pair_counts$n_pairs)] <- 0

  lts_df <- pair_counts %>%
    filter(group == LTS) %>%
    select(celltype, LTS_pairs = n_pairs)

  hts_df <- pair_counts %>%
    filter(group == HTS) %>%
    select(celltype, HTS_pairs = n_pairs)

  fold_df <- lts_df %>%
    left_join(hts_df, by = "celltype")

  fold_df <- fold_df[
    fold_df$LTS_pairs > 0 | fold_df$HTS_pairs > 0,
    ,
    drop = FALSE
  ]

  # LTS is set to 1; HTS is represented as HTS_pairs / LTS_pairs.
  fold_df$fold_change <- fold_df$HTS_pairs / (fold_df$LTS_pairs + epsilon)
  fold_df <- fold_df[
    order(fold_df$fold_change, decreasing = TRUE),
    ,
    drop = FALSE
  ]

  long_df <- rbind(
    data.frame(
      celltype = fold_df$celltype,
      group = LTS,
      value = -1,
      n_pairs = fold_df$LTS_pairs
    ),
    data.frame(
      celltype = fold_df$celltype,
      group = HTS,
      value = fold_df$fold_change,
      n_pairs = fold_df$HTS_pairs
    )
  )

  long_df$direction <- direction
  long_df
}

figure4d_df <- make_pair_fold_long("incoming")
figure4e_df <- make_pair_fold_long("outgoing")

common_levels <- unique(as.character(figure4d_df$celltype))
extra_levels <- setdiff(
  unique(as.character(figure4e_df$celltype)),
  common_levels
)
common_levels <- c(common_levels, extra_levels)

figure4d_df$celltype <- factor(
  as.character(figure4d_df$celltype),
  levels = common_levels
)
figure4e_df$celltype <- factor(
  as.character(figure4e_df$celltype),
  levels = common_levels
)

fill_map <- c(LTS = "#4DBBD5", HTS = "#E64B35")

plot_diverging_lr <- function(data, panel_label, subtitle_text) {
  ggplot(data, aes(x = value, y = celltype, fill = group)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = fill_map) +
    scale_x_continuous(labels = function(x) abs(x)) +
    scale_y_discrete(position = "right") +
    labs(
      x = "Fold change in ligand-receptor pair number (HTS/LTS; LTS = 1)",
      y = NULL,
      title = panel_label,
      subtitle = subtitle_text
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.margin = margin(5.5, 30, 5.5, 5.5),
      axis.text.y.right = element_text(size = 10)
    )
}

p4d <- plot_diverging_lr(
  figure4d_df,
  "D",
  "Cell types interacting with HTS/LTS as receptor providers [incoming]"
)

p4e <- plot_diverging_lr(
  figure4e_df,
  "E",
  "Cell types interacting with HTS/LTS as ligand providers [outgoing]"
)

print(p4d)
print(p4e)

ggsave(
  filename = "Figure4D_LR_pair_fold_incoming.svg",
  plot = p4d,
  device = grDevices::svg,
  width = 7.5,
  height = 5.0
)

ggsave(
  filename = "Figure4E_LR_pair_fold_outgoing.svg",
  plot = p4e,
  device = grDevices::svg,
  width = 7.5,
  height = 5.0
)


# ============================================================================== 
# Figure 4F | Top ligand-receptor interactions targeting HTS/MTS/LTS
# ============================================================================== 

target_groups <- c(HTS, MTS, LTS)
source_groups <- setdiff(levels(cellchat@idents), target_groups)

bubble_df <- subsetCommunication(
  cellchat,
  sources.use = source_groups,
  targets.use = target_groups
)

bubble_df_sig <- subset(bubble_df, pval < 0.05)
bubble_df_top <- bubble_df_sig[order(bubble_df_sig$prob, decreasing = TRUE), ]
bubble_df_top <- head(bubble_df_top, 50)

pair_use <- data.frame(
  interaction_name = unique(bubble_df_top$interaction_name)
)

grDevices::svg(
  "Figure4F_top50_LR_bubble.svg",
  width = 10,
  height = 8
)
netVisual_bubble(
  cellchat,
  sources.use = source_groups,
  targets.use = target_groups,
  pairLR.use = pair_use,
  remove.isolate = FALSE
)
dev.off()


# ============================================================================== 
# Figure 4G | Hallmark GSVA: HTS versus LTS
# ============================================================================== 

hallmark_sets <- getGmt("h.all.v2026.1.Hs.symbols.gmt")

selected_cells <- names(micro_obj[[GROUP_COL, drop = TRUE]])[
  micro_obj[[GROUP_COL, drop = TRUE]] %in% c(HTS, LTS)
]

expr_matrix <- GetAssayData(micro_obj, layer = "counts")[, selected_cells]
message(
  sprintf(
    "Expression matrix: %d genes x %d cells",
    nrow(expr_matrix),
    ncol(expr_matrix)
  )
)

sample_groups <- data.frame(
  Sample = selected_cells,
  Group = micro_obj[[GROUP_COL, drop = TRUE]][selected_cells],
  stringsAsFactors = FALSE
)

message("Group counts:")
print(table(sample_groups$Group))

expression_range <- range(expr_matrix, na.rm = TRUE)
message(
  sprintf(
    "Expression range: %.2f to %.2f",
    expression_range[1],
    expression_range[2]
  )
)

# Log2-transform when the matrix appears to contain raw counts.
if (max(expr_matrix, na.rm = TRUE) > 50) {
  expr_matrix <- log2(expr_matrix + 1)
  transformed_range <- range(expr_matrix, na.rm = TRUE)
  message(
    sprintf(
      "After log2 transformation: %.2f to %.2f",
      transformed_range[1],
      transformed_range[2]
    )
  )
}

# Retain genes detected in at least 5% of selected cells.
keep_genes <- rowSums(expr_matrix > 0) >= ncol(expr_matrix) * 0.05
expr_matrix <- expr_matrix[keep_genes, , drop = FALSE]
message(sprintf("Genes retained after filtering: %d", nrow(expr_matrix)))

gsva_parameters <- gsvaParam(
  exprData = as.matrix(expr_matrix),
  geneSets = hallmark_sets,
  kcdf = "Gaussian",
  minSize = 10,
  maxSize = 500
)

gsva_results <- gsva(gsva_parameters)
message(sprintf("GSVA completed for %d Hallmark pathways", nrow(gsva_results)))

design <- model.matrix(~0 + Group, data = sample_groups)
contrast_matrix <- makeContrasts(
  HTSvsLTS = GroupHTS - GroupLTS,
  levels = design
)

fit <- lmFit(gsva_results, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit)

diff_pathways <- topTable(fit, n = Inf, adjust.method = "BH")

# Retain the original nominal P-value summary for compatibility with the
# previous analysis; the figure below uses FDR-adjusted P < 0.05.
sig_pathways_nominal <- diff_pathways[diff_pathways$P.Value < 0.05, , drop = FALSE]
message(sprintf("Nominally significant pathways: %d", nrow(sig_pathways_nominal)))
message(sprintf("Positive logFC pathways: %d", sum(sig_pathways_nominal$logFC > 0)))
message(sprintf("Negative logFC pathways: %d", sum(sig_pathways_nominal$logFC < 0)))

diff_pathways$Pathway <- rownames(diff_pathways)
diff_pathways$Pathway <- gsub("HALLMARK_", "", diff_pathways$Pathway)

plot_df <- diff_pathways %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(Regulation = ifelse(t > 0, "Up", "Down"))

up_df <- plot_df %>%
  filter(t > 0) %>%
  arrange(desc(t))

down_df <- plot_df %>%
  filter(t < 0) %>%
  arrange(desc(t))

plot_df_ordered <- bind_rows(up_df, down_df)
plot_df_ordered$Pathway <- factor(
  plot_df_ordered$Pathway,
  levels = rev(plot_df_ordered$Pathway)
)

p4g <- ggplot(
  plot_df_ordered,
  aes(x = t, y = Pathway, fill = Regulation)
) +
  geom_col(width = 0.8) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  scale_fill_manual(
    values = c(
      Up = "#E64B35",
      Down = "#4DBBD5"
    )
  ) +
  labs(
    x = "t value of GSVA score (HTS vs LTS)",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11)
  )

n_pathways <- nrow(plot_df_ordered)
print(p4g)

ggsave(
  filename = "Figure4G_GSVA_HTS_vs_LTS.svg",
  plot = p4g,
  device = grDevices::svg,
  width = 5,
  height = max(6, n_pathways * 0.25)
)

write.csv(
  diff_pathways,
  "Figure4G_GSVA_HTS_vs_LTS_all_pathways.csv",
  row.names = FALSE
)


# ============================================================================== 
# Figure 4H | Microglial metabolic pathway activity
# ============================================================================== 

Idents(micro_obj) <- GROUP_COL

# Convert the RNA assay for compatibility with scMetabolism when required.
micro_obj[["RNA"]] <- as(object = micro_obj[["RNA"]], Class = "Assay")

metabolism_obj <- sc.metabolism.Seurat(
  obj = micro_obj,
  method = "AUCell",
  imputation = FALSE,
  metabolism.type = "REACTOME"
)

metabolism_scores <- metabolism_obj@assays$METABOLISM$score
metabolism_obj@meta.data <- cbind(
  metabolism_obj@meta.data,
  t(metabolism_scores)
)

all_metabolic_pathways <- rownames(metabolism_obj@assays$METABOLISM$score)

excluded_pathways <- c(
  "Vitamin D calciferol metabolism",
  "Vitamin C ascorbate metabolism",
  "Vitamin B2 riboflavin metabolism",
  "Vitamin B1 thiamin metabolism",
  "Pp2a mediated dephosphorylation of key metabolic factors",
  "PKA mediated phosphorylation of key metabolic factors",
  "Phenylalanine and tyrosine metabolism",
  "Phenylalanine metabolism",
  "Metabolism of steroid hormones",
  "Metabolism of ingested semet sec mesec into H2SE",
  "Metabolism of angiotensinogen to angiotensins",
  "Metabolism of amine derived hormones",
  "Ketone body metabolism",
  "Hyaluronan metabolism",
  "Diseases associated with surfactant metabolism",
  "Defective csf2rb causes pulmonary surfactant metabolism dysfunction 5 smdp5 ",
  "Class C3 metabotropic glutamate pheromone receptors ",
  "Chrebp activates metabolic gene expression",
  "Abacavir transport and metabolism",
  "Abacavir metabolism",
  "Creatine metabolism",
  "Fructose metabolism"
)

metabolic_pathways_to_plot <- setdiff(
  all_metabolic_pathways,
  excluded_pathways
)

p4h <- DotPlot.metabolism(
  obj = metabolism_obj,
  pathway = metabolic_pathways_to_plot,
  phenotype = GROUP_COL,
  norm = "y"
)

print(p4h)

ggsave(
  filename = "Figure4H_metabolism_dotplot.svg",
  plot = p4h,
  device = grDevices::svg,
  width = 10,
  height = 16
)

# ============================================================================== 
# End of script
# ============================================================================== 

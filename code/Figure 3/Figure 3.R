# ==============================================================================
# Figure 3 | Microglia single-cell analysis
# Publication-ready, panel-separated script
#
# Required object in the current R environment:
#   scRNA_harmony
#
# Required metadata columns:
#   celltype_global, Scoring, diagnosis, orig.ident,
#   CytoTRACE2_Score
#
# Figure panels:
#   Figure 3A | Scoring density on UMAP
#   Figure 3B | Scoring distribution and HTS/MTS/LTS thresholds
#   Figure 3C | HTS/MTS/LTS UMAP
#   Figure 3D | Observed/expected (Ro/e) enrichment
#   Figure 3E | Milo differential abundance analysis
#   Figure 3F | CytoTRACE2 percentile across HTS/MTS/LTS
#   Figure 3G | Correlation between Scoring and CytoTRACE2 percentile
#
# All figure outputs are saved as SVG files in the current working directory.
# ==============================================================================


# ==============================================================================
# 0. Packages
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(Nebulosa)
  library(miloR)
  library(scater)
  library(scran)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gghalves)
  library(cowplot)
  library(patchwork)
  library(grid)
})


# ==============================================================================
# Figure 3A | Scoring density on UMAP
# ==============================================================================

# ------------------------------
# A.1 Prepare microglia
# ------------------------------

micro_4a <- subset(
  scRNA_harmony,
  subset = celltype_global == "Micro"
)


# ------------------------------
# A.2 Plot Scoring density
# ------------------------------

p3a <- Nebulosa::plot_density(
  micro_4a,
  features = "Scoring",
  reduction = "umap"
) +
  scale_color_viridis_c(
    option = "viridis",
    name = "Density"
  ) +
  labs(
    title = "Microglia: Scoring activity"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0,
      face = "bold",
      size = 14
    ),
    axis.title = element_text(size = 12)
  )

print(p3a)


# ------------------------------
# A.3 Save Figure 3A
# ------------------------------

ggsave(
  filename = "Figure3A_Microglia_Scoring_density_UMAP.svg",
  plot = p3a,
  width = 6,
  height = 5,
  units = "in",
  device = grDevices::svg,
  bg = "white"
)


# ==============================================================================
# Figure 3B | Scoring distribution and HTS/MTS/LTS thresholds
# ==============================================================================

# ------------------------------
# B.1 Prepare microglia and Scoring values
# ------------------------------

micro_4b <- subset(
  scRNA_harmony,
  subset = celltype_global == "Micro"
)

score_df_4b <- micro_4b@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell = cell,
    Score = as.numeric(Scoring)
  ) %>%
  filter(!is.na(Score))


# ------------------------------
# B.2 Define HTS/MTS/LTS using quartile thresholds
# ------------------------------

score_thresholds_4b <- quantile(
  score_df_4b$Score,
  probs = c(0.25, 0.75),
  na.rm = TRUE
)

low_thr_4b  <- unname(score_thresholds_4b[1])
high_thr_4b <- unname(score_thresholds_4b[2])

score_df_4b <- score_df_4b %>%
  mutate(
    activity_group = case_when(
      Score < low_thr_4b  ~ "LTS",
      Score > high_thr_4b ~ "HTS",
      TRUE                ~ "MTS"
    ),
    activity_group = factor(
      activity_group,
      levels = c("LTS", "MTS", "HTS")
    )
  )

# Write the 3B-derived classification back to the main Seurat object.
scRNA_harmony$IC_class <- NA_character_
micro_cells_4b <- WhichCells(
  scRNA_harmony,
  expression = celltype_global == "Micro"
)

micro_class_map_4b <- setNames(
  as.character(score_df_4b$activity_group),
  score_df_4b$cell
)

scRNA_harmony$IC_class[micro_cells_4b] <- unname(
  micro_class_map_4b[micro_cells_4b]
)

scRNA_harmony$IC_class <- factor(
  scRNA_harmony$IC_class,
  levels = c("HTS", "MTS", "LTS")
)

group_counts_4b <- table(score_df_4b$activity_group)

x_range_4b <- quantile(
  score_df_4b$Score,
  probs = c(0.005, 0.995),
  na.rm = TRUE
)


# ------------------------------
# B.3 Plot Scoring distribution
# ------------------------------

p3b <- ggplot(
  score_df_4b,
  aes(x = Score)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 90,
    fill = "#D3D3D3",
    color = "white",
    alpha = 0.8
  ) +
  geom_density(
    color = "#FF00FF",
    linetype = "dashed",
    linewidth = 1.0,
    adjust = 1.5
  ) +
  geom_vline(
    xintercept = c(low_thr_4b, high_thr_4b),
    color = "#8B0000",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = (x_range_4b[1] + low_thr_4b) / 2,
    y = Inf,
    label = paste0("LTS\nn = ", group_counts_4b["LTS"]),
    vjust = 2,
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = (low_thr_4b + high_thr_4b) / 2,
    y = Inf,
    label = paste0("MTS\nn = ", group_counts_4b["MTS"]),
    vjust = 2,
    size = 4,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = (high_thr_4b + x_range_4b[2]) / 2,
    y = Inf,
    label = paste0("HTS\nn = ", group_counts_4b["HTS"]),
    vjust = 2,
    size = 4,
    fontface = "bold"
  ) +
  coord_cartesian(
    xlim = x_range_4b
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "TWAS activity distribution in microglia",
    x = "Scoring",
    y = "Density"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    axis.line = element_line(linewidth = 0.6),
    axis.text = element_text(color = "black")
  )

print(p3b)


# ------------------------------
# B.4 Save Figure 3B
# ------------------------------

ggsave(
  filename = "Figure3B_Microglia_Scoring_distribution.svg",
  plot = p3b,
  width = 8,
  height = 5,
  units = "in",
  device = grDevices::svg,
  bg = "white"
)

# Save the updated Seurat object with the 3B-derived IC_class.
saveRDS(
  scRNA_harmony,
  file = "scRNA_harmony_with_IC_class.rds"
)

# Save the thresholds used to define LTS/MTS/HTS.
write.csv(
  data.frame(
    threshold = c("LTS_MTS", "MTS_HTS"),
    value = c(low_thr_4b, high_thr_4b)
  ),
  file = "Figure3B_HTS_MTS_LTS_thresholds.csv",
  row.names = FALSE
)


# ==============================================================================
# Figure 3C | HTS/MTS/LTS UMAP
# ==============================================================================

# ------------------------------
# C.1 Prepare microglia and group labels
# ------------------------------

micro_4c <- subset(
  scRNA_harmony,
  subset = celltype_global == "Micro"
)

micro_4c$IC_class <- factor(
  micro_4c$IC_class,
  levels = c("HTS", "MTS", "LTS")
)

ic_colors_4c <- c(
  HTS = "#D88F91",
  MTS = "#BDBDBD",
  LTS = "#80BCC8"
)


# ------------------------------
# C.2 Plot HTS/MTS/LTS UMAP
# ------------------------------

p3c <- DimPlot(
  micro_4c,
  reduction = "umap",
  group.by = "IC_class",
  cols = ic_colors_4c,
  pt.size = 0.6,
  order = c("LTS", "MTS", "HTS")
) +
  theme_void() +
  labs(title = "Microglia") +
  theme(
    legend.position = "right",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    plot.margin = margin(20, 20, 20, 20)
  )


# ------------------------------
# C.3 Add compact UMAP axis arrows
# ------------------------------

umap_coordinates_4c <- Embeddings(
  micro_4c,
  reduction = "umap"
)

xmin_4c <- min(umap_coordinates_4c[, 1])
xmax_4c <- max(umap_coordinates_4c[, 1])
ymin_4c <- min(umap_coordinates_4c[, 2])
ymax_4c <- max(umap_coordinates_4c[, 2])

x_axis_length_4c <- 0.10 * (xmax_4c - xmin_4c)
y_axis_length_4c <- 0.10 * (ymax_4c - ymin_4c)

p3c <- p3c +
  annotate(
    "segment",
    x = xmin_4c,
    xend = xmin_4c + x_axis_length_4c,
    y = ymin_4c,
    yend = ymin_4c,
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.7
  ) +
  annotate(
    "segment",
    x = xmin_4c,
    xend = xmin_4c,
    y = ymin_4c,
    yend = ymin_4c + y_axis_length_4c,
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = xmin_4c + x_axis_length_4c / 2,
    y = ymin_4c - 0.04 * (ymax_4c - ymin_4c),
    label = "UMAP_1",
    size = 3
  ) +
  annotate(
    "text",
    x = xmin_4c - 0.04 * (xmax_4c - xmin_4c),
    y = ymin_4c + y_axis_length_4c / 2,
    label = "UMAP_2",
    size = 3,
    angle = 90
  )

print(p3c)


# ------------------------------
# C.4 Save Figure 3C
# ------------------------------

ggsave(
  filename = "Figure3C_Microglia_HTS_MTS_LTS_UMAP.svg",
  plot = p3c,
  width = 6,
  height = 5,
  units = "in",
  device = grDevices::svg,
  bg = "white"
)


# ==============================================================================
# Figure 3D | Observed/expected (Ro/e) enrichment
# ==============================================================================

# ------------------------------
# D.1 Prepare microglia metadata
# ------------------------------

micro_4d <- subset(
  scRNA_harmony,
  subset = celltype_global == "Micro"
)

roe_meta_4d <- micro_4d@meta.data %>%
  filter(
    IC_class %in% c("HTS", "MTS", "LTS"),
    !is.na(diagnosis)
  )


# ------------------------------
# D.2 Calculate observed/expected enrichment
# ------------------------------

observed_4d <- table(
  roe_meta_4d$IC_class,
  roe_meta_4d$diagnosis
)

expected_4d <- outer(
  rowSums(observed_4d),
  colSums(observed_4d)
) / sum(observed_4d)

roe_matrix_4d <- observed_4d / expected_4d

roe_df_4d <- as.data.frame(
  as.table(roe_matrix_4d)
) %>%
  setNames(c("IC_class", "diagnosis", "Ro_e")) %>%
  mutate(
    diagnosis = factor(
      diagnosis,
      levels = c("Control", "AD")
    ),
    IC_class = factor(
      IC_class,
      levels = c("LTS", "MTS", "HTS")
    ),
    score_label = case_when(
      Ro_e > 1.0                 ~ "+++",
      Ro_e > 0.8 & Ro_e <= 1.0  ~ "++",
      Ro_e >= 0.2 & Ro_e <= 0.8 ~ "+",
      Ro_e > 0.0 & Ro_e < 0.2   ~ "+/-",
      Ro_e == 0                  ~ "-",
      TRUE                       ~ ""
    )
  )


# ------------------------------
# D.3 Plot Ro/e line graph
# ------------------------------

diagnosis_colors_4d <- c(
  Control = "#4DBBD5",
  AD = "#E64B35"
)

p3d_line <- ggplot(
  roe_df_4d,
  aes(
    x = factor(IC_class, levels = c("HTS", "MTS", "LTS")),
    y = Ro_e,
    color = diagnosis,
    group = diagnosis
  )
) +
  geom_hline(
    yintercept = 1,
    linewidth = 0.5,
    linetype = "dotted",
    color = "grey40"
  ) +
  geom_line(
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_point(size = 4) +
  scale_color_manual(values = diagnosis_colors_4d) +
  labs(
    x = NULL,
    y = "Ro/e index",
    color = "Diagnosis"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black")
  )


# ------------------------------
# D.4 Plot Ro/e heatmap
# ------------------------------

p3d_heatmap <- ggplot(
  roe_df_4d,
  aes(
    x = diagnosis,
    y = IC_class,
    fill = Ro_e
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = score_label),
    size = 4
  ) +
  scale_fill_gradient(
    low = "#FFE5CC",
    high = "#F08A4B"
  ) +
  scale_y_discrete(position = "right") +
  labs(
    x = NULL,
    y = NULL,
    fill = "Ro/e"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    panel.grid = element_blank()
  )


# ------------------------------
# D.5 Combine and save Figure 3D
# ------------------------------

p3d <- p3d_line + p3d_heatmap +
  plot_layout(widths = c(2.2, 1)) +
  plot_annotation(
    title = "IC-class enrichment in control and AD microglia"
  )

print(p3d)

ggsave(
  filename = "Figure3D_Microglia_RoE_enrichment.svg",
  plot = p3d,
  width = 8,
  height = 5,
  units = "in",
  device = grDevices::svg,
  bg = "white"
)


# ==============================================================================
# Figure 3E | Milo differential abundance analysis
# ==============================================================================

# ------------------------------
# E.1 Prepare microglia and SingleCellExperiment object
# ------------------------------

set.seed(1234)

micro_4e <- subset(
  scRNA_harmony,
  subset = celltype_global == "Micro"
)

micro_sce_4e <- as.SingleCellExperiment(micro_4e)

micro_sce_4e <- runPCA(
  micro_sce_4e,
  ncomponents = 30
)

milo_4e <- Milo(micro_sce_4e)

if ("UMAP" %in% reducedDimNames(micro_sce_4e)) {
  reducedDim(milo_4e, "UMAP") <- reducedDim(micro_sce_4e, "UMAP")
}


# ------------------------------
# E.2 Build Milo graph and neighborhoods
# ------------------------------

milo_4e <- buildGraph(
  milo_4e,
  k = 30,
  d = 15,
  reduced.dim = "PCA"
)

milo_4e <- makeNhoods(
  milo_4e,
  prop = 0.20,
  k = 30,
  d = 15,
  refined = TRUE,
  reduced_dims = "PCA"
)

milo_4e <- countCells(
  milo_4e,
  meta.data = as.data.frame(colData(milo_4e)),
  sample = "orig.ident"
)


# ------------------------------
# E.3 Prepare design matrix and test differential abundance
# ------------------------------

milo_design_4e <- as.data.frame(
  colData(milo_4e)
) %>%
  select(orig.ident, diagnosis) %>%
  distinct()

rownames(milo_design_4e) <- milo_design_4e$orig.ident

milo_design_4e <- milo_design_4e[
  colnames(nhoodCounts(milo_4e)),
  ,
  drop = FALSE
]

milo_design_4e$diagnosis <- factor(
  milo_design_4e$diagnosis,
  levels = c("Control", "AD")
)

milo_4e <- calcNhoodDistance(
  milo_4e,
  d = 15,
  reduced.dim = "PCA"
)

da_results_4e <- testNhoods(
  milo_4e,
  design = ~ diagnosis,
  design.df = milo_design_4e
)

milo_4e <- buildNhoodGraph(milo_4e)

da_results_4e <- annotateNhoods(
  milo_4e,
  da_results_4e,
  coldata_col = "IC_class"
)


# ------------------------------
# E.4 Plot Milo differential abundance on UMAP
# ------------------------------

p3e_umap <- plotReducedDim(
  milo_4e,
  dimred = "UMAP",
  colour_by = "diagnosis",
  text_by = "IC_class",
  text_size = 3,
  point_size = 0.5
)

p3e_graph <- plotNhoodGraphDA(
  milo_4e,
  da_results_4e,
  layout = "UMAP",
  alpha = 0.90
)

p3e <- p3e_umap + p3e_graph +
  plot_layout(guides = "collect")

print(p3e)


# ------------------------------
# E.5 Plot Milo differential abundance by HTS/MTS/LTS
# ------------------------------

p3e_beeswarm <- plotDAbeeswarm(
  da_results_4e,
  group.by = "IC_class",
  alpha = 0.90
)

print(p3e_beeswarm)


# ------------------------------
# E.6 Save Figure 3E and Milo results
# ------------------------------

ggsave(
  filename = "Figure3E_Milo_DA_UMAP.svg",
  plot = p3e,
  width = 10,
  height = 5,
  units = "in",
  device = grDevices::svg,
  bg = "white"
)

ggsave(
  filename = "Figure3E_Milo_DA_beeswarm.svg",
  plot = p3e_beeswarm,
  width = 8,
  height = 6,
  units = "in",
  device = grDevices::svg,
  bg = "white"
)

write.csv(
  da_results_4e,
  file = "Figure3E_Milo_differential_abundance_results.csv",
  row.names = FALSE
)


# ==============================================================================
# Figure 3F | CytoTRACE2 percentile across HTS/MTS/LTS
# ==============================================================================

# ------------------------------
# F.1 Prepare CytoTRACE2 data
# ------------------------------

micro_4f <- subset(
  scRNA_harmony,
  subset = celltype_global == "Micro"
)

cytotrace_df_4f <- micro_4f@meta.data %>%
  transmute(
    IC_class,
    CytoTRACE2_Score = as.numeric(CytoTRACE2_Score)
  ) %>%
  filter(
    IC_class %in% c("HTS", "MTS", "LTS"),
    !is.na(CytoTRACE2_Score)
  ) %>%
  mutate(
    IC_class = factor(
      IC_class,
      levels = c("HTS", "MTS", "LTS")
    ),
    CytoTRACE2_percentile = percent_rank(CytoTRACE2_Score)
  )

# Pairwise Wilcoxon rank-sum tests with Benjamini-Hochberg correction.
wilcox_comparisons_4f <- list(
  c("HTS", "MTS"),
  c("HTS", "LTS"),
  c("MTS", "LTS")
)

wilcox_results_4f <- bind_rows(
  lapply(
    wilcox_comparisons_4f,
    function(pair) {

      test_df <- cytotrace_df_4f %>%
        filter(
          IC_class %in% pair
        ) %>%
        droplevels()

      test_result <- wilcox.test(
        CytoTRACE2_percentile ~ IC_class,
        data = test_df,
        exact = FALSE
      )

      data.frame(
        group1 = pair[1],
        group2 = pair[2],
        p_value = test_result$p.value
      )
    }
  )
) %>%
  mutate(
    p_adj = p.adjust(
      p_value,
      method = "BH"
    ),
    significance = case_when(
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001  ~ "***",
      p_adj < 0.01   ~ "**",
      p_adj < 0.05   ~ "*",
      TRUE           ~ "ns"
    )
  )

write.csv(
  wilcox_results_4f,
  "Figure3F_CytoTRACE2_Wilcoxon_results.csv",
  row.names = FALSE
)

print(wilcox_results_4f)


# ------------------------------
# F.2 Calculate distribution summary
# ------------------------------

cytotrace_summary_4f <- cytotrace_df_4f %>%
  group_by(IC_class) %>%
  summarise(
    q1 = quantile(
      CytoTRACE2_percentile,
      0.25,
      na.rm = TRUE
    ),
    median = median(
      CytoTRACE2_percentile,
      na.rm = TRUE
    ),
    q3 = quantile(
      CytoTRACE2_percentile,
      0.75,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

rect_width_4f <- 0.02
rect_shift_4f <- 0.10

cytotrace_summary_4f <- cytotrace_summary_4f %>%
  mutate(
    x = as.numeric(IC_class),
    xmin = x - rect_shift_4f - rect_width_4f,
    xmax = x - rect_shift_4f + rect_width_4f,
    xmid = x - rect_shift_4f
  )

ic_colors_4f <- c(
  HTS = "#D88F91",
  MTS = "#BDBDBD",
  LTS = "#80BCC8"
)


# ------------------------------
# F.3 Plot CytoTRACE2 percentile distribution
# ------------------------------

sig_positions_4f <- data.frame(
  group1 = c("HTS", "HTS", "MTS"),
  group2 = c("MTS", "LTS", "LTS"),
  y = c(1.00, 1.04, 1.08)
) %>%
  left_join(
    wilcox_results_4f,
    by = c("group1", "group2")
  ) %>%
  mutate(
    x1 = match(group1, c("HTS", "MTS", "LTS")),
    x2 = match(group2, c("HTS", "MTS", "LTS")),
    xm = (x1 + x2) / 2
  )

p3f <- ggplot(
  cytotrace_df_4f,
  aes(
    x = IC_class,
    y = CytoTRACE2_percentile,
    fill = IC_class
  )
) +
  gghalves::geom_half_violin(
    side = "r",
    trim = FALSE,
    color = NA,
    alpha = 0.85,
    width = 0.75
  ) +
  geom_rect(
    data = cytotrace_summary_4f,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = q1,
      ymax = q3,
      fill = IC_class
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.0
  ) +
  geom_linerange(
    data = cytotrace_summary_4f,
    aes(
      x = xmid,
      ymin = 0,
      ymax = q1
    ),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  geom_linerange(
    data = cytotrace_summary_4f,
    aes(
      x = xmid,
      ymin = q3,
      ymax = 1
    ),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  geom_segment(
    data = cytotrace_summary_4f,
    aes(
      x = xmin,
      xend = xmax,
      y = median,
      yend = median
    ),
    inherit.aes = FALSE,
    linewidth = 1.0
  ) +
  scale_fill_manual(values = ic_colors_4f) +
  scale_y_continuous(
    limits = c(0, 1.12),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "CytoTRACE2 score percentile"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  )

p3f <- p3f +
  geom_segment(
    data = sig_positions_4f,
    aes(
      x = x1,
      xend = x2,
      y = y,
      yend = y
    ),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_segment(
    data = sig_positions_4f,
    aes(
      x = x1,
      xend = x1,
      y = y - 0.01,
      yend = y
    ),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_segment(
    data = sig_positions_4f,
    aes(
      x = x2,
      xend = x2,
      y = y - 0.01,
      yend = y
    ),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_text(
    data = sig_positions_4f,
    aes(
      x = xm,
      y = y + 0.012,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  )

print(p3f)


# ------------------------------
# F.4 Save Figure 3F
# ------------------------------

ggsave(
  filename = "Figure3F_CytoTRACE2_by_HTS_MTS_LTS.svg",
  plot = p3f,
  width = 6,
  height = 4,
  units = "in",
  device = grDevices::svg,
  bg = "white"
)


# ==============================================================================
# Figure 3G | Scoring-CytoTRACE2 correlation
# ==============================================================================

# ------------------------------
# G.1 Prepare correlation data
# ------------------------------

micro_4g <- subset(
  scRNA_harmony,
  subset = celltype_global == "Micro"
)

cor_df_4g <- FetchData(
  micro_4g,
  vars = c("Scoring", "CytoTRACE2_Score")
) %>%
  transmute(
    Scoring = as.numeric(Scoring),
    CytoTRACE2_Score = as.numeric(CytoTRACE2_Score)
  ) %>%
  drop_na() %>%
  mutate(
    CytoTRACE2_percentile = percent_rank(CytoTRACE2_Score)
  )


# ------------------------------
# G.2 Pearson correlation analysis
# ------------------------------

pearson_test_4g <- cor.test(
  cor_df_4g$Scoring,
  cor_df_4g$CytoTRACE2_percentile,
  method = "pearson"
)

cor_annotation_4g <- sprintf(
  paste0(
    "Microglia, Pearson: r = %.2f, 95%% CI [%.2f, %.2f], ",
    "p = %.2e, n = %d"
  ),
  unname(pearson_test_4g$estimate),
  pearson_test_4g$conf.int[1],
  pearson_test_4g$conf.int[2],
  pearson_test_4g$p.value,
  nrow(cor_df_4g)
)

zero_margin_4g <- theme(
  plot.margin = margin(0, 0, 0, 0)
)


# ------------------------------
# G.3 Top marginal histogram
# ------------------------------

p3g_top <- ggplot(
  cor_df_4g,
  aes(x = Scoring)
) +
  geom_histogram(
    bins = 40,
    fill = "#2CA25F",
    color = "black",
    linewidth = 0.25
  ) +
  labs(
    x = NULL,
    y = "Count"
  ) +
  theme_classic(base_size = 9) +
  zero_margin_4g +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  )


# ------------------------------
# G.4 Correlation scatter plot
# ------------------------------

p3g_scatter <- ggplot(
  cor_df_4g,
  aes(
    x = Scoring,
    y = CytoTRACE2_percentile
  )
) +
  geom_point(
    size = 0.5,
    alpha = 0.35
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "blue",
    se = FALSE,
    linewidth = 0.8
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Scoring",
    y = "CytoTRACE2 score percentile"
  ) +
  theme_classic(base_size = 10) +
  zero_margin_4g


# ------------------------------
# G.5 Right marginal histogram
# ------------------------------

p3g_right <- ggplot(
  cor_df_4g,
  aes(x = CytoTRACE2_percentile)
) +
  geom_histogram(
    bins = 40,
    fill = "#FDAE6B",
    color = "black",
    linewidth = 0.25
  ) +
  coord_flip() +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Count"
  ) +
  theme_classic(base_size = 9) +
  zero_margin_4g +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 7)
  )


# ------------------------------
# G.6 Combine correlation panel
# ------------------------------

p3g_blank <- ggplot() +
  theme_void() +
  zero_margin_4g

p3g_bottom <- cowplot::plot_grid(
  p3g_scatter,
  p3g_right,
  nrow = 1,
  rel_widths = c(1, 0.18),
  align = "h",
  axis = "tb"
)

p3g_toprow <- cowplot::plot_grid(
  p3g_top,
  p3g_blank,
  nrow = 1,
  rel_widths = c(1, 0.18),
  align = "h"
)

p3g_main <- cowplot::plot_grid(
  p3g_toprow,
  p3g_bottom,
  ncol = 1,
  rel_heights = c(0.18, 1),
  align = "v"
)

p3g <- ggdraw(p3g_main) +
  draw_label(
    cor_annotation_4g,
    x = 0.02,
    y = 0.99,
    hjust = 0,
    vjust = 1,
    size = 8.5
  )

print(p3g)


# ------------------------------
# G.7 Save Figure 3G
# ------------------------------

ggsave(
  filename = "Figure3G_Scoring_CytoTRACE2_correlation.svg",
  plot = p3g,
  width = 4.2,
  height = 3.8,
  units = "in",
  device = grDevices::svg,
  bg = "white"
)


# ==============================================================================
# Optional reproducibility information
# ==============================================================================

capture.output(
  sessionInfo(),
  file = "Figure3_sessionInfo.txt"
)

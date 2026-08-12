# ==============================================================================
# Figure 2 — Single-cell overview, marker visualization, and TWAS activity
# Publication-ready analysis script
#
# Required input:
#   scRNA_harmony : Seurat object containing:
#     - clusters
#     - diagnosis
#     - celltype_global
#     - orig.ident
#     - RNA assay with counts and normalized data
#     - UMAP reduction
#
# Figure panels:
#   2A  Cluster-level circular embedding with diagnosis composition ring
#   2B  Cell-type circular embedding with diagnosis track
#   2C  Canonical marker-gene dot plot
#   2D  Nebulosa marker-density maps
#   2E  AD-specific five-method TWAS activity dot plot
#   2F  AD-specific five-method TWAS activity distributions
#   2G  AD-specific composite TWAS activity density landscape
#   2H  Control versus AD composite TWAS activity at donor level
#
# Scoring design for 2E–2H:
#   Positive/Negative gene sets are imported directly from Figure 1.
#   The five algorithms are calculated on ALL cells first.
#   A single composite Scoring value is then calculated on ALL cells.
#   Panels 2E–2G subset diagnosis == "AD" for visualization.
#   Panel 2H compares Control versus AD using the same global Scoring scale.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages and global settings
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(plot1cell)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(MASS)
  library(scales)
  library(wordcloud)
  library(ggplot2)
  library(svglite)
  library(Nebulosa)
  library(AUCell)
  library(UCell)
  library(GSVA)
  library(singscore)
})

set.seed(1234)

output_dir <- "Figure2"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

diagnosis_colors <- c(
  Control = "#80BCC8",
  AD      = "#D88F91"
)

celltype_colors <- c(
  Oligo = "#D9B44A",
  OPC   = "#E57C5B",
  Astro = "#B24B5A",
  Exc   = "#8CC1E6",
  Inh   = "#9CCB7A",
  Micro = "#7A55A3",
  Vasc  = "#3B78C2"
)

celltype_order <- c(
  "Exc",
  "Inh",
  "Astro",
  "Micro",
  "OPC",
  "Oligo",
  "Vasc"
)

required_meta <- c(
  "clusters",
  "diagnosis",
  "celltype_global",
  "orig.ident"
)

missing_meta <- setdiff(
  required_meta,
  colnames(scRNA_harmony@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

unexpected_diagnosis <- setdiff(
  unique(na.omit(as.character(scRNA_harmony$diagnosis))),
  c("Control", "AD")
)

if (length(unexpected_diagnosis) > 0) {
  stop(
    "Unexpected diagnosis levels: ",
    paste(unexpected_diagnosis, collapse = ", ")
  )
}


# ==============================================================================
# 2A. Cluster-level circular embedding with diagnosis composition ring
# ==============================================================================

cluster_col <- "clusters"
diagnosis_col <- "diagnosis"

cluster_numeric <- suppressWarnings(
  as.integer(as.character(scRNA_harmony@meta.data[[cluster_col]]))
)

cluster_levels <- if (all(!is.na(cluster_numeric))) {
  as.character(sort(unique(cluster_numeric)))
} else {
  sort(unique(as.character(scRNA_harmony@meta.data[[cluster_col]])))
}

scRNA_harmony[[cluster_col]] <- factor(
  as.character(scRNA_harmony@meta.data[[cluster_col]]),
  levels = cluster_levels
)

Idents(scRNA_harmony) <- cluster_col

cluster_circ_data <- prepare_circlize_data(
  scRNA_harmony,
  scale = 0.68
)

if (!"Cluster" %in% names(cluster_circ_data)) {
  cluster_circ_data$Cluster <- as.character(Idents(scRNA_harmony))
}

cluster_circ_data$Cluster <- factor(
  as.character(cluster_circ_data$Cluster),
  levels = cluster_levels
)

if (!diagnosis_col %in% names(cluster_circ_data)) {
  cluster_circ_data[[diagnosis_col]] <-
    scRNA_harmony@meta.data[[diagnosis_col]]
}

max_cells_for_plot <- 80000L

if (nrow(cluster_circ_data) > max_cells_for_plot) {
  cluster_circ_data <- cluster_circ_data[
    sample.int(nrow(cluster_circ_data), max_cells_for_plot),
    ,
    drop = FALSE
  ]
}

cluster_colors <- setNames(
  grDevices::hcl.colors(
    length(cluster_levels),
    palette = "Dark 3"
  ),
  cluster_levels
)

cluster_meta <- scRNA_harmony@meta.data
cluster_meta[[cluster_col]] <- factor(
  as.character(cluster_meta[[cluster_col]]),
  levels = cluster_levels
)
cluster_meta[[diagnosis_col]] <- as.character(
  cluster_meta[[diagnosis_col]]
)

diagnosis_table <- table(
  cluster_meta[[cluster_col]],
  cluster_meta[[diagnosis_col]]
)

diagnosis_prop <- as.data.frame.matrix(
  prop.table(diagnosis_table, margin = 1)
)

if (!"Control" %in% colnames(diagnosis_prop)) {
  diagnosis_prop$Control <- 0
}
if (!"AD" %in% colnames(diagnosis_prop)) {
  diagnosis_prop$AD <- 0
}

diagnosis_prop <- diagnosis_prop[
  cluster_levels,
  c("Control", "AD"),
  drop = FALSE
]

plot_cluster_circlize <- function(
    data_plot,
    cluster_colors,
    show_center_labels = FALSE,
    contour_levels = c(0.2, 0.3),
    point_size = 0.002,
    kde_n = 120
) {

  centers <- data_plot %>%
    group_by(Cluster) %>%
    summarise(
      x = median(x, na.rm = TRUE),
      y = median(y, na.rm = TRUE),
      .groups = "drop"
    )

  density_est <- MASS::kde2d(
    data_plot$x,
    data_plot$y,
    n = kde_n
  )

  cluster_names <- levels(data_plot$Cluster)

  data_plot$Colors <- unname(
    cluster_colors[as.character(data_plot$Cluster)]
  )

  circos.clear()
  par(bg = "white", mar = c(0, 0, 0, 0), xpd = FALSE)

  circos.par(
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.01, 0),
    track.height = 0.01,
    gap.degree = c(rep(2, length(cluster_names) - 1), 12),
    points.overflow.warning = FALSE
  )

  circos.initialize(
    sectors = data_plot$Cluster,
    x = data_plot$x_polar2
  )

  circos.track(
    factors = data_plot$Cluster,
    x = data_plot$x_polar2,
    y = data_plot$dim2,
    bg.border = NA,
    panel.fun = function(x, y) {
      circos.axis(
        major.tick = TRUE,
        labels = FALSE,
        col = "black",
        lwd = 0.8
      )
    }
  )

  for (ct in cluster_names) {
    subset_ct <- data_plot[
      data_plot$Cluster == ct,
      ,
      drop = FALSE
    ]

    if (nrow(subset_ct) == 0) next

    circos.segments(
      x0 = min(subset_ct$x_polar2, na.rm = TRUE),
      y0 = 0,
      x1 = max(subset_ct$x_polar2, na.rm = TRUE),
      y1 = 0,
      col = cluster_colors[ct],
      lwd = 3,
      sector.index = ct
    )
  }

  graphics::points(
    data_plot$x,
    data_plot$y,
    pch = 19,
    col = scales::alpha(data_plot$Colors, 0.20),
    cex = point_size
  )

  graphics::contour(
    density_est,
    drawlabels = FALSE,
    levels = contour_levels,
    col = "#AE9C76",
    add = TRUE
  )

  if (show_center_labels) {
    wordcloud::textplot(
      x = centers$x,
      y = centers$y,
      words = centers$Cluster,
      cex = 0.6,
      new = FALSE,
      show.lines = FALSE
    )
  }
}

draw_cluster_panel <- function(
    filename,
    show_center_labels = FALSE
) {

  png(
    filename = filename,
    width = 8,
    height = 8,
    units = "in",
    res = 600,
    bg = "white"
  )

  plot_cluster_circlize(
    data_plot = cluster_circ_data,
    cluster_colors = cluster_colors,
    show_center_labels = show_center_labels
  )

  circos.trackPlotRegion(
    factors = cluster_levels,
    ylim = c(0, 1),
    track.height = 0.012,
    bg.border = NA,
    panel.fun = function(x, y) {

      sector <- as.character(
        get.cell.meta.data("sector.index")
      )
      xlim <- get.cell.meta.data("xlim")
      sector_width <- diff(xlim)

      p_control <- diagnosis_prop[sector, "Control"]
      p_ad <- diagnosis_prop[sector, "AD"]

      total <- p_control + p_ad

      if (is.na(total) || total == 0) {
        p_control <- 0
        p_ad <- 0
      } else {
        p_control <- p_control / total
        p_ad <- p_ad / total
      }

      circos.rect(
        xlim[1],
        0,
        xlim[1] + p_control * sector_width,
        1,
        col = diagnosis_colors["Control"],
        border = NA
      )

      circos.rect(
        xlim[1] + p_control * sector_width,
        0,
        xlim[2],
        1,
        col = diagnosis_colors["AD"],
        border = NA
      )
    }
  )

  dev.off()
  circos.clear()
}

draw_cluster_panel(
  file.path(output_dir, "Fig_2A_cluster_circlize_no_labels.png"),
  FALSE
)

draw_cluster_panel(
  file.path(output_dir, "Fig_2A_cluster_circlize_with_labels.png"),
  TRUE
)


# ==============================================================================
# 2B. Cell-type circular embedding with diagnosis track
# ==============================================================================

celltype_obj <- scRNA_harmony
Idents(celltype_obj) <- celltype_obj$celltype_global

celltype_levels <- levels(Idents(celltype_obj))

missing_celltype_colors <- setdiff(
  celltype_levels,
  names(celltype_colors)
)

if (length(missing_celltype_colors) > 0) {
  stop(
    "Missing cell-type colors for: ",
    paste(missing_celltype_colors, collapse = ", ")
  )
}

celltype_circ_data <- prepare_circlize_data(
  celltype_obj,
  scale = 0.75
)

plot_celltype_circlize <- function(
    data_plot,
    colors,
    contour_levels = c(0.2, 0.3),
    point_size = 0.01,
    kde_n = 200
) {

  density_est <- MASS::kde2d(
    data_plot$x,
    data_plot$y,
    n = kde_n
  )

  celltypes <- names(table(data_plot$Cluster))

  color_df <- data.frame(
    Cluster = celltypes,
    color = unname(colors[celltypes]),
    stringsAsFactors = FALSE
  )

  original_order <- rownames(data_plot)

  data_plot <- merge(
    data_plot,
    color_df,
    by = "Cluster"
  )

  rownames(data_plot) <- data_plot$cells
  data_plot <- data_plot[original_order, , drop = FALSE]
  data_plot$Colors <- data_plot$color

  circos.clear()
  par(bg = "white", mar = c(0, 0, 0, 0), xpd = FALSE)

  circos.par(
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.01, 0),
    track.height = 0.01,
    gap.degree = c(rep(2, length(celltypes) - 1), 12),
    points.overflow.warning = FALSE
  )

  circos.initialize(
    sectors = data_plot$Cluster,
    x = data_plot$x_polar2
  )

  circos.track(
    data_plot$Cluster,
    data_plot$x_polar2,
    y = data_plot$dim2,
    bg.border = NA,
    panel.fun = function(x, y) {
      circos.axis(
        labels = FALSE,
        major.tick = FALSE
      )
    }
  )

  for (ct in celltypes) {
    subset_ct <- data_plot[
      data_plot$Cluster == ct,
      ,
      drop = FALSE
    ]

    if (nrow(subset_ct) == 0) next

    circos.segments(
      x0 = min(subset_ct$x_polar2, na.rm = TRUE),
      y0 = 0,
      x1 = max(subset_ct$x_polar2, na.rm = TRUE),
      y1 = 0,
      col = colors[ct],
      lwd = 3,
      sector.index = ct
    )
  }

  graphics::points(
    data_plot$x,
    data_plot$y,
    pch = 19,
    col = scales::alpha(data_plot$Colors, 0.20),
    cex = point_size
  )

  graphics::contour(
    density_est,
    drawlabels = FALSE,
    levels = contour_levels,
    col = "#AE9C76",
    add = TRUE
  )
}

add_diagnosis_track <- function(
    data_plot,
    colors,
    track_lwd = 3
) {

  circos.track(
    data_plot$Cluster,
    data_plot$x_polar2,
    y = data_plot$dim2,
    bg.border = NA
  )

  celltypes <- names(table(data_plot$Cluster))

  for (ct in celltypes) {

    subset_ct <- data_plot[
      data_plot$Cluster == ct,
      ,
      drop = FALSE
    ]

    diagnosis_names <- names(table(subset_ct$diagnosis))

    missing_colors <- setdiff(
      diagnosis_names,
      names(colors)
    )

    if (length(missing_colors) > 0) {
      stop(
        "Missing diagnosis colors for: ",
        paste(missing_colors, collapse = ", ")
      )
    }

    segment_colors <- unname(
      colors[diagnosis_names]
    )

    segment_start <- get_segment(
      subset_ct,
      group = "diagnosis"
    )

    segment_end <- c(
      segment_start[-1] - 1,
      nrow(subset_ct)
    )

    scale_factor <-
      max(subset_ct$x_polar2, na.rm = TRUE) /
      nrow(subset_ct)

    circos.segments(
      x0 = segment_start * scale_factor,
      y0 = 0,
      x1 = segment_end * scale_factor,
      y1 = 0,
      col = segment_colors,
      sector.index = ct,
      lwd = track_lwd
    )
  }
}

png(
  filename = file.path(
    output_dir,
    "Fig_2B_celltype_circlize.png"
  ),
  width = 7,
  height = 7,
  units = "in",
  res = 600,
  bg = "white"
)

plot_celltype_circlize(
  data_plot = celltype_circ_data,
  colors = celltype_colors
)

add_diagnosis_track(
  data_plot = celltype_circ_data,
  colors = diagnosis_colors
)

dev.off()
circos.clear()


# ==============================================================================
# 2C. Canonical marker-gene dot plot
# ==============================================================================

DefaultAssay(scRNA_harmony) <- "RNA"

marker_genes <- c(
  "SLC17A7",
  "GAD1",
  "GAD2",
  "AQP4",
  "CX3CR1",
  "CSF1R",
  "PDGFRA",
  "PLP1",
  "CLDN5"
)

missing_marker_genes <- setdiff(
  marker_genes,
  rownames(scRNA_harmony[["RNA"]])
)

if (length(missing_marker_genes) > 0) {
  stop(
    "Marker genes not found in RNA assay: ",
    paste(missing_marker_genes, collapse = ", ")
  )
}

expr_marker <- GetAssayData(
  scRNA_harmony,
  assay = "RNA",
  layer = "data"
)[marker_genes, , drop = FALSE]

celltype_group <- factor(
  scRNA_harmony$celltype_global,
  levels = celltype_order
)

marker_summary <- lapply(
  marker_genes,
  function(gene) {

    values <- as.numeric(expr_marker[gene, ])
    by_group <- split(values, celltype_group)

    avg_exp <- vapply(
      by_group,
      function(x) mean(expm1(x), na.rm = TRUE),
      numeric(1)
    )

    pct_exp <- vapply(
      by_group,
      function(x) mean(x > 0, na.rm = TRUE) * 100,
      numeric(1)
    )

    data.frame(
      gene = gene,
      celltype = names(avg_exp),
      avg_exp = as.numeric(avg_exp),
      pct_exp = as.numeric(pct_exp)
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    gene = factor(gene, levels = marker_genes),
    celltype = factor(celltype, levels = celltype_order),
    avg_exp_log1p = log1p(avg_exp),
    pct_exp = pmin(pct_exp, 75)
  )

p_2C <- ggplot(
  marker_summary,
  aes(
    x = celltype,
    y = gene
  )
) +
  geom_point(
    aes(
      size = pct_exp,
      color = avg_exp_log1p
    )
  ) +
  scale_size(
    range = c(2, 9),
    limits = c(0, 75),
    name = "Cells expressing (%)"
  ) +
  scale_color_gradientn(
    colours = c(
      "#39489F",
      "#39BBEC",
      "#F9ED36",
      "#F38466",
      "#B81F25"
    ),
    values = scales::rescale(
      c(0, 0.5, 1.0, 1.5, 2.0)
    ),
    limits = c(0, 2),
    oob = scales::squish,
    name = "Average expression\n(log1p)"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.grid.major = element_line(
      color = "grey95"
    ),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(
    output_dir,
    "Fig_2C_marker_dotplot.svg"
  ),
  plot = p_2C,
  width = 5,
  height = 5,
  units = "in"
)


# ==============================================================================
# 2D. Nebulosa density maps for canonical marker genes
# ==============================================================================

nebulosa_dir <- file.path(
  output_dir,
  "Fig_2D_Nebulosa_markers"
)

dir.create(
  nebulosa_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

for (gene in marker_genes) {

  p_gene_list <- Nebulosa::plot_density(
    scRNA_harmony,
    features = gene,
    slot = "data",
    reduction = "umap",
    size = 0.3,
    combine = FALSE
  )

  p_gene <- p_gene_list[[1]]

  ggsave(
    filename = file.path(
      nebulosa_dir,
      paste0("Nebulosa_", gene, ".png")
    ),
    plot = p_gene,
    width = 5,
    height = 4.5,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}



# ==============================================================================
# PREPARATION BEFORE FIGURE 2E
# Five-method TWAS-DEG gene-set scoring across ALL cells
#
# Positive and negative gene sets are taken directly from Figure 1 output:
#   TWAS_DEG_gene_sets.csv
#
# Figure 1 defines these genes by combining:
#   TWAS P-value threshold + TWAS direction + cell-type-specific DEG support.
#
# The same Positive/Negative gene sets are then scored across ALL cells using:
#   AUCell, UCell, ssGSEA, AddModuleScore, and singscore.
# Panels 2E–2G subsequently subset diagnosis == "AD".
# Panel 2H uses the same global composite Scoring for Control versus AD.
# ===============================================================================

TWAS_DEG_FILE <- "TWAS_DEG_gene_sets.csv"

if (!file.exists(TWAS_DEG_FILE)) {
  stop(
    "Figure 1 gene-set file not found: ",
    TWAS_DEG_FILE,
    ". Run Figure 1 first to generate TWAS_DEG_gene_sets.csv."
  )
}

twas_deg_table <- read.csv(
  TWAS_DEG_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_gene_set_columns <- c(
  "direction",
  "gene"
)

missing_gene_set_columns <- setdiff(
  required_gene_set_columns,
  colnames(twas_deg_table)
)

if (length(missing_gene_set_columns) > 0) {
  stop(
    "Figure 1 gene-set file is missing required columns: ",
    paste(
      missing_gene_set_columns,
      collapse = ", "
    )
  )
}

twas_gene_sets <- list(
  Positive = unique(
    twas_deg_table$gene[
      twas_deg_table$direction == "Positive"
    ]
  ),
  Negative = unique(
    twas_deg_table$gene[
      twas_deg_table$direction == "Negative"
    ]
  )
)

twas_gene_sets <- lapply(
  twas_gene_sets,
  function(x) {
    x <- trimws(as.character(x))
    unique(
      x[
        !is.na(x) & nzchar(x)
      ]
    )
  }
)

if (
  length(twas_gene_sets$Positive) == 0 ||
  length(twas_gene_sets$Negative) == 0
) {
  stop(
    "Figure 1 gene-set file contains an empty Positive or Negative gene set."
  )
}

DefaultAssay(scRNA_harmony) <- "RNA"

expr_counts <- GetAssayData(
  scRNA_harmony,
  assay = "RNA",
  layer = "counts"
)

expr_data <- GetAssayData(
  scRNA_harmony,
  assay = "RNA",
  layer = "data"
)

genes_present <- rownames(expr_counts)

matched_counts <- vapply(
  twas_gene_sets,
  function(x) sum(x %in% genes_present),
  integer(1)
)

message(
  "Figure 1 TWAS-DEG genes before assay matching: Positive = ",
  length(twas_gene_sets$Positive),
  "; Negative = ",
  length(twas_gene_sets$Negative)
)

message(
  "Figure 1 TWAS-DEG genes matched to RNA assay: Positive = ",
  matched_counts["Positive"],
  "; Negative = ",
  matched_counts["Negative"]
)

twas_gene_sets <- lapply(
  twas_gene_sets,
  function(x) intersect(x, genes_present)
)

if (
  length(twas_gene_sets$Positive) < 5 ||
  length(twas_gene_sets$Negative) < 5
) {
  stop(
    "Too few Figure 1 TWAS-DEG genes matched the RNA assay. ",
    "Positive = ",
    length(twas_gene_sets$Positive),
    "; Negative = ",
    length(twas_gene_sets$Negative)
  )
}

write.csv(
  data.frame(
    direction = c(
      rep("Positive", length(twas_gene_sets$Positive)),
      rep("Negative", length(twas_gene_sets$Negative))
    ),
    gene = c(
      twas_gene_sets$Positive,
      twas_gene_sets$Negative
    )
  ),
  file.path(
    output_dir,
    "Figure1_TWAS_DEG_gene_sets_used_for_Figure2.csv"
  ),
  row.names = FALSE
)

# 2. Helper functions
# ------------------------------------------------------------------------------

minmax01 <- function(x) {

  x <- as.numeric(x)
  valid <- is.finite(x)

  out <- rep(
    NA_real_,
    length(x)
  )

  if (sum(valid) == 0) {
    return(out)
  }

  x_range <- range(
    x[valid]
  )

  if (
    x_range[2] ==
    x_range[1]
  ) {
    out[valid] <- 0
    return(out)
  }

  out[valid] <- (
    x[valid] -
    x_range[1]
  ) / (
    x_range[2] -
    x_range[1]
  )

  out
}


safe_zscore <- function(x) {

  x <- as.numeric(x)
  valid <- is.finite(x)

  out <- rep(
    NA_real_,
    length(x)
  )

  if (sum(valid) == 0) {
    return(out)
  }

  x_sd <- sd(
    x[valid],
    na.rm = TRUE
  )

  if (
    !is.finite(x_sd) ||
    x_sd == 0
  ) {
    out[valid] <- 0
    return(out)
  }

  out[valid] <- (
    x[valid] -
    mean(
      x[valid],
      na.rm = TRUE
    )
  ) / x_sd

  out
}


safe_row_sum <- function(df) {

  mat <- as.matrix(df)

  sums <- rowSums(
    mat,
    na.rm = TRUE
  )

  n_nonmissing <- rowSums(
    !is.na(mat)
  )

  sums[
    n_nonmissing == 0
  ] <- NA_real_

  sums
}


safe_scale_vector <- function(x) {

  x <- as.numeric(x)
  valid <- is.finite(x)

  out <- rep(
    NA_real_,
    length(x)
  )

  if (sum(valid) == 0) {
    return(out)
  }

  x_sd <- sd(
    x[valid],
    na.rm = TRUE
  )

  if (
    !is.finite(x_sd) ||
    x_sd == 0
  ) {
    out[valid] <- 0
    return(out)
  }

  out[valid] <- as.numeric(
    scale(
      x[valid]
    )
  )

  out
}


# ------------------------------------------------------------------------------
# 3. AUCell
# ------------------------------------------------------------------------------

auc_rankings <- AUCell::AUCell_buildRankings(
  expr_counts,
  plotStats = FALSE
)

auc_scores <- AUCell::AUCell_calcAUC(
  twas_gene_sets,
  auc_rankings,
  aucMaxRank = ceiling(
    nrow(
      auc_rankings
    ) * 0.05
  )
)

auc_matrix <- AUCell::getAUC(
  auc_scores
)

scRNA_harmony$AUCell_Positive <- as.numeric(
  auc_matrix[
    "Positive",
    colnames(
      scRNA_harmony
    )
  ]
)

scRNA_harmony$AUCell_Negative <- as.numeric(
  auc_matrix[
    "Negative",
    colnames(
      scRNA_harmony
    )
  ]
)

scRNA_harmony$AUCell_Positive_z <- safe_zscore(
  scRNA_harmony$AUCell_Positive
)

scRNA_harmony$AUCell_Negative_z <- safe_zscore(
  scRNA_harmony$AUCell_Negative
)

scRNA_harmony$AUCell_diff_z <-
  scRNA_harmony$AUCell_Positive_z -
  scRNA_harmony$AUCell_Negative_z

scRNA_harmony$AUCell_diff_norm <- minmax01(
  scRNA_harmony$AUCell_diff_z
)


# ------------------------------------------------------------------------------
# 4. UCell
# ------------------------------------------------------------------------------

metadata_before_ucell <- colnames(
  scRNA_harmony@meta.data
)

scRNA_harmony <- UCell::AddModuleScore_UCell(
  scRNA_harmony,
  features = twas_gene_sets,
  name = "TWAS_UCell"
)

new_ucell_cols <- setdiff(
  colnames(
    scRNA_harmony@meta.data
  ),
  metadata_before_ucell
)

positive_ucell_col <- new_ucell_cols[
  grepl(
    "Positive",
    new_ucell_cols,
    ignore.case = TRUE
  )
]

negative_ucell_col <- new_ucell_cols[
  grepl(
    "Negative",
    new_ucell_cols,
    ignore.case = TRUE
  )
]

if (
  length(positive_ucell_col) < 1 ||
  length(negative_ucell_col) < 1
) {
  stop(
    "Could not identify Positive/Negative UCell columns."
  )
}

scRNA_harmony$UCell_Positive <- as.numeric(
  scRNA_harmony@meta.data[
    ,
    positive_ucell_col[1]
  ]
)

scRNA_harmony$UCell_Negative <- as.numeric(
  scRNA_harmony@meta.data[
    ,
    negative_ucell_col[1]
  ]
)

scRNA_harmony$UCell_Positive_z <- safe_zscore(
  scRNA_harmony$UCell_Positive
)

scRNA_harmony$UCell_Negative_z <- safe_zscore(
  scRNA_harmony$UCell_Negative
)

scRNA_harmony$UCell_diff_z <-
  scRNA_harmony$UCell_Positive_z -
  scRNA_harmony$UCell_Negative_z

scRNA_harmony$UCell_diff_norm <- minmax01(
  scRNA_harmony$UCell_diff_z
)


# ------------------------------------------------------------------------------
# 5. ssGSEA
# ------------------------------------------------------------------------------

ssgsea_parameters <- GSVA::ssgseaParam(
  exprData = expr_data,
  geneSets = twas_gene_sets,
  minSize = 5,
  maxSize = 500,
  normalize = FALSE
)

ssgsea_matrix <- GSVA::gsva(
  ssgsea_parameters
)

common_cells_ssgsea <- intersect(
  colnames(
    scRNA_harmony
  ),
  colnames(
    ssgsea_matrix
  )
)

scRNA_harmony$ssgsea_Positive <- NA_real_
scRNA_harmony$ssgsea_Negative <- NA_real_

scRNA_harmony@meta.data[
  common_cells_ssgsea,
  "ssgsea_Positive"
] <- as.numeric(
  ssgsea_matrix[
    "Positive",
    common_cells_ssgsea
  ]
)

scRNA_harmony@meta.data[
  common_cells_ssgsea,
  "ssgsea_Negative"
] <- as.numeric(
  ssgsea_matrix[
    "Negative",
    common_cells_ssgsea
  ]
)

scRNA_harmony$ssgsea_Positive_z <- safe_zscore(
  scRNA_harmony$ssgsea_Positive
)

scRNA_harmony$ssgsea_Negative_z <- safe_zscore(
  scRNA_harmony$ssgsea_Negative
)

scRNA_harmony$ssgsea_diff_z <-
  scRNA_harmony$ssgsea_Positive_z -
  scRNA_harmony$ssgsea_Negative_z

scRNA_harmony$ssgsea_diff_norm <- minmax01(
  scRNA_harmony$ssgsea_diff_z
)


# ------------------------------------------------------------------------------
# 6. AddModuleScore
# ------------------------------------------------------------------------------

metadata_before_addmodule <- colnames(
  scRNA_harmony@meta.data
)

scRNA_harmony <- Seurat::AddModuleScore(
  object = scRNA_harmony,
  features = twas_gene_sets,
  name = "TWAS_AddModuleScore",
  seed = 1234
)

new_addmodule_cols <- setdiff(
  colnames(
    scRNA_harmony@meta.data
  ),
  metadata_before_addmodule
)

addmodule_cols <- new_addmodule_cols[
  grepl(
    "^TWAS_AddModuleScore",
    new_addmodule_cols
  )
]

if (length(addmodule_cols) < 2) {
  stop(
    "AddModuleScore did not generate two score columns."
  )
}

# The feature-list order is Positive first and Negative second.
scRNA_harmony$AddModuleScore_Positive <- as.numeric(
  scRNA_harmony@meta.data[
    ,
    addmodule_cols[1]
  ]
)

scRNA_harmony$AddModuleScore_Negative <- as.numeric(
  scRNA_harmony@meta.data[
    ,
    addmodule_cols[2]
  ]
)

scRNA_harmony$AddModuleScore_Positive_z <- safe_zscore(
  scRNA_harmony$AddModuleScore_Positive
)

scRNA_harmony$AddModuleScore_Negative_z <- safe_zscore(
  scRNA_harmony$AddModuleScore_Negative
)

scRNA_harmony$AddModuleScore_diff_z <-
  scRNA_harmony$AddModuleScore_Positive_z -
  scRNA_harmony$AddModuleScore_Negative_z

scRNA_harmony$AddModuleScore_diff_norm <- minmax01(
  scRNA_harmony$AddModuleScore_diff_z
)


# ------------------------------------------------------------------------------
# 7. singscore
# ------------------------------------------------------------------------------

singscore_chunk_size <- 2000L

all_cells <- colnames(
  scRNA_harmony
)

singscore_positive <- setNames(
  rep(
    NA_real_,
    length(all_cells)
  ),
  all_cells
)

singscore_negative <- setNames(
  rep(
    NA_real_,
    length(all_cells)
  ),
  all_cells
)

for (
  start_index in seq(
    1,
    length(all_cells),
    by = singscore_chunk_size
  )
) {

  index <- start_index:min(
    start_index +
      singscore_chunk_size -
      1,
    length(all_cells)
  )

  cells_chunk <- all_cells[
    index
  ]

  expr_chunk <- as.matrix(
    expr_data[
      ,
      cells_chunk,
      drop = FALSE
    ]
  )

  rank_data <- singscore::rankGenes(
    expr_chunk
  )

  positive_result <- singscore::simpleScore(
    rank_data,
    upSet = twas_gene_sets$Positive
  )

  negative_result <- singscore::simpleScore(
    rank_data,
    upSet = twas_gene_sets$Negative
  )

  singscore_positive[
    cells_chunk
  ] <- as.numeric(
    positive_result$TotalScore
  )

  singscore_negative[
    cells_chunk
  ] <- as.numeric(
    negative_result$TotalScore
  )
}

scRNA_harmony$singscore_Positive <- as.numeric(
  singscore_positive[
    colnames(
      scRNA_harmony
    )
  ]
)

scRNA_harmony$singscore_Negative <- as.numeric(
  singscore_negative[
    colnames(
      scRNA_harmony
    )
  ]
)

scRNA_harmony$singscore_Positive_z <- safe_zscore(
  scRNA_harmony$singscore_Positive
)

scRNA_harmony$singscore_Negative_z <- safe_zscore(
  scRNA_harmony$singscore_Negative
)

scRNA_harmony$singscore_diff_z <-
  scRNA_harmony$singscore_Positive_z -
  scRNA_harmony$singscore_Negative_z

scRNA_harmony$singscore_diff_norm <- minmax01(
  scRNA_harmony$singscore_diff_z
)


# ------------------------------------------------------------------------------
# 8. Columns used directly by Figures 2E–2G
# ------------------------------------------------------------------------------

twas_z_methods <- c(
  "AUCell_diff_z",
  "UCell_diff_z",
  "AddModuleScore_diff_z",
  "ssgsea_diff_z",
  "singscore_diff_z"
)

twas_norm_methods <- c(
  "AUCell_diff_norm",
  "UCell_diff_norm",
  "AddModuleScore_diff_norm",
  "ssgsea_diff_norm",
  "singscore_diff_norm"
)

required_score_columns <- c(
  twas_z_methods,
  twas_norm_methods
)

missing_score_columns <- setdiff(
  required_score_columns,
  colnames(
    scRNA_harmony@meta.data
  )
)

if (length(missing_score_columns) > 0) {
  stop(
    "Five-method scoring failed to generate: ",
    paste(
      missing_score_columns,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------------------------
# Global composite TWAS activity score
# ------------------------------------------------------------------------------

# Each method-specific diff_norm score is already scaled to 0-1 across ALL cells.
# To keep Panels 2E-2H directly comparable, the five-method composite score is
# also calculated ONCE across all cells before any diagnosis-specific subsetting.

score_matrix_global <- as.matrix(
  scRNA_harmony@meta.data[
    ,
    twas_norm_methods,
    drop = FALSE
  ]
)

complete_score_rows <- rowSums(
  is.finite(score_matrix_global)
) == length(twas_norm_methods)

scoring_raw_global <- rep(
  NA_real_,
  nrow(score_matrix_global)
)

scoring_raw_global[complete_score_rows] <- rowSums(
  score_matrix_global[
    complete_score_rows,
    ,
    drop = FALSE
  ]
)

scRNA_harmony$Scoring_raw <- scoring_raw_global
scRNA_harmony$Scoring <- minmax01(
  scRNA_harmony$Scoring_raw
)

# Save the reusable object only after all five method scores and the global
# composite Scoring column have been created.
saveRDS(
  scRNA_harmony,
  file.path(
    output_dir,
    "scRNA_harmony_with_5method_TWAS_scores.rds"
  )
)


# ==============================================================================
# 2E. AD-specific Figure 1 TWAS-DEG activity across all major cell types
# ==============================================================================

ad_twas_meta <- scRNA_harmony@meta.data %>%
  as.data.frame() %>%
  filter(
    diagnosis == "AD",
    celltype_global %in% celltype_order
  )

twas_methods_all <- c(
  twas_z_methods,
  "Scoring"
)

ad_twas_long <- ad_twas_meta %>%
  rownames_to_column(
    "cell"
  ) %>%
  select(
    cell,
    celltype_global,
    all_of(
      twas_methods_all
    )
  ) %>%
  pivot_longer(
    cols = all_of(
      twas_methods_all
    ),
    names_to = "method",
    values_to = "score"
  )

# For the five method-specific Z-scores, activity-positive cells are defined as
# score > 0. For the global 0-1 composite Scoring value, the display threshold
# is Scoring > 0.5. These thresholds are used only for dot size in Panel 2E.

ad_twas_summary <- ad_twas_long %>%
  group_by(
    celltype_global,
    method
  ) %>%
  summarise(
    avg = mean(
      score,
      na.rm = TRUE
    ),
    pct = if (
      first(method) == "Scoring"
    ) {
      mean(
        score > 0.5,
        na.rm = TRUE
      ) * 100
    } else {
      mean(
        score > 0,
        na.rm = TRUE
      ) * 100
    },
    .groups = "drop"
  ) %>%
  group_by(
    method
  ) %>%
  mutate(
    avg_scaled = safe_scale_vector(
      avg
    )
  ) %>%
  ungroup() %>%
  mutate(
    avg_scaled = pmax(
      pmin(
        avg_scaled,
        2
      ),
      -1
    ),
    pct = ifelse(
      is.na(pct),
      0,
      pct
    ),
    method = factor(
      method,
      levels = twas_methods_all
    ),
    celltype_global = factor(
      celltype_global,
      levels = celltype_order
    )
  )

p_2E <- ggplot(
  ad_twas_summary,
  aes(
    x = method,
    y = celltype_global
  )
) +
  geom_point(
    aes(
      size = pct,
      color = avg_scaled
    ),
    shape = 16,
    stroke = 0
  ) +
  scale_color_gradientn(
    colors = c(
      "#39489F",
      "#39BBEC",
      "#F9ED36",
      "#F38466",
      "#B81F25"
    ),
    name = "Scaled activity",
    limits = c(-1, 2),
    breaks = c(-1, 0, 1, 2),
    oob = scales::squish
  ) +
  scale_size_continuous(
    range = c(1.5, 4),
    name = "Cells above threshold (%)"
  ) +
  labs(
    x = "Scoring method",
    y = "Cell type"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    panel.grid.major = element_line(
      color = "grey95",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    legend.position = "right"
  )

ggsave(
  filename = file.path(
    output_dir,
    "Fig_2E_AD_TWAS_activity_dotplot.svg"
  ),
  plot = p_2E,
  width = 4.5,
  height = 3.5,
  units = "in"
)


# ==============================================================================
# 2F. AD-specific five-method TWAS activity distributions
# ==============================================================================

twas_score_columns <- c(
  "AUCell_diff_z",
  "UCell_diff_z",
  "AddModuleScore_diff_z",
  "ssgsea_diff_z",
  "singscore_diff_z",
  "Scoring"
)

twas_method_labels <- c(
  AUCell_diff_z = "AUCell",
  UCell_diff_z = "UCell",
  AddModuleScore_diff_z = "AddModuleScore",
  ssgsea_diff_z = "ssGSEA",
  singscore_diff_z = "singscore",
  Scoring = "Scoring"
)

ad_twas_violin <- ad_twas_meta %>%
  select(
    celltype_global,
    all_of(
      twas_score_columns
    )
  ) %>%
  pivot_longer(
    cols = all_of(
      twas_score_columns
    ),
    names_to = "method",
    values_to = "score"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = rev(
        twas_score_columns
      ),
      labels = rev(
        unname(
          twas_method_labels[
            twas_score_columns
          ]
        )
      )
    ),
    score = as.numeric(
      score
    ),
    celltype_global = factor(
      celltype_global,
      levels = celltype_order
    )
  ) %>%
  filter(
    !is.na(score),
    !is.na(celltype_global)
  )

p_2F <- ggplot(
  ad_twas_violin,
  aes(
    x = celltype_global,
    y = score,
    fill = celltype_global
  )
) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    adjust = 1.2,
    color = "black",
    linewidth = 0.3
  ) +
  geom_boxplot(
    width = 0.10,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.3
  ) +
  facet_grid(
    method ~ .,
    scales = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = celltype_colors
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    panel.spacing.y = grid::unit(
      0.10,
      "lines"
    ),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      hjust = 1
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    legend.position = "none"
  )

ggsave(
  filename = file.path(
    output_dir,
    "Fig_2F_AD_TWAS_activity_violin.svg"
  ),
  plot = p_2F,
  width = 8,
  height = 8,
  units = "in"
)


# ==============================================================================
# 2G. AD-only composite TWAS activity density landscape
# ==============================================================================

ad_scoring_df <- ad_twas_meta %>%
  transmute(
    celltype = factor(
      celltype_global,
      levels = celltype_order
    ),
    score = as.numeric(
      Scoring
    )
  ) %>%
  filter(
    !is.na(celltype),
    !is.na(score)
  )

scoring_quantiles <- ad_scoring_df %>%
  group_by(
    celltype
  ) %>%
  summarise(
    p10 = quantile(
      score,
      0.10,
      na.rm = TRUE
    ),
    p25 = quantile(
      score,
      0.25,
      na.rm = TRUE
    ),
    p50 = quantile(
      score,
      0.50,
      na.rm = TRUE
    ),
    p75 = quantile(
      score,
      0.75,
      na.rm = TRUE
    ),
    p90 = quantile(
      score,
      0.90,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

n_density_grid <- 260

score_min <- min(
  ad_scoring_df$score,
  na.rm = TRUE
)

score_max <- max(
  ad_scoring_df$score,
  na.rm = TRUE
)

score_grid <- seq(
  score_min,
  score_max,
  length.out = n_density_grid
)

scoring_density_df <- bind_rows(
  lapply(
    levels(
      ad_scoring_df$celltype
    ),
    function(ct) {

      values <- ad_scoring_df$score[
        ad_scoring_df$celltype == ct
      ]

      values <- values[
        is.finite(values)
      ]

      if (length(values) < 2) {
        return(
          data.frame(
            celltype = ct,
            score = score_grid,
            density_scaled = 0
          )
        )
      }

      density_est <- density(
        values,
        from = score_min,
        to = score_max,
        n = n_density_grid,
        na.rm = TRUE
      )

      density_scaled <- density_est$y

      if (
        max(
          density_scaled
        ) > 0
      ) {
        density_scaled <-
          density_scaled /
          max(
            density_scaled
          )
      } else {
        density_scaled <- rep(
          0,
          length(
            density_scaled
          )
        )
      }

      data.frame(
        celltype = ct,
        score = density_est$x,
        density_scaled = density_scaled
      )
    }
  )
) %>%
  mutate(
    celltype = factor(
      celltype,
      levels = celltype_order
    )
  )

p_2G <- ggplot() +
  geom_raster(
    data = scoring_density_df,
    aes(
      x = celltype,
      y = score,
      fill = density_scaled
    )
  ) +
  geom_line(
    data = scoring_quantiles,
    aes(
      x = celltype,
      y = p50,
      group = 1
    ),
    linewidth = 1.0,
    color = "#B2182B"
  ) +
  geom_line(
    data = scoring_quantiles,
    aes(
      x = celltype,
      y = p25,
      group = 1
    ),
    linetype = "dashed",
    linewidth = 0.7,
    color = "#B2182B"
  ) +
  geom_line(
    data = scoring_quantiles,
    aes(
      x = celltype,
      y = p75,
      group = 1
    ),
    linetype = "dashed",
    linewidth = 0.7,
    color = "#B2182B"
  ) +
  geom_line(
    data = scoring_quantiles,
    aes(
      x = celltype,
      y = p10,
      group = 1
    ),
    linetype = "dotted",
    linewidth = 0.6,
    color = "#B2182B"
  ) +
  geom_line(
    data = scoring_quantiles,
    aes(
      x = celltype,
      y = p90,
      group = 1
    ),
    linetype = "dotted",
    linewidth = 0.6,
    color = "#B2182B"
  ) +
  scale_fill_gradientn(
    colours = c(
      "#39489F",
      "#39BBEC",
      "#F9ED36",
      "#F38466",
      "#B81F25"
    ),
    values = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),
    limits = c(0, 1),
    name = "Relative density"
  ) +
  labs(
    x = NULL,
    y = "TWAS activity score"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.line = element_line(
      linewidth = 0.6
    ),
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(
    output_dir,
    "Fig_2G_AD_Scoring_density.svg"
  ),
  plot = p_2G,
  width = 9,
  height = 5,
  units = "in"
)

# ==============================================================================
# 2H. Control versus AD composite TWAS activity at donor level
# ==============================================================================

# Statistical unit: donor/sample (orig.ident), not individual cells.
# The same global Scoring column used in Panels 2E-2G is used here.

control_ad_cell_scores <- scRNA_harmony@meta.data %>%
  as.data.frame() %>%
  transmute(
    sample = as.character(orig.ident),
    diagnosis = factor(
      diagnosis,
      levels = c("Control", "AD")
    ),
    celltype = factor(
      celltype_global,
      levels = celltype_order
    ),
    score = as.numeric(Scoring)
  ) %>%
  filter(
    !is.na(sample),
    nzchar(sample),
    !is.na(diagnosis),
    !is.na(celltype),
    is.finite(score)
  )

# Collapse cells to one donor-level value per cell type. Median is used to reduce
# sensitivity to extreme single-cell scores.
control_ad_donor_scores <- control_ad_cell_scores %>%
  group_by(
    sample,
    diagnosis,
    celltype
  ) %>%
  summarise(
    score = median(
      score,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  )

# Ensure that a sample is not assigned to more than one diagnosis.
ambiguous_samples <- control_ad_donor_scores %>%
  distinct(
    sample,
    diagnosis
  ) %>%
  count(sample) %>%
  filter(n > 1)

if (nrow(ambiguous_samples) > 0) {
  stop(
    "At least one orig.ident is associated with more than one diagnosis: ",
    paste(
      ambiguous_samples$sample,
      collapse = ", "
    )
  )
}

wilcoxon_results <- control_ad_donor_scores %>%
  group_by(celltype) %>%
  summarise(
    n_control = sum(
      diagnosis == "Control"
    ),
    n_ad = sum(
      diagnosis == "AD"
    ),
    p_value = if (
      n_distinct(diagnosis) == 2 &&
      sum(diagnosis == "Control") >= 2 &&
      sum(diagnosis == "AD") >= 2
    ) {
      wilcox.test(
        score ~ droplevels(diagnosis),
        exact = FALSE
      )$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(
      p_value,
      method = "BH"
    ),
    significance = case_when(
      is.na(p_adj) ~ "",
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

annotation_positions <- control_ad_donor_scores %>%
  group_by(celltype) %>%
  summarise(
    ymin = min(
      score,
      na.rm = TRUE
    ),
    ymax = max(
      score,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    y = pmin(
      1.05,
      ymax + 0.08 * pmax(
        ymax - ymin,
        0.05
      )
    )
  ) %>%
  left_join(
    wilcoxon_results,
    by = "celltype"
  )

p_2H <- ggplot(
  control_ad_donor_scores,
  aes(
    x = celltype,
    y = score,
    fill = diagnosis
  )
) +
  geom_violin(
    position = position_dodge(
      width = 0.85
    ),
    trim = FALSE,
    scale = "width",
    adjust = 1.2,
    color = "black",
    linewidth = 0.3,
    alpha = 0.70
  ) +
  geom_boxplot(
    position = position_dodge(
      width = 0.85
    ),
    width = 0.16,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.3,
    alpha = 0.85
  ) +
  geom_point(
    aes(
      color = diagnosis,
      group = diagnosis
    ),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.85
    ),
    size = 1.4,
    alpha = 0.85,
    show.legend = FALSE
  ) +
  geom_text(
    data = annotation_positions,
    aes(
      x = celltype,
      y = y,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = diagnosis_colors
  ) +
  scale_color_manual(
    values = diagnosis_colors
  ) +
  scale_y_continuous(
    limits = c(0, 1.08),
    expand = expansion(
      mult = c(0.01, 0.02)
    )
  ) +
  labs(
    x = NULL,
    y = "TWAS activity score",
    fill = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    legend.position = "top",
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(
    output_dir,
    "Fig_2H_Control_vs_AD_Scoring_donor_level.svg"
  ),
  plot = p_2H,
  width = 10,
  height = 5,
  units = "in"
)

write.csv(
  control_ad_donor_scores,
  file = file.path(
    output_dir,
    "Fig_2H_donor_level_scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  wilcoxon_results,
  file = file.path(
    output_dir,
    "Fig_2H_Control_vs_AD_Wilcoxon_BH.csv"
  ),
  row.names = FALSE
)


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

message("Figure 2A-2H analysis completed successfully using Figure 1 TWAS-DEG gene sets.")

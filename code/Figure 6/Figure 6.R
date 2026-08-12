# ==============================================================================
# Figure 6 | Expression of LCORL, ZFP36L2 and ZNF652
#
# Figure 6A: Gene expression across major cell types
# Figure 6B: Microglia expression in Control vs AD + differential-expression test
# Figure 6C: Gene-expression dot plot in Control vs AD
#
# Notes:
#   - All three panels use the same genes: LCORL, ZFP36L2 and ZNF652.
#   - Cells with gene expression equal to 0 are retained in all analyses.
#   - Figure 6B contains violin plots only; statistical testing is exported
#     separately and is not annotated on the figure.
#   - No local/absolute paths are used.
#   - All figures are saved as SVG files.
# ==============================================================================


# ==============================================================================
# 0. Packages and shared settings
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# Genes displayed in Figure 6A-C
genes <- c("LCORL", "ZFP36L2", "ZNF652")

# Use the RNA assay explicitly.
if (!"RNA" %in% Assays(scRNA_harmony)) {
  stop("RNA assay was not found in scRNA_harmony.")
}

DefaultAssay(scRNA_harmony) <- "RNA"

# Ensure the normalized RNA data layer is available.
rna_data <- tryCatch(
  GetAssayData(
    scRNA_harmony,
    assay = "RNA",
    layer = "data"
  ),
  error = function(e) NULL
)

if (
  is.null(rna_data) ||
  nrow(rna_data) == 0 ||
  ncol(rna_data) == 0
) {
  scRNA_harmony <- NormalizeData(
    scRNA_harmony,
    verbose = FALSE
  )
}

# Confirm that all requested genes are present.
genes_present <- intersect(
  genes,
  rownames(scRNA_harmony[["RNA"]])
)

genes_missing <- setdiff(
  genes,
  genes_present
)

if (length(genes_missing) > 0) {
  warning(
    "The following genes were not found and will be skipped: ",
    paste(genes_missing, collapse = ", ")
  )
}

if (length(genes_present) == 0) {
  stop("None of the requested genes were found in scRNA_harmony.")
}


required_meta_6 <- c(
  "celltype_global",
  "diagnosis",
  "orig.ident"
)

missing_meta_6 <- setdiff(
  required_meta_6,
  colnames(scRNA_harmony@meta.data)
)

if (length(missing_meta_6) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(
      missing_meta_6,
      collapse = ", "
    )
  )
}

cluster_colors <- c(
  Oligo = "#D9B44A",
  OPC   = "#E57C5B",
  Astro = "#B24B5A",
  Exc   = "#8CC1E6",
  Inh   = "#9CCB7A",
  Micro = "#7A55A3",
  Vasc  = "#3B78C2"
)

diagnosis_colors <- c(
  Control = "#80BCC8",
  AD      = "#D88F91"
)


# ==============================================================================
# Figure 6A | Expression across major cell types
# ==============================================================================

# Figure 6A uses all cells in scRNA_harmony.
# IMPORTANT: cells with expression = 0 are retained.

p6a <- VlnPlot(
  object = scRNA_harmony,
  features = genes_present,
  group.by = "celltype_global",
  pt.size = 0,
  cols = cluster_colors,
  ncol = length(genes_present),
  combine = TRUE
) &
  theme_classic(base_size = 11) &
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_text(
      face = "italic",
      hjust = 0.5,
      size = 12
    )
  )

print(p6a)

ggsave(
  filename = "Figure6A_LCORL_ZFP36L2_ZNF652_by_celltype.svg",
  plot = p6a,
  width = 12,
  height = 5,
  units = "in",
  device = grDevices::svg
)


# ==============================================================================
# Figure 6B | Microglia Control vs AD violin plots
# ==============================================================================

# ----------------------------------------------------------------------------
# 6B.1 Prepare microglia
# ----------------------------------------------------------------------------

micro_6b <- subset(
  scRNA_harmony,
  subset = celltype_global == "Micro" & diagnosis %in% c("Control", "AD")
)

micro_6b$diagnosis <- factor(
  micro_6b$diagnosis,
  levels = c("Control", "AD")
)

# IMPORTANT:
# Do NOT remove cells with LCORL/ZFP36L2/ZNF652 expression equal to 0.
# Every Control and AD microglial cell is retained.


# ----------------------------------------------------------------------------
# 6B.2 Violin plots only
# ----------------------------------------------------------------------------

p6b <- VlnPlot(
  object = micro_6b,
  features = genes_present,
  group.by = "diagnosis",
  pt.size = 0,
  cols = diagnosis_colors,
  ncol = length(genes_present),
  combine = TRUE
) &
  theme_classic(base_size = 11) &
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 10
    ),
    axis.text.y = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_text(
      face = "italic",
      hjust = 0.5,
      size = 12
    )
  )

print(p6b)

ggsave(
  filename = "Figure6B_Microglia_Control_vs_AD_violin.svg",
  plot = p6b,
  width = 10,
  height = 5,
  units = "in",
  device = grDevices::svg
)


# ----------------------------------------------------------------------------
# 6B.3 Donor-level Control vs AD differential-expression test
# ----------------------------------------------------------------------------

# The violin plot remains cell-level for visualization.
# Formal statistical testing is performed at the donor/sample level using
# orig.ident to avoid treating cells from the same donor as independent
# biological replicates.

expr_6b <- FetchData(
  micro_6b,
  vars = c(
    "diagnosis",
    "orig.ident",
    genes_present
  )
)

donor_level_6b <- do.call(
  rbind,
  lapply(
    genes_present,
    function(gene) {

      gene_df <- expr_6b %>%
        dplyr::select(
          diagnosis,
          orig.ident,
          all_of(gene)
        )

      colnames(gene_df)[3] <- "expression"

      gene_df <- gene_df %>%
        group_by(
          orig.ident,
          diagnosis
        ) %>%
        summarise(
          expression = median(
            expression,
            na.rm = TRUE
          ),
          .groups = "drop"
        )

      if (
        n_distinct(gene_df$diagnosis) < 2
      ) {
        return(
          data.frame(
            gene = gene,
            n_AD = sum(
              gene_df$diagnosis == "AD"
            ),
            n_Control = sum(
              gene_df$diagnosis == "Control"
            ),
            median_AD = NA_real_,
            median_Control = NA_real_,
            p_value = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }

      test_result <- wilcox.test(
        expression ~ diagnosis,
        data = gene_df,
        exact = FALSE
      )

      data.frame(
        gene = gene,
        n_AD = sum(
          gene_df$diagnosis == "AD"
        ),
        n_Control = sum(
          gene_df$diagnosis == "Control"
        ),
        median_AD = median(
          gene_df$expression[
            gene_df$diagnosis == "AD"
          ],
          na.rm = TRUE
        ),
        median_Control = median(
          gene_df$expression[
            gene_df$diagnosis == "Control"
          ],
          na.rm = TRUE
        ),
        p_value = test_result$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
)

donor_level_6b$p_adj_BH <- p.adjust(
  donor_level_6b$p_value,
  method = "BH"
)

donor_level_6b$avg_log2FC <- NA_real_

for (i in seq_len(nrow(donor_level_6b))) {

  if (
    is.finite(donor_level_6b$median_AD[i]) &&
    is.finite(donor_level_6b$median_Control[i])
  ) {

    donor_level_6b$avg_log2FC[i] <-
      donor_level_6b$median_AD[i] -
      donor_level_6b$median_Control[i]
  }
}

print(donor_level_6b)

write.csv(
  donor_level_6b,
  file = "Figure6B_Microglia_AD_vs_Control_donor_level_Wilcoxon.csv",
  row.names = FALSE
)


# ==============================================================================
# Figure 6C | Control vs AD dot plot
# ==============================================================================

# Figure 6C follows the original DotPlot logic and uses the full Seurat object.
# DotPlot calculations include cells with expression = 0; no nonzero filtering
# is performed before plotting.

obj_6c <- scRNA_harmony

obj_6c$diagnosis <- factor(
  obj_6c$diagnosis,
  levels = c("Control", "AD")
)

p6c <- DotPlot(
  object = obj_6c,
  features = genes_present,
  group.by = "diagnosis",
  dot.scale = 6
) +
  RotatedAxis() +
  scale_color_gradientn(
    colors = c(
      "#39489F",
      "#39BBEC",
      "#F9ED36",
      "#F38466",
      "#B81F25"
    )
  ) +
  scale_size_continuous(
    name = "Percent expressed"
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Average expression"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      family = "Arial",
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(
      family = "Arial",
      size = 11,
      color = "black"
    ),
    axis.ticks = element_blank()
  )

print(p6c)

ggsave(
  filename = "Figure6C_LCORL_ZFP36L2_ZNF652_by_diagnosis_DotPlot.svg",
  plot = p6c,
  width = 7,
  height = 4.5,
  units = "in",
  device = grDevices::svg
)

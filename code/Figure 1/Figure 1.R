# ==============================================================================
# Figure 1 | TWAS association, cell-type-specific DEGs, and functional enrichment
#
# Figure 1A: TWAS Manhattan plot
# Figure 1B: GO enrichment analysis and dot plot
# Figure 1C: KEGG pathway enrichment analysis and dot plot
#
# Species: Homo sapiens
# Output format: SVG
#
# Analysis overview:
#   1. Read and annotate TWAS results
#   2. Generate the TWAS Manhattan plot
#   3. Perform cell-type-specific AD vs Control differential expression
#   4. Intersect TWAS-associated genes with DEGs
#   5. Perform GO and KEGG enrichment on the combined TWAS-DEG gene set
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages and analysis settings
# ------------------------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "ggplot2",
  "biomaRt",
  "clusterProfiler",
  "org.Hs.eg.db",
  "stringr",
  "Seurat"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "The following packages are not installed:\n",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(biomaRt)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  library(Seurat)
})

set.seed(1234)

# Input files should be placed in the working directory or supplied as
# project-relative paths.
TWAS_FILE <- "SMulTiXcan.tsv"
DEG_OUTPUT_FILE <- "DEG_AD_vs_Control_by_celltype_minpct0.25.csv"

TWAS_P_CUTOFF <- 0.01
DEG_PADJ_CUTOFF <- 0.05

ENRICHMENT_P_CUTOFF <- 0.05
ENRICHMENT_Q_CUTOFF <- 0.20

MANHATTAN_OUTPUT <- "Figure1A_TWAS_Manhattan.svg"
GO_OUTPUT <- "Figure1B_GO_enrichment.svg"
KEGG_OUTPUT <- "Figure1C_KEGG_enrichment.svg"

if (!exists("scRNA_harmony")) {
  stop(
    "The Seurat object 'scRNA_harmony' is not available in the R environment."
  )
}

if (!file.exists(TWAS_FILE)) {
  stop(
    "TWAS input file not found: ",
    TWAS_FILE
  )
}


# ==============================================================================
# 1. Common TWAS data preparation
# ==============================================================================

twas_data <- read.delim(
  TWAS_FILE,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_twas_columns <- c(
  "gene",
  "pvalue",
  "z_mean"
)

missing_twas_columns <- setdiff(
  required_twas_columns,
  colnames(twas_data)
)

if (length(missing_twas_columns) > 0) {
  stop(
    "TWAS data is missing required columns:\n",
    paste(
      missing_twas_columns,
      collapse = ", "
    )
  )
}

# Remove Ensembl version suffixes.
twas_data$gene <- sub(
  "\\..*$",
  "",
  trimws(
    as.character(
      twas_data$gene
    )
  )
)

ensembl_ids <- unique(
  na.omit(
    twas_data$gene
  )
)

ensembl_ids <- ensembl_ids[
  nzchar(
    ensembl_ids
  )
]

if (length(ensembl_ids) == 0) {
  stop(
    "No valid Ensembl gene IDs were found."
  )
}


# ------------------------------------------------------------------------------
# 1.1 Gene annotation
# ------------------------------------------------------------------------------

required_annotation_columns <- c(
  "gene_name",
  "chromosome_name",
  "start_position",
  "end_position"
)

annotation_available <- all(
  required_annotation_columns %in%
    colnames(twas_data)
)

if (!annotation_available) {

  ensembl_mart <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = "hsapiens_gene_ensembl"
  )

  gene_annotation <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "hgnc_symbol",
      "chromosome_name",
      "start_position",
      "end_position"
    ),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl_mart
  ) %>%
    mutate(
      hgnc_symbol = trimws(
        as.character(
          hgnc_symbol
        )
      ),
      canonical_chr = chromosome_name %in% c(
        as.character(
          1:22
        ),
        "X",
        "Y"
      ),
      has_hgnc = !is.na(hgnc_symbol) &
        nzchar(hgnc_symbol)
    ) %>%
    arrange(
      ensembl_gene_id,
      desc(canonical_chr),
      desc(has_hgnc)
    ) %>%
    distinct(
      ensembl_gene_id,
      .keep_all = TRUE
    ) %>%
    transmute(
      ensembl_gene_id,
      hgnc_symbol_mart = hgnc_symbol,
      chromosome_name_mart = chromosome_name,
      start_position_mart = start_position,
      end_position_mart = end_position
    )

  twas_data <- twas_data %>%
    left_join(
      gene_annotation,
      by = c(
        "gene" = "ensembl_gene_id"
      )
    )

  if (!"gene_name" %in% colnames(twas_data)) {
    twas_data$gene_name <- NA_character_
  }

  if (!"chromosome_name" %in% colnames(twas_data)) {
    twas_data$chromosome_name <- NA_character_
  }

  if (!"start_position" %in% colnames(twas_data)) {
    twas_data$start_position <- NA_real_
  }

  if (!"end_position" %in% colnames(twas_data)) {
    twas_data$end_position <- NA_real_
  }

  twas_data <- twas_data %>%
    mutate(
      gene_name = coalesce(
        na_if(
          trimws(
            as.character(
              gene_name
            )
          ),
          ""
        ),
        na_if(
          trimws(
            as.character(
              hgnc_symbol_mart
            )
          ),
          ""
        ),
        gene
      ),
      chromosome_name = coalesce(
        na_if(
          trimws(
            as.character(
              chromosome_name
            )
          ),
          ""
        ),
        na_if(
          trimws(
            as.character(
              chromosome_name_mart
            )
          ),
          ""
        )
      ),
      start_position = coalesce(
        suppressWarnings(
          as.numeric(
            start_position
          )
        ),
        suppressWarnings(
          as.numeric(
            start_position_mart
          )
        )
      ),
      end_position = coalesce(
        suppressWarnings(
          as.numeric(
            end_position
          )
        ),
        suppressWarnings(
          as.numeric(
            end_position_mart
          )
        )
      )
    ) %>%
    select(
      -any_of(
        c(
          "hgnc_symbol_mart",
          "chromosome_name_mart",
          "start_position_mart",
          "end_position_mart"
        )
      )
    )
}

# Gene midpoint for Manhattan plotting.
twas_data <- twas_data %>%
  mutate(
    BP = (
      as.numeric(
        start_position
      ) +
        as.numeric(
          end_position
        )
    ) / 2
  )

write.csv(
  twas_data,
  "SMulTiXcan_annotated.csv",
  row.names = FALSE
)


# ==============================================================================
# 2. Figure 1A | TWAS Manhattan plot
# ==============================================================================

manhattan_data <- twas_data %>%
  transmute(
    Gene = as.character(
      gene
    ),
    Gene_name = as.character(
      gene_name
    ),
    Chromosome = suppressWarnings(
      as.integer(
        chromosome_name
      )
    ),
    Position = as.numeric(
      BP
    ),
    P_value = as.numeric(
      pvalue
    )
  ) %>%
  filter(
    Chromosome %in% 1:22,
    !is.na(Position),
    is.finite(Position),
    !is.na(P_value),
    is.finite(P_value),
    P_value > 0,
    P_value <= 1
  ) %>%
  arrange(
    Chromosome,
    Position
  )

if (nrow(manhattan_data) == 0) {
  stop(
    "No valid autosomal TWAS records were available for Manhattan plotting."
  )
}

chromosome_offsets <- manhattan_data %>%
  group_by(
    Chromosome
  ) %>%
  summarise(
    chromosome_max_position = max(
      Position,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    Chromosome
  ) %>%
  mutate(
    offset = cumsum(
      chromosome_max_position
    ) -
      chromosome_max_position
  )

manhattan_data <- manhattan_data %>%
  left_join(
    chromosome_offsets,
    by = "Chromosome"
  ) %>%
  mutate(
    cumulative_position =
      Position +
      offset
  )

chromosome_axis <- manhattan_data %>%
  group_by(
    Chromosome
  ) %>%
  summarise(
    center = (
      min(
        cumulative_position,
        na.rm = TRUE
      ) +
        max(
          cumulative_position,
          na.rm = TRUE
        )
    ) / 2,
    .groups = "drop"
  ) %>%
  arrange(
    Chromosome
  )

nominal_twas_points <- manhattan_data %>%
  filter(
    P_value < TWAS_P_CUTOFF
  )

other_twas_points <- manhattan_data %>%
  filter(
    P_value >= TWAS_P_CUTOFF
  )

p1a <- ggplot() +
  geom_point(
    data = other_twas_points,
    aes(
      x = cumulative_position,
      y = -log10(P_value),
      color = factor(Chromosome)
    ),
    size = 0.8,
    alpha = 0.80
  ) +
  geom_point(
    data = nominal_twas_points,
    aes(
      x = cumulative_position,
      y = -log10(P_value)
    ),
    color = "#D62728",
    size = 1.0,
    alpha = 0.95
  ) +
  geom_hline(
    yintercept = -log10(
      TWAS_P_CUTOFF
    ),
    linetype = 2,
    linewidth = 0.6,
    color = "grey40"
  ) +
  scale_x_continuous(
    breaks = chromosome_axis$center,
    labels = chromosome_axis$Chromosome,
    expand = expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +
  guides(
    color = "none"
  ) +
  labs(
    x = "Chromosome",
    y = expression(
      -log[10](P)
    )
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    axis.text = element_text(
      color = "black"
    )
  )

print(p1a)

ggsave(
  filename = MANHATTAN_OUTPUT,
  plot = p1a,
  device = grDevices::svg,
  width = 15,
  height = 5,
  bg = "white"
)


# ==============================================================================
# 3. Cell-type-specific AD versus Control differential expression
# ==============================================================================

DefaultAssay(
  scRNA_harmony
) <- "RNA"

scRNA_harmony <- NormalizeData(
  scRNA_harmony
)

if (!"diagnosis" %in% colnames(scRNA_harmony@meta.data)) {
  stop(
    "Metadata column 'diagnosis' is required."
  )
}

if (!"celltype_global" %in% colnames(scRNA_harmony@meta.data)) {
  stop(
    "Metadata column 'celltype_global' is required."
  )
}

diagnosis_values <- unique(
  na.omit(
    as.character(
      scRNA_harmony$diagnosis
    )
  )
)

if (
  !"AD" %in% diagnosis_values ||
  !"Control" %in% diagnosis_values
) {
  stop(
    "The 'diagnosis' metadata must contain both 'AD' and 'Control'."
  )
}

Idents(
  scRNA_harmony
) <- "diagnosis"

celltypes <- sort(
  unique(
    na.omit(
      as.character(
        scRNA_harmony$celltype_global
      )
    )
  )
)

degs <- vector(
  mode = "list",
  length = length(celltypes)
)

names(degs) <- celltypes

for (ct in celltypes) {

  subset_obj <- subset(
    scRNA_harmony,
    subset =
      celltype_global == ct
  )

  degs[[ct]] <- FindMarkers(
    subset_obj,
    ident.1 = "AD",
    ident.2 = "Control",
    test.use = "wilcox",
    min.pct = 0.25,
    logfc.threshold = 0,
    verbose = FALSE
  )
}

deg_total <- bind_rows(
  lapply(
    names(degs),
    function(ct) {

      deg_df <- as.data.frame(
        degs[[ct]]
      )

      if (nrow(deg_df) == 0) {
        return(
          data.frame(
            gene = character(),
            cluster = character()
          )
        )
      }

      deg_df$gene <- rownames(
        deg_df
      )

      deg_df$cluster <- ct

      deg_df
    }
  )
)

write.csv(
  deg_total,
  DEG_OUTPUT_FILE,
  row.names = FALSE
)


# ==============================================================================
# 4. TWAS-positive and TWAS-negative gene sets
# ==============================================================================

genes_pos_twas <- twas_data %>%
  filter(
    !is.na(pvalue),
    as.numeric(
      pvalue
    ) <= TWAS_P_CUTOFF,
    !is.na(z_mean),
    as.numeric(
      z_mean
    ) > 0,
    !is.na(gene_name),
    nzchar(
      trimws(
        as.character(
          gene_name
        )
      )
    )
  ) %>%
  pull(
    gene_name
  ) %>%
  as.character() %>%
  unique()

genes_neg_twas <- twas_data %>%
  filter(
    !is.na(pvalue),
    as.numeric(
      pvalue
    ) < TWAS_P_CUTOFF,
    !is.na(z_mean),
    as.numeric(
      z_mean
    ) < 0,
    !is.na(gene_name),
    nzchar(
      trimws(
        as.character(
          gene_name
        )
      )
    )
  ) %>%
  pull(
    gene_name
  ) %>%
  as.character() %>%
  unique()


# ==============================================================================
# 5. Intersect TWAS-associated genes with cell-type-specific DEGs
#
# Final criteria:
#   Positive:
#     TWAS P <= 0.01
#     TWAS z_mean > 0
#     DEG adjusted P <= 0.05
#
#   Negative:
#     TWAS P < 0.01
#     TWAS z_mean < 0
#     DEG adjusted P <= 0.05
#
# Positive and negative genes are retained separately here, then combined for
# the downstream GO and KEGG enrichment analyses.
# ==============================================================================

genes_deg_pos <- unique(
  deg_total$gene[
    deg_total$gene %in% genes_pos_twas &
      deg_total$p_val_adj <= DEG_PADJ_CUTOFF
  ]
)

genes_deg_pos <- genes_deg_pos[
  !is.na(
    genes_deg_pos
  )
]

genes_deg_neg <- unique(
  deg_total$gene[
    deg_total$gene %in% genes_neg_twas &
      deg_total$p_val_adj <= DEG_PADJ_CUTOFF
  ]
)

genes_deg_neg <- genes_deg_neg[
  !is.na(
    genes_deg_neg
  )
]

gene_list <- list(
  TWAS_DEG_POS = genes_deg_pos,
  TWAS_DEG_NEG = genes_deg_neg
)

write.table(
  gene_list$TWAS_DEG_POS,
  "TWAS_DEG_positive_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  gene_list$TWAS_DEG_NEG,
  "TWAS_DEG_negative_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.csv(
  data.frame(
    direction = c(
      rep(
        "Positive",
        length(
          gene_list$TWAS_DEG_POS
        )
      ),
      rep(
        "Negative",
        length(
          gene_list$TWAS_DEG_NEG
        )
      )
    ),
    gene = c(
      gene_list$TWAS_DEG_POS,
      gene_list$TWAS_DEG_NEG
    )
  ),
  "TWAS_DEG_gene_sets.csv",
  row.names = FALSE
)


# ==============================================================================
# 6. Helper for mixed SYMBOL / ENSEMBL to Entrez conversion
# ==============================================================================

map_mixed_gene_ids_to_entrez <- function(
    gene_ids,
    orgdb
) {

  gene_ids <- unique(
    trimws(
      as.character(
        gene_ids
      )
    )
  )

  gene_ids <- gene_ids[
    !is.na(
      gene_ids
    ) &
      nzchar(
        gene_ids
      )
  ]

  is_ensembl <- grepl(
    "^ENSG[0-9]+",
    gene_ids,
    ignore.case = TRUE
  )

  ensembl_ids_local <- sub(
    "\\..*$",
    "",
    gene_ids[
      is_ensembl
    ]
  )

  symbol_ids_local <- toupper(
    gene_ids[
      !is_ensembl
    ]
  )

  mappings <- list()

  if (length(ensembl_ids_local) > 0) {

    map_ensembl <- suppressMessages(
      clusterProfiler::bitr(
        unique(
          ensembl_ids_local
        ),
        fromType = "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = orgdb
      )
    )

    if (
      !is.null(
        map_ensembl
      ) &&
      nrow(
        map_ensembl
      ) > 0
    ) {

      mappings[[length(mappings) + 1]] <- map_ensembl %>%
        transmute(
          input_gene = ENSEMBL,
          input_type = "ENSEMBL",
          ENTREZID = ENTREZID
        )
    }
  }

  if (length(symbol_ids_local) > 0) {

    map_symbol <- suppressMessages(
      clusterProfiler::bitr(
        unique(
          symbol_ids_local
        ),
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = orgdb
      )
    )

    if (
      !is.null(
        map_symbol
      ) &&
      nrow(
        map_symbol
      ) > 0
    ) {

      mappings[[length(mappings) + 1]] <- map_symbol %>%
        transmute(
          input_gene = SYMBOL,
          input_type = "SYMBOL",
          ENTREZID = ENTREZID
        )
    }
  }

  if (length(mappings) == 0) {

    return(
      tibble(
        input_gene = character(),
        input_type = character(),
        ENTREZID = character()
      )
    )
  }

  bind_rows(
    mappings
  ) %>%
    filter(
      !is.na(
        ENTREZID
      ),
      nzchar(
        as.character(
          ENTREZID
        )
      )
    ) %>%
    distinct(
      input_gene,
      ENTREZID,
      .keep_all = TRUE
    )
}


# ==============================================================================
# 7. Figure 1B | GO enrichment
# ==============================================================================

# Combine positive and negative TWAS-DEG genes as requested.
go_genes <- unique(
  c(
    gene_list$TWAS_DEG_POS,
    gene_list$TWAS_DEG_NEG
  )
)

go_genes <- go_genes[
  !is.na(
    go_genes
  ) &
    nzchar(
      go_genes
    )
]

if (length(go_genes) == 0) {
  stop(
    "No TWAS-DEG genes were available for GO enrichment."
  )
}

write.table(
  go_genes,
  "Figure1B_GO_input_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

go_mapping <- map_mixed_gene_ids_to_entrez(
  go_genes,
  org.Hs.eg.db
)

if (
  is.null(
    go_mapping
  ) ||
  nrow(
    go_mapping
  ) == 0
) {
  stop(
    "No GO input genes could be converted to Entrez IDs."
  )
}

go_entrez_ids <- unique(
  as.character(
    go_mapping$ENTREZID
  )
)

if (length(go_entrez_ids) < 5) {
  stop(
    "Fewer than five genes were available for GO enrichment."
  )
}

write.csv(
  go_mapping,
  "Figure1B_GO_gene_to_Entrez_mapping.csv",
  row.names = FALSE
)

go_bp <- clusterProfiler::enrichGO(
  gene = go_entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = ENRICHMENT_P_CUTOFF,
  qvalueCutoff = ENRICHMENT_Q_CUTOFF,
  readable = TRUE
)

go_cc <- clusterProfiler::enrichGO(
  gene = go_entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = ENRICHMENT_P_CUTOFF,
  qvalueCutoff = ENRICHMENT_Q_CUTOFF,
  readable = TRUE
)

go_mf <- clusterProfiler::enrichGO(
  gene = go_entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = ENRICHMENT_P_CUTOFF,
  qvalueCutoff = ENRICHMENT_Q_CUTOFF,
  readable = TRUE
)

go_bp_result <- as.data.frame(
  go_bp
)

go_cc_result <- as.data.frame(
  go_cc
)

go_mf_result <- as.data.frame(
  go_mf
)

write.csv(
  go_bp_result,
  "Figure1B_GO_BP_complete.csv",
  row.names = FALSE
)

write.csv(
  go_cc_result,
  "Figure1B_GO_CC_complete.csv",
  row.names = FALSE
)

write.csv(
  go_mf_result,
  "Figure1B_GO_MF_complete.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7.1 Representative GO terms for display
# ------------------------------------------------------------------------------

selected_go_terms <- c(
  "endosomal transport",
  "regulation of amyloid-beta clearance",
  "intracellular calcium ion homeostasis",
  "gliogenesis",
  "amyloid precursor protein catabolic process",
  "amyloid-beta formation",
  "regulation of synapse organization",
  "regulation of synaptic plasticity",
  "regulation of reactive oxygen species metabolic process",
  "response to oxidative stress",
  "synaptic vesicle cycle",
  "astrocyte activation"
)

go_bp_selected <- go_bp_result %>%
  filter(
    Description %in%
      selected_go_terms
  ) %>%
  mutate(
    Ontology = "BP"
  )

go_cc_selected <- go_cc_result %>%
  filter(
    Description %in%
      selected_go_terms
  ) %>%
  mutate(
    Ontology = "CC"
  )

go_mf_selected <- go_mf_result %>%
  filter(
    Description %in%
      selected_go_terms
  ) %>%
  mutate(
    Ontology = "MF"
  )

go_selected <- bind_rows(
  go_bp_selected,
  go_cc_selected,
  go_mf_selected
) %>%
  distinct(
    Description,
    .keep_all = TRUE
  )

if (nrow(go_selected) == 0) {
  stop(
    "None of the selected GO terms were found in the enrichment results."
  )
}

missing_go_terms <- setdiff(
  selected_go_terms,
  as.character(
    go_selected$Description
  )
)

if (length(missing_go_terms) > 0) {
  message(
    "Selected GO terms not found: ",
    paste(
      missing_go_terms,
      collapse = "; "
    )
  )
}

write.csv(
  go_selected,
  "Figure1B_GO_selected.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7.2 Prepare and draw GO dot plot
# ------------------------------------------------------------------------------

go_plot_data <- go_selected %>%
  filter(
    !is.na(
      p.adjust
    ),
    p.adjust <= ENRICHMENT_P_CUTOFF,
    !is.na(
      GeneRatio
    ),
    !is.na(
      Count
    )
  )

go_plot_data$GeneRatio_num <- vapply(
  strsplit(
    as.character(
      go_plot_data$GeneRatio
    ),
    "/",
    fixed = TRUE
  ),
  function(parts) {

    if (length(parts) != 2) {
      return(
        NA_real_
      )
    }

    numerator <- suppressWarnings(
      as.numeric(
        parts[1]
      )
    )

    denominator <- suppressWarnings(
      as.numeric(
        parts[2]
      )
    )

    if (
      is.na(
        numerator
      ) ||
      is.na(
        denominator
      ) ||
      denominator == 0
    ) {
      return(
        NA_real_
      )
    }

    numerator /
      denominator
  },
  FUN.VALUE = numeric(1)
)

go_plot_data <- go_plot_data %>%
  filter(
    !is.na(
      GeneRatio_num
    )
  )

go_plot_data$Description <- factor(
  as.character(
    go_plot_data$Description
  ),
  levels = rev(
    selected_go_terms[
      selected_go_terms %in%
        as.character(
          go_plot_data$Description
        )
    ]
  )
)

p1b <- ggplot(
  go_plot_data,
  aes(
    x = GeneRatio_num,
    y = Description,
    fill = -log10(
      p.adjust
    ),
    size = Count
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 1.0,
    alpha = 1
  ) +
  scale_fill_gradientn(
    colors = c(
      "#39489f",
      "#39bbec",
      "#f9ed36",
      "#f38466",
      "#b81f25"
    ),
    name = expression(
      -log[10](
        adjusted~p
      )
    )
  ) +
  scale_size_continuous(
    name = "Count",
    range = c(
      4,
      11
    )
  ) +
  scale_y_discrete(
    labels = function(x) {
      stringr::str_wrap(
        x,
        width = 55
      )
    }
  ) +
  scale_x_continuous(
    name = "Gene ratio",
    expand = expansion(
      mult = c(
        0.02,
        0.08
      )
    )
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.2
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(
      colour = "black"
    ),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      size = 12,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    legend.title = element_text(
      size = 12
    ),
    legend.text = element_text(
      size = 11
    ),
    plot.margin = margin(
      t = 8,
      r = 15,
      b = 8,
      l = 18
    )
  ) +
  guides(
    size = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 1
      )
    ),
    fill = guide_colorbar(
      order = 2
    )
  )

print(p1b)

ggsave(
  filename = GO_OUTPUT,
  plot = p1b,
  device = grDevices::svg,
  width = 9,
  height = 8,
  bg = "white"
)


# ==============================================================================
# 8. Figure 1C | KEGG enrichment
# ==============================================================================

# Use the same combined TWAS-DEG gene set as GO enrichment.
kegg_genes <- unique(
  c(
    gene_list$TWAS_DEG_POS,
    gene_list$TWAS_DEG_NEG
  )
)

kegg_genes <- kegg_genes[
  !is.na(
    kegg_genes
  ) &
    nzchar(
      kegg_genes
    )
]

if (length(kegg_genes) == 0) {
  stop(
    "No TWAS-DEG genes were available for KEGG enrichment."
  )
}

write.table(
  kegg_genes,
  "Figure1C_KEGG_input_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

kegg_mapping <- map_mixed_gene_ids_to_entrez(
  kegg_genes,
  org.Hs.eg.db
)

if (
  is.null(
    kegg_mapping
  ) ||
  nrow(
    kegg_mapping
  ) == 0
) {
  stop(
    "No KEGG input genes could be converted to Entrez IDs."
  )
}

kegg_entrez_ids <- unique(
  as.character(
    kegg_mapping$ENTREZID
  )
)

if (length(kegg_entrez_ids) < 5) {
  stop(
    "Fewer than five genes were available for KEGG enrichment."
  )
}

write.csv(
  kegg_mapping,
  "Figure1C_KEGG_gene_to_Entrez_mapping.csv",
  row.names = FALSE
)

kegg_enrichment <- clusterProfiler::enrichKEGG(
  gene = kegg_entrez_ids,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = ENRICHMENT_P_CUTOFF,
  qvalueCutoff = ENRICHMENT_Q_CUTOFF
)

kegg_enrichment <- clusterProfiler::setReadable(
  kegg_enrichment,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

kegg_result <- as.data.frame(
  kegg_enrichment
)

if (nrow(kegg_result) == 0) {
  stop(
    "KEGG enrichment returned no pathways."
  )
}

write.csv(
  kegg_result,
  "Figure1C_KEGG_complete.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8.1 Representative KEGG pathways for display
# ------------------------------------------------------------------------------

selected_kegg_terms <- c(
  "Axon guidance",
  "Phagosome",
  "Gap junction",
  "Amyotrophic lateral sclerosis",
  "Neurotrophin signaling pathway",
  "Pathways of neurodegeneration - multiple diseases",
  "Regulation of actin cytoskeleton",
  "Adherens junction",
  "Calcium signaling pathway",
  "Cellular senescence"
)

kegg_selected <- kegg_result %>%
  filter(
    Description %in%
      selected_kegg_terms
  )

if (nrow(kegg_selected) == 0) {
  stop(
    "None of the selected KEGG pathways were found in the enrichment results."
  )
}

missing_kegg_terms <- setdiff(
  selected_kegg_terms,
  as.character(
    kegg_selected$Description
  )
)

if (length(missing_kegg_terms) > 0) {
  message(
    "Selected KEGG pathways not found: ",
    paste(
      missing_kegg_terms,
      collapse = "; "
    )
  )
}

write.csv(
  kegg_selected,
  "Figure1C_KEGG_selected.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8.2 Prepare and draw KEGG dot plot
# ------------------------------------------------------------------------------

kegg_plot_data <- kegg_selected %>%
  filter(
    !is.na(
      p.adjust
    ),
    p.adjust <= ENRICHMENT_P_CUTOFF,
    !is.na(
      GeneRatio
    ),
    !is.na(
      Count
    )
  )

kegg_plot_data$GeneRatio_num <- vapply(
  strsplit(
    as.character(
      kegg_plot_data$GeneRatio
    ),
    "/",
    fixed = TRUE
  ),
  function(parts) {

    if (length(parts) != 2) {
      return(
        NA_real_
      )
    }

    numerator <- suppressWarnings(
      as.numeric(
        parts[1]
      )
    )

    denominator <- suppressWarnings(
      as.numeric(
        parts[2]
      )
    )

    if (
      is.na(
        numerator
      ) ||
      is.na(
        denominator
      ) ||
      denominator == 0
    ) {
      return(
        NA_real_
      )
    }

    numerator /
      denominator
  },
  FUN.VALUE = numeric(1)
)

kegg_plot_data <- kegg_plot_data %>%
  filter(
    !is.na(
      GeneRatio_num
    )
  )

kegg_plot_data$Description <- factor(
  as.character(
    kegg_plot_data$Description
  ),
  levels = rev(
    selected_kegg_terms[
      selected_kegg_terms %in%
        as.character(
          kegg_plot_data$Description
        )
    ]
  )
)

p1c <- ggplot(
  kegg_plot_data,
  aes(
    x = GeneRatio_num,
    y = Description,
    fill = -log10(
      p.adjust
    ),
    size = Count
  )
) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 1.0,
    alpha = 1
  ) +
  scale_fill_gradientn(
    colors = c(
      "#39489f",
      "#39bbec",
      "#f9ed36",
      "#f38466",
      "#b81f25"
    ),
    name = expression(
      -log[10](
        adjusted~p
      )
    )
  ) +
  scale_size_continuous(
    name = "Count",
    range = c(
      4,
      11
    )
  ) +
  scale_y_discrete(
    labels = function(x) {
      stringr::str_wrap(
        x,
        width = 50
      )
    }
  ) +
  scale_x_continuous(
    name = "Gene ratio",
    expand = expansion(
      mult = c(
        0.02,
        0.08
      )
    )
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.2
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(
      colour = "black"
    ),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      size = 12,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    legend.title = element_text(
      size = 12
    ),
    legend.text = element_text(
      size = 11
    ),
    plot.margin = margin(
      t = 8,
      r = 15,
      b = 8,
      l = 18
    )
  ) +
  guides(
    size = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 1
      )
    ),
    fill = guide_colorbar(
      order = 2
    )
  )

print(p1c)

ggsave(
  filename = KEGG_OUTPUT,
  plot = p1c,
  device = grDevices::svg,
  width = 9,
  height = 7,
  bg = "white"
)


# ==============================================================================
# 9. Reproducibility information
# ==============================================================================

capture.output(
  sessionInfo(),
  file = "Figure1_sessionInfo.txt"
)

message(
  "Figure 1A-1C analysis completed successfully."
)

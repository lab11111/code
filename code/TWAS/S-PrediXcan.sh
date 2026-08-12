#!/usr/bin/env bash

# =========================
# Paths
# =========================
DATA="$HOME/lab/S-PrediXcan/data"
METAXCAN="$HOME/lab/S-PrediXcan/MetaXcan/software"
OUTPUT="$HOME/lab/S-PrediXcan/output"
MODEL="$HOME/lab/S-PrediXcan/eqtl/mashr"

GWAS_FILE="$DATA/df.txt.gz"
OUTDIR="$OUTPUT/spredixcan/eqtl/utmost_11"

# Create output directory
mkdir -p "$OUTDIR"


# =========================
# Tissues
# =========================
TISSUES=(
    Brain_Cortex
    Brain_Anterior_cingulate_cortex_BA24
    Brain_Hippocampus
    Brain_Amygdala
    Brain_Caudate_basal_ganglia
    Brain_Nucleus_accumbens_basal_ganglia
    Brain_Putamen_basal_ganglia
    Brain_Substantia_nigra
    Brain_Hypothalamus
    Pituitary
    Spleen
)


# =========================
# Run S-PrediXcan
# =========================
for TISSUE in "${TISSUES[@]}"; do

    echo "========================================"
    echo "Running S-PrediXcan: ${TISSUE}"
    echo "========================================"

    python "$METAXCAN/SPrediXcan.py" \
        --gwas_file "$GWAS_FILE" \
        --snp_column panel_variant_id \
        --effect_allele_column effect_allele \
        --non_effect_allele_column other_allele \
        --zscore_column zscore \
        --model_db_path "$MODEL/mashr_${TISSUE}.db" \
        --covariance "$MODEL/mashr_${TISSUE}.txt.gz" \
        --model_db_snp_key varID \
        --keep_non_rsid \
        --additional_output \
        --throw \
        --output_file "$OUTDIR/${TISSUE}.csv"

    echo "Finished: ${TISSUE}"
    echo
done

echo "All tissues completed."

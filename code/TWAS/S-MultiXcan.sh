DATA=$HOME/lab/S-PrediXcan/data
METAXCAN=$HOME/lab/S-PrediXcan/MetaXcan/software
OUTPUT=$HOME/lab/S-PrediXcan/output

# 11个组织的 mashr 模型 db 所在目录（你原来用的那个）
MODEL_DB=$HOME/lab/S-PrediXcan/eqtl/mashr

# 你刚确认存在的 smultixcan 跨组织 SNP covariance 文件所在目录
MODEL_COV=$HOME/lab/S-PrediXcan/eqtl/models/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz

TISSUES_RE='(Brain_Cortex|Brain_Anterior_cingulate_cortex_BA24|Brain_Hippocampus|Brain_Amygdala|Brain_Caudate_basal_ganglia|Brain_Nucleus_accumbens_basal_ganglia|Brain_Putamen_basal_ganglia|Brain_Substantia_nigra|Brain_Hypothalamus|Pituitary|Spleen)'

mkdir -p $OUTPUT/smultixcan/eqtl/utmost_11

python $METAXCAN/SMulTiXcan.py \
  --models_folder $MODEL_DB \
  --models_name_filter \
    mashr_Brain_Cortex.db \
    mashr_Brain_Anterior_cingulate_cortex_BA24.db \
    mashr_Brain_Hippocampus.db \
    mashr_Brain_Amygdala.db \
    mashr_Brain_Caudate_basal_ganglia.db \
    mashr_Brain_Nucleus_accumbens_basal_ganglia.db \
    mashr_Brain_Putamen_basal_ganglia.db \
    mashr_Brain_Substantia_nigra.db \
    mashr_Brain_Hypothalamus.db \
    mashr_Pituitary.db \
    mashr_Spleen.db \
  --snp_covariance $MODEL_COV \
  --metaxcan_folder $OUTPUT/spredixcan/eqtl/utmost_11 \
  --metaxcan_filter "mashr_.*\.csv" \
  --metaxcan_file_name_parse_pattern "(.*)\.csv" \
  --gwas_file  $DATA/df.txt.gz \
  --snp_column panel_variant_id \
  --effect_allele_column effect_allele \
  --non_effect_allele_column other_allele \
  --zscore_column zscore \
  --model_db_snp_key varID \
  --keep_non_rsid \
  --cutoff_threshold 30 \
  --throw \
  --output $OUTPUT/smultixcan/eqtl/utmost_11/SMulTiXcan.tsv


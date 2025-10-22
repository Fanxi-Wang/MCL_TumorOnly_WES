#!/bin/bash
#SBATCH --job-name=merge_vcfs
#SBATCH --output=logs/merge_vcfs_%A.out
#SBATCH --error=logs/merge_vcfs_%A.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-35

# ==============================================================
# Merge Mutect2 and VarDict VCFs
# Author: Fanxi Wang
# Description:
#   Merges tumor-only Mutect2 and VarDict VCF files using 
#   a custom Python2 script adapted from the MC3 pipeline.
# ==============================================================

# Step 1. Load dependencies
module load anaconda3
source ~/.bashrc
conda activate mc3_env   # contains Python 2 + hgsc_vcf dependency

# Step 2. Define directories and script
MERGE_SCRIPT="/path/to/2_variant_calling/merge_vcfs.py"
INPUT_DIR_MUTECT="/path/to/output_variantcalling/mutect2"
INPUT_DIR_VARDICT="/path/to/output_variant_integration/split_vcfs"
OUTPUT_DIR="/path/to/output_variant_integration/merged"

mkdir -p "$OUTPUT_DIR" logs

# Step 3. Define sample list
SAMPLES=("Sample1" "Sample2" "Sample3" "Sample4" "Sample5"
         "Sample6" "Sample7" "Sample8" "Sample9" "Sample10"
         "Sample11" "Sample12" "Sample13" "Sample14" "Sample15"
         "Sample16" "Sample17" "Sample18" "Sample19" "Sample20"
         "Sample21" "Sample22" "Sample23" "Sample24" "Sample25"
         "Sample26" "Sample27" "Sample28" "Sample29" "Sample30"
         "Sample31" "Sample32" "Sample33" "Sample34" "Sample35")

# Step 4. Assign current sample
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}
VCF_MUTECT="${INPUT_DIR_MUTECT}/mutect2_${SAMPLE}_tumor.vcf.gz"
VCF_VARDICT="${INPUT_DIR_VARDICT}/${SAMPLE}_nonSV.vcf"
OUT_VCF="${OUTPUT_DIR}/${SAMPLE}_merged.vcf"

# Step 5. Run merge
echo "[INFO] Merging VCFs for ${SAMPLE}..."
python2 "$MERGE_SCRIPT" \
  --input "$VCF_MUTECT" "$VCF_VARDICT" \
  --output "$OUT_VCF" \
  --keys Mutect2 VarDict

echo "[DONE] Merged VCF generated for ${SAMPLE}: ${OUT_VCF}"

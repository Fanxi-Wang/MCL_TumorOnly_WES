#!/bin/bash
#SBATCH --job-name=split_vcfs
#SBATCH --output=logs/split_vcfs_%A.out
#SBATCH --error=logs/split_vcfs_%A.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-35

# ==============================================================
# Split VarDict VCF into Structural Variant (SV) and non-SV files
# Author: Fanxi Wang
# Description:
#   Uses bcftools to separate VarDict VCFs into:
#     (1) nonSV.vcf  – SNVs and INDELs
#     (2) SV.vcf     – Structural Variants (DEL/DUP/INV/INS)
# ==============================================================

# Step 1. Load dependencies
module load bcftools

# Step 2. Define directories
INPUT_DIR="/path/to/output_variantcalling/vardict"       # Input: VarDict VCFs
OUTPUT_DIR="/path/to/output_variant_integration/split_vcfs"  # Output: split VCFs
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
VCF_IN="${INPUT_DIR}/vardict_${SAMPLE}_tumor.vcf"
VCF_NONSV="${OUTPUT_DIR}/${SAMPLE}_nonSV.vcf"
VCF_SV="${OUTPUT_DIR}/${SAMPLE}_SV.vcf"

echo "[INFO] Splitting VarDict VCF for ${SAMPLE}..."
echo "       Input:  $VCF_IN"
echo "       Output: $VCF_NONSV"
echo "       Output: $VCF_SV"

# Step 5. Split into SV and non-SV
bcftools view -i 'INFO/SVTYPE!=""' "$VCF_IN" -Ov -o "$VCF_SV"
bcftools view -e 'INFO/SVTYPE!=""' "$VCF_IN" -Ov -o "$VCF_NONSV"

echo "[DONE] Completed splitting for ${SAMPLE}"


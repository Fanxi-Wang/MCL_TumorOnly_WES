#!/bin/bash
#SBATCH --job-name=vardict
#SBATCH --output=logs/vardict_%A.out
#SBATCH --error=logs/vardict_%A.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-35

# ==============================================================
# Tumor-only Somatic Variant Calling with VarDict
# Author: Fanxi Wang
# Description:
#   Runs VarDict in tumor-only mode for multiple WES samples
#   using Slurm array jobs.
# ==============================================================

# Step 1. Load dependencies
module load anaconda3
source ~/.bashrc
conda activate my_env   # contains VarDict, bgzip, tabix

# Step 2. Define directories
REF_DIR="/path/to/reference/hg38"
INPUT_DIR="/path/to/output_postalign"
OUT_DIR="/path/to/output_variantcalling/vardict"
REFERENCE="${REF_DIR}/Homo_sapiens_assembly38.fasta"
BED_FILE="${REF_DIR}/hg38_exome_targets.bed"

mkdir -p "$OUT_DIR" logs

# Step 3. Define 35 samples
SAMPLES=("Sample1" "Sample2" "Sample3" "Sample4" "Sample5"
         "Sample6" "Sample7" "Sample8" "Sample9" "Sample10"
         "Sample11" "Sample12" "Sample13" "Sample14" "Sample15"
         "Sample16" "Sample17" "Sample18" "Sample19" "Sample20"
         "Sample21" "Sample22" "Sample23" "Sample24" "Sample25"
         "Sample26" "Sample27" "Sample28" "Sample29" "Sample30"
         "Sample31" "Sample32" "Sample33" "Sample34" "Sample35")

# Step 4. Assign current sample based on SLURM_ARRAY_TASK_ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}
TUMOR_BAM="${INPUT_DIR}/recalibrated_${SAMPLE}.bam"
OUTPUT_VCF="${OUT_DIR}/vardict_${SAMPLE}_tumor.vcf"
AF_THR="0.01"

# Step 5. Run VarDict
echo "[INFO] Running VarDict for sample ${SAMPLE}..."

vardict-java \
  -G "$REFERENCE" \
  -f "$AF_THR" \
  -N "$SAMPLE" \
  -b "$TUMOR_BAM" \
  -c 1 -S 2 -E 3 \
  -Q 30 \
  -th 8 \
  "$BED_FILE" \
  | teststrandbias.R \
  | var2vcf_valid.pl -N "$SAMPLE" -E -f "$AF_THR" > "$OUTPUT_VCF"

# Step 6. Compress & Index
bgzip -c "$OUTPUT_VCF" > "${OUTPUT_VCF}.gz"
tabix -p vcf "${OUTPUT_VCF}.gz"

echo "[DONE] VarDict tumor-only calling completed for sample ${SAMPLE}"

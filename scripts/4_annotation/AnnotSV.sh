#!/bin/bash
#SBATCH --job-name=annotsv
#SBATCH --output=logs/annotsv_%A_%a.out
#SBATCH --error=logs/annotsv_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-35

# ==============================================================
# Structural Variant Annotation using AnnotSV
# Author: Fanxi Wang
# Description:
#   Annotates VarDict-derived SV VCF files using AnnotSV.
# ==============================================================

# Step 1. Load environment
module load anaconda3
source ~/.bashrc
conda activate annotsv_env  # contains AnnotSV installation

# Step 2. Define directories
ANNOTATION_DIR="$CONDA_PREFIX/share/AnnotSV/AnnotSV_annotations"   # automatic path from environment
INPUT_DIR="/path/to/output_variantintegration/split_vcfs"          # input SV-only VCFs
OUTPUT_DIR="/path/to/output_annotation/annotsv"                    # output TSVs
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
VCF_FILE="${INPUT_DIR}/${SAMPLE}_SV.vcf"
TMP_DIR="/scratch/$USER/annotsv_${SAMPLE}_$SLURM_JOB_ID"

echo "[INFO] Running AnnotSV for ${SAMPLE}..."
echo "       Input VCF:  $VCF_FILE"
echo "       Output Dir: $OUTPUT_DIR"

# Step 5. Run AnnotSV
mkdir -p "$TMP_DIR"
cd "$TMP_DIR"

AnnotSV \
  -annotationsDir "$ANNOTATION_DIR" \
  -SVinputFile "$VCF_FILE" \
  -genomeBuild GRCh38 \
  -outputFile "${SAMPLE}_SVs_annotated.tsv" \
  -outputDir "$TMP_DIR" \
  -overwrite 1

# Step 6. Copy results back to output directory
cp "${TMP_DIR}/${SAMPLE}_SVs_annotated.tsv" "$OUTPUT_DIR/"

# Optional: clean temporary workspace
rm -rf "$TMP_DIR"

echo "[DONE] AnnotSV annotation completed for ${SAMPLE}"

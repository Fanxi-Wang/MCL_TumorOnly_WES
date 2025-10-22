#!/bin/bash
#SBATCH --job-name=annovar_annotation
#SBATCH --output=logs/annovar_%A.out
#SBATCH --error=logs/annovar_%A.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-35

# ==============================================================
# Functional Annotation of SNV/INDEL Variants using ANNOVAR
# Author: Fanxi Wang
# Description:
#   Annotates merged Mutect2 + VarDict SNV VCFs using ANNOVAR
#   with multiple databases (refGene, dbNSFP, gnomAD, ClinVar, COSMIC).
# ==============================================================

# Step 1. Load environment
module load spack/latest
module load perl/5.40.0-b3brrq5  # required for ANNOVAR

# Step 2. Define directories
ANNOVAR_DIR="/path/to/tools/annovar"                    # ANNOVAR installation path
DB_DIR="${ANNOVAR_DIR}/humandb"                         # ANNOVAR database folder
INPUT_DIR="/path/to/output_variant_integration/merged"  # merged Mutect2 + VarDict SNV VCFs
OUTPUT_DIR="/path/to/output_annotation/annovar"         # annotation output folder
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
VCF_FILE="${INPUT_DIR}/${SAMPLE}_merged.vcf"
OUT_PREFIX="${OUTPUT_DIR}/${SAMPLE}_merged_annotated"

echo "[INFO] Running ANNOVAR for ${SAMPLE}..."
echo "       Input VCF:  $VCF_FILE"
echo "       Output Dir: $OUTPUT_DIR"

# Step 5. Run ANNOVAR annotation
cd "$ANNOVAR_DIR"

perl table_annovar.pl "$VCF_FILE" "$DB_DIR" \
  --buildver hg38 \
  --out "$OUT_PREFIX" \
  --remove \
  --protocol refGene,dbnsfp47a,gnomad41_exome,clinvar_20240917,avsnp151,cosmic102 \
  --operation g,f,f,f,f,f \
  --nastring . \
  --vcfinput \
  --polish

echo "[DONE] ANNOVAR annotation completed for ${SAMPLE}"

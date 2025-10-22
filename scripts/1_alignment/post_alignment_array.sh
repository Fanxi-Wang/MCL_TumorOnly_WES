#!/bin/bash
#SBATCH --job-name=post_alignment
#SBATCH --output=logs/post_alignment_%A.out
#SBATCH --error=logs/post_alignment_%A.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-35

# Step 1. Load required modules / environment
module load gatk/4.6.1
module load samtools
module load picard

# (Optional) Activate conda env if used
module load anaconda3
source ~/.bashrc
conda activate my_env

# Step 2. Define input / output directories
REFERENCE="/path/to/reference/Homo_sapiens_assembly38.fasta"
DBSNP_VCF="/path/to/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
INPUT_DIR="path/to/output_alignment"
OUTPUT_DIR="/path/to/output_postalign"  
REF_DIR="/path/to/reference"

# Step 3. Define 35 samples
SAMPLES=("Sample_01" "Sample_02" "Sample_03" "Sample_04" "Sample_05"
         "Sample_06" "Sample_07" "Sample_08" "Sample_09" "Sample_10"
         "Sample_11" "Sample_12" "Sample_13" "Sample_14" "Sample_15"
         "Sample_16" "Sample_17" "Sample_18" "Sample_19" "Sample_20"
         "Sample_21" "Sample_22" "Sample_23" "Sample_24" "Sample_25"
         "Sample_26" "Sample_27" "Sample_28" "Sample_29" "Sample_30"
         "Sample_31" "Sample_32" "Sample_33" "Sample_34" "Sample_35")

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

INPUT_BAM="${INPUT_DIR}/sorted_${SAMPLE}.bam"
DEDUP_BAM="${OUTPUT_DIR}/dedup_${SAMPLE}.bam"
RECAL_TABLE="${OUTPUT_DIR}/recal_data_${SAMPLE}.table"
RECAL_BAM="${OUTPUT_DIR}/recalibrated_${SAMPLE}.bam"

# Step 4. Generate required reference indices
echo "Preparing reference indices..."
samtools faidx "$REFERENCE"

if [ ! -f "${REF_DIR}/Homo_sapiens_assembly38.dict" ]; then
  gatk CreateSequenceDictionary -R "$REFERENCE" -O "${REF_DIR}/Homo_sapiens_assembly38.dict"
fi

# Step 5. Mark Duplicates
echo "Running MarkDuplicates for $SAMPLE..."
gatk MarkDuplicates \
  -I "$INPUT_BAM" \
  -O "$DEDUP_BAM" \
  -M "${OUTPUT_DIR}/dedup_metrics_${SAMPLE}.txt" \
  --REMOVE_DUPLICATES true

samtools index "$DEDUP_BAM"

# Step 6. Base Quality Score Recalibration (BQSR)
echo "Running BaseRecalibrator for $SAMPLE..."
gatk BaseRecalibrator \
  -R "$REFERENCE" \
  -I "$DEDUP_BAM" \
  --known-sites "$DBSNP_VCF" \
  -O "$RECAL_TABLE"

gatk ApplyBQSR \
  -R "$REFERENCE" \
  -I "$DEDUP_BAM" \
  --bqsr-recal-file "$RECAL_TABLE" \
  -O "$RECAL_BAM"

samtools index "$RECAL_BAM"

# Step 7. Completion message
echo "Post-alignment processing completed for $SAMPLE!"
echo "Output files saved to: $OUTPUT_DIR"


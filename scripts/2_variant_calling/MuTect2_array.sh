#!/bin/bash
#SBATCH --job-name=mutect2
#SBATCH --output=logs/mutect2_%A.out
#SBATCH --error=logs/mutect2_%A.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-35   

# ==============================================================
# Tumor-only Somatic Variant Calling with GATK Mutect2
# Author: Fanxi Wang
# Description:
#   Runs GATK Mutect2 in tumor-only mode for multiple samples
#   using Slurm array jobs.
# ==============================================================

# Step 1. Load dependencies
module load gatk/4.6.1
module load samtools
module load anaconda3
source ~/.bashrc
conda activate my_env 


# Step 2. Define directories
REF_DIR="/path/to/reference/hg38"
INPUT_DIR="/path/to/output_postalign"   
OUT_DIR="/path/to/output_variantcalling/mutect"
REFERENCE="${REF_DIR}/Homo_sapiens_assembly38.fasta"
GERMLINE_RESOURCE="${REF_DIR}/af-only-gnomad.hg38.vcf.gz"

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
INPUT_BAM="${INPUT_DIR}/recalibrated_${SAMPLE}.bam"
OUTPUT_VCF="${OUT_DIR}/mutect2_${SAMPLE}_tumor.vcf.gz"

echo ">>> [$(date)] Starting Mutect2 variant calling for ${SAMPLE}"

# Step 5. Run Mutect2
gatk Mutect2 \
  -R "$REFERENCE" \
  -I "$INPUT_BAM" \
  -tumor "$SAMPLE" \
  --germline-resource "$GERMLINE_RESOURCE" \
  -O "$OUTPUT_VCF"

# Step 6. Filter variants
gatk FilterMutectCalls \
  -R "$REFERENCE" \
  -V "$OUTPUT_VCF" \
  -O "${OUT_DIR}/mutect2_${SAMPLE}_tumor_filtered.vcf.gz"

echo ">>> [$(date)] Mutect2 tumor-only variant calling completed for ${SAMPLE}"

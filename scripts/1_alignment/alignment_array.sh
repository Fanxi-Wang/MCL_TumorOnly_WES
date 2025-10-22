#!/bin/bash
#SBATCH --job-name=bwa_align
#SBATCH --output=logs/bwa_align_%A.out
#SBATCH --error=logs/bwa_align_%A.err
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-35

# Step 1. Load required modules
module load bwa
module load samtools

# Step 2. Define directories
REFERENCE="/path/to/reference/Homo_sapiens_assembly38.fasta"
FASTQ_DIR="/path/to/fastq"
OUTPUT_DIR="path/to/output_alignment"
SCRATCH_DIR="/scratch/$USER/$SLURM_JOB_ID"

# Step 3. Define 35 samples
SAMPLES=("Sample_01" "Sample_02" "Sample_03" "Sample_04" "Sample_05"
         "Sample_06" "Sample_07" "Sample_08" "Sample_09" "Sample_10"
         "Sample_11" "Sample_12" "Sample_13" "Sample_14" "Sample_15"
         "Sample_16" "Sample_17" "Sample_18" "Sample_19" "Sample_20"
         "Sample_21" "Sample_22" "Sample_23" "Sample_24" "Sample_25"
         "Sample_26" "Sample_27" "Sample_28" "Sample_29" "Sample_30"
         "Sample_31" "Sample_32" "Sample_33" "Sample_34" "Sample_35")

# Get the sample name based on SLURM array ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

FASTQ_R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
FASTQ_R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

# Step 4. Run alignment
echo "Starting BWA alignment for ${SAMPLE}"

mkdir -p $SCRATCH_DIR
mkdir -p $OUTPUT_DIR
cp $FASTQ_R1 $FASTQ_R2 $REFERENCE* $SCRATCH_DIR
cd $SCRATCH_DIR

bwa mem -t 8 Homo_sapiens_assembly38.fasta $FASTQ_R1 $FASTQ_R2 | \
  samtools view -Sb - > aligned_${SAMPLE}.bam

samtools sort -o sorted_${SAMPLE}.bam aligned_${SAMPLE}.bam
samtools index sorted_${SAMPLE}.bam

cp sorted_${SAMPLE}.bam* $OUTPUT_DIR/

rm -rf $SCRATCH_DIR

echo "Alignment completed for ${SAMPLE}! Output in $OUTPUT_DIR"


#!/bin/bash
# ==============================================
# Script: submit_all.sh
# Description:
#   Submit all Slurm jobs for each pipeline stage
#   (Alignment → Variant Calling → Variant Integration → Annotation → Visualization)
# Author: Fanxi Wang
# ==============================================

# === 1. Alignment stage ===
echo "Submitting alignment jobs..."
sbatch ../1_alignment/alignment_array.sh
sbatch ../1_alignment/post_alignment_array.sh

# === 2. Variant Calling stage ===
echo "Submitting variant calling jobs..."
sbatch ../2_variant_calling/MuTect2_array.sh
sbatch ../2_variant_calling/VarDict_array.sh

# === 3. Variant Integration stage ===
echo "Submitting variant integration jobs..."
sbatch ../3_variant_integration/vcf_SV_split.sh
sbatch ../3_variant_integration/merge_vcfs.sh


# === 4. Annotation stage ===
echo "Submitting annotation jobs..."
sbatch ../4_annotation/ANNOVAR_annotation.sh
sbatch ../4_annotation/AnnotSV.sh

# === 5. Visualization (if applicable) ===
echo "Visualization step can be run locally (e.g., in RStudio)."


# MCL_TumorOnly_WES

This repository contains the **tumor-only whole exome sequencing (WES) analysis pipeline** for **mantle cell lymphoma (MCL)** samples.  
The workflow integrates alignment, variant calling, variant integration, annotation, and visualization, implemented on the HPC cluster.

---

## Workflow Overview

1. **Alignment**
   - `alignment_array.sh`: Performs BWA-MEM alignment of paired-end FASTQ reads.
   - `post_alignment_array.sh`: Marks duplicates and applies base quality score recalibration (BQSR).

2. **Variant Calling**
   - `MuTect2_array.sh`: Runs GATK Mutect2 for somatic SNV/indel calling.
   - `VarDict_array.sh`: Runs VarDict for complementary variant detection.

3. **Variant Integration**
   - `vcf_SV_split.sh`: Splits structural variants (SVs) and SNVs.
   - `merge_vcfs.py` & `merge_vcfs.sh`: Merges Mutect2 and VarDict VCFs for SNV integration.

4. **Annotation**
   - `ANNOVAR_annotation.sh`: Functional annotation of SNVs/indels with ANNOVAR.
   - `AnnotSV.sh`: Annotation of structural variants (SVs) with AnnotSV.

5. **Visualization**
   - `SNV_Annotation.R`: Post-processed SNV filtering and statistics.
   - `SVs_Annotation.R`: Aggregates AnnotSV results and computes quality metrics.
   - `Downstream_Oncoprint.R`: Generates mutation heatmaps (oncoprints) for key genes.

---

## Environment Setup

Conda environment configuration files are provided in the [`env/`](./env) folder:

| Environment | Description |
|--------------|-------------|
| `my_env.yml` | Core pipeline: GATK, Samtools, Mutect2, VarDict, bcftools |
| `mc3_env.yml` | Legacy VCF merging (Python 2 compatibility) |
| `annotsv_env.yml` | Structural variant annotation |

To recreate an environment:
```bash
conda env create -f env/my_env.yml

```
## HPC Execution

All .sh scripts are designed for Slurm-based HPC scheduling.  
Example job submission:

```bash
sbatch scripts/2_variant_calling/MuTect2_array.sh
```

## Output Summary

| Step | Description | Example Outputs |
|------|--------------|------------------|
| **1. Alignment** | Raw alignment and recalibration | `sorted_*.bam`, `recalibrated_*.bam` |
| **2. Variant Calling** | Somatic SNV and indel calling using Mutect2 and VarDict | `*_mutect2_tumor.vcf.gz`, `*_vardict_tumor.vcf` |
| **3. Variant Integration** | Merge Mutect2 + VarDict results; split SNV/SV | `*_merged.vcf`, `*_nonSV.vcf`, `*_SV.vcf` |
| **4. Annotation** | Functional and structural variant annotation | `*_hg38_multianno.txt`, `*_SVs_annotated.tsv` |
| **5. Visualization** | Oncoprint plots and variant summary tables | HTML reports, annotated figures |

---

**Author:** Fanxi Wang  
**Institution:** Weill Cornell Medicine  
**Lab:** Elemento Lab, Englander Institute for Precision Medicine

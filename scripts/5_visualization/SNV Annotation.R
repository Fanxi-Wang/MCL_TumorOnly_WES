library(maftools)
library(dplyr)

# Step 1: Set working directory
setwd("/path/to/output_annotation/annovar")

# Step 2: Get list of all 35 annotation files
anno_files <- list.files(pattern = "*.hg38_multianno.txt$")

# Step 3: Prepare empty list to store processed data
updated_list <- list()

# Step 4: Loop through each file
for (file in anno_files) {
  # Extract sample_id: string between "qc_" and ".hg38"
  sample_id <- sub(".*qc_", "", file)  # remove everything before qc_
  sample_id <- sub(".hg38_multianno.txt", "", sample_id)  # remove suffix
  
  # Read file
  this_df <- read.delim(file, sep = "\t", stringsAsFactors = FALSE)
  
  # Add Tumor_Sample_Barcode column
  this_df$Tumor_Sample_Barcode <- sample_id
  
  # Save to list (do not overwrite files!)
  updated_list[[sample_id]] <- this_df
}

# Step 5: Combine all sample dataframes into one big table
combined_df <- bind_rows(updated_list)

# Step 6: Write this combined file to use as input for annovarToMaf
write.table(combined_df, file = "combined_with_barcode.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Step 7: Convert to MAF
maf <- annovarToMaf(
  annovar = "combined_with_barcode.txt",
  Center = "WCM",
  refBuild = "hg38",
  table = "refGene",
  MAFobj = FALSE,
  tsbCol = "Tumor_Sample_Barcode"
)

# Step 8: Save MAF table
write.table(maf, file = "new_final_merged.maf", sep = "\t", quote = FALSE, row.names = FALSE) 

colnames(maf)
 
combined_df <- read.delim(
  "combined_with_barcode.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

merged_maf <- read.delim("new_final_merged.maf", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Step 0: Load required libraries
library(dplyr)
library(stringr)

# Step 1: Extract AF from Otherinfo13
fields <- do.call(rbind, strsplit(maf$Otherinfo13, ":", fixed = TRUE))
maf$GT <- fields[, 1]
maf$DP <- as.numeric(fields[, 2])
maf$AD <- fields[, 3]
maf$AF <- as.numeric(fields[, 4])


# Step 2: Extract COSMIC occurrence count from cosmic102 
# Count the number of occurrences by counting bracket pairs (each sample)
maf$COSMIC_CNT <- str_count(maf$cosmic102, "\\([^)]+\\)")

# Step 3: Build GnomAD + COSMIC filter logic 
maf$gnomAD_high <- ifelse(is.na(maf$gnomad41_exome_AF), FALSE, maf$gnomad41_exome_AF >= 0.01)
maf$cosmic_keep <- maf$COSMIC_CNT >= 2
maf$gnomAD_cosmic_pass <- ifelse(maf$gnomAD_high & !maf$cosmic_keep, FALSE, TRUE)


# Step 4: Apply full filtering
filtered_maf <- maf %>%
  filter(
    # GnomAD high-frequency variants are removed unless COSMIC count is high
    gnomAD_cosmic_pass,

    # Sequencing depth filter
    DP >= 10,

    # Allele frequency filter
    AF >= 0.02,

    # Functional classification filter
    !(Variant_Classification %in% c(
      "Intron", "5'Flank", "3'Flank", "IGR",
      "Targeted_Region", "RNA", "Silent"
    )),

    # Deleteriousness filter: at least one satisfied
    REVEL_score >= 0.5 |
    MetaSVM_pred == "D" |
    CADD_phred >= 20
  )


benign_labels <- c(
  "Benign", "Likely_benign", 
  "Benign/Likely_benign", 
  "Benign/Likely_benign|other"
)

filtered_maf <- filtered_maf %>%
  filter(!(CLNSIG %in% benign_labels))


# Step 6: Annotate caller type
filtered_maf$caller_type <- case_when(
  filtered_maf$Otherinfo11 == "CENTERS=VarDict" ~ "vardict-only",
  filtered_maf$Otherinfo11 == "CENTERS=Mutect2" ~ "mutect2-only",
  str_detect(filtered_maf$Otherinfo11, "CENTERS=Mutect2\\|VarDict") ~ "merged",
  TRUE ~ "unknown"
)

# Step 7: Output summary
cat("total variants after filtering: ", nrow(filtered_maf), "\n")
print(table(filtered_maf$caller_type))

# Step 8: Post-hoc ClinSig analysis 
cat("ClinSig distribution (not used for filtering):\n")
print(table(filtered_maf$CLIN_SIG))

# Step 9: Save result
write.table(filtered_maf, "filtered_maf.tsv", sep = "\t", quote = FALSE, row.names = FALSE)



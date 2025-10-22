library(dplyr)
library(OncoScience)

# ===============================
# Step 1. Load merged annotation TSV
# ===============================
merged_df <- read.delim("merged_all.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ===============================
# Step 2. Parse the "PRIMARY" field
# The field structure looks like "GT:DP:VD:AD", e.g., "0/1:251:3:248,3"
# ===============================

# Split each PRIMARY field by ":"
fields <- strsplit(merged_df$PRIMARY, ":")

# Extract AD (4th element, e.g., "248,3")
merged_df$AD <- sapply(fields, function(x) x[4])

# Extract VD (3rd element, e.g., "3")
merged_df$VD <- sapply(fields, function(x) x[3])

# ===============================
# Step 3. Calculate DP (depth) from AD
# AD field contains two comma-separated numbers, e.g. "248,3"
# ===============================
merged_df$DP <- sapply(merged_df$AD, function(ad_str) {
  ad_values <- as.numeric(unlist(strsplit(ad_str, ",")))
  sum(ad_values)
})

# Convert VD to numeric
merged_df$VD <- as.numeric(merged_df$VD)

# ===============================
# Step 4. Compute allele frequency (AF = VD / DP)
# ===============================
merged_df$AF <- merged_df$VD / merged_df$DP

# ===============================
# Step 5. Apply filtering criteria
# Keep variants with sufficient depth (DP ≥ 10) and allele frequency (AF ≥ 0.02)
# ===============================
merged_df_filtered <- merged_df %>%
  filter(AF >= 0.02 & DP >= 10)

# ===============================
# Step 6. Keep only rows with annotation mode = "split"
# ===============================
merged_df_filtered_split <- merged_df_filtered %>%
  filter(Annotation_mode == "split")

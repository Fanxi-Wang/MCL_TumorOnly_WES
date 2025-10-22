library(dplyr)

# Set the folder path
folder_path <- "/path/to/output_annotation/annotsv"

# Get all TSV filenames 
file_list <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)

# Initialize a list to store data frames
merged_list <- list()

# Loop over each file
for (i in seq_along(file_list)) {
  file_path <- file_list[i]
  file_name <- basename(file_path)

  # Extract sample ID
  sample_id <- sub("_SVs_annotated\\.tsv$", "", file_name)

  # Read the file 
  df <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Overwrite the Samples_ID column with the correct sample name
  df$Samples_ID <- sample_id

  # Append to the list
  merged_list[[i]] <- df
}

# Combine all data frames into one big data frame
merged_df <- do.call(rbind, merged_list)

# Save to a merged output file
write.table(merged_df,
            file = file.path(folder_path, "merged_all.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)



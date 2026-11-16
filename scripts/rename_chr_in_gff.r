# read gff file at /mnt/hc-storage/users/esharma/rice/ms/test_por/Input_fasta_and_gff/oryza_sativa_nipponbaremerged.gff 
# in the first column, replace the values like 1, 2, 3 etc with chr01, chr02, chr03 etc.. upto chr12


#!/usr/bin/env Rscript

# Read the GFF file
gff_file <- "/mnt/hc-storage/users/esharma/rice/ms/test_por/Input_fasta_and_gff/oryza_sativa_nipponbaremerged.gff"
output_file <- "/mnt/hc-storage/users/esharma/rice/ms/oryza_sativa_nipponbaremerged_modified.gff"

# Read the file
gff_data <- readLines(gff_file)

# Process each line
modified_data <- sapply(gff_data, function(line) {
  # Skip comment lines
  if (grepl("^#", line)) {
    return(line)
  }
  
  # Split the line by tabs
  fields <- strsplit(line, "\t")[[1]]
  
  # Check if the first field is a number between 1 and 12
  if (length(fields) > 0 && grepl("^[0-9]+$", fields[1])) {
    chr_num <- as.integer(fields[1])
    if (chr_num >= 1 && chr_num <= 12) {
      # Replace with chr01, chr02, etc.
      fields[1] <- sprintf("chr%02d", chr_num)
    }
  }
  
  # Reconstruct the line
  return(paste(fields, collapse = "\t"))
})

# Write the modified data to output file
writeLines(modified_data, output_file)

cat("Successfully processed GFF file.\n")
cat("Output written to:", output_file, "\n")

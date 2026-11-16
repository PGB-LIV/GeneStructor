#!/usr/bin/env Rscript
# Optimized script to process 100k+ protein structures (.CIF files) using multiple cores
# Extracts pLDDT metrics 
# Handles naming mismatch: CIF files are shortened (e.g., G69055.cif)

# Load required libraries
library(data.table)
library(tidyverse)
library(bio3d)
library(parallel)

# ============================================================================
# CONFIGURATION
# ============================================================================

setwd("./")

base_dir_cif <- "./consolidated_nip/cif"

# Detect cores (leave 2 free for system)
n_cores <- max(1, detectCores() - 2)
cat(sprintf("Using %d cores for parallel processing\n\n", n_cores))

# ============================================================================
# PROCESS CIF FILES
# ============================================================================

cat("========================================\n")
cat("STEP: Processing CIF files\n")
cat("========================================\n")

pattern_cif <- "\\.cif$"
cif_files <- list.files(path = base_dir_cif, pattern = pattern_cif, 
                         recursive = FALSE, full.names = TRUE)

cat(sprintf("Found %d CIF files\n", length(cif_files)))

if (length(cif_files) == 0) {
  stop("No CIF files found in: ", base_dir_cif)
}

# Process CIF files in parallel
cat("Parsing CIF files and calculating pLDDT statistics (parallel processing)...\n")
cat("This may take several minutes for 100k+ files...\n")

system.time({
  cif_records <- mclapply(cif_files, function(f) {
    # Extract ID from filename
    short_id <- sub(pattern_cif, "", basename(f))
    
    # Read and parse CIF file
    cif_data <- tryCatch({
      cif_obj <- read.cif(f)
      as.data.frame(cif_obj$atom)
    }, error = function(e) {
      warning("Failed to read CIF: ", basename(f), " - ", e$message)
      return(NULL)
    })
    
    if (is.null(cif_data)) {
      return(data.frame(
        short_id = short_id,
        mean_plddt_cif = NA_real_,
        total_residues = NA_integer_,
        plddt_50 = NA_integer_,
        plddt_60 = NA_integer_,
        plddt_70 = NA_integer_,
        plddt_80 = NA_integer_,
        plddt_90 = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }
    
    # Calculate per-residue pLDDT (mean b-factor per residue)
    cif_data$b <- as.numeric(cif_data$b)
    total_residues <- length(unique(cif_data$resno))
    
    # Group by residue and calculate statistics
    residue_stats <- cif_data %>%
      group_by(resno) %>%
      summarise(plddt = mean(b, na.rm = TRUE), .groups = 'drop') %>%
      summarise(
        mean_plddt_cif = mean(plddt, na.rm = TRUE),
        plddt_50 = sum(plddt > 50),
        plddt_60 = sum(plddt > 60),
        plddt_70 = sum(plddt > 70),
        plddt_80 = sum(plddt > 80),
        plddt_90 = sum(plddt > 90),
        .groups = 'drop'
      )
    
    data.frame(
      short_id = short_id,
      mean_plddt_cif = residue_stats$mean_plddt_cif,
      total_residues = total_residues,
      plddt_50 = residue_stats$plddt_50,
      plddt_60 = residue_stats$plddt_60,
      plddt_70 = residue_stats$plddt_70,
      plddt_80 = residue_stats$plddt_80,
      plddt_90 = residue_stats$plddt_90,
      stringsAsFactors = FALSE
    )
  }, mc.cores = n_cores)
})

# Combine results
cif_df <- rbindlist(cif_records, fill = TRUE)

cat(sprintf("Successfully parsed %d CIF files\n", nrow(cif_df)))
cat(sprintf("Preview of CIF data:\n"))
print(head(cif_df, 5))
cat("\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("========================================\n")
cat("STEP 4: Saving results\n")
cat("========================================\n")

# Save as RDS for faster loading in R
output_rds <- "./protenix_combined_stats_all_nip.rds"
saveRDS(cif_df, output_rds)
cat(sprintf("Saved RDS to: %s\n", output_rds))

# Save separate CIF results as backups

fwrite(cif_df, "./RE_protenix_cif_stats.csv")
cat("Saved CIF results \n")

cat("\n========================================\n")
cat("PROCESSING COMPLETE!\n")
cat("========================================\n")
cat(sprintf("Total time: %.2f minutes\n", proc.time()[3] / 60))

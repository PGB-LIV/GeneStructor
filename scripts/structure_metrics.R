#!/usr/bin/env Rscript
# pLDDT metrics from predicted protein structures (.pdb from AlphaFold2/ColabFold, .cif from Protenix).
# pLDDT is read from the B-factor column and averaged per residue.
#
# Usage:
#   Rscript scripts/structure_metrics.R --input-dir structures/ --out structure_metrics.csv [--cores 4]
#
# Output columns: idx (file name without extension), mean_plddt, total_residues and the number of residues with
# pLDDT above 50, 60, 70, 80 and 90 (plddt_50 ... plddt_90). GeneStructor uses plddt_80 for AlphaFold2 and
# plddt_70 for Protenix-mini.

suppressPackageStartupMessages({
  library(data.table)
  library(bio3d)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
if (length(args) == 0 || "--help" %in% args) {
  cat("See the header of this script or the README for usage.\n")
  quit(status = 0)
}

input_dir <- get_arg("input-dir")
out_file <- get_arg("out")
if (is.null(input_dir) || is.null(out_file)) stop("--input-dir and --out are required")
# forked parallelism is not available on Windows
n_cores <- if (.Platform$OS.type == "windows") 1L else as.integer(get_arg("cores", "1"))

files <- list.files(input_dir, pattern = "\\.(pdb|cif)$", full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0) stop("No .pdb or .cif files found in ", input_dir)
cat(sprintf("Processing %d structure files with %d core(s)\n", length(files), n_cores))

structure_stats <- function(f) {
  idx <- sub("\\.(pdb|cif)$", "", basename(f), ignore.case = TRUE)
  atoms <- tryCatch({
    s <- if (grepl("\\.cif$", f, ignore.case = TRUE)) read.cif(f) else read.pdb(f, verbose = FALSE)
    as.data.frame(s$atom)
  }, error = function(e) NULL)
  if (is.null(atoms) || nrow(atoms) == 0) {
    return(data.table(idx = idx, mean_plddt = NA_real_, total_residues = NA_integer_,
                      plddt_50 = NA_integer_, plddt_60 = NA_integer_, plddt_70 = NA_integer_,
                      plddt_80 = NA_integer_, plddt_90 = NA_integer_))
  }
  # one pLDDT value per residue (chain + residue number)
  plddt <- tapply(as.numeric(atoms$b), paste(atoms$chain, atoms$resno), mean, na.rm = TRUE)
  data.table(idx = idx, mean_plddt = mean(plddt), total_residues = length(plddt),
             plddt_50 = sum(plddt > 50), plddt_60 = sum(plddt > 60), plddt_70 = sum(plddt > 70),
             plddt_80 = sum(plddt > 80), plddt_90 = sum(plddt > 90))
}

stats <- rbindlist(mclapply(files, structure_stats, mc.cores = n_cores))
fwrite(stats, out_file)
cat(sprintf("Wrote pLDDT metrics for %d structures to %s (%d could not be read)\n",
            nrow(stats), out_file, sum(is.na(stats$total_residues))))

#!/usr/bin/env Rscript
# Plot the ranked gene models of one gene: exon structure from a GFF3, coloured by source
# annotation, with the GeneStructor rank and score of each gene model.
#
# Usage:
#   Rscript scripts/plot_ranked_gene.R \
#     --ranked ranked.tsv \
#     --gff example/Os06g0113150.gff3 \
#     --gene Os06g0113150 \
#     --out Os06g0113150.png
#
# Arguments:
#   --ranked   output of rank_isoforms.R
#   --gff      GFF3 with CDS (and optionally UTR) features of the gene models; features are matched to
#              gene models by their Parent attribute
#   --gene     parent gene to plot (column Gene in ranked.tsv)
#   --out      output image (.png, .pdf, ...); default <gene>.png
#   --width    plot width in inches (default 10)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

ranked_file <- get_arg("ranked")
gff_file    <- get_arg("gff")
gene_id     <- get_arg("gene")
if (is.null(ranked_file) || is.null(gff_file) || is.null(gene_id)) {
  stop("--ranked, --gff and --gene are required")
}
out_file <- get_arg("out", paste0(str_replace_all(gene_id, "[^[:alnum:]_.-]", "_"), ".png"))
width    <- as.numeric(get_arg("width", 10))

# plot_gene_with_ml_rank() lives next to this script
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(sub("^--file=", "", script_arg[1])) else "scripts"
source(file.path(script_dir, "plot_single_gene_with_labels.R"))

gff_lines <- grep("^#", readLines(gff_file), value = TRUE, invert = TRUE)
gff <- fread(text = gff_lines, sep = "\t", header = FALSE, quote = "",
             col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

# one row per gene model; gene models with identical metrics share a row in ranked.tsv
gene_transcripts <- fread(ranked_file) %>%
  filter(Gene == gene_id) %>%
  select(transcript_id_collapsed, ML_rank, median_pred_CT) %>%
  separate_rows(transcript_id_collapsed, sep = ";|,") %>%
  mutate(
    transcript_id  = str_trim(transcript_id_collapsed),
    parent_gene_id = gene_id,
    .pred_CT       = median_pred_CT,
    database = case_when(
      str_detect(transcript_id, "^Os(?!Nip)") ~ "RAPDB",
      str_starts(transcript_id, "LOC")        ~ "RGAP",
      str_starts(transcript_id, "OsNip")      ~ "Gramene",
      TRUE                                    ~ "Other"
    )
  ) %>%
  select(-transcript_id_collapsed) %>%
  filter(transcript_id != "") %>%
  distinct()

if (nrow(gene_transcripts) == 0) stop("Gene ", gene_id, " not found in ", ranked_file)

plot_gene_with_ml_rank(gff, gene_transcripts, output_file = out_file, width = width)
cat("Plot written to", out_file, "\n")

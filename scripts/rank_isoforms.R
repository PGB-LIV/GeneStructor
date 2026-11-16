#!/usr/bin/env Rscript
# Rank the gene models (protein isoforms) at each gene locus with the GeneStructor ensemble of 25 logistic regression models.
#
# Usage:
#   Rscript scripts/rank_isoforms.R --metrics metrics.tsv --source alphafold --out ranked.tsv [--models models]
#
#   --metrics  Table from build_metrics_table.R, or any table with columns Gene, transcript_id, bits.max,
#              total_domain_length, in_frame_score and plddt_80 (AlphaFold2) or plddt_70 (Protenix-mini)
#   --source   alphafold (uses plddt_80) or protenix (uses plddt_70 in place of plddt_80)
#   --models   Folder with model_1.rds ... model_25.rds (default: models/ in this repository)
#   --out      Output table (tab-separated)
#   --digits   Decimal places used when ranking probabilities (default 2)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidymodels)
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

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
metrics_file <- get_arg("metrics")
out_file <- get_arg("out")
source_type <- tolower(get_arg("source", "alphafold"))
models_dir <- get_arg("models", file.path(script_dir, "..", "models"))
digits <- as.integer(get_arg("digits", "2"))
if (is.null(metrics_file) || is.null(out_file)) stop("--metrics and --out are required")
if (!source_type %in% c("alphafold", "protenix")) stop("--source must be 'alphafold' or 'protenix'")

plddt_col <- if (source_type == "protenix") "plddt_70" else "plddt_80"
predictors <- c(plddt_col, "bits.max", "total_domain_length", "in_frame_score")

metrics <- fread(metrics_file, integer64 = "double")
missing_cols <- setdiff(c("Gene", "transcript_id", predictors), names(metrics))
if (length(missing_cols) > 0) stop("Missing columns in metrics table: ", paste(missing_cols, collapse = ", "))

# A missing metric means no structure, hit, domain or score: treat it as 0, as in the manuscript
metrics <- metrics %>%
  mutate(across(all_of(predictors), ~ ifelse(is.na(.x), 0, as.numeric(.x)))) %>%
  mutate(bits.max = pmax(bits.max, 0))

# Gene models with identical metrics within a gene (usually identical proteins) are scored once
metrics <- metrics %>%
  group_by(Gene, across(all_of(predictors))) %>%
  mutate(transcript_id_collapsed = paste(transcript_id, collapse = ",")) %>%
  slice(1) %>%
  ungroup()

# the models were trained with plddt_80; predict_with_glm_ensemble() copies plddt_70 into it for Protenix
if (!"plddt_80" %in% names(metrics)) metrics$plddt_80 <- NA_real_
if (!"plddt_70" %in% names(metrics)) metrics$plddt_70 <- NA_real_

source(file.path(script_dir, "predict_with_glm_ensemble.r"))
n_models <- length(list.files(models_dir, pattern = "^model_[0-9]+\\.rds$"))
if (n_models == 0) stop("No model_*.rds files found in ", models_dir)

predictions <- predict_with_glm_ensemble(new_data = metrics, models_dir = models_dir,
                                         n_models = n_models, source = source_type)

ranked <- predictions %>%
  group_by(Gene) %>%
  mutate(
    ML_rank = dense_rank(desc(round(median_pred_CT, digits))),
    prob_diff_to_rank_1 = round(round(max(median_pred_CT), digits) - round(median_pred_CT, digits), digits)
  ) %>%
  ungroup() %>%
  left_join(metrics %>% select(Gene, transcript_id, transcript_id_collapsed, all_of(predictors)),
            by = c("Gene", "transcript_id")) %>%
  select(Gene, ML_rank, transcript_id, transcript_id_collapsed, median_pred_CT, prob_diff_to_rank_1,
         sd_pred_CT, n_models, all_of(predictors)) %>%
  arrange(Gene, ML_rank, desc(median_pred_CT))

fwrite(ranked, out_file, sep = "\t")

multi <- ranked %>% group_by(Gene) %>% filter(n() > 1) %>% ungroup()
gap <- multi %>% group_by(Gene) %>% filter(any(ML_rank == 2)) %>%
  summarise(diff = min(prob_diff_to_rank_1[ML_rank == 2]), .groups = "drop")
cat(sprintf("Ranked %d distinct gene models in %d genes (%d genes with more than one model)\n",
            nrow(ranked), n_distinct(ranked$Gene), n_distinct(multi$Gene)))
cat(sprintf("Genes where rank 1 beats rank 2 by >= 0.1: %d of %d\n", sum(gap$diff >= 0.1), nrow(gap)))
cat("Results written to", out_file, "\n")

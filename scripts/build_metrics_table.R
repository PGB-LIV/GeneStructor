#!/usr/bin/env Rscript
# Build the GeneStructor metrics table: one row per transcript with the predictors used by the models
# (pLDDT residue counts, maximum Foldseek bitscore, total InterPro/Pfam domain length, PSAURON in-frame score).
#
# Usage:
#   Rscript scripts/build_metrics_table.R \
#     --ids transcript_gene.tab \
#     --structures structure_metrics.csv \
#     --foldseek-pdb fs_pdb.m8 --foldseek-afsp fs_afsp.m8 \
#     --interpro interproscan.tsv \
#     --psauron psauron_scores.csv \
#     --out metrics.tsv
#
# Required:
#   --ids            Tab-separated table with columns ID (transcript) and Parent (gene), e.g. from Simplified_GFF.rmd
#   --out            Output table (tab-separated)
# Optional (missing metrics are set to 0 when ranking):
#   --structure-map  Two columns, transcript_id and structure_id, when transcripts with identical proteins share one
#                    structure. Without it, structures are matched to transcripts by ID.
#   --structures     Output of process_colabfold_all.py (column idx) or process_cif_for_pLDDT.R (column short_id)
#   --foldseek-pdb   Foldseek easy-search output (m8) against PDB, run with the --format-output columns given in the README
#   --foldseek-afsp  Foldseek easy-search output (m8) against AlphaFold/Swiss-Prot
#   --foldseek-csv   Alternative to the two m8 files: output of process_fs_protenix_all_nip.py
#   --exclude-taxon  Drop AlphaFold/Swiss-Prot hits to this taxon name (e.g. "Oryza sativa Japonica Group")
#   --interpro       InterProScan TSV output (-f TSV)
#   --psauron        PSAURON output (CSV with description and in-frame_score columns)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
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
ids_file <- get_arg("ids")
out_file <- get_arg("out")
if (is.null(ids_file) || is.null(out_file)) stop("--ids and --out are required")

# strip structure file extensions and predictor suffixes so IDs match transcript IDs, e.g.
# Os01t0100100-01.pdb, Os01t0100100-01_scores_rank_001_alphafold2_ptm_model_1_seed_000 (ColabFold),
# Os01t0100100-01_sample_0 (Protenix) -> Os01t0100100-01
clean_id <- function(x) {
  x <- sub("\\.(cif|pdb|json)(\\.gz)?$", "", basename(x), ignore.case = TRUE)
  x <- sub("_(scores|unrelaxed|relaxed)_rank_[0-9]+.*$", "", x)
  sub("_sample_[0-9]+$", "", x)
}

# ---- transcript to gene ----
metrics <- fread(ids_file) %>%
  transmute(transcript_id = ID, Gene = Parent) %>%
  distinct()

map_file <- get_arg("structure-map")
if (!is.null(map_file)) {
  structure_map <- fread(map_file) %>% select(transcript_id = 1, structure_id = 2) %>% distinct(transcript_id, .keep_all = TRUE)
  metrics <- metrics %>% left_join(structure_map, by = "transcript_id")
} else {
  metrics <- metrics %>% mutate(structure_id = transcript_id)
}

# ---- structure metrics (pLDDT counts) ----
structures_file <- get_arg("structures")
if (!is.null(structures_file)) {
  structures <- fread(structures_file)
  if ("short_id" %in% names(structures)) structures <- rename(structures, idx = short_id)
  if ("mean_plddt_cif" %in% names(structures)) structures <- rename(structures, mean_plddt = mean_plddt_cif)
  # ColabFold writes several ranked models per protein: keep the top-ranked one, preferring its JSON (it carries pTM)
  if ("file_name" %in% names(structures) && any(grepl("_rank_[0-9]+", structures$file_name))) {
    structures <- structures %>% filter(!grepl("_rank_[0-9]+", file_name) | grepl("_rank_0*1([^0-9]|$)", file_name))
  }
  if ("file_type" %in% names(structures)) structures <- structures %>% arrange(file_type != "json")
  structures <- structures %>%
    mutate(structure_id = clean_id(idx)) %>%
    select(structure_id, any_of(c("ptm", "mean_plddt", "total_residues", "plddt_70", "plddt_80"))) %>%
    distinct(structure_id, .keep_all = TRUE)
  metrics <- metrics %>% left_join(structures, by = "structure_id")
}

# ---- Foldseek: maximum bitscore of the top hit (lowest E-value) against PDB and AlphaFold/Swiss-Prot ----
fs_cols <- c("query", "target", "alntmscore", "lddt", "fident", "alnlen", "mismatch", "gapopen",
             "qstart", "qend", "tstart", "tend", "evalue", "bits", "taxid", "taxname", "taxlineage")

top_hit_bits <- function(m8_file, exclude_taxon = NULL) {
  hits <- fread(m8_file, header = FALSE, sep = "\t", quote = "", integer64 = "double", fill = TRUE)
  setnames(hits, fs_cols[seq_len(ncol(hits))])
  hits <- hits %>% filter(evalue < 1e-5)
  if (!is.null(exclude_taxon) && "taxname" %in% names(hits)) hits <- hits %>% filter(taxname != exclude_taxon)
  hits %>%
    mutate(structure_id = clean_id(query)) %>%
    group_by(structure_id) %>%
    slice_min(evalue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(structure_id, bits)
}

exclude_taxon <- get_arg("exclude-taxon")
fs_pdb <- get_arg("foldseek-pdb")
fs_afsp <- get_arg("foldseek-afsp")
fs_csv <- get_arg("foldseek-csv")
foldseek <- NULL
if (!is.null(fs_csv)) {
  foldseek <- fread(fs_csv, integer64 = "double") %>%
    transmute(structure_id = clean_id(idx), bits.pdb = as.numeric(bits.pdb), bits.sp = as.numeric(bits.sp))
} else if (!is.null(fs_pdb) || !is.null(fs_afsp)) {
  empty <- tibble(structure_id = character(), bits = numeric())
  pdb <- if (!is.null(fs_pdb)) top_hit_bits(fs_pdb) else empty
  afsp <- if (!is.null(fs_afsp)) top_hit_bits(fs_afsp, exclude_taxon) else empty
  foldseek <- full_join(rename(pdb, bits.pdb = bits), rename(afsp, bits.sp = bits), by = "structure_id")
}
if (!is.null(foldseek)) {
  foldseek <- foldseek %>%
    # 64-bit integer columns can surface as the 32-bit NA sentinel; treat negative bitscores as missing
    mutate(across(c(bits.pdb, bits.sp), ~ ifelse(.x < 0, NA, .x))) %>%
    mutate(bits.max = pmax(bits.pdb, bits.sp, na.rm = TRUE)) %>%
    distinct(structure_id, .keep_all = TRUE)
  metrics <- metrics %>% left_join(foldseek, by = "structure_id")
}

# Sequence-based results can be keyed by transcript ID or, for deduplicated FASTA files, by structure ID
join_by_either <- function(metrics, tbl, value_col) {
  metrics %>%
    left_join(rename(tbl, v1 = all_of(value_col)), by = c("transcript_id" = "id")) %>%
    left_join(rename(tbl, v2 = all_of(value_col)), by = c("structure_id" = "id")) %>%
    mutate(!!value_col := coalesce(v1, v2)) %>%
    select(-v1, -v2)
}

# ---- InterProScan: total domain length (larger of InterPro and Pfam-only coverage, overlaps merged) ----
interpro_file <- get_arg("interpro")
if (!is.null(interpro_file)) {
  Interpro_results <- fread(interpro_file, header = FALSE, sep = "\t", quote = "", fill = TRUE)
  source(file.path(script_dir, "process_interproscan_results.r"))
  domains <- bind_rows(total_domain_length_df, total_domain_length_df_Pfam) %>%
    group_by(id = transcript_id) %>%
    summarise(total_domain_length = max(total_covered_length), .groups = "drop")
  metrics <- join_by_either(metrics, domains, "total_domain_length")
}

# ---- PSAURON in-frame score ----
psauron_file <- get_arg("psauron")
if (!is.null(psauron_file)) {
  # PSAURON writes the command and an overall score above the header line
  psauron <- fread(psauron_file, skip = "description") %>%
    group_by(id = description) %>%
    summarise(in_frame_score = max(`in-frame_score`, na.rm = TRUE), .groups = "drop")
  metrics <- join_by_either(metrics, psauron, "in_frame_score")
}

fwrite(metrics %>% arrange(Gene, transcript_id), out_file, sep = "\t", na = "NA")
cat(sprintf("Wrote metrics for %d transcripts in %d genes to %s\n",
            nrow(metrics), n_distinct(metrics$Gene), out_file))
for (col in intersect(c("plddt_80", "plddt_70", "bits.max", "total_domain_length", "in_frame_score"), names(metrics))) {
  cat(sprintf("  %-20s missing for %d transcripts\n", col, sum(is.na(metrics[[col]]))))
}

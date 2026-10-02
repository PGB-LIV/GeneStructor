# Final ranking of isoforms by GeneStructor based on ensemble LR model

```r

## Load libraries
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(cowplot)

options(width = 150)
```

## 1. Create a simplified version of merged GFF with only ids as described in [Simplified_GFF.Rmd](Simplified_GFF.Rmd)
```r
#Set working directory
data_dir <- "./data"

# Read the simplified version of Merged GFF
new_gff_data_subset_ids <- fread(file.path(data_dir, "new_gff_data_subset_ids.tab"))

```

## 2. Protein structure results (e.g., Alphafold2 or Protenix.)
Protein structure results were processed using process_cif_for_pLDDT.R. Read in the processed results in long format
In our analysis, Identical transcripts from different databases for same merged gene locus were grouped under same sequence group.
Select relevant columns from the protenix results such as sequence groups, identifiers and number of residues with pLDDT above 70 or 80 (pLDDT70 for Protenix and pLDDT80 for Alphafold).
```r
# Need plddt_70 / plddt_80 values for the model
protenix_res_selected <- fread(file.path(data_dir, "protenix/source_df_long_protenix_cif_stats.csv")) %>%
  select(Source_group, transcript_id, total_residues,plddt_70, plddt_80 , mean_plddt_cif)

```

## 3. Foldseek results
Foldseek results were processed using [process_fs_protenix_all_nip.py](process_fs_protenix_all_nip.py). Read in the foldseek results and process them to select relevant columns and compute the maximum bit score.

```r
foldseek_res <- fread(file.path(data_dir, "foldseek_results_processed.csv"))

# remove .cif from idx column
foldseek_res <- foldseek_res %>%
  mutate(idx = gsub(".cif", "", idx))

# Need maximum bit score for the model

foldseek_res_selected <- foldseek_res %>% 
  select(idx,bits.pdb, bits.sp) %>% 
  mutate(
    bits.pdb = as.numeric(bits.pdb),
    bits.sp = as.numeric(bits.sp)
  ) %>%
  mutate(bits.max = pmax(bits.pdb, bits.sp, na.rm = TRUE))
```

## 4. Combine Protenix and Foldseek results with the simplified GFF data
```r
new_gff_protenix_foldseek_df <- new_gff_data_subset_ids %>% 
    left_join(protenix_res_selected, by = c("ID" = "transcript_id")) %>%
    left_join(foldseek_res_selected, by = c("Source_group" = "idx")) %>%
    select(Parent, ID, Source_group, total_residues, plddt_70, plddt_80, mean_plddt_cif, bits.max ) 
```

## 5. Add interpro total domain length and psauron scores
Total interpro domain length was calculated using [process_interproscan_results.R](process_interproscan_results.r) from Interproscan results.
```r
### Add interpro and Pfam total domain length and keep the maximum value
total_covered_length_df <- fread(file.path(data_dir, "total_domain_length_df.tsv"))

total_covered_length_df_Pfam <- fread(file.path(data_dir, "total_domain_length_df_Pfam.tsv"))


total_covered_length_df_combined <- total_covered_length_df %>%
  bind_rows(total_covered_length_df_Pfam) %>%
  group_by(transcript_id) %>%
  summarise(total_domain_length = max(total_covered_length)) %>%
  ungroup()


### Add psauron scores
psauron_scores <- fread(file.path(data_dir, "protenix/psauron_scores_all3_nip_pep.tsv")) %>% distinct() %>% 
  dplyr::rename(in_frame_score = `in-frame_score` , transcript_id = description)


# combined psauron scores and keep max in_frame_score for each transcript_id
psauron_scores <- psauron_scores %>%
  group_by(transcript_id) %>%
  summarise(in_frame_score = max(in_frame_score, na.rm = TRUE)) %>%
  ungroup()


# Optional but useful: protein length
# protein length

source_groups_protein_length <- fread(file.path(data_dir, "protenix/protein_length_df_groups.csv")) %>% 
select(Source_group, Width) %>% distinct() %>%
rename(protein_length = Width)


# add interpro total domain length and psauron scores to new_gff_protenix_foldseek_df
new_gff_protenix_foldseek_df <- new_gff_protenix_foldseek_df %>%
rename(transcript_id = ID , Gene = Parent) %>%
  left_join(total_covered_length_df_combined, by = c("transcript_id" = "transcript_id")) %>%
  left_join(psauron_scores, by = c("transcript_id" = "transcript_id")) %>%
  left_join(source_groups_protein_length, by = c("Source_group" = "Source_group"))

head(new_gff_protenix_foldseek_df)
```

## 6. Get counts of transcripts per gene
```r
# Get counts of transcripts per gene
transcript_counts <- new_gff_protenix_foldseek_df %>%
  group_by(Gene) %>%
  summarise(transcript_count = n_distinct(transcript_id)) %>%
  ungroup()

transcript_counts %>% group_by(transcript_count) %>% summarise(n_distinct(Gene))

# Separate out those with 1 transcript in another dataframe for comparison later
single_transcript_genes <- new_gff_protenix_foldseek_df %>%
  left_join(transcript_counts, by = "Gene") %>%
  filter(transcript_count == 1) %>%
  select(-transcript_count)

# Filtered for genes with more than 1 transcript 
multiple_transcripts_df <- new_gff_protenix_foldseek_df %>%
  left_join(transcript_counts, by = "Gene") %>%
  filter(transcript_count > 1) %>%
  select(-transcript_count)

# Add transcript_counts to original dataframe
new_gff_protenix_foldseek_df <- new_gff_protenix_foldseek_df %>%
left_join(transcript_counts , by = "Gene")

```

## 7. Reduce redundant transcripts within the same gene and replace NA values with 0
```r
# replace NA with 0 for plddt_80, bits.max, total_domain_length, in_frame_score
new_gff_protenix_foldseek_df <- new_gff_protenix_foldseek_df %>%
  mutate(
    plddt_70 = ifelse(is.na(plddt_70), 0, plddt_70),
    plddt_80 = ifelse(is.na(plddt_80), 0, plddt_80),
    bits.max = ifelse(is.na(bits.max), 0, bits.max),
    total_domain_length = ifelse(is.na(total_domain_length), 0, total_domain_length),
    in_frame_score = ifelse(is.na(in_frame_score), 0, in_frame_score)
  )

# collapse those transcripts within same gene that have same plddt_70, plddt_80, bits.max, total_domain_length, in_frame_score
new_gff_protenix_foldseek_df <- new_gff_protenix_foldseek_df %>%
  group_by(Gene, plddt_70, plddt_80, bits.max, total_domain_length, in_frame_score) %>%
  mutate(transcript_id_collapsed = paste(transcript_id, collapse = ",")) %>%
  slice(1) %>% # keep only the first row of each group
  ungroup()
```

## 8. Apply the ensemble prediction models to the new data and rank isoforms
```r

models_dir <- file.path(data_dir, "models") # Folder with all the saved LR models

source(file.path(data_dir, "predict_with_glm_ensemble.r"))
```

## 9. Apply ensemble prediction and rank isoforms
```r
#options(digits = 2)
# Specify the source of the new data as "protenix" ; The path to the directory with the saved models is given by `models_dir`
all_predictions <- predict_with_glm_ensemble(new_data = new_gff_protenix_foldseek_df, 
                                               models_dir = models_dir, 
                                               n_models = 25 ,
                                               source = "protenix")

# round to 2 decimal places for a sensible ranking
all_predictions <- all_predictions %>%
    group_by(Gene) %>%
    mutate(ML_rank = dense_rank(desc(round(median_pred_CT , 2))))  %>%
    ungroup() %>% 
    arrange(Gene, ML_rank)

# Calculate the difference in predicted probability to the top-ranked isoform within each gene
all_predictions <- all_predictions %>%
    group_by(Gene) %>%
    arrange(ML_rank) %>%
    mutate(prob_diff_to_rank_1 = round(median_pred_CT[ML_rank == 1][1], 2) - round(median_pred_CT, 2)) %>%
    ungroup()

all_predictions <- all_predictions %>% 
left_join(new_gff_protenix_foldseek_df %>% select(Gene, transcript_id, protein_length, plddt_70, plddt_80, bits.max, total_domain_length, in_frame_score , transcript_count, transcript_id_collapsed), 
          by = c("Gene", "transcript_id"))

# Save the predictions to a CSV file
fwrite(all_predictions, file.path(data_dir, "all_predictions.csv"))
         
```
## 10. GeneStructor results for a Sample Gene 

![Gene models for Os06g0113150](https://github.com/user-attachments/assets/df33a3e2-ed68-4c8f-a7d1-14d6aecf3378)

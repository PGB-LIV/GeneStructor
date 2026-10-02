# function_calc_total_IPR_covered_length.R
# Function to merge intervals and calculate domain length from Interproscan results
calc_total_domain_length <- function(df) {
  if(nrow(df) == 0) return(0L)
  merged_intervals <- list()
  current_start <- df$start[1]
  current_end <- df$end[1]
  if(nrow(df) > 1){
    for(i in 2:nrow(df)){
      s <- df$start[i]
      e <- df$end[i]
      if(s <= current_end){
        current_end <- max(current_end, e)
      } else {
        merged_intervals <- append(merged_intervals, list(c(current_start, current_end)))
        current_start <- s
        current_end <- e
      }
    }
  }
  merged_intervals <- append(merged_intervals, list(c(current_start, current_end)))
  sum(sapply(merged_intervals, function(x) x[2] - x[1] + 1L))
}

# For each transcript_id in Interpro results filter for IPR domains, rename start/end, and calculate total_domain_length
Interpro_results_onlyIPR <- Interpro_results %>%
  filter(grepl("^IPR", V12)) %>%
  rename(start = V7, end = V8 , transcript_id = V1) %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end)
  ) %>%
  filter(!is.na(start) & !is.na(end)) %>%
  arrange(transcript_id, start)

Interpro_results_onlyPfam <- Interpro_results %>%
  filter(grepl("Pfam", V4)) %>%
  rename(start = V7, end = V8 , transcript_id = V1) %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end)
  ) %>%
  filter(!is.na(start) & !is.na(end)) %>%
  arrange(transcript_id, start) %>%
    filter(!(grepl("^IPR", V12 )))


# Apply per transcript_id
total_domain_length_df <- Interpro_results_onlyIPR %>%
  group_by(transcript_id) %>%
  arrange(start) %>%
  group_modify(~ tibble(total_covered_length = calc_total_domain_length(.x))) %>%
  ungroup()

total_domain_length_df_Pfam <- Interpro_results_onlyPfam %>%
  group_by(transcript_id) %>%
  arrange(start) %>%
  group_modify(~ tibble(total_covered_length = calc_total_domain_length(.x))) %>%
  ungroup()

  # These data frames contain the total domain length for interpro (IPR) and Pfam domains per protein (transcript_id)
predict_with_glm_ensemble <- function(new_data, models_dir, n_models = n_splits , source) {
  
  library(tidymodels)
  library(dplyr)
  
  cat("=== Ensemble Prediction ===\n")
  
  all_predictions <- list()
  
  for (i in 1:n_models) {
    model_path <- file.path(models_dir, paste0("model_", i, ".rds"))
    
    if (!file.exists(model_path)) {
      warning(sprintf("Model %d not found at %s", i, model_path))
      next
    }
    
    # Load model
    log_fit <- readRDS(model_path)
    
    # if source = "protenix", then consider plddt_70 same as plddt_80 whcih is the default in ensembl models.
    # if source = "alphafold", the make no changes to plddt_80 and use the plddt_80 values as is.
    if (source == "protenix") {
      new_data <- new_data %>%
        mutate(plddt_80 = plddt_70)
    }
    
    # Predict on new data
    preds <- predict(log_fit, new_data = new_data, type = "prob") %>%
      bind_cols(predict(log_fit, new_data)) %>%
      bind_cols(new_data %>% select(transcript_id))
    
    all_predictions[[i]] <- preds
  }
  
  cat(sprintf("Successfully loaded %d models\n", length(all_predictions)))
  
  # Calculate median across all models
  final_scores <- bind_rows(all_predictions, .id = "model_id") %>%
    group_by(transcript_id) %>%
    summarise(
      median_pred_CT = median(.pred_CT),
      median_pred_NCT = median(.pred_NCT),
      mean_pred_CT = mean(.pred_CT),
      sd_pred_CT = sd(.pred_CT),
      min_pred_CT = min(.pred_CT),
      max_pred_CT = max(.pred_CT),
      n_models = n(),
      .groups = "drop"
    )
  
  final_scores <- final_scores %>% 
    left_join(new_data %>% select(Gene, transcript_id), by = "transcript_id")
  cat(sprintf("Generated ensemble predictions for %d transcripts\n", nrow(final_scores)))
  return(final_scores)
}
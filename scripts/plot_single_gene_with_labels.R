# Plot single gene with transcript labels
plot_single_gene_with_labels <- function(gff,
                                         gene_transcripts,
                                         output_file = NULL,
                                         width = 10,
                                         height = NULL) {
  transcript_ids <- unique(gene_transcripts$transcript_id)
  parent_gene <- unique(gene_transcripts$parent_gene_id)[1]
  

  track_spacing <- 0.6
  rect_half_height <- 0.15

  all_features <- gff %>%
    filter(type %in% c("CDS", "five_prime_UTR", "three_prime_UTR")) %>%
    mutate(
      feature = type,
      Parent = str_extract(attributes, "(?<=Parent=)[^;]+")
    ) %>%
    filter(Parent %in% transcript_ids) %>%
    inner_join(gene_transcripts, by = c("Parent" = "transcript_id")) %>%
    arrange(Parent, start)

  db_col <- if ("database" %in% names(all_features)) {
    all_features$database
  } else if ("source" %in% names(all_features)) {
    all_features$source
  } else {
    rep(NA_character_, nrow(all_features))
  }
  
  # Calculate transcript lengths and order by length (descending - longest at top)
  transcript_order <- all_features %>%
    group_by(Parent) %>%
    summarise(
      transcript_start = min(start),
      transcript_end = max(end),
      .groups = "drop"
    ) %>%
    mutate(transcript_length = transcript_end - transcript_start) %>%
    arrange(desc(transcript_length)) %>%
    pull(Parent)
  
  # Calculate dynamic arrow parameters based on gene range
  gene_start <- min(all_features$start)
  gene_end <- max(all_features$end)
  gene_range <- gene_end - gene_start
  arrow_gap <- max(100, gene_range * 0.01)  # 1% of gene range, minimum 100 bp
  arrow_len <- max(20, gene_range * 0.005)  # 0.5% of gene range, minimum 20 bp
  
  # Assign tracks based on length order and determine source database
  db_col <- if ("database" %in% names(all_features)) {
    all_features$database
  } else if ("source" %in% names(all_features)) {
    all_features$source
  } else {
    rep(NA_character_, nrow(all_features))
  }

  all_features <- all_features %>%
    mutate(
      track = as.numeric(factor(Parent, levels = transcript_order)) * track_spacing,
      source_db = case_when(
        !is.na(db_col) & str_detect(as.character(db_col), "(?i)rap") ~ "RAPDB",
        !is.na(db_col) & str_detect(as.character(db_col), "(?i)rgap") ~ "RGAP",
        !is.na(db_col) & str_detect(as.character(db_col), "(?i)gramene") ~ "Gramene",
        str_detect(Parent, "^LOC") ~ "RGAP",
        str_detect(Parent, "^OsNip") ~ "Gramene",
        str_detect(Parent, "^Os(?!Nip)") ~ "RAPDB",
        TRUE ~ ifelse(!is.na(db_col), as.character(db_col), "Other")
      ),
      rect_h = ifelse(feature == "CDS", rect_half_height, rect_half_height * 0.5)
    ) %>%
    group_by(Parent) %>%
    mutate(transcript_label = dplyr::first(Parent)) %>%
    ungroup()

  gene_ranges <- all_features %>%
    group_by(Parent) %>%
    summarise(
      start = min(start),
      end = max(end),
      strand = dplyr::first(strand),
      track = dplyr::first(track),
      transcript_label = dplyr::first(transcript_label),
      source_db = dplyr::first(source_db),
      .groups = "drop"
    )

  if (is.null(height)) {
    height <- max(6, 2 + 0.35 * (nrow(gene_ranges) + 1))
  }

  # Database colors
  db_colors <- c(
    "RAPDB"   = "#E69F00",  # Orange
    "RGAP"     = "#56B4E9",  # Sky Blue
    "Gramene" = "#009E73",  # Green
    "Other"   = "#999999"   # Gray
  )

  # Calculate x-axis limits with padding for labels and arrows
  x_min <- min(gene_ranges$start)
  x_max <- max(gene_ranges$end)
  x_range <- x_max - x_min
  x_padding_left <- x_range * 0.25  # 25% padding on the left for labels
  x_padding_right <- arrow_gap + arrow_len + 500  # Extra space for arrows on the right
  x_padding_left <- max(x_padding_left, arrow_gap + arrow_len + 500)  # Ensure enough space for - strand arrows

  gene_plot <- ggplot() +
    geom_segment(
      data = gene_ranges,
      aes(x = start, xend = end, y = track, yend = track),
      color = "black",
      linewidth = 0.8
    ) +
    geom_rect(
      data = all_features,
      aes(xmin = start, xmax = end, ymin = track - rect_h, ymax = track + rect_h, fill = source_db)
    ) +
    geom_segment(
      data = gene_ranges,
      aes(
        x = ifelse(strand == "+", end + arrow_gap, start - arrow_gap),
        xend = ifelse(strand == "+", end + arrow_gap + arrow_len, start - arrow_gap - arrow_len),
        y = track,
        yend = track
      ),
      arrow = grid::arrow(length = grid::unit(0.1, "inches"), type = "closed"),
      color = "black",
      linewidth = 0.5
    ) +
    geom_text(
      data = gene_ranges,
      aes(x = start, y = track, label = transcript_label, color = source_db),
      hjust = 1.1,
      vjust = 0.5,
      size = 2.5,
      fontface = "bold"
    ) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 5),
      labels = scales::label_number(scale = 1e-6, suffix = " Mb", accuracy = 0.01),
      limits = c(x_min - x_padding_left, x_max + x_padding_right),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(
      values = db_colors,
      name = "Source"
    ) +
    scale_color_manual(
      values = db_colors,
      name = "Source"
    ) +
    xlab("Genomic Coordinates") +
    ylab(NULL) +
    ggtitle(paste("Gene Models for", parent_gene)) +
    theme_cowplot() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.title = element_text(size = 11, hjust = 0.5, margin = margin(b = 4))
    )
  
  if (!is.null(output_file)) {
    ggsave(output_file, plot = gene_plot, width = width, height = height, units = "in", bg = "white")
  }
  
  return(gene_plot)
}


# Enhanced function to plot gene models with ML rank information
plot_gene_with_ml_rank <- function(gff,
                                   gene_transcripts,
                                   output_file = NULL,
                                   width = 10,
                                   height = NULL,
                                   merged_gene_ids = character(),
                                   false_positive_transcript_ids = character()) {
  transcript_ids <- unique(gene_transcripts$transcript_id)
  parent_gene <- unique(gene_transcripts$parent_gene_id)[1]
  merged_gene_ids <- unique(na.omit(str_trim(as.character(merged_gene_ids))))
  false_positive_transcript_ids <- unique(na.omit(str_trim(as.character(false_positive_transcript_ids))))
  is_merged_gene <- parent_gene %in% merged_gene_ids
  
  track_spacing <- 0.6
  
  all_features <- gff %>%
    filter(type %in% c("CDS", "five_prime_UTR", "three_prime_UTR")) %>%
    mutate(
      feature = type,
      Parent = str_extract(attributes, "(?<=Parent=)[^;]+")
    ) %>%
    filter(Parent %in% transcript_ids) %>%
    inner_join(gene_transcripts, by = c("Parent" = "transcript_id")) %>%
    arrange(Parent, start)

  db_col <- if ("database" %in% names(all_features)) {
    all_features$database
  } else if ("source" %in% names(all_features)) {
    all_features$source
  } else {
    rep(NA_character_, nrow(all_features))
  }
  
  # Calculate transcript lengths and order by ML_rank (if available), then by length
  transcript_order <- all_features %>%
    group_by(Parent) %>%
    summarise(
      transcript_start = min(start),
      transcript_end = max(end),
      ML_rank = dplyr::first(ML_rank),
      .pred_CT = dplyr::first(.pred_CT),
      database = dplyr::first(database),
      .groups = "drop"
    ) %>%
    mutate(transcript_length = transcript_end - transcript_start) %>%
    arrange(ML_rank, desc(transcript_length)) %>%
    pull(Parent)

  rect_half_height <- min(0.06, 0.02 * length(transcript_order))
  
  # Calculate dynamic arrow parameters based on gene range
  gene_start <- min(all_features$start)
  gene_end <- max(all_features$end)
  gene_range <- gene_end - gene_start
  arrow_gap <- max(100, gene_range * 0.001)
  arrow_len <- max(20, gene_range * 0.005)
  
  # Assign tracks based on order and create enhanced labels with ML rank
  all_features <- all_features %>%
    mutate(
      track = as.numeric(factor(Parent, levels = transcript_order)) * track_spacing,
      source_db = case_when(
        !is.na(db_col) & str_detect(as.character(db_col), "(?i)rap") ~ "RAPDB",
        !is.na(db_col) & str_detect(as.character(db_col), "(?i)rgap") ~ "RGAP",
        !is.na(db_col) & str_detect(as.character(db_col), "(?i)gramene") ~ "Gramene",
        str_detect(Parent, "^LOC") ~ "RGAP",
        str_detect(Parent, "^OsNip") ~ "Gramene",
        str_detect(Parent, "^Os(?!Nip)") ~ "RAPDB",
        TRUE ~ ifelse(!is.na(db_col), as.character(db_col), "Other")
      )
    ) %>%
    group_by(Parent) %>%
    mutate(
      transcript_label = paste0(
        dplyr::first(Parent),
        "\n(Rank ", dplyr::first(ML_rank),
        ", p=", round(dplyr::first(.pred_CT), 2), ")"
      ),
      false_positive_label = ifelse(
        dplyr::first(Parent) %in% false_positive_transcript_ids,
        "FALSE POSITIVE",
        NA_character_
      )
    ) %>%
    ungroup()
  
  gene_ranges <- all_features %>%
    group_by(Parent) %>%
    summarise(
      start = min(start),
      end = max(end),
      strand = dplyr::first(strand),
      track = dplyr::first(track),
      transcript_label = dplyr::first(transcript_label),
      false_positive_label = dplyr::first(false_positive_label),
      source_db = dplyr::first(source_db),
      .groups = "drop"
    )

  if (is.null(height)) {
    height <- max(6, 2 + 0.35 * (nrow(gene_ranges) + 1))
  }
  
  # Calculate x-axis limits with padding for labels and arrows
  x_min <- min(gene_ranges$start)
  x_max <- max(gene_ranges$end)
  x_range <- x_max - x_min
  label_padding <- x_range * 0.20
  arrow_padding <- arrow_gap + arrow_len + 500
  x_padding_left <- max(
    if (any(gene_ranges$strand == "+", na.rm = TRUE)) label_padding else 0,
    arrow_padding
  )
  x_padding_right <- max(
    if (any(gene_ranges$strand != "+", na.rm = TRUE)) label_padding else 0,
    arrow_padding
  )

  db_colors <- c(
    "RAPDB" = "#E69F00",
    "RGAP" = "#56B4E9",
    "Gramene" = "#009E73",
    "Other" = "#999999"
  )
  
  gene_plot <- ggplot() +
    geom_segment(
      data = gene_ranges,
      aes(x = start, xend = end, y = track, yend = track),
      color = "black",
      linewidth = 1.2
    ) +
    geom_rect(
      data = all_features,
      aes(xmin = start, xmax = end, ymin = track - rect_half_height, ymax = track + rect_half_height, fill = source_db),
      color = "white",
      linewidth = 0.3
    ) +
    geom_segment(
      data = gene_ranges,
      aes(
        x = ifelse(strand == "+", end + arrow_gap, start - arrow_gap),
        xend = ifelse(strand == "+", end + arrow_gap + arrow_len, start - arrow_gap - arrow_len),
        y = track,
        yend = track
      ),
      arrow = grid::arrow(length = grid::unit(0.1, "inches"), type = "closed"),
      color = "black",
      linewidth = 0.8
    ) +
    geom_text(
      data = gene_ranges,
      aes(
        x = ifelse(strand == "+", start, end),
        y = track,
        label = transcript_label,
        hjust = ifelse(strand == "+", 1.05, -0.05)
      ),
      color = "black",
      vjust = 0.5,
      size = 2.8,
      fontface = "bold"
    ) +
    geom_text(
      data = gene_ranges %>% filter(!is.na(false_positive_label)),
      aes(
        x = ifelse(strand == "+", start, end),
        y = track - track_spacing / 2,
        label = false_positive_label,
        hjust = ifelse(strand == "+", 1.05, -0.05)
      ),
      color = "#e90f28",
      vjust = 0.5,
      size = 2.8,
      fontface = "bold"
    ) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 5),
      labels = scales::label_number(scale = 1e-6, suffix = " Mb", accuracy = 0.01),
      limits = c(x_min - x_padding_left, x_max + x_padding_right),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(
      values = db_colors,
      name = "Source"
    ) +
    scale_y_continuous(
  limits = c(
    min(gene_ranges$track) - track_spacing,
    max(gene_ranges$track) + track_spacing
  ),
  expand = expansion(mult = c(0, 0))
) +
    xlab("Genomic Coordinates") +
    ylab(NULL) +
    ggtitle(
      paste("Gene Models for", parent_gene, "- Colored by Database"),
      subtitle = if (is_merged_gene) "MERGED GENE DETECTED" else NULL
    ) +
    theme_cowplot(font_size = 11) +
    theme(
      legend.position = "bottom",
      plot.subtitle = element_text(size = 12, hjust = 0.5, face = "bold", color = "#e90f28"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA)
    )
  
  if (!is.null(output_file)) {
    ggsave(output_file, plot = gene_plot, width = width, height = height, units = "in", bg = "white")
  }
  
  return(gene_plot)
}

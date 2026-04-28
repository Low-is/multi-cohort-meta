
  prepare_pathway_df <- function(up_list, down_list, n_show = 5) {
    
    # Function to process each pathway enrichment dataframe
    process_pathway <- function(df, direction, n_show) {
      if (is.null(df)) return(NULL)  # skip if NULL
      df %>%
        arrange(p.adjust) %>%
        slice_head(n = n_show) %>%
        mutate(direction = direction,
               minus_log10_p = -log10(p.adjust),
               term = Description)
    }
    
    # Process up-regulated pathways
    up_df <- imap(up_list, ~ process_pathway(.x, direction = 1, n_show = n_show)) %>%
      compact() %>%  # remove NULLs
      bind_rows(.id = "study")
    
    # Process down-regulated pathways
    down_df <- imap(down_list, ~ process_pathway(.x, direction = -1, n_show = n_show)) %>%
      compact() %>%  # remove NULLs
      bind_rows(.id = "study")
    
    # Combine up and down
    combined_df <- bind_rows(up_df, down_df) %>%
      mutate(count_dir = Count * direction,
             signed_logp = minus_log10_p * direction)
    
    return(combined_df)
  }
  
  
  
  plot_pathways <- function(pathway_df, combine = TRUE, title = NULL, save_dir = NULL) {
    
    # Skip if pathway_df is empty
    if (nrow(pathway_df) == 0) {
      message("No pathways to plot. Exiting.")
      return(NULL)
    }
    
    # Determine save directory
    if (isTRUE(save_dir)) save_dir <- getwd()
    
    # Add signed logp and count_dir
    pathway_df <- pathway_df %>%
      dplyr::mutate(
        signed_logp = minus_log10_p * direction,
        count_dir = Count * direction
      ) %>%
      # Create term_ordered for plotting
      dplyr::arrange(direction, desc(abs(count_dir))) %>%
      dplyr::mutate(term_ordered = factor(term, levels = unique(term)))
    
    if (combine) {
      p <- ggplot(pathway_df, aes(x = term_ordered,
                                  y = count_dir,
                                  fill = signed_logp)) +
        geom_col(width = 0.6, colour = "black") +
        scale_fill_gradient2(low = "#7CAE00", mid = "white", high = "#C77CFF",
                             midpoint = 0, name = "-logP") +
        geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
        coord_flip() +
        labs(x = NULL, y = "Gene counts", title = title) +
        theme_classic(base_size = 12) +
        theme(axis.text.y = element_text(hjust = 1),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              legend.title = element_text(size = 9),
              legend.position = "right",
              plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
      
      # Save if directory is specified
      if (!is.null(save_dir)) {
        ggsave(filename = file.path(save_dir, "combined_pathways.png"),
               plot = p,
               width = 10,
               height = 6,
               dpi = 300,
               units = "in")
      }
      
      return(p)
      
    } else {
      # Split per study, skip studies with 0 rows
      plots <- pathway_df %>%
        split(.$study) %>%
        keep(~ nrow(.x) > 0) %>%  # skip empty studies
        map(~ ggplot(.x, aes(x = term_ordered,
                             y = count_dir,
                             fill = signed_logp)) +
              geom_col(width = 0.6, colour = "black") +
              scale_fill_gradient2(low = "#F68B1F", mid = "white", high = "#2166AC",
                                   midpoint = 0, name = "-logP") +
              geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
              coord_flip() +
              labs(x = NULL, y = "Gene counts", title = unique(.x$study)) +
              theme_classic(base_size = 12) +
              theme(axis.text.y = element_text(hjust = 1),
                    axis.ticks.y = element_blank(),
                    axis.line.y = element_blank(),
                    legend.title = element_text(size = 9),
                    legend.position = "right",
                    plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
        )
      
      # Save each plot if save_dir specified
      if (!is.null(save_dir)) {
        for (study_name in names(plots)) {
          ggsave(filename = file.path(save_dir, paste0(study_name, "_pathways.png")),
                 plot = plots[[study_name]],
                 width = 8,
                 height = 6,
                 dpi = 300,
                 units = "in")
        }
      }
      
      return(plots)
    }
  }

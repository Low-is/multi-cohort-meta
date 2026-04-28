 plot_immune_deconv_list <- function(
    dna_list,
    rna_list,
    dna_pData,
    rna_pData,
    dna_method = "mcp_counter",
    rna_method = "quantiseq",
    cell_types = NULL,
    genes = NULL,                # <<---- NEW ARGUMENT
    save_path = NULL,
    dpi = 600
  ) {
    
    all_plots <- list()
    
    # ---- Function to process a single dataset ----
    process_dataset <- function(expr_matrix, pData, method) {
      
      # ---- OPTIONAL GENE FILTER ----
      if (!is.null(genes)) {
        keep <- intersect(genes, rownames(expr_matrix))
        if (length(keep) == 0) {
          stop("None of the provided genes were found in this expression matrix.")
        }
        expr_matrix <- expr_matrix[keep, , drop = FALSE]
      }
      # ----------------------------------------------
      
      # Deconvolution
      deconv_res <- immunedeconv::deconvolute(expr_matrix, method)
      
      # Remove rows with all NAs
      deconv_res <- deconv_res[rowSums(is.na(deconv_res[,-1])) < ncol(deconv_res[,-1]), ]
      
      # Pivot longer
      long_df <- deconv_res %>%
        tidyr::pivot_longer(
          cols = -cell_type,
          names_to = "gsm",
          values_to = "fraction"
        ) %>%
        dplyr::left_join(pData %>% dplyr::select(gsm, condition), by = "gsm")
      
      
      # ---- Normalize DNA MCP-counter scores to 0-1 ----
      if (method == "mcp_counter") {
        
        long_df <- long_df %>%
          dplyr::filter(!grepl("cytotoxic", cell_type, ignore.case = TRUE))
        
        long_df <- long_df %>%
          dplyr::group_by(cell_type) %>%
          dplyr::mutate(
            fraction = (fraction - min(fraction, na.rm = TRUE)) /
              (max(fraction, na.rm = TRUE) - min(fraction, na.rm = TRUE))
          ) %>%
          dplyr::ungroup()
      }
      
      # Summary
      long_df_summary <- long_df %>%
        dplyr::group_by(cell_type, condition) %>%
        dplyr::summarize(avg_fraction = mean(fraction), .groups = "drop")
      
      # P-values
      cell_pvals <- long_df %>%
        dplyr::group_by(cell_type) %>%
        dplyr::summarize(p_value = wilcox.test(fraction ~ condition)$p.value,
                         .groups = "drop") %>%
        dplyr::mutate(
          signif = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01  ~ "**",
            p_value < 0.05  ~ "*",
            TRUE            ~ "ns"
          ),
          label = paste0("p = ", signif(p_value, 3), " ", signif)
        )
      
      long_df_summary <- long_df_summary %>% left_join(cell_pvals, by = "cell_type")
      
      # Filter cell types
      if (!is.null(cell_types)) {
        filtered_cell_types <- intersect(cell_types, unique(long_df_summary$cell_type))
      } else {
        filtered_cell_types <- long_df_summary %>%
          dplyr::filter(avg_fraction > 0 & p_value <= 0.05) %>%
          dplyr::pull(cell_type) %>% unique()
      }
      
      long_df_summary <- long_df_summary %>% filter(cell_type %in% filtered_cell_types)
      
      pval_labels <- long_df_summary %>%
        dplyr::distinct(cell_type, label) %>%
        dplyr::mutate(y_pos = max(long_df_summary$avg_fraction) + 0.1)
      
      # Plot
      p <- ggplot2::ggplot(long_df_summary, aes(x = condition, y = avg_fraction, fill = condition)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), colour = "black", linewidth = 1.2) +
        ggplot2::geom_text(
          data = pval_labels,
          aes(x = 1.5, y = y_pos, label = label),
          inherit.aes = FALSE,
          size = 6
        ) +
        ggplot2::facet_wrap(~cell_type, scales = "free_y") +
        ggplot2::labs(
          x = "", y = "Average Cell Fraction", fill = "Condition"
        ) +
        ggplot2::scale_y_continuous(limits = c(0.0, 1)) +
        ggplot2::theme_classic() +
        ggplot2::scale_fill_manual(values = c("NoNEC" = "#2CA02C",
                                              "NEC"   = "#9467BD")) +
        ggplot2::theme(plot.margin = unit(c(.75,.75,.75,.75), "inches"),
                       strip.text = ggplot2::element_text(size = 14, face = "bold"),
                       axis.text = ggplot2::element_text(size = 16),
                       axis.title = ggplot2::element_text(size = 16),
                       legend.position = "none")
      
      return(list(plot = p, long_df = long_df, summary = long_df_summary, pvals = cell_pvals))
    }
    
    # ---- Updated workflow ----
    all_plots <- list()
    get_plot_name <- function(lst, prefix, i) {
      if (!is.null(names(lst)) && names(lst)[i] != "") names(lst)[i] else paste0(prefix, "_", i)
    }
    
    # DNA datasets
    for (i in seq_along(dna_list)) {
      expr_matrix <- dna_list[[i]]
      pData <- dna_pData[[i]]
      plot_name <- get_plot_name(dna_list, "DNA", i)
      all_plots[[plot_name]] <- process_dataset(expr_matrix, pData, dna_method)
    }
    
    # RNA datasets
    for (i in seq_along(rna_list)) {
      expr_matrix <- rna_list[[i]]
      pData <- rna_pData[[i]]
      plot_name <- get_plot_name(rna_list, "RNA", i)
      all_plots[[plot_name]] <- process_dataset(expr_matrix, pData, rna_method)
    }
    
    # Save images
    if (!is.null(save_path)) {
      if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
      
      for (name in names(all_plots)) {
        ggsave(
          filename = file.path(save_path, paste0(name, ".png")),
          plot = all_plots[[name]]$plot,
          dpi = dpi, width = 10, height = 8, units = "in"
        )
      }
    }
    
    return(all_plots)
  }

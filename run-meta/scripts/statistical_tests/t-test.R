# student's t-test and plotting 
  
  plot_ttest_genes <- function(list_of_mtx,
                               list_of_pData,
                               genes_of_interest = NULL,
                               save_plots = FALSE,
                               prefix = "dataset_boxplot",
                               width = 8,
                               height = 6,
                               dpi = 600) {
    
    results_list <- list()
    
    # Loop through each dataset separately
    for (i in seq_along(list_of_mtx)) {
      
      mtx <- list_of_mtx[[i]]
      pdata <- list_of_pData[[i]]
      dataset_name <- names(list_of_mtx)[i]
      if (is.null(dataset_name)) dataset_name <- paste0("Dataset_", i)
      
      # Ensure condition factor
      condition <- factor(pdata$condition)
      if (nlevels(condition) != 2)
        stop(paste0("Dataset ", dataset_name, " does not have exactly 2 condition groups."))
      
      # ------------------------------
      # 1. Z-score normalization
      # ------------------------------
      # mtx_z <- t(scale(t(mtx)))  # z-score per gene
      
      # ------------------------------
      # 2. Select genes
      # ------------------------------
      if (!is.null(genes_of_interest)) {
        genes_keep <- intersect(genes_of_interest, rownames(mtx))
        if (length(genes_keep) == 0)
          stop(paste0("No requested genes found in ", dataset_name))
      } else {
        genes_keep <- rownames(mtx)  # use all genes
      }
      
      mtx_sub <- mtx[genes_keep, , drop = FALSE]
      
      # ------------------------------
      # 3. Build long dataframe
      # ------------------------------
      df <- data.frame(
        gene       = rep(rownames(mtx_sub), each = ncol(mtx_sub)),
        expression = as.vector(as.matrix(mtx_sub)),
        condition  = rep(condition, times = nrow(mtx_sub))
      )
      
      # ------------------------------
      # 4. Compute t-test p-values per gene
      # ------------------------------
      pvals <- df %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(
          p_value = t.test(expression ~ condition, var.equal = FALSE)$p.value, # changed 'var.equal = TRUE' to FALSE
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          sig = dplyr::case_when(
            p_value <= 0.001 ~ "***",
            p_value <= 0.01  ~ "**",
            p_value <= 0.05  ~ "*",
            TRUE             ~ "ns"
          ),
          label = paste0("p=", signif(p_value, 3), " ", sig)
        )
      
      # ------------------------------
      # 5. Annotation positions
      # ------------------------------
      y_max <- df %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(y = max(expression, na.rm = TRUE) * 1.15)
      
      annotation_df <- dplyr::left_join(pvals, y_max, by = "gene") %>%
        dplyr::rename(p = label)
      
      # ------------------------------
      # 6. Create boxplots for this dataset
      # ------------------------------
      p <- ggplot(df, aes(x = condition, y = expression, fill = condition)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
        facet_wrap(~ gene, scales = "free_y") +
        geom_text(
          data = annotation_df,
          aes(x = 1.5, y = y, label = p),
          inherit.aes = FALSE,
          size = 3.5
        ) +
        labs(
          title = paste0("Expression Boxplots: ", dataset_name),
          x = NULL,
          y = "Normalized Expression"  # changed 'Z-scored Expression' to 'Normalized Expression'
        ) +
        scale_fill_manual(values = c("NoNEC" = "#2CA02C", # chaning NoSepsis to NoNEC
                                     "NEC"   = "#9467BD")) + # changing Sepsis to NEC
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank(),
              legend.position = "none",
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"))
      
      # ------------------------------
      # 7. Optionally save plot
      # ------------------------------
      if (save_plots) {
        fname <- paste0(prefix, "_", dataset_name, ".png")
        ggsave(fname, plot = p, dpi = dpi, width = width, height = height, units = "in")
      }
      
      # ------------------------------
      # 8. Store results
      # ------------------------------
      results_list[[dataset_name]] <- list(
        plot = p,
        p_values = annotation_df,
        data_used = df
      )
    }
    
    return(results_list)
  }

# student's t-test and plotting

### Plots wilcox of combined matrices ####
plot_wilcox_genes <- function(combined_mtx_res,
                              genes_of_interest,
                              save_plot = FALSE,
                              filename  = "wilcox_boxplots.png",
                              width     = 8,
                              height    = 6,
                              dpi       = 600) {

    # Extract components
    expr_mat  <- combined_mtx_res$expr_matrix
    condition <- combined_mtx_res$condition

    # Keep only valid genes
    genes_of_interest <- intersect(genes_of_interest, rownames(expr_mat))
    if (length(genes_of_interest) == 0) {
      stop("No valid genes found.")
    }

    # Subset matrix
    sub_expr <- expr_mat[genes_of_interest, , drop = FALSE]

    # Long-format dataframe
    df <- data.frame(
      gene       = rep(rownames(sub_expr), each = ncol(sub_expr)),
      expression = as.vector(as.matrix(sub_expr)),
      condition  = rep(condition, times = nrow(sub_expr))
    )

    # Compute Wilcoxon p-values per gene
    pvals <- df %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(
        p_value = wilcox.test(expression ~ condition)$p.value,
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        sig = dplyr::case_when(
          p_value <= 0.001 ~ "***",
          p_value <= 0.01  ~ "**",
          p_value <= 0.05  ~ "*",
          TRUE             ~ "ns"
        ),
        label = paste0("p = ", signif(p_value, 4), " ", sig)
      )

    # Build annotation dataframe like your example
    y_max <- df %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(y = max(expression, na.rm = TRUE) * 1.15)

    custom_p_L <- dplyr::left_join(pvals, y_max, by = "gene") %>%
      dplyr::rename(p = label)

    # Create plot
    p <- ggplot(df, aes(x = condition, y = expression, fill = condition)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
      facet_wrap(~ gene, scales = "free_y") +
      geom_text(
        data = custom_p_L,
        aes(x = 1.5, y = y, label = p),
        inherit.aes = FALSE,
        size = 3.5
      ) +
      labs(x = NULL, y = "Scaled Expression") +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank()) +
      theme(
        legend.position  = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold")
      )

    # Save if requested
    if (save_plot) {
      ggsave(
        filename = filename,
        plot     = p,
        dpi      = dpi,
        width    = width,
        height   = height,
        units    = "in"
      )
    }

    # Return both plot and p-values table
    return(list(
      plot = p,
      p_values = custom_p_L
    ))
  }
### Plots wilcox of combined matrices ####

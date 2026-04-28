library(shiny)
library(DT)
library(bslib)





# --------------------------
# Shiny Server
# --------------------------
server <- function(input, output, session) {
  
  pdata_reactive <- reactiveVal(list())
  expr_reactive  <- reactiveVal(list())
  
  # --------------------------
  # Load GEO Studies + Expression Matrices
  # --------------------------
  observeEvent(input$load_data, {
    
    geo_vec <- str_trim(unlist(strsplit(input$geo_ids, ",\\s*")))
    
    if(length(geo_vec) == 0){
      showNotification("Please enter at least one GEO accession.", type = "error")
      return(NULL)
    }
    
    # Build study lists properly
    dna_studies <- if("DNA" %in% input$study_type) {
      stats::setNames(as.list(geo_vec), geo_vec)
    } else {
      list()
    }
    
    rna_studies <- if("RNA" %in% input$study_type) {
      stats::setNames(as.list(geo_vec), geo_vec)
    } else {
      list()
    }
    
    if(length(dna_studies) == 0 && length(rna_studies) == 0){
      showNotification("Select at least one study type (DNA or RNA).", type = "error")
      return(NULL)
    }
    
    withProgress(message = "Loading Studies...", value = 0, {
      
      incProgress(0.2, detail = "Downloading pData...")
      
      # Combine studies for pData retrieval
      combined_studies <- c(dna_studies, rna_studies)
      pdata_list <- get_pData(combined_studies)
      pdata_reactive(pdata_list)
      
      updateSelectInput(session, "study_to_inspect", 
                        choices = names(pdata_list))
      updateSelectInput(session, "study_for_column", 
                        choices = names(pdata_list))
      
      incProgress(0.4, detail = "Generating expression matrices...")
      
      expr_list <- generate_exprs_mtx(
        DNA = "DNA" %in% input$study_type,
        RNA = "RNA" %in% input$study_type,
        dna_studies = dna_studies,
        rna_studies = rna_studies
      )
      
      expr_reactive(expr_list)
      
      updateSelectInput(session, "study_to_inspect_expr",
                        choices = names(expr_list))
      
      updateSelectInput(session, "study_to_clean_single",
                        choices = names(expr_list))
      
      incProgress(0.4, detail = "Done.")
    })
    
    showNotification("Studies and expression matrices loaded successfully!")
  })
  
  
  observeEvent(input$load_data, {
    
    req(expr_reactive())
    
    updateCheckboxGroupInput(
      session,
      "studies_to_clean",
      choices = names(expr_reactive()),
      selected = NULL
    )
  })
  
  
  # --------------------------
  # Inspect pData
  # --------------------------
  output$pdata_structure_table <- renderDT({
    req(pdata_reactive(), input$study_to_inspect)
    datatable(
      head(pdata_reactive()[[input$study_to_inspect]], 50),
      options = list(scrollX = TRUE)
    )
  })
  
  
  # --------------------------
  # Inspect Expression Matrix
  # --------------------------
  output$expr_structure_table <- renderDT({
    req(expr_reactive(), input$study_to_inspect_expr)
    
    expr_mat <- expr_reactive()[[input$study_to_inspect_expr]]
    
    # Convert matrix to data.frame for DT
    expr_df <- as.data.frame(expr_mat)
    
    datatable(
      head(expr_df, 50),
      options = list(scrollX = TRUE)
    )
  })
  
  
  # --------------------------
  # Add Condition Column
  # --------------------------
  observeEvent(input$study_for_column, {
    req(pdata_reactive(), input$study_for_column)
    
    study_name <- input$study_for_column
    cols <- colnames(pdata_reactive()[[study_name]])
    
    updateSelectInput(session, "column_name_input", choices = cols)
  })
  
  
  observeEvent(input$convert_column, {
    
    req(pdata_reactive(), input$study_for_column, input$column_name_input)
    
    pdata_list  <- pdata_reactive()
    study_name  <- input$study_for_column
    column_name <- input$column_name_input
    
    if (!(column_name %in% colnames(pdata_list[[study_name]]))) {
      showNotification(
        paste("Column", column_name, "not found."),
        type = "error"
      )
      return(NULL)
    }
    
    # Parse user regex patterns (GENERALIZABLE)
    case_patterns <- str_trim(unlist(strsplit(input$case_patterns_text, ",")))
    control_patterns <- str_trim(unlist(strsplit(input$control_patterns_text, ",")))
    
    # Remove empty strings
    case_patterns <- case_patterns[case_patterns != ""]
    control_patterns <- control_patterns[control_patterns != ""]
    
    if (length(case_patterns) == 0 || length(control_patterns) == 0) {
      showNotification(
        "Please enter at least one pattern for Case and Control.",
        type = "error"
      )
      return(NULL)
    }
    
    # Build regex
    case_regex <- paste(case_patterns, collapse = "|")
    control_regex <- paste(control_patterns, collapse = "|")
    
    column_values <- as.character(pdata_list[[study_name]][[column_name]])
    
    # Assign condition using regex (case insensitive)
    new_condition <- ifelse(
      grepl(case_regex, column_values, ignore.case = TRUE),
      "Case",
      ifelse(
        grepl(control_regex, column_values, ignore.case = TRUE),
        "Control",
        NA
      )
    )
    
    pdata_list[[study_name]]$condition <- factor(
      new_condition,
      levels = c("Control", "Case")
    )
    
    pdata_reactive(pdata_list)
    
    # Summary message
    case_n <- sum(new_condition == "Case", na.rm = TRUE)
    control_n <- sum(new_condition == "Control", na.rm = TRUE)
    na_n <- sum(is.na(new_condition))
    
    showNotification(
      paste0(
        "Condition column added to ", study_name,
        " | Case: ", case_n,
        " | Control: ", control_n,
        " | Unmatched: ", na_n
      ),
      type = "message"
    )
  })
  
  
  # --------------------------
  # Preview Selected Column
  # --------------------------
  output$column_preview_table <- renderDT({
    req(pdata_reactive(), input$study_for_column, input$column_name_input)
    
    pdata_list  <- pdata_reactive()
    study_name  <- input$study_for_column
    column_name <- input$column_name_input
    
    if (column_name %in% colnames(pdata_list[[study_name]])) {
      datatable(
        pdata_list[[study_name]][, column_name, drop = FALSE],
        options = list(
          scrollX = TRUE,
          pageLength = 25,
          lengthMenu = c(25, 50, 100, 500)
        ),
        caption = paste("Preview of column:", column_name)
      )
    }
  })
  
  
  
  # --------------------------
  # Show pre-cleaning dimensions 
  # --------------------------
  output$clean_dims_before <- renderPrint({
    req(expr_reactive(), pdata_reactive(), input$study_to_clean_single)
    
    expr_list  <- expr_reactive()
    pdata_list <- pdata_reactive()
    study      <- input$study_to_clean_single
    
    if (!study %in% names(expr_list) || !study %in% names(pdata_list)) {
      cat("Selected study not found in current data.\n")
      return(invisible(NULL))
    }
    
    expr_mat <- expr_list[[study]]
    pdata    <- pdata_list[[study]]
    
    cat("Before cleaning:\n")
    cat("  Expression matrix dimensions: ",
        paste(dim(expr_mat), collapse = " x "), "\n", sep = "")
    cat("  pData dimensions: ",
        paste(dim(pdata), collapse = " x "), "\n", sep = "")
    if ("condition" %in% colnames(pdata)) {
      cat("  Condition table:\n")
      print(table(pdata$condition, useNA = "ifany"))
    } else {
      cat("  'condition' column not found in pData.\n")
    }
  })
  
  
  
  observeEvent(input$apply_cleaning_single, {
    req(expr_reactive(), pdata_reactive(), input$study_to_clean_single)
    
    expr_list  <- expr_reactive()
    pdata_list <- pdata_reactive()
    study      <- input$study_to_clean_single
    
    if (!study %in% names(expr_list) || !study %in% names(pdata_list)) {
      showNotification("Selected study not found.", type = "error")
      return(NULL)
    }
    
    pdata <- pdata_list[[study]]
    expr_mat <- expr_list[[study]]
    
    if (!"condition" %in% colnames(pdata)) {
      showNotification("Condition column not found for this study. Please create it first in 'Add Column' tab.", type = "error")
      return(NULL)
    }
    
    cond <- pdata$condition
    
    # build keep mask
    keep <- rep(TRUE, length(cond))
    
    # keep only selected condition levels
    if (!is.null(input$condition_levels_keep) && length(input$condition_levels_keep) > 0) {
      keep <- keep & cond %in% input$condition_levels_keep
    }
    
    # option to drop NAs in condition
    if (isTRUE(input$drop_na_condition)) {
      keep <- keep & !is.na(cond)
    }
    
    if (!any(keep)) {
      showNotification("No samples left after cleaning with current settings.", type = "error")
      return(NULL)
    }
    
    # subset pData and expression matrix
    pdata_clean <- pdata[keep, , drop = FALSE]
    
    # assume columns of expr_mat correspond to samples in pData (your workflow)
    if (ncol(expr_mat) != nrow(pdata)) {
      showNotification("Mismatch between number of samples in expression matrix and pData before cleaning.", type = "error")
      return(NULL)
    }
    
    expr_clean <- expr_mat[, keep, drop = FALSE]
    
    # store back
    pdata_list[[study]] <- pdata_clean
    expr_list[[study]]  <- expr_clean
    
    pdata_reactive(pdata_list)
    expr_reactive(expr_list)
    
    showNotification(
      paste0(
        "Cleaning complete for ", study,
        ". Remaining samples: ", nrow(pdata_clean),
        " (", paste(table(pdata_clean$condition), collapse = ", "), ")."
      ),
      type = "message"
    )
  })
  
  
  # --------------------------
  # Apply cleaning and show post dimensions
  # --------------------------
  output$clean_dims_after <- renderPrint({
    req(expr_reactive(), pdata_reactive(), input$study_to_clean_single)
    
    expr_list  <- expr_reactive()
    pdata_list <- pdata_reactive()
    study      <- input$study_to_clean_single
    
    if (!study %in% names(expr_list) || !study %in% names(pdata_list)) {
      cat("Selected study not found in current data.\n")
      return(invisible(NULL))
    }
    
    expr_mat <- expr_list[[study]]
    pdata    <- pdata_list[[study]]
    
    cat("After cleaning (current in-memory objects):\n")
    cat("  Expression matrix dimensions: ",
        paste(dim(expr_mat), collapse = " x "), "\n", sep = "")
    cat("  pData dimensions: ",
        paste(dim(pdata), collapse = " x "), "\n", sep = "")
    if ("condition" %in% colnames(pdata)) {
      cat("  Condition table:\n")
      print(table(pdata$condition, useNA = "ifany"))
    } else {
      cat("  'condition' column not found in pData.\n")
    }
    
    if (ncol(expr_mat) != nrow(pdata)) {
      cat("\nWARNING: ncol(expr) != nrow(pData); meta-analysis will fail.\n")
    }
  })
  
  
  # --------------------------
  # Reactive to get the gene list from upload, or NULL if nothing uploaded
  # --------------------------
  selected_genes <- reactive({
    if (is.null(input$gene_file)) return(NULL)  # No file uploaded → auto-detect
    
    # Read uploaded file (support CSV or TXT)
    tryCatch({
      ext <- tools::file_ext(input$gene_file$name)
      
      if (ext == "txt") {
        genes <- scan(input$gene_file$datapath, what = character(), sep = "\n")
      } else if (ext == "csv") {
        genes <- read.csv(input$gene_file$datapath, header = FALSE, stringsAsFactors = FALSE)[,1]
      } else {
        stop("Unsupported file type")
      }
      
      unique(genes)
      
    }, error = function(e) {
      showNotification("Failed to read uploaded gene list. Using auto-detection.", type = "warning")
      return(NULL)
    })
  })
  
  
  # --------------------------
  # Run Meta-Analysis
  # --------------------------
  results <- eventReactive(input$run_analysis, {
    
    req(expr_reactive(), pdata_reactive())
    
    # Grab current reactive data
    expr_list  <- expr_reactive()
    pdata_list <- pdata_reactive()
    
    # --------------------------
    # CLEAN EACH STUDY
    # --------------------------
    cleaned <- clean_study_for_meta(expr_list, pdata_list)
    expr_list  <- cleaned$expr_list
    pdata_list <- cleaned$pdata_list
    
    # --------------------------
    # CHECK THAT EACH STUDY HAS BOTH CASE & CONTROL
    # --------------------------
    for (study in names(expr_list)) {
      g <- pdata_list[[study]]$condition
      if (length(unique(g)) < 2) {
        showNotification(
          paste("Study", study, "does not contain both Case and Control after cleaning."),
          type = "error"
        )
        return(NULL)
      }
    }
    
    genes_to_use <- selected_genes()
    
    # --------------------------
    # RUN META-ANALYSIS
    # --------------------------
    meta_res <- generate_list_for_meta_analysis(
      DNA = TRUE,
      list_of_dna_mtx = expr_list,
      list_of_pData  = pdata_list,
      study = names(expr_list),
      common_genes = genes_to_use
    )
    
    # Extract full meta-analysis table
    meta_res$all_genes <- rownames(meta_res$g)

    
    # Store back cleaned data into reactives
    expr_reactive(expr_list)
    pdata_reactive(pdata_list)
    
    showNotification("Meta-analysis complete!")
    
    # --------------------------
    # RETURN FULL TABLE AND PLOT DATA (no thresholds yet)
    # --------------------------
    list(
      table = meta_res$pooled_estimates,  # full table, thresholds applied later
      plot_data = meta_res,
      all_genes = rownames(meta_res$summary$g)
    )
  })
  
  # --------------------------
  # Reactive: dynamically filter genes by user thresholds
  # --------------------------
  filtered_genes <- reactive({
    req(results())
    
    meta_obj <- results()$plot_data
    
    # Use top-level vectors/matrices from meta_res
    pool      <- meta_obj$pooled$summary
    se.pool   <- meta_obj$pooled$se.summary
    p.het     <- meta_obj$pooled$pval.het
    fdr       <- meta_obj$pooled$FDR
    g         <- meta_obj$summary$g
    lower_CI  <- meta_obj$pooled$lower_CI
    upper_CI  <- meta_obj$pooled$upper_CI
    
    # Apply thresholds dynamically using identify_robust_genes()
    identify_robust_genes(
      pool = pool,
      fdr = fdr,
      p.het = p.het,
      g = g,
      lower_CI = lower_CI,
      upper_CI = upper_CI,
      effect_thresh = input$effect_thresh,
      fdr_thresh = input$fdr_thresh,
      het_thresh = input$het_thresh
    )
  })
  
  # --------------------------
  # Update gene selectors based on filtered genes
  # --------------------------
  observe({
    req(results())
    
    genes <- filtered_genes()
    
    # Fallback if no robust genes
    if (length(genes) == 0) {
      genes <- results()$plot_data$all_genes
    }
    
    updateSelectInput(session,
                      "selected_gene",
                      choices = genes,
                      selected = genes[1])
    
    updateSelectInput(session,
                      "selected_gene_download",
                      choices = genes,
                      selected = genes[1])
  })
  
  
  # --------------------------
  # Display filtered results table
  # --------------------------
  output$results_table <- renderDT({
    req(results())
    
    table <- results()$table
    genes <- filtered_genes()
    
    if (length(genes) > 0) {
      table <- table[table$genes %in% genes, , drop = FALSE]
    }
    
    datatable(table, options = list(scrollX = TRUE))
  })
  
  
  # --------------------------
  # Dynamic plot UI
  # --------------------------
  output$gene_selector <- renderUI({
    req(results())
    
    genes <- filtered_genes()
    
    # Fallback if no robust genes
    if (length(genes) == 0) {
      genes <- results()$plot_data$all_genes
    }
    
    selectInput("selected_gene", "Select Gene:", choices = genes)
  })
  
  
  
  # --------------------------
  # Dynamic forest plot based on selected gene
  # --------------------------
  output$meta_plot_dynamic <- renderPlot({
    req(results())
    req(input$selected_gene)
    
    meta_obj <- results()$plot_data
    meta_res <- meta_obj$meta
    gene <- input$selected_gene
    
    # Only plot if gene is in the robust set
    if (!(gene %in% filtered_genes())) return(NULL)
    
    # Extract per-gene data
    g_gene <- meta_res$summary$g[gene, ]
    se.g_gene <- meta_res$summary$se.g[gene, ]
    pool <- meta_res$summary$pooled.estimates$summary[gene]
    se.pool <- meta_res$summary$pooled.estimates$se.summary[gene]
    study_names <- gsub("_g", "", names(g_gene))
    
    # Safety check
    if (is.null(g_gene) || length(g_gene) == 0) return(NULL)
    
    par(cex = 0.65)
    metaplot(
      g_gene, se.g_gene,
      labels = study_names,
      summn = pool,
      sumse = se.pool,
      sumnn = 1 / se.pool^2,
      summlabel = "Pooled Effect",
      xlab = "Standardized Mean Difference",
      ylab = "",
      xlim = c(-1.5, 1.5),
      main = gene,
      colors = meta.colors(
        box = "violetred",
        lines = "plum",
        summary = "mediumpurple",
        text = "black",
        axes = "black",
        zero = "black"
      ),
      cex = 0.75,
      boxsize = 1,
      lty.random = 1,
      lwd.random = 2,
      zero = 0,
      col.zero = "black",
      lty.zero = 3
    )
  }, width = 650, height = 550, res = 250)
  
  
  
  # --------------------------
  # Download Handler for Forest Plots
  # --------------------------
  output$download_plot_final <- downloadHandler(
    filename = function() {
      ext <- input$download_format
      paste0("meta_plot_", input$selected_gene_download, "_", Sys.Date(), ".", ext)
    },
    content = function(file) {
      
      req(results(), input$selected_gene_download)
      
      meta_obj <- results()$plot_data
      meta_res <- meta_obj$meta
      gene <- input$selected_gene_download
      
      # Extract per-gene data
      g_gene <- meta_res$summary$g[gene, ]
      se.g_gene <- meta_res$summary$se.g[gene, ]
      pool <- meta_res$summary$pooled.estimates$summary[gene]
      se.pool <- meta_res$summary$pooled.estimates$se.summary[gene]
      study_names <- gsub("_g", "", names(g_gene))
      
      # Safety check
      if (is.null(g_gene) || length(g_gene) == 0) return(NULL)
      
      # --------------------------
      # Generate plot depending on format
      # --------------------------
      if (input$download_format == "pdf") {
        pdf(file = file,
            width = input$download_width,
            height = input$download_height)
      } else if (input$download_format == "png") {
        png(file = file,
            width = input$download_width,
            height = input$download_height,
            units = "in",
            res = input$png_download_res)
      }
      
      # --------------------------
      # Plot using metaplot
      # --------------------------
      metaplot(
        g_gene, se.g_gene,
        labels = study_names,
        summn = pool,
        sumse = se.pool,
        sumnn = 1 / se.pool^2,
        summlabel = "Pooled Effect",
        xlab = "Standardized Mean Difference",
        ylab = "",
        xlim = c(-3, 3),
        main = gene,
        colors = meta.colors(
          box = "violetred",
          lines = "plum",
          summary = "mediumpurple",
          text = "black",
          axes = "black",
          zero = "black"
        ),
        cex = 0.75,
        boxsize = 1,
        lty.random = 1,
        lwd.random = 2,
        zero = 0,
        col.zero = "black",
        lty.zero = 3
      )
      
      dev.off()
    }
  )
  
  # --------------------------
  # Optional: Show info about the downloaded plot
  # --------------------------
  output$download_info <- renderPrint({
    req(input$selected_gene_download)
    cat("Gene:", input$selected_gene_download, "\n",
        "Format:", input$download_format, "\n",
        "Width:", input$download_width, "in\n",
        "Height:", input$download_height, "in\n",
        if (input$download_format == "png") paste0("DPI: ", input$png_download_res))
  })
}
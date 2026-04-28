find_common_genes <- function(
    DNA = FALSE, RNA = FALSE,
    list_of_dna_mtx = NULL,
    list_of_rna_mtx = NULL,
    use_DEG = FALSE,
    dna_deg_list = NULL,
    rna_deg_list = NULL,
    sig_label = "Significant"
  ) {
    
    # --- Helper: intersect genes across a list of objects ---
    intersection_from_lists <- function(lst) {
      Reduce(intersect, lst)
    }
    
    # ============================================================
    #       MODE 1: Use ONLY rownames (original behavior)
    # ============================================================
    if (!use_DEG) {
      
      if (DNA && !RNA) {
        if (!is.list(list_of_dna_mtx))
          stop("DNA matrices need to be in a list")
        
        return(intersection_from_lists(lapply(list_of_dna_mtx, rownames)))
      }
      
      if (RNA && !DNA) {
        if (!is.list(list_of_rna_mtx))
          stop("RNA matrices need to be in a list")
        
        return(intersection_from_lists(lapply(list_of_rna_mtx, rownames)))
      }
      
      if (DNA && RNA) {
        if (!is.list(list_of_dna_mtx) || !is.list(list_of_rna_mtx))
          stop("Both list_of_dna_mtx and list_of_rna_mtx must be lists of matrices.")
        
        list_of_mtx <- c(list_of_dna_mtx, list_of_rna_mtx)
        return(intersection_from_lists(lapply(list_of_mtx, rownames)))
      }
      
      return(NULL)
    }
    
    # ============================================================
    #       MODE 2: Use DE results ("Significance" == sig_label)
    # ============================================================
    if (use_DEG) {
      
      # --- DNA DEGs ---
      dna_sig_genes <- NULL
      if (DNA) {
        if (!is.list(dna_deg_list))
          stop("dna_deg_list must be a list of DEG data frames.")
        
        dna_sig_genes <- lapply(dna_deg_list, function(df) {
          if (!"Significance" %in% colnames(df))
            stop("DNA DEG result missing a 'Significance' column.")
          
          rownames(df)[df$Significance == sig_label]
        })
        dna_sig_genes <- intersection_from_lists(dna_sig_genes)
      }
      
      # --- RNA DEGs ---
      rna_sig_genes <- NULL
      if (RNA) {
        if (!is.list(rna_deg_list))
          stop("rna_deg_list must be a list of DEG data frames.")
        
        rna_sig_genes <- lapply(rna_deg_list, function(df) {
          if (!"Significance" %in% colnames(df))
            stop("RNA DEG result missing a 'Significance' column.")
          
          rownames(df)[df$Significance == sig_label]
        })
        rna_sig_genes <- intersection_from_lists(rna_sig_genes)
      }
      
      # Combine DNA + RNA intersections if both requested
      if (DNA && RNA) {
        return(intersect(dna_sig_genes, rna_sig_genes))
      }
      
      if (DNA) return(dna_sig_genes)
      if (RNA) return(rna_sig_genes)
    }
    
    return(NULL)
  }
  

library(shiny)
library(shinythemes)
library(DT)

ui <- fluidPage(
  bootstrapLib(bs_theme(preset = "sketchy")),  
  
  # --------------------------
  # Title
  # --------------------------
  titlePanel(
    div(
      "Meta-Multicohort Pipeline",
      style = "
      text-align: center; 
      font-size: 36px; 
      font-weight: bold; 
      color: #0d6efd; 
      text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
      margin-bottom: 20px;"
    )
  ),
  
  # --------------------------
  # Custom CSS for tabs
  # --------------------------
  tags$head(
    tags$style(HTML("
    /* Tab container underline */
    .nav-tabs {
      border-bottom: 3px solid black !important;  
      position: relative;
      z-index: 1;
    }

    /* Individual tab styling with black outline */
    .nav.nav-tabs > li a {
      border: 2px solid black;      
      border-radius: 0.5rem;        
      font-weight: bold;
      color: white !important;
      padding: 8px 15px;
      margin-right: 5px;
      position: relative;
      z-index: 2;
      background-clip: padding-box;
    }

    /* Individual tab colors */
    .nav.nav-tabs > li:nth-child(1) a { background-color: #0d6efd; }
    .nav.nav-tabs > li:nth-child(2) a { background-color: #FF69B4; }
    .nav.nav-tabs > li:nth-child(3) a { background-color: #198754; }
    .nav.nav-tabs > li:nth-child(4) a { background-color: #0dcaf0; color: black; }
    .nav.nav-tabs > li:nth-child(5) a { background-color: #ffc107; color: black; }

    /* Active tab slightly darker and shadow */
    .nav.nav-tabs > li.active a {
      box-shadow: 0 0 5px rgba(0,0,0,0.3);
    }

    /* Hover effect */
    .nav.nav-tabs > li a:hover {
      opacity: 0.85;
      text-decoration: none;
    }

    /* Tab content background */
    .tab-content > div {
      padding: 15px;
      border-radius: 8px;
      margin-top: 5px;
      background-color: #f8f9fa;
      position: relative;
      z-index: 0;
    }
  "))
  ),
  
  # --------------------------
  # Sidebar Layout
  # --------------------------
  sidebarLayout(
    
    # --------------------------
    # Sidebar Panel: Inputs, Welcome Box, Download
    # --------------------------
    sidebarPanel(
      
      # --------------------------
      # Welcome / Instructions Box
      # --------------------------
      div(
        style = "background-color: #fff3cd; border: 2px solid #ffc107; border-radius: 10px; padding: 15px; margin-bottom: 15px; box-shadow: 2px 2px 5px rgba(0,0,0,0.1);",
        h4("👋 Welcome!", style="color:#856404;"),
        p("This app allows you to load GEO datasets, process them, and run meta-analyses."),
        p("Follow these steps:"),
        tags$ol(
          tags$li("Enter GEO IDs and select study type (DNA/RNA)."),
          tags$li("Click 'Load Studies' to fetch and process data."),
          tags$li("Inspect and clean datasets if needed."),
          tags$li("Run Meta-Analysis and download plots."),
          tags$li("Explore genes using the dynamic plots.")
        ),
        p(strong("Tip:"), " Ensure your data columns are correctly labeled for Case/Control.")
      ),
      
      # --------------------------
      # Main Inputs
      # --------------------------
      textInput("geo_ids", "Enter GEO IDs (comma-separated):", value = ""),
      checkboxGroupInput("study_type", "Select Study Type:", choices=c("DNA","RNA"), selected="DNA"),
      actionButton("load_data", "Load Studies", class = "btn btn-primary"),
      br(), hr(),
      
      
      # --------------------------
      # Upload predefined gene set
      # --------------------------
      fileInput("gene_file", "Upload Gene List (optional)", accept = c(".txt", ".csv")),
      br(), hr(),
      
      
      # --------------------------
      # Robust gene thresholds
      # --------------------------
      div(
        style = "background-color: #f2f2f2; padding: 15px; border-radius: 8px; margin-bottom: 15px;",
        h4("⚙️ Robust Gene Thresholds", style = "color: #004080;"),
        
        numericInput(
          "effect_thresh", 
          "Summary effect size threshold (|pool| > ...)", 
          value = 0.5, min = 0, max = 5, step = 0.05
        ),
        
        numericInput(
          "fdr_thresh", 
          "FDR threshold (fdr < ...)", 
          value = 0.05, min = 0, max = 1, step = 0.01
        ),
        
        numericInput(
          "het_thresh", 
          "Heterogeneity p-value threshold (p.het > ...)", 
          value = 0.05, min = 0, max = 1, step = 0.01
        ),
        
        br(), hr()
      ),
    
      # --------------------------
      # Plot Download Section
      # --------------------------
      div(
        style = "background-color: #e6f2ff; padding: 15px; border-radius: 8px; box-shadow: 2px 2px 5px rgba(0,0,0,0.1); margin-top: 15px;",
        
        h4("📥 Plot Downloads", style = "color: #0059b3;"),
        
        # Gene selector (dynamic choices updated in server)
        selectInput("selected_gene_download", "Select Gene:", choices = character(0)),
        
        # Format selector
        radioButtons(
          "download_format", 
          "Format:", 
          choices = c("PDF" = "pdf", "PNG" = "png"),
          selected = "pdf",
          inline = TRUE
        ),
        
        # Width & Height
        numericInput("download_width", "Width (in)", value = 8, min = 3, max = 20, step = 0.5),
        numericInput("download_height", "Height (in)", value = 6, min = 3, max = 15, step = 0.5),
        
        # PNG DPI, only show if PNG is selected
        conditionalPanel(
          condition = "input.download_format == 'png'",
          numericInput("png_download_res", "DPI", value = 300, min = 72, max = 600, step = 50)
        ),
        
        # Download button
        downloadButton(
          "download_plot_final", 
          "⬇️ Download Plot", 
          class = "btn-primary btn-block", 
          style = "margin-top: 10px;"
        ),
        
        hr(),
        
        # Info about plot
        h5("Plot Info", style = "color: #004080;"),
        verbatimTextOutput("download_info", placeholder = TRUE)
      )
      
    ),
    
    # --------------------------
    # Main Panel: Tabs
    # --------------------------
    mainPanel(
      tabsetPanel(
        tabPanel("Inspect Expression Matrix",
                 selectInput("study_to_inspect_expr", "Select Study:", choices = NULL),
                 DTOutput("expr_structure_table")),
        
        tabPanel("Inspect pData",
                 selectInput("study_to_inspect", "Select Study:", choices = NULL),
                 DTOutput("pdata_structure_table")),
        
        tabPanel("Add Column",
                 selectInput("study_for_column","Select Study:", choices=NULL),
                 selectInput("column_name_input","Column to Convert to Binary:", choices=NULL),
                 DTOutput("column_preview_table"),
                 textInput("case_patterns_text","Values to assign as Case (regex, comma-separated):", value=""),
                 textInput("control_patterns_text","Values to assign as Control (regex, comma-separated):", value=""),
                 actionButton("convert_column","Add Condition Column")),
        
        tabPanel("Clean Studies",
                 selectInput("study_to_clean_single","Select Study:", choices=NULL),
                 verbatimTextOutput("clean_dims_before"),
                 checkboxGroupInput("condition_levels_keep","Condition levels to keep:", choices=c("Control","Case"), selected=c("Control","Case")),
                 checkboxInput("drop_na_condition","Remove samples with missing condition (NA)?", value=TRUE),
                 actionButton("apply_cleaning_single","Apply Cleaning"),
                 br(), hr(),
                 verbatimTextOutput("clean_dims_after")),
        
        tabPanel("Meta-Analysis",
                 actionButton("run_analysis","Run Meta-Analysis"),
                 br(), br(),
                 DTOutput("results_table"),
                 br(),
                 uiOutput("gene_selector"),
                 br(),
                 div(style="display:flex; justify-content:center;",
                     plotOutput("meta_plot_dynamic", width="400px", height="200px"))
        )
      )
    )
  )
)
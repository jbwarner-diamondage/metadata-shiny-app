library(shiny)
library(shinyjs)
library(shinydashboard)
library(readxl)
library(rhandsontable)
library(dplyr)
library(DT)

is_mainData_loaded <- reactiveVal(FALSE)  # Flag to track whether mainData is fully loaded
mainData <- reactiveVal(NULL)

HEADERS <- c(
  "ngs_sample_name",               
  "viral_plasmid_library_metadata",
  "comparison_group",
  "ngs_number",                    
  "sample_id",
  "sample_number",                 
  "sample_name",                   
  "cell_or_tissue_or_viral_or_plasmid",
  "barcode_or_nnk",
  "dna_or_rna",                    
  "insertion_size",
  "sample_type",
  "ligandscan_tissue_target",
  "sample_tissue",
  "capsid",
  "insertion_site",
  "restriction_site_desc",         
  "restriction_site_sequence",     
  "plasmid_library_associated_with_sample",
  "viral_library_associated_with_sample",
  "previously_sequenced_in",
  "amplicon",                      
  "guide",                         
  "i7_index_id",                   
  "index",                         
  "i5_index_id",                   
  "index2",
  "primer_scheme",
  "pcr_1_primer",
  "pcr_2_primer",                  
  "animal_id",
  "species"
)




initialize_editable_data <- function() {
  df <- data.frame(matrix("", ncol=length(HEADERS), nrow=30), stringsAsFactors=FALSE)  # Start with 20 rows
  colnames(df) <- HEADERS
  
  df$ngs_number <- "NGS_1"
  df$sample_id   <- ""
  df$sample_number <- paste0("S", seq_len(nrow(df)))  # Assign S1, S2, ..., S20
  df$sample_name <- paste0(df$sample_id, "_", df$sample_number)
  df$ngs_sample_name <- paste0(df$ngs_number, "_", df$sample_name)
  
  rownames(df) <- seq_len(nrow(df))  # Properly set row names
  
  df
}

search_folder <- "/cloud-data/digitalrnd-projects-ireland/BOWD/Magellan/CH_AAV_LAB/ngs_automation/samplesheets/historical_metadata_csvs"


SHEETS_TO_CLEAN <- c("AAV Preps","Tissues","Oligos","Plasmids")

search_in_files <- function(directory, search_text, session) {
  file_extensions <- c("csv", "csv")  # Include both CSV & XLSX
  file_pattern <- paste0("\\.(", paste(file_extensions, collapse = "|"), ")$")
  
  matching_files <- list.files(directory, pattern = file_pattern, full.names = TRUE)
  
  result_files <- c()  # Store matching files
  total_files <- length(matching_files)  # Get total number of files
  
  if (total_files == 0) {
    return(character(0))  # No files found
  }
  
  # Start progress
  withProgress(message = "Searching files...", value = 0, {
    for (i in seq_along(matching_files)) {
      file <- matching_files[i]
      
      incProgress(1 / total_files, detail = paste("Checking:", basename(file)))
      
      tryCatch({
        if (grepl("\\.csv$", file, ignore.case = TRUE)) {
          file_content <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
        } else if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
          file_content <- readxl::read_excel(file)
        }
        
        if (any(grepl(search_text, file_content, ignore.case = TRUE))) {
          result_files <- c(result_files, basename(file))
        }
        
      }, error = function(e) {
        message("Error reading file: ", file, " - ", e$message)
      })
    }
  })
  
  return(result_files)
}

clean_table <- function(df) {
  if (ncol(df)==0 || nrow(df)==0) return(df)
  col1 <- colnames(df)[1]
  df[[col1]] <- gsub("^CH_", "", df[[col1]])
  df[[col1]] <- gsub("[_ ]","-", df[[col1]])
  df <- df[df[[col1]]!="" & !is.na(df[[col1]]),,drop=FALSE]
  #df <- df[!duplicated(df[[col1]]),,drop=FALSE]
  #df <- df[order(df[[col1]]),,drop=FALSE]
  df
}

clean_primer_table <- function(df) {
  if (ncol(df)==0 || nrow(df)==0) return(df)
  col1 <- colnames(df)[1]
  df <- df %>%
    add_column(add_column = "constant_value")
  df <- df[!duplicated(df[[col1]]),,drop=FALSE]
  df <- df[order(df[[col1]]),,drop=FALSE]
  df
}

read_and_clean_sheet <- function(file_path, sheet_name) {
  tmp <- read_excel(file_path, sheet=sheet_name)
  if (sheet_name %in% SHEETS_TO_CLEAN) {
    tmp <- clean_table(tmp)
  }
  # Add 'primer_scheme' for 'Primers' table
  if (sheet_name == "Primers") {
    # Ensure required columns exist
    required_cols <- c("primer 1", "primer 2", "AAV cap", "insertion AA", "restriction site")
    missing_cols <- setdiff(required_cols, colnames(tmp))
    if (length(missing_cols) == 0) {
      # Create primer_scheme
      tmp$primer_scheme <- paste(
        tmp[["primer 1"]],
        "/",
        tmp[["primer 2"]],
        "_",
        tmp[["AAV cap"]],
        "_",
        tmp[["insertion AA"]],
        "_",
        tmp[["restriction site"]],
        sep = ""
      )
      
      # Move primer_scheme to the first column
      tmp <- tmp[, c("primer_scheme", setdiff(colnames(tmp), "primer_scheme"))]
      tmp <- tmp[order(tmp[["primer_scheme"]]),,drop=FALSE]
    } else {
      warning(paste("Missing columns in 'Primers' table:", paste(missing_cols, collapse = ", ")))
    }
  }
  
  tmp
}

# Function to load cached data if available, otherwise initialize
load_cached_data <- function() {
  shared_folder_path <- "/cloud-data/digitalrnd-projects-ireland/BOWD/Magellan/CH_AAV_LAB/ngs_automation/samplesheets/historical_metadata_csvs"
  
  cache_file <- file.path(shared_folder_path, "cached_metadata_do_not_delete.csv")
  
  
  if (file.exists(cache_file)) {
    message("🔄 Loading cached data from 'cached_metadata_do_not_delete.csv'...")
    
    df <- read.csv(cache_file, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
    
    # Replace all NA values with empty strings
    df[is.na(df)] <- ""
    
    # Ensure all expected columns exist
    missing_cols <- setdiff(HEADERS, colnames(df))
    for (col in missing_cols) {
      df[[col]] <- ""
      message(paste("⚠ Added missing column:", col))
    }
    
    # Reorder columns to match HEADERS
    df <- df[, HEADERS, drop = FALSE]
    return(df)
  } else {
    message("🚀 No cached data found. Initializing with default table...")
    return(initialize_editable_data())  # Default if no cache
  }
}

initialize_main_data <- function() {
  shared_folder_path <- "/cloud-data/digitalrnd-projects-ireland/BOWD/Magellan/CH_AAV_LAB/ngs_automation/samplesheets/historical_metadata_csvs"
  
  cache_file <- file.path(shared_folder_path, "cached_metadata_do_not_delete.csv")
  #"cached_metadata_test.csv"
  
  if (file.exists(cache_file)) {
    message("🔄 Loading cached data from 'cached_metadata_do_not_delete.csv'...")
    
    tryCatch({
      df <- read.csv(cache_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))
      df[is.na(df)] <- ""  # Convert NAs to empty strings
      
      message("✅ Successfully loaded mainData from cache.")
      
      # Mark that the data has been successfully loaded
      is_mainData_loaded(TRUE)
      
      return(df)
    }, error = function(e) {
      message("❌ Error loading cached data: ", e$message)
      is_mainData_loaded(TRUE)  # Mark as loaded even if failed, to prevent infinite loops
    })
  } 
  
  message("🆕 No cache found. Initializing fresh data...")
  
  # Ensure fresh data is initialized properly
  is_mainData_loaded(TRUE)  # Mark as loaded
  return(initialize_editable_data())
}


# ✅ Ensure mainData is set only **after** data is fully loaded
mainData <- reactiveVal(initialize_main_data())



# ✅ Load cached data **only once at startup**
observe({
  req(is_mainData_loaded())  # 🚀 Ensure data is loaded before auto-saving
  req(mainData())  # Ensure mainData exists before saving
  
  df <- isolate(mainData())  # Get current mainData
  df[is.na(df)] <- ""  # Replace NA values with empty strings
  
  message("💾 Auto-saving mainData to cache...")  # Debugging message
  
  tryCatch({
    shared_folder_path <- "/cloud-data/digitalrnd-projects-ireland/BOWD/Magellan/CH_AAV_LAB/ngs_automation/samplesheets/historical_metadata_csvs"
    
    write.csv(df, file.path(shared_folder_path, "cached_metadata_do_not_delete.csv"), row.names = FALSE)
    message("✅ Successfully saved mainData to 'cached_metadata_do_not_delete'.")
  }, error = function(e) {
    message("❌ Error saving cached data: ", e$message)
  })
})





server <- function(input, output, session) {
  
  shared_folder_path <- "/cloud-data/digitalrnd-projects-ireland/BOWD/Magellan/CH_AAV_LAB/ngs_automation/samplesheets/historical_metadata_csvs"
  #mainData <- reactiveVal(initialize_editable_data())
  mainData <- reactiveVal(load_cached_data())
  excel_data <- reactiveVal(NULL)
  
  rv_timer <- reactiveValues(
    prev = NULL, current = NULL
  )
  
  # Initialize file save system
  shinyFileSave(input, "save_as", roots = c(wd = "."), session = session)
  
  # Download handler for saving CSV
  output$save_as <- downloadHandler(
    filename = function() {
      paste("metadata_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(mainData())
      df <- mainData()
      #df <- df[, !colnames(df) %in% c("primer_scheme"), drop = FALSE]
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  predetermined_file <- "AAV_Key.xlsx"
  if (file.exists(predetermined_file)) {
    sh <- excel_sheets(predetermined_file)
    df_list <- setNames(
      lapply(sh, function(s) read_and_clean_sheet(predetermined_file, s)),
      sh
    )
    excel_data(df_list)
  }
  
  observeEvent(input$file1, {
    req(input$file1)
    file_path <- input$file1$datapath
    sh <- excel_sheets(file_path)
    df_list <- setNames(
      lapply(sh, function(s) read_and_clean_sheet(file_path, s)),
      sh
    )
    excel_data(df_list)
  })
  
  # Populate dropdown with available CSV files
  observe({
    csv_files <- list.files(shared_folder_path, pattern = "\\.csv$", full.names = FALSE)
    updateSelectInput(session, "selected_csv", choices = csv_files)
  })
  
  # Load selected CSV file when button is clicked
  observeEvent(input$load_server_csv, {
    req(input$selected_csv)
    
    file_path <- file.path(shared_folder_path, input$selected_csv)
    
    tryCatch({
      df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NAA"))
      df[is.na(df)] <- ""  # Replace NA with empty strings
      missing_cols <- setdiff(HEADERS, colnames(df))
      for (col in missing_cols) {
        df[[col]] <- ""
      }
      df <- df[, HEADERS, drop = FALSE]
      mainData(df)
      showModal(modalDialog(
        title = "File Loaded",
        paste("Successfully loaded:", input$selected_csv),
        easyClose = TRUE
      ))
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error Loading File",
        paste("Error loading file:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  output$file_info <- renderText({
    if (is.null(input$file1)) {
      if (!is.null(excel_data())) {
        "Oligo, Plasmid, AAV Key file loaded or data in memory."
      } else {
        "No file uploaded."
      }
    } else {
      paste("File uploaded:", input$file1$name)
    }
  })
  
  # Auto-save every 30 seconds (only after `mainData` is loaded)
  auto_save_interval <- 30  # Interval in seconds
  
  
  observe({
    req(is_mainData_loaded())  # 🚀 Prevent auto-saving before `mainData` is fully loaded
    req(mainData())  # Ensure `mainData` exists
    df <- mainData()  # Get current mainData
    invalidateLater(auto_save_interval * 1000)  # Convert seconds to milliseconds
    
    message("💾 Auto-saving 'cached_metadata_do_not_delete.csv'...")  # Debugging log
    shared_folder_path <- "/cloud-data/digitalrnd-projects-ireland/BOWD/Magellan/CH_AAV_LAB/ngs_automation/samplesheets/historical_metadata_csvs"
    
    cache_file <- file.path(shared_folder_path, "cached_metadata_do_not_delete.csv")
    write.csv(df, cache_file, row.names = FALSE)
    message("✅ Auto-saved 'cached_metadata_do_not_delete.csv' at", Sys.time())  # Debugging message
    
    #auto_save_data(df)  # Call the auto-save function
  })
  
  observeEvent(input$save_csv, {
    req(mainData())  # Ensure data exists
    req(input$save_filename)  # Ensure filename is entered
    
    # Define server-side folder where CSVs will be saved
    server_save_path <- "/cloud-data/digitalrnd-projects-ireland/BOWD/Magellan/CH_AAV_LAB/ngs_automation/samplesheets/historical_metadata_csvs"
    
    # Sanitize filename: remove special characters and ensure it ends with .csv
    filename <- gsub("[^A-Za-z0-9_-]", "_", input$save_filename)  
    filename <- paste0(filename, ".csv")
    
    # Full file path on server
    full_path <- file.path(server_save_path, filename)
    
    tryCatch({
      # Save the current mainData() to the specified file
      df <- mainData()
      write.csv(df, full_path, row.names = FALSE)
      
      showModal(modalDialog(
        title = "Success",
        paste0("✅ File saved successfully: ", server_save_path,"/",filename),
        easyClose = TRUE
      ))
      
      message("💾 File saved on server: ", full_path)  # Debugging message
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("❌ Failed to save file:", e$message),
        easyClose = TRUE
      ))
      message("❌ Error saving file: ", e$message)  # Debugging message
    })
  })
  
  output$save_as_info <- renderText("  Save metadata as csv")
  
 # observe({
 #   req(mainData())  # Ensure mainData is available
    
  #  message("📊 Updating UI with cached mainData...")
    
    # Force a refresh of Handsontable UI
 ##   session$sendCustomMessage("refreshHandsontable", list())
 # })
 # 
  # A dictionary for column-specific help
  columnHelp <- list(
    "ngs_sample_name"="Automatically created from <b>sample_name</b> and <b>ngs_number</b>.",
    "viral_plasmid_library_metadata"="Optional. Automatically created from sample tables, if possible. Edittable.",
    "comparison_group"="Optional. Entered by user.",
    "ngs_number"="Editable only in row 1. Must match NGS_# and be the same for entire table.",
    "sample_id"="<span style='color:red;'>Required.</span> Must match a sample in the Plasmids, AAV, or Tissues table or contain the words 'water' or 'Cell'.",
    "sample_number"="Automatically created from row number. <b>S1</b> for the sample in row 1, etc. ",
    "sample_name"="Automatically created from <b>sample_id</b> and <b>ngs_number</b>.",
    "cell_or_tissue_or_viral_or_plasmid"="Automatically created if <b>sample_id</b> in Plasmids/AAV/Tissues tables or sample_id contains the text 'Cell'.",
    "barcode_or_nnk"="<span style='color:red;'>Required.</span> Must be 'BARCODE', 'NNK', 'PROMOTER', or 'immunosilencing' ",
    "dna_or_rna"="Must be 'DNA' or 'RNA'.",
    "insertion_size"="<span style='color:red;'>Required for NNK or immunosilencing.</span> Entered by user. Must be a whole number.",
    "sample_type"="<span style='color:red;'>Required for NNK or immunosilencing.</span> For NNK samples, must be LigandScan or RandomScan.",
    "ligandscan_tissue_target"="<span style='color:red;'>Required for LigandScan</span> For LigandScan samples, must be </b>ocular</b>, <b>neuro</b>, <b>myo</b> or a valid <b>csv file</b>.",
    "sample_tissue"="Automatically created from <i>Tissues</i> table, if this is a Tissues sample. Otherwise, it is manually entered.",
    "capsid"="Automatically created from <b>primer_scheme</b>",
    "insertion_site"="Automatically created from <b>primer_scheme</b>",
    "restriction_site_desc"="Automatically created from <b>primer_scheme</b>",
    "restriction_site_sequence"="Automatically created from <b>primer_scheme</b>",
    "plasmid_library_associated_with_sample"="Optional. Entered by user.",
    "viral_library_associated_with_sample"="Optional. Entered by user.",
    "previously_sequenced_in"="Optional. Entered by user.",
    "amplicon"="Automatically created from <b>primer_scheme</b>",
    "guide"="Automatically created from <b>primer_scheme</b>",
    "i7_index_id"="<span style='color:red;'>Required.</span> Must match a known index on the <b>I7_index</b> table.",
    "index"="Automatically created from <b>i7_index_id</b> based on the <i><b>I7_Index</i></b> table.",
    "i5_index_id"="<span style='color:red;'>Required.</span> Must match a known index on the <b>I5_index</b> table.",
    "index2"="Automatically created from <b>i5_index_id</b> based on the <i><b>I5_Index</i></b> table.",
    "primer_scheme"="<span style='color:red;'>Required.</span> Select from a dropdown menu, based on the <i>Primers</i> table.",
    "pcr_1_primer"="Automatically created from <b>primer_scheme</b> based on the <i><b>Primers</i></b> table.",	
    "pcr_2_primer"="Automatically created from <b>primer_scheme</b> based on the <i><b>Primers</i></b> table.",	
    "animal_id"="Optional. Entered by user.",
    "species"="Optional. Entered by user."
    # ...
  )
  
  output$main_table <- renderRHandsontable({
    df <- mainData()
    if (is.null(df)) return(NULL)
    
    # Non-editable columns
    readOnlyCols <- c(
      "ngs_sample_name",
      "sample_number",
      "sample_name",
      "cell_or_tissue_or_viral_or_plasmid",
      "restriction_site_desc",
      "restriction_site_sequence",
      "capsid",
      "insertion_site",
      "amplicon",
      "guide",
      "index",
      "index2",
      "pcr_1_primer",
      "pcr_2_primer"
    )
    
    df_list <- excel_data()
    i7Choices <- c("")
    i5Choices <- c("")
    if (!is.null(df_list)) {
      if ("I7_index" %in% names(df_list)) {
        i7DF <- df_list[["I7_index"]]
        if (!is.null(i7DF) && "I7_Index_ID" %in% colnames(i7DF)) {
          i7Choices <- c("", unique(as.character(i7DF$I7_Index_ID)))
        }
      }
      if ("I5_index" %in% names(df_list)) {
        i5DF <- df_list[["I5_index"]]
        if (!is.null(i5DF) && "I5_Index_ID" %in% colnames(i5DF)) {
          i5Choices <- c("", unique(as.character(i5DF$I5_Index_ID)))
        }
      }
      if ("Primers" %in% names(df_list)) {
        primerDF <- df_list[["Primers"]]
        if (!is.null(primerDF) && "primer_scheme" %in% colnames(primerDF)) {
          primerChoices <- c("", unique(as.character(primerDF$primer_scheme)))
        }
      }
    }
    
    rh <- rhandsontable(df, rowHeaders=FALSE, height = 800) %>%
      hot_cols(manualColumnResize = TRUE, manualColumnMove = TRUE) %>%
      hot_table(contextMenu = TRUE, stretchH = "all") %>%
      hot_table(
        fixedColumnsLeft = 1,
        afterSelection = htmlwidgets::JS("
        function(r, c, r2, c2, preventScrolling) {
          Shiny.setInputValue('main_table_col_selected', c, {priority: 'event'});
          Shiny.setInputValue('main_table_rows_selected', {start: r + 1, end: r2 + 1}, {priority: 'event'});
          Shiny.setInputValue('main_table_row_selected', r + 1, {priority: 'event'});  // Convert to 1-based index
       
        }
      "),
        afterCreateRow = htmlwidgets::JS("
        function(index, amount) {
          Shiny.setInputValue('insert_row_index', {start: index, count: amount}, {priority: 'event'});
        }
      "),
        afterRemoveRow = htmlwidgets::JS("
        function(index, amount) {
          Shiny.setInputValue('delete_row_index', {start: index, count: amount}, {priority: 'event'});
        }
      ")
      )
    
 
  
    
    # Set custom width for primer_scheme
    if ("primer_scheme" %in% colnames(df)) {
      rh <- hot_col(rh, "primer_scheme", width = 250)  # Adjust width as needed
    }
    
    # Add conditional formatting to "sample_id"
    rh <- hot_col(rh, c("sample_id"), renderer = htmlwidgets::JS("
    function(instance, td, row, col, prop, value, cellProperties) {
      Handsontable.renderers.TextRenderer.apply(this, arguments);
      if (value === null || value === '') {
        td.style.background = '#FFC0CB'; // Pink
      }
    }
  "))
    
    rh <- hot_col(rh, c("i7_index_id","i5_index_id","primer_scheme","barcode_or_nnk"), renderer = htmlwidgets::JS("
    function(instance, td, row, col, prop, value, cellProperties) {
      Handsontable.renderers.DropdownRenderer.apply(this, arguments); // Retain dropdown functionality
      if (value === null || value === '') {
        td.style.background = '#FFC0CB'; // Pink for blank cells
      }
    }
  "))
    
    
    # Mark entire readOnly columns with myReadOnlyColumn class
    for (colName in readOnlyCols) {
      if (colName %in% colnames(df)) {
        rh <- hot_col(rh, colName, readOnly=TRUE, className="myReadOnlyColumn")
      }
    }
    
    # "ngs_number": mark entire column readOnly + color,
    # then re-enable row=1
    if ("ngs_number" %in% colnames(df)) {
      # Mark entire column readOnly
      rh <- hot_col(rh, "ngs_number", readOnly=TRUE, className="myReadOnlyColumn")
      # Re-enable row=1 => hot_cell does NOT accept className, so we omit
      if (nrow(df)>0) {
        rh <- hot_cell(rh, row=1, col="ngs_number", readOnly=FALSE)
      }
    }
    
    # dna_or_rna => must be '','DNA','RNA'
    rh <- hot_col(rh, "dna_or_rna",
                  type="dropdown",
                  source=c("","DNA","RNA"),
                  #className="dropdownColumn",
                  strict=TRUE,
                  allowInvalid=FALSE)
    # nnk_or_barcode => must be '','NNK','BARCODE','PROMOTER','immunosilencing'
    rh <- hot_col(rh, "barcode_or_nnk",
                  type="dropdown",
                  source=c("","NNK","BARCODE","PROMOTER","immunosilencing"),
                  #className="dropdownColumn",
                  strict=TRUE,
                  allowInvalid=FALSE)
    
    # Conditional formatting for "insertion_size"
    rh <- hot_col(rh, "insertion_size", 
                  renderer = htmlwidgets::JS("
                  function(instance, td, row, col, prop, value, cellProperties) {
                  Handsontable.renderers.TextRenderer.apply(this, arguments);
                   // var barcodeValue = instance.getDataAtCell(row, instance.propToCol('barcode_or_nnk'));
                  var barcodeValue = instance.getDataAtCell(row, col - 2); // Assuming barcode_or_nnk is one column left
                  if ((barcodeValue === 'NNK' || barcodeValue === 'immunosilencing') && (!value || value === '')) {
                    td.style.background = '#FFC0CB'; // Pink for empty insertion_size
                    } 
                  }
               "))
    
    # Conditional formatting for "ligandscan_tissue_target"
    rh <- hot_col(rh, "ligandscan_tissue_target", renderer = htmlwidgets::JS("
         function(instance, td, row, col, prop, value, cellProperties) {
         Handsontable.renderers.TextRenderer.apply(this, arguments);
         var sampleType = instance.getDataAtCell(row, col - 1); // Assuming sample_type is 1 column to the left
         if (sampleType === 'LigandScan') {
            var validValues = ['ocular', 'myo', 'neuro'];
            var isValidCSV = value && value.endsWith('.csv');
            if (!validValues.includes(value) && !isValidCSV) {
                td.style.background = '#FFC0CB'; // Highlight invalid cells in pink
              }
           }
         }
      "))
    
    rh <- hot_col(rh, "sample_type",
                  renderer = htmlwidgets::JS("
                function(instance, td, row, col, prop, value, cellProperties) {
                  Handsontable.renderers.DropdownRenderer.apply(this, arguments);
                  var barcodeValue = instance.getDataAtCell(row, col - 3); // Assuming barcode_or_nnk is three columns to the left
                  if ((barcodeValue === 'NNK' || barcodeValue === 'immunosilencing')) {
                    if (value !== 'LigandScan' && value !== 'RandomScan') {
                      td.style.background = '#FFC0CB'; // Pink for invalid values
                    }
                  }
                }
              "),
                  type = "dropdown",
                  source = c("", "LigandScan", "RandomScan"))
    
    # I7_Index_ID => i7Choices
    rh <- hot_col(rh, 
                  "i7_index_id",
                  type="dropdown",
                  source=i7Choices,
                  #className="dropdownColumn",
                  strict=TRUE,
                  allowInvalid=FALSE)
    
    # I5_Index_ID => i5Choices
    rh <- hot_col(rh, "i5_index_id",
                  type="dropdown",
                  source=i5Choices,
                  # className="dropdownColumn",
                  strict=TRUE,
                  allowInvalid=FALSE)
    
    # Primers => primerChoices
    rh <- hot_col(rh, "primer_scheme",
                  type="dropdown",
                  source=primerChoices,
                  # className="dropdownColumn",
                  strict=TRUE,
                  allowInvalid=FALSE)
    
    rh
  })
  
  output$info_main <- renderUI({
    HTML("Read-only columns are shaded <span style='color:gray;'>gray</span>, <br>Mandatory values are <span style='color:red;'>red</span> <br> Note: only Row 1 of ngs_number is editable.")
  })
  
  output$download_file <- downloadHandler(
    filename = function() {
      input$file_to_download  # Keep original filename
    },
    content = function(file) {
      req(input$file_to_download)  # Ensure a file is selected
      
      server_file_path <- file.path(search_folder, input$file_to_download)
      
      if (file.exists(server_file_path)) {
        file.copy(server_file_path, file)  # Copy the file to the user's download location
      } else {
        showModal(modalDialog(
          title = "Error",
          "Selected file does not exist on the server.",
          easyClose = TRUE
        ))
      }
    }
  )
  
  
  
  observeEvent(input$insert_row_index, {
    req(mainData())
    df <- mainData()
    
    # Parse the index of the inserted row
    insert_info <- input$insert_row_index
    start_index <- insert_info$start + 1  # Convert from 0-based to 1-based
    count <- insert_info$count
    
    # Validate the start index
    if (is.null(start_index) || start_index < 1 || start_index > (nrow(df) + 1)) {
      showModal(modalDialog(
        title = "Invalid Row Insertion",
        "The row index for insertion is invalid.",
        easyClose = TRUE
      ))
      return()
    }
    
    # Create new rows
    new_rows <- data.frame(matrix("", ncol = ncol(df), nrow = count), stringsAsFactors = FALSE)
    colnames(new_rows) <- colnames(df)
    
    # Handle special cases for inserting at the start
    if (start_index == 1) {
      # Insert rows at the very beginning
      df <- rbind(new_rows, df)
    } else if (start_index > nrow(df)) {
      # Append rows at the end
      df <- rbind(df, new_rows)
    } else {
      # Insert rows in the middle
      df <- rbind(
        df[1:(start_index - 1), , drop = FALSE],
        new_rows,
        df[start_index:nrow(df), , drop = FALSE]
      )
    }
    
    # Update sample_number and dependent columns
    df$sample_number <- paste0("S", seq_len(nrow(df)))
    if (nrow(df) > 1) {
      for (r in 2:nrow(df)) {
        df[r, "ngs_number"] <- df[1, "ngs_number"]
      }
    }
    for (r in seq_len(nrow(df))) {
      df[r, "sample_name"] <- paste0(df[r, "sample_id"], "_", df[r, "sample_number"])
      df[r, "ngs_sample_name"] <- paste0(df[r, "ngs_number"], "_", df[r, "sample_name"])
    }
    
    # Update the reactive value
    mainData(df)
  })
  
  
  
  # Observe column selection -> show help
  observeEvent(input$main_table_col_selected, {
    colIndex <- input$main_table_col_selected  # zero-based
    df <- mainData()
    if (!is.null(df)) {
      if (colIndex>=0 && colIndex<ncol(df)) {
        colName <- colnames(df)[colIndex+1]
        if (colName %in% names(columnHelp)) {
          output$column_help <- renderUI({
            HTML("<b>", colName, "</b>: ", columnHelp[[colName]])
          })
        } else {
          output$column_help <- renderText({
            HTML("No specific help for column <b>", colName, "</b>")
          })
        }
      } else {
        output$column_help <- renderText("No column selected.")
      }
    }
  })
  
  observeEvent(input$delete_row_index, {
    req(mainData())
    df <- mainData()
    
    # Parse the index of the removed row
    delete_info <- input$delete_row_index
    start_index <- delete_info$start + 1  # Convert from 0-based to 1-based
    count <- delete_info$count
    
    # Remove the rows from the data
    if (start_index > 0 && start_index <= nrow(df)) {
      rows_to_remove <- seq(start_index, length.out = count)
      df <- df[-rows_to_remove, , drop = FALSE]
      
      # Update sample_number and dependent columns
      df$sample_number <- paste0("S", seq_len(nrow(df)))
      if (nrow(df) > 1) {
        for (r in 2:nrow(df)) {
          df[r, "ngs_number"] <- df[1, "ngs_number"]
        }
      }
      for (r in seq_len(nrow(df))) {
        df[r, "sample_name"] <- paste0(df[r, "sample_id"], "_", df[r, "sample_number"])
        df[r, "ngs_sample_name"] <- paste0(df[r, "ngs_number"], "_", df[r, "sample_name"])
      }
      
      # Update the reactive value
      mainData(df)
    }
  })
  
  observeEvent(input$load_csv, {
    req(input$load_csv)  # Ensure a file is uploaded
    
    message("🚀 Starting CSV Load Process...")  # Debugging message
    
    file_path <- input$load_csv$datapath
    
    df <- tryCatch({
      read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error Loading CSV",
        paste("Could not read the CSV file:", e$message),
        easyClose = TRUE
      ))
      return(NULL)  # Stop execution
    })
    
    req(df)  # Ensure data was successfully read
    
    df[is.na(df)] <- ""  # Replace all NA values with empty strings
    
    # Ensure necessary columns exist
    missing_cols <- setdiff(HEADERS, colnames(df))
    for (col in missing_cols) {
      df[[col]] <- ""
      message(paste("⚠ Added missing column:", col))  # Debugging message
    }
    
    # Ensure only expected columns
    df <- df[, HEADERS, drop = FALSE]
    
    message("✅ CSV file loaded and validated.")  # Debugging message
    
    mainData(df)  # Update reactive variable
    
    showModal(modalDialog(
      title = "CSV Loaded Successfully",
      "The main table has been replaced with the new CSV data.",
      easyClose = TRUE
    ))
  })
  
  
  
  observeEvent(input$main_table, {
    df <- hot_to_r(input$main_table)
    if (is.null(df)) return()
    
    # If multiple rows, copy row1's ngs_number
    if (nrow(df)>1) {
      for (r in 2:nrow(df)) {
        df[r,"ngs_number"] <- df[1,"ngs_number"]
      }
    }
    
    # sample_id logic
    df_list <- excel_data()
    oligoIDs <- c()
    plasmidIDs <- c()
    aavIDs <- c()
    tissueIDs <- c()
    if (!is.null(df_list)) {
      if ("Oligos" %in% names(df_list)) {
        tmp <- df_list[["Oligos"]]
        if (ncol(tmp)>0) {
          col1 <- colnames(tmp)[1]
          oligoIDs <- unique(as.character(tmp[[col1]]))
        }
      }
      if ("Plasmids" %in% names(df_list)) {
        tmp <- df_list[["Plasmids"]]
        plasmidDF <- df_list[["Plasmids"]]
        if (ncol(tmp)>0) {
          col1 <- colnames(tmp)[1]
          plasmidIDs <- unique(as.character(tmp[[col1]]))
        }
      }
      if ("AAV Preps" %in% names(df_list)) {
        tmp <- df_list[["AAV Preps"]]
        aavDF <- df_list[["AAV Preps"]]
        if (ncol(tmp)>0) {
          col1 <- colnames(tmp)[1]
          aavIDs <- unique(as.character(tmp[[col1]]))
        }
      }
      if ("Tissues" %in% names(df_list)) {
        tmp <- df_list[["Tissues"]]
        tissueDF <- df_list[["Tissues"]]
        if (ncol(tmp)>0) {
          col1 <- colnames(tmp)[1]
          tissueIDs <- unique(as.character(tmp[[col1]]))
        }
      }
    }
    
    oldDF <- isolate(mainData())
    
    # Validation for ligandscan_tissue_target based on sample_type
    for (r in seq_len(nrow(df))) {
      if (df[r, "sample_type"] == "LigandScan") {
        valid_values <- c("ocular", "myo", "neuro", "")
        valid_csv <- grepl("\\.csv$", df[r, "ligandscan_tissue_target"], ignore.case = TRUE)
        if (!(df[r, "ligandscan_tissue_target"] %in% valid_values || valid_csv)) {
          # Revert to the previous value and show an error message
          df[r, "ligandscan_tissue_target"] <- oldDF[r, "ligandscan_tissue_target"]
          showModal(modalDialog(
            title = "Invalid Value",
            paste0(
              "Invalid value '", df[r, "ligandscan_tissue_target"],
              "' in 'ligandscan_tissue_target'. When 'sample_type' is 'LigandScan', ",
              "'ligandscan_tissue_target' must be 'ocular', 'myo', 'neuro', or end with '.csv'. ",
              "Reverting to the previous value."
            ),
            easyClose = TRUE
          ))
        }
      }
    }
    
    # Validation for ligandscan_tissue_target based on barcode_or_nnk
    for (r in seq_len(nrow(df))) {
      if (df[r, "barcode_or_nnk"] == "NNK") {
        valid_values <- c("LigandScan", "RandomScan","")
        if (!df[r, "sample_type"] %in% valid_values) {
          # Revert to the previous value and show an error message
          df[r, "sample_type"] <- oldDF[r, "sample_type"]
          showModal(modalDialog(
            title = "Invalid Value",
            paste0(
              "Invalid value '", df[r, "sample_type"],
              "' in 'sample_type'. When 'barcode_or_nnk' is 'NNK', ",
              "'sample_type' must be 'LigandScan' or 'RandomScan'. Reverting to the previous value."
            ),
            easyClose = TRUE
          ))
        }
      }
    }
    
    for (r in seq_len(nrow(df))) {
      sid <- df[r,"sample_id"]
      inOligos   <- sid %in% oligoIDs
      inPlasmids <- sid %in% plasmidIDs
      inAAV      <- sid %in% aavIDs
      inTissues  <- sid %in% tissueIDs
      
      changedSid <- !identical(sid, oldDF[r,"sample_id"])
      validSampleID <- (inOligos || inPlasmids || inAAV || inTissues ||
                          grepl("water", sid, ignore.case=TRUE) ||
                          grepl("Cell", sid, ignore.case=TRUE))
      
      if (changedSid && !validSampleID) {
        # revert
        df[r,"sample_id"] <- oldDF[r,"sample_id"]
        showModal(modalDialog(
          title="Invalid sample_id",
          paste("Value '", sid, "' not found in Oligos/Plasmids/AAV/Tissues nor has 'water'/'Cell'. Reverting."),
          easyClose=TRUE
        ))
        df[r,"cell_or_tissue_or_viral_or_plasmid"] <- oldDF[r,"cell_or_tissue_or_viral_or_plasmid"]
      } else {
        # auto-set cell_or_tissue_or_viral_or_plasmid
        if (inPlasmids) {
          df[r,"cell_or_tissue_or_viral_or_plasmid"] <- "Plasmid"
          
          df[r,"viral_plasmid_library_metadata"] <- plasmidDF[plasmidDF$`CH Lab Name`==df[r,"sample_id"],"Plasmid Description"]
        } else if (inAAV) {
          df[r,"cell_or_tissue_or_viral_or_plasmid"] <- "Viral"
          df[r,"viral_plasmid_library_metadata"] <- aavDF[aavDF$`CH Lab Name`==df[r,"sample_id"],"AAV Prep Description"]
        } else if (inTissues) {
          df[r,"cell_or_tissue_or_viral_or_plasmid"] <- "Tissue"
          df[r,"sample_tissue"] <- tissueDF[tissueDF$CH_TI_ID==df[r,"sample_id"],"Tissue received"]
          df[r,"viral_plasmid_library_metadata"] <- tissueDF[tissueDF$CH_TI_ID==df[r,"sample_id"],"Project"]
        } else if (grepl("Cell", sid, ignore.case=TRUE)) {
          df[r,"cell_or_tissue_or_viral_or_plasmid"] <- "Cell"
        }
      }
    }
    
    # Recompute sample_name, ngs_sample_name
    for (r in seq_len(nrow(df))) {
      df[r,"sample_name"] <- paste0(df[r,"sample_id"], "_", df[r,"sample_number"])
      df[r,"ngs_sample_name"] <- paste0(df[r,"ngs_number"], "_", df[r,"sample_name"])
    }
    
    # Auto-fill i7/i5 index
    if (!is.null(df_list)) {
      if ("I7_index" %in% names(df_list)) {
        i7DF <- df_list[["I7_index"]]
        if (!is.null(i7DF) && all(c("I7_Index_ID","index") %in% colnames(i7DF))) {
          for (r in seq_len(nrow(df))) {
            chosenID <- df[r,"i7_index_id"]
            rowMatch <- which(i7DF$I7_Index_ID==chosenID)
            if (length(rowMatch)>0) {
              df[r,"index"] <- as.character(i7DF[rowMatch[1],"index"])
            } else {
              df[r,"index"] <- ""
            }
          }
        }
      }
      if ("I5_index" %in% names(df_list)) {
        i5DF <- df_list[["I5_index"]]
        if (!is.null(i5DF) && all(c("I5_Index_ID","index2") %in% colnames(i5DF))) {
          for (r in seq_len(nrow(df))) {
            chosenID <- df[r,"i5_index_id"]
            rowMatch <- which(i5DF$I5_Index_ID==chosenID)
            if (length(rowMatch)>0) {
              df[r,"index2"] <- as.character(i5DF[rowMatch[1],"index2"])
            } else {
              df[r,"index2"] <- ""
            }
          }
        }
      }
      if ("Primers" %in% names(df_list)) {
        primerDF <- df_list[["Primers"]]
        if (!is.null(primerDF) && all(c("primer_scheme","primer 1","primer 2") %in% colnames(primerDF))) {
          for (r in seq_len(nrow(df))) {
            chosenID <- df[r,"primer_scheme"]
            print(chosenID)
            rowMatch <- which(primerDF$primer_scheme==chosenID)
            if (length(rowMatch)>0) {
              df[r,"pcr_1_primer"] <- as.character(primerDF[rowMatch[1],"primer 1"])
              df[r,"pcr_2_primer"] <- as.character(primerDF[rowMatch[1],"primer 2"])
              df[r,"amplicon"] <- as.character(primerDF[rowMatch[1],"amplicon"])
              df[r,"guide"] <- as.character(primerDF[rowMatch[1],"guide"])
              df[r,"restriction_site_desc"] <- as.character(primerDF[rowMatch[1],"restriction site"])   
              df[r,"restriction_site_sequence"] <- as.character(primerDF[rowMatch[1],"restriction site sequence"]) 
              df[r,"capsid"] <- as.character(primerDF[rowMatch[1],"AAV cap"]) 
              df[r,"insertion_site"] <- as.character(primerDF[rowMatch[1],"insertion AA"]) 
            } else {
              df[r,"pcr_1_primer"] <- ""
              df[r,"pcr_2_primer"] <- ""
              df[r,"amplicon"] <- ""
              df[r,"guide"] <- ""
              df[r,"restriction_site_desc"] <- ""
              df[r,"restriction_site_sequence"] <- "" 
              df[r,"capsid"] <- ""
              df[r,"insertion_site"] <- ""
            }
          }
        }
      }
    }
    
    mainData(df)
  })
  
  ############################
  # Add row
  ############################
  # Reactive trigger for adding rows
  # Reactive trigger for adding rows
  add_row_trigger <- reactiveVal(NULL)
  
  observeEvent(input$add_row, {
    if (is.null(rv_timer$prev)) {
      rv_timer$prev <- Sys.time()
      #return(NULL)
    }
    
    # if it's not the first time to edit table, get current clock time:
    rv_timer$current <- Sys.time()
    
    # if the difference btwn prev recorded time and current time is less 
    # than 1second, don't do anything, just return:
    if ((rv_timer$current - rv_timer$prev) < 1.0) {
      return(NULL)
    }
    
    
    
    # finally set current clock time as `rv_timer$prev` for use in the next 
    # invalidation:
    
    req(mainData())  # Ensure table exists
    shinyjs::disable("add_row")  # Disable button to prevent spam clicks
    
    isolate({
      df <- mainData()  # Get the current table
      
      # Determine the selected row; default to last row if none selected
      selected_row <- input$main_table_row_selected
      if (is.null(selected_row) || selected_row < 1 || selected_row > nrow(df)) {
        selected_row <- nrow(df)
      }
      
      # Create a new empty row
      new_row <- data.frame(matrix("", ncol = ncol(df), nrow = 1), stringsAsFactors = FALSE)
      colnames(new_row) <- colnames(df)
      
      # **Insert the new row BELOW the selected row**
      if (selected_row == nrow(df)) {
        df <- rbind(df, new_row)  # Append at the end
      } else {
        df <- rbind(
          df[1:selected_row, , drop = FALSE],  # Rows up to and including the selected row
          new_row,                             # Insert new row
          df[(selected_row + 1):nrow(df), , drop = FALSE]  # Remaining rows
        )
      }
      
      # **Ensure new row is empty before validation runs**
      df[selected_row + 1, ] <- ""
      
      # **Update sample_number for all rows**
      df$sample_number <- paste0("S", seq_len(nrow(df)))
      
      # **Only update sample_name & ngs_sample_name for non-empty sample_id rows**
      for (r in seq_len(nrow(df))) {
        if (nzchar(df[r, "sample_id"])) {  # Skip blank rows
          df[r, "sample_name"] <- paste0(df[r, "sample_id"], "_", df[r, "sample_number"])
          df[r, "ngs_sample_name"] <- paste0(df[r, "ngs_number"], "_", df[r, "sample_name"])
        }
      }
      
      # **Freeze row selection to prevent flickering**
      freezeReactiveValue(input, "main_table_row_selected")
      updateNumericInput(session, "main_table_row_selected", value = selected_row + 1)
      
      rv_timer$prev <- Sys.time()
      # **Prevent Handsontable from resetting focus**
      isolate({
        mainData(df)
      })
    })
    
    # **Debounce re-enabling of button**
    delay(300, { shinyjs::enable("add_row") })
  })
  
  observe({
    file_extensions <- c("csv", "xlsx")  # Include both file types
    file_pattern <- paste0("\\.(", paste(file_extensions, collapse = "|"), ")$")
    
    available_files <- list.files(search_folder, pattern = file_pattern, full.names = FALSE)
    
    updateSelectInput(session, "file_to_download", choices = available_files)
    # Enable button if files exist
    if (length(available_files) > 0) {
      shinyjs::enable("download_file")
    } else {
      shinyjs::disable("download_file")
    }
  })
  
  observeEvent(input$search_button, {
    req(input$search_text)
    
    output$search_progress <- renderText("🔍 Searching... Please wait.")
    
    found_files <- search_in_files(search_folder, input$search_text, session)
    
    output$search_results <- renderUI({
      if (length(found_files) > 0) {
        HTML(paste("<SPAN STYLE='color:#000000'><b>✅ Found in files:</b><br>", paste(found_files, collapse = "<br>")))
      } else {
        HTML("<SPAN STYLE='color:#000000'>❌ No files found with the given search text.")
      }
    })
    
    # Clear progress message once complete
    output$search_progress <- renderText("")
  })
  
  
  
  # Handle row addition logic (debounced ob)
  observeEvent(add_row_trigger(), {
    req(add_row_trigger())
    
    df <- mainData()  # Get the current data
    
    # Determine the selected row to add a new row below
    selected_row <- input$main_table_row_selected
    if (is.null(selected_row) || selected_row < 1 || selected_row > nrow(df)) {
      selected_row <- nrow(df)  # Default to the last row if no valid selection
    }
    
    # Create the new row
    new_row <- data.frame(matrix("", ncol = ncol(df), nrow = 1), stringsAsFactors = FALSE)
    colnames(new_row) <- colnames(df)
    
    # Insert the new row below the selected row
    if (selected_row == nrow(df)) {
      df <- rbind(df, new_row)  # Append at the end
    } else {
      df <- rbind(
        df[1:selected_row, , drop = FALSE],  # Rows up to and including the selected row
        new_row,                             # New row inserted after the selected row
        df[(selected_row + 1):nrow(df), , drop = FALSE]  # Remaining rows
      )
    }
    
    # Update dependent columns
    df$sample_number <- paste0("S", seq_len(nrow(df)))
    if (nrow(df) > 1) {
      for (r in 2:nrow(df)) {
        df[r, "ngs_number"] <- df[1, "ngs_number"]
      }
    }
    for (r in seq_len(nrow(df))) {
      df[r, "sample_name"] <- paste0(df[r, "sample_id"], "_", df[r, "sample_number"])
      df[r, "ngs_sample_name"] <- paste0(df[r, "ngs_number"], "_", df[r, "sample_name"])
    }
    
    mainData(df)  # Update the reactive data
    
    shinyjs::enable("add_row")  # Re-enable the button
  }, ignoreInit = TRUE)  # Prevent immediate triggering on app load
  
  
  
  
  ############################
  # Delete row 
  ############################
  observeEvent(input$delete_row, {
    if (is.null(rv_timer$prev)) {
      rv_timer$prev <- Sys.time()-2
      #return(NULL)
    }
    
    # if it's not the first time to edit table, get current clock time:
    rv_timer$current <- Sys.time()
    
    # if the difference btwn prev recorded time and current time is less 
    # than 1second, don't do anything, just return:
    if ((rv_timer$current - rv_timer$prev) < 1.0) {
      return(NULL)
    }
    shinyjs::disable("delete_row")  # Disable button to prevent spam clicks
    req(mainData())
    df <- mainData()
    
    
    selected_rows <- input$main_table_rows_selected
    
    # Validate selected_rows
    if (is.null(selected_rows) || selected_rows$start < 1 || selected_rows$end > nrow(df)) {
      showModal(modalDialog(
        title = "Invalid Row Selection",
        "Please select valid rows to delete.",
        easyClose = TRUE
      ))
      return()
    }
    
    # Get range of rows to delete
    rows_to_delete <- seq(selected_rows$start, selected_rows$end)
    
    # Remove the selected rows
    df <- df[-rows_to_delete, , drop = FALSE]
    
    # Update dependent columns
    if (nrow(df) > 0) {
      df$sample_number <- paste0("S", seq_len(nrow(df)))
      if (nrow(df) > 1) {
        for (r in 2:nrow(df)) {
          df[r, "ngs_number"] <- df[1, "ngs_number"]
        }
      }
      for (r in seq_len(nrow(df))) {
        df[r, "sample_name"] <- paste0(df[r, "sample_id"], "_", df[r, "sample_number"])
        df[r, "ngs_sample_name"] <- paste0(df[r, "ngs_number"], "_", df[r, "sample_name"])
      }
    }
    # **Freeze row selection to prevent flickering**
    freezeReactiveValue(input, "main_table_rows_selected")
    updateNumericInput(session, "main_table_rows_selected", value = selected_rows$end + 1)
    
    rv_timer$prev <- Sys.time()
    # **Prevent Handsontable from resetting focus**
    isolate({
      mainData(df)
    })
    # **Debounce re-enabling of button**
    delay(300, { shinyjs::enable("delete_row") })

  })
  
  
  
  
  
  ############################
  # TABS UI
  ############################
  output$tabs_ui <- renderUI({
    df_list <- excel_data()
    if (is.null(df_list)) {
      return(span("No Excel data loaded."))
    }
    sheet_names <- names(df_list)
    tabs <- lapply(sheet_names, function(shName){
      tabPanel(
        shName,
        rHandsontableOutput(paste0("sheet_", shName))
      )
    })
    do.call(tabsetPanel, tabs)
  })
  
  observeEvent(excel_data(), {
    df_list <- excel_data()
    if (is.null(df_list)) return()
    sheet_names <- names(df_list)
    for (shName in sheet_names) {
      local({
        sheetTitle <- shName
        output[[paste0("sheet_", sheetTitle)]] <- renderRHandsontable({
          df_sheet <- df_list[[sheetTitle]]
          if (is.null(df_sheet)) return(NULL)
          df_sheet[] <- lapply(df_sheet, as.character)
          rh <- rhandsontable(df_sheet, readOnly=TRUE, rowHeaders=FALSE)
          rh %>% hot_cols(manualColumnMove=TRUE, manualColumnResize=TRUE)
        })
      })
    }
  })
}



library(shiny)
library(shinyjs)
library(shinydashboard)
library(readxl)
library(rhandsontable)
library(dplyr)
library(shinyFiles)

#/cloud-data/digitalrnd-labdata-virginia/MA-Lab/Genomics/Illumina/NextSeq2000_

useShinyjs()  # Enable shinyjs for managing UI elements

ui <- dashboardPage(
  dashboardHeader(title="Metadata Maker v0.103",titleWidth=350),
  dashboardSidebar(
    
    
    width=350,
    fileInput("file1","Choose AAV Key File",accept=".xlsx"),
    textOutput("file_info"),
    fileInput("load_csv", "Load a Local Metadata CSV", accept = c(".csv")),  # NEW CSV UPLOAD
    selectInput("selected_csv","Choose a Shared Metadata CSV",choices=NULL),
    actionButton("load_server_csv", "Load Shared Metadata CSV"),
    br(),
    
    tags$style(HTML("
      #file_info {
        margin-left: 15px; /* Move the text 15 pixels to the right */
      }
    ")),
    
    br(),
    actionButton("add_row", "Add Row", class = "btn-primary", disabled = FALSE),
    actionButton("delete_row","Delete Row",class = "btn-primary", disabled = FALSE),
    br(),
    
    tags$style(HTML("
      #save_as_info {
        margin-left: 15px; /* Move the button 15 pixels to the right */
      }
    ")),
    textOutput("save_as_info"),
    br(),
    
    tags$style(HTML("
      #save_as {
        margin-left: 15px; /* Move the button 15 pixels to the right */
      }
    ")),
    
    textInput("save_filename", "Enter Filename for CSV to Save:", value = paste0("metadata_", Sys.Date())),
    actionButton("save_csv", "Save Metadata CSV to Server"),
    box(
      title = "🔎 Search for Text in Metadata Files",
      status = "warning",
      solidHeader = TRUE,
      width = 12,
      
      textInput("search_text", "Enter Search Term:", ""),
      actionButton("search_button", "Search Files"),
      
      br(),
      tags$hr(),
      verbatimTextOutput("search_progress"),  # Show progress
      uiOutput("search_results")  # Show results
    ),
    
    box(
      title = "📥 Download File from Server", status = "info", solidHeader = TRUE, width = 12,
      selectInput("file_to_download", "Select a file:", choices = NULL),
      downloadButton("download_file", label="Download Selected File", class = "butt"),
      tags$head(tags$style(".butt{background-color:#FFFFFF;} .butt{color: black;} .skin-blue .sidebar a { color: #444; }")), # background color and font color
      
    )
    
    
    
    
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        /* Enforce consistent row height for Handsontable */
        .handsontable td, .handsontable th {
          white-space: nowrap !important;  /* Prevent text from wrapping */
          overflow: hidden !important;     /* Hide overflowing text */
          text-overflow: ellipsis !important;  /* Add '...' if text overflows */
          height: 30px !important;         /* Set a fixed row height */
          line-height: 30px !important;     /* Ensure text is vertically centered */
        }
        
        /* Prevent column height mismatch when scrolling */
        .ht_master table, .ht_clone_left table {
          table-layout: fixed !important;  /* Prevents auto row height expansion */
        }
        .dropdownColumn {
          background-color: #ffe0b3 !important; /* peach for dropdowns */
        }
        .myReadOnlyColumn {
          background-color: #eeeeee !important; /* light gray for readOnly columns */
        }
        /* Fix overflowing text in selectInput */
        .selectize-input {
          text-overflow: ellipsis !important;  /* Adds '...' for long text */
          overflow: hidden !important;         /* Hides overflowing text */
          white-space: nowrap !important;      /* Prevents wrapping */
          display: block !important;           /* Ensures proper behavior */
          max-width: 100% !important;          /* Ensures it does not stretch */
        }
      "))
    ),
    
    fluidRow(
      box(
        title="Metadata Table",
        status="primary",
        solidHeader=TRUE,
        width=12,
        rHandsontableOutput("main_table", width = "100%", height = "800px"),
        br(),
        uiOutput("info_main"),
        br(),
        uiOutput("column_help")
      )
    ),
    
    fluidRow(
      box(
        title="Oligo, Plasmid, AAV Key Tables",
        status="info",
        solidHeader=TRUE,
        width=12,
        uiOutput("tabs_ui")
      )
    )
  )
)

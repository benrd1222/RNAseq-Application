library(shiny)
library(bslib)
library(DESeq2)

# This needs to be carefully considered if you want to host this app online,
# but for local use it is dependent on your computers available memory
options(shiny.maxRequestSize = 4000 * 1024^2) #set to 4 gigabytes

# NEXT STEPS:
# Address the various TODOs and get a version 0.0.1 working

# UI ----
ui <- page_fluid(
  input_dark_mode(id = "mode"),
  navset_pill(
    # Panel for data set up and initial analysis ----
    nav_panel("Set Up",
              layout_columns(
                card(
                  card_header("Data Import"),
                  
                  fileInput("data",
                            "Select the CSV file of the raw counts for this analysis",
                            accept = c("text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")
                  ), # may want to support other file types
                  fileInput("meta",
                            "Select the CSV file of the meta data for this analysis",
                            accept = c("text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")
                  ),
                  "Note: Make sure your metadata only includes the sample names and
                the treatment variables you list in your formula. The example text 
                in the formula field below would be valid for a metadata structure 
                of one column labeled condition that notes the treatment of each 
                sample as control vs. experimental"
                ),
                card(
                  card_header("Analysis"),
                  textInput("formula", 
                            "Please enter a formula valid for DESeq2", 
                            value = "~ condition"
                  ),
                  # do I make them write the tilde or just imply it behind the scenes
                  # TODO: more helpfully I could update the available variable names 
                  #when they upload the metadata
                  input_switch("relevel_switch", "
                               Define reference levels"),
                  uiOutput("reveal_relevel"),
                  input_switch("pre_filter_switch", 
                               "Apply Pre-filtering"),
                  uiOutput("reveal_pre_filter"),
                  input_switch("shrink_switch", 
                               "Apply LFCshrinkage"),
                  # TODO: this switch also needs to bring up a drop down of the
                  # different coefficients that can be shrunk
                  tags$hr(),
                  actionButton("run_dds", 
                               "Run Analysis")
                )
              ),
              
              card(
                card_header("Analysis Export"),
                downloadButton("export_dds", 
                             "Download DESeq Datastructure")
              )
              ),
    
    # Overview of the data ----
    # accordion panels with a few different tools collapsed inside them
    nav_panel("Basics", 
              "Here we allow the user to choose their reference condition, and
              explore the various contrasts available as well as export data. Given
              a reference level and the listed contrasts of interest we can also 
              generate boxplots of the LFC for a list of up to 5 genes of interest.
              We would also like to make an option for them to look at the boxplots
              for the raw counts per million."
              ), 
    
    # Ontology ----
    nav_panel("Ontology", 
              "Here we can output a visualization of the ontology for the chosen
              contrast from a drop down menu. There should also be a few button
              to export. This could also include KEGG pathways"
              ),
    # Extras ----
    nav_menu( 
      "Other links", 
      nav_panel("Venn Diagrams *Beta*", 
                "May break"
                ), 
      "----",
      "Help",
      nav_item( 
        a("DESeq2 Tutorial", 
          href = "https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html",
          target = "_blank") 
      ),
      nav_item( 
        a("Valid Formulas", 
          href = "https://www.datacamp.com/tutorial/r-formula-tutorial",
          target = "_blank") 
      )
    ), 
  ), 
  id = "tab" 
)

# Server ----
server <- function(input, output) {
  # Conditional Switch inputs ----
  # Allowing user re-leveling
  output$reveal_relevel <- renderUI({
    meta<-metaData()
    if (input$relevel_switch && !is.null(meta)) {
      textInput("Placeholder","Placeholder") 
      # TODO: Dynamic UI element, based on the amount of treatment variables
      # listed in the metadata allows for the selection of reference level for each
    } else {
      NULL
    }
  })
  
  # Allowing user pre-filtering preferences 
  output$reveal_pre_filter <- renderUI({
    if(input$pre_filter_switch){
      tagList(
        numericInput( 
          "smallestCount", 
          "Smallest raw counts to keep", 
          value = 10, 
          min = 1, 
          max = 100 
        ),
        numericInput( 
          "smallestGroup", 
          "Smallest group cutoff", 
          value = 1, 
          min = 1, 
          max = 10 
        )
      )
    }
    else{
      NULL
    }
  })
  # Data Input ------------------------------------------------
  # Reactive observation of the data or metadata:
  rawData <- reactive({
    req(input$data) # Requires that a file has been uploaded
    
    df <- tryCatch(
      {
        read.csv(input$data$datapath,
                 header = input$header,
                 sep = input$sep,
                 stringsAsFactors = FALSE)
      },
      error = function(e) {
        NULL
      }
    )
    return(df)
  })
  
  metaData <- reactive({
    req(input$meta) # Requires that a file has been uploaded
    
    df <- tryCatch(
      {
        read.csv(input$data$datapath,
                 header = input$header,
                 sep = input$sep,
                 stringsAsFactors = TRUE)
      },
      error = function(e) {
        NULL
      }
    )
    return(df)
  })
  
  usrFormula <- reactive({
    req(input$formula)
    
    f <- tryCatch(
      {
        formula(f)
      },
      error = function(e){
        NULL
      }
    )
    return(f)
  })
  # Initial Analysis Run ------------------------------------------------
  analysisOutput <- eventReactive(input$run_dds, {

    counts <- rawData()
    meta <- metaData()
    f <- usrFormula()
    req(counts,meta,f)
    
    # TODO: design validation, check scratch.R for info on formula validation
    validate(
      need(
        !is.null(counts) && nrow(counts) > 0,
        "No data has been uploaded."
      ),
      need(
        !is.null(meta) && nrow(meta) > 0,
        "No metadata has been uploaded."
      ),
      need(
        !is.null(f),
        "Invalid formula provided."
      ),
      need(
        "The first columns of the data to be a specfic name and be of type character"
           ),
      need(
        "colnames to match the first column of the metadata"
      ),
      need(
        "check that meta[,1]==colnames(counts[,-c(1,2)])"
      ),
      need(
        "columns in metadata to conform with the pieces of valid formula"
      )
    )
    
    # Validation should have us proceeding with the following datastructures
    # counts <-
    # Gene_id | External_gene_name | Sample1_name | ... | SampleN_name
    # <chr>   | <chr>              | <numeric>    | ... | <numeric>
    
    # meta <- 
    # Sample_names | Treatment1 | ... | TreatmentN
    # <factor>     | <factor>   | ... | <factor>

    # f <- valid formula
    
    # TODO: something is broken with the progress run, or my validation is not
    # informative enough. I would guess it is an issue with how I am passing my 
    # forumula and the reactive check. I believ we are failing on line 205 calling
    # to the function on line 187
    
    withProgress(meassage = "Running DESeq2 Analysis", value=0,{
      # body of initial analysis
      genes<-cbind(counts[,1],counts[,2])

      row.names(counts)<-counts[,1]
      counts<-counts[,-c(1,2)]
      
      row.names(meta)<-meta[,1]
      meta<-meta[,-1]
      
      if(input$relevel_switch){
        #Do some re leveling based off of the dynamic re leveling solution
      }
      
      dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta, design=f)
      
      if(input$pre_filter_switch){
        smallestGroupSize <- input$smallestGroup
        smallestCounts <- input$smallestCount
        dds <- dds[rowSums(counts(dds) >= smallestCount) >= smallestGroupSize, ]
      }
      
      if(input$shrink_switch){
        # DOES NOTHING FOR NOW
        #lfcShrink(dds, coef="fix needed", type="apeglm")
      }

      dds <- DESeq(dds)
    })
    
    dds
  }) # end analysisOutput
  
  # Clear analysis output if file changes after an analysis was run ----
  observeEvent(c(input$data,input$meta,input$formula), {
    output$analysisOutput <- renderPrint("")
  }, ignoreInit = TRUE)

  # dds download handler -----
  output$export_dds <- downloadHandler(
    filename = function() {
      paste0("DESeq_app_dds_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(analysisOutput) # this statement may need to point to reactive funciton
      
      saveRDS(analysisOuput, file = file)
    }
  ) #TODO: check that this works properly with small example dataset
  
  
}

shinyApp(ui = ui, server = server)









library(shiny)
library(bslib)

# tags$hr(), provides a line for distinciton within a section of the webpage

ui <- page_fluid(
  navset_pill( 
    nav_panel("Set Up",
              card(
                card_header("Data Import:"),
                fileInput("data",
                          "Select the CSV file of the raw counts for this analysis",
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
                          ), # may need to support other file types
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
                sample as control vs. experimental",
                textInput("formula", 
                          "Please enter a formula valid for DESeq2", 
                          value = "~ condition" 
                          # do I make them write the tilde or just imply it behind the scenes
                          ),
                tags$hr(),
                actionButton("run", "Run Analysis")
              ),
              ),
    
    nav_panel("Basics", 
              "Page B content"
              ), 
    nav_panel("Ontology", 
              "Page C content"
              ),
    nav_menu( 
      "Other links", 
      nav_panel("Venn Diagrams *Beta*", 
                "May break"
                ), 
      "----", 
      nav_item( 
        a("Tutorial", 
          href = "https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html",
          target = "_blank") 
      ), 
    ), 
  ), 
  id = "tab" 
)

server <- function(input, output) {
  # Data Input ------------------------------------------------
  # Observe the changing of the data or metadata:
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
  dds <- eventReactive(input$runAnalysis, {

    data <- rawData()
    meta <- metaData()
    f <- usrFormula()
    req(data,meta,f)
    
    validate(
      need(
        !is.null(data) && nrow(data) > 0,
        "No Data has been uploaded."
      ),
      need(
        !is.null(data) && nrow(data) > 0,
        "No Data has been uploaded."
      ),
      need(
        "The first columns of the data to be a specfic name and be of type character"
           ),
      need(
        "colnames to match the first column of the metadata"
      ),
      need(
        "some check on metadata"
      ),
      need(
        "columns in metadata to conform with the pieces of valid formula"
      )
    )
    
    # Example Analysis:
    mean_val <- mean(data$Value, na.rm = TRUE)
    sd_val <- sd(data$Value, na.rm = TRUE)
    n_obs <- nrow(data)
    
    # return a list of useful things on complete analysis. I believe I will really
    # only need to return the dds object
    list(
      mean = mean_val,
      sd = sd_val,
      n = n_obs,
      message = "Analysis completed successfully!"
    )
  })
  
  # withProgress(meassage = "Running DESeq2 Analysis", value=0,{
  #   # body of initial analysis
  #   genes<-cbind(counts[,1],counts[,2])
  #   
  #   row.names(counts)<-counts[,1]
  #   counts<-counts[,-c(1,2)]
  #   row.names(counts)<-strip(row.names(counts))
  #   
  #   
  # })

  
  # Optional: Clear analysis output if file changes after an analysis was run
  observeEvent(c(input$data,input$meta,input$formula), {
    output$analysisOutput <- renderPrint("")
    output$validationMessage <- renderText("")
  }, ignoreInit = TRUE) # Don't run on initial app load
  
  
}

shinyApp(ui = ui, server = server)


# useful bits and bobs

# If you include a message with the output of a validated function in the server
# side, you can render the message for the user to inform them if there is an error
# output$validationMessage <- renderText({
#   results <- analysisResults()
#   if (!is.null(results)) {
#     return(results$message)
#   }
#   return("")
# })





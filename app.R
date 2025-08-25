library(shiny)
library(bslib)
library(DESeq2)
library(stringr)
library(ggplot2)
library(ggrepel)

# TODO: think about how to install these packages for the user, or make that part
# of the markdown tutorial that will go along with this app

# This needs to be carefully considered if you want to host this app online,
# but for local use it is dependent on your computers available memory
options(shiny.maxRequestSize = 4000 * 1024^2) #set to 4 gigabytes

# TODO: MED- replace the horrible cat errors and conditional UI error
# with the showNotification, just found out this is an option, and it is exactly
# the option I was looking for

# showNotification(
#   "Error: The uploaded file is not a valid 'DESeqDataSet' object.",
#   type = "error"
# )
#
# The premise here is to have a tryCatch around the code of interest, then catch
# errors and return NULL. Have an if(NULL){} check that leads to the showNotification

# TODO: LOW- set it up so the cards are different widths for Set Up and Analysis
# make the SetUp 1/3 of the page and analysis taking up the rest
# could also look to get the dynamic listing of covariates to choose reference levels
# for renders them all side by side rather than vertically

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
                            value = "~ treatment"
                  ),
                  
                  input_switch("relevel_switch", "
                               Define reference levels"), #switch to trigger multiselect option
                  uiOutput("reveal_relevel"), # multiselect option and commit button
                  uiOutput("reveal_relevel_choice"), # commit triggers final input lines
                  
                  input_switch("pre_filter_switch", 
                               "Apply Pre-filtering",
                               value = TRUE),
                  uiOutput("reveal_pre_filter"),
                  
                  input_switch("shrink_switch", 
                               "Apply LFCshrinkage"), #TODO: MED think of how to apply this easily
                  
                  tags$hr(),
                  
                  actionButton("run_dds", 
                               "Run Analysis"),
                  uiOutput("analysis_error_message")
                )
              ),
              
              card(
                card_header("Analysis Export"),
                downloadButton("export_dds", 
                             "Download DESeq Datastructure")
              )
              ),
    
    # Overview of the data ----
    # accordion panels with a few different tools collapsed inside them would
    # be a better set up to not be overwhelming
    nav_panel("Basics",
              accordion(
                accordion_panel(
                  "PCA",
                  plotOutput("PCA")
                ),
                accordion_panel(
                  "Volcano Plot",
                  plotOutput("Volcano")
                )
                )
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
      
      nav_panel("Load previous analysis from .rds", 
                card(
                  card_header("Load previous analysis"),
                  
                  fileInput("prev_model",
                             "Select a valid .rds where you saved a previous DESeq output",
                             accept = ".rds"
                  ),
                  
                  actionButton("load_prev_model",
                               "Load Previous Model"
                  )
                     )
      ),
      
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
  
  # Allowing user re-leveling: observing input$relevel_switch
  # Actually there are just conditional panels, so I don't need this and shoul recode the other reveal switch
  output$reveal_relevel <- renderUI({
    meta<-metaData()
    if (input$relevel_switch && !is.null(meta)) {
      tagList(
        selectInput(
          "col_relevel",
          "Select Covariate(s) to relevel:",
          choices = names(meta)[-1],
          multiple = TRUE,
          selected = names(df)[1]
        ),
        actionButton(
          "go_relevel",
          "Select"
        )
      )
    }else{NULL}
  })
  
  output$reveal_relevel_choice <- renderUI({
    meta<-metaData()
    
    req(meta, input$col_relevel)
    
    
    if(input$go_relevel && length(input$col_relevel) > 0){
      
      cov_to_re<-input$col_relevel
      
      elements<-list()
      #need to make a number of elements equal to the number of variables in meta
      for(i in 1:length(cov_to_re)){
        
        name<-cov_to_re[i]
        
        current_choices <- unique(meta[[name]])
        
        elements[[i]]<-selectInput(
          paste0("relevel_",name),
          paste0("Select reference condition for ",name,":"),
          choices = current_choices,
          selected = current_choices[1]
        )
      }
      
      tagList(elements)
      
    }else{NULL}
  })
  
  # Allowing user pre-filtering preferences: observing input input$pre_filter_switch
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
    else{NULL}
  })
  
  
  # Finally the last conditional switch is to do with applying LFC shrinkage
  # which will be a bonus feature to deal with at the very end
  
  # Data Input ------------------------------------------------
  
  # Reactive read in of the selected user files
  rawData <- reactive({
    req(input$data)
    
    df <- tryCatch(
      {
        read.csv(input$data$datapath,
                 header = TRUE)
      },
      error = function(e){NULL}
    )
    return(df)
  })
  
  metaData <- reactive({
    req(input$meta)
    
    df <- tryCatch(
      {
        read.csv(input$meta$datapath,
                 header = TRUE,
                 stringsAsFactors = TRUE)
      },
      error = function(e){NULL}
    )
    return(df)
  })
  
  # Initial Analysis Run ------------------------------------------------
  
  # variables we need to hang on to for exploratory analysis
  dds_object <- reactiveVal(NULL)
  meta_out <- reactiveVal(NULL)
  genes_out <- reactiveVal(NULL)
  
  observeEvent(input$run_dds, {
    
    cat("Starting analysis...\n")
    
    counts <- rawData()
    meta<-metaData()
    if(length(input$formula)>0){f_usr <- input$formula}
    
    req(counts,meta,f_usr) # silently errors: could add notifications to rawData() or metData()
    
    cat("Data read in...\n")
    
    # pre-validation specific to formula inputs
    f<-gsub(" ","",f_usr)
    if(grepl("~",f)){f<-str_extract(f,"(?<=~).*")} 
    # if the user made a formula with a twiddle, just take everything after the twiddle
    
    f_fac<-unlist(str_split(f,"[+,*,:,|,-]"))
    # these are the pieces we need to match against the meta colnames()
    f_fac<-f_fac[f_fac!=1] 
    #-1 is a common indicator to exclude the intercept in formulas which we don't need to check against the metadata
    
    if(!(ncol(counts) >= 2 && colnames(counts)[1]=="gene.id" && 
         colnames(counts)[2]=="gene.name")){
      showNotification("Error: The first two columns must be vectors of characters named gene.id and gene.name respectively",
                       type = "error",
                       duration = NULL)
      shiny:::reactiveStop()
    }
    
    if(!(length(colnames(counts[,-c(1,2)])) > 0 && length(meta[,1]) > 0 && 
       all(sort(colnames(counts[,-c(1,2)])) == sort(as.character(meta[,1]))))){
    showNotification("Error: Sample names in count data columns and 
           metadata do not match exactly when sorted. Please check for discrepancies.",
                     type = "error",
                     duration = NULL)
    shiny:::reactiveStop()
    }
    
    if(f!="~."){
      if(!all(f_fac %in% colnames(meta))){
        showNotification("Error: The pieces of your formula need to exactly match the column names of the metadata provided",
                         type = "error",
                         duration = NULL)
        shiny:::reactiveStop()
      }
    }
    
    cat("Passed Validation!\n")
    
    # Post validation data description:
    # Validation should have us proceeding with the following datastructures
    # counts <-
    # Gene_id | External_gene_name | Sample1_name | ... | SampleN_name
    # <chr>   | <chr>              | <numeric>    | ... | <numeric>
    
    # meta <- 
    # Sample_names | Treatment1 | ... | TreatmentN
    # <factor>     | <factor>   | ... | <factor>
    
    # f <- valid pieces to a formula
    
    
    # below I would like to enclose the analysis in a withProgress
    # I believe I can avoid a tryCatch statement because the above validation
    # should ensure that the following should run
    
    genes_out(cbind(counts[,1],counts[,2]))
    
    row.names(counts)<-counts[,1]
    
    counts<-counts[,-c(1,2)]
    
    row.names(meta)<-meta[,1]
    meta<-meta[,-1,drop=FALSE]
    meta_out(meta)
    
    f<-paste0("~",f)
    f<-formula(f)
    
    if(input$relevel_switch){
      # here is where releveling the selected columns actually happens
      # IN: User desired columns to relevel, and reference levels
      # input$col_relevel and input$relevel_`input$col_relevel[i]`
      
      relevel_names<-input$col_relevel
      
      for(i in relevel_names){
        meta[[i]]<-relevel(meta[[i]],ref=input[[paste0("relevel_",i)]])
      }
      
      cat("The user requested releveling of the following columns: ", relevel_names,"\n")
      for(i in 1:length(relevel_names)){
        cat("The reference level of ", relevel_names[i], " is: ", levels(meta[[relevel_names]])[i],"\n")
        }
    }
    
    # An okay application of withProgress at least lets the user see some form of
    # progress in the browser
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a minute', value = 0, min = 0, max = 100, {
                   
                   
                   dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta, design=f)
                   
                   setProgress(value= 10,message="Made DESeq matrix...")
                   cat(file=stderr(),"Made our matrix...\n")
                   
                   if(input$pre_filter_switch){
                     dds <- dds[rowSums(counts(dds) >= input$smallestCount) >= input$smallestGroup, ]
                   }
                   
                   if(input$shrink_switch){
                     # DOES NOTHING FOR NOW
                     #lfcShrink(dds, coef="fix needed", type="apeglm")
                   }
                   setProgress(value= 30,message="Running DESeq() analysis...")
                   cat(file=stderr(),"Running differential expression GLM...\n")
                   dds <- DESeq(dds)
                   
                   
                   setProgress(value= 90,message="Finished DESeq computation")
                   
                   dds_object(dds)
                   cat("Initial analysis complete... try out the downloader\n")
                   #return(dds)
                 })
    
    showNotification("DESeq() has been run, try out the downloader to save your analysis", 
                     duration = 8, 
                     type = "message")
    
    
  }) # end analysisOutput
  

  # TODO: VERY LOW- consider the benefits of clearing analyses when changing data
  # 
  # # Clear analysis output if file changes after an analysis was run: observing key changed inputs to indicate a new analysis is being run
  # observeEvent(c(input$data,input$meta,input$formula), {
  #   dds_object(NULL)
  # }, ignoreInit = TRUE)

  # DESeq Dataset download handler -----

  # TO DO: LOW =The best way to solve the future issue, is to maybe have the app export the
  # entire R environment, then when you reload it it just resets the app to a
  # saved state
  
  # Output processed into an .rds to be passed to the download handler
  output$export_dds <- downloadHandler(
    filename = function() {
      paste0("DESeq_app_dds_", Sys.Date(), ".rds")
    },
    content = function(file) {
      dds_to_download <- dds_object()
      
      req(dds_to_download)
      
      saveRDS(dds_to_download, file = file)
    }
  )
  
  
  
  # Load Previous Model (ISSUE: UNCOMPATIBLE WITH DOWNSTREAM ANALYSIS- see DESeq Dataset download handler todo) ----
  # Looks for trigger from "load_prev_model"
  # requires a model has been uploaded, readRDS(), validate then
  # store in dds_object(dds), to overwrite any analysis currently being stored
  
  # I didn't think of this but can we even load a previous model,
  # a lot of this is contingent on having the other pieces of the data available
  # to snag from, but the dds should be built of the same pieces
  # actually we can rebuild out of dds@colData, every factor vector in this dataframe
  # is one of the covariates and dds@design delivers us the formula
  # the only thing missing would be the gene names associated with the metadata
  # which we could map using some BioConductor function
  
  # or is it just easier to require them to reupload the metadata
  
  # No idea if this will work for all downstream applications
  
  observeEvent(input$load_prev_model,{
    req(input$prev_model)
    
    tmp <- tryCatch({
      readRDS(input$prev_model$datapath)
    }, error = function(e){NULL}
      )
    
    if(is.null(tmp) || class(tmp)!="DESeqDataSet"){
      showNotification(
        "Error: The uploaded file is not a valid 'DESeqDataSet' object.",
        type = "error"
      )
      return()
    }
    
    dds_object(tmp)
      
    cat("Analysis uploaded successfully with results loaded in as\n", 
        resultsNames(dds_object()),
        "\nIf this looks incorrect please validate the model was downloaded correctly\n")
      
  })
  
  # Post Analysis Functionality ----
  
  # Expects to have a DESeq Dataset object to call results on
  # as well as the metadata associated with the DESeq Object
  
  # All other features actually using the post-processing data here
  
  # PCA and initial volcano plots.
  # This working properly depends on setting the reference level set properly
  
  # PCA is easy just run the rda_all with a switch that toggle the BLIND parameter

  
  output$PCA <- renderPlot({
    req(dds_object())
    
    dds<-dds_object()
    
    cat(file = stderr(),"PCA triggered...\n")
    
    # for PCAs we need normalized counts and not LFC
    rld_all <- rlog(dds, blind = FALSE)
    rld_all_df <- as.data.frame(assay(rld_all))
    
    cat(file = stderr(),"PCA needed calculation done...\n")
    
    
    # TODO: NEXT- figure out how to make the PCA customizable
    
    PCA_plot_all <- plotPCA(rld_all, intgroup = "treatment", returnData=TRUE,ntop=1000)
    percentVar <- round(100 * attr(PCA_plot_all, "percentVar"))
    
    PCA_plot_all$treatment<-sub("\\."," ",PCA_plot_all$treatment)
    
    # need to make a dropdown that allows the user to specify intgroup (groups of interest e.g. one of the columns),color, and label
    # as well as a color palette for the PCA
    
    PCA_group <- ggplot(PCA_plot_all, aes(x=PC1, y=PC2, color=treatment,label=treatment)) +
      geom_point(size=5) +
      scale_x_continuous(name=paste0("PC1: ",percentVar[1],"% variance"), limits=c(-20,20))+
      scale_y_continuous(name=paste0("PC2: ",percentVar[2],"% variance"), limits=c(-20,20))+
      coord_fixed()+
      ggforce::geom_mark_ellipse(aes(label=NULL,color = treatment))+
      geom_text_repel(size=7,vjust = "inward", hjust = "inward",
                      show.legend = FALSE, point.padding=10)+
      ggtitle ("PCA")+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(axis.text=element_text(size=20), 
            axis.title=element_text(size=25), 
            title=element_text(size=25,face="bold"), 
            legend.text=element_text(size=20))
    PCA_group
  })
  
  # some option to save the PCA
  
  
  
  # Volcano plots are a little harder, depending on the formula there could be
  # multiple different ways that the results could be stored or displayed in a DESeqDataSet object
  
  # - the simplest is one intercept with the treatment listed with the different comparison to the reference level
  # - a more detailed GLM formula where displaying the differences becomes a matter of, picking the
  # the relevant covariate and the contrasts within
  # - These options are made more complicated by the different calls you can make
  # to DESeq to display the results # TODO: investigate the vignette to determine options
  # the two options are name and contrast
  # 
  # One exception to the equivalence of these two commands (name and contrast), is that, using contrast 
  # will additionally set to 0 the estimated LFC in a comparison of two groups, where 
  # all of the counts in the two groups are equal to 0 (while other groups have positive counts). 
  # As this may be a desired feature to have the LFC in these cases set to 0, one 
  # can use contrast to build these results tables. More information about extracting 
  # specific coefficients from a fitted DESeqDataSet object can be found in the help page 
  
  
  
  
  
}

shinyApp(ui = ui, server = server)









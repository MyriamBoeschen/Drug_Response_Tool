library(shiny)
library(stringr)
library(epitools)
library(ggplot2)
library(gridExtra)
library(grid)
#library(RColorBrewer)
library(writexl)
library(shinycssloaders)
#library(dplyr)
library(shinyjs) 
library(DT)

ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
        .shiny-output-error-validation {
          color: red;
          font-weight: bold;
          font-size:20px;
        }
      "))
  ),
  
  #titlePanel(h3("GDSC Drug Response Tool")),
  fluidRow(
    useShinyjs(),
    
    #column(4,
     #      tags$img(src="GDSC_logo.png", title="GDSC_logo", width="315", height="100")
    #),
    column(3,
           h2("Drug Response Tool")
    ),
    column(4,
           h4("")
    ),
    column(1,
           h4(a("Home", href='https://www.hornlab.org', target='_blank'))
    ),
    column(1,
           tags$a(href="https://www.uni-due.de/phenotime/", tags$img(src="logo_phenotime.png", title="Phenotime", width="120", height="40"))
    ),
    column(1,
           tags$a(href="https://biochemie.medizin.uni-leipzig.de/bch_cms/index.php/en/", tags$img(src="logo_RSI.png", title="RSI",
                                                                                                  width="130", height="50") )
    ),
    column(1,
           tags$a(href="https://biochemie.medizin.uni-leipzig.de/bch_cms/index.php/en/", tags$img(src="logo_UL.png", title="UL",
                                                                                                  width="130", height="50"))
    ),
    
    column(12,
           h4("Test the effect of drugs on cancer cell lines that are genetically changed compared to wild type cancer cell lines. The grouping of cell lines for comparison will be based on genetic data at www.cBioportal.org and prepared as a table manually or using this tool. The drug data and genetic data provided for analysis here were ascertained in Yang et al. 2013, PMID:23180760; Iorio et al. 2016, PMID:27397505 and Garnett et al. 2012, PMID:27397505."),
           #h5("The drug data and genetic data provided here were ascertained in Yang et al. 2013, PMID:23180760; Iorio et al. 2016, PMID:27397505 and Garnett et al. 2012, PMID:27397505."),
           h5(a("Step 1 - ", strong("Obtain a file"),  "'alterations_across_samples.tsv' for cell lines", strong("from www.cBioportal.org"), href='Manual.pdf', target='_blank')),
           h5("Step 2 - Use this tool to create and analyze a groups.txt file from the .tsv table with up to three genes queried. By default we define supposedly activating changes (DNA amplifications or mRNA high) and deactivating changes (DNA homozygous deletions or mRNA low). You can also include point mutations in the altered group, leave them out or analyze point mutations only."),
           h5("Alternatively, upload custom groups.txt file", a(" - see example", href='groups.txt', target='_blank')),
           #h5("NOTE: Make sure to upload a .tsv file from cbioportal.org", strong("OR"), "a selfmade groups.txt file.")
    )
  ),
  hr(), ### line ###
  sidebarLayout(
    sidebarPanel(
           splitLayout(cellWidths = c("85%", "15%"),
                       fileInput("alterations","Upload .tsv file retrieved from cbioportal.org"),
                       actionButton('reset', 'reset', style="color: black;
                       background-color:  white; 
                       margin-top: 25px;
                       text-align:center;
                       text-indent: -2px")
           ),
           splitLayout(cellWidths = c("85%", "15%"),
                       fileInput("groupstxt","Upload custom groups.txt file"),
                       actionButton('RESET', 'reset', style="color: black;
                       background-color:  white;
                       margin-top: 25px;
                       text-align:center;
                       text-indent: -2px")
           ),
           selectInput("tissue",label = "Analyze specific tissue origins only?", choices = list("across entities"="EN",
                                                                                                "lung_NSCLC"="lung_NSCLC",
                                                                                                "urogenital_system"= "urogenital_system",
                                                                                                "digestive_system"="digestive_system",
                                                                                                "nervous_system"="nervous_system",
                                                                                                "skin"="skin",
                                                                                                "leukemia"="leukemia",
                                                                                                "kidney"="kidney",
                                                                                                "thyroid"="thyroid",
                                                                                                "soft_tissue"="soft_tissue",
                                                                                                "aero_dig_tract"="aero_dig_tract",
                                                                                                "lymphoma"="lymphoma",
                                                                                                "myeloma"="myeloma",
                                                                                                "pancreas"="pancreas",
                                                                                                "breast"="breast",
                                                                                                "neuroblastoma"="neuroblastoma",
                                                                                                "large_intestine"="large_intestine",
                                                                                                "bone"="bone",
                                                                                                "lung_SCLC"="lung_SCLC",
                                                                                                "lung"="lung"), multiple = TRUE, selected = list("across entities"="EN")),
           hr(),
           uiOutput("test"),
           selectInput("amountgenes","Number of genes analyzed (1 to 3 possible)",choices = list("one"="One", 
                                                                                                 "two"="Two",
                                                                                                 "three"="Three")),
           selectInput("notprofiled",label = "Exclude cell lines not profiled in:", choices = list("exclude no data"="EN",
                                                                                                   "copy number variations"="..AMP",
                                                                                                   "mRNA alterations"= "..EXP..2",
                                                                                                   "point mutations"="..MUT",
                                                                                                   "fusions"="..FUSION"), multiple = TRUE, selected = list("copy number variations"="..AMP", "mRNA alterations"= "..EXP..2")),
           hr(),
           selectizeInput("oncgenes", "Select gene name", choices=NULL, multiple = FALSE ),
           #textInput("oncgenes","Enter gene name"),
           tags$div(id="txt",h5(strong("Choose alteration type to test against wildtype"), style="margin-bottom: 0px; margin-top: 15px")),
           fluidRow(
             column(6,
                    h5(radioButtons("mutationkind",
                                    label = "",
                                    choices = list("presumably activating changes"="am",
                                                   "presumably deactivating changes"="dm"),
                                    selected = list("presumably activating changes"="am")), 
                       style="margin-top: 0px")),
                    column(6,
                             conditionalPanel(
                               condition="input.mutationkind=='am'",
                               h5(checkboxGroupInput("pointmut",
                                                     label = "",
                                                     choices = list("CNV up"="up",
                                                                    "CNV down"="down",
                                                                    "mRNA high"="high",
                                                                    "mRNA low"="low",
                                                                    "missense"="aaa",
                                                                    "frameshift"="fs",
                                                                    "nonsense"= "tm",
                                                                    "splice"="sv",
                                                                    "fusions"="fn"), 
                                                     selected = list("CNV up"="up", "mRNA high"="high")),
                                  style="margin-top: 0px")),
                             
                             conditionalPanel(
                               condition="input.mutationkind=='dm'",
                               h5(checkboxGroupInput("pointmutdm",
                                                     label = "",
                                                     choices = list("CNV up"="up",
                                                                    "CNV down"="down",
                                                                    "mRNA high"="high",
                                                                    "mRNA low"="low",
                                                                    "missense"="aaa",
                                                                    "frameshift"="fs",
                                                                    "nonsense"= "tm",
                                                                    "splice"="sv",
                                                                    "fusions"="fn"), selected = list("CNV down"="down",
                                                                                                    "mRNA low"="low",
                                                                                                    "frameshift"="fs",
                                                                                                    "nonsense"= "tm")),
                                  style="margin-top: 0px"))
                             )#column
           ), #fluidrow
           conditionalPanel(
             condition="input.amountgenes=='Two'",
             selectizeInput("secondgene", "Select second gene", choices=NULL, multiple = FALSE ),
             #textInput("secondgene","Enter second gene name"),
             tags$div(id="txt",h5(strong("Choose alteration type to test against wildtype"), style="margin-bottom: 0px; margin-top: 15px")),
             fluidRow(
               column(6,
                      h5(radioButtons("mutationkind_second",
                          label = "",
                          choices = list("presumably activating changes"="amm",
                                         "presumably deactivating changes"="dmm"),
                          selected = list("presumably activating changes"="amm")),
                         style="margin-top: 0px")),
               column(6,
                    conditionalPanel(
                    condition="input.mutationkind_second=='amm'",
                    h5(checkboxGroupInput("secondpointmut",
                                label = "",
                                choices = list("CNV up"="up",
                                               "CNV down"="down",
                                               "mRNA high"="high",
                                               "mRNA low"="low",
                                               "missense"="aaa",
                                               "frameshift"="fs",
                                               "nonsense"= "tm",
                                               "splice"="sv",
                                               "fusions"="fn"),
                                selected = list("CNV up"="up",
                                                "mRNA high"="high")),
                       style="margin-top: 0px")),
                    conditionalPanel(
                      condition="input.mutationkind_second=='dmm'",
                      h5(checkboxGroupInput("secondpointmutdm",
                                            label = "",
                                            choices = list("CNV up"="up",
                                                           "CNV down"="down",
                                                           "mRNA high"="high",
                                                           "mRNA low"="low",
                                                           "missense"="aaa",
                                                           "frameshift"="fs",
                                                           "nonsense"= "tm",
                                                           "splice"="sv",
                                                           "fusions"="fn"),
                                            selected = list("CNV down"="down",
                                                            "mRNA low"="low",
                                                            "frameshift"="fs",
                                                            "nonsense"= "tm")),
                         style="margin-top: 0px"))
                    ) #column
               ), #fluidrow
             radioButtons("mutgroup",
                          label = "Define cell line as mutated if",
                          choices = list("at least one of all specified genes is mutated"="onlyone",
                                         "all specified genes are co-mutated only"="comut"),
                          selected = "onlyone")
           ), #conditionalpanel if 2 genes selected
           conditionalPanel(
             condition="input.amountgenes=='Three'",
             selectizeInput("second3gene", "Select second gene", choices=NULL, multiple = FALSE ),
             #textInput("second3gene","Enter second gene name"),
             tags$div(id="txt",h5(strong("Choose alteration type to test against wildtype"), style="margin-bottom: 0px; margin-top: 15px")),
             fluidRow(
               column(6,
                      h5(radioButtons("mutationkind_second3",
                                      label = "",
                                      choices = list("presumably activating changes"="amm",
                                                     "presumably deactivating changes"="dmm"),
                                      selected = list("presumably activating changes"="amm")),
                         style="margin-top: 0px")),
               column(6,
                      conditionalPanel(
                        condition="input.mutationkind_second3=='amm'",
                        h5(checkboxGroupInput("second3pointmut",
                                              label = "",
                                              choices = list("CNV up"="up",
                                                             "CNV down"="down",
                                                             "mRNA high"="high",
                                                             "mRNA low"="low",
                                                             "missense"="aaa",
                                                             "frameshift"="fs",
                                                             "nonsense"= "tm",
                                                             "splice"="sv",
                                                             "fusions"="fn"),
                                              selected = list("CNV up"="up",
                                                              "mRNA high"="high")),
                           style="margin-top: 0px")),
                      conditionalPanel(
                        condition="input.mutationkind_second3=='dmm'",
                        h5(checkboxGroupInput("second3pointmutdm",
                                              label = "",
                                              choices = list("CNV up"="up",
                                                             "CNV down"="down",
                                                             "mRNA high"="high",
                                                             "mRNA low"="low",
                                                             "missense"="aaa",
                                                             "frameshift"="fs",
                                                             "nonsense"= "tm",
                                                             "splice"="sv",
                                                             "fusions"="fn"),
                                              selected = list("CNV down"="down",
                                                              "mRNA low"="low",
                                                              "frameshift"="fs",
                                                              "nonsense"= "tm")),
                           style="margin-top: 0px"))
               ) #column
             ), #fluidrow
             selectizeInput("third3gene", "Select third gene", choices=NULL, multiple = FALSE ),
             #textInput("third3gene","Enter third gene name"),
             tags$div(id="txt",h5(strong("Choose alteration type to test against wildtype"), style="margin-bottom: 0px; margin-top: 15px")),
             fluidRow(
               column(6,
                      h5(radioButtons("mutationkind_third3",
                                      label = "",
                                      choices = list("presumably activating changes"="ammm",
                                                     "presumably deactivating changes"="dmmm"),
                                      selected = list("presumably activating changes"="ammm")),
                         style="margin-top: 0px")),
               column(6,
                      conditionalPanel(
                        condition="input.mutationkind_third3=='ammm'",
                        h5(checkboxGroupInput("third3pointmut",
                                              label = "",
                                              choices = list("CNV up"="up",
                                                             "CNV down"="down",
                                                             "mRNA high"="high",
                                                             "mRNA low"="low",
                                                             "missense"="aaa",
                                                             "frameshift"="fs",
                                                             "nonsense"= "tm",
                                                             "splice"="sv",
                                                             "fusions"="fn"),
                                              selected = list("CNV up"="up",
                                                              "mRNA high"="high")),
                           style="margin-top: 0px")),
                      conditionalPanel(
                        condition="input.mutationkind_third3=='dmmm'",
                        h5(checkboxGroupInput("third3pointmutdm",
                                              label = "",
                                              choices = list("CNV up"="up",
                                                             "CNV down"="down",
                                                             "mRNA high"="high",
                                                             "mRNA low"="low",
                                                             "missense"="aaa",
                                                             "frameshift"="fs",
                                                             "nonsense"= "tm",
                                                             "splice"="sv",
                                                             "fusions"="fn"),
                                              selected = list("CNV down"="down",
                                                              "mRNA low"="low",
                                                              "frameshift"="fs",
                                                              "nonsense"= "tm")),
                           style="margin-top: 0px"))
               ) #column
             ), #fluidrow
             radioButtons("mutgroup3",
                          label = "Define cell line as mutated if",
                          choices = list("at least one of all specified genes is mutated"="onlyone",
                                         "all specified genes are co-mutated only"="comut"),
                          selected = "onlyone")
           ),
           actionButton(inputId = "submit",label = "submit", style="color: black; 
                       background-color:   #7fb3d5; 
                       position: relative; 
                       #left: 3%;
                       height: 35px;
                       width: 100px;
                       text-align:center;
                       text-indent: -2px;
                       border-radius: 6px;
                       border-width: 2px"),
           hr(),
           downloadButton("groupsfiles", "Download groups.txt file", style="color: black; 
                       background-color:  white; 
                       position: relative; 
                       #left: 3%;
                       height: 35px;
                       width: 200px;
                       text-align:center;
                       text-indent: -2px;
                       border-radius: 5px;
                       border-width: 1px"),
           downloadButton("downloadData", "Download tables", style="color: black; 
                       background-color:  white; 
                       position: relative; 
                       #left: 3%;
                       height: 35px;
                       width: 200px;
                       text-align:center;
                       text-indent: -2px;
                       border-radius: 5px;
                       border-width: 1px"),
           downloadButton("downloadData_plot", "Download barplot",style="color: black; 
                       background-color:  white; 
                       position: relative; 
                       #left: 3%;
                       height: 35px;
                       width: 200px;
                       text-align:center;
                       text-indent: -2px;
                       border-radius: 5px;
                       border-width: 1px"),
           downloadButton("downloadData_forest", "Download forestplot",style="color: black; 
                       background-color:  white; 
                       position: relative; 
                       #left: 3%;
                       height: 35px;
                       width: 200px;
                       text-align:center;
                       text-indent: -2px;
                       border-radius: 5px;
                       border-width: 1px"),
           downloadButton("downloadData_boxplot", "Download boxplots",style="color: black; 
                       background-color:  white; 
                       position: relative; 
                       #left: 3%;
                       height: 35px;
                       width: 200px;
                       text-align:center;
                       text-indent: -2px;
                       border-radius: 5px;
                       border-width: 1px"),
           downloadButton("downloadData_wilcoxtable", "Download boxplot data",style="color: black; 
                       background-color:  white; 
                       position: relative; 
                       #left: 3%;
                       height: 35px;
                       width: 200px;
                       text-align:center;
                       text-indent: -2px;
                       border-radius: 5px;
                       border-width: 1px")
    
    ), #sidebarpanel
    mainPanel(
      tabsetPanel(type="tabs", id = "tabs",
                  tabPanel("Intro",
                           fluidRow(
                             column(12,
                                    #htmlOutput("title_example"), 
                                    uiOutput("pdfview")),
                                    #htmlOutput("title_example_forest"), 
                                    #uiOutput("pdfview_forest")),
                           )
                  ),
                  tabPanel("Table",
                           fluidRow(
                             column(12,
                                    htmlOutput("title_table"), 
                                    withSpinner(DT::dataTableOutput("Table_druggroups"),proxy.height = "400px")),
                             #DT::dataTableOutput("Table_druggroups")),
                           )
                  ),
                  tabPanel("Barplot",
                           fluidRow(
                             column(12,
                                    htmlOutput("title_plot"), 
                                    withSpinner(plotOutput("plotop",height = "2000px", width = "800px"),proxy.height = "400px")),
                           )
                  ),
                  tabPanel("Forestplot",
                           fluidRow(
                             column(12,
                                    htmlOutput("title_forestplot"), 
                                    withSpinner(plotOutput("Forest",height = "auto"),proxy.height = "400px")),
                           )
                  ),
                  tabPanel("Boxplot",
                           fluidRow(
                             column(12,
                                    htmlOutput("title_boxplot"), 
                                    withSpinner(plotOutput("BoxPlot",height = "auto", width = "1000px"),proxy.height = "400px"))
                           )
                  )
                  ) #tabsetpanel
      ) #mainpanel
    ) #sidebarlayout
  ) #fluidpage

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

server <- function(input, output,session) {
  
  output$pdfview <- renderUI({
    tags$iframe(style="height: 1000px; width:100%; border:0", src="Intro.png")
  })
  
  output$pdfview_forest <- renderUI({
    tags$iframe(style="height:1067px; width:100%; border:0", src="BRAFV600_forestplot_800px.png")
  })
  
  values <- reactiveValues(
    upload_state = NULL
  )
  
  values2 <- reactiveValues(
    upload_state = NULL
  )
  
  observeEvent(input$alterations, {
    values$upload_state <- 'uploaded'
  })
  
  observeEvent(input$groupstxt, {
    values2$upload_state <- 'uploaded2'
  })
  
  observeEvent(input$reset, {
    values$upload_state <- 'reset'
  })
  
  observeEvent(input$RESET, {
    values2$upload_state <- 'reset2'
  })
  
  file_input <- reactive({
    if (is.null(values$upload_state)) {
      return(NULL)
    } else if (values$upload_state == 'uploaded') {
      return(input$alterations)
    } else if (values$upload_state == 'reset') {
      return(NULL)
    }
  })
  
  file_input_groups <- reactive({
    if (is.null(values2$upload_state)) {
      return(NULL)
    } else if (values2$upload_state == 'uploaded2') {
      return(input$groupstxt)
    } else if (values2$upload_state == 'reset2') {
      return(NULL)
    }
  })
  
  
  observeEvent(input$groupstxt, {
    removeUI(selector = sprintf('.shiny-input-container:has(#%s)','oncgenes'),multiple = TRUE, immediate = TRUE)
    removeUI(selector = sprintf('.shiny-input-container:has(#%s)','amountgenes'),multiple = TRUE, immediate = TRUE)
    removeUI(selector = sprintf('.shiny-input-container:has(#%s)','notprofiled'),multiple = TRUE, immediate = TRUE)
    removeUI(selector = "#mutationkind", multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#pointmut', multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#pointmutdm', multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#mutationkind_second', multiple = TRUE, immediate = TRUE)
    removeUI(selector = sprintf('.shiny-input-container:has(#%s)','secondgene'),multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#secondpointmut', multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#secondpointmutdm', multiple = TRUE, immediate = TRUE)
    removeUI(selector = sprintf('.shiny-input-container:has(#%s)','second3gene'),multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#mutationkind_second3', multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#second3pointmut', multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#second3pointmutdm', multiple = TRUE, immediate = TRUE)
    removeUI(selector = sprintf('.shiny-input-container:has(#%s)','third3gene'),multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#mutationkind_third3', multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#third3pointmut', multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#third3pointmutdm', multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#mutgroup', multiple = TRUE, immediate = TRUE)
    removeUI(selector = '#mutgroup3', multiple = TRUE, immediate = TRUE)
    removeUI(selector = 'div#txt', multiple = TRUE, immediate = TRUE)
  })
  
  output$test <- renderUI({
    
    if (is.null(values2$upload_state)) {
      return(NULL)
    } else if (values2$upload_state == 'reset2') {
      tagList(
        selectInput("amountgenes","Number of genes analyzed (1 to 3 possible)",choices = list("one"="One", 
                                                                                              "two"="Two",
                                                                                              "three"="Three")),
        selectInput("notprofiled",label = "Exclude cell lines not profiled in:", choices = list("exclude no data"="EN",
                                                                                                "copy number variations"="..AMP",
                                                                                                "mRNA alterations"= "..EXP..2",
                                                                                                "point mutations"="..MUT",
                                                                                                "fusions"="..FUSION"), multiple = TRUE, selected = list("copy number variations"="..AMP", "mRNA alterations"= "..EXP..2")),
        hr(),
        selectizeInput("oncgenes", "Select gene name", choices=NULL, multiple = FALSE ),
        #textInput("oncgenes","Enter gene name"),
        tags$div(id="txt",h5(strong("Choose alteration type to test against wildtype"), style="margin-bottom: 0px; margin-top: 15px")),
        fluidRow(
          column(6,
                 h5(radioButtons("mutationkind",
                                 label = "",
                                 choices = list("presumably activating changes"="am",
                                                "presumably deactivating changes"="dm"),
                                 selected = list("presumably activating changes"="am")), 
                    style="margin-top: 0px")),
          column(6,
                 conditionalPanel(
                   condition="input.mutationkind=='am'",
                   h5(checkboxGroupInput("pointmut",
                                         label = "",
                                         choices = list("CNV up"="up",
                                                        "CNV down"="down",
                                                        "mRNA high"="high",
                                                        "mRNA low"="low",
                                                        "missense"="aaa",
                                                        "frameshift"="fs",
                                                        "nonsense"= "tm",
                                                        "splice"="sv",
                                                        "fusions"="fn"), 
                                         selected = list("CNV up"="up", "mRNA high"="high")),
                      style="margin-top: 0px")),
                 
                 conditionalPanel(
                   condition="input.mutationkind=='dm'",
                   h5(checkboxGroupInput("pointmutdm",
                                         label = "",
                                         choices = list("CNV up"="up",
                                                        "CNV down"="down",
                                                        "mRNA high"="high",
                                                        "mRNA low"="low",
                                                        "missense"="aaa",
                                                        "frameshift"="fs",
                                                        "nonsense"= "tm",
                                                        "splice"="sv",
                                                        "fusions"="fn"), selected = list("CNV down"="down",
                                                                                         "mRNA low"="low",
                                                                                         "frameshift"="fs",
                                                                                         "nonsense"= "tm")),
                      style="margin-top: 0px"))
          )#column
        ), #fluidrow
        conditionalPanel(
          condition="input.amountgenes=='Two'",
          selectizeInput("secondgene", "Select second gene", choices=NULL, multiple = FALSE ),
          #textInput("secondgene","Enter second gene name"),
          tags$div(id="txt",h5(strong("Choose alteration type to test against wildtype"), style="margin-bottom: 0px; margin-top: 15px")),
          fluidRow(
            column(6,
                   h5(radioButtons("mutationkind_second",
                                   label = "",
                                   choices = list("presumably activating changes"="amm",
                                                  "presumably deactivating changes"="dmm"),
                                   selected = list("presumably activating changes"="amm")),
                      style="margin-top: 0px")),
            column(6,
                   conditionalPanel(
                     condition="input.mutationkind_second=='amm'",
                     h5(checkboxGroupInput("secondpointmut",
                                           label = "",
                                           choices = list("CNV up"="up",
                                                          "CNV down"="down",
                                                          "mRNA high"="high",
                                                          "mRNA low"="low",
                                                          "missense"="aaa",
                                                          "frameshift"="fs",
                                                          "nonsense"= "tm",
                                                          "splice"="sv",
                                                          "fusions"="fn"),
                                           selected = list("CNV up"="up",
                                                           "mRNA high"="high")),
                        style="margin-top: 0px")),
                   conditionalPanel(
                     condition="input.mutationkind_second=='dmm'",
                     h5(checkboxGroupInput("secondpointmutdm",
                                           label = "",
                                           choices = list("CNV up"="up",
                                                          "CNV down"="down",
                                                          "mRNA high"="high",
                                                          "mRNA low"="low",
                                                          "missense"="aaa",
                                                          "frameshift"="fs",
                                                          "nonsense"= "tm",
                                                          "splice"="sv",
                                                          "fusions"="fn"),
                                           selected = list("CNV down"="down",
                                                           "mRNA low"="low",
                                                           "frameshift"="fs",
                                                           "nonsense"= "tm")),
                        style="margin-top: 0px"))
            ) #column
          ), #fluidrow
          radioButtons("mutgroup",
                       label = "Define cell line as mutated if",
                       choices = list("at least one of all specified genes is mutated"="onlyone",
                                      "all specified genes are co-mutated only"="comut"),
                       selected = "onlyone")
        ), #conditionalpanel if 2 genes selected
        conditionalPanel(
          condition="input.amountgenes=='Three'",
          selectizeInput("second3gene", "Select second gene", choices=NULL, multiple = FALSE ),
          #textInput("second3gene","Enter second gene name"),
          tags$div(id="txt",h5(strong("Choose alteration type to test against wildtype"), style="margin-bottom: 0px; margin-top: 15px")),
          fluidRow(
            column(6,
                   h5(radioButtons("mutationkind_second3",
                                   label = "",
                                   choices = list("presumably activating changes"="amm",
                                                  "presumably deactivating changes"="dmm"),
                                   selected = list("presumably activating changes"="amm")),
                      style="margin-top: 0px")),
            column(6,
                   conditionalPanel(
                     condition="input.mutationkind_second3=='amm'",
                     h5(checkboxGroupInput("second3pointmut",
                                           label = "",
                                           choices = list("CNV up"="up",
                                                          "CNV down"="down",
                                                          "mRNA high"="high",
                                                          "mRNA low"="low",
                                                          "missense"="aaa",
                                                          "frameshift"="fs",
                                                          "nonsense"= "tm",
                                                          "splice"="sv",
                                                          "fusions"="fn"),
                                           selected = list("CNV up"="up",
                                                           "mRNA high"="high")),
                        style="margin-top: 0px")),
                   conditionalPanel(
                     condition="input.mutationkind_second3=='dmm'",
                     h5(checkboxGroupInput("second3pointmutdm",
                                           label = "",
                                           choices = list("CNV up"="up",
                                                          "CNV down"="down",
                                                          "mRNA high"="high",
                                                          "mRNA low"="low",
                                                          "missense"="aaa",
                                                          "frameshift"="fs",
                                                          "nonsense"= "tm",
                                                          "splice"="sv",
                                                          "fusions"="fn"),
                                           selected = list("CNV down"="down",
                                                           "mRNA low"="low",
                                                           "frameshift"="fs",
                                                           "nonsense"= "tm")),
                        style="margin-top: 0px"))
            ) #column
          ), #fluidrow
          selectizeInput("third3gene", "Select third gene", choices=NULL, multiple = FALSE ),
          #textInput("third3gene","Enter third gene name"),
          tags$div(id="txt",h5(strong("Choose alteration type to test against wildtype"), style="margin-bottom: 0px; margin-top: 15px")),
          fluidRow(
            column(6,
                   h5(radioButtons("mutationkind_third3",
                                   label = "",
                                   choices = list("presumably activating changes"="ammm",
                                                  "presumably deactivating changes"="dmmm"),
                                   selected = list("presumably activating changes"="ammm")),
                      style="margin-top: 0px")),
            column(6,
                   conditionalPanel(
                     condition="input.mutationkind_third3=='ammm'",
                     h5(checkboxGroupInput("third3pointmut",
                                           label = "",
                                           choices = list("CNV up"="up",
                                                          "CNV down"="down",
                                                          "mRNA high"="high",
                                                          "mRNA low"="low",
                                                          "missense"="aaa",
                                                          "frameshift"="fs",
                                                          "nonsense"= "tm",
                                                          "splice"="sv",
                                                          "fusions"="fn"),
                                           selected = list("CNV up"="up",
                                                           "mRNA high"="high")),
                        style="margin-top: 0px")),
                   conditionalPanel(
                     condition="input.mutationkind_third3=='dmmm'",
                     h5(checkboxGroupInput("third3pointmutdm",
                                           label = "",
                                           choices = list("CNV up"="up",
                                                          "CNV down"="down",
                                                          "mRNA high"="high",
                                                          "mRNA low"="low",
                                                          "missense"="aaa",
                                                          "frameshift"="fs",
                                                          "nonsense"= "tm",
                                                          "splice"="sv",
                                                          "fusions"="fn"),
                                           selected = list("CNV down"="down",
                                                           "mRNA low"="low",
                                                           "frameshift"="fs",
                                                           "nonsense"= "tm")),
                        style="margin-top: 0px"))
            ) #column
          ), #fluidrow
          radioButtons("mutgroup3",
                       label = "Define cell line as mutated if",
                       choices = list("at least one of all specified genes is mutated"="onlyone",
                                      "all specified genes are co-mutated only"="comut"),
                       selected = "onlyone")
        )
    )}
  })
  
  output$title_example <-renderUI({        ## output text intro on Example
    HTML(paste("<br/><p style='font-size:14pt;'>Example output shows BRAF V600 mutation analysis across entities</p>Genetically altered samples (mut) are compared to wild type samples (wt) with respect to their resistance (red) or sensitivity (white) to each group of drugs. Numbers of cell lines in each category are indicated in the bars: altered and resistant, altered and sensitive, wild type and resistant, wild type and sensitive. Odds ratios for resistance in the altered cell lines are computed and shown above each barplot. An OR of >1 indicates increased resistance of altered cell lines to drugs in this drug pathway (Fisher tests). An OR<1 indicates sensitivity. Here, the output for BRAF V600 mutations are shown<br/><br/>"))
  })
  output$title_table <-renderUI({        ## output text intro on Table
    HTML(paste("<br/>Genetically altered samples (mut) are compared to wild type samples (wt) with respect to their resistance or sensitivity to each group of drugs. Odds ratios for resistance in the altered cell lines are computed. An OR of >1 indicates increased resistance of altered cell lines to drugs in this drug pathway and an OR<1 indicates sensitivity. Significance is calculated using fisher's tests and chi-squared tests. Low p-values are indicative of pathways that the altered cell lines may be resistant/sensitive to. Full data, including group sizes, are provided in tables. After Bonferroni correction for testing 22 pathways, p-values <0.002 should be considered significant. <br/><br/>"))
  })
  output$title_plot <-renderUI({        ## output text intro on BarPlot
    HTML(paste("<br/>Genetically altered samples (mut) are compared to wild type samples (wt) with respect to their resistance (red) or sensitivity (white) to each group of drugs. Numbers of cell lines in each category are indicated in the bars: altered and resistant, altered and sensitive, wild type and resistant, wild type and sensitive. Odds ratios for resistance in the altered cell lines are computed and shown above each barplot. An OR of >1 indicates increased resistance of altered cell lines to drugs in this drug pathway (Fisher tests). An OR<1 indicates sensitivity. After Bonferroni correction for testing 22 pathways, p-values <0.002 should be considered significant. <br/><br/>"))
  })
  output$title_forestplot <-renderUI({        ## output text intro on ForestPlot
    HTML(paste("<br/>The forestplot indicates individual drugs that are associated with resistance and sensitivity in altered cell lines compared to the wild type cell lines (Fisher's p<0.05; OR>1 for resistance, OR<1 for sensitivity). If odds ratios can not be calculated due to a group size of 0, these significant compounds are shown at OR=0 (sensitive) or OR=7 (resistant) without a CI. Full data are given in tables provided for download in the lower left panel. After Bonferroni correction for testing 250 compounds, p-values <0.0002 should be considered significant. <br/><br/>"))
  })
  output$title_example_forest <-renderUI({        ## output text intro on Example ForestPlot
    HTML(paste("<br/>The forestplot indicates individual drugs that are associated with resistance and sensitivity in altered cell lines compared to the wild type cell lines (Fisher's p<0.05; OR>1 for resistance, OR<1 for sensitivity). If odds ratios can not be calculated due to a group size of 0, these significant compounds are shown at OR=0 (sensitive) or OR=7 (resistant) without a CI. Full data are given in tables provided for download in the lower left panel. After Bonferroni correction for testing 250 compounds, p-values <0.0002 should be considered significant. <br/><br/>")) 
    })
  
  output$title_boxplot <-renderUI({        ## output text intro on BoxPlot
    HTML(paste("<br/>IC50 values are compared using two-sided Wilcoxon tests. Full boxplot data are provided for download in the lower left panel. After Bonferroni correction for testing 250 compounds, p-values <0.0002 should be considered significant. <br/><br/>")) 
  })
  
  observeEvent(input$submit, {
    hideTab(inputId = "tabs", target = "Example Output")
  })
  
  data_input <- reactive({
    
    if (is.null(file_input())) {
      mydata_input=list('', '', '', '', '')
    
    } else {
      validate(need(input$alterations, ''))
      infile <- file_input()
      
      mydata_input=read.table(infile$datapath,
                              header=TRUE, 
                              sep="\t", 
                              na.strings = "")
    }
    
    return(mydata_input)
  })
  
  names_input <- reactive({
    
    data=data_input()
    nam<-names(data[5:length(data)])
    
    
    names_input=nam[!grepl('\\..', nam)]
    
    return(names_input)
  })
  
  observe({
    updateSelectizeInput(session, 'oncgenes',
                         choices=names_input(),
                         server=TRUE #,
                         #selected=placeh_gene1
    )
  })
  
  observe({
    updateSelectizeInput(session, 'secondgene',
                         choices=names_input(),
                         server=TRUE #,
                         #selected=placeh_gene1
    )
  })
  
  observe({
    updateSelectizeInput(session, 'second3gene',
                         choices=names_input(),
                         server=TRUE #,
                         #selected=placeh_gene1
    )
  })
  
  observe({
    updateSelectizeInput(session, 'third3gene',
                         choices=names_input(),
                         server=TRUE #,
                         #selected=placeh_gene1
    )
  })
  
  mydata_bigfile <- eventReactive(input$submit, {
    
    inFile <- file_input_groups()
    inFile2 <- file_input()

    
    if (!is.null(inFile) & is.null(inFile2)) {
      
      groups <- read.table(inFile$datapath, header=TRUE, sep="\t", na.strings = "")
      
      ##### ACHTUNG COLUMN REIHENFOLGE #####
      #colnames(groups)<-c("Sample.ID", "Patient.ID", "changed")
      #groups<-groups[grep("COLO829_MATCHED_NORMAL_TISSUE|HCC1143_MATCHED_NORMAL_TISSUE|HCC1187_MATCHED_NORMAL_TISSUE|HCC1395_MATCHED_NORMAL_TISSUE|HCC1599_MATCHED_NORMAL_TISSUE|HCC1937_MATCHED_NORMAL_TISSUE|HCC1954_MATCHED_NORMAL_TISSUE|HCC2157_MATCHED_NORMAL_TISSUE|HCC2218_MATCHED_NORMAL_TISSUE|HCC38_MATCHED_NORMAL_TISSUE|HS940T_FIBROBLAST|J82_MATCHED_NORMAL_TISSUE|KMH2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE|LS1034_MATCHED_NORMAL_TISSUE|MS1_SKIN|TT_OESOPHAGUS|SIHA_CERVIX|SW626_LARGE_INTESTINE|SW954_CERVIX", groups$Patient.ID,invert = TRUE),]
      
      groups<-groups[grep("MATCHED_NORMAL_TISSUE", groups$Sample.ID, invert = TRUE),]
      groups$Patient.ID =str_replace(groups$Patient.ID, pattern = "KM-H2", replacement =  "KM-H2_l") #gibt es sonst doppelt, wenn man alle Zeichen raus nimmt und alles klein schreibt
      groups$Patient.ID =str_replace(groups$Patient.ID, pattern = "T-T", replacement =  "T-T_adt")
      
      groups$Patient.ID<-gsub("_.*", "", groups$Patient.ID)
      
      #remove all duplicates
      fff<-groups[duplicated(groups$Patient.ID),]
      groups<-subset(groups, groups$Patient.ID%in%fff$Patient.ID==FALSE)
      
      groupsx<-groups
      groups<-subset.data.frame(x = groups,select = c("Patient.ID", "changed"))
      colnames(groups)<-c("cellline","mut")
      groups[,1]=tolower(groups[,1])
      groups[,1]=gsub('[[:punct:] ]+','',groups[,1])
      groups<-groups[!is.na(groups$mut), ]
      
    } else if (!is.null(inFile2) & is.null(inFile)){
      
      #inFile2 <- input$alterations
      
      #if (is.null(inFile2))
      # return(NULL)
      
      #isolate({
      
      input$alterations
      
      #oncgene<-unique(input$oncgenes)
      
      #for (o in oncgene) {
      tbl <- read.csv(inFile2$datapath, header=TRUE, sep="\t")
      ### delete cell lines not in CCLE19 datasat -> otherways we would have dupplicates 
      #tbl<-tbl[grep("COLO829_MATCHED_NORMAL_TISSUE|HCC1143_MATCHED_NORMAL_TISSUE|HCC1187_MATCHED_NORMAL_TISSUE|HCC1395_MATCHED_NORMAL_TISSUE|HCC1599_MATCHED_NORMAL_TISSUE|HCC1937_MATCHED_NORMAL_TISSUE|HCC1954_MATCHED_NORMAL_TISSUE|HCC2157_MATCHED_NORMAL_TISSUE|HCC2218_MATCHED_NORMAL_TISSUE|HCC38_MATCHED_NORMAL_TISSUE|HS940T_FIBROBLAST|J82_MATCHED_NORMAL_TISSUE|KMH2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE|LS1034_MATCHED_NORMAL_TISSUE|MS1_SKIN|TT_OESOPHAGUS|SIHA_CERVIX|SW626_LARGE_INTESTINE|SW954_CERVIX", tbl$Patient.ID,invert = TRUE),]
      tbl<-tbl[grep("MATCHED_NORMAL_TISSUE", tbl$Sample.ID, invert = TRUE),]
      
      # delete "(driver)"
      tbl <- data.frame(lapply(tbl, gsub, pattern="\\(driver\\)", replacement=""))
      
      
      ###################################################################################################
      ######################################## Create groups files ######################################
      ###################################################################################################
      
      
      ####################################### groups file one gene ###########################################
      
      if (input$amountgenes=="One") {
        
        validate(
          need(input$oncgenes%in%colnames(tbl)==TRUE, "The entered gene name does not exist in your .tsv file.")
        )
        
        if (length(input$notprofiled)>1) {
          validate(
            need(!"EN"%in%input$notprofiled==TRUE, "Please check your settings in <Exclude cell lines not profiled in:>")
          )
        }
        
        if (length(input$notprofiled)==1 & !input$notprofiled[1]=="EN") {
          spalte<-paste0("*",input$notprofiled[1])
          z=grep(spalte, names(tbl),value = TRUE)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled"),]
        } else if (length(input$notprofiled)==2 & !input$notprofiled[1]=="EN" & !input$notprofiled[2]=="EN") {
          spalte1<-paste0("*",input$notprofiled[1])
          z=grep(spalte1, names(tbl),value = TRUE)
          spalte2<-paste0("*",input$notprofiled[2])
          zz=grep(spalte2, names(tbl),value = TRUE)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zz[1], names(tbl))]=="not profiled"),]
        } else if (length(input$notprofiled)==3 & !input$notprofiled[1]=="EN" & !input$notprofiled[2]=="EN" & !input$notprofiled[3]=="EN") {
          spalte1<-paste0("*",input$notprofiled[1])
          z=grep(spalte1, names(tbl),value = TRUE)
          spalte2<-paste0("*",input$notprofiled[2])
          zz=grep(spalte2, names(tbl),value = TRUE)
          spalte3<-paste0("*",input$notprofiled[3])
          zzz=grep(spalte3, names(tbl),value = TRUE)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zz[1], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zzz[1], names(tbl))]=="not profiled"),]
        } else if (length(input$notprofiled)==4 & !input$notprofiled[1]=="EN" & !input$notprofiled[2]=="EN" & !input$notprofiled[3]=="EN" & !input$notprofiled[4]=="EN") {
          spalte1<-paste0("*",input$notprofiled[1])
          z=grep(spalte1, names(tbl),value = TRUE)
          spalte2<-paste0("*",input$notprofiled[2])
          zz=grep(spalte2, names(tbl),value = TRUE)
          spalte3<-paste0("*",input$notprofiled[3])
          zzz=grep(spalte3, names(tbl),value = TRUE)
          spalte4<-paste0("*",input$notprofiled[4])
          zzzz=grep(spalte4, names(tbl),value = TRUE)
          print(zzzz)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zz[1], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zzz[1], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zzzz[1], names(tbl))]=="not profiled"),]
        } 
        
        if (length(input$notprofiled)==1 & !input$notprofiled[1]=="EN") {
          spalte<-paste0(input$oncgenes,input$notprofiled[1])
          tbl=subset(tbl, !tbl[spalte]=="not profiled")
        } else if (length(input$notprofiled)==2 & !input$notprofiled[1]=="EN" & !input$notprofiled[2]=="EN") {
          spalte1<-paste0(input$oncgenes,input$notprofiled[1])
          spalte2<-paste0(input$oncgenes,input$notprofiled[2])
          tbl=subset(tbl, !tbl[spalte1]=="not profiled")
          tbl=subset(tbl, !tbl[spalte2]=="not profiled")
        }
        
        tbl<-subset.data.frame(x = tbl,select = c("Sample.ID","Patient.ID", input$oncgenes))
        colnames(tbl) <- c("Sample.ID","Patient.ID", "mut")
        tbl$mut<-str_replace(string = tbl$mut,pattern = "no alteration",replacement = "wt")
        tbl<-tbl[grep("not profiled", tbl$mut,invert = TRUE),]
        
        #################################################################################################################
        #activating changes#
        
        if (input$mutationkind=="am") {
          
          tbl<-tbl[grep("LOW|HOMDEL", tbl$mut,invert = TRUE),]
          
          if ("up"%in%input$pointmut==TRUE) {
            tbl$mut<-str_replace(string = tbl$mut,pattern = "AMP",replacement = "mut")
          } else {
            tbl<-tbl[grep("AMP", tbl$mut,invert = TRUE),]
          }
          
          if ("high"%in%input$pointmut==TRUE) {
            tbl$mut<-str_replace(string = tbl$mut,pattern = "HIGH",replacement = "mut")
          } else {
            tbl<-tbl[grep("HIGH", tbl$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$pointmut==TRUE) {
            tbl$mut<-str_replace(string = tbl$mut,pattern = "_splice",replacement = "_mut@")
            tbl$mut<-gsub("@.*", "", tbl$mut)
            tbl$mut<-gsub(".*_", "", tbl$mut)
          } else {
            tbl<-tbl[grep("_splice", tbl$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$pointmut==TRUE) {
            x="*"
            tbl.2<-subset(tbl,(str_detect(string = tbl$mut,fixed(x))==TRUE))
            tbl.2$mut<-gsub(".*", "mut", tbl.2$mut)
            tbl.2$mut<-replace(tbl.2$mut,list = (grep(tbl.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl.2$Sample.ID<- as.character(tbl.2$Sample.ID)
            y<-tbl.2[["Sample.ID"]]
            tbl3<- subset(tbl, tbl$Sample.ID%in%y==FALSE)
            tbl<-rbind(tbl3, tbl.2)
          } else {
            x="*"
            tbl<-subset(tbl,(str_detect(string = tbl$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$pointmut==TRUE) {
            tbl$mut<-str_replace(string = tbl$mut,pattern = "fs",replacement = "_mut@")
            tbl$mut<-gsub("@.*", "", tbl$mut)
            tbl$mut<-gsub(".*_", "", tbl$mut)
          } else {
            tbl<-tbl[grep(".fs,.|.fs", tbl$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$pointmut==TRUE) { 
            tbl$mut<-str_replace(string = tbl$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl$mut<-gsub("@.*", "", tbl$mut)
            tbl$mut<-gsub(".*_", "", tbl$mut)
          } else {
            tbl<-tbl[grep(input$oncgenes, tbl$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$pointmut==TRUE) {
            tbl$mut<-replace(tbl$mut,list = (grep(tbl$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl<-tbl[grep("wt|mut", tbl$mut),]
            tbl$mut<-replace(tbl$mut,list = (grep(tbl$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        } 
        
        if (input$mutationkind=="dm") {
          
          tbl<-tbl[grep("AMP|HIGH", tbl$mut,invert = TRUE),]
          
          if ("down"%in%input$pointmutdm==TRUE) {
            tbl$mut<-str_replace(string = tbl$mut,pattern = "HOMDEL",replacement = "mut")
          } else {
            tbl<-tbl[grep("HOMDEL", tbl$mut,invert = TRUE),]
          }
          
          if ("low"%in%input$pointmutdm==TRUE) {
            tbl$mut<-str_replace(string = tbl$mut,pattern = "LOW",replacement = "mut")
          } else {
            tbl<-tbl[grep("LOW", tbl$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$pointmutdm==TRUE) {
            tbl$mut<-str_replace(string = tbl$mut,pattern = "_splice",replacement = "_mut@")
            tbl$mut<-gsub("@.*", "", tbl$mut)
            tbl$mut<-gsub(".*_", "", tbl$mut)
          } else {
            tbl<-tbl[grep("_splice", tbl$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$pointmutdm==TRUE) {
            x="*"
            tbl.2<-subset(tbl,(str_detect(string = tbl$mut,fixed(x))==TRUE))
            tbl.2$mut<-gsub(".*", "mut", tbl.2$mut)
            tbl.2$mut<-replace(tbl.2$mut,list = (grep(tbl.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl.2$Sample.ID<- as.character(tbl.2$Sample.ID)
            y<-tbl.2[["Sample.ID"]]
            tbl3<- subset(tbl, tbl$Sample.ID%in%y==FALSE)
            tbl<-rbind(tbl3, tbl.2)
          } else {
            x="*"
            tbl<-subset(tbl,(str_detect(string = tbl$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$pointmutdm==TRUE) {
            tbl$mut<-str_replace(string = tbl$mut,pattern = "fs",replacement = "_mut@")
            tbl$mut<-gsub("@.*", "", tbl$mut)
            tbl$mut<-gsub(".*_", "", tbl$mut)
          } else {
            tbl<-tbl[grep(".fs,.|.fs", tbl$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$pointmutdm==TRUE) { 
            tbl$mut<-str_replace(string = tbl$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl$mut<-gsub("@.*", "", tbl$mut)
            tbl$mut<-gsub(".*_", "", tbl$mut)
          } else {
            tbl<-tbl[grep(input$oncgenes, tbl$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$pointmutdm==TRUE) {
            tbl$mut<-replace(tbl$mut,list = (grep(tbl$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl<-tbl[grep("wt|mut", tbl$mut),]
            tbl$mut<-replace(tbl$mut,list = (grep(tbl$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        }
          
        tbl$Patient.ID<-gsub("_.*", "", tbl$Patient.ID)
          
        tbl$Patient.ID =str_replace(tbl$Patient.ID, pattern = "KM-H2", replacement =  "KM-H2_l") #gibt es sonst doppelt, wenn man alle Zeichen raus nimmt und alles klein schreibt
        tbl$Patient.ID =str_replace(tbl$Patient.ID, pattern = "T-T", replacement =  "T-T_adt")
        #ermove all duplicates
        fff<-tbl[duplicated(tbl$Patient.ID),]
        tbl<-subset(tbl, tbl$Patient.ID%in%fff$Patient.ID==FALSE)
        #save data frame for groups file download
        groupsx<-tbl
        #groupsfile to work with
        colnames(groupsx) <- c("Sample.ID","Patient.ID", "changed")
        tbl<-subset.data.frame(x = tbl,select = c("Patient.ID", "mut")) 
        colnames(tbl) <- c("cellline", "mut")
        tbl[,1]=tolower(tbl[,1])
        tbl[,1]=gsub('[[:punct:] ]+','',tbl[,1])
        groups<-tbl
      }
      
      ####################################### groups file two genes ###########################################
      
      if (input$amountgenes=="Two") {
        
        validate(
          need(input$oncgenes%in%colnames(tbl)==TRUE, "One of your entered gene names does not exist in your .tsv file.")
        )
        
        validate(
          need(input$secondgene%in%colnames(tbl)==TRUE, "One of your entered gene name does not exist in your .tsv file.")
        )
        
        if (length(input$notprofiled)>1) {
          validate(
            need(!"EN"%in%input$notprofiled==TRUE, "Please check your settings in <Exclude cell lines not profiled in:>")
          )
        }
        
        if (length(input$notprofiled)==1 & !input$notprofiled[1]=="EN") {
          spalte<-paste0("*",input$notprofiled[1])
          z=grep(spalte, names(tbl),value = TRUE)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled" & tbl[grep(z[2], names(tbl))]=="not profiled"),]
        } else if (length(input$notprofiled)==2 & !input$notprofiled[1]=="EN" & !input$notprofiled[2]=="EN") {
          spalte1<-paste0("*",input$notprofiled[1])
          z=grep(spalte1, names(tbl),value = TRUE)
          spalte2<-paste0("*",input$notprofiled[2])
          zz=grep(spalte2, names(tbl),value = TRUE)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled" & tbl[grep(z[2], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zz[1], names(tbl))]=="not profiled" & tbl[grep(zz[2], names(tbl))]=="not profiled"),]
        } else if (length(input$notprofiled)==3 & !input$notprofiled[1]=="EN" & !input$notprofiled[2]=="EN" & !input$notprofiled[3]=="EN") {
          spalte1<-paste0("*",input$notprofiled[1])
          z=grep(spalte1, names(tbl),value = TRUE)
          spalte2<-paste0("*",input$notprofiled[2])
          zz=grep(spalte2, names(tbl),value = TRUE)
          spalte3<-paste0("*",input$notprofiled[3])
          zzz=grep(spalte3, names(tbl),value = TRUE)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled" & tbl[grep(z[2], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zz[1], names(tbl))]=="not profiled" & tbl[grep(zz[2], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zzz[1], names(tbl))]=="not profiled" & tbl[grep(zzz[2], names(tbl))]=="not profiled"),]
        } else if (length(input$notprofiled)==4 & !input$notprofiled[1]=="EN" & !input$notprofiled[2]=="EN" & !input$notprofiled[3]=="EN" & !input$notprofiled[4]=="EN") {
          spalte1<-paste0("*",input$notprofiled[1])
          z=grep(spalte1, names(tbl),value = TRUE)
          spalte2<-paste0("*",input$notprofiled[2])
          zz=grep(spalte2, names(tbl),value = TRUE)
          spalte3<-paste0("*",input$notprofiled[3])
          zzz=grep(spalte3, names(tbl),value = TRUE)
          spalte4<-paste0("*",input$notprofiled[4])
          zzzz=grep(spalte4, names(tbl),value = TRUE)
          print(zzzz)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled" & tbl[grep(z[2], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zz[1], names(tbl))]=="not profiled" & tbl[grep(zz[2], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zzz[1], names(tbl))]=="not profiled" & tbl[grep(zzz[2], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zzzz[1], names(tbl))]=="not profiled" & tbl[grep(zzzz[2], names(tbl))]=="not profiled"),]
        } 
        
        tbl1<-subset.data.frame(x = tbl,select = c("Sample.ID","Patient.ID", input$oncgenes))
        colnames(tbl1) <- c("Sample.ID","Patient.ID", "mut")
        tbl1$mut<-str_replace(string = tbl1$mut,pattern = "no alteration",replacement = "wt")
        tbl1<-tbl1[grep("not profiled", tbl1$mut,invert = TRUE),]
        
        #################################################################################################################
        #### activating changes #### first gene ####
        
        if (input$mutationkind=="am") {
          
          tbl1<-tbl1[grep("LOW|HOMDEL", tbl1$mut,invert = TRUE),]
          
          if ("up"%in%input$pointmut==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "AMP",replacement = "mut")
          } else {
            tbl1<-tbl1[grep("AMP", tbl1$mut,invert = TRUE),]
          }
          
          if ("high"%in%input$pointmut==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "HIGH",replacement = "mut")
          } else {
            tbl1<-tbl1[grep("HIGH", tbl1$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$pointmut==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "_splice",replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep("_splice", tbl1$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$pointmut==TRUE) {
            x="*"
            tbl1.2<-subset(tbl1,(str_detect(string = tbl1$mut,fixed(x))==TRUE))
            tbl1.2$mut<-gsub(".*", "mut", tbl1.2$mut)
            tbl1.2$mut<-replace(tbl1.2$mut,list = (grep(tbl1.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl1.2$Sample.ID<- as.character(tbl1.2$Sample.ID)
            y<-tbl1.2[["Sample.ID"]]
            tbl13<- subset(tbl1, tbl1$Sample.ID%in%y==FALSE)
            tbl1<-rbind(tbl13, tbl1.2)
          } else {
            x="*"
            tbl1<-subset(tbl1,(str_detect(string = tbl1$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$pointmut==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "fs",replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep(".fs,.|.fs", tbl1$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$pointmut==TRUE) { 
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep(input$oncgenes, tbl1$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$pointmut==TRUE) {
            tbl1$mut<-replace(tbl1$mut,list = (grep(tbl1$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl1<-tbl1[grep("wt|mut", tbl1$mut),]
            tbl1$mut<-replace(tbl1$mut,list = (grep(tbl1$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        }
        
        #################################################################################################################
        #### deactivating changes #### first gene ####
        
        if (input$mutationkind=="dm") {
          
          tbl1<-tbl1[grep("AMP|HIGH", tbl1$mut,invert = TRUE),]
          
          if ("down"%in%input$pointmutdm==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "HOMDEL",replacement = "mut")
          } else {
            tbl1<-tbl1[grep("HOMDEL", tbl1$mut,invert = TRUE),]
          }
          
          if ("low"%in%input$pointmutdm==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "LOW",replacement = "mut")
          } else {
            tbl1<-tbl1[grep("LOW", tbl1$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$pointmutdm==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "_splice",replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep("_splice", tbl1$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$pointmutdm==TRUE) {
            x="*"
            tbl1.2<-subset(tbl1,(str_detect(string = tbl1$mut,fixed(x))==TRUE))
            tbl1.2$mut<-gsub(".*", "mut", tbl1.2$mut)
            tbl1.2$mut<-replace(tbl1.2$mut,list = (grep(tbl1.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl1.2$Sample.ID<- as.character(tbl1.2$Sample.ID)
            y<-tbl1.2[["Sample.ID"]]
            tbl13<- subset(tbl1, tbl1$Sample.ID%in%y==FALSE)
            tbl1<-rbind(tbl13, tbl1.2)
          } else {
            x="*"
            tbl1<-subset(tbl1,(str_detect(string = tbl1$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$pointmutdm==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "fs",replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep(".fs,.|.fs", tbl1$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$pointmutdm==TRUE) { 
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep(input$oncgenes, tbl1$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$pointmutdm==TRUE) {
            tbl1$mut<-replace(tbl1$mut,list = (grep(tbl1$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl1<-tbl1[grep("wt|mut", tbl1$mut),]
            tbl1$mut<-replace(tbl1$mut,list = (grep(tbl1$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        }
        
        tbl1$Patient.ID<-gsub("_.*", "", tbl1$Patient.ID)
        
        ############################ second gene ##############################
        
        tbl2<-subset.data.frame(x = tbl,select = c("Sample.ID","Patient.ID", input$secondgene))
        colnames(tbl2) <- c("Sample.ID","Patient.ID", "mut")
        tbl2$mut<-str_replace(string = tbl2$mut,pattern = "no alteration",replacement = "wt")
        tbl2<-tbl2[grep("not profiled", tbl2$mut,invert = TRUE),]
        
        #################################################################################################################
        #### activating changes #### second gene ####
        
        if (input$mutationkind_second=="amm") {
          
          tbl2<-tbl2[grep("LOW|HOMDEL", tbl2$mut,invert = TRUE),]
          
          if ("up"%in%input$secondpointmut==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "AMP",replacement = "mut")
          } else {
            tbl2<-tbl2[grep("AMP", tbl2$mut,invert = TRUE),]
          }
          
          if ("high"%in%input$secondpointmut==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "HIGH",replacement = "mut")
          } else {
            tbl2<-tbl2[grep("HIGH", tbl2$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$secondpointmut==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "_splice",replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep("_splice", tbl2$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$secondpointmut==TRUE) {
            x="*"
            tbl2.2<-subset(tbl2,(str_detect(string = tbl2$mut,fixed(x))==TRUE))
            tbl2.2$mut<-gsub(".*", "mut", tbl2.2$mut)
            tbl2.2$mut<-replace(tbl2.2$mut,list = (grep(tbl2.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl2.2$Sample.ID<- as.character(tbl2.2$Sample.ID)
            y<-tbl2.2[["Sample.ID"]]
            tbl23<- subset(tbl2, tbl2$Sample.ID%in%y==FALSE)
            tbl2<-rbind(tbl23, tbl2.2)
          } else {
            x="*"
            tbl2<-subset(tbl2,(str_detect(string = tbl2$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$secondpointmut==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "fs",replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep(".fs,.|.fs", tbl2$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$secondpointmut==TRUE) { 
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep(input$secondgene, tbl2$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$secondpointmut==TRUE) {
            tbl2$mut<-replace(tbl2$mut,list = (grep(tbl2$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl2<-tbl2[grep("wt|mut", tbl2$mut),]
            tbl2$mut<-replace(tbl2$mut,list = (grep(tbl2$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        }
        
        #################################################################################################################
        #### deactivating changes #### second gene ####
        
        if (input$mutationkind_second=="dmm") {
          
          tbl2<-tbl2[grep("AMP|HIGH", tbl2$mut,invert = TRUE),]
          
          if ("down"%in%input$secondpointmutdm==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "HOMDEL",replacement = "mut")
          } else {
            tbl2<-tbl2[grep("HOMDEL", tbl2$mut,invert = TRUE),]
          }
          
          if ("low"%in%input$secondpointmutdm==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "LOW",replacement = "mut")
          } else {
            tbl2<-tbl2[grep("LOW", tbl2$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$secondpointmutdm==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "_splice",replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep("_splice", tbl2$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$secondpointmutdm==TRUE) {
            x="*"
            tbl2.2<-subset(tbl2,(str_detect(string = tbl2$mut,fixed(x))==TRUE))
            tbl2.2$mut<-gsub(".*", "mut", tbl2.2$mut)
            tbl2.2$mut<-replace(tbl2.2$mut,list = (grep(tbl2.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl2.2$Sample.ID<- as.character(tbl2.2$Sample.ID)
            y<-tbl2.2[["Sample.ID"]]
            tbl23<- subset(tbl2, tbl2$Sample.ID%in%y==FALSE)
            tbl2<-rbind(tbl23, tbl2.2)
          } else {
            x="*"
            tbl2<-subset(tbl2,(str_detect(string = tbl2$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$secondpointmutdm==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "fs",replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep(".fs,.|.fs", tbl2$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$secondpointmutdm==TRUE) { 
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep(input$secondgene, tbl2$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$secondpointmutdm==TRUE) {
            tbl2$mut<-replace(tbl2$mut,list = (grep(tbl2$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl2<-tbl2[grep("wt|mut", tbl2$mut),]
            tbl2$mut<-replace(tbl2$mut,list = (grep(tbl2$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        }
        
        tbl2$Patient.ID<-gsub("_.*", "", tbl2$Patient.ID)
        
        ########### merge first and second gene ##############
        
        tbl1$Patient.ID =str_replace(tbl1$Patient.ID, pattern = "KM-H2", replacement =  "KM-H2_l")
        tbl1$Patient.ID =str_replace(tbl1$Patient.ID, pattern = "T-T", replacement =  "T-T_adt")
        #ermove all duplicates
        fff1<-tbl1[duplicated(tbl1$Patient.ID),]
        tbl1<-subset(tbl1, tbl1$Patient.ID%in%fff1$Patient.ID==FALSE)
        
        tbl2$Patient.ID =str_replace(tbl2$Patient.ID, pattern = "KM-H2", replacement =  "KM-H2_l")
        tbl2$Patient.ID =str_replace(tbl2$Patient.ID, pattern = "T-T", replacement =  "T-T_adt")
        #ermove all duplicates
        fff2<-tbl2[duplicated(tbl2$Patient.ID),]
        tbl2<-subset(tbl2, tbl2$Patient.ID%in%fff2$Patient.ID==FALSE)
        #print(tbl1)
        #print(tbl2)
        gx1<-tbl1
        gx2<-tbl2
        groupsx<-merge(gx1, gx2, by="Patient.ID", all = TRUE)
        groupsx$mut.x <- paste(groupsx$mut.x,groupsx$mut.y,sep = "_")
        groupsx$mut.y<-NULL
        
        q <- sapply(groupsx, is.factor)
        groupsx[q] <- lapply(groupsx[q], as.character)
        
        for(i in 1:length(groupsx$Sample.ID.x)) {
          if (is.na(groupsx$Sample.ID.x[i])) {
            groupsx$Sample.ID.x[i]<-(groupsx$Sample.ID.y[i])
          }
        }
        
        groupsx$Sample.ID.y<-NULL
        
        if (input$mutgroup=="comut") {
          groupsx$mut.x<-replace(groupsx$mut.x,list = (grep(groupsx$mut.x,pattern = "*wt*")),values = "wt")
          groupsx$mut.x<-replace(groupsx$mut.x,list = (grep(groupsx$mut.x,pattern = "*mut*")),values = "mut")
        } else if (input$mutgroup=="onlyone") {
          groupsx$mut.x<-replace(groupsx$mut.x,list = (grep(groupsx$mut.x,pattern = "*mut*")),values = "mut")
          groupsx$mut.x<-replace(groupsx$mut.x,list = (grep(groupsx$mut.x,pattern = "*wt*")),values = "wt")
        }
        colnames(groupsx) <- c("Patient.ID", "Sample.ID", "changed")
        groupsx<-data.frame(groupsx,stringsAsFactors = FALSE)
        groupsx<-groupsx[c(2,1,3)]
        #print(groupsx)
        
        groups<-groupsx
        
        groups<-subset.data.frame(x = groups,select = c("Patient.ID", "changed")) 
        colnames(groups) <- c("cellline", "mut")
        groups[,1]=tolower(groups[,1])
        groups[,1]=gsub('[[:punct:] ]+','',groups[,1])
        groups<-data.frame(groups,stringsAsFactors = FALSE)
      }
      
      
      ####################################### groups three genes ###########################################
      
      if (input$amountgenes=="Three") {
        
        validate(
          need(input$oncgenes%in%colnames(tbl)==TRUE, "One of your entered gene names does not exist in your .tsv file.")
        )
        
        validate(
          need(input$second3gene%in%colnames(tbl)==TRUE, "One of your entered gene name does not exist in your .tsv file.")
        )
        
        validate(
          need(input$third3gene%in%colnames(tbl)==TRUE, "One of your entered gene names does not exist in your .tsv file.")
        )
        
        if (length(input$notprofiled)>1) {
          validate(
            need(!"EN"%in%input$notprofiled==TRUE, "Please check your settings in <Exclude cell lines not profiled in:>")
          )
        }
        
        if (length(input$notprofiled)==1 & !input$notprofiled[1]=="EN") {
          spalte<-paste0("*",input$notprofiled[1])
          z=grep(spalte, names(tbl),value = TRUE)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled" & tbl[grep(z[2], names(tbl))]=="not profiled" & tbl[grep(z[3], names(tbl))]=="not profiled"),]
        } else if (length(input$notprofiled)==2 & !input$notprofiled[1]=="EN" & !input$notprofiled[2]=="EN") {
          spalte1<-paste0("*",input$notprofiled[1])
          z=grep(spalte1, names(tbl),value = TRUE)
          spalte2<-paste0("*",input$notprofiled[2])
          zz=grep(spalte2, names(tbl),value = TRUE)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled" & tbl[grep(z[2], names(tbl))]=="not profiled" & tbl[grep(z[3], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zz[1], names(tbl))]=="not profiled" & tbl[grep(zz[2], names(tbl))]=="not profiled"& tbl[grep(zz[3], names(tbl))]=="not profiled"),]
        } else if (length(input$notprofiled)==3 & !input$notprofiled[1]=="EN" & !input$notprofiled[2]=="EN" & !input$notprofiled[3]=="EN") {
          spalte1<-paste0("*",input$notprofiled[1])
          z=grep(spalte1, names(tbl),value = TRUE)
          spalte2<-paste0("*",input$notprofiled[2])
          zz=grep(spalte2, names(tbl),value = TRUE)
          spalte3<-paste0("*",input$notprofiled[3])
          zzz=grep(spalte3, names(tbl),value = TRUE)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled" & tbl[grep(z[2], names(tbl))]=="not profiled" & tbl[grep(z[3], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zz[1], names(tbl))]=="not profiled" & tbl[grep(zz[2], names(tbl))]=="not profiled"& tbl[grep(zz[3], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zzz[1], names(tbl))]=="not profiled" & tbl[grep(zzz[2], names(tbl))]=="not profiled" & tbl[grep(zzz[3], names(tbl))]=="not profiled"),]
        } else if (length(input$notprofiled)==4 & !input$notprofiled[1]=="EN" & !input$notprofiled[2]=="EN" & !input$notprofiled[3]=="EN" & !input$notprofiled[4]=="EN") {
          spalte1<-paste0("*",input$notprofiled[1])
          z=grep(spalte1, names(tbl),value = TRUE)
          spalte2<-paste0("*",input$notprofiled[2])
          zz=grep(spalte2, names(tbl),value = TRUE)
          spalte3<-paste0("*",input$notprofiled[3])
          zzz=grep(spalte3, names(tbl),value = TRUE)
          spalte4<-paste0("*",input$notprofiled[4])
          zzzz=grep(spalte4, names(tbl),value = TRUE)
          print(zzzz)
          tbl=tbl[!(tbl[grep(z[1], names(tbl))]=="not profiled" & tbl[grep(z[2], names(tbl))]=="not profiled" & tbl[grep(z[3], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zz[1], names(tbl))]=="not profiled" & tbl[grep(zz[2], names(tbl))]=="not profiled"& tbl[grep(zz[3], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zzz[1], names(tbl))]=="not profiled" & tbl[grep(zzz[2], names(tbl))]=="not profiled" & tbl[grep(zzz[3], names(tbl))]=="not profiled"),]
          tbl=tbl[!(tbl[grep(zzzz[1], names(tbl))]=="not profiled" & tbl[grep(zzzz[2], names(tbl))]=="not profiled" & tbl[grep(zzzz[3], names(tbl))]=="not profiled"),]
        }
        
        tbl1<-subset.data.frame(x = tbl,select = c("Sample.ID","Patient.ID", input$oncgenes))
        colnames(tbl1) <- c("Sample.ID","Patient.ID", "mut")
        tbl1$mut<-str_replace(string = tbl1$mut,pattern = "no alteration",replacement = "wt")
        tbl1<-tbl1[grep("not profiled", tbl1$mut,invert = TRUE),]
        
        #################################################################################################################
        #### activating changes #### first gene ####
        
        if (input$mutationkind=="am") {
          
          tbl1<-tbl1[grep("LOW|HOMDEL", tbl1$mut,invert = TRUE),]
          
          if ("up"%in%input$pointmut==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "AMP",replacement = "mut")
          } else {
            tbl1<-tbl1[grep("AMP", tbl1$mut,invert = TRUE),]
          }
          
          if ("high"%in%input$pointmut==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "HIGH",replacement = "mut")
          } else {
            tbl1<-tbl1[grep("HIGH", tbl1$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$pointmut==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "_splice",replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep("_splice", tbl1$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$pointmut==TRUE) {
            x="*"
            tbl1.2<-subset(tbl1,(str_detect(string = tbl1$mut,fixed(x))==TRUE))
            tbl1.2$mut<-gsub(".*", "mut", tbl1.2$mut)
            tbl1.2$mut<-replace(tbl1.2$mut,list = (grep(tbl1.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl1.2$Sample.ID<- as.character(tbl1.2$Sample.ID)
            y<-tbl1.2[["Sample.ID"]]
            tbl13<- subset(tbl1, tbl1$Sample.ID%in%y==FALSE)
            tbl1<-rbind(tbl13, tbl1.2)
          } else {
            x="*"
            tbl1<-subset(tbl1,(str_detect(string = tbl1$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$pointmut==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "fs",replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep(".fs,.|.fs", tbl1$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$pointmut==TRUE) { 
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep(input$oncgenes, tbl1$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$pointmut==TRUE) {
            tbl1$mut<-replace(tbl1$mut,list = (grep(tbl1$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl1<-tbl1[grep("*wt*|*mut*", tbl1$mut),]
            tbl1$mut<-replace(tbl1$mut,list = (grep(tbl1$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        }
        
        #################################################################################################################
        #### deactivating changes #### first gene ####
        
        if (input$mutationkind=="dm") {
          
          tbl1<-tbl1[grep("AMP|HIGH", tbl1$mut,invert = TRUE),]
          
          if ("down"%in%input$pointmutdm==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "HOMDEL",replacement = "mut")
          } else {
            tbl1<-tbl1[grep("HOMDEL", tbl1$mut,invert = TRUE),]
          }
          
          if ("low"%in%input$pointmutdm==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "LOW",replacement = "mut")
          } else {
            tbl1<-tbl1[grep("LOW", tbl1$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$pointmutdm==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "_splice",replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep("_splice", tbl1$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$pointmutdm==TRUE) {
            x="*"
            tbl1.2<-subset(tbl1,(str_detect(string = tbl1$mut,fixed(x))==TRUE))
            tbl1.2$mut<-gsub(".*", "mut", tbl1.2$mut)
            tbl1.2$mut<-replace(tbl1.2$mut,list = (grep(tbl1.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl1.2$Sample.ID<- as.character(tbl1.2$Sample.ID)
            y<-tbl1.2[["Sample.ID"]]
            tbl13<- subset(tbl1, tbl1$Sample.ID%in%y==FALSE)
            tbl1<-rbind(tbl13, tbl1.2)
          } else {
            x="*"
            tbl1<-subset(tbl1,(str_detect(string = tbl1$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$pointmutdm==TRUE) {
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = "fs",replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep(".fs,.|.fs", tbl1$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$pointmutdm==TRUE) { 
            tbl1$mut<-str_replace(string = tbl1$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl1$mut<-gsub("@.*", "", tbl1$mut)
            tbl1$mut<-gsub(".*_", "", tbl1$mut)
          } else {
            tbl1<-tbl1[grep(input$oncgenes, tbl1$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$pointmutdm==TRUE) {
            tbl1$mut<-replace(tbl1$mut,list = (grep(tbl1$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl1<-tbl1[grep("wt|mut", tbl1$mut),]
            tbl1$mut<-replace(tbl1$mut,list = (grep(tbl1$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        }
        
        tbl1$Patient.ID<-gsub("_.*", "", tbl1$Patient.ID)
        
        #################### second gene ####################
        
        tbl2<-subset.data.frame(x = tbl,select = c("Sample.ID","Patient.ID", input$second3gene))
        colnames(tbl2) <- c("Sample.ID","Patient.ID", "mut")
        tbl2$mut<-str_replace(string = tbl2$mut,pattern = "no alteration",replacement = "wt")
        tbl2<-tbl2[grep("not profiled", tbl2$mut,invert = TRUE),]
        
        #################################################################################################################
        #### activating changes #### second gene ####
        
        if (input$mutationkind_second3=="amm") {
          
          tbl2<-tbl2[grep("LOW|HOMDEL", tbl2$mut,invert = TRUE),]
          
          if ("up"%in%input$second3pointmut==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "AMP",replacement = "mut")
          } else {
            tbl2<-tbl2[grep("AMP", tbl2$mut,invert = TRUE),]
          }
          
          if ("high"%in%input$second3pointmut==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "HIGH",replacement = "mut")
          } else {
            tbl2<-tbl2[grep("HIGH", tbl2$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$second3pointmut==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "_splice",replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep("_splice", tbl2$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$second3pointmut==TRUE) {
            x="*"
            tbl2.2<-subset(tbl2,(str_detect(string = tbl2$mut,fixed(x))==TRUE))
            tbl2.2$mut<-gsub(".*", "mut", tbl2.2$mut)
            tbl2.2$mut<-replace(tbl2.2$mut,list = (grep(tbl2.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl2.2$Sample.ID<- as.character(tbl2.2$Sample.ID)
            y<-tbl2.2[["Sample.ID"]]
            tbl23<- subset(tbl2, tbl2$Sample.ID%in%y==FALSE)
            tbl2<-rbind(tbl23, tbl2.2)
          } else {
            x="*"
            tbl2<-subset(tbl2,(str_detect(string = tbl2$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$second3pointmut==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "fs",replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep(".fs,.|.fs", tbl2$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$second3pointmut==TRUE) { 
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep(input$second3gene, tbl2$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$second3pointmut==TRUE) {
            tbl2$mut<-replace(tbl2$mut,list = (grep(tbl2$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl2<-tbl2[grep("wt|mut", tbl2$mut),]
            tbl2$mut<-replace(tbl2$mut,list = (grep(tbl2$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        }
        
        #################################################################################################################
        #### deactivating changes #### second gene ####
        
        if (input$mutationkind_second3=="dmm") {
          
          tbl2<-tbl2[grep("AMP|HIGH", tbl2$mut,invert = TRUE),]
          
          if ("down"%in%input$second3pointmutdm==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "HOMDEL",replacement = "mut")
          } else {
            tbl2<-tbl2[grep("HOMDEL", tbl2$mut,invert = TRUE),]
          }
          
          if ("low"%in%input$second3pointmutdm==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "LOW",replacement = "mut")
          } else {
            tbl2<-tbl2[grep("LOW", tbl2$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$second3pointmutdm==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "_splice",replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep("_splice", tbl2$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$second3pointmutdm==TRUE) {
            x="*"
            tbl2.2<-subset(tbl2,(str_detect(string = tbl2$mut,fixed(x))==TRUE))
            tbl2.2$mut<-gsub(".*", "mut", tbl2.2$mut)
            tbl2.2$mut<-replace(tbl2.2$mut,list = (grep(tbl2.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl2.2$Sample.ID<- as.character(tbl2.2$Sample.ID)
            y<-tbl2.2[["Sample.ID"]]
            tbl23<- subset(tbl2, tbl2$Sample.ID%in%y==FALSE)
            tbl2<-rbind(tbl23, tbl2.2)
          } else {
            x="*"
            tbl2<-subset(tbl2,(str_detect(string = tbl2$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$second3pointmutdm==TRUE) {
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = "fs",replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep(".fs,.|.fs", tbl2$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$second3pointmutdm==TRUE) { 
            tbl2$mut<-str_replace(string = tbl2$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl2$mut<-gsub("@.*", "", tbl2$mut)
            tbl2$mut<-gsub(".*_", "", tbl2$mut)
          } else {
            tbl2<-tbl2[grep(input$second3gene, tbl2$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$second3pointmutdm==TRUE) {
            tbl2$mut<-replace(tbl2$mut,list = (grep(tbl2$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl2<-tbl2[grep("wt|mut", tbl2$mut),]
            tbl2$mut<-replace(tbl2$mut,list = (grep(tbl2$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        }
        
        tbl2$Patient.ID<-gsub("_.*", "", tbl2$Patient.ID)
        
        ############################ third gene ##############################
        
        tbl3<-subset.data.frame(x = tbl,select = c("Sample.ID","Patient.ID", input$third3gene))
        colnames(tbl3) <- c("Sample.ID","Patient.ID", "mut")
        tbl3$mut<-str_replace(string = tbl3$mut,pattern = "no alteration",replacement = "wt")
        tbl3<-tbl3[grep("not profiled", tbl3$mut,invert = TRUE),]
        
        #################################################################################################################
        #### activating changes #### third gene ####
        
        if (input$mutationkind_third3=="ammm") {
          
          tbl3<-tbl3[grep("LOW|HOMDEL", tbl3$mut,invert = TRUE),]
          
          if ("up"%in%input$third3pointmut==TRUE) {
            tbl3$mut<-str_replace(string = tbl3$mut,pattern = "AMP",replacement = "mut")
          } else {
            tbl3<-tbl3[grep("AMP", tbl3$mut,invert = TRUE),]
          }
          
          if ("high"%in%input$third3pointmut==TRUE) {
            tbl3$mut<-str_replace(string = tbl3$mut,pattern = "HIGH",replacement = "mut")
          } else {
            tbl3<-tbl3[grep("HIGH", tbl3$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$third3pointmut==TRUE) {
            tbl3$mut<-str_replace(string = tbl3$mut,pattern = "_splice",replacement = "_mut@")
            tbl3$mut<-gsub("@.*", "", tbl3$mut)
            tbl3$mut<-gsub(".*_", "", tbl3$mut)
          } else {
            tbl3<-tbl3[grep("_splice", tbl3$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$third3pointmut==TRUE) {
            x="*"
            tbl3.2<-subset(tbl3,(str_detect(string = tbl3$mut,fixed(x))==TRUE))
            tbl3.2$mut<-gsub(".*", "mut", tbl3.2$mut)
            tbl3.2$mut<-replace(tbl3.2$mut,list = (grep(tbl3.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl3.2$Sample.ID<- as.character(tbl3.2$Sample.ID)
            y<-tbl3.2[["Sample.ID"]]
            tbl33<- subset(tbl3, tbl3$Sample.ID%in%y==FALSE)
            tbl3<-rbind(tbl33, tbl3.2)
          } else {
            x="*"
            tbl3<-subset(tbl3,(str_detect(string = tbl3$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$third3pointmut==TRUE) {
            tbl3$mut<-str_replace(string = tbl3$mut,pattern = "fs",replacement = "_mut@")
            tbl3$mut<-gsub("@.*", "", tbl3$mut)
            tbl3$mut<-gsub(".*_", "", tbl3$mut)
          } else {
            tbl3<-tbl3[grep(".fs,.|.fs", tbl3$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$third3pointmut==TRUE) { 
            tbl3$mut<-str_replace(string = tbl3$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl3$mut<-gsub("@.*", "", tbl3$mut)
            tbl3$mut<-gsub(".*_", "", tbl3$mut)
          } else {
            tbl3<-tbl3[grep(input$third3gene, tbl3$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$third3pointmut==TRUE) {
            tbl3$mut<-replace(tbl3$mut,list = (grep(tbl3$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl3<-tbl3[grep("wt|mut", tbl3$mut),]
            tbl3$mut<-replace(tbl3$mut,list = (grep(tbl3$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        }
        
        #################################################################################################################
        #### deactivating changes #### third gene ####
        
        if (input$mutationkind_third3=="dmmm") {
          
          tbl3<-tbl3[grep("AMP|HIGH", tbl3$mut,invert = TRUE),]
          
          if ("down"%in%input$third3pointmutdm==TRUE) {
            tbl3$mut<-str_replace(string = tbl3$mut,pattern = "HOMDEL",replacement = "mut")
          } else {
            tbl3<-tbl3[grep("HOMDEL", tbl3$mut,invert = TRUE),]
          }
          
          if ("low"%in%input$third3pointmutdm==TRUE) {
            tbl3$mut<-str_replace(string = tbl3$mut,pattern = "LOW",replacement = "mut")
          } else {
            tbl3<-tbl3[grep("LOW", tbl3$mut,invert = TRUE),]
          }
          
          if ("sv"%in%input$third3pointmutdm==TRUE) {
            tbl3$mut<-str_replace(string = tbl3$mut,pattern = "_splice",replacement = "_mut@")
            tbl3$mut<-gsub("@.*", "", tbl3$mut)
            tbl3$mut<-gsub(".*_", "", tbl3$mut)
          } else {
            tbl3<-tbl3[grep("_splice", tbl3$mut,invert = TRUE),]
          }
          
          if ("tm"%in%input$third3pointmutdm==TRUE) {
            x="*"
            tbl3.2<-subset(tbl3,(str_detect(string = tbl3$mut,fixed(x))==TRUE))
            tbl3.2$mut<-gsub(".*", "mut", tbl3.2$mut)
            tbl3.2$mut<-replace(tbl3.2$mut,list = (grep(tbl3.2$mut,pattern = "wt",invert = TRUE)),values = "mut")
            tbl3.2$Sample.ID<- as.character(tbl3.2$Sample.ID)
            y<-tbl3.2[["Sample.ID"]]
            tbl33<- subset(tbl3, tbl3$Sample.ID%in%y==FALSE)
            tbl3<-rbind(tbl33, tbl3.2)
          } else {
            x="*"
            tbl3<-subset(tbl3,(str_detect(string = tbl3$mut,fixed(x))==FALSE))
          }
          
          if ("fs"%in%input$third3pointmutdm==TRUE) {
            tbl3$mut<-str_replace(string = tbl3$mut,pattern = "fs",replacement = "_mut@")
            tbl3$mut<-gsub("@.*", "", tbl3$mut)
            tbl3$mut<-gsub(".*_", "", tbl3$mut)
          } else {
            tbl3<-tbl3[grep(".fs,.|.fs", tbl3$mut,invert = TRUE),]
          }
          
          if ("fn"%in%input$third3pointmutdm==TRUE) { 
            tbl3$mut<-str_replace(string = tbl3$mut,pattern = input$oncgenes,replacement = "_mut@")
            tbl3$mut<-gsub("@.*", "", tbl3$mut)
            tbl3$mut<-gsub(".*_", "", tbl3$mut)
          } else {
            tbl3<-tbl3[grep(input$third3gene, tbl3$mut,invert = TRUE),]
          }
          
          if ("aaa"%in%input$third3pointmutdm==TRUE) {
            tbl3$mut<-replace(tbl3$mut,list = (grep(tbl3$mut,pattern = "wt",invert = TRUE)),values = "mut")
          } else {
            tbl3<-tbl3[grep("wt|mut", tbl3$mut),]
            tbl3$mut<-replace(tbl3$mut,list = (grep(tbl3$mut,pattern = "wt",invert = TRUE)),values = "mut")
          }
        }
        
        tbl3$Patient.ID<-gsub("_.*", "", tbl3$Patient.ID)

        ##################### merge genes #####################
        
        tbl1$Patient.ID =str_replace(tbl1$Patient.ID, pattern = "KM-H2", replacement =  "KM-H2_l")
        tbl1$Patient.ID =str_replace(tbl1$Patient.ID, pattern = "T-T", replacement =  "T-T_adt")
        #ermove all duplicates
        fff1<-tbl1[duplicated(tbl1$Patient.ID),]
        tbl1<-subset(tbl1, tbl1$Patient.ID%in%fff1$Patient.ID==FALSE)
        
        tbl2$Patient.ID =str_replace(tbl2$Patient.ID, pattern = "KM-H2", replacement =  "KM-H2_l")
        tbl2$Patient.ID =str_replace(tbl2$Patient.ID, pattern = "T-T", replacement =  "T-T_adt")
        #ermove all duplicates
        fff2<-tbl2[duplicated(tbl2$Patient.ID),]
        tbl2<-subset(tbl2, tbl2$Patient.ID%in%fff2$Patient.ID==FALSE)
        
        tbl3$Patient.ID =str_replace(tbl3$Patient.ID, pattern = "KM-H2", replacement =  "KM-H2_l")
        tbl3$Patient.ID =str_replace(tbl3$Patient.ID, pattern = "T-T", replacement =  "T-T_adt")
        #ermove all duplicates
        fff3<-tbl3[duplicated(tbl3$Patient.ID),]
        tbl3<-subset(tbl3, tbl3$Patient.ID%in%fff3$Patient.ID==FALSE)
        
        #merge data frames
        #merge first 2
        gx1<-tbl1
        gx2<-tbl2
        groupsx<-merge(gx1, gx2, by="Patient.ID", all = TRUE)
        groupsx$mut.x <- paste(groupsx$mut.x,groupsx$mut.y,sep = "_")
        groupsx$mut.y<-NULL
        
        q <- sapply(groupsx, is.factor)
        groupsx[q] <- lapply(groupsx[q], as.character)
        
        for(i in 1:length(groupsx$Sample.ID.x)) {
          if (is.na(groupsx$Sample.ID.x[i])) {
            groupsx$Sample.ID.x[i]<-(groupsx$Sample.ID.y[i])
          }
        }
        
        groupsx$Sample.ID.y<-NULL
        
        if (input$mutgroup3=="comut") {
          groupsx$mut.x<-replace(groupsx$mut.x,list = (grep(groupsx$mut.x,pattern = "*wt*")),values = "wt")
          groupsx$mut.x<-replace(groupsx$mut.x,list = (grep(groupsx$mut.x,pattern = "*mut*")),values = "mut")
        } else if (input$mutgroup3=="onlyone") {
          groupsx$mut.x<-replace(groupsx$mut.x,list = (grep(groupsx$mut.x,pattern = "*mut*")),values = "mut")
          groupsx$mut.x<-replace(groupsx$mut.x,list = (grep(groupsx$mut.x,pattern = "*wt*")),values = "wt")
        }
        colnames(groupsx) <- c("Patient.ID", "Sample.ID", "changed")
        groupsx<-data.frame(groupsx,stringsAsFactors = FALSE)
        groupsx<-groupsx[c(2,1,3)]
        #print(groupsx)
        
        #merge 2 and 3 
        gx3<-groupsx
        colnames(gx3)<-c("Sample.ID", "Patient.ID", "mut")
        gx4<-tbl3
        groupsx2<-merge(gx3, gx4, by="Patient.ID", all = TRUE)
        groupsx2$mut.x <- paste(groupsx2$mut.x,groupsx2$mut.y,sep = "_")
        groupsx2$mut.y<-NULL
        
        q <- sapply(groupsx2, is.factor)
        groupsx2[q] <- lapply(groupsx2[q], as.character)
        
        for(i in 1:length(groupsx2$Sample.ID.x)) {
          if (is.na(groupsx2$Sample.ID.x[i])) {
            groupsx2$Sample.ID.x[i]<-(groupsx2$Sample.ID.y[i])
          }
        }
        
        groupsx2$Sample.ID.y<-NULL
        
        if (input$mutgroup3=="comut") {
          groupsx2$mut.x<-replace(groupsx2$mut.x,list = (grep(groupsx2$mut.x,pattern = "*wt*")),values = "wt")
          groupsx2$mut.x<-replace(groupsx2$mut.x,list = (grep(groupsx2$mut.x,pattern = "*mut*")),values = "mut")
        } else if (input$mutgroup3=="onlyone") {
          groupsx2$mut.x<-replace(groupsx2$mut.x,list = (grep(groupsx2$mut.x,pattern = "*mut*")),values = "mut")
          groupsx2$mut.x<-replace(groupsx2$mut.x,list = (grep(groupsx2$mut.x,pattern = "*wt*")),values = "wt")
        }
        colnames(groupsx2) <- c("Patient.ID", "Sample.ID", "changed")
        groupsx2<-data.frame(groupsx2,stringsAsFactors = FALSE)
        groupsx2<-groupsx2[c(2,1,3)]
        groupsx<-groupsx2
        #print(groupsx)
        
        groups<-groupsx
        
        groups<-subset.data.frame(x = groups,select = c("Patient.ID", "changed")) 
        colnames(groups) <- c("cellline", "mut")
        groups[,1]=tolower(groups[,1])
        groups[,1]=gsub('[[:punct:] ]+','',groups[,1])
        groups<-data.frame(groups,stringsAsFactors = FALSE)
        
      }
    }
    
    print(unique(groups$mut))
    print(groups)
    validate(
      need(exists("groups")==TRUE, "Please check that ONE file is uploaded.")
    )
    
    validate(
      need(length(unique(groups$mut))==2, "Please check column <changed> in your groups.txt file. It consists of more words than <mut>, <wt> and empty spaces for NA values."),
      need("wt"%in%groups$mut==TRUE, "Please check column <changed> in your groups.txt file. It consists of more words than <mut>, <wt> and empty spaces for NA values."),
      need("mut"%in%groups$mut==TRUE, "Please check column <changed> in your groups.txt file. It consists of more words than <mut>, <wt> and empty spaces for NA values.")
    )
    
    ############################################################################################################################
    ####################################### Expanded Drug Response #############################################################
    ############################################################################################################################
    
    drug_groups=read.table(file="20200810_drug_groups.txt",header=TRUE,sep="\t",fill=TRUE)
    unique_drug_groupe<-unique(drug_groups$target_pathway)
    
    tab=read.table(na.strings = "", file="20200408_expanded_drug_data.csv", header=TRUE,sep=",",dec=".",fill=TRUE,check.names = FALSE)
    
    origin=unique(tab$GDSC_Tissue_descriptor_1)
    
    ###### include for analyzing melanoma celllines only ######
    
    if (length(input$tissue)>1) {
      validate(
        need(!"EN"%in%input$tissue==TRUE, "Please check your settings in <Analyze specific tissue origins only?>")
      )
    }
    
    if (!input$tissue=="EN") {
      #skin_subgroup<-read.table("skin_subgroup.txt", header=TRUE,sep = "\t",fill=TRUE)
      #skin_subgroup[,1]=tolower(skin_subgroup[,1])
      #skin_subgroup[,1]=gsub('[[:punct:] ]+','',skin_subgroup[,1])
      tab=subset(tab, tab$GDSC_Tissue_descriptor_1%in%input$tissue==TRUE)
    }
    ###### general infos for infofile #######
    
    n=length(tab[1,])
    print(paste(n,"/2 samples"))
    
    ### Total amount of celllines ####
    num<-subset(tab, tab$Screened_Compounds%in%groups$cellline==TRUE)
    n<-as.numeric(length(num[,1]))
    ### Amount of celllines wt ###
    numwt<-subset(groups, groups$mut=="wt")
    num_wt<-subset(tab, tab$Screened_Compounds%in%numwt$cellline==TRUE)
    m<-length(num_wt[,1])
    ### Amount of celllines mut ###
    nummut<-subset(groups, groups$mut=="mut")
    num_mut<-subset(tab, tab$Screened_Compounds%in%nummut$cellline==TRUE)
    o<-length(num_mut[,1])
    
    ### perdent mutated ###
    overallmut<-(o/n*100)
    overallmut<-round(overallmut, digits = 2)
    
    ### percentage distribution of cell lines origins and percentage mutated per origin ###
    percent_list<-list() #overall percentage origin
    percent_plot_list<-list() #overall percentage origin for plotting
    mut_perc_list<-list() #percentage mutated per origin
    
    for (orig in origin) {
      
      ### percent origin ###
      nochwas<-subset(num, num$GDSC_Tissue_descriptor_1==orig)
      p=as.numeric(length(nochwas[,1]))
      prozent=(p/n*100)
      
      ### n(wt) ### 
      nochwaswt<-subset(num_wt, num_wt$GDSC_Tissue_descriptor_1==orig)
      q<-as.numeric(length(nochwaswt[,1]))
      
      ### n(mut) ###
      nochwasmut<-subset(num_mut, num_mut$GDSC_Tissue_descriptor_1==orig)
      r<-as.numeric(length(nochwasmut[,1]))
      
      ### percent cell lines mutated per origin ###
      prozent_mut<-(r/p*100)
      
      print(paste(orig, p, prozent))
      
      ### overall percentage origin ###
      dat<-c(orig, p, prozent)
      ### overall percentage origin for plotting ###
      dat0<-(c("overall", orig, "NA", p, prozent))
      ### Amount of cellines mutated per origin ### 
      dat1<-c(orig, q,r,prozent_mut)
      
      percent_plot_list[[orig]]<- dat0
      percent_list[[orig]] <- dat
      mut_perc_list[[orig]]<- dat1
    }
    
    ### create ordered dataframe overall percentage distribution per origin ###
    percent = do.call(rbind, percent_list)
    percent<-data.frame(percent, stringsAsFactors = FALSE)
    percent$X2<-as.numeric(percent$X2)
    percent<-percent[order(percent$X2, decreasing = TRUE),]
    
    ### create dataframe overall percentage distribution per origin for plotting ###
    percent_plot = do.call(rbind, percent_plot_list)
    percent_plot<-data.frame(percent_plot, stringsAsFactors = FALSE)
    
    ### create ordered dataframe percent mutated per origin ###
    mut_perc=do.call(rbind, mut_perc_list)
    mut_perc<-data.frame(mut_perc, stringsAsFactors = FALSE)
    mut_perc$X4<-as.numeric(mut_perc$X4)
    mut_perc<-mut_perc[order(mut_perc$X4, decreasing = TRUE),]
    
    
    ######## looping through all compounds #########
    
    for (i in colnames(tab[2:251])) {
      
      compound<-i
      print(i) 
      
      ### subset all sensitive per compound ###
      temp_tab_S=subset(tab, tab[[i]]=="S")	
      nS=length(temp_tab_S[,1])
      
      ### subset all resistant per compound ###
      temp_tab_R=subset(tab, tab[[i]]=="R")	
      nR=length(temp_tab_R[,1]) 
      
      mutorwt=unique(groups$mut)
      print(mutorwt)
      
      ### add up amound of wt/mut r/s celllines per compound ###
      for (type in mutorwt) {
        type_group=subset(groups, groups$mut==type)
        type_group=type_group[,1]
        tab_type_S=subset(temp_tab_S, temp_tab_S$Screened_Compounds%in%type_group==TRUE)
        n_type_S=length(tab_type_S[,1])
        
        tab_type_R=subset(temp_tab_R, temp_tab_R$Screened_Compounds%in%type_group==TRUE)
        n_type_R=length(tab_type_R[,1])
        
        #output=print(paste(compound, type, n_type_S, n_type_R))
        ###
        screened_compounds<-data.frame(compound, type, n_type_S, n_type_R,stringsAsFactors = FALSE)
        #screened_compounds<-matrix(screened_compounds=list(compound, type, n_type_S, n_type_R),nrow=1,ncol=2)
        #screened_compounds <- as.data.frame(screened_compounds,stringsAsFactors = FALSE)
        colnames(screened_compounds) <- c("Screened_Compounds","wt_or_mut", "sensitive", "resistant")
        if (exists("Screened_Compounds")==TRUE) {
          Screened_Compounds <- rbind(Screened_Compounds, screened_compounds)
        } else {
          Screened_Compounds <- screened_compounds
        }
        ###
        #write.table(output,file = "Screened_Compounds.txt",append = TRUE,col.names = FALSE, row.names = FALSE, sep = "\t",quote = FALSE)
        
        if (type=="wt") {
          wt_S<-n_type_S
          data_wt_S<-tab_type_S
          wt_R<-n_type_R 
          data_wt_R<-tab_type_R
        }
        
        if (type=="mut") {
          mut_S<-n_type_S
          data_mut_S<-tab_type_S
          mut_R<-n_type_R
          data_mut_R<-tab_type_R
        }
      }
      
      
      #i matchen mit target pathway
      tp<-subset(drug_groups, drug_groups$Screened_Compounds==i)
      tp<-tp$target_pathway[1]
      
      ### calculate odds ratios and fisher's exact test ###
      if (wt_S>0 & wt_R>0 & mut_S>0 & mut_R>0) {
        
        ### twosided 95% ###
        odds_ratio=oddsratio( c(mut_R, wt_R, mut_S, wt_S), conf.level=0.95) #turne OR #mut_r, wt_r, mut_s, wt_s
        odds_ratio_p=odds_ratio$p.value[6] #chi.square p value
        if (odds_ratio_p<0.1) {
          odds_ratio_p<-formatC(odds_ratio_p, format="e", digits = 3) 
        } else {
          odds_ratio_p<-round(odds_ratio_p, digits = 4)}
        fisher_pvalue=odds_ratio$p.value[4] #fisher p value
        if (fisher_pvalue<0.1) {
          fisher_pvalue<-formatC(fisher_pvalue, format="e", digits = 3) 
        } else {
          fisher_pvalue<-round(fisher_pvalue, digits = 4)}
        Odds_Ratio=odds_ratio$measure[2] #odds ratio
        if (Odds_Ratio<0.1) {
          Odds_Ratio<-formatC(Odds_Ratio, format="e", digits = 3) 
        } else {
          Odds_Ratio<-round(Odds_Ratio, digits = 4)}
        Lower=odds_ratio$measure[4] #lower
        if (Lower<0.1) {
          Lower<-formatC(Lower, format="e", digits = 3) 
        } else {
          Lower<-round(Lower, digits = 4)}
        Upper=odds_ratio$measure[6] #upper
        if (Upper<0.1) {
          Upper<-formatC(Upper, format="e", digits = 3) 
        } else {
          Upper<-round(Upper, digits = 4)}
        compound_odds_ratios<-data.frame(compound, tp, Odds_Ratio, fisher_pvalue, odds_ratio_p, Lower, Upper,stringsAsFactors = FALSE)
        colnames(compound_odds_ratios) <- c("Screened_Compounds","target_pathway","OR","fisher_p","chi_p", "lower", "upper")
        if (exists("Compound_odds_ratios")==TRUE) {
          Compound_odds_ratios <- rbind(Compound_odds_ratios, compound_odds_ratios)
        } else {
          Compound_odds_ratios <- compound_odds_ratios
        }
      } else { 
        ### if a group of 0 exists -> OR not possible -> only calculating fisher's exact test ###
        tafel<-matrix(c(mut_S, mut_R, wt_S, wt_R), nrow = 2)
        fis=fisher.test(tafel, conf.level = 0.95)
        fisher_p=fis$p.value
        if (mut_R<1) {
          compound_odds_ratios<-data.frame(compound, tp, 0, fisher_p, "NA", 0, 0,stringsAsFactors = FALSE)
        } else if (mut_S<1) {
          compound_odds_ratios<-data.frame(compound, tp, 7, fisher_p, "NA", 7, 7,stringsAsFactors = FALSE)
        } else {
          compound_odds_ratios<-data.frame(compound, tp, "NA", fisher_p, "NA", "NA", "NA",stringsAsFactors = FALSE)
        }
        colnames(compound_odds_ratios) <- c("Screened_Compounds","target_pathway","OR","fisher_p","chi_p", "lower", "upper")
        if (exists("Compound_odds_ratios")==TRUE) {
          Compound_odds_ratios <- rbind(Compound_odds_ratios, compound_odds_ratios)
        } else {
          Compound_odds_ratios <- compound_odds_ratios
        }
      }
    } #compound
    
    
    ####### same for drug groups #######
    
    ### read data ###
    print(unique_drug_groupe)
    result_data=Screened_Compounds
    
    ### odds ratios per drug group ###
    for (druggroup in unique_drug_groupe) {
      drug_subset<-subset(drug_groups, drug_groups$target_pathway==druggroup) #alle compouds einer Medikamentengruppe
      drug_subset<-drug_subset[,1]
      tab_subset<-subset(result_data, result_data$Screened_Compounds%in%drug_subset==TRUE)
      #print(tab_subset)
      subset_wt_mut<-unique(tab_subset$wt_or_mut)
      
      for (j in subset_wt_mut) {
        tab_subset_wt_mut=subset(tab_subset, tab_subset$wt_or_mut==j)
        sum_S<-sum(as.numeric(tab_subset_wt_mut$sensitive))
        sum_R<-sum(as.numeric(tab_subset_wt_mut$resistant))
        percent_S<-(sum_S/(sum_S+sum_R)*100)
        percent_R<-(sum_R/(sum_S+sum_R)*100)
        druggroups<-data.frame(druggroup, j, "sensitive",sum_S, percent_S, stringsAsFactors = FALSE)
        colnames(druggroups) <- c("drug_group", "wt_or_mut", "response", "sum", "percent")
        if (exists("Druggroups")==TRUE) {
          Druggroups <- rbind(Druggroups, druggroups)
        } else {
          Druggroups <- druggroups
        }
        druggroups<-data.frame(druggroup, j, "resistant", sum_R, percent_R, stringsAsFactors = FALSE)
        colnames(druggroups) <- c("drug_group", "wt_or_mut", "response", "sum", "percent")
        if (exists("Druggroups")==TRUE) {
          Druggroups <- rbind(Druggroups, druggroups)
        } else {
          Druggroups <- druggroups
        }
        
        if (j=="wt") {
          wt_s<-sum_S
          wt_r<-sum_R }
        
        if (j=="mut") {
          mut_s<-sum_S
          mut_r<-sum_R }
      }
      
      if (wt_s>0 & wt_r>0 & mut_s>0 & mut_r>0) {
        ## Calculate ORs #twosided
        or=oddsratio(c(mut_r, wt_r, mut_s, wt_s), conf.level=0.95)
        #	  print(or)
        OR_p=or$p.value[6] #chi.square p value
        if (OR_p<0.1) {
          OR_p<-formatC(OR_p, format="e", digits = 3) 
        } else {
          OR_p<-round(OR_p, digits = 4)}
        fisher_pv=or$p.value[4] #fisher p value
        if (fisher_pv<0.1) {
          fisher_pv<-formatC(fisher_pv, format="e", digits = 3) 
        } else {
          fisher_pv<-round(fisher_pv, digits = 4)}
        OR=or$measure[2]
        if (OR<0.1) {
          OR<-formatC(OR, format="e", digits = 3) 
        } else {
          OR<-round(OR, digits = 4)}
        lower=or$measure[4]
        if (lower<0.1) {
          lower<-formatC(lower, format="e", digits = 3) 
        } else {
          lower<-round(lower, digits = 4)}
        upper=or$measure[6]
        if (upper<0.1) {
          upper<-formatC(upper, format="e", digits = 3) 
        } else {
          upper<-round(upper, digits = 4)}
        odds_ratios<-data.frame(druggroup, OR, fisher_pv, OR_p, lower, upper, stringsAsFactors = FALSE)
        colnames(odds_ratios) <- c("drug_group","OR","fisher_p","chi_p","lower", "upper")
        if (exists("Odds_ratios")==TRUE) {
          Odds_ratios <- rbind(Odds_ratios, odds_ratios)
        } else {
          Odds_ratios <- odds_ratios
        }
        
      } else { 
        
        tafel1<-matrix(c(mut_r, mut_s, wt_r, wt_s), nrow = 2)
        fis=fisher.test(tafel1, conf.level = 0.95)
        fisherp=fis$p.value
        if (fisherp<0.1) {
          fisherp<-formatC(fisherp, format="e", digits = 3) 
        } else {
          fisherp<-round(fisherp, digits = 4)}
        odds_ratios<-data.frame(druggroup, "NA",fisherp , "NA", "NA", "NA", stringsAsFactors = FALSE)
        colnames(odds_ratios) <- c("drug_group","OR","fisher_p","chi_p","lower", "upper")
        if (exists("Odds_ratios")==TRUE) {
          Odds_ratios <- rbind(Odds_ratios, odds_ratios)
        } else {
          Odds_ratios <- odds_ratios
        }
      }
    } #drug group 
    
    #ddd <- data.frame(lapply(ddd, as.numeric), stringsAsFactors=FALSE)
    
    ### crating list of data frames ###
    Screened_Compounds<-data.frame(Screened_Compounds,stringsAsFactors = FALSE)
    Druggroups<-data.frame(Druggroups,stringsAsFactors = FALSE)
    Compound_odds_ratios<-data.frame(Compound_odds_ratios,stringsAsFactors = FALSE)
    Compound_odds_ratios$fisher_p<-as.numeric(Compound_odds_ratios$fisher_p)
    Compound_odds_ratios$chi_p<-as.numeric(Compound_odds_ratios$chi_p)
    Compound_odds_ratios$OR<-as.numeric(Compound_odds_ratios$OR)
    Compound_odds_ratios$lower<-as.numeric(Compound_odds_ratios$lower)
    Compound_odds_ratios$upper<-as.numeric(Compound_odds_ratios$upper)
    Odds_ratios<-data.frame(Odds_ratios,stringsAsFactors = FALSE)
    Odds_ratios$fisher_p<-as.numeric(Odds_ratios$fisher_p)
    Odds_ratios$chi_p<-as.numeric(Odds_ratios$chi_p)
    Odds_ratios$OR<-as.numeric(Odds_ratios$OR)
    Odds_ratios$lower<-as.numeric(Odds_ratios$lower)
    Odds_ratios$upper<-as.numeric(Odds_ratios$upper)
    str(Screened_Compounds)
    str(Druggroups)
    str(Compound_odds_ratios)
    str(Odds_ratios)
    plot_data_new_mydata_bigfile <- list("OR_druggroup"=Odds_ratios, "OR_compounds"=Compound_odds_ratios, "drug_groups" = Druggroups, "screened_compounds" = Screened_Compounds, "groups file"=groupsx) 
    
    #plot_data_new_mydata_bigfile <- list("groups file"=groupsx, "screened_compounds" = Screened_Compounds, "drug_groups" = Druggroups, "OR_compounds"=Compound_odds_ratios, "OR_druggroup"=Odds_ratios) 
    #write_xlsx(sheets, "odds_ratio_data.xlsx")
    
    
    #})
    
    plot_data_new_mydata_bigfile
  })
  
  observeEvent(input$reset, {
    reset("alterations")
  })
  
  observeEvent(input$RESET, {
    reset('groupstxt')
  })
  
  drug_results<-function(){
    drug<- mydata_bigfile()
    drug<-data.frame(Reduce(rbind, drug[1]))
    
    drug
    
  }
  
  output$Table_druggroups<-DT::renderDataTable(drug_results(), options = list(pageLength = 25))
  
  groupsfile<-function(){
    groups<- mydata_bigfile()
    groups<-data.frame(Reduce(rbind, groups[5]))
    
    groups
    
  }
  
  
  #########################################################################################################################
  ############################################## Plot hight Boxplot #######################################################
  #########################################################################################################################
  
  bpdata<-eventReactive(input$submit, {
    neededdata<- mydata_bigfile()
    neededdata<-data.frame(Reduce(rbind, neededdata[5]))
    
    groups<-subset.data.frame(x = neededdata,select = c("Patient.ID", "changed")) 
    
    groups<-groups[grep("COLO829_MATCHED_NORMAL_TISSUE|HCC1143_MATCHED_NORMAL_TISSUE|HCC1187_MATCHED_NORMAL_TISSUE|HCC1395_MATCHED_NORMAL_TISSUE|HCC1599_MATCHED_NORMAL_TISSUE|HCC1937_MATCHED_NORMAL_TISSUE|HCC1954_MATCHED_NORMAL_TISSUE|HCC2157_MATCHED_NORMAL_TISSUE|HCC2218_MATCHED_NORMAL_TISSUE|HCC38_MATCHED_NORMAL_TISSUE|HS940T_FIBROBLAST|J82_MATCHED_NORMAL_TISSUE|KMH2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE|LS1034_MATCHED_NORMAL_TISSUE|MS1_SKIN|TT_OESOPHAGUS|SIHA_CERVIX|SW626_LARGE_INTESTINE|SW954_CERVIX", groups$Patient.ID,invert = TRUE),]
    groups$Patient.ID =str_replace(groups$Patient.ID, pattern = "KM-H2", replacement =  "KM-H2_l") #gibt es sonst doppelt, wenn man alle Zeichen raus nimmt und alles klein schreibt
    groups$Patient.ID =str_replace(groups$Patient.ID, pattern = "T-T", replacement =  "T-T_adt")
    groups[,1]=tolower(groups[,1])
    groups[,1]=gsub('[[:punct:] ]+','',groups[,1])
    ##
    colnames(groups)<-c("Screened_Compounds", "mut")
    
    IC50<-read.csv("20200408_IC50.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    
    IC50<-IC50[grep("COLO829_MATCHED_NORMAL_TISSUE|HCC1143_MATCHED_NORMAL_TISSUE|HCC1187_MATCHED_NORMAL_TISSUE|HCC1395_MATCHED_NORMAL_TISSUE|HCC1599_MATCHED_NORMAL_TISSUE|HCC1937_MATCHED_NORMAL_TISSUE|HCC1954_MATCHED_NORMAL_TISSUE|HCC2157_MATCHED_NORMAL_TISSUE|HCC2218_MATCHED_NORMAL_TISSUE|HCC38_MATCHED_NORMAL_TISSUE|HS940T_FIBROBLAST|J82_MATCHED_NORMAL_TISSUE|KMH2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE|LS1034_MATCHED_NORMAL_TISSUE|MS1_SKIN|TT_OESOPHAGUS|SIHA_CERVIX|SW626_LARGE_INTESTINE|SW954_CERVIX", IC50$Screened_Compounds,invert = TRUE),]
    IC50$Screened_Compounds =str_replace(IC50$Screened_Compounds, pattern = "KM-H2", replacement =  "KM-H2_l") #gibt es sonst doppelt, wenn man alle Zeichen raus nimmt und alles klein schreibt
    IC50$Screened_Compounds =str_replace(IC50$Screened_Compounds, pattern = "T-T", replacement =  "T-T_adt")
    IC50[,1]=tolower(IC50[,1])
    IC50[,1]=gsub('[[:punct:] ]+','',IC50[,1])
    
    tab=read.table(na.strings = "", file="20200408_expanded_drug_data.csv", header=TRUE,sep=",",dec=".",fill=TRUE,check.names = FALSE)
    x<-subset.data.frame(tab,select = c("Screened_Compounds", "GDSC_Tissue_descriptor_1"))
    IC50<-merge(IC50, x, by="Screened_Compounds")
    
    if (length(input$tissue)>1) {
      validate(
        need(!"EN"%in%input$tissue==TRUE, "Please check your settings in <Analyze specific tissue origins only?>")
      )
    }
    
    if (!input$tissue=="EN") {
      IC50=subset(IC50, IC50$GDSC_Tissue_descriptor_1%in%input$tissue==TRUE)
    }
    
    boxplotdata<-merge(groups, IC50, by="Screened_Compounds")
  })
  
  plothight_boxplot<- function(){
    
    boxplotdata<-bpdata()
    #sum(is.na(boxplotdata$Crizotinib)) #Anzahl an NAs gecheckt -> enspricht Anzahl an removed rows in warnmeldung beim ploten
    colnames(boxplotdata)=gsub('[[:punct:] ]+','',colnames(boxplotdata))
    #print(boxplotdata)
    p <- list()
    colNames <- names(boxplotdata)[3:252] #nicht notwendig, auch direkt mit colnames im loop möglich
    for(i in colNames) {
      colnames(boxplotdata)[colnames(boxplotdata) == i] <- paste0("IC50_", i)
      j<-paste0("IC50_", i)
      
      tab_wt<-subset(boxplotdata, boxplotdata$mut=="wt")
      n_wt<- sum(is.na(tab_wt[[j]])==FALSE)
      #print(n_wt)
      tab_mut<-subset(boxplotdata, boxplotdata$mut=="mut")
      n_mut<- sum(is.na(tab_mut[[j]])==FALSE)
      #print(n_mut)
      
      w=tryCatch({wilcox.test(tab_wt[[j]],tab_mut[[j]], alternative="two.sided")}, error=function(err) { print(paste("Wilcox ERROR: ",err)); w=NA } )  
      
      if(is.na(w)!=TRUE) {
        p_val=round(w$p.value,5) }
      
      Significant<-p_val
      
      if (exists("significant")==TRUE) {
        significant <- rbind(significant, Significant)
      } else {
        significant <- Significant
      }
    }
    
    significant<- subset(significant,significant[,1]<0.05)
    n=as.numeric(length(significant[,1]))
    print(n)
    if (n>6) {
      plothight_boxplot<-(n*50)
    } else {
      plothight_boxplot<-250
    }
  }
  
  downloadplothight_boxplot<- function(){
    x<-plothight_boxplot()
    downloadplothight<-(x/40)
  }
  
  #########################################################################################################################
  ################################################## Output Boxplot #######################################################
  #########################################################################################################################
  
  
  output$BoxPlot <- renderPlot({
    
    boxplotdata<-bpdata()
    #sum(is.na(boxplotdata$Crizotinib)) #Anzahl an NAs gecheckt -> enspricht Anzahl an removed rows in warnmeldung beim ploten
    colnames(boxplotdata)=gsub('[[:punct:] ]+','',colnames(boxplotdata))
    #print(boxplotdata)
    p <- list()
    colNames <- names(boxplotdata)[3:252] #nicht notwendig, auch direkt mit colnames im loop möglich
    for(i in colNames) {
      colnames(boxplotdata)[colnames(boxplotdata) == i] <- paste0("IC50_", i)
      j<-paste0("IC50_", i)
      
      tab_wt<-subset(boxplotdata, boxplotdata$mut=="wt")
      n_wt<- sum(is.na(tab_wt[[j]])==FALSE)
      #print(n_wt)
      tab_mut<-subset(boxplotdata, boxplotdata$mut=="mut")
      n_mut<- sum(is.na(tab_mut[[j]])==FALSE)
      #print(n_mut)
      
      w=tryCatch({wilcox.test(tab_wt[[j]],tab_mut[[j]], alternative="two.sided")}, error=function(err) { print(paste("Wilcox ERROR: ",err)); w=NA } )  
      
      if (is.na(w)!=TRUE) {
        p_val=round(w$p.value,5)
        #print(paste(j,":",p_val))
        
        if (p_val<0.05) {
          p[[i]] <- ggplot(boxplotdata, aes_string(x="mut", y = j, fill="mut")) +
            geom_boxplot()+
            scale_x_discrete(name = paste("p =", p_val), limits=c("wt","mut"))+
            scale_fill_manual(values=c("red", "gray"))+
            labs(caption=paste0("n=", n_wt, " vs. n=", n_mut)) +
            theme(legend.position = "none",
                  plot.caption=element_text(size=15, hjust=0),
                  axis.text=element_text(size=13),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 15))
        }
      }
    }
    #test<-do.call(grid.arrange,p) #geht auch
    test<-grid.arrange(grobs = p, ncol = 6)
    test
  }, height = plothight_boxplot)
  
  #########################################################################################################################
  ############################################## Boxplot for download #####################################################
  #########################################################################################################################
  
  boxPlot <- function() {
    boxplotdata<-bpdata()
    #sum(is.na(boxplotdata$Crizotinib)) #Anzahl an NAs gecheckt -> enspricht Anzahl an removed rows in warnmeldung beim ploten
    colnames(boxplotdata)=gsub('[[:punct:] ]+','',colnames(boxplotdata))
    #print(boxplotdata)
    p <- list()
    colNames <- names(boxplotdata)[3:252] #nicht notwendig, auch direkt mit colnames im loop möglich
    for(i in colNames) {
      colnames(boxplotdata)[colnames(boxplotdata) == i] <- paste0("IC50_", i)
      j<-paste0("IC50_", i)
      
      tab_wt<-subset(boxplotdata, boxplotdata$mut=="wt")
      n_wt<- sum(is.na(tab_wt[[j]])==FALSE)
      #print(n_wt)
      tab_mut<-subset(boxplotdata, boxplotdata$mut=="mut")
      n_mut<- sum(is.na(tab_mut[[j]])==FALSE)
      #print(n_mut)
      
      w=tryCatch({wilcox.test(tab_wt[[j]],tab_mut[[j]], alternative="two.sided")}, error=function(err) { print(paste("Wilcox ERROR: ",err)); w=NA } )  
      
      if (is.na(w)!=TRUE) {
        p_val=round(w$p.value,5)
        #print(paste(j,":",p_val))
        
        if (p_val<0.05) {
          p[[i]] <- ggplot(boxplotdata, aes_string(x="mut", y = j, fill="mut")) +
            geom_boxplot()+
            scale_x_discrete(name = paste("p =", p_val), limits=c("wt","mut"))+
            scale_fill_manual(values=c("red", "gray"))+
            labs(caption=paste0("n=", n_wt, " vs. n=", n_mut)) +
            theme(legend.position = "none",
                  plot.caption=element_text(size=11, hjust=0))
        }
      }
    }
    #test<-do.call(grid.arrange,p)
    test<-grid.arrange(grobs = p, ncol = 6)
    
    test
  }
  
  #########################################################################################################################
  ############################################## Boxplot Wilcox Table #####################################################
  ######################################################################################################################### 
  
  Wilcoxtable<-function(){
    
    
    boxplotdata<-bpdata()
    #sum(is.na(boxplotdata$Crizotinib)) #Anzahl an NAs gecheckt -> enspricht Anzahl an removed rows in warnmeldung beim ploten
    colnames(boxplotdata)=gsub('[[:punct:] ]+','',colnames(boxplotdata))
    #print(boxplotdata)
    #p <- list()
    
    for(i in colnames(boxplotdata)[3:252]) {
      colnames(boxplotdata)[colnames(boxplotdata) == i] <- paste0("IC50_", i)
      j<-paste0("IC50_", i)
      
      tab_wt<-subset(boxplotdata, boxplotdata$mut=="wt")
      n_wt<- sum(is.na(tab_wt[[j]])==FALSE)
      #print(n_wt)
      tab_mut<-subset(boxplotdata, boxplotdata$mut=="mut")
      n_mut<- sum(is.na(tab_mut[[j]])==FALSE)
      #print(n_mut)
            #upper_wt<-quantile(tab_wt, 0.25)
      median_wt<-quantile(tab_wt[[j]], 0.5, na.rm = TRUE)
      median_mut<-quantile(tab_mut[[j]], 0.5, na.rm = TRUE)
      
      w=tryCatch({wilcox.test(tab_wt[[j]],tab_mut[[j]], alternative="two.sided")}, error=function(err) { print(paste("Wilcox ERROR: ",err)); w=NA } )  
      
      if (is.na(w)!=TRUE) {
        p_val=round(w$p.value,5)
      } else {
        p_val="NA"
      }
      
      Significant<-data.frame(i, n_wt, n_mut, median_wt, median_mut, p_val, stringsAsFactors = FALSE)
      colnames(Significant) <- c("compound","n_wt","n_mut", "median_wt", "median_mut", "wilcox_p-value")
      
      if (exists("significant")==TRUE) {
        significant <- rbind(significant, Significant)
      } else {
        significant <- Significant
      }
    }
    significant
  }
  
  #########################################################################################################################
  ################################# plot Expanded Drug Response for download ##############################################
  #########################################################################################################################
  
  myplot<- function() {
    
    plot_data_new_list <- mydata_bigfile()
    dg_or<-data.frame(Reduce(rbind, plot_data_new_list[1]))
    dg<-data.frame(Reduce(rbind, plot_data_new_list[3]))
    merge_data <- merge(dg, dg_or, by ="drug_group")
    plot_data_new<-merge_data
    
    if (is.na(plot_data_new$OR[1])==FALSE) {
      uno<-paste0("OR=",plot_data_new$OR[1])
    } else {
      uno<-paste0(round(plot_data_new$percent[4], digits = 2),"%"," vs ",round(plot_data_new$percent[2], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[5])==FALSE) {
      dos<-paste0("OR=",plot_data_new$OR[5])
    } else {
      dos<-paste0(round(plot_data_new$percent[8], digits = 2),"%"," vs ",round(plot_data_new$percent[6], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[9])==FALSE) {
      tres<-paste0("OR=",plot_data_new$OR[9])
    } else {
      tres<-paste0(round(plot_data_new$percent[12], digits = 2),"%"," vs ",round(plot_data_new$percent[10], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[13])==FALSE) {
      quatro<-paste0("OR=",plot_data_new$OR[13])
    } else {
      quatro<-paste0(round(plot_data_new$percent[16], digits = 2),"%"," vs ",round(plot_data_new$percent[14], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[17])==FALSE) {
      cinco<-paste0("OR=",plot_data_new$OR[17])
    } else {
      cinco<-paste0(round(plot_data_new$percent[20], digits = 2),"%"," vs ",round(plot_data_new$percent[18], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[21])==FALSE) {
      seis<-paste0("OR=",plot_data_new$OR[21])
    } else {
      seis<-paste0(round(plot_data_new$percent[24], digits = 2),"%"," vs ",round(plot_data_new$percent[22], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[25])==FALSE) {
      siete<-paste0("OR=",plot_data_new$OR[25])
    } else {
      siete<-paste0(round(plot_data_new$percent[28], digits = 2),"%"," vs ",round(plot_data_new$percent[26], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[29])==FALSE) {
      ocho<-paste0("OR=",plot_data_new$OR[29])
    } else {
      ocho<-paste0(round(plot_data_new$percent[32], digits = 2),"%"," vs ",round(plot_data_new$percent[30], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[33])==FALSE) {
      nueve<-paste0("OR=",plot_data_new$OR[33])
    } else {
      nueve<-paste0(round(plot_data_new$percent[36], digits = 2),"%"," vs ",round(plot_data_new$percent[34], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[37])==FALSE) {
      diez<-paste0("OR=",plot_data_new$OR[37])
    } else {
      diez<-paste0(round(plot_data_new$percent[40], digits = 2),"%"," vs ",round(plot_data_new$percent[38], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[41])==FALSE) {
      once<-paste0("OR=",plot_data_new$OR[41])
    } else {
      once<-paste0(round(plot_data_new$percent[44], digits = 2),"%"," vs ",round(plot_data_new$percent[42], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[45])==FALSE) {
      doce<-paste0("OR=",plot_data_new$OR[45])
    } else {
      doce<-paste0(round(plot_data_new$percent[48], digits = 2),"%"," vs ",round(plot_data_new$percent[46], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[49])==FALSE) {
      trece<-paste0("OR=",plot_data_new$OR[49])
    } else {
      trece<-paste0(round(plot_data_new$percent[52], digits = 2),"%"," vs ",round(plot_data_new$percent[50], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[53])==FALSE) {
      catorce<-paste0("OR=",plot_data_new$OR[53])
    } else {
      catorce<-paste0(round(plot_data_new$percent[56], digits = 2),"%"," vs ",round(plot_data_new$percent[54], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[57])==FALSE) {
      quince<-paste0("OR=",plot_data_new$OR[57])
    } else {
      quince<-paste0(round(plot_data_new$percent[60], digits = 2),"%"," vs ",round(plot_data_new$percent[58], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[61])==FALSE) {
      dieciseis<-paste0("OR=",plot_data_new$OR[61])
    } else {
      dieciseis<-paste0(round(plot_data_new$percent[64], digits = 2),"%"," vs ",round(plot_data_new$percent[62], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[65])==FALSE) {
      diecisiete<-paste0("OR=",plot_data_new$OR[65])
    } else {
      diecisiete<-paste0(round(plot_data_new$percent[68], digits = 2),"%"," vs ",round(plot_data_new$percent[66], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[69])==FALSE) {
      dieciocho<-paste0("OR=",plot_data_new$OR[69])
    } else {
      dieciocho<-paste0(round(plot_data_new$percent[72], digits = 2),"%"," vs ",round(plot_data_new$percent[70], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[73])==FALSE) {
      diecinueve<-paste0("OR=",plot_data_new$OR[73])
    } else {
      diecinueve<-paste0(round(plot_data_new$percent[76], digits = 2),"%"," vs ",round(plot_data_new$percent[74], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[77])==FALSE) {
      veinte<-paste0("OR=",plot_data_new$OR[77])
    } else {
      veinte<-paste0(round(plot_data_new$percent[80], digits = 2),"%"," vs ",round(plot_data_new$percent[78], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[81])==FALSE) {
      veintiuno<-paste0("OR=",plot_data_new$OR[81])
    } else {
      veintiuno<-paste0(round(plot_data_new$percent[84], digits = 2),"%"," vs ",round(plot_data_new$percent[82], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[85])==FALSE) {
      veintiuno2<-paste0("OR=",plot_data_new$OR[85])
    } else {
      veintiuno2<-paste0(round(plot_data_new$percent[88], digits = 2),"%"," vs ",round(plot_data_new$percent[86], digits = 2),"%") }
    
    un<-print(paste0("p=",plot_data_new$fisher_p[1]))
    deux<-paste0("p=",plot_data_new$fisher_p[5])
    trois<-paste0("p=",plot_data_new$fisher_p[9])
    quatre<-paste0("p=",plot_data_new$fisher_p[13])
    cinq<-paste0("p=",plot_data_new$fisher_p[17])
    six<-paste0("p=",plot_data_new$fisher_p[21])
    sept<-paste0("p=",plot_data_new$fisher_p[25])
    huit<-paste0("p=",plot_data_new$fisher_p[29])
    neuf<-paste0("p=",plot_data_new$fisher_p[33])
    dix<-paste0("p=",plot_data_new$fisher_p[37])
    onze<-paste0("p=",plot_data_new$fisher_p[41])
    douze<-paste0("p=",plot_data_new$fisher_p[45])
    treize<-paste0("p=",plot_data_new$fisher_p[49])
    quatorze<-paste0("p=",plot_data_new$fisher_p[53])
    quinze<-paste0("p=",plot_data_new$fisher_p[57])
    seize<-paste0("p=",plot_data_new$fisher_p[61])
    dixsept<-paste0("p=",plot_data_new$fisher_p[65])
    dixhuit<-paste0("p=",plot_data_new$fisher_p[69])
    dixneuf<-paste0("p=",plot_data_new$fisher_p[73])
    vingt<-paste0("p=",plot_data_new$fisher_p[77])
    vingtetun<-paste0("p=",plot_data_new$fisher_p[81])
    vingtetun2<-paste0("p=",plot_data_new$fisher_p[85])
    
    eins<-paste(plot_data_new$drug_group[1], uno, un, sep = "\n")
    zwei<-paste(plot_data_new$drug_group[5], dos, deux, sep = "\n")
    drei<-paste(plot_data_new$drug_group[9], tres, trois, sep = "\n")
    vier<-paste(plot_data_new$drug_group[13], quatro, quatre, sep = "\n")
    fuenf<-paste(plot_data_new$drug_group[17], cinco, cinq, sep = "\n")
    sechs<-paste(plot_data_new$drug_group[21], seis, six, sep = "\n")
    sieben<-paste(plot_data_new$drug_group[25], siete, sept, sep = "\n")
    acht<-paste(plot_data_new$drug_group[29], ocho, huit, sep = "\n")
    neun<-paste(plot_data_new$drug_group[33], nueve, neuf, sep = "\n")
    zehn<-paste(plot_data_new$drug_group[37], diez, dix, sep = "\n")
    elf<-paste(plot_data_new$drug_group[41], once, onze, sep = "\n")
    zwoelf<-paste(plot_data_new$drug_group[45], doce, douze, sep = "\n")
    dreizehn<-paste(plot_data_new$drug_group[49], trece, treize, sep = "\n")
    vierzehn<-paste(plot_data_new$drug_group[53], catorce, quatorze, sep = "\n")
    fuenzehn<-paste(plot_data_new$drug_group[57], quince, quinze, sep = "\n")
    sechtszehn<-paste(plot_data_new$drug_group[61], dieciseis ,seize, sep = "\n")
    siebzehn<-paste(plot_data_new$drug_group[65], diecisiete, dixsept, sep = "\n")
    achtzehn<-paste(plot_data_new$drug_group[69], dieciocho, dixhuit, sep = "\n")
    neunzehn<-paste(plot_data_new$drug_group[73], diecinueve, dixneuf, sep = "\n")
    zwanzig<-paste(plot_data_new$drug_group[77], veinte, vingt, sep = "\n")
    einundzwanzig<-paste(plot_data_new$drug_group[81], veintiuno, vingtetun, sep = "\n")
    einundzwanzig2<-paste(plot_data_new$drug_group[85], veintiuno2, vingtetun2, sep = "\n")
    
    drug_names <- c(
      "ABL_signaling" = eins,
      "Apoptosis_regulation" = zwei,
      "Cell_cycle" = drei,
      "Chromatin_histone_acetylation" = vier,
      "Chromatin_histone_methylation" = fuenf,
      "Chromatin_other" = sechs,
      "Cytoskeleton" = sieben,
      "DNA_replication" = acht,
      "EGFR_signaling" = neun,
      "ERK_MAPK_signaling" = zehn,
      "Genome_integrity" = elf,
      "IGF1R_signaling" = zwoelf,
      "JNK_and_p38_signaling" = dreizehn,
      "Metabolism" = vierzehn,
      "Mitosis" = fuenzehn,
      "Other" = sechtszehn,
      "Other_kinases" = siebzehn,
      "p53_pathway" = achtzehn,
      "PI3K/MTOR_signaling" = neunzehn,
      "Protein_stab./deg." = zwanzig,
      "RTK_signaling" = einundzwanzig,
      "WNT_signaling" = einundzwanzig2
    )
    
    a=ggplot(data=plot_data_new[1:16,], 
             aes(fill=response, y=percent, x=wt_or_mut)) + 
      geom_bar(width=0.8, position="stack", stat="identity", color="gray") + 
      facet_grid(. ~ drug_group,labeller = as_labeller(drug_names))+
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid = element_line(size = 0.5, linetype = 'solid',
                                      colour = "grey90"),
            axis.ticks = element_line(color = "gray90"),
            axis.text = element_text(size=10, colour = "black"),
            strip.text.x = element_text(size = 8))+
      scale_x_discrete(name="")+
      scale_y_continuous(expand = expand_scale(mult = c(0.01, .03)))+
      scale_fill_manual("legend", values = c("resistant" = "indianred1", "sensitive" = "white"))+
      geom_text(data=plot_data_new[1:16,], aes(x = wt_or_mut, y = percent,
                                               label = paste0("n=",sum)), size=3, position = position_stack(vjust = 0.5))
    
    b=ggplot(data=plot_data_new[17:32,], 
             aes(fill=response, y=percent, x=wt_or_mut)) + 
      geom_bar(width=0.8, position="stack", stat="identity", color="gray") + 
      facet_grid(. ~ drug_group, labeller = as_labeller(drug_names))+
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid = element_line(size = 0.5, linetype = 'solid',
                                      colour = "grey90"),
            axis.ticks = element_line(color = "gray90"),
            axis.text = element_text(size=10, colour = "black"),
            strip.text.x = element_text(size = 8))+
      scale_x_discrete(name="")+
      scale_y_continuous(expand = expand_scale(mult = c(0.01, .03)))+
      scale_fill_manual("legend", values = c("resistant" = "indianred1", "sensitive" = "white"))+
      geom_text(data=plot_data_new[17:32,], aes(x = wt_or_mut, y = percent,
                                                label = paste0("n=",sum)), size=3, position = position_stack(vjust = 0.5))
    
    c=ggplot(data=plot_data_new[33:48,],
             aes(fill=response, y=percent, x=wt_or_mut)) + 
      geom_bar(width=0.8, position="stack", stat="identity", color="gray") + 
      facet_grid(. ~ drug_group, labeller = as_labeller(drug_names))+
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid = element_line(size = 0.5, linetype = 'solid',
                                      colour = "grey90"),
            axis.ticks = element_line(color = "gray90"),
            axis.text = element_text(size=10, colour = "black"),
            strip.text.x = element_text(size = 8))+
      scale_x_discrete(name="")+
      scale_y_continuous(expand = expand_scale(mult = c(0.01, .03)))+
      scale_fill_manual("legend", values = c("resistant" = "indianred1", "sensitive" = "white"))+
      geom_text(data=plot_data_new[33:48,], aes(x = wt_or_mut, y = percent,
                                                label = paste0("n=",sum)), size=3, position = position_stack(vjust = 0.5))
    d=ggplot(data=plot_data_new[49:68,], 
             aes(fill=response, y=percent, x=wt_or_mut)) + 
      geom_bar(width=0.9, position="stack", stat="identity", color="gray") + 
      facet_grid(. ~ drug_group, labeller = as_labeller(drug_names))+
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid = element_line(size = 0.5, linetype = 'solid',
                                      colour = "grey90"),
            axis.ticks = element_line(color = "gray90"),
            axis.text = element_text(size=10, colour = "black"),
            strip.text.x = element_text(size = 8))+
      scale_x_discrete(name="")+
      scale_y_continuous(expand = expand_scale(mult = c(0.01, .03)))+
      scale_fill_manual("legend", values = c("resistant" = "indianred1", "sensitive" = "white"))+
      geom_text(data=plot_data_new[49:68,], aes(x = wt_or_mut, y = percent,
                                                label = paste0("n=",sum)), size=3, position = position_stack(vjust = 0.5))
    
    e=ggplot(data=plot_data_new[69:nrow(plot_data_new),], 
             aes(fill=response, y=percent, x=wt_or_mut)) + 
      geom_bar(width=0.9, position="stack", stat="identity", color="gray") + 
      facet_grid(. ~ drug_group, labeller = as_labeller(drug_names))+
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid = element_line(size = 0.5, linetype = 'solid',
                                      colour = "grey90"),
            axis.ticks = element_line(color = "gray90"),
            axis.text = element_text(size=10, colour = "black"),
            strip.text.x = element_text(size = 8))+
      scale_x_discrete(name="")+
      scale_y_continuous(expand = expand_scale(mult = c(0.01, .03)))+
      scale_fill_manual("legend", values = c("resistant" = "indianred1", "sensitive" = "white"))+
      geom_text(data=plot_data_new[69:nrow(plot_data_new),], aes(x = wt_or_mut, y = percent,
                                                                 label = paste0("n=",sum)), size=3, position = position_stack(vjust = 0.5))
    
    R=grid.arrange(a,b,c,d,e, nrow=5)
    
  }
  
  #########################################################################################################################
  ################################# define variable plothight for forestplot ##############################################
  #########################################################################################################################
  
  plothight<- function(){
    plot_data_forest_list <- mydata_bigfile()
    forest<-data.frame(Reduce(rbind, plot_data_forest_list[2]),stringsAsFactors = FALSE)
    #forest$fisher_p<-str_replace(string = forest$fisher_p,pattern = "'",replacement = "")
    #forest$fisher_p <- as.numeric(forest$fisher_p)
    data=subset(forest, forest$fisher_p<0.05)
    data=data.frame(data,stringsAsFactors = FALSE)
    
    n=as.numeric(length(data$OR))
    if (n>7) {
      plothight<-(n*25)
    } else if (n>3) {
      plothight<-(n*35)
    } else {
      plothight<-100
    }
  }
  
  downloadplothight<- function(){
    x<-plothight()
    downloadplothight<-(x/50)
  }
  
  #########################################################################################################################
  ################################################# forestplot for download ###############################################
  #########################################################################################################################
  
  myforest<-function(){
    plot_data_forest_list <- mydata_bigfile()
    forest<-data.frame(Reduce(rbind, plot_data_forest_list[2]),stringsAsFactors = FALSE)
    #forest$fisher_p<-str_replace(string = forest$fisher_p,pattern = "'",replacement = "")
    #forest$fisher_p <- as.numeric(forest$fisher_p)
    data=subset(forest, forest$fisher_p<0.05)
    #data$OR<- as.numeric(data$OR)
    data <- data[!is.na(data$OR),]
    data=data.frame(data,stringsAsFactors = FALSE)
    #data$fisher_p <- as.numeric(data$fisher_p)
    #data$lower<- as.numeric(data$lower)
    #data$upper<- as.numeric(data$upper)
    
    print(data)
    str(data)
    
    f<-ggplot(data=data, aes(x=Screened_Compounds, y=OR, ymin=lower, ymax=upper)) +
      geom_pointrange() + 
      facet_grid(rows=vars(target_pathway),scales = "free_y",space="free_y")+
      geom_hline(yintercept=1, lty=2,color="red") +  # add a dotted line at x=1 after flip
      coord_flip(ylim = c(0, 7)) +  # flip coordinates (puts labels on y axis)
      #facet_grid(~target_pathway, strip.position="left")
      xlab("Compounds") + ylab("Mean (95% CI)") +
      scale_y_continuous(breaks = c(seq(0, 7,by=1)),label = c("0\nsensitive","1","2","3", "4", "5", "6","7\nresistant")) +
      #theme_bw()+  # use a white background+
      theme(strip.text.y = element_text(angle = 0))
    return(f)
  }
  
  #########################################################################################################################
  ########################################## forestplot for output at page ################################################
  #########################################################################################################################
  
  output$Forest <- renderPlot({
    plot_data_forest_list <- mydata_bigfile()
    forest<-data.frame(Reduce(rbind, plot_data_forest_list[2]),stringsAsFactors = FALSE)
    #forest$fisher_p<-str_replace(string = forest$fisher_p,pattern = "'",replacement = "")
    #forest$fisher_p <- as.numeric(forest$fisher_p)
    #forest$OR <- as.numeric(forest$OR)
    data_xy=subset(forest, forest$fisher_p<0.05)
    #data$OR<- as.numeric(data$OR)
    data_xy <- data_xy[!is.na(data_xy$OR),]
    #data=data.frame(data,stringsAsFactors = FALSE)
    #data$fisher_p <- as.numeric(data$fisher_p)
    #data$OR<- as.numeric(data$OR)
    #data$lower<- as.numeric(data$lower)
    #data$upper<- as.numeric(data$upper)
    print(data_xy)
    str(data_xy)
    #n=lenght(data)
    #plothight<-(n*50)
    
    #This warning occures 6 times only for output on page - can't find problem
    #Input to asJSON(keep_vec_names=TRUE) is a named vector. In a future version of jsonlite, this option will not be supported, 
    #and named vectors will be translated into arrays instead of objects. If you want JSON object output, please use a named list instead. See ?toJSON.
    
    f<-ggplot(data=data_xy, aes(x=Screened_Compounds, y=OR, ymin=lower, ymax=upper)) +
      geom_pointrange() + 
      facet_grid(rows=vars(target_pathway),scales = "free_y",space="free_y")+
      geom_hline(yintercept=1, lty=2,color="red") +  # add a dotted line at x=1 after flip
      coord_flip(ylim = c(0, 7)) +  # flip coordinates (puts labels on y axis)
      #facet_grid(~target_pathway, strip.position="left")
      xlab("Compounds") + ylab("OR (95% CI)") +
      scale_y_continuous(breaks = c(seq(0,7, by=1)),label = c("0\nsensitive","1","2","3", "4", "5", "6","7\nresistant")) +
      #theme_bw()+  # use a white background+
      theme(strip.text.y = element_text(angle = 0),
            axis.text = element_text(size = 13),
            axis.title = element_text(size = 15))
    
    #f2<-grid.arrange(f, bottom= textGrob("sensitive", y = 0, gp = gpar(fontsize = 12)))
    
    return(f)
  }, height = plothight)
  
  #########################################################################################################################
  ################################# plot Expanded Drug Response for output at page  ########################################
  #########################################################################################################################
  
  output$plotop <- renderPlot({
    plot_data_new_list <- mydata_bigfile()
    dg_or<-data.frame(Reduce(rbind, plot_data_new_list[1]))
    dg<-data.frame(Reduce(rbind, plot_data_new_list[3]))
    merge_data <- merge(dg, dg_or, by ="drug_group")
    plot_data_new<-merge_data
    
    if (is.na(plot_data_new$OR[1])==FALSE) {
      uno<-paste0("OR=",plot_data_new$OR[1])
    } else {
      uno<-paste0(round(plot_data_new$percent[4], digits = 2),"%"," vs ",round(plot_data_new$percent[2], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[5])==FALSE) {
      dos<-paste0("OR=",plot_data_new$OR[5])
    } else {
      dos<-paste0(round(plot_data_new$percent[8], digits = 2),"%"," vs ",round(plot_data_new$percent[6], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[9])==FALSE) {
      tres<-paste0("OR=",plot_data_new$OR[9])
    } else {
      tres<-paste0(round(plot_data_new$percent[12], digits = 2),"%"," vs ",round(plot_data_new$percent[10], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[13])==FALSE) {
      quatro<-paste0("OR=",plot_data_new$OR[13])
    } else {
      quatro<-paste0(round(plot_data_new$percent[16], digits = 2),"%"," vs ",round(plot_data_new$percent[14], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[17])==FALSE) {
      cinco<-paste0("OR=",plot_data_new$OR[17])
    } else {
      cinco<-paste0(round(plot_data_new$percent[20], digits = 2),"%"," vs ",round(plot_data_new$percent[18], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[21])==FALSE) {
      seis<-paste0("OR=",plot_data_new$OR[21])
    } else {
      seis<-paste0(round(plot_data_new$percent[24], digits = 2),"%"," vs ",round(plot_data_new$percent[22], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[25])==FALSE) {
      siete<-paste0("OR=",plot_data_new$OR[25])
    } else {
      siete<-paste0(round(plot_data_new$percent[28], digits = 2),"%"," vs ",round(plot_data_new$percent[26], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[29])==FALSE) {
      ocho<-paste0("OR=",plot_data_new$OR[29])
    } else {
      ocho<-paste0(round(plot_data_new$percent[32], digits = 2),"%"," vs ",round(plot_data_new$percent[30], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[33])==FALSE) {
      nueve<-paste0("OR=",plot_data_new$OR[33])
    } else {
      nueve<-paste0(round(plot_data_new$percent[36], digits = 2),"%"," vs ",round(plot_data_new$percent[34], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[37])==FALSE) {
      diez<-paste0("OR=",plot_data_new$OR[37])
    } else {
      diez<-paste0(round(plot_data_new$percent[40], digits = 2),"%"," vs ",round(plot_data_new$percent[38], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[41])==FALSE) {
      once<-paste0("OR=",plot_data_new$OR[41])
    } else {
      once<-paste0(round(plot_data_new$percent[44], digits = 2),"%"," vs ",round(plot_data_new$percent[42], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[45])==FALSE) {
      doce<-paste0("OR=",plot_data_new$OR[45])
    } else {
      doce<-paste0(round(plot_data_new$percent[48], digits = 2),"%"," vs ",round(plot_data_new$percent[46], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[49])==FALSE) {
      trece<-paste0("OR=",plot_data_new$OR[49])
    } else {
      trece<-paste0(round(plot_data_new$percent[52], digits = 2),"%"," vs ",round(plot_data_new$percent[50], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[53])==FALSE) {
      catorce<-paste0("OR=",plot_data_new$OR[53])
    } else {
      catorce<-paste0(round(plot_data_new$percent[56], digits = 2),"%"," vs ",round(plot_data_new$percent[54], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[57])==FALSE) {
      quince<-paste0("OR=",plot_data_new$OR[57])
    } else {
      quince<-paste0(round(plot_data_new$percent[60], digits = 2),"%"," vs ",round(plot_data_new$percent[58], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[61])==FALSE) {
      dieciseis<-paste0("OR=",plot_data_new$OR[61])
    } else {
      dieciseis<-paste0(round(plot_data_new$percent[64], digits = 2),"%"," vs ",round(plot_data_new$percent[62], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[65])==FALSE) {
      diecisiete<-paste0("OR=",plot_data_new$OR[65])
    } else {
      diecisiete<-paste0(round(plot_data_new$percent[68], digits = 2),"%"," vs ",round(plot_data_new$percent[66], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[69])==FALSE) {
      dieciocho<-paste0("OR=",plot_data_new$OR[69])
    } else {
      dieciocho<-paste0(round(plot_data_new$percent[72], digits = 2),"%"," vs ",round(plot_data_new$percent[70], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[73])==FALSE) {
      diecinueve<-paste0("OR=",plot_data_new$OR[73])
    } else {
      diecinueve<-paste0(round(plot_data_new$percent[76], digits = 2),"%"," vs ",round(plot_data_new$percent[74], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[77])==FALSE) {
      veinte<-paste0("OR=",plot_data_new$OR[77])
    } else {
      veinte<-paste0(round(plot_data_new$percent[80], digits = 2),"%"," vs ",round(plot_data_new$percent[78], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[81])==FALSE) {
      veintiuno<-paste0("OR=",plot_data_new$OR[81])
    } else {
      veintiuno<-paste0(round(plot_data_new$percent[84], digits = 2),"%"," vs ",round(plot_data_new$percent[82], digits = 2),"%") }
    
    if (is.na(plot_data_new$OR[85])==FALSE) {
      veintiuno2<-paste0("OR=",plot_data_new$OR[85])
    } else {
      veintiuno2<-paste0(round(plot_data_new$percent[88], digits = 2),"%"," vs ",round(plot_data_new$percent[86], digits = 2),"%") }
    
    un<-print(paste0("p=",plot_data_new$fisher_p[1]))
    deux<-paste0("p=",plot_data_new$fisher_p[5])
    trois<-paste0("p=",plot_data_new$fisher_p[9])
    quatre<-paste0("p=",plot_data_new$fisher_p[13])
    cinq<-paste0("p=",plot_data_new$fisher_p[17])
    six<-paste0("p=",plot_data_new$fisher_p[21])
    sept<-paste0("p=",plot_data_new$fisher_p[25])
    huit<-paste0("p=",plot_data_new$fisher_p[29])
    neuf<-paste0("p=",plot_data_new$fisher_p[33])
    dix<-paste0("p=",plot_data_new$fisher_p[37])
    onze<-paste0("p=",plot_data_new$fisher_p[41])
    douze<-paste0("p=",plot_data_new$fisher_p[45])
    treize<-paste0("p=",plot_data_new$fisher_p[49])
    quatorze<-paste0("p=",plot_data_new$fisher_p[53])
    quinze<-paste0("p=",plot_data_new$fisher_p[57])
    seize<-paste0("p=",plot_data_new$fisher_p[61])
    dixsept<-paste0("p=",plot_data_new$fisher_p[65])
    dixhuit<-paste0("p=",plot_data_new$fisher_p[69])
    dixneuf<-paste0("p=",plot_data_new$fisher_p[73])
    vingt<-paste0("p=",plot_data_new$fisher_p[77])
    vingtetun<-paste0("p=",plot_data_new$fisher_p[81])
    vingtetun2<-paste0("p=",plot_data_new$fisher_p[85])
    
    eins<-paste(plot_data_new$drug_group[1], uno, un, sep = "\n")
    zwei<-paste(plot_data_new$drug_group[5], dos, deux, sep = "\n")
    drei<-paste(plot_data_new$drug_group[9], tres, trois, sep = "\n")
    vier<-paste(plot_data_new$drug_group[13], quatro, quatre, sep = "\n")
    fuenf<-paste(plot_data_new$drug_group[17], cinco, cinq, sep = "\n")
    sechs<-paste(plot_data_new$drug_group[21], seis, six, sep = "\n")
    sieben<-paste(plot_data_new$drug_group[25], siete, sept, sep = "\n")
    acht<-paste(plot_data_new$drug_group[29], ocho, huit, sep = "\n")
    neun<-paste(plot_data_new$drug_group[33], nueve, neuf, sep = "\n")
    zehn<-paste(plot_data_new$drug_group[37], diez, dix, sep = "\n")
    elf<-paste(plot_data_new$drug_group[41], once, onze, sep = "\n")
    zwoelf<-paste(plot_data_new$drug_group[45], doce, douze, sep = "\n")
    dreizehn<-paste(plot_data_new$drug_group[49], trece, treize, sep = "\n")
    vierzehn<-paste(plot_data_new$drug_group[53], catorce, quatorze, sep = "\n")
    fuenzehn<-paste(plot_data_new$drug_group[57], quince, quinze, sep = "\n")
    sechtszehn<-paste(plot_data_new$drug_group[61], dieciseis ,seize, sep = "\n")
    siebzehn<-paste(plot_data_new$drug_group[65], diecisiete, dixsept, sep = "\n")
    achtzehn<-paste(plot_data_new$drug_group[69], dieciocho, dixhuit, sep = "\n")
    neunzehn<-paste(plot_data_new$drug_group[73], diecinueve, dixneuf, sep = "\n")
    zwanzig<-paste(plot_data_new$drug_group[77], veinte, vingt, sep = "\n")
    einundzwanzig<-paste(plot_data_new$drug_group[81], veintiuno, vingtetun, sep = "\n")
    einundzwanzig2<-paste(plot_data_new$drug_group[85], veintiuno2, vingtetun2, sep = "\n")
    
    drug_names <- c(
      "ABL_signaling" = eins,
      "Apoptosis_regulation" = zwei,
      "Cell_cycle" = drei,
      "Chromatin_histone_acetylation" = vier,
      "Chromatin_histone_methylation" = fuenf,
      "Chromatin_other" = sechs,
      "Cytoskeleton" = sieben,
      "DNA_replication" = acht,
      "EGFR_signaling" = neun,
      "ERK_MAPK_signaling" = zehn,
      "Genome_integrity" = elf,
      "IGF1R_signaling" = zwoelf,
      "JNK_and_p38_signaling" = dreizehn,
      "Metabolism" = vierzehn,
      "Mitosis" = fuenzehn,
      "Other" = sechtszehn,
      "Other_kinases" = siebzehn,
      "p53_pathway" = achtzehn,
      "PI3K/MTOR_signaling" = neunzehn,
      "Protein_stab./deg." = zwanzig,
      "RTK_signaling" = einundzwanzig,
      "WNT_signaling" = einundzwanzig2
    )
    
    a=ggplot(data=plot_data_new[1:16,], 
             aes(fill=response, y=percent, x=wt_or_mut)) + 
      geom_bar(width=0.8, position="stack", stat="identity", color="gray") + 
      facet_grid(. ~ drug_group,labeller = as_labeller(drug_names))+
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid = element_line(size = 0.5, linetype = 'solid',
                                      colour = "grey90"),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 11),
            axis.ticks = element_line(color = "gray90"),
            axis.text = element_text(size=13, colour = "black"),
            strip.text.x = element_text(size = 11),
            axis.title.y = element_text(size = 14))+
      scale_x_discrete(name="")+
      scale_y_continuous(expand = expand_scale(mult = c(0.01, .03)))+
      scale_fill_manual("legend", values = c("resistant" = "indianred1", "sensitive" = "white"))+
      geom_text(data=plot_data_new[1:16,], aes(x = wt_or_mut, y = percent,
                                               label = paste0("n=",sum)), size=4, position = position_stack(vjust = 0.5))
    
    b=ggplot(data=plot_data_new[17:32,], 
             aes(fill=response, y=percent, x=wt_or_mut)) + 
      geom_bar(width=0.8, position="stack", stat="identity", color="gray") + 
      facet_grid(. ~ drug_group, labeller = as_labeller(drug_names))+
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid = element_line(size = 0.5, linetype = 'solid',
                                      colour = "grey90"),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 11),
            axis.ticks = element_line(color = "gray90"),
            axis.text = element_text(size=13, colour = "black"),
            strip.text.x = element_text(size = 11),
            axis.title.y = element_text(size = 14))+
      scale_x_discrete(name="")+
      scale_y_continuous(expand = expand_scale(mult = c(0.01, .03)))+
      scale_fill_manual("legend", values = c("resistant" = "indianred1", "sensitive" = "white"))+
      geom_text(data=plot_data_new[17:32,], aes(x = wt_or_mut, y = percent,
                                                label = paste0("n=",sum)), size=4, position = position_stack(vjust = 0.5))
    
    c=ggplot(data=plot_data_new[33:48,],
             aes(fill=response, y=percent, x=wt_or_mut)) + 
      geom_bar(width=0.8, position="stack", stat="identity", color="gray") + 
      facet_grid(. ~ drug_group, labeller = as_labeller(drug_names))+
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid = element_line(size = 0.5, linetype = 'solid',
                                      colour = "grey90"),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 11),
            axis.ticks = element_line(color = "gray90"),
            axis.text = element_text(size=13, colour = "black"),
            strip.text.x = element_text(size = 11),
            axis.title.y = element_text(size = 14))+
      scale_x_discrete(name="")+
      scale_y_continuous(expand = expand_scale(mult = c(0.01, .03)))+
      scale_fill_manual("legend", values = c("resistant" = "indianred1", "sensitive" = "white"))+
      geom_text(data=plot_data_new[33:48,], aes(x = wt_or_mut, y = percent,
                                                label = paste0("n=",sum)), size=4, position = position_stack(vjust = 0.5))
    d=ggplot(data=plot_data_new[49:68,], 
             aes(fill=response, y=percent, x=wt_or_mut)) + 
      geom_bar(width=0.9, position="stack", stat="identity", color="gray") + 
      facet_grid(. ~ drug_group, labeller = as_labeller(drug_names))+
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid = element_line(size = 0.5, linetype = 'solid',
                                      colour = "grey90"),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 11),
            axis.ticks = element_line(color = "gray90"),
            axis.text = element_text(size=13, colour = "black"),
            strip.text.x = element_text(size = 11),
            axis.title.y = element_text(size = 14))+
      scale_x_discrete(name="")+
      scale_y_continuous(expand = expand_scale(mult = c(0.01, .03)))+
      scale_fill_manual("legend", values = c("resistant" = "indianred1", "sensitive" = "white"))+
      geom_text(data=plot_data_new[49:68,], aes(x = wt_or_mut, y = percent,
                                                label = paste0("n=",sum)), size=4, position = position_stack(vjust = 0.5))
    
    e=ggplot(data=plot_data_new[69:nrow(plot_data_new),], 
             aes(fill=response, y=percent, x=wt_or_mut)) + 
      geom_bar(width=0.9, position="stack", stat="identity", color="gray") + 
      facet_grid(. ~ drug_group, labeller = as_labeller(drug_names))+
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            panel.grid = element_line(size = 0.5, linetype = 'solid',
                                      colour = "grey90"),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 11),
            axis.ticks = element_line(color = "gray90"),
            axis.text = element_text(size=13, colour = "black"),
            strip.text.x = element_text(size = 11),
            axis.title.y = element_text(size = 14))+
      scale_x_discrete(name="")+
      scale_y_continuous(expand = expand_scale(mult = c(0.01, .03)))+
      scale_fill_manual("legend", values = c("resistant" = "indianred1", "sensitive" = "white"))+
      geom_text(data=plot_data_new[69:nrow(plot_data_new),], aes(x = wt_or_mut, y = percent,
                                                                 label = paste0("n=",sum)), size=4, position = position_stack(vjust = 0.5))
    
    R=grid.arrange(a,b,c,d,e, nrow=5)
  })
  
  #########################################################################################################################
  ################################################## downloadHandler ######################################################
  #########################################################################################################################
  
  output$groupsfiles <- downloadHandler(
    filename = function() {
      if (input$amountgenes=="One") {
        paste(input$oncgenes, "_groups.txt", sep = "")
      } else if (input$amountgenes=="Two") {
        paste(input$oncgenes, "_", input$secondgene, "_groups.txt", sep = "")
      } else
        paste(input$oncgenes, "_", input$second3gene, "_", input$third3gene, "_groups.txt", sep = "")
    },
    content = function(file) {
      write.table(groupsfile(), file, sep = "\t",
                  row.names = FALSE,quote = FALSE)
    }
  )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      if (input$amountgenes=="One") {
        paste(input$oncgenes, "_tables.xlsx", sep = "")
      } else if (input$amountgenes=="Two") {
        paste(input$oncgenes, "_", input$secondgene, "_tables.xlsx", sep = "")
      } else
        paste(input$oncgenes, "_", input$second3gene, "_", input$third3gene, "_tables.xlsx", sep = "")
    },
    content = function(file) {
      write_xlsx(mydata_bigfile(),path = file)
    }
  )
  
  output$downloadData_plot <- downloadHandler(
    filename = function() {
      if (input$amountgenes=="One") {
        paste(input$oncgenes, "_plot.png", sep = "")
      } else if (input$amountgenes=="Two") {
        paste(input$oncgenes, "_", input$secondgene, "_plot.png", sep = "")
      } else
        paste(input$oncgenes, "_", input$second3gene, "_", input$third3gene, "_plot.png", sep = "")
    },
    content = function(file) {
      ggsave(file,plot = myplot(),width = 21.0,height = 48, units = "cm",device = "png")
    }
  )
  
  output$downloadData_forest <- downloadHandler(
    filename = function() {
      if (input$amountgenes=="One") {
        paste(input$oncgenes, "_forestplot.png", sep = "")
      } else if (input$amountgenes=="Two") {
        paste(input$oncgenes, "_", input$secondgene, "_forestplot.png", sep = "")
      } else
        paste(input$oncgenes, "_", input$second3gene, "_", input$third3gene, "_forestplot.png", sep = "")
    },
    content = function(file) {
      ggsave(file,plot = myforest(),width = 21.0,height = downloadplothight(), units = "cm",device = "png",limitsize = FALSE)
    }
  )
  
  output$downloadData_boxplot <- downloadHandler(
    filename = function() {
      if (input$amountgenes=="One") {
        paste(input$oncgenes, "_boxplot.png", sep = "")
      } else if (input$amountgenes=="Two") {
        paste(input$oncgenes, "_", input$secondgene, "_boxplot.png", sep = "")
      } else
        paste(input$oncgenes, "_", input$second3gene, "_", input$third3gene, "_boxplot.png", sep = "")
    },
    content = function(file) {
      ggsave(file,plot = boxPlot(),width = 26.0,height = downloadplothight_boxplot(), units = "cm",device = "png",limitsize = FALSE)
    }
  )
  
  output$downloadData_wilcoxtable<-downloadHandler(
    filename = function() {
      if (input$amountgenes=="One") {
        paste(input$oncgenes, "_boxplot_data.csv", sep = "")
      } else if (input$amountgenes=="Two") {
        paste(input$oncgenes, "_", input$secondgene, "_boxplot_data.csv", sep = "")
      } else
        paste(input$oncgenes, "_", input$second3gene, "_", input$third3gene, "_boxplot_data.csv", sep = "")
    },
    content = function(file) {
      write.csv(Wilcoxtable(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui,server = server)

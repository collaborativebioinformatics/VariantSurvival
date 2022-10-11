source("global.R", local = TRUE)

ui <- fluidPage( 
  dashboardPage(
    dashboardHeader(title = "geneTarget"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Import Data", tabName = "DataImport", icon = icon("dna")),
        menuItem("Survival Analysis", tabName = "Survival_Analysis", icon = icon("dna"))
      )
    ),
    dashboardBody(use_waiter(),
                  
                  tabItems(
                    #### ------------------ Import data dashboard -------------------###############
                    tabItem(tabName = "DataImport",
                            #start fluidrow csv file load
                            fluidRow(
                              # start box right (contain input) ==> read csv file
                              box(status = "info", width = 4,
                                  solidHeader = TRUE,
                                  collapsible = TRUE,
                                  title = "Import File",
                                  br(),
                                  "Import the variant data vcf file here",
                                  br(),
                                  fileInput("input$vcf_file", label = "Upload vcf file"),
                                  
                                  #, shinyDirButton("dir", "Input directory", "Upload")

                              ), # end box right
                              box(status = "info", width = 4,
                                  solidHeader = TRUE,
                                  collapsible = TRUE,
                                  title = "Import File",
                                  br(),
                                  "Import the metadata file here",
                                  br(),
                                  fileInput("input$meta_file", label = "Upload metadata file"),
                              ),
                              box(status = "info", width = 4,
                                  solidHeader = TRUE,
                                  collapsible = TRUE,
                                  title = "Chose disease",
                                  br(),
                                  "Chose the disease from the list below",
                                  br(),
                                  selectizeInput("disease", "diseases", choices = c("Alzheimer's disease"=1,"Parkinson's disease"=2, selected = 1)),
                                  )),  # end fluid row 1
                              fluidRow(
                               box(width = 6,
                                   title = "Genes",
                                  h1("Based on literature The following genes are associated with the disease mechanism"),
                                  br(),
                                  DT::dataTableOutput("geneslist"),
                                  br(),
                                  textInput("targetGene", label = h3("Select your gene of interest"), value = "Gene Symbol")
                                  #selectizeInput('gene', 'Genes', selected =NULL, choices = levels(geneslist))
                                  ),
                               box(width = 6, height = 700,
                                   title = "Structural Variant in selected gene",
                                 h1("barplot count of SV in the selected gene, in each sample + color mark the samples with same group")
                                 #,plotOutput("barplot")
                               )
                              ),
                            fluidRow(
                              box(title = "Variants", width = 12, height = 600,
                                  br(),
                                DT::dataTableOutput("data_matrix_vcf") 
                              )
                              
                            ),
                              fluidRow(
                                box(actionButton("analyze", "Start Analysis"))
                              )
                    ), #import tabitem

                    #
                    tabItem(tabName = "Survival_Analysis", h2("Survival Analysis Results"),
                            box(width = 12,
                              "survival plot here", 
                              br(),
                              downloadButton("plot.png", "dowload plot"),
                              plotOutput("plot1", height = 500))
                              
                            )
                  )#tabitems

                  )#dashbody
    )#dashpage
)#fluidpage


server <- function(input, output) {

  ## read data vcf file input$vcf_file
  file_data <- reactive({
    file1 <- input$vcf_file
    if(!is.null(file1)){read.csv(file1$datapath)}
  })
  
  # 
  output$data_matrix_vcf <- DT::renderDataTable({
    #req(file_data())
    req(input$vcf_file)
    file_data()
  }) 
  
  #
  output$geneslist <- DT::renderDataTable({
    if(as.numeric(input$disease == 1)){
      geneslist  <-read_csv("disease_gene/AD/genes_list.txt")
    } #else if = the other diseases
  }) 
  #source("server/survival.R", local = TRUE)

}


shinyApp(ui = ui, server = server)

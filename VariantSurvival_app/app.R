source("global.R", local = TRUE)

ui <- fluidPage(
  dashboardPage(
    dashboardHeader(title = "VariantSurvival"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Import Data", tabName = "DataImport", icon = icon("dna")),
        menuItem("Survival Analysis", tabName = "Survival_Analysis", icon = icon("dna"))
      )
    ),
    dashboardBody(
                  
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
                                  selectizeInput("disease", "diseases", choices = c("Amyotrophic lateral sclerosis"=1,"Parkinson's disease"=2,"Alzheimer's disease"=3,
                                                                                    "Friedreich ataxia"=4, "Huntington's disease"=5, "Lewy body disease"=6, 
                                                                                    "Spinal muscular atrophy"=7, selected = 1)),
                                  )),  # end fluid row 1
                              fluidRow(
                               box(width = 6,
                                   title = "Genes",
                                   h2("Note : It is recommended to use the Illumina ExpansionHunter tool for the SV calling."),
                                  h1("Based on literature The following genes are associated with the disease mechanism"),
                                  br(),
                                  DT::dataTableOutput("geneslist"),
                                  br(),
                                  textInput("targetGene", label = h3("Select your gene of interest")#, value = "SETX"
                                            )
                                  
                                  #selectizeInput('gene', 'Genes', selected =NULL, choices = levels(geneslist))
                                  ),
                               box(width = 6, height = 800,
                                   title = "Structural Variant in selected gene",
                                   br(),
                                   br(),
                                   br(),
                                   br(),
                                plotOutput("barplot")
                               )
                              
                                #,box(DT::dataTableOutput("svtable"))
                    )), #import tabitem

                    #
                    tabItem(tabName = "Survival_Analysis", h2("Survival Analysis Results"),
                            box(title = "Survival Plot according to existing or not of the SVs",
                                width = 9, height = 600,
                                downloadButton(outputId = "plotsr", label = "Download the plot"),
                                plotOutput("plot4")
                                ),
                            box(width = 3, height = 300,
                                title = "Competing risks regression",
                                h4("Factor = Variant"),
                                DT::dataTableOutput("table")),
                            box(width = 3, height = 300,
                                title = "Competing risks regression", h4("Factor = Phenotype"),
                                DT::dataTableOutput("table2")),
                            box(width = 12, height = 600,
                              title = "Survival Plot according to SV count",
                              plotOutput("plot5")
                            )
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
  source("server/step1.R", local = TRUE)
  source("server/step2.R", local = TRUE)
  


}


shinyApp(ui = ui, server = server)

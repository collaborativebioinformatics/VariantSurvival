#' VariantSurvival
#'
#' @param vcffile path to the vcf file containing the Structural variant data
#' @param metadatafile path to the txt file containing the samples metadata
#'
#' @return
#' @export
#'

source("global.R") # install/load requirements
source("server.R") # main functions

VariantSurvival <- function(vcffile, metadatafile){
  # parse inputs
  vcf <- vcfR::read.vcfR(vcffile, verbose = FALSE)
  metadata <- readxl::read_excel(metadatafile)
  # remove empty extra lines
  metadata <- na.omit(metadata)
  gene_ids_table <- read.csv(file = 'ensembleTogenes.csv')
  rownames(gene_ids_table) <- gene_ids_table$ensembleID
  # create user interface layout
  ui <- fluidPage(
    dashboardPage(
      dashboardHeader(title = "VariantSurvival"),
      dashboardSidebar(
        selectizeInput(inputId ="disease_n",
                       label = "Diseases",
                       choices = c("Amyotrophic lateral sclerosis"=1,
                                   "Parkinson's disease"=2,
                                   "Alzheimer's disease"=3,
                                   "Friedreich ataxia"=4,
                                   "Huntington's disease"=5,
                                   "Lewy body disease"=6,
                                   "Spinal muscular atrophy"=7),
                       selected = FALSE
        ),
        selectInput(inputId = "target_gene",
                    label = "Gene of interest:",
                    choices = NULL,
                    selected = "SETX"
        ),
        selectInput(inputId = "time",
                    label = "Select the time factor:",
                    choices = colnames(metadata),
                    selected = FALSE
        ),
        selectInput(inputId = "event",
                    label = "Select the event factor:",
                    choices = colnames(metadata),
                    selected = FALSE
        ),
        selectInput(inputId = "phenotype",
                    label = "Select the phenotype factor:",
                    choices = colnames(metadata),
                    selected = FALSE
        )
      ),
      dashboardBody(
        fluidRow(
          box(width = 7, #height = 600,
              title = "TBC",
              br(), #linebreak
              plotOutput("sc_no_sv")
          ),
          box(width = 5, #height = 600,
              title = "Structural Variant in selected gene",
              br(),
              plotOutput("barplot")
          )
        ),
        fluidRow(
          tabBox(width = 6,# height = 600,
                 tabPanel(title = "Survival Plot",
                          h6("Phenotype 0 = Placebo ; 1= Treatment"),
                          plotOutput("plot2"),
                          downloadButton("download2plot", "Download as PNG")
                 ),
                 tabPanel(title = "Survival Plot according to SV count",
                          h6("starta 0 = Placebo ; 1= Treatment"),
                          plotOutput("plot3"),
                          downloadButton("download3plot", "Download as PNG")
                 ),
                 tabPanel( title = "Competing risks regression",
                           DT::dataTableOutput("table3"),
                           "HR = Hazard Ratio, CI = Confidence Interval")
          ),
          tabBox(width = 6, #height = 700,
                 tabPanel(title = "Kaplan-Meier plot",
                          plotOutput("plot1"),
                          downloadButton("download1plot", "Download as PNG")
                 ),
                 tabPanel(title = "x-year survival time time",
                          DT::dataTableOutput("table1"),
                          "CI = Confidence Interval"),
                 tabPanel(title = "median survival time",
                          DT::dataTableOutput("table2"),
                          "CI = Confidence Interval")
          )
        )
      )
    )
  )
  # Run the application
  shinyApp(ui = ui, server = server)
}


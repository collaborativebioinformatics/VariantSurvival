#' VariantSurvival
#'
#' @param vcffile path to the vcf file containing the Structural variant data
#' @param metadatafile path to the txt file containing the samples metadata
#'
#' @return
#' @export
#'

VariantSurvival <- function(vcffile, metadatafile){
  install_load_requirements()
  # parse inputs
  vcf <- vcfR::read.vcfR(vcffile, verbose = FALSE)
  metadata <- readxl::read_excel(metadatafile)
  # remove empty extra lines
  metadata <- na.omit(metadata)
  # create user interface layout


  ui <- bootstrapPage(
    navbarPage(theme = shinytheme("flatly"),
               collapsible = TRUE,
               HTML('<a style="text-decoration:none;
               cursor:default;
                    color:#FFFFFF;
                    " class="active" href="#">VariantSurvival</a>'),
               id="nav",
               windowTitle ="VariantSurvival",
               sidebarLayout(
                 sidebarPanel(
                   pickerInput(inputId ="disease_n",
                               label = "Diseases",
                               choices = colnames(disease_gene),
                               selected = FALSE
                               ),
                               span(shiny::tags$i(
                               h4("Based on literature the following genes are
                               associated with the disease mechanism")
                               ),
                               style="color:#045a8d"),
                   selectInput(inputId = "target_gene",
                               label = "Gene of interest:",
                               choices = NULL,
                               selected = FALSE
                               ),
                   selectInput("time",
                               label = "Select the time factor:",
                               choices = colnames(metadata)
                               ),
                   selectInput("phenotype",
                               label = "Select the study group / phenotype factor:",
                               choices = colnames(metadata)
                               ),
                   selectInput("event",
                               label = "Select the alive/dead factor:",
                               choices = colnames(metadata)
                               )
                   ),
                 mainPanel(
                   span(shiny::tags$i(h2("Structural Variants Distribution")),
                   shinycssloaders::withSpinner(plotOutput("histogram"))
                   )
                 ),
               sidebarLayout(
                 sidebarPanel(
                   span(shiny::tags$i(h3("1-year survival time")),
                        style="color:#045a8d"),
                   DT::dataTableOutput("table1"),
                   br(),
                   br(),
                   span(shiny::tags$i(h3("Median survival time")),
                        style="color:#045a8d"),
                   DT::dataTableOutput("table2"),
                   br(),
                   br(),
                   span(shiny::tags$i(h3("Cox regression table")),
                        style="color:#045a8d"),
                   DT::dataTableOutput("table3")
                   ),
                 mainPanel(
                   span(shiny::tags$i(h2("Kaplan–Meier")), 
                   shinycssloaders::withSpinner(plotOutput("plot_km"))
                   )
                 )
               )
    )

  server <- function(input, output, session) {
    gene_ids_table <- read.csv(file = 'ensembleTogenes.csv')
    rownames(gene_ids_table) <- gene_ids_table$ensembleID
    # get disease_n input and update the genes list accordingly
    observeEvent(input$disease_n,
                 {
                   genes_list = c(get_disease_gene_list(input$disease_n))
                   updateSelectizeInput(session,
                                        input = "target_gene",
                                        choices = genes_list,
                                        selected = NULL)
                 }
    )
    # Update genes drop-down after disease input is given
    reactive_gene_list <- reactive({get_disease_gene_list(input$disease_n)})
    disease_genes_names = gene_ids_table$GeneName
    sample_names = colnames(vcf@gt)[-1] # VCF genotype information

    # genes are repeated since a single gene can have more than one SV.
    genes_with_svs_in_sample <- apply(vcf@fix, 1, getGeneName, gene_ids_table)

    # count dataframe with patient_ids in rows and gene_ids in columns
    count_df <- CountSVsDf(length(disease_genes_names),
                           length(sample_names),
                           disease_genes_names,
                           sample_names,
                           genes_with_svs_in_sample,
                           vcf
    )

    reactive_no_NAs_metadata <- reactive({
      new_md <- RemoveNAs(metadata, input$time)
      new_df <- subset(count_df,
                       patient_ID %in% new_md$patient_ID)
      no_na_df <- merge(new_df,
                        #-- need to fix this so it's not hardcoded!!
                        new_md[c("patient_ID",
                                 "Phenotype",
                                 input$event,
                                 input$time)],
                        on = "patient_ID") %>%  rename(time = input$time,
                                                       event = input$event)
      no_na_df <- transform(no_na_df,
                            time = as.numeric(time),
                            event = as.numeric(event)
      )
      no_na_df %>% mutate(SV_binary = ifelse(new_df[input$target_gene]>0, 1, 0))
    })

    ## output - histogram ##
    output$histogram <- renderPlot(
      {
        new_df <- reactive_no_NAs_metadata()
        # get a df with counts
        svs_gene_input_df <- hist_df(new_df, input)
        ggplot(svs_gene_input_df,
               aes(SV_count_per_gene,
                   fill=Phenotype)) +
          geom_histogram(binwidth=1) +
          stat_bin(binwidth=1,
                   geom='text',
                   color='white',
                   aes(label=after_stat(count)),
                   position=position_stack(vjust = 0.5)) +
          xlab("Number of SVs in target gene") + ylab("Frequency")+
        theme(                                             # All font sizes
            text = element_text(size = 20),
            panel.background = element_rect(fill = "white",
                                            colour = "black"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()#change font size of legend title   
            
          )+
          scale_fill_manual(values=c("#8B1D4D", "#5275A6")) #  bins color
      }
    )

    output$plot_km <- renderPlot(
      {
        svs_gene_input_df <- reactive_no_NAs_metadata()
        # subset those patients with and without the SV
        without_sv <- svs_gene_input_df[svs_gene_input_df$SV_binary == 0,]
        with_sv <- svs_gene_input_df[svs_gene_input_df$SV_binary == 1,]
        # generate survival curve objects for each group
        sc_without <- survfit2(Surv(time, event)~Phenotype, data = with_sv)
        sv_with <- survfit2(Surv(time, event)~Phenotype, data = without_sv)
        # create a list and plot
        surv_fit_list <- list("with SV"=sv_with,
                              "without SV"=sc_without)
        ggsurvplot_combine(surv_fit_list,
                           data=svs_gene_input_df,
                           risk.table=TRUE,
                           conf.int = FALSE,
                           conf.int.style = "step",
                           risk.table.y.text = FALSE,
                           ggtheme = theme_light(),
                           legend.labs =
                             c(paste("with", input$target_gene, "- placebo"),
                               paste("with", input$target_gene, "- treatment"),
                               paste("without", input$target_gene, "- placebo"),
                               paste("without", input$target_gene, "- treatment")
                             ),
                           palette = c("royalblue4",
                                                   "steelblue",
                                                   "seagreen3",
                                                   "turquoise3")
        )
      })
  }

  # Run the application
  shinyApp(ui = ui, server = server)
}


############# Helper functions ###################
install_load_requirements<- function() {
  if (!require("shiny")) install.packages("shiny")
  if (!require("shinyWidgets")) install.packages("shinyWidgets")
  if (!require("shinythemes")) install.packages("shinythemes")
  if (!require("shinycssloaders")) install.packages("shinycssloaders")
  if (!require("shinydashboard")) install.packages("shinydashboard")
  if (!require("DT")) install.packages("DT")
  if (!require("dplyr")) install.packages("dplyr")
  if (!require("tidyverse")) install.packages("tidyverse")
  if (!require("vcfR")) install.packages("vcfR")
  if (!require("readr")) install.packages("readr")
  if (!require("ggsurvfit")) install.packages("ggsurvfit")
  if (!require("survival")) install.packages("survival")
  if (!require("ggplot2")) install.packages("ggplot2")
  if (!require("survminer")) install.packages("survminer")
  if (!require("lubridate")) install.packages("lubridate")
  if (!require("gtsummary")) install.packages("gtsummary")
  library(shiny)
  library(shinydashboard)
  library(shinythemes)
  library(shinyWidgets)
  library(shinycssloaders)
  library(DT)
  library(dplyr)
  library(tidyverse)
  library(vcfR)
  library(readr)
  library(readxl)
  library(ggsurvfit)
  library(survival)
  library(ggplot2)
  library(survminer)
  library(lubridate)
  library(gtsummary)
}


#' `get_disease_gene_list` parses an existing excel 
#' file containing all known target genes associated
#' with the selected disease.
#' @param input_disease: input disease selection
#' @return  genes_list
get_disease_gene_list <- function(input_disease) {
  disease_gene_select <- disease_gene[,input$disease_n]
  colnames(disease_gene_select)[1] = "gene"
  genes_list <- disease_gene_select[1]
  genes_list <- na.omit(genes_list)
  return(genes_list)
}

#' `getGeneName`
#'
#' @param info
#' @return
getGeneName <- function(info, geneIDS) {
  x <- str_extract(info['INFO'], "(?<=ensembl_gene_id=)[^;]+")
  return(geneIDS[x,]$GeneName)
}


#' `CountSVsDf`
#'
#' @param ncol
#' @param nrow
#' @param disease_genes_names
#' @param  sample_names
#' @param genes_with_svs_in_sample
#' @param vcf
#' @return
CountSVsDf <- function(ncol, nrow,
                       disease_genes_names,
                       sample_names,
                       genes_with_svs_in_sample,
                       vcf){
  sample_disease_gene_df <- data.frame(
    matrix(0,
           ncol = ncol,
           nrow = nrow
    )
  )
  colnames(sample_disease_gene_df) = disease_genes_names
  rownames(sample_disease_gene_df) = sample_names

  for (sv_idx in 1:length(genes_with_svs_in_sample)) {
    # if the sv is in a gene of interest
    sv_gene_name = genes_with_svs_in_sample[sv_idx]
    if(is.na(sv_gene_name))
      next;
    for(individual_idx in 1:length(sample_names)){
      indiv_id = sample_names[individual_idx]
      gt=vcf@gt[sv_idx, indiv_id]
      if(!is.na(gt))
      {
        sample_disease_gene_df[indiv_id, sv_gene_name] %+=% 1
      }
    }
  }
  sample_disease_gene_df <- rownames_to_column(sample_disease_gene_df,
                                               "patient_ID")
  return(sample_disease_gene_df)
}


#' `hist_df`
#'
#' @param df
#' @param input
#' @return
hist_df <- function(df,input){
  svs_gene_input_df <- df[c("patient_ID",
                            input$target_gene,
                            input$phenotype)
  ]
  colnames(svs_gene_input_df) <- c('patient_ID',
                                   'SV_count_per_gene',
                                   'Phenotype')

  svs_gene_input_df <- svs_gene_input_df %>%
    mutate(Phenotype = ifelse(Phenotype=="0",
                              "Placebo",
                              "Treatment" )
    )
  return(svs_gene_input_df)
}


#' `RemoveNAs` removes all rows with Na
#' entries in the time factor column.
#' @param df: metadata data frame
#' @param time_col: time factor input
#' @return
RemoveNAs <- function(df, time_col) {
  entries <- df[[time_col]]
  nas_entries <- is.na(entries)
  NAs_entries <- entries=="NA"
  if (sum(nas_entries) > 0){
    return(df[!nas_entries,])
  } else if (sum(NAs_entries) > 0){
    return(df[!NAs_entries,])
  } else {
    return(df)
  }
}


#' implementation of += operator
#' https://stackoverflow.com/questions/5738831/
`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))



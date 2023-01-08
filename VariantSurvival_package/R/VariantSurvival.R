#' VariantSurvival
#'
#' @param vcffile path to the vcf file containing the Structural variant data
#' @param metadatafile path to the txt file containing the samples metadata
#'
#' @return
#' @export
#'

function(vcffile, metadatafile){
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
                               choices = c("Amyotrophic lateral sclerosis"=1,
                                           "Parkinson's disease"=2,
                                           "Alzheimer's disease"=3,
                                           "Friedreich ataxia"=4,
                                           "Huntington's disease"=5,
                                           "Lewy body disease"=6,
                                           "Spinal muscular atrophy"=7),
                               selected = FALSE
                   ),
                   selectInput(inputId = "ids",
                               label = "Select the identifier factor:",
                               choices = colnames(metadata)
                   ),
                   span(shiny::tags$i(
                     h6("Based on literature the following genes are
                      associated with the disease mechanism")
                   ),
                   style="color:#045a8d"),
                   selectInput(inputId = "target_gene",
                               label = "Gene of interest:",
                               choices = NULL,
                               selected = FALSE
                   ),
                   selectInput(inputId = "time",
                               label = "Select the time factor:",
                               choices = colnames(metadata)
                   ),
                   span(shiny::tags$i(
                     h6("The following selections must refer to binary columns")
                   ),
                   style="color:#045a8d"),
                   selectInput(inputId = "group",
                               label = "Select the clinical trial groups factor:",
                               choices = colnames(metadata)
                   ),
                   selectInput(inputId = "event",
                               label = "Select the alive/dead factor:",
                               choices = colnames(metadata)
                   )
                 ),
                 mainPanel(
                   shinycssloaders::withSpinner(plotOutput("histogram"))
                 )
               ),
               sidebarLayout(
                 sidebarPanel(
                   span(shiny::tags$i(h6("add text here")),
                        style="color:#045a8d"),
                   DT::dataTableOutput("table1"),
                   br(),
                   br(),
                   span(shiny::tags$i( h6("add text here")),
                        style="color:#045a8d"),
                   DT::dataTableOutput("table2"),
                   br(),
                   br(),
                   span(shiny::tags$i(h6("add text here")),
                        style="color:#045a8d"),
                   DT::dataTableOutput("table3")
                 ),
                 mainPanel(
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

    reactive_no_NAs_metadata <- reactive({
      disease_genes_names = gene_ids_table$GeneName
      sample_names = colnames(vcf@gt)[-1] # ignore VCF genotype information

      # genes are repeated since a single gene can have more than one SV.
      genes_with_svs_in_sample <- apply(vcf@fix,
                                        1,
                                        getGeneName,
                                        gene_ids_table)
      # count dataframe with patient_ids in rows and gene_ids in columns
      count_df <- CountSVsDf(length(disease_genes_names),
                             length(sample_names),
                             disease_genes_names,
                             sample_names,
                             genes_with_svs_in_sample,
                             vcf,
                             input$ids)
      new_md <- (RemoveNAs(metadata, input$time)
                 %>% rename(ids = input$ids,
                            trial_group = input$group,
                            time = input$time,
                            event = input$event)
                 )
      new_df <- (subset(count_df, ids %in% new_md$ids)
                 %>% rename(SV_count_per_gene = input$target_gene)
                 )
      no_na_df <- merge(new_df, new_md[c("ids",
                                         "trial_group",
                                         "event",
                                         "time")], on = "ids")
      no_na_df <- transform(no_na_df,
                            time = as.numeric(time),
                            event = as.numeric(event),
                            SV_count_per_gene = as.numeric(SV_count_per_gene)
      )
      no_na_df <- (no_na_df
                   %>% mutate(SV_binary = ifelse(new_df$SV_count_per_gene>0, 1, 0))
                   )
      })

    ## output - histogram ##
    output$histogram <- renderPlot(
      {
        new_df <- reactive_no_NAs_metadata()
        # get a df with counts
        svs_gene_input_df <- hist_df(new_df, input)
        ggplot(svs_gene_input_df,
               aes(SV_count_per_gene,
                   fill=trial_group)) +
          geom_histogram(binwidth=1) +
          stat_bin(binwidth=1,
                   geom='text',
                   color='white',
                   aes(label=after_stat(count)),
                   position=position_stack(vjust = 0.5)) +
          xlab("Number of SVs in target gene") + ylab("Frequency")
      }
    )

    output$plot_km <- renderPlot({
      svs_gene_input_df <- reactive_no_NAs_metadata()
      survfit2(Surv(time, event)~ trial_group + SV_binary,
               data = svs_gene_input_df) %>%
        ggsurvfit() +
        labs(
          x = "years",
          y = "Overall survival probability"
        ) +
        theme(legend.position = c(0.85, 0.85))+
        add_confidence_interval() +
        add_risktable()
    })


    # output$plot_km <- renderPlot(
    #   {
    #     svs_gene_input_df <- reactive_no_NAs_metadata()
    #     # subset those patients with and without the SV
    #     without_sv <- svs_gene_input_df[svs_gene_input_df$SV_binary == 0,]
    #     with_sv <- svs_gene_input_df[svs_gene_input_df$SV_binary == 1,]
    #     # generate survival curve objects for each group
    #     sc_without <- survfit2(Surv(time, event)~Phenotype, data = with_sv)
    #     sv_with <- survfit2(Surv(time, event)~Phenotype, data = without_sv)
    #     # create a list and plot
    #     surv_fit_list <- list("with SV"=sv_with,
    #                           "without SV"=sc_without)
    #     ggsurvplot_combine(surv_fit_list,
    #                        data=svs_gene_input_df,
    #                        risk.table=TRUE,
    #                        conf.int = FALSE,
    #                        conf.int.style = "step",
    #                        risk.table.y.text = FALSE,
    #                        ggtheme = theme_light(),
    #                        legend.labs =
    #                          c(paste("with", input$target_gene, "- placebo"),
    #                            paste("with", input$target_gene, "- treatment"),
    #                            paste("without", input$target_gene, "- placebo"),
    #                            paste("without", input$target_gene, "- treatment")
    #                          ),
    #                        palette = c("royalblue4",
    #                                                "steelblue",
    #                                                "seagreen3",
    #                                                "turquoise3")
    #     )
    #   })
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


#' `get_disease_gene_list` parses an existing .txt
#' file containing all known target genes associated
#' with the selected disease.
#' @param input_disease: input disease selection
#' @return  genes_list
get_disease_gene_list <- function(input_disease) {
  if(as.numeric(input_disease == 1)){
    genes_list  <-read_csv("disease_gene/ALS.txt")
  }
  else if(as.numeric(input_disease == 2)){
    genes_list  <-read_csv("disease_gene/PD.txt")
  }
  else if(as.numeric(input_disease == 3)) {
    genes_list  <-read_csv("disease_gene/AD.txt")
  }
  else if(as.numeric(input_disease == 4)) {
    genes_list  <-read_csv("disease_gene/FD.txt")
  }
  else {
    genes_list  <-read_csv("disease_gene/DLB.txt")
  }
  return(genes_list)
}

#' `getGeneName`
#' @param info
#' @return
getGeneName <- function(info, geneIDS) {
  x <- str_extract(info['INFO'], "(?<=ensembl_gene_id=)[^;]+")
  return(geneIDS[x,]$GeneName)
}


#' `CountSVsDf`
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
                       vcf,
                       ids_col){
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
  sample_disease_gene_df <- (rownames_to_column(sample_disease_gene_df,ids_col)
                             %>%  rename(ids = ids_col))
  return(sample_disease_gene_df)
}


#' `hist_df`
#' @param df
#' @param input
#' @return
hist_df <- function(df,input){
  svs_gene_input_df <- df[c("ids",
                            "trial_group",
                            "SV_count_per_gene")]
  svs_gene_input_df <- (svs_gene_input_df
                        %>% mutate(trial_group = ifelse(trial_group=="0",
                                                        "Placebo",
                                                        "Treatment" )
                        )
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



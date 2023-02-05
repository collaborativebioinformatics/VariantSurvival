#' VariantSurvival
#'
#' @param vcffile path to the vcf file containing the Structural variant data
#' @param metadatafile path to the txt file containing the samples metadata
#' @param demo true or false
#'
#' @return
#' @export
#'

VariantSurvival <- function(vcffile, metadatafile,demo=FALSE){
  install_load_requirements()
  #demo or input files
  if (demo==TRUE){
    vcffile_demo <- "merged.filtered.vcf"
    vcf <- vcfR::read.vcfR(vcffile_demo, verbose = FALSE)
    metadata_demo <-"metadata.xlsx"
    metadata <- readxl::read_excel(metadata_demo)
    metadata <- na.omit(metadata)
  } else if (demo==FALSE){
    vcf <- vcfR::read.vcfR(vcffile, verbose = FALSE)
    metadata <- readxl::read_excel(metadatafile)
    metadata <- na.omit(metadata)

  }
  disease_gene <- read_excel("disease_gene.xlsx")
  days_year <- 365.25

  ui <- bootstrapPage(
    navbarPage(theme = shinytheme("flatly"),
               collapsible = TRUE,
               HTML('<a style="text-decoration:none;
               cursor:default;
                    color:#FFFFFF;
                    " class="active" href="#">VariantSurvival</a>'),
               id="nav",
               windowTitle ="VariantSurvival",
               tabPanel("Select Target Gene",
               sidebarLayout(
                 sidebarPanel(
                   pickerInput(inputId ="disease_n",
                               label = "Select the disease of interest:",
                               choices = c(colnames(disease_gene), "N/A"),
                               selected = "N/A"
                               ),
                   selectInput(inputId = "ids",
                               label = "Select the participant ID:",
                               choices = c(colnames(metadata), "N/A"),
                               selected = "N/A"
                               ),
                   selectInput(inputId = "time",
                               label = "Select the time factor:",
                               choices = c(colnames(metadata), "N/A"),
                               selected = "N/A"
                               ),
                   radioButtons(inputId = "time_unit",
                                label = "Time factor unit:",
                                choices =  c("years", "days"),
                                inline = TRUE
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
                   span(shiny::tags$i(
                     h6("The following selections must refer to binary factors")),
                     style="color:#045a8d"),
                   selectInput(inputId = "group",
                               label = "Select the clinical trial groups factor:",
                               choices =c(colnames(metadata), "N/A"),
                               selected = "N/A"
                               ),
                   selectInput(inputId = "event",
                               label = "Select the alive/dead factor:",
                               choices = c(colnames(metadata), "N/A"),
                               selected = "N/A"
                               )
                   ),
                 mainPanel(
                   tabBox(
                     selected = "Summary",
                     tabPanel("Summary",
                              span(DT::dataTableOutput("summ_table"))
                              ),
                     tabPanel("Histogram",
                              span(shiny::tags$i(
                                h2("Structural Variants Distribution"))),
                              shinycssloaders::withSpinner(
                                plotOutput(outputId = "histogram"))
                     ),
                     tabPanel("Table",
                              selectInput(inputId = "table_cols",
                                          label = "Include columns (optional)",
                                          choices = NULL,
                                          selected = FALSE,
                                          multiple = TRUE),
                              span(DT::dataTableOutput("table"))
                              )
                     )
                   )
                 )
               ),
         tabPanel("Survival Analysis",
                  mainPanel(
                    span(shiny::tags$i(
                      h2("Kaplan–Meier")),
                      shinyjs::useShinyjs(),
                      checkboxGroupInput("km_feat",
                                         "Modify plot",
                                         choices = c("confidence interval" = "conf_itv",
                                                     "risk table" = "risk_table",
                                                     "y grid line" = "grid_line")),
                      checkboxInput("all_n_svs",
                                    "Include all counts",
                                    value = TRUE),
                      selectInput(inputId = "n_svs_min",
                                  label = "min",
                                  choices = NULL,
                                  selected = FALSE
                                  ),
                      selectInput(inputId = "n_svs_max",
                                  label = "max",
                                  choices = NULL,
                                  selected = FALSE
                      ),
                      shinycssloaders::withSpinner(
                        plotOutput(outputId = "plot_km",width = "100%")
                        )
                    )
                  ) 
                  ),
         tabPanel("Cox regression",
                  span(shiny::tags$i(
                    h3("Cox regression table")),
                    style="color:#045a8d"),
                  selectizeInput(inputId = "sel_cov",
                                 label = "Select binary covariates",
                                 # SV_bin is added by us, 0/1 without/with SV
                                 choices = NULL,
                                 selected = FALSE,
                                 multiple = TRUE
                  ),
                  selectInput(inputId = "sel_strata",
                              label = "Select strata covariate (optional)",
                              choices = NULL),
                  DT::dataTableOutput("table3")
                  )
         )
    )
  
  
  server <- function(input, output, session) {
    gene_ids_table <- read.csv(file = 'ensembleTogenes.csv')
    rownames(gene_ids_table) <- gene_ids_table$ensembleID
    # get disease_n input and update the genes list accordingly
    observeEvent(input$disease_n, {
                   if(input$disease_n!= "N/A"){
                     genes_list <- c(get_disease_gene_list(disease_gene,input$disease_n))
                     disease_genes_names <- gene_ids_table$GeneName
                     updateSelectizeInput(session,
                                          input = "target_gene",
                                          choices = genes_list,
                                          selected = NULL)
                     }
      })
  
    observe({ if(input$disease_n != "N/A" 
                 & input$ids != "N/A" 
                 &input$time != "N/A"){
      genes_with_svs_in_sample <- apply(vcf@fix, 1, getGeneName,gene_ids_table)
      sample_names <- colnames(vcf@gt)[-1]
      disease_genes_names <- gene_ids_table$GeneName
      output$summ_table <- DT::renderDataTable({
        # count dataframe with patient_ids in rows and gene_ids in columns
        count_df <- CountSVsDf(length(disease_genes_names),
                               length(sample_names),
                               disease_genes_names,
                               sample_names,
                               genes_with_svs_in_sample,
                               vcf,
                               input$ids)
        new_md <- RemoveNAs(metadata, input$time) %>% rename(ids = input$ids)
        # what happen if the input$target_gene is not in the vcf file?
        new_df <- subset(count_df, ids %in% new_md$ids)
        count <- as.data.frame(apply(new_df[,-1], 2, function(c)sum(c!=0)))
        colnames(count) <- "Count of patients with SVs"
        count <- count %>% arrange(desc(count))
        count
        })
      }
    })
    
    reactive_no_NAs_metadata <- reactive({
      if(checkInput(input)){
        disease_genes_names <- gene_ids_table$GeneName
        sample_names <- colnames(vcf@gt)[-1] # ignore VCF genotype information
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
                              trial_group_bin = input$group,
                              time = input$time,
                              event = input$event)
        )
        # what happen if the input$target_gene is not in the vcf file?
        new_df <- (subset(count_df, ids %in% new_md$ids)
                   %>% rename(SV_count_per_gene = input$target_gene)
        )
        no_na_df <- merge(new_df, new_md, by = "ids")
        no_na_df <- transform(no_na_df,
                              time = as.numeric(time),
                              event = as.numeric(event),
                              SV_count_per_gene = as.numeric(SV_count_per_gene)
        )
        no_na_df <- (no_na_df
                     %>% mutate(SV_bin = ifelse(new_df$SV_count_per_gene>0, 1, 0)
                     )
                     %>% mutate(trial_group = ifelse(trial_group_bin=="0",
                                                     "Placebo",
                                                     "Treatment")
                     )
                     %>% mutate(SV = ifelse(SV_bin=="0", "with", "without")
                     )
        )
        if(input$time_unit == 'years'){
          no_na_df$time_days <- floor(no_na_df[["time"]] * days_year)
        } else{
          no_na_df$time_years <- no_na_df[["time"]] / days_year
        }
        no_na_df
        }
      }
      )

    ## output - histogram ##
    output$histogram <- renderPlot(
      {
        if(checkInput(input)){
          new_df <- reactive_no_NAs_metadata()
          # get a dataframe with counts
          svs_gene_input_df <-  new_df[c("ids", "trial_group", "SV_count_per_gene")]
          legend_title <- "Group"
          ggplot(svs_gene_input_df, aes(SV_count_per_gene, fill=trial_group)) +
            geom_histogram(binwidth=1) +
            stat_bin(binwidth=1,
                     geom='text',
                     color='white',
                     aes(label=after_stat(count)),
                     position=position_stack(vjust = 0.5)) +
            xlab("Number of SVs in target gene") + ylab("Frequency") +
            theme(text = element_text(size = 20),
                  panel.background = element_rect(fill = "white",
                                                  colour = "white"),
                  axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position="top"
            ) +
            scale_fill_manual(legend_title, values=c("#8B1D4D", "#5275A6"))
          }
        },
      height = 500,
      width = 700
      )
    
    observe({
             if(checkInput(input) & !is.null(input$target_gene)){
               svs_gene_input_df <- reactive_no_NAs_metadata()
               # renaming back to original column names, so it's not confusing for the user.
               # maybe not the best approach
               rename_cols <- c("ids", "trial_group_bin", "time", "event")
               rename_cols_with <- c(input$ids, input$group, input$time, input$event)
               names(svs_gene_input_df)[names(svs_gene_input_df) %in% rename_cols]<- rename_cols_with
               coln_list <- colnames(svs_gene_input_df)
               coln_list <- coln_list[!(coln_list %in% c(input$ids,"SV_count_per_gene"))]
               updateSelectizeInput(session,
                                    input = "table_cols",
                                    choices = coln_list,
                                    selected = NULL)
             }
      })
    
    observe({
      if(checkInput(input)){
        svs_gene_input_df <- reactive_no_NAs_metadata()
        # renaming back to original column names, so it's not confusing for the user.
        # maybe not the best approach
        rename_cols <- c("ids", "trial_group_bin", "time", "event")
        rename_cols_with <- c(input$ids, input$group, input$time, input$event)
        names(svs_gene_input_df)[names(svs_gene_input_df) %in% rename_cols]<- rename_cols_with
        if(any(!is.na(input$table_cols))){
          output$table <- DT::renderDataTable({
            col_names_list <- c(c(input$ids, "SV_count_per_gene"),
                                input$table_cols)
            svs_gene_input_df <- as_tibble(svs_gene_input_df[col_names_list])
            svs_gene_input_df
            })
          }
        else{
          output$table <- DT::renderDataTable({
          svs_gene_input_df <- as_tibble(svs_gene_input_df[c(input$ids, "SV_count_per_gene")])
          colnames(svs_gene_input_df) <- c("ID", "SV count per gene")
          svs_gene_input_df
          })
          }
        }
      })
        

    observe({
      if(checkInput(input)){
        risk_table = FALSE
        conf_itv = FALSE
        grid_line = "white"
          if(!is.null(input$km_feat)){
            if("conf_itv" %in% input$km_feat){
              conf_itv = TRUE}
            if("risk_table" %in% input$km_feat){
              risk_table = TRUE}
            if("grid_line" %in% input$km_feat){
              grid_line = "grey"}
          }
        svs_gene_input_df <- reactive_no_NAs_metadata()
        if (input$all_n_svs == FALSE 
            & input$n_svs_min!= "" 
            & input$n_svs_max!= ""){
          n_svs_min = as.numeric(input$n_svs_min)
          n_svs_max = as.numeric(input$n_svs_max)
          svs_gene_input_df$SV_count_per_gene
          sub <- (svs_gene_input_df$SV_count_per_gene >= n_svs_min 
                  & svs_gene_input_df$SV_count_per_gene <= n_svs_max)
          svs_gene_input_df <- svs_gene_input_df[sub,]
        }
        n_sv_groups <- unique(svs_gene_input_df[["SV_bin"]])
        n_groups <- unique(svs_gene_input_df[["trial_group_bin"]])
        # can it happen that there are more than 1 patient groups?
        if(length(n_sv_groups) == 2){
          # generate survival curve objects for each group
          sc_without <- survfit2(Surv(time, event)~trial_group_bin,
                                 data = svs_gene_input_df, subset=(SV_bin==0))
          sv_with <- survfit2(Surv(time, event)~trial_group_bin,
                              data = svs_gene_input_df, subset=(SV_bin==1))
          # create a list and plot
          surv_fit_list <- list("with SV"=sv_with, "without SV"=sc_without)
          lables <- c("with variant - placebo",
                     "with variant - treatment",
                     "without variant - placebo",
                     "without variant - treatment")
          cols <- c("Violetred4", "steelblue",
                   "Violetred2", "turquoise3")
        } 
        else if (length(n_sv_groups) == 1){
          if (n_sv_groups == 0){
            type <- "without"
            cols <- c("Violetred4", "steelblue")
          }
          else if (n_sv_groups == 1){
            type <- "with"
            cols <- c("Violetred2", "turquoise3")
          }
          sc <- survfit2(Surv(time, event)~trial_group_bin, 
                         data = svs_gene_input_df)
          # create a list and plot
          temp <- paste(type , " SV") 
          surv_fit_list <- list(temp = sc)
          if (length(n_groups) == 2){
            lables = c(paste(type, " variant - placebo"),
                       paste(type, " variant - treatment"))
          }
          else if (length(n_groups) < 2){
            if (n_groups == 0){
              lables = c(paste(type, " variant - placebo"))
              # for consistency with the coloring of placebo group
              cols <- cols[1] 
            }
            else if (n_groups == 1){
              lables = c(paste(type, " variant - treatment"))
              # for consistency with the coloring of treatment group
              cols <- cols[2] 
            }
          }
        }
        output$plot_km <- renderPlot({
          ggsurvplot_combine(surv_fit_list,
                             data=svs_gene_input_df,
                             risk.table = risk_table,
                             conf.int = conf_itv,
                             conf.int.style = "step",
                             risk.table.y.text = FALSE,
                             xlab = toTitleCase(input$time_unit),
                             ggtheme = theme(
                               text = element_text(size = 20),
                               panel.background = element_rect(fill = "white",
                                                               colour = "white"),
                               axis.line = element_line(colour = "black"),
                               panel.grid.major.y = element_line(colour=grid_line),
                               panel.grid.minor.y = element_line(colour=grid_line),
                               legend.position = c(0.2, 0.5)
                             ),
                             legend.title = "Group",
                             legend.labs =lables,
                             palette = cols
                             )},
          height = 600,
          width = 1000)
        }
        else{
          output$plot_km <- renderText({ "Missing input data"})
        }
      })
    
    inputs_react <- reactive({list(input$event,
                                   input$time,
                                   input$ids,
                                   input$target_gene
                                   )})
    # hide the event, time and ids covariates from the covariates  and strata
    # drop-down
    observeEvent(inputs_react(),
                 {if(checkInput(input)){
                   svs_gene_input_df <- reactive_no_NAs_metadata()
                   svs_levels <- unique(svs_gene_input_df["SV_count_per_gene"])$SV_count_per_gene
                   cov_list <- c(get_cov_list(metadata, input), "SV_bin")
                   # update the covariate drop down
                   updateSelectizeInput(session,
                                        input = "sel_cov",
                                        choices = cov_list,
                                        options = list(create = TRUE))
                   # update the strata drop down, this field is optional
                   updateSelectizeInput(session,
                                        input = "sel_strata",
                                        choices = c(cov_list, c("SV_bin","N/A")),
                                        selected = "N/A")
                   }
                }
              )
    
    observeEvent(input$all_n_svs, {
      if (checkInput(input)){
        if (input$all_n_svs == FALSE){
          shinyjs::enable("n_svs_min")
          shinyjs::enable("n_svs_max")
          svs_gene_input_df <- reactive_no_NAs_metadata()
          svs_levels <- unique(svs_gene_input_df["SV_count_per_gene"])$SV_count_per_gene
          updateSelectizeInput(session,
                               input = "n_svs_min",
                               choices = sort(svs_levels))
        }
        else {
          shinyjs::disable("n_svs_min")
          shinyjs::disable("n_svs_max")
        }
      }
      }
      )
    
    observeEvent(input$n_svs_min, {
      if (checkInput(input)){
        if (input$all_n_svs == FALSE & !is.null(input$n_svs_min) & input$n_svs_min!=""){
          svs_gene_input_df <- reactive_no_NAs_metadata()
          svs_levels <- unique(svs_gene_input_df["SV_count_per_gene"])$SV_count_per_gene
          svs_levels <- svs_levels[svs_levels >= input$n_svs_min]
          updateSelectizeInput(session,
                               input = "n_svs_max",
                               choices = sort(svs_levels))
        }
      }
      }
      )
    
    
    #regression table
    output$table3 <- DT::renderDataTable({
      if(checkInput(input) & any(!is.na(input$sel_cov))){
        input_cov_cox <- input$sel_cov
        # mapping original column names to the new ones
        if(input$event %in% input_cov_cox){
          input_cov_cox <- str_replace(input_cov_cox,
                                       input$event,
                                       "event")}
        if(input$group %in% input_cov_cox){
          input_cov_cox <- str_replace(input_cov_cox,
                                       input$group,
                                       "trial_group_bin")
          }
        svs_gene_input_df <- reactive_no_NAs_metadata()
        formulaString <- paste("Surv(time, event) ~", 
                               paste(input_cov_cox, collapse="+"))
        x3 <- (coxph(as.formula(formulaString), data=svs_gene_input_df) 
               %>% tbl_regression(exp = TRUE))
        t3 <-as_tibble(x3)
        t3
      }
    })
  }

  # Run the application
  shinyApp(ui = ui, server = server)
}


############# Helper functions ###################
install_load_requirements<- function() {
  if (!require("shiny")) install.packages("shiny")
  if (!require("shiny")) install.packages("shinyjs")
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
  if (!require("tools")) install.packages("tools")
  library(shiny)
  library(shinyjs)
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
get_disease_gene_list <- function(disease_gene, input_disease) {
  genes_list <-( na.omit(disease_gene[,input_disease])
                 %>% rename(genes = input_disease))
  return(genes_list$genes)
}


#' `checkInput`
#' @param info
#' @return
checkInput <- function(input) {
  if (input$ids != "N/A"
      & input$time != "N/A"
      & input$group != "N/A"
      & input$event != "N/A"
      & input$target_gene != "N/A"){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


#' `getGeneName`
#' @param info
#' @return
getGeneName <- function(info, geneIDS) {
  x <- str_extract(info['INFO'], "(?<=ensembl_gene_id=)[^;]+")
  return(geneIDS[x,]$GeneName)
}


get_cov_list <- function(metadata, input) {
  inputs_list <- c(input$event, input$time, input$ids)
  cov_list <- colnames(metadata)
  cov_list <- cov_list[!(cov_list %in% inputs_list)]
  return(cov_list)
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
  sample_disease_gene_df <- data.frame(matrix(0, ncol = ncol, nrow = nrow))
  colnames(sample_disease_gene_df) <- disease_genes_names
  rownames(sample_disease_gene_df) <- sample_names

  for (sv_idx in 1:length(genes_with_svs_in_sample)) {
    # if the sv is in a gene of interest
    sv_gene_name <- genes_with_svs_in_sample[sv_idx]
    if(is.na(sv_gene_name))
      next
    for(individual_idx in 1:length(sample_names)){
      indiv_id <- sample_names[individual_idx]
      gt <- vcf@gt[sv_idx, indiv_id]
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
`%+=%` <- function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

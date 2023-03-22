

#' VariantSurvival
#'
#' @param vcffile path to the vcf file containing the Structural variant data
#' @param metadatafile path to the txt file containing the samples metadata
#' @param demo
#'
#' @return
#' @export
#'
#' @examples
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
  disease_type_gene <- read_csv("disease_type_gene.csv")
  #
  gene_ids_table <- read.csv(file = 'ensembleTogenes.csv')
  gene_ids_table <- unique(gene_ids_table)
  rownames(gene_ids_table) <- gene_ids_table$ensembleID


  ui <- bootstrapPage(
    navbarPage(theme = shinytheme("flatly"),
               collapsible = TRUE,
               HTML('<a style="text-decoration:none;
               cursor:default;
                    color:#FFFFFF;
                    " class="active" href="#">VariantSurvival</a>'),
               id="nav",
               windowTitle ="VariantSurvival",
               #########################" tab 1 ################################################
               tabPanel("Select Target Gene",
                        sidebarLayout(
                          sidebarPanel(
                            pickerInput(inputId ="disease_n",
                                        label = "Select the disease of interest:",
                                        choices = c(colnames(disease_gene), "N/A"),
                                        selected = "N/A"
                            ),

                            #
                            h2("Annotate metadata :"),
                            selectInput(inputId = "ids",
                                        label = "Select the participant ID:",
                                        choices = c(colnames(metadata), "N/A"),
                                        selected = "N/A"
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
                           )
                          ),
                          mainPanel(
                            h2("The following table resume the  Biomarkers associated
                               with the diseases based on the ClinGen database"),
                            span(DT::dataTableOutput("table")),
                            #
                            br(),br(),br(),br(),br(),br(),
                            h2("Select Target Gene"),
                           box(width = 6,
                               selectInput(inputId = "target_gene",
                                           label = "Gene of interest:",
                                           choices = NULL,
                                           selected = FALSE
                               ),
                               br(),
                               DT::dataTableOutput("summ_table") ),
                            box(width = 6,plotOutput(outputId = "histogram") )
                          )
                        )
               ),
               #######################""" tab2 ###############################
               tabPanel("Kaplan–Meier",

                          # to enable/disable the n_svs_min & n_svs_max fields,
                          #this should be added at the panel level
                          shinyjs::useShinyjs(),
                          tabBox(
                            tabPanel("Null model",
                                     checkboxInput("life_table_null_model","Display life table",value = FALSE),
                                     fixedPanel(title = "", draggable = TRUE, left=50,width="50%",

                                                shinycssloaders::withSpinner(plotOutput(outputId = "null_model_km"))
                                     ),

                                     fixedPanel(title = "",  draggable = TRUE, right  = 50,id = "myBox", width="50%",
                                                span(DT::dataTableOutput("null_model_life_table"))
                                     )
                            ),
                            tabPanel("Multiple model",
                                     checkboxInput("life_table_multiple_model","Display life table",value = FALSE),
                                     fixedPanel(title = "",  draggable = TRUE, left  = 50, #id = "",width="50%",
                                                dropdownButton(
                                                  checkboxInput("all_n_svs",
                                                                "Include all counts",
                                                                value = TRUE
                                                  ),
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
                                                  checkboxGroupInput("km_feat",
                                                                     "Plot layout options:",
                                                                     choices = c("confidence interval" = "conf_itv",
                                                                                 "risk table" = "risk_table",
                                                                                 "y grid line" = "grid_line")),
                                                  circle = TRUE,
                                                  status = "danger",
                                                  icon = icon("gear"), width = "200",
                                                  tooltip = tooltipOptions(title = "Click to see inputs !")
                                                ),#button
                                                shinycssloaders::withSpinner(plotOutput(outputId = "plot_km"#,width = "100%"
                                                ))
                                     ),
                                     fixedPanel(title = "",  draggable = TRUE, right  = 50,id = "myBox",width="50%",
                                                span(DT::dataTableOutput("multiple_model_life_table"))
                                     )

                            )
                          )),
               ###########################################"""
               tabPanel("Cox regression",
                        span(shiny::tags$i(
                          h3("Cox regression table")),
                          style="color:#045a8d",
                          shinyjs::useShinyjs(),
                          checkboxInput("cox_reg_td",
                                        "With time-dependent covariates",
                                        value = TRUE),
                          selectizeInput(inputId = "sel_cov",
                                         label = "Select categorical covariates",
                                         # SV_bin is added by us, 0/1 without/with SV
                                         choices = NULL,
                                         selected = FALSE,
                                         multiple = TRUE
                          ),
                          selectizeInput(inputId = "sel_cov_cont",
                                         label = "Select continuous covariates",
                                         choices = NULL,
                                         selected = FALSE,
                                         multiple = TRUE
                          ),
                          selectInput(inputId = "sel_strata",
                                      label = "Select strata covariate (optional)",
                                      choices = NULL),
                         box(span(DT::dataTableOutput("table3")))

                        )
               )
                          )
                        )

  server <- function(input, output, session) {

# all gens
output$table <- DT::renderDataTable({
  disease_type_gene <- disease_type_gene %>% filter(disease_type_gene$GCEP ==input$disease_n )
  disease_type_gene <-unique( disease_type_gene[,-c(2,4,6,9,10)])
  disease_type_gene %>% arrange(DISEASE_LABEL)
})
##
# add % to table
##
# target gene
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

################################################################################
#extract svs from vcf
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
    combine_newdf_metadata <-merge ( metadata,new_df,  by.x = c(1), by.y ='ids')
    #patients
    genes_in_patients <- combine_newdf_metadata %>% filter(Phenotype =="1")
    genes_in_patients <- genes_in_patients[,-c(2,3,4,5,6,7)]
    #invert rows and columns
    genes_in_patients2 <- data.frame(t(genes_in_patients[-1]))
    colnames(genes_in_patients2) <- genes_in_patients[, 1]

    #max sv in each gene in patients
    genes_in_patients2$max_patient <- apply(genes_in_patients2, 1, max, na.rm=TRUE)
    #min
    genes_in_patients2$min_patient <- apply(genes_in_patients2, 1, min, na.rm=TRUE)
    #control
    genes_in_controls <- combine_newdf_metadata %>% filter(Phenotype =="0")
    genes_in_controls <- genes_in_controls[,-c(2,3,4,5,6,7)]
    #invert rows and columns
    genes_in_controls2 <- data.frame(t(genes_in_controls[-1]))
    colnames(genes_in_controls2) <- genes_in_controls[, 1]

    #max sv in each gene in patients
    genes_in_controls2$max_controls <- apply(genes_in_controls2, 1, max, na.rm=TRUE)
    #min
    genes_in_controls2$min_controls <- apply(genes_in_controls2, 1, min, na.rm=TRUE)
    #combine all in a table
    genes_in_patients2 <- genes_in_patients2 %>%
      rownames_to_column("gene")


    genes_in_patients2  <- genes_in_patients2[c("gene", "min_patient", "max_patient")]
    #
    genes_in_controls2 <- genes_in_controls2 %>%
      rownames_to_column("gene")
    genes_in_controls2 <- genes_in_controls2[c("gene", "min_controls", "max_controls")]
    #% of patients samples with SV
    genes_in_patients3 <-data.frame(t(genes_in_patients[-1]))
    colnames(genes_in_patients3) <- genes_in_patients[, 1]
    genes_in_patients3$percentage_sample_patients_with_sv <- (as.numeric(rowSums(genes_in_patients3!=0))
                                                              / as.numeric(ncol(genes_in_patients3)))*100
    genes_in_patients3 <- genes_in_patients3 %>%
      rownames_to_column("gene")
    genes_in_patients3 <- genes_in_patients3[c("gene", "percentage_sample_patients_with_sv")]

    #% of control samples with SV
    genes_in_controls3 <-data.frame(t(genes_in_controls[-1]))
    colnames(genes_in_controls3) <- genes_in_controls[, 1]
    genes_in_controls3$percentage_sample_controls_with_sv <- (as.numeric(rowSums(genes_in_controls3!=0))
                                                              / as.numeric(ncol(genes_in_controls3)))*100
    genes_in_controls3 <- genes_in_controls3 %>%
      rownames_to_column("gene")
    genes_in_controls3 <- genes_in_controls3[c("gene", "percentage_sample_controls_with_sv")]
    #merge
    output$summ_table <- DT::renderDataTable({
      all_tabble = list( genes_in_controls2,genes_in_patients2,genes_in_patients3,genes_in_controls3)
      statistics_table <- all_tabble %>% reduce(inner_join, by='gene')
      statistics_table <- statistics_table %>% filter(statistics_table$gene ==input$target_gene)
      svgene <- data.frame(t(statistics_table[-1]))
      colnames(svgene) <-statistics_table[, 1]
      options(digits=3)
      svgene
    })
  })
}
  })
####
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
                 aes(label=after_stat(if_else(condition = count>0,
                                              as.character(count), ""))),
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

############################""" tab 2 km #################################
observe({
  if(checkInput(input)){
    svs_gene_input_df <- reactive_no_NAs_metadata()
    null_model <- survfit(Surv(time, event)~1,
                          data = svs_gene_input_df,
                          type = "kaplan-meier")
    output$null_model_km <- renderPlot({
      ggsurvplot_combine(list(null_model),
                         data=no_na_df,
                         risk.table = FALSE,
                         conf.int = FALSE,
                         conf.int.style = "step",
                         risk.table.y.text = FALSE,
                         xlab = "Years",
                         ylab = "Overall survival probability",
                         legend = "none",
                         ggtheme = theme(
                           text = element_text(size = 20),
                           panel.background = element_rect(fill = "white",
                                                           colour = "white"),
                           axis.line = element_line(colour = "black"),
                           panel.grid.major.y = element_line(colour='white'),
                           panel.grid.minor.y = element_line(colour='white'),
                           # legend.position = c(0.2, 0.5)
                         ))
    },
    height = 500,
    width = 700)

    observeEvent(input$life_table_null_model, {
      if(input$life_table_null_model == FALSE){
        shinyjs::hide(id = "myBox")
      }
      else if(input$life_table_null_model==TRUE){
        shinyjs::show(id = "myBox")
        output$null_model_life_table <- DT::renderDataTable({
          as_tibble(round_df(as.data.frame(surv_summary(null_model)), 3))
        })
      }
    })
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
      if (input$n_svs_max > input$n_svs_min){
        n_svs_min = as.numeric(input$n_svs_min)
        n_svs_max = as.numeric(input$n_svs_max)
        sub <- (svs_gene_input_df$SV_count_per_gene >= n_svs_min
                & svs_gene_input_df$SV_count_per_gene <= n_svs_max)
        svs_gene_input_df <- svs_gene_input_df[sub,]
      }
    }
    groups_df <- (unique(svs_gene_input_df[c("SV_bin", "trial_group_bin")])
                  %>% arrange(SV_bin, trial_group_bin))
    lables <- c()
    cols <- c()
    for (row in 1:nrow(groups_df)){
      sv_bin <- groups_df[row, "SV_bin"]
      trial_g <- groups_df[row, "trial_group_bin"]
      if (sv_bin == 0){
        if (trial_g == 0){
          lables <- c(lables, "without sv - placebo")
          cols <- c(cols, "Violetred2")
        }
        else if (trial_g == 1){
          lables <- c(lables, "without sv - treatment")
          cols <- c(cols, "turquoise3")
        }
      }
      if (sv_bin == 1){
        if (trial_g == 0){
          lables <- c(lables, "with sv - placebo")
          cols <- c(cols, "Violetred4")
        }
        else if (trial_g == 1){
          lables <- c(lables, "with sv - treatment")
          cols <- c(cols, "steelblue")
        }
      }
    }
    n_sv_groups <- unique(svs_gene_input_df[["SV_bin"]])
    # can it happen that there are more than 1 patient groups?
    if(length(n_sv_groups) == 2){
      # generate survival curve objects for each group
      sc_without <- survfit2(Surv(time, event)~trial_group_bin,
                             data = svs_gene_input_df,
                             subset=(SV_bin==0),
                             type = "kaplan-meier")
      sv_with <- survfit2(Surv(time, event)~trial_group_bin,
                          data = svs_gene_input_df,
                          subset=(SV_bin==1),
                          type = "kaplan-meier")
      surv_fit_list <- list("with SV" = sv_with, "without SV" = sc_without)
    }
    else if(length(n_sv_groups) ==1){
      if (n_sv_groups == 0){
        type = "without"}
      else if (n_sv_groups == 1){
        type = "with"}
      temp <- paste(type , " SV")
      sc <- survfit2(Surv(time, event)~trial_group_bin,
                     data = svs_gene_input_df,
                     type = "kaplan-meier")
      surv_fit_list <- list(temp = sc)
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
                           # legend.position = c(0.2, 0.5)
                         ),
                         legend = "right",
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
                                    input = "sel_cov_cont",
                                    choices = cov_list,
                                    options = list(create = TRUE))
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
    if (input$all_n_svs == FALSE & input$n_svs_min!=""){
      svs_gene_input_df <- reactive_no_NAs_metadata()
      svs_levels <- unique(svs_gene_input_df["SV_count_per_gene"])$SV_count_per_gene
      max_levels <- sort(svs_levels[svs_levels >= as.numeric(input$n_svs_min)])
      updateSelectizeInput(session,
                           input = "n_svs_max",
                           choices = max_levels)
    }
  }
}
)

#############################" cox

observe({
  if(checkInput(input) & (any(!is.na(input$sel_cov)) | any(!is.na(input$sel_cov_cont)))){
    covariates <- c()
    input_df <- reactive_no_NAs_metadata()
    if (any(!is.na(input$sel_cov_cont))){
      input_cov_cont <- map_col_names(input, input$sel_cov_cont)
      input_df[input_cov_cont] <- sapply(input_df[input_cov_cont],as.numeric)
      covariates <- c(covariates, input_cov_cont)
    }
    if(any(!is.na(input$sel_cov))){
      input_cov_cat <- map_col_names(input, input$sel_cov)
      input_df[input_cov_cat] <- sapply(input_df[input_cov_cat],as.character)
      covariates <- c(covariates, input_cov_cat)
    }
    if(input$sel_strata != "N/A"){
      covariates <- c(covariates[covariates != input$sel_strata], sprintf("strata(%s)", input$sel_strata))
    }
    # Standard model
    formulaString <- paste("Surv(time, event) ~", paste(covariates, collapse="+"))
    cox_reg.std <- coxph(as.formula(formulaString), data=input_df)
    res.std <- cox.zph(cox_reg.std)
  }


  output$table3 <- DT::renderDataTable({
    if(checkInput(input) & (any(!is.na(input$sel_cov)) | any(!is.na(input$sel_cov_cont)))){
      covariates <- c()
      input_df <- reactive_no_NAs_metadata()
      if (any(!is.na(input$sel_cov_cont))){
        input_cov_cont <- map_col_names(input, input$sel_cov_cont)
        input_df[input_cov_cont] <- sapply(input_df[input_cov_cont],as.numeric)
        covariates <- c(covariates, input_cov_cont)
      }
      if(any(!is.na(input$sel_cov))){
        input_cov_cat <- map_col_names(input, input$sel_cov)
        input_df[input_cov_cat] <- sapply(input_df[input_cov_cat],as.character)
        covariates <- c(covariates, input_cov_cat)
      }
      formulaString <- paste("Surv(time, event) ~", paste(covariates, collapse="+"))
      x3 <- (coxph(as.formula(formulaString), data=input_df)
             %>% tbl_regression(exp = TRUE))
      proport_hazard_assump <- cox.zph(cox_reg.std)
      t3 <-as_tibble(x3)
      t3
    }
  })

})
  }

  # Run the application
  shinyApp(ui = ui, server = server)
}


############# Helper functions ###################
install_load_requirements<- function() {
  if (!require("shiny")) install.packages("shiny")
  if (!require("shinyjs")) install.packages("shinyjs")
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
################################
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

round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}


map_col_names <- function(input, cov_list){
  if(input$event %in% cov_list){
    cov_list <- str_replace(cov_list,
                            input$event,
                            "event")}
  if(input$group %in% cov_list){
    cov_list <- str_replace(cov_list,
                            input$group,
                            "trial_group_bin")
  }
  return(cov_list)
}


#' implementation of += operator
#' https://stackoverflow.com/questions/5738831/
`%+=%` <- function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

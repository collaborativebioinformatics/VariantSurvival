#' VariantSurvival
#'
#' @param vcf_file path to the vcf file containing the Structural variant data
#' @param metadata_file path to the txt file containing the samples metadata
#' @param demo true or false
#' @importFrom dplyr %>%
#' @export

VariantSurvival <- function(vcf_file, metadata_file, demo = FALSE) {
  days_year <- 365.25
  #demo or input files
  if (demo == TRUE) {
    vcf_file_demo <- system.file(
      "extdata",
      "merged.filtered.vcf",
      package ="VariantSurvival"
      )
    metadata_demo <- system.file(
      "extdata",
      "metadata.xlsx",
      package ="VariantSurvival"
      )
    vcf <- vcfR::read.vcfR(
      vcf_file_demo,
      verbose = FALSE
      )
    metadata <- readxl::read_excel(metadata_demo)
    metadata <- stats::na.omit(metadata)
  } else if (demo == FALSE) {
    vcf <- vcfR::read.vcfR(
      vcf_file, 
      verbose = FALSE
      )
    metadata <- readxl::read_excel(metadata_file)
    metadata <- stats::na.omit(metadata)
  }
  
  disease_genes_path <- system.file(
    "extdata",
    "disease_gene.xlsx",
    package ="VariantSurvival"
    )
  disease_types_path <- system.file(
    "extdata",
    "disease_type_gene.csv",
    package ="VariantSurvival"
    )
  disease_type_gene <- readr::read_csv(
    disease_types_path,
    show_col_types = FALSE
    )
  
  disease_gene <- readxl::read_excel(disease_genes_path)

  ui <- shiny::bootstrapPage(
    shiny::navbarPage(
      theme = shinythemes::shinytheme("flatly"),
      collapsible = TRUE,
      htmltools::HTML(
        '<a style="text-decoration:none;
        cursor:default;
        color:#FFFFFF;
        " class="active" href="#">VariantSurvival</a>'
        ),
      id = "nav",
      windowTitle = "VariantSurvival",
      shiny::tabPanel(
        "Select Target Gene",
        shiny::sidebarLayout(
          shiny::sidebarPanel(
            shinyWidgets::pickerInput(
              inputId = "disease_n",
              label = "Select the disease of interest:",
              choices = c(colnames(disease_gene), "N/A"),
              selected = "N/A"
              ),
            shiny::selectInput(
              inputId = "ids",
              label = "Select the participant ID:",
              choices = c(colnames(metadata), "N/A"),
              selected = "N/A"
              ),
            shiny::selectInput(
              inputId = "time",
              label = "Select the time factor:",
              choices = c(colnames(metadata), "N/A"),
              selected = "N/A"
              ),
            shiny::radioButtons(
              inputId = "time_unit",
              label = "Time factor unit:",
              choices =  c("years", "days"),
              inline = TRUE
              ),
            htmltools::span(
              shiny::tags$i(
              htmltools::h6(
                "The following selections must refer to binary factors"
                )
              ),
              style = "color:#045a8d"
              ),
            shiny::selectInput(
              inputId = "group",
              label = "Select the clinical trial groups factor:",
              choices = c(colnames(metadata), "N/A"),
              selected = "N/A"
              ),
            shiny::selectInput(
              inputId = "event",
              label = "Select the alive/dead factor:",
              choices = c(colnames(metadata), "N/A"),
              selected = "N/A"
              )
            ),
          shiny::mainPanel(
            shiny::fluidRow(
              shiny::column(
                12,
                htmltools::h4(
                "The Biomarkers table resumes the Biomarkers associated
                with the diseases based on the ClinGen database"
                ),
                shinydashboard::tabBox(
                  selected = "Biomarkers Table",
                  shiny::tabPanel(
                    "Biomarkers Table",
                    htmltools::span(
                      DT::dataTableOutput("biomarkers_table")
                      )
                    ),
                  shiny::tabPanel(
                    "Summary",
                    htmltools::span(
                      DT::dataTableOutput("summ_table")
                      )
                    )
                  )
                )
              ),
            lapply(1:2, function(x) shiny::br()),
            htmltools::h4("Select Target Gene"),
            shinydashboard::box(
              width = 6,
              shiny::selectInput(
                inputId = "target_gene",
                label = "Gene of interest:",
                choices = NULL,
                selected = FALSE
                )
              ),
            lapply(
              1:4,
              function(x) shiny::br()
              ),
            shiny::fluidRow(
              shiny::column(
                12,
                shinydashboard::tabBox(
                  selected = "Information Summary",
                  shiny::tabPanel(
                    "Information Summary",
                    DT::dataTableOutput("gene_summary_table_i")
                    ),
                  shiny::tabPanel(
                    "Structural Variants Distribution",
                    shinycssloaders::withSpinner(
                      shiny::plotOutput(outputId = "histogram")
                      )
                    ),
                  shiny::tabPanel(
                    "Participants table",
                    shiny::selectInput(
                      inputId = "table_cols",
                      label = "Include columns (optional)",
                      choices = NULL,
                      selected = FALSE,
                      multiple = TRUE
                      ),
                    shiny::span(
                      DT::dataTableOutput("table")
                      )
                    )
                  )
                )
              
              ),
            lapply(
              1:4,
              function(x) shiny::br()
            ),
            )
          )
        ),
      shiny::tabPanel(
        "Kaplan-Meier",
        shiny::mainPanel(
          shiny::span(
            # to enable/disable the n_svs_min & n_svs_max fields,
            #this should be added at the panel level
            shinyjs::useShinyjs(),
            shinydashboard::tabBox(
              shiny::tabPanel(
                "Multiple model",
                shiny::fluidRow(
                  shiny::column(
                    12,
                    shiny::div(
                      style = "height:70px; Topleft"
                      ),
                    shinyWidgets::dropdownButton(
                      shiny::checkboxInput(
                        inputId = "all_n_svs",
                        label = "Include all counts",
                        value = TRUE
                        ),
                      shiny::selectInput(
                        inputId = "n_svs_min",
                        label = "min",
                        choices = NULL,
                        selected = FALSE
                        ),
                      shiny::selectInput(
                        inputId = "n_svs_max",
                        label = "max",
                        choices = NULL,
                        selected = FALSE
                        ),
                      shiny::checkboxGroupInput(
                        "km_feat",
                        "Plot layout options:",
                        choices = c(
                          "confidence interval" = "conf_itv",
                          "risk table" = "risk_table",
                          "y grid line" = "grid_line",
                          "life table" = "life_table"
                          )
                        ),
                      circle = TRUE,
                      status = "danger",
                      icon = shiny::icon("gear"),
                      width = "300px",
                      tooltip = shinyWidgets::tooltipOptions(
                        title = "Click to see inputs!"
                        )
                      ),
                    shinycssloaders::withSpinner(
                      shiny::plotOutput(
                        outputId = "plot_km",
                        width = "100%"
                        )
                      ),
                    )
                  ),
                shiny::fluidRow(
                  shiny::column(
                  12,
                  shiny::div(
                    style = "height:200px; Bottomleft"
                    ),
                  shiny::span(
                    DT::dataTableOutput("p_values_km")
                    )
                  )
                  ),
                shiny::conditionalPanel(
                  condition = "input.km_feat.includes('life_table')",
                  shiny::fluidRow(
                    shiny::column(
                      12,
                      shiny::div(
                        style = "height:200px; Bottomleft"
                      ),
                      shinydashboard::tabBox(
                        id = "myBox_mm",
                        title = "",
                        width = 12,
                        shiny::tabPanel(
                          "with SV",
                          shiny::span(
                            DT::dataTableOutput("lt_mm_0")
                          )
                        ),
                        shiny::tabPanel(
                          "without SV",
                          shiny::span(
                            DT::dataTableOutput("lt_mm_1")
                          )
                        )
                      )
                    )
                  )
                )
                ),
              shiny::tabPanel(
                "Null model",
                shiny::fluidRow(
                  shiny::column(
                    12,
                    shiny::div(
                      style = "height:70px; Topleft"
                      ),
                    shinycssloaders::withSpinner(
                      shiny::plotOutput(
                        outputId = "null_model_km"
                        )
                      )
                    )
                  ),
                shiny::fluidRow(
                  shiny::column(
                    12,
                    shiny::div(
                      style = "height:200px; Bottomleft"
                      ),
                    shinydashboard::box(
                      id = "myBox",
                      title = "",
                      width = '200px',
                      shiny::span(
                        DT::dataTableOutput(
                          "null_model_life_table"
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        ),
      shiny::tabPanel(
        "Cox regression",
        shiny::span(
          shiny::tags$i(
            htmltools::h3("Cox regression table")
            ),
          style = "color:#045a8d",
          shinyjs::useShinyjs(),
          shiny::selectizeInput(
            inputId = "sel_cov",
            label = "Select categorical covariate(s)",
            # SV_bin is added by us, 0/1 without/with SV
            choices = NULL,
            selected = FALSE,
            multiple = TRUE
            ),
          shiny::selectizeInput(
            inputId = "sel_cov_cont",
            label = "Select numerical covariate(s)",
            choices = NULL,
            selected = FALSE,
            multiple = TRUE
            ),
          shiny::selectInput(
            inputId = "sel_strata",
            label = "Select strata covariate (optional)",
            choices = NULL
            ),
          shinydashboard::tabBox(
            shiny::tabPanel(
              "Standard model",
              shiny::span(
                DT::dataTableOutput("summ_std")
                ),
              lapply(
                1:2,
                function(x) shiny::br()
                ),
              shiny::span(
                DT::dataTableOutput("prop_h_std")
                ),
              shiny::br(),
              shiny::span(
                shiny::tags$i(
                  htmltools::h2("Residuals")
                  )
                ),
              shinycssloaders::withSpinner(
                shiny::plotOutput(outputId = "residues_std")
                )
              ),
            shiny::tabPanel(
              "Multiple model",
              shinydashboard::tabBox(
                width = 12,
                shiny::tabPanel(
                  "with SV",
                  shiny::span(
                    DT::dataTableOutput("summ_mm_1")
                    ),
                  lapply(
                    1:2,
                    function(x) shiny::br()
                    ),
                  shiny::tabPanel(
                    "with",
                    shiny::span(
                      DT::dataTableOutput("prop_h_mm0")
                      ),
                    shiny::br(),
                    shiny::span(
                      shiny::tags$i(
                        htmltools::h2("Residuals")
                        )
                      ),
                    shinycssloaders::withSpinner(
                      shiny::plotOutput(
                        outputId = "residues_std_mm0"
                        )
                      )
                    ),
                  shiny::br()
                  ),
                shiny::tabPanel(
                  "without SV",
                  shiny::span(
                    DT::dataTableOutput("summ_mm_0")
                    ),
                  lapply(
                    1:2, 
                    function(x) shiny::br()
                    ),
                  shiny::tabPanel(
                    "without",
                    shiny::span(
                      DT::dataTableOutput("prop_h_mm1")
                      ),
                    shiny::br(),
                    shiny::span(
                      shiny::tags$i(
                        htmltools::h2("Residuals")
                        )
                      ),
                    shinycssloaders::withSpinner(
                      shiny::plotOutput(
                        outputId = "residues_std_mm1"
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  
  server <- function(input,
                     output,
                     session
                     ) {
    ensmbl_gene_ids <- system.file(
      "extdata",
      'ensemblTogenes.csv',
      package ="VariantSurvival"
      )
    gene_ids_table <- utils::read.csv(file=ensmbl_gene_ids)
    rownames(gene_ids_table) <- gene_ids_table$ensembleID
    # get disease_n input and update the genes list accordingly
    shiny::observeEvent(
      input$disease_n, {
        if (input$disease_n != "N/A") {
        genes_list <- c(get_disease_gene_list(disease_gene,input$disease_n))
        shiny::updateSelectizeInput(
          session,
          input = "target_gene",
          choices = c(genes_list, 'N/A'),
          selected = 'N/A'
          )
        disease_type_gene <- disease_type_gene %>% 
          dplyr::filter(
            GCEP == input$disease_n
            )
        disease_type_gene <- unique(
          disease_type_gene[, -c(2, 4, 6, 9, 10)]
          )
        # Convert URLs to HTML links
        disease_type_gene$ONLINE_REPORT <- sapply(
          disease_type_gene$ONLINE_REPORT, function(x) {
          paste0('<a href="', x, '" target="_blank">', x, '</a>')
        })
        
        output$biomarkers_table <- DT::renderDataTable(
          {
            disease_type_gene %>% 
              dplyr::arrange(DISEASE_LABEL)
            },
          escape = FALSE,  # Important to render HTML
          options = list(
            autoWidth = TRUE,
            columnDefs = list(
              list(width = '200px', targets = "_all")  # Adjust column width as needed
            )
          )
        )
        }
        }
      )
    observe(
      {
        if (input$disease_n != "N/A"
            & input$ids != "N/A"
            & input$time != "N/A") {
          count_table <- genesCountTable(
            vcf,
            metadata, 
            input,
            gene_ids_table
            )
          output$summ_table <- DT::renderDataTable(
            {
              # count data frame with patient_ids in rows and gene_ids in columns
              count_table
              }
          )
          }
        }
      )

    reactive_no_NAs_metadata <- reactive(
      {
        if (checkInput(input)) {
          disease_genes_names <- gene_ids_table$GeneName
          if (input$target_gene %in% disease_genes_names) {
            sample_names <- colnames(vcf@gt)[-1] # ignore VCF genotype information
            # genes are repeated since a single gene can have more than one SV.
            genes_with_svs_in_sample <- apply(
              vcf@fix,
              1,
              getGeneName, 
              gene_ids_table
              )
          # count dataframe with patient_ids in rows and gene_ids in columns
          count_df <- CountSVsDf(
            length(disease_genes_names),
            length(sample_names),
            disease_genes_names,
            sample_names,
            genes_with_svs_in_sample,
            vcf,
            input$ids
            )
          new_md <- (
            Remove_NA(metadata, input$time) %>%
              dplyr::rename(
                ids = input$ids,
                trial_group_bin = input$group,
                time = input$time,
                event = input$event
                )
            )
          # what happen if the input$target_gene is not in the vcf file?
          new_df <- (
            subset(count_df, ids %in% new_md$ids) %>%
              dplyr::rename(SV_count_per_gene = input$target_gene)
            )
          no_na_df <- merge(new_df, new_md, by = "ids")
          no_na_df <- transform(
            no_na_df,
            time = as.numeric(time),
            event = as.numeric(event),
            SV_count_per_gene = as.numeric(SV_count_per_gene)
            )
          no_na_df <- (no_na_df %>%
                         dplyr::mutate(
                           SV_bin = ifelse(
                             new_df$SV_count_per_gene > 0, 1, 0)
                           ) %>%
                         dplyr::mutate(
                           trial_group = ifelse(
                             trial_group_bin == "0", "Placebo", "Treatment")
                           ) %>% 
                         dplyr::mutate(
                           SV = ifelse(
                             SV_bin == "0", "with", "without"
                             )
                           )
                       )
          if (input$time_unit == 'years') {
            no_na_df$time_days <- floor(no_na_df[["time"]] * days_year)
            } else{
              no_na_df$time_years <- no_na_df[["time"]] / days_year
              }
          no_na_df
        }
        else{
          print(
            paste(
              'Please include the gene ID ',
              input$target_gene,
              'in the ensemblTogenes.csvfile'
            )
          )
          genes_list <- c(get_disease_gene_list(disease_gene, input$disease_n))
          shiny::updateSelectizeInput(
            session,
            input = "target_gene",
            choices = c(genes_list, 'N/A'),
            selected = 'N/A'
          )
        }
      }
    }
    )

    observe(
      {
      if (checkInput(input)) {
        new_df <- reactive_no_NAs_metadata()
        count_tab <- genesCountTable(
          vcf, 
          metadata,
          input,
          gene_ids_table
          )
        target_gene_df <- new_df[, c("SV_count_per_gene", "trial_group_bin")]
        # subset the count of structural variants in the control/treatment patients
        control <- target_gene_df$SV_count_per_gene[target_gene_df$trial_group_bin == 0]
        number_of_patients <- count_tab %>%
          dplyr::filter(count_tab$Gene_ID == input$target_gene) %>%
          dplyr::select(`Count of patients with SVs`)

        if (number_of_patients$`Count of patients with SVs` == 0) {
          output$gene_summary_table <- DT::renderDataTable(
            {
              print(
                DT::datatable("No patient carries structural
                                 variants in this gene")
              )
              }
            )
          output$gene_summary_table_i <- DT::renderDataTable({
            print(
              dplyr::tibble("No patient carries structural
                            variants in this gene")
              )
            }
          )
          }
        else{
          output$gene_summary_table_i <- DT::renderDataTable({
            sketch <- htmltools::withTags(
              table(class = 'display',
                    thead(
                      tr(
                        th(rowspan = 2, 'Group'),
                        th(colspan = 2, 'Patients with SVs in target gene')),
                      tr(
                        lapply(
                          rep(c('Total',  'Percentage over entire sample'), 1),
                          th
                          )
                        )
                      )
                    )
              )
            rows_names <- c('control', 'treatment')
            length(control[control > 0])
            control <- target_gene_df$SV_count_per_gene[target_gene_df$trial_group_bin == 0]
            treatment <- target_gene_df$SV_count_per_gene[target_gene_df$trial_group_bin == 1]
            n_control <- length(control[control > 0])
            n_trearment <- length(treatment[treatment > 0])
            n_indiv <- nrow(target_gene_df)
            first_column <- c(n_control, n_trearment)
            second_column <- c(round(n_control / n_indiv, 2),
                               round(n_trearment / n_indiv, 2)
                               )
            df <- data.frame(
              rows_names,
              first_column, 
              second_column
              )
            DT::datatable(
              df,
              container = sketch,
              rownames = FALSE,
              options = list(dom = 't')
            )
            }
          )
        }
        }
      }
      )
    ## output - histogram ##
    output$histogram <- renderPlot({
      if (checkInput(input)) {
        new_df <- reactive_no_NAs_metadata()
        # get a dataframe with counts
        svs_gene_input_df <- new_df[c("ids", "trial_group", "SV_count_per_gene")]
        legend_title <- "Group"
        ggplot2::ggplot(
          svs_gene_input_df,
          ggplot2::aes(SV_count_per_gene, fill = trial_group)) +
          ggplot2::geom_histogram(binwidth = 1) +
          ggplot2::stat_bin(
            binwidth=1,
            geom='text',
            color='white',
            ggplot2::aes(
              label=ggplot2::after_stat(
                dplyr::if_else(
                  condition=count > 0,
                  as.character(count),
                  ""
                  )
                )
              ),
            position=ggplot2::position_stack(vjust = 0.5)
            ) +
          ggplot2::xlab("Number of SVs in target gene") + 
          ggplot2::ylab("Number of patients") +
          ggplot2::theme(
            text = ggplot2::element_text(size = 20),
            panel.background = ggplot2::element_rect(
              fill = "white",
              colour = "white"
              ),
            axis.line = ggplot2::element_line(colour = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "top"
          ) +
          ggplot2::scale_fill_manual(
            legend_title,
            values = c("#8B1D4D", "#5275A6")
            )
        }
      },
    height = 500,
    width = 700
    )

    observe({
      if (checkInput(input) & !is.null(input$target_gene)) {
        count_tab <- genesCountTable(vcf, metadata, input, gene_ids_table)
        non_zero <- count_tab[count_tab$`Count of patients with SVs`!=0,]
        shiny::updateSelectizeInput(
          session,
          input = "table_cols",
          choices = c(non_zero$Gene_ID,colnames(metadata)),
          selected = NULL
        )
        }
    })

    observe({
      if (checkInput(input)) {
        svs_gene_input_df <- reactive_no_NAs_metadata()
        # renaming back to original column names, so it's not confusing for the user.
        # maybe not the best approach
        rename_cols <- c("ids", "trial_group_bin", "time", "event")
        rename_cols_with <- c(input$ids, input$group, input$time, input$event)
        names(svs_gene_input_df)[names(svs_gene_input_df) %in% rename_cols] <-
          rename_cols_with
        if (any(!is.na(input$table_cols))) {
          output$table <- DT::renderDataTable({
            col_names_list <- c(c(input$ids, "SV_count_per_gene"), input$table_cols)
            svs_gene_input_df <- tibble::as_tibble(svs_gene_input_df[col_names_list])
            svs_gene_input_df
          })
        }
        else{
          output$table <- DT::renderDataTable({
            svs_gene_input_df <- tibble::as_tibble(svs_gene_input_df[c(input$ids, "SV_count_per_gene")])
            colnames(svs_gene_input_df) <- c("ID", "SV count in gene of interest")
            svs_gene_input_df
          })
        }
      }
    })

    observe({
      if (checkInput(input)) {
        svs_gene_input_df <- reactive_no_NAs_metadata()
        null_model <- survival::survfit(
          survival::Surv(time, event) ~ 1,
          data = svs_gene_input_df,
          type = "kaplan-meier"
          )
        output$null_model_km <- renderPlot({
          survminer::ggsurvplot_combine(
            list(null_model),
            data=no_na_df,
            risk.table=FALSE,
            conf.int=FALSE,
            conf.int.style="step",
            risk.table.y.text=FALSE,
            xlab="Years",
            ylab="Overall survival probability",
            legend="none",
            ggtheme=ggplot2::theme(
              text=ggplot2::element_text(size = 20),
              panel.background=ggplot2::element_rect(
                fill="white",
                colour="white"
                ),
              axis.line=ggplot2::element_line(colour = "black"),
              panel.grid.major.y=ggplot2::element_line(
                colour='white'
                ),
              panel.grid.minor.y=ggplot2::element_line(
                colour='white'
                ),
              # legend.position = c(0.2, 0.5)
            )
          )
        },
        height = 500,
        width = 700)

        shiny::observeEvent(input$life_table_null_model, {
          if (input$life_table_null_model == FALSE) {
            shinyjs::hide(id = "myBox")
          }
          else if (input$life_table_null_model == TRUE) {
            shinyjs::show(id = "myBox")
            output$null_model_life_table <- DT::renderDataTable({
              tibble::as_tibble(round_df(as.data.frame(
                survminer::surv_summary(null_model)
                ),
                3)
                )
            })
          }
        })
      }
    })
    
    observe({
      if (checkInput(input)) {
        risk_table <- FALSE
        conf_itv <- FALSE
        life_table <- FALSE
        grid_line <- "white"
        if (!is.null(input$km_feat)) {
          if ("conf_itv" %in% input$km_feat) {
            conf_itv <- TRUE
          }
          if ("risk_table" %in% input$km_feat) {
            risk_table <- TRUE
          }
          if ("grid_line" %in% input$km_feat) {
            grid_line <- "grey"
          }
          if ("life_table" %in% input$km_feat) {
            life_table <- TRUE
          }
        }
        
        svs_gene_input_df <- reactive_no_NAs_metadata()
        if (input$all_n_svs == FALSE
            & input$n_svs_min != ""
            & input$n_svs_max != "") {
          if (input$n_svs_max > input$n_svs_min) {
            n_svs_min <- as.numeric(input$n_svs_min)
            n_svs_max <- as.numeric(input$n_svs_max)
            sub <- (svs_gene_input_df$SV_count_per_gene >= n_svs_min
              & svs_gene_input_df$SV_count_per_gene <= n_svs_max
              )
            svs_gene_input_df <- svs_gene_input_df[sub, ]
          }
        }
        groups_df <- (unique(
          svs_gene_input_df[c("SV_bin", "trial_group_bin")]
          ) %>% 
            dplyr::arrange(SV_bin, trial_group_bin)
          )
        lables <- c()
        cols <- c()
        for (row in 1:nrow(groups_df)) {
          sv_bin <- groups_df[row, "SV_bin"]
          trial_g <- groups_df[row, "trial_group_bin"]
          if (sv_bin == 0) {
            if (trial_g == 0) {
              lables <- c(lables, "without sv - placebo")
              cols <- c(cols, "Violetred2")
            }
            else if (trial_g == 1) {
              lables <- c(lables, "without sv - treatment")
              cols <- c(cols, "turquoise3")
            }
          }
          if (sv_bin == 1) {
            if (trial_g == 0) {
              lables <- c(lables, "with sv - placebo")
              cols <- c(cols, "Violetred4")
            }
            else if (trial_g == 1) {
              lables <- c(lables, "with sv - treatment")
              cols <- c(cols, "steelblue")
            }
          }
        }
        n_sv_groups <- unique(svs_gene_input_df[["SV_bin"]])
        # can it happen that there are more than 1 patient groups?
        if (length(n_sv_groups) == 2) {
          # generate survival curve objects for each group
          sc_without <- ggsurvfit::survfit2(
            survival::Surv(time, event) ~ trial_group_bin,
            data = svs_gene_input_df,
            subset = (SV_bin == 0),
            type = "kaplan-meier"
          )
          sv_with <- ggsurvfit::survfit2(
            survival::Surv(time, event) ~ trial_group_bin,
            data = svs_gene_input_df,
            subset = (SV_bin == 1),
            type = "kaplan-meier"
          )
          surv_fit_list <- list(
            "with" = sv_with,
            "without" = sc_without
            )
        }
        else if (length(n_sv_groups) == 1) {
          if (n_sv_groups == 0) {
            type <- "without"
          }
          else if (n_sv_groups == 1) {
            type <- "with"
          }
          temp <- paste(type , " SV")
          sc <- ggsurvfit::survfit2(
            survival::Surv(time, event) ~ trial_group_bin,
            data = svs_gene_input_df,
            type = "kaplan-meier"
            )
          surv_fit_list <- list(temp = sc)
        }
        output$plot_km <- renderPlot({
          survminer::ggsurvplot_combine(
            surv_fit_list,
            data = svs_gene_input_df,
            risk.table = risk_table,
            conf.int = conf_itv,
            conf.int.style = "step",
            risk.table.y.text = FALSE,
            xlab = tools::toTitleCase(input$time_unit),
            ggtheme = ggplot2::theme(
              text = ggplot2::element_text(size = 20),
              panel.background = ggplot2::element_rect(
                fill = "white",
                colour = "white"
                ),
              axis.line = ggplot2::element_line(
                colour = "black"
                ),
              panel.grid.major.y = ggplot2::element_line(
                colour=grid_line
                ),
              panel.grid.minor.y = ggplot2::element_line(
                colour=grid_line
                ),
              # legend.position = c(0.2, 0.5)
            ),
            legend = "right",
            legend.title = "Group",
            legend.labs = lables,
            palette = cols
          )
        },
        height = 600,
        width = 1000
        )

        pval_tab <- survminer::surv_pvalue(
          surv_fit_list,
          data=svs_gene_input_df,
          combine = TRUE
          )
        pval_tab <- pval_tab[c("id","method","pval.txt")]
        colnames(pval_tab) <- c("Structural Variant", "Method", "p-value")

        output$p_values_km <- DT::renderDataTable(
          
          DT::datatable(
            tibble::as_tibble(pval_tab),
            options = list(
              paging = FALSE, 
              searching = FALSE,
              lengthChange = FALSE,
              info = FALSE  # Hides the 'Showing 1 to N of N entries' information
            )
          )
        )
        
        if (life_table == TRUE){
          shinyjs::show(id = "myBox_mm")
          output$lt_mm_0 <- DT::renderDataTable({
            tl_sv_without <- data.frame(
              lapply(
                surv_summary(
                  sc_without,
                  data = svs_gene_input_df
                  ),
                function(x) if (is.numeric(x)) round(x, 3) else x
                )
              )
            tl_sv_without <- tl_sv_without[,!names(tl_sv_without) %in% c("strata")]
            tibble::as_tibble(tl_sv_without)
          })
          output$lt_mm_1 <- DT::renderDataTable({
            tl_sv_with <- data.frame(
              lapply(
                surv_summary(
                  sv_with,
                  data = svs_gene_input_df
                  ),
                function(x) if (is.numeric(x)) round(x, 3) else x
                )
              )
            tl_sv_with <- tl_sv_with[,!names(tl_sv_with) %in% c("strata")]
            tibble::as_tibble(tl_sv_with)
          })
        }
        }
      else{
        output$plot_km <- shiny::renderText({
          "Missing input data"
        })
      }
    })

    inputs_react <- reactive({
      list(input$event,
           input$time,
           input$ids,
           input$target_gene)
    })
    # hide the event, time and ids covariates from the covariates  and strata
    # drop-down
    shiny::observeEvent(inputs_react(), {
      if (checkInput(input)) {
        cov_list <- c(get_cov_list(metadata, input), "SV_bin")
        # update the covariate drop down
        shiny::updateSelectizeInput(
          session,
          input = "sel_cov",
          choices = cov_list,
          options = list(create = TRUE)
        )
        # update the strata drop down, this field is optional
        shiny::updateSelectizeInput(
          session,
          input = "sel_cov_cont",
          choices = cov_list,
          options = list(create = TRUE)
        )
        shiny::updateSelectizeInput(
          session,
          input = "sel_strata",
          choices = c(cov_list, c("SV_bin", "N/A")),
          selected = "N/A"
        )
      }
    })

    shiny::observeEvent(input$all_n_svs, {
      if (checkInput(input)) {
        if (input$all_n_svs == FALSE) {
          shinyjs::enable("n_svs_min")
          shinyjs::enable("n_svs_max")
          svs_gene_input_df <- reactive_no_NAs_metadata()
          svs_levels <- unique(svs_gene_input_df["SV_count_per_gene"])$SV_count_per_gene
          shiny::updateSelectizeInput(
            session,
            input = "n_svs_min",
            choices = sort(svs_levels)
            )
        }
        else {
          shinyjs::disable("n_svs_min")
          shinyjs::disable("n_svs_max")
        }
      }
    })

    shiny::observeEvent(input$n_svs_min, {
      if (checkInput(input)) {
        if (input$all_n_svs == FALSE & input$n_svs_min != "") {
          svs_gene_input_df <- reactive_no_NAs_metadata()
          svs_levels <- unique(svs_gene_input_df["SV_count_per_gene"])$SV_count_per_gene
          max_levels <- sort(svs_levels[svs_levels >= as.numeric(input$n_svs_min)])
          shiny::updateSelectizeInput(
            session,
            input = "n_svs_max",
            choices = max_levels
            )
        }
      }
    })

    observe({
      if (checkInput(input) &
          (any(!is.na(input$sel_cov)) | any(!is.na(input$sel_cov_cont)))) {
        covariates <- c()
        input_df <- reactive_no_NAs_metadata()
        if (any(!is.na(input$sel_cov_cont))) {
          input_cov_cont <- map_col_names(input, input$sel_cov_cont)
          input_df[input_cov_cont] <-
            sapply(input_df[input_cov_cont], as.numeric)
          covariates <- c(covariates, input_cov_cont)
        }
        if (any(!is.na(input$sel_cov))) {
          input_cov_cat <- map_col_names(input, input$sel_cov)
          input_df[input_cov_cat] <- sapply(input_df[input_cov_cat], as.character)
          covariates <- c(covariates, input_cov_cat)
        }
        if (input$sel_strata != "N/A") {
          covariates <- c(covariates[covariates != input$sel_strata],
                          sprintf("strata(%s)", input$sel_strata))
        }
        # Standard model
        formulaString <- paste(
          "survival::Surv(time, event) ~", 
          paste(covariates, collapse = "+")
          )
        cox_reg.std <- survival::coxph(
          stats::as.formula(formulaString),
          data = input_df
          )
        ## Check for violation of proportional hazard
        res.std <- survival::cox.zph(cox_reg.std)


        # Multiple model
        input_df_mm <- data.frame(input_df)
        input_df_mm["SV_bin"] <- sapply(input_df_mm["SV_bin"], as.numeric)
        formulaString_mm <- paste(
          "survival::Surv(time, event) ~",
          paste(covariates[!(covariates %in% c("SV_bin","strata(SV_bin)"))],
                collapse = "+")
          )
        cox_reg.mul0 <- survival::coxph(
          stats::as.formula(formulaString_mm),
          data = input_df_mm,
          subset = (SV_bin == 0)
          )
        cox_reg.mul1 <- survival::coxph(
          stats::as.formula(formulaString_mm),
          data = input_df_mm,
          subset = (SV_bin == 1)
          )
        res.std_mul0 <- survival::cox.zph(cox_reg.mul0)
        res.std_mul1 <- survival::cox.zph(cox_reg.mul1)

        output$summ_std <- DT::renderDataTable({
          tibble::as_tibble(
            cox_reg.std %>% 
              gtsummary::tbl_regression(exp = TRUE)
            )
        })
        output$prop_h_std <- DT::renderDataTable({
          DT::datatable(tibble::rownames_to_column(
            round_df(as.data.frame(res.std$table), 3), " "),
            options = list(dom = 't')
            )
        })
        output$residues_std <- renderPlot({
          graphics::par(mfrow = c(length(covariates), 1))
          plot(res.std)
        })
        output$summ_mm_0 <- DT::renderDataTable({
          tibble::as_tibble(
            cox_reg.mul0 %>% 
              gtsummary::tbl_regression(exp = TRUE)
            )
        })

        output$summ_mm_1 <- DT::renderDataTable({
          tibble::as_tibble(
            cox_reg.mul1 %>% 
              gtsummary::tbl_regression(exp = TRUE)
            )
        })
        output$prop_h_mm0 <- DT::renderDataTable({
          DT::datatable(tibble::rownames_to_column(
            round_df(as.data.frame(res.std_mul0$table), 3), " "),
            options = list(dom = 't')
            )
        })
        output$residues_std_mm0 <- renderPlot({
          graphics::par(mfrow = c(length(covariates), 1))
          plot(res.std_mul0)
        })
        output$prop_h_mm1 <- DT::renderDataTable({
          DT::datatable(tibble::rownames_to_column(
            round_df(as.data.frame(res.std_mul1$table), 3), " "),
            options = list(dom = 't')
            )
        })
        output$residues_std_mm1 <- renderPlot({
          graphics::par(mfrow = c(length(covariates), 1))
          plot(res.std_mul1)
        })
      }
    })
  }
  # Run the application
  # We launch the app in a browser to optimize the user experience
  shiny::shinyApp(
    ui=ui,
    server=server,
    options=list(launch.browser = TRUE),
    )
}


############# Helper functions ###################

#' `get_disease_gene_list` parses an existing .txt
#' file containing all known target genes associated
#' with the selected disease.
#' @param input_disease: input disease selection
#' @return  genes_list
get_disease_gene_list <- function(
    disease_gene,
    input_disease
    ) {
  genes_list <- stats::na.omit(disease_gene[, input_disease]) %>%
    dplyr::rename(genes = input_disease)
  return(genes_list$genes)
}


#' `checkInput`
#' @param input A list containing input values.
#' @return A logical value indicating whether the input is valid.
checkInput <- function(
    input
    ) {
  return (
    input$ids != "N/A"
    & input$time != "N/A"
    & input$group != "N/A"
    & input$event != "N/A"
    & input$target_gene != "N/A"
  )
}

#' `getID`
#' @param x
#' @param geneIDS
#' @return 
getID <- function(
    x,
    geneIDS
    ) {
  gene_name <- geneIDS[x, ]$GeneName
  if (is.na(gene_name)) {
    print(paste(
      'Please include the gene ID ',
      gene_name,
      'in the ensemblTogenes.csvfile'
    ))
  }
  return(gene_name)
}


getGeneName <- function(
    info, 
    geneIDS
    ) {
  x <- stringr::str_extract(
    info['INFO'], 
    "(?<=ensembl_gene_id=)[^;]+"
    )
  if (grepl(",", x,  fixed = TRUE)) {
    grep_id <- base::strsplit(x, split = ",")[[1]]
    return(sapply(
      grep_id,
      getID,
      geneIDS,
      USE.NAMES=FALSE)
      )
  }
  return(getID(x, geneIDS))
}


get_cov_list <- function(
    metadata,
    input
    ) {
  inputs_list <- c(input$event, input$time, input$ids)
  cov_list <- colnames(metadata)
  cov_list <- cov_list[!(cov_list %in% inputs_list)]
  return(cov_list)
}

genesCountTable <- function(
    vcf,
    metadata,
    input,
    gene_ids_table
    ) {
  genes_with_svs_in_sample <- apply(vcf@fix, 1, getGeneName, gene_ids_table)
  sample_names <- colnames(vcf@gt)[-1]
  disease_genes_names <- gene_ids_table$GeneName
  count_df <- CountSVsDf( 
    length(disease_genes_names),
    length(sample_names),
    disease_genes_names,
    sample_names,
    genes_with_svs_in_sample,
    vcf,
    input$ids
    )
  new_md <- Remove_NA(metadata, input$time) %>% 
    dplyr::rename(ids = input$ids)
  # what happen if the input$target_gene is not in the vcf file?
  new_df <- subset(count_df, ids %in% new_md$ids)
  count_table <- as.data.frame(apply(new_df[, -1], 2, function(c) sum(c != 0)))
  colnames(count_table) <- "Count of patients with SVs"
  count_table <- tibble::rownames_to_column(count_table, "Gene_ID")
  count_table <- count_table %>% 
    dplyr::arrange(
    desc(count_table["Count of patients with SVs"])
    )
  return(count_table)
  }

#' `CountSVsDf`
#' @param ncol
#' @param nrow
#' @param disease_genes_names
#' @param sample_names
#' @param genes_with_svs_in_sample
#' @param vcf
#' @param ids_col
#' @return
CountSVsDf <- function(
    ncol,
    nrow,
    disease_genes_names,
    sample_names,
    genes_with_svs_in_sample,
    vcf,
    ids_col
    ) {
  sample_disease_gene_df <- data.frame(
    matrix(0, ncol = ncol, nrow = nrow)
    )
  colnames(sample_disease_gene_df) <- disease_genes_names
  rownames(sample_disease_gene_df) <- sample_names
  for (sv_idx in 1:length(genes_with_svs_in_sample)) {
    # if the sv is in a gene of interest
    for (individual_idx in 1:length(sample_names)) {
      indiv_id <- sample_names[individual_idx]
      gt <- vcf@gt[sv_idx, indiv_id]
      if (!is.na(gt))
      {
        sv_gene_name <- genes_with_svs_in_sample[sv_idx]
        for (gene_id in unlist(sv_gene_name)) {
          sample_disease_gene_df[indiv_id, gene_id] %+=% 1
        }
      }
    }
  }
  sample_disease_gene_df <- tibble::rownames_to_column(
    sample_disease_gene_df,
    ids_col
    )
  sample_disease_gene_df <- sample_disease_gene_df %>% 
    dplyr::rename(ids = ids_col)
  return(sample_disease_gene_df)
}


#' `Remove_NA` Removes rows from a data frame where the specified
#'  time column has missing values or values equal to "NA".
#' @param df: A data frame containing the dataset.
#' @param time_col: A character string specifying the name 
#' of the time column to check.
#' @return A modified data frame with rows removed
#'  based on the specified time column.
Remove_NA <- function(df, time_col) {
  entries <- df[[time_col]]
  nas_entries <- is.na(entries)
  NAs_entries <- entries == "NA"
  if (sum(nas_entries) > 0) {
    return(df[!nas_entries, ])
  } else if (sum(NAs_entries) > 0) {
    return(df[!NAs_entries, ])
  } else {
    return(df)
  }
}

#' `round_df` Round Numeric Columns in a Data Frame
#' @param x A data frame to be modified.
#' @param digits The number of decimal places to round the 
#' numeric columns to.
round_df <- function(x, digits) {
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <- round(x[numeric_columns], digits)
  return(x)
}


map_col_names <- function(input, cov_list) {
  if (input$event %in% cov_list) {
    cov_list <- str_replace(cov_list, input$event, "event")
  }
  if (input$group %in% cov_list) {
    cov_list <- str_replace(cov_list, input$group, "trial_group_bin")
  }
  return(cov_list)
}

`%+=%` <- function(e1, e2)
  eval.parent(substitute(e1 <- e1 + e2))

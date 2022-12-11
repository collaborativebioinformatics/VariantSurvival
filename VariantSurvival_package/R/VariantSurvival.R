#' VariantSurvival
#'
#' @param vcffile path to the vcf file containing the Structural variant data
#' @param metadatafile path to the txt file containing the samples metadata
#'
#' @return
#' @export
#'

source("global.R") # install/load requirements
source("fun.R") # main functions

VariantSurvival <- function(vcffile, metadatafile){
  # parse inputs
  vcf <- vcfR::read.vcfR(vcffile, verbose = FALSE)
  metadata <- readxl::read_excel(metadatafile)
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
                    label = "Gene of interest",
                    choices = NULL,
                    selected = FALSE
        ),
        selectInput(inputId = "time",
                    label = "Select the time factor:",
                    choices = colnames(metadata),
                    selected = FALSE
        ),
        br(),
        selectInput(inputId = "phenotype",
                    label = "Select the phenotype factor:",
                    choices = colnames(metadata),
                    selected = FALSE
        )
      ),
      dashboardBody(
        fluidRow(
          box(width = 6, #height = 600,
              title = "Genes",
              "Based on literature The following genes are associated
              with the disease mechanism",
              br(), #linebreak
              DT::dataTableOutput("geneslist")
          ),
          box(width = 6, #height = 600,
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
  
  # Define server logic required to draw a histogram
  server <- function(input, output, session) {
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
    output$geneslist <- DT::renderDataTable(
      datatable(reactive_gene_list(),
                selection = 'none')
    )
    
    gene_ids_table <- read.csv(file = 'ensembleTogenes.csv')
    rownames(gene_ids_table) <- gene_ids_table$ensembleID
    gene_names = gene_ids_table$GeneName
    sample_names = colnames(vcf@gt)[-1] # VCF genotype information
    
    getGeneName <- function(info) {
      x <- str_extract(info['INFO'], "(?<=ensembl_gene_id=)[^;]+")
      return(geneIDS[x,]$GeneName)
    }
    
    sv_gene <- apply(vcf@fix,1,getGeneName)
    #     allgenes_sv <- data.frame(matrix(0,
    #                                      ncol = length(genes),
    #                                      nrow = length(samples)
    #                                      )
    #                               )
    #     colnames(allgenes_sv) = genes
    #     rownames(allgenes_sv) = samples
    #
    #     for (i in 1:length(sv_gene)) {
    #       if(is.na(sv_gene[i]))
    #         next;
    #       for(j in 1:length(samples)){
    #         gt=vcf@gt[i,samples[j]]
    #         if(!is.na(gt))
    #         {
    #           allgenes_sv[samples[j], sv_gene[[i]]] = allgenes_sv[samples[j],
    #                                                               sv_gene[[i]] ] + 1
    #         }
    #       }
    #     }
    #     allgenes_sv<- rownames_to_column(allgenes_sv, "patient_ID")
    #
    #     metadata2 = data.frame(metadata)
    #     #test
    #     metadata3 = metadata2[c(1,2)]
    #     #merge with metadata
    #     dx <- merge(allgenes_sv, metadata3, by=0, all=TRUE)
    #
    #     dx2 <- dx[-(1)]
    #     dx2 <- dx2[-(22)]
    #     #
    #     output$barplot <- renderPlot({
    #
    #       dx3 <- as.data.frame(
    #         c(dx2["patient_ID.x"],
    #           dx2[input$targetGene],
    #           dx2[input$phenotype]
    #           )
    #         )
    #
    #       colnames(dx3) <- c('patient_ID', 'gene', 'Phenotype')
    #       dx3 <- dx3 %>%
    #         mutate(Phenotype= ifelse(Phenotype=="0", "Placebo","treatment" ))
    #       ggplot(data=dx3,
    #              aes(x=patient_ID,
    #                  y=gene,
    #                  fill=Phenotype)
    #              ) + geom_bar(stat="identity") + theme_classic()
    #       }
    #       )
    #
    #     #reactive output
    #     output$svtable <-DT::renderDataTable({
    #       gene_sv <- as.data.frame(
    #         c(allgenes_sv["patient_ID"],
    #           allgenes_sv[input$targetGene]
    #           )
    #         )
    #       }
    #       )
    #     #step 1 kaplan-meier
    #     output$plot1 <- renderPlot({
    #       gene_sv2 <- as.data.frame(
    #         c(allgenes_sv["patient_ID"],
    #           allgenes_sv[input$targetGene]
    #           )
    #         )
    #       colnames(gene_sv2) <- c('patient_ID',
    #                               'Structural_Variants_count')
    #       gene_sv2$variant <- "No"
    #       gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
    #       de2 <- merge(gene_sv2,
    #                    metadata2,
    #                    by=0,
    #                    all=TRUE)
    #       de2 <- de2[-(1)]
    #       de2 <- de2[-(4)]
    #
    #       df <- de2 %>%
    #         mutate(variant= ifelse(variant=="Yes", 1, 0))
    #       df <- df %>% rename(Time = input$time)
    #       df <- df %>% rename(Phenotype = input$phenotype)
    #       #
    #       survfit2(Surv(Time,variant)~ 1, data = df) %>%
    #         ggsurvfit() +
    #         labs(
    #           x = "Days",
    #           y = "Overall survival probability"
    #         ) +
    #         add_confidence_interval() +
    #         add_risktable()
    #       }
    #       )
    #
    #     # x-year survival time
    #     output$table1 <- DT::renderDataTable({
    #       x<- survfit(Surv(Time,variant)~ 1, data = df) %>%
    #         tbl_survfit(
    #           times = 365.25,
    #           label_header = "**1-year survival (95% CI)**"
    #           )
    #       t <- as_data_frame(x)
    #       t
    #       }
    #       )
    #
    #     # Median survival time
    #     output$table2 <- DT::renderDataTable({
    #       x2 <-  survfit(Surv(Time,variant)~ 1, data = df) %>%
    #         gtsummary::tbl_survfit(
    #           probs = 0.5,
    #           label_header = "**Median survival (95% CI)**"
    #         )
    #       t2 <- as_data_frame(x2)
    #       t2
    #     })
    #
    #     # step2 survival curve
    #     #plot2
    #     output$plot2 <- renderPlot({
    #       gene_sv2 <- as.data.frame(
    #         c(allgenes_sv["patient_ID"],
    #           allgenes_sv[input$targetGene]
    #           )
    #         )
    #       colnames(gene_sv2) <- c('patient_ID',
    #                               'Structural_Variants_count')
    #       gene_sv2$variant <- "No"
    #       gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
    #       de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)
    #       de2 <- de2[-(1)]
    #       de2 <- de2[-(4)]
    #       df <- de2 %>%
    #         mutate(variant= ifelse(variant=="Yes", 1, 0))
    #       #input factors
    #       #df <- df %>% rename(Time = Time_to_death_or_last_followup_days)
    #       df <- df %>% rename(Time = input$time)
    #       df <- df %>% rename(Phenotype = input$phenotype)
    #       #
    #       s <- survfit(Surv(Time,variant)~ Phenotype, data = df)
    #       g <- ggsurvplot(
    #         s,
    #         conf.int = TRUE,
    #         data = df,
    #         risk.table = TRUE
    #       )
    #       g
    #
    #     })
    #
    #
    #     #plot3
    #     output$plot3 <- renderPlot({
    #       gene_sv2 <- as.data.frame(
    #         c(allgenes_sv["patient_ID"],
    #           allgenes_sv[input$targetGene]
    #           )
    #         )
    #       colnames(gene_sv2) <- c('patient_ID',
    #                               'Structural_Variants_count')
    #       gene_sv2$variant <- "No"
    #       gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
    #       de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)
    #       de2 <- de2[-(1)]
    #       de2 <- de2[-(4)]
    #       df <- de2 %>%
    #         mutate(variant= ifelse(variant=="Yes", 1, 0))
    #       #input factors
    #       #• df <- df %>% rename(Time = Time_to_death_or_last_followup_days)
    #       df <- df %>% rename(Time = input$time)
    #
    #       #df <- df %>% rename(Phenotype = "Phenotype")
    #       df <- df %>% rename(Phenotype = input$phenotype)
    #       #
    #       s <- survfit(Surv(Time,variant)~ Phenotype, data = df)
    #       g2 <- ggsurvplot_facet(
    #         s,
    #         conf.int = TRUE,
    #         data = df,
    #         facet.by = "Structural_Variants_count",
    #       )
    #       g2
    #
    #     })
    #
    #     #regression table
    #
    #     output$table3 <- DT::renderDataTable({
    #       x3 <- coxph(Surv(Time,variant)~ Phenotype, data = df) %>%
    #         tbl_regression(exp = TRUE)
    #
    #       t3 <- as_data_frame(x3)
    #
    #       t3
    #     })
    # #d plot1
    #     output$download1plot <- downloadHandler(
    #       filename = function(){
    #         paste("kaplan_meier_plot", "png", sep=".")
    #       },
    #       content = function(file){
    #         png(file)
    #         gene_sv2 <- as.data.frame(
    #           c(allgenes_sv["patient_ID"],
    #             allgenes_sv[input$targetGene]
    #             )
    #           )
    #         colnames(gene_sv2) <- c('patient_ID',
    #                                 'Structural_Variants_count')
    #         gene_sv2$variant <- "No"
    #         gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
    #         de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)
    #         de2 <- de2[-(1)]
    #         de2 <- de2[-(4)]
    #         df <- de2 %>%
    #           mutate(variant= ifelse(variant=="Yes", 1, 0))
    #         #input factors
    #         #df <- df %>% rename(Time = Time_to_death_or_last_followup_days)
    #         df <- df %>% rename(Time = input$time)
    #         df <- df %>% rename(Phenotype = input$phenotype)
    #         #
    #         g1 <-survfit2(Surv(Time,variant)~ 1, data = df) %>%
    #           ggsurvfit() +
    #           labs(
    #             x = "Days",
    #             y = "Overall survival probability"
    #           ) +
    #           add_confidence_interval() +
    #           add_risktable()
    #         g1
    #         dev.off()
    #       })
    #     #d plot2
    #     output$download2plot <- downloadHandler(
    #       filename = function(){
    #         paste("survival_plot", "png", sep=".")
    #       },
    #       content = function(file){
    #         png(file)
    #         gene_sv2 <- as.data.frame(
    #           c(allgenes_sv["patient_ID"],
    #             allgenes_sv[input$targetGene]
    #             )
    #           )
    #         colnames(gene_sv2) <- c('patient_ID',
    #                                 'Structural_Variants_count')
    #         gene_sv2$variant <- "No"
    #         gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
    #         de2 <- merge(gene_sv2,
    #                      metadata2,
    #                      by=0,
    #                      all=TRUE)
    #         de2 <- de2[-(1)]
    #         de2 <- de2[-(4)]
    #
    #         df <- de2 %>%
    #           mutate(variant= ifelse(variant=="Yes", 1, 0))
    #         #input factors
    #         #df <- df %>% rename(Time = Time_to_death_or_last_followup_days)
    #         df <- df %>% rename(Time = input$time)
    #         df <- df %>% rename(Phenotype = input$phenotype)
    #         #
    #         s <- survfit(Surv(Time,variant)~ Phenotype, data = df)
    #         g <- ggsurvplot(
    #           s,
    #           conf.int = TRUE,
    #           data = df,
    #           risk.table = TRUE
    #         )
    #         g
    #         dev.off()
    #       })
    #     #d plot
    #     output$download3plot <- downloadHandler(
    #       filename = function(){
    #         paste("survival_plot2", "png", sep=".")
    #       },
    #       content = function(file){
    #         png(file)
    #         gene_sv2 <- as.data.frame(
    #           c(allgenes_sv["patient_ID"],
    #             allgenes_sv[input$targetGene]
    #             )
    #           )
    #         colnames(gene_sv2) <- c('patient_ID',
    #                                 'Structural_Variants_count')
    #         gene_sv2$variant <- "No"
    #         gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
    #         de2 <- merge(gene_sv2,
    #                      metadata2,
    #                      by=0,
    #                      all=TRUE)
    #         de2 <- de2[-(1)]
    #         de2 <- de2[-(4)]
    #         df <- de2 %>%
    #           mutate(variant= ifelse(variant=="Yes", 1, 0))
    #         df <- df %>% rename(Time = input$time)
    #         df <- df %>% rename(Phenotype = input$phenotype)
    #         #
    #         s <- survfit(Surv(Time,variant)~ Phenotype, data = df)
    #         g2 <-ggsurvplot_facet(
    #           s,
    #           conf.int = TRUE,
    #           data = df,
    #           facet.by = "Structural_Variants_count")
    #         g2
    #         dev.off()
    #       })
  }
  
  # Run the application
  shinyApp(ui = ui, server = server)
}


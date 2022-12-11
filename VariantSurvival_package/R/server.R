source("fun.R") # main functions


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
  ## output ##
  output$geneslist <- DT::renderDataTable(
    datatable(reactive_gene_list(),
              selection = 'none')
  )
  rownames(gene_ids_table) <- gene_ids_table$ensembleID
  disease_genes_names = gene_ids_table$GeneName
  sample_names = colnames(vcf@gt)[-1] # VCF genotype information
  
  # genes are repeated since a single gene can have more than one SV.
  genes_with_svs_in_sample <- apply(vcf@fix,1,getGeneName)
  
  # count dataframe with patient_ids in rows and gene_ids in columns
  sample_disease_gene_df <- CountSVsDf(length(disease_genes_names),
                                       length(sample_names),
                                       disease_genes_names,
                                       sample_names,
                                       genes_with_svs_in_sample,
                                       vcf)
  
  #incorporate phenotype column (placebo/treatment)
  new_sample_disease_gene_df <- merge(sample_disease_gene_df, 
                                      metadata[, 1:2],
                                      on = "patient_ID")
  
  ## output ##
  output$barplot <- renderPlot(
    {
    svs_gene_input_df <- as.data.frame(
      c(new_sample_disease_gene_df["patient_ID"] ,
        new_sample_disease_gene_df[input$target_gene],
        new_sample_disease_gene_df[input$phenotype]
        )
      )
    colnames(svs_gene_input_df) <- c('patient_ID',
                                     'SVs_number_per_gene',
                                     'Phenotype')

    svs_gene_input_df <- svs_gene_input_df %>%
      mutate(Phenotype = ifelse(Phenotype=="0",
                                "Placebo", 
                                "Treatment" )
             )
    ggplot(svs_gene_input_df,
           aes(SVs_number_per_gene, 
               fill=Phenotype)) +
      geom_histogram(binwidth=1) +
      stat_bin(binwidth=1,
               geom='text',
               color='white',
               aes(label=after_stat(count)),
               position=position_stack(vjust = 0.5)) + 
      xlab("Number of SVs in target gene") + ylab("Frequency")
    }
    )
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

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
  # output$geneslist <- DT::renderDataTable(
  #   datatable(reactive_gene_list(),
  #             selection = 'none')
  # )
  
  disease_genes_names = gene_ids_table$GeneName
  sample_names = colnames(vcf@gt)[-1] # VCF genotype information
  
  # genes are repeated since a single gene can have more than one SV.
  genes_with_svs_in_sample <- apply(vcf@fix,1,getGeneName)
  
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
  output$barplot <- renderPlot(
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
        xlab("Number of SVs in target gene") + ylab("Frequency")
      }
    )
  
  output$sc_no_sv <- renderPlot(
    {
      svs_gene_input_df <- reactive_no_NAs_metadata()
      # subset those patients with no sv
      no_sv <- svs_gene_input_df[svs_gene_input_df$SV_binary == 0,]
      survfit2(Surv(time, event)~ Phenotype, data = no_sv) %>%
      ggsurvfit() + labs(x = "Days", y = "Overall survival probability") +
        add_confidence_interval() + add_risktable()
      }
    )
    
  # survival_ <- apply(vcf@fix,1,getGeneName)
  #       
  # ## output - Kaplan - Meier ##
  # output$plot1 <- renderPlot({
  #   # 0 if there's no variant in the target gene , 1 otherwise
  #   svs_gene_input_df <- svs_gene_input_df %>%
  #     mutate(variant= ifelse(svs_gene_input_df$SVs_count_per_gene > 0,
  #                            "Has Variant", 
  #                            "No Variant"))
  #   
  #   svs_gene_input_df <- merge(svs_gene_input_df,
  #                              metadata[,c("patient_ID",input$time)],
  #                              on = "patient_ID")
  #   # rename columns
  #   svs_gene_input_df <- svs_gene_input_df %>%
  #     rename(time = input$time) %>%
  #     rename(phenotype = input$phenotype)
  # #   
  #   survfit2(Surv(time, new_list)~ variant, data = svs_gene_input_df) %>%
  #     ggsurvfit() +
  #     labs(
  #       x = "Days",
  #       y = "Overall survival probability"
  #     ) +
  #     add_confidence_interval() +
  #     add_risktable()
  #   
  # 
  #   survfit2(Surv(time, new_list)~ new_sample_disease_gene_df$Phenotype, 
  #            data = svs_gene_input_df) %>%
  #     ggsurvfit() +
  #     labs(
  #       x = "Days",
  #       y = "Overall survival probability"
  #     ) +
  #     add_confidence_interval() +
  #     add_risktable()
  # })
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

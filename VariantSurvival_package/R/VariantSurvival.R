#' VariantSurvival
#'
#' @param vcffile path to the vcf file containing the Structural variant data
#' @param metadatafile path to the txt file containing the samples metadata
#'
#' @return
#' @export
#'


VariantSurvival <- function(vcffile, metadatafile){
  if (!require("shiny")) install.packages("shiny")
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

  vcf <- vcfR::read.vcfR(vcffile, verbose = FALSE)
  metadata <- readxl::read_excel(metadatafile)


  ui <- fluidPage(
    dashboardPage(
      dashboardHeader(title = "VariantSurvival"),
      dashboardSidebar(
        selectizeInput("disease", "Diseases", choices = c("Amyotrophic lateral sclerosis"=1,"Parkinson's disease"=2,"Alzheimer's disease"=3,
                                                          "Friedreich ataxia"=4, "Huntington's disease"=5, "Lewy body disease"=6,
                                                          "Spinal muscular atrophy"=7, selected = 1)),
        textInput("targetGene", label = "Gene of interest"
        ),
        selectInput("time", label = "Select the time factor:",
                    choices = colnames(metadata)),
        br(),

        selectInput("phenotype", label = "Select the phenotype factor:",
                    choices = colnames(metadata)
        )
      ),

      dashboardBody(
        fluidRow(
          box(width = 6, #height = 600,
              title = "Genes",
              "Based on literature The following genes are associated with the disease mechanism",
              br(),
              DT::dataTableOutput("geneslist")
              #renderDataTable("geneslist")
          ),
          box(width = 6, #height = 600,
              title = "Structural Variant in selected gene",
              br(),

              plotOutput("barplot")
          )),

        fluidRow(

          tabBox(width = 6,# height = 600,
                 tabPanel(title = "Survival Plot",
                          h6("Phenotype 0 = Placebo ; 1= Treatment"),
                          plotOutput("plot2"),downloadButton("download2plot", "Download as PNG")
                 ),

                 tabPanel(title = "Survival Plot according to SV count",                     
                          h6("starta 0 = Placebo ; 1= Treatment"),
                          plotOutput("plot3"),downloadButton("download3plot", "Download as PNG")
                 ),
                 tabPanel(
                   title = "Cox regression model",
                   DT::dataTableOutput("table3")
                   ,"HR = Hazard Ratio, CI = Confidence Interval")

          ),

          tabBox(width = 6, #height = 700,
                 tabPanel(title = "Kaplan-Meier plot",

                          plotOutput("plot1"),downloadButton("download1plot", "Download as PNG")
                 ),
                 tabPanel(
                   title = "1-year survival time time",
                   DT::dataTableOutput("table1")
                   ,"CI = Confidence Interval"),
                 tabPanel(
                   title = "median survival time",
                   DT::dataTableOutput("table2")
                   ,"CI = Confidence Interval"))

        )
      )))



  # Define server logic required to draw a histogram
  server <- function(input, output) {

    output$geneslist <- DT::renderDataTable({

      if(as.numeric(input$disease == 1)){
        geneslist  <-read_csv("disease_gene/ALS/genes_list.txt")
      } else if(as.numeric(input$disease == 2)){
        geneslist  <-read_csv("disease_gene/PD/genes_list.txt")
      } else if(as.numeric(input$disease == 3)) {
        geneslist  <-read_csv("disease_gene/AD/genes_list.txt")
      } else if(as.numeric(input$disease == 4)) {
        geneslist  <-read_csv("disease_gene/FD/genes_list.txt")
      } else {
        geneslist  <-read_csv("disease_gene/DLB/genes_list.txt")
      }

    })

    geneIDS<- read.csv(file = 'ensembleTogenes.csv')
    rownames(geneIDS) <- geneIDS$ensembleID


    samples=colnames(vcf@gt)
    samples=samples[2:length(samples)]

    genes=geneIDS$GeneName


    getGeneName <- function(info) {
      s <- str_extract(info["INFO"], "ensembl_gene_id=[^;]*")
      s2=str_split(s,'=')[[1]][2]
      s3=as.array(str_split(s2,',')[[1]])
      apply(s3,1,function(x) {geneIDS[x,]$GeneName} )
    }


    sv_gene <- apply(vcf@fix,1,getGeneName)


    allgenes_sv <- data.frame(matrix(0,ncol = length(genes), nrow = length(samples)))
    colnames(allgenes_sv) = genes
    rownames(allgenes_sv) = samples

    for (i in 1:length(sv_gene)) {
      if(is.na(sv_gene[i]))
        next;
      for(j in 1:length(samples)){
        gt=vcf@gt[i,samples[j]]
        if(!is.na(gt))
        {
          allgenes_sv[samples[j], sv_gene[[i]] ]=allgenes_sv[samples[j], sv_gene[[i]] ]+1
        }
      }
    }
    allgenes_sv<- rownames_to_column(allgenes_sv, "patient_ID")


    metadata2 = data.frame(metadata)
    #test
    metadata3 = metadata2[c(1,2)]
    #merge with metadata
    dx <- merge(allgenes_sv, metadata3, by=0, all=TRUE)

    dx2 <- dx[-(1)]
    dx2 <- dx2[-(22)]
    #
    output$barplot <- renderPlot({

      dx3 <- as.data.frame(c(dx2["patient_ID.x"] ,dx2[input$targetGene], dx2[input$phenotype]))

      #dx3 <- as.data.frame(c(dx2["patient_ID.x"] ,dx2["SETX"], dx2["Phenotype"]))

      colnames(dx3) <- c('patient_ID','gene', 'Phenotype')
      dx3 <- dx3 %>%
        mutate(Phenotype= ifelse(Phenotype=="0", "Placebo","treatment" ))
      ggplot(data=dx3, aes(x=patient_ID, y=gene,fill=Phenotype)) +
        labs(y = "Structural variant count", x = "Sample")+
        geom_bar(stat="identity")+
        theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())

    })

    #reactive output
    output$svtable <-DT::renderDataTable({

      gene_sv <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))

    })

    #step 1 kaplan-meier
    output$plot1 <- renderPlot({
      gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))
      colnames(gene_sv2) <- c('patient_ID','Structural_Variants_count')
      gene_sv2$variant <- "No"
      gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
      de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)
      de2 <- de2[-(1)]
      de2 <- de2[-(4)]
      df <- de2 %>%
        mutate(variant= ifelse(variant=="Yes", 1, 0))
      #input factors
      #df <- df %>% rename(Time = Time_to_death_or_last_followup_days)
      df <- df %>% rename(Time = input$time)
      df <- df %>% rename(Phenotype = input$phenotype)
      #
      survfit2(Surv(Time,variant)~ 1, data = df) %>%
        ggsurvfit() +
        labs(
          x = "Days",
          y = "Overall survival probability"
        ) +
        add_confidence_interval() +
        add_risktable()

    })


    # x-year survival time

    output$table1 <- DT::renderDataTable({
      x<- survfit(Surv(Time,variant)~ 1, data = df) %>%
        tbl_survfit(
          times = 365.25,
          label_header = "**1-year survival (95% CI)**"
        )
      t <- as_tibble(x)

    })

    # Median survival time
    output$table2 <- DT::renderDataTable({
      x2 <-  survfit(Surv(Time,variant)~ 1, data = df) %>%
        gtsummary::tbl_survfit(
          probs = 0.5,
          label_header = "**Median survival (95% CI)**"
        )
      t2 <- as_tibble(x2)
 
    })

    # step2 survival curve
    #plot2
    output$plot2 <- renderPlot({
      gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))
      colnames(gene_sv2) <- c('patient_ID','Structural_Variants_count')
      gene_sv2$variant <- "No"
      gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
      de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)
      de2 <- de2[-(1)]
      de2 <- de2[-(4)]
      df <- de2 %>%
        mutate(variant= ifelse(variant=="Yes", 1, 0))
      #input factors
      #df <- df %>% rename(Time = Time_to_death_or_last_followup_days)
      df <- df %>% rename(Time = input$time)
      df <- df %>% rename(Phenotype = input$phenotype)
      #
      s <- survfit(Surv(Time,variant)~ Phenotype, data = df)
      g <- ggsurvplot(
        s,
        conf.int = TRUE,
        data = df,
        risk.table = TRUE
      )
      g

    })


    #plot3
    output$plot3 <- renderPlot({
      gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))
      colnames(gene_sv2) <- c('patient_ID','Structural_Variants_count')
      gene_sv2$variant <- "No"
      gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
      de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)
      de2 <- de2[-(1)]
      de2 <- de2[-(4)]
      df <- de2 %>%
        mutate(variant= ifelse(variant=="Yes", 1, 0))
      #input factors
      #• df <- df %>% rename(Time = Time_to_death_or_last_followup_days)
      df <- df %>% rename(Time = input$time)

      #df <- df %>% rename(Phenotype = "Phenotype")
      df <- df %>% rename(Phenotype = input$phenotype)
      #
      s <- survfit(Surv(Time,variant)~ Phenotype, data = df)
      g2 <- ggsurvplot_facet(
        s,
        conf.int = TRUE,
        data = df,
        facet.by = "Structural_Variants_count",
      )
      g2

    })

    #regression table

    output$table3 <- DT::renderDataTable({
      x3 <- coxph(Surv(Time,variant)~ Phenotype, data = df) %>%
        tbl_regression(exp = TRUE)

      t3 <- as_tibble(x3)

    })
# plot1
    output$download1plot <- downloadHandler(
      filename = function(){
        paste("kaplan_meier_plot", "png", sep=".")
      },
      content = function(file){
        png(file)
        gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))
        colnames(gene_sv2) <- c('patient_ID','Structural_Variants_count')
        gene_sv2$variant <- "No"
        gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
        de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)
        de2 <- de2[-(1)]
        de2 <- de2[-(4)]
        df <- de2 %>%
          mutate(variant= ifelse(variant=="Yes", 1, 0))
        #input factors
        #df <- df %>% rename(Time = Time_to_death_or_last_followup_days)
        df <- df %>% rename(Time = input$time)
        df <- df %>% rename(Phenotype = input$phenotype)
        #
        g1 <-survfit2(Surv(Time,variant)~ 1, data = df) %>%
          ggsurvfit() +
          labs(
            x = "Days",
            y = "Overall survival probability"
          ) +
          add_confidence_interval() +
          add_risktable()
        g1
        dev.off()
      })
    #d plot2
    output$download2plot <- downloadHandler(
      filename = function(){
        paste("survival_plot", "png", sep=".")
      },
      content = function(file){
        png(file)
        gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))
        colnames(gene_sv2) <- c('patient_ID','Structural_Variants_count')
        gene_sv2$variant <- "No"
        gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
        de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)
        de2 <- de2[-(1)]
        de2 <- de2[-(4)]
        df <- de2 %>%
          mutate(variant= ifelse(variant=="Yes", 1, 0))
        #input factors
        #df <- df %>% rename(Time = Time_to_death_or_last_followup_days)
        df <- df %>% rename(Time = input$time)
        df <- df %>% rename(Phenotype = input$phenotype)
        #
        s <- survfit(Surv(Time,variant)~ Phenotype, data = df)
        g <- ggsurvplot(
          s,
          conf.int = TRUE,
          data = df,
          risk.table = TRUE
        )
        g
        dev.off()
      })
    #d plot
    output$download3plot <- downloadHandler(
      filename = function(){
        paste("survival_plot2", "png", sep=".")
      },
      content = function(file){
        png(file)
        gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))
        colnames(gene_sv2) <- c('patient_ID','Structural_Variants_count')
        gene_sv2$variant <- "No"
        gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
        de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)
        de2 <- de2[-(1)]
        de2 <- de2[-(4)]
        df <- de2 %>%
          mutate(variant= ifelse(variant=="Yes", 1, 0))
        #input factors
        #• df <- df %>% rename(Time = Time_to_death_or_last_followup_days)
        df <- df %>% rename(Time = input$time)

        #df <- df %>% rename(Phenotype = "Phenotype")
        df <- df %>% rename(Phenotype = input$phenotype)
        #
        s <- survfit(Surv(Time,variant)~ Phenotype, data = df)
        g2 <-ggsurvplot_facet(
          s,
          conf.int = TRUE,
          data = df,
          facet.by = "Structural_Variants_count")
        g2
        dev.off()
      })
  }

  # Run the application
  shinyApp(ui = ui, server = server)
}

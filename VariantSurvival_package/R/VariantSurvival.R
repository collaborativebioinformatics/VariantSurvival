#' Title
#'
#' @param vcffile path to the vcf file containing the Structural variant data
#' @param metadatafile path to the txt file containing the samples metadata
#'
#' @return
#' @export
#'
#' @examples
VariantSurvival <- function(vcffile, metadatafile){
  #vcffile <- "merged.filtered.vcf"
  #metadatafile <- "metadata.xlsx"
  vcf <- read.vcfR(vcffile, verbose = FALSE)
  metadata <- read_excel(metadatafile)

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

  ui <- fluidPage(
    dashboardPage(
      dashboardHeader(title = "VariantSurvival"),
      dashboardSidebar(
        selectizeInput("disease", "Diseases", choices = c("Amyotrophic lateral sclerosis"=1,"Parkinson's disease"=2,"Alzheimer's disease"=3,
                                                          "Friedreich ataxia"=4, "Huntington's disease"=5, "Lewy body disease"=6,
                                                          "Spinal muscular atrophy"=7, selected = 1)),
        textInput("targetGene", label = h3("Gene of interest")#, value ="SETX"
        )
      ),

      dashboardBody(
        fluidRow(
          box(width = 6, height = 600,
              title = "Genes",
              "Based on literature The following genes are associated with the disease mechanism",
              br(),
              DT::dataTableOutput("geneslist"),
          ),
          box(width = 6, height = 600,
              title = "Structural Variant in selected gene",
              br(),
              br(),
              br(),
              br(),
              plotOutput("barplot")
          )

        ),
        #part1
        fluidRow(box(title = "about the survival plot",
                     width = 4, height = 400,
                     "The survival analysis is perfomed as follows: survfit(Surv(input$factor1, Sructural variant Count)~ input$factor2,, data = df)",
                     "Defualt setting are",
                     br(),
                     "Factor 1 = Time to death or last followup",
                     br(),
                     "Factor 2 = Phenotype"),

                 box(width = 4 , height = 400,
                     "Please define the factors column name in the metadata file.",
                     br(),
                     textInput("factor1", label = h3("Select the first factor.")),
                     br(),
                     textInput("factor2", label = h3("Select the second factor.")),

                 ),
                 box(width = 4, height = 400,
                     title = "Competing risks regression",
                     DT::dataTableOutput("table")
                     ,"HR = Hazard Ratio, CI = Confidence Interval")
        ),
        fluidRow(
          box(title = "Survival Plot",
              width = 6, height = 720,
              plotOutput("plot4")
          ),

          box(title = "Survival Plot according to SV count",
              "set the second factor as the structural variant count, to visualize the count effect.",
              width = 6, height = 720,
              plotOutput("plot5"),
          )


        )

      )
    ))

  server <- function(input, output) {

    output$geneslist <- DT::renderDataTable({
      if(as.numeric(input$disease == 1)){
        geneslist  <-read_csv("disease_gene/ALS/genes_list.txt")
      } else if(as.numeric(input$disease == 2)){
        geneslist  <-read_csv("disease_gene/PD/genes_list.txt")
      } else if(as.numeric(input$disease == 3)) {
        geneslist  <-read_csv("disease_gene/AD/genes_list.txt")
      } else if(as.numeric(input$disease == 4)) {
        geneslist  <-read_csv("disease_gene/FA/genes_list.txt")
      } else if(as.numeric(input$disease == 5)) {
        geneslist  <-read_csv("disease_gene/HD/genes_list.txt")
      } else if(as.numeric(input$disease == 6)) {
        geneslist  <-read_csv("disease_gene/LD/genes_list.txt")
      } else {
        geneslist  <-read_csv("disease_gene/SMA/genes_list.txt")
      }
    })
    #step1

    geneIDS<- read.csv(file = 'ensembleTogenes.csv')
    rownames(geneIDS) <- geneIDS$ensembleID


    samples=colnames(vcf@gt)
    samples=samples[2:length(samples)]

    genes=geneIDS$GeneName


    getGeneName <- function(info) {
      s <- str_extract(info["INFO"], "ensembl_gene_id=[^;]*")
      s2=str_split(s,'=')[[1]][2]
      s3=as.array(str_split(s2,',')[[1]])
      apply(s3,1,(\(x) geneIDS[x,]$GeneName))
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

    #metadata
    metadata2 = data.frame(metadata)
    metadata3 = metadata2[c(1,2)]
    #merge with metadata
    dx <- merge(allgenes_sv, metadata3, by=0, all=TRUE)

    dx2 <- dx[-(1)]
    dx2 <- dx2[-(22)]

    # reactive output
    output$barplot <- renderPlot({

      dx3 <- as.data.frame(c(dx2["patient_ID.x"] ,dx2[input$targetGene], dx2["Phenotype"]))

      #dx3 <- as.data.frame(c(dx2["patient_ID.x"] ,dx2["SETX"], dx2["Phenotype"]))
      colnames(dx3) <- c('patient_ID','gene', 'Phenotype')
      dx3 <- dx3 %>%
        mutate(Phenotype= ifelse(Phenotype=="0", "Placebo","treatment" ))
      ggplot(data=dx3, aes(x=patient_ID, y=gene,fill=Phenotype)) +
        geom_bar(stat="identity")+
        theme_classic()
    })

    #reactive output
    output$svtable <-DT::renderDataTable({

      gene_sv <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))
      #gene_sv <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv["SETX"]))
    })

    output$detable <-DT::renderDataTable({
      gene_sv <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))
      #gene_sv <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv["SETX"]))
      colnames(gene_sv) <- c('patient_ID','Structural_Variants_count')
      gene_sv$variant <- "No"
      gene_sv$variant[gene_sv$Structural_Variants_count > 1 ] <- "Yes"
      de <- merge(gene_sv, metadata2, by=0, all=TRUE)  # merge by row names (by=0 or by="row.names")
      de <- de[-(1)]
      de <- de[-(4)]
    })

    # step2


    output$plot4 <- renderPlot({
      gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))
      #gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv["SETX"])) #remove
      colnames(gene_sv2) <- c('patient_ID','Structural_Variants_count')
      gene_sv2$variant <- "No"
      gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
      de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)  # merge by row names
      de2 <- de2[-(1)]
      de2 <- de2[-(4)]
      df <- de2 %>%
        mutate(variant= ifelse(variant=="Yes", 1, 0))
      #input factors
      df <- df %>% rename(factor1 = input$factor1)
      df <- df %>% rename(factor2 = input$factor2)

      #
      s <- survfit(Surv(factor1, variant)~ factor2, data = df)
      g <- ggsurvplot(
        s,
        conf.int = TRUE,
        data = df,
        risk.table = TRUE
      )
      g


    })

    #regression table
    output$table0 <- DT::renderDataTable({
      gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))
      #gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv["SETX"])) #remove
      colnames(gene_sv2) <- c('patient_ID','Structural_Variants_count')
      gene_sv2$variant <- "No"
      gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
      de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)  # merge by row names
      de2 <- de2[-(1)]
      de2 <- de2[-(4)]
      df <- de2 %>%
        mutate(variant= ifelse(variant=="Yes", 1, 0))
      #input factors
      df <- df %>% rename(factor1 = input$factor1)
      df <- df %>% rename(factor2 = input$factor2)

      #
      s <- survfit(Surv(factor1, variant)~ factor2, data = df)
      x3 <- coxph(Surv(factor1, variant)~ factor2, data = df) %>%
        tbl_regression(exp = TRUE)
      t3 <- as_data_frame(x3)
      t3[t3 == "factor2"] <- input$factor2
    })
    output$table <- DT::renderDataTable({
      t3

    })

    #plot2
    output$plot5 <- renderPlot({
      gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv[input$targetGene]))
      #gene_sv2 <- as.data.frame(c(allgenes_sv["patient_ID"], allgenes_sv["SETX"])) #remove
      colnames(gene_sv2) <- c('patient_ID','Structural_Variants_count')
      gene_sv2$variant <- "No"
      gene_sv2$variant[gene_sv2$Structural_Variants_count > 1 ] <- "Yes"
      de2 <- merge(gene_sv2, metadata2, by=0, all=TRUE)  # merge by row names
      de2 <- de2[-(1)]
      de2 <- de2[-(4)]
      df <- de2 %>%
        mutate(variant= ifelse(variant=="Yes", 1, 0))
      #input factors
      df <- df %>% rename(factor1 = input$factor1)
      df <- df %>% rename(factor2 = input$factor2)
      #
      s <- survfit(Surv(factor1, variant)~ factor2, data = df)
      g2 <- ggsurvplot_facet(
        s,
        conf.int = TRUE,
        data = df,
        facet.by = "Structural_Variants_count",
      )
      g2

    })


  }


  shinyApp(ui = ui, server = server)


}

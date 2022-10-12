
vcf <- read.vcfR("server/merged.filtered.vcf", verbose = FALSE )

geneIDS<- read.csv(file = 'server/ensembleTogenes.csv')
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
#allgenes_sv<- rownames_to_column(allgenes_sv, "patient_ID")


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
#allgenes_sv$patient_id <- rownames(allgenes_sv)
#read metadata
metadata <- read_excel("server/metadata.xlsx")   #input
metadata2 = data.frame(metadata)
#test
metadata3 = metadata2[c(1,2)]
#merge with metadata
dx <- merge(allgenes_sv, metadata3, by=0, all=TRUE) 

dx2 <- dx[-(1)]
dx2 <- dx2[-(22)]

#

# reactive output
output$barplot <- renderPlot({

  dx3 <- as.data.frame(c(dx2["patient_ID.x"] ,dx2[input$targetGene], dx2["Phenotype"]))
  colnames(dx3) <- c('patient_ID','gene', 'Phenotype')
  dx3 <- dx3 %>%
    mutate(Phenotype= ifelse(Phenotype=="0", "Placebo","treatment" ))
  ggplot(data=dx3, aes(x=patient_ID, y=gene,fill=Phenotype)) +
    geom_bar(stat="identity")+
    theme_classic()
})
#test
#8gene_sv <-  as.data.frame(allgenes_sv [c("patient_ID","SETX")])
#reactive output
output$svtable <-DT::renderDataTable({
  gene_sv <- as.data.frame(c( allgenes_sv["patient_ID"] ,allgenes_sv[input$targetGene]))
  gene_sv
}) 

gene_sv2 <- gene_sv 

colnames(gene_sv2) <- c('patient_ID','SV_count')
gene_sv2$variant <- "No"
gene_sv2$variant[gene_sv2$SV_count > 1 ] <- "Yes"

#merge metadata file + variant file
# metadata = input$metafile

#metadata <- read_excel("server/metadata.xlsx")  
#metadata2 = data.frame(metadata)

# variant =gene_sv2
de <- merge(gene_sv2, metadata2, by=0, all=TRUE)  # merge by row names (by=0 or by="row.names")

de2 <- de[-(1)]
de2 <- de2[-(4)]




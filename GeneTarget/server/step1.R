#unzip gff3 folder + read
#gzip -d gencode.v19.annotation.gff3.gz

#gene = input$targetGene
#generate genes regions

#%%bash
#parallel --gnu " ./GeneTarget/getRegions.sh {} > {}.regions  " ::: ALS2 ANG CHMP2B DAO DCTN1 FIG4 FUS NEFH OPTN PFN1 PON1 PON2   PON3 PRPH SETX  SOD1  SQST
#%%bash
#wc -l *regions

# read inputVCF files 
#%%bash
#dx download SV/*

#%%bash

#echo "ALS2 ANG CHMP2B DAO DCTN1 FIG4 FUS NEFH OPTN PFN1 PON1 PON2   PON3 PRPH SETX  SOD1  SQSTM1 TARDBP TREM2 UBQLN2  VAPB VCPVEGFA"| tr -s ' ' $'\n' > genes
#ls GeneTarget/data/ | grep -P "gz$" > samples
#parallel --gnu 'bcftools view -R {1}.regions GeneTarget/data/{2} > events/{1}_{2} ' :::: genes :::: samples


# create dataframe
#reticulate::repl_python()

#source(table.py)




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
rownames(allgenes_sv) = samples
colnames(allgenes_sv) = genes

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

allgenes_sv$patient_id <- rownames(allgenes_sv)

gene_sv <- as.data.frame(c( allgenes_sv["patient_id"] ,allgenes_sv["FIG4"]) )


library(readr)
#sv_events <- read_csv("server/sv_events.csv")
#allgenes_sv <- as.data.frame(sv_events)
#gene_sv <-  as.data.frame(allgenes_sv [c(1,5)])



print(gene_sv) 
gene_sv2 <- gene_sv 
#change name column
# consider variant 
colnames(gene_sv2) <- c('patient_ID','SV_count')
gene_sv2$variant <- "No"
gene_sv2$variant[gene_sv2$SV_count  > 1 ] <- "Yes"

#merge metadata file + variant file
# metadata = input$metafile
library(readxl)
metadata <- read_excel("server/metadata.xlsx")                                                                          
metadata2 = data.frame(metadata)

# variant =gene_sv2
de <- merge(gene_sv2, metadata2, by=0, all=TRUE)  # merge by row names (by=0 or by="row.names")

de2 <- de[-(1)]
de2 <- de2[-(4)]




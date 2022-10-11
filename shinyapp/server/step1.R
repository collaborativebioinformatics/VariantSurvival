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
reticulate::repl_python()

source(table.py)

# ouput = sv_events.csv as dataframe for selected gene

allgenes_sv <- as.data.frame(sv_events)
gene_sv <- as.data.frame(allgenes_sv [c(1,5)])
gene_sv 
gene_sv2 <- gene_sv 
#change name column
# consider variant 
colnames(gene_sv2) <- c('patient_ID','gene')
gene_sv2$significant <- "No"
gene_sv2$significant[gene_sv2$gene  > 1 ] <- "Yes"
#










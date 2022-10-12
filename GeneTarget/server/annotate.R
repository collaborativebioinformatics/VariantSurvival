library(StructuralVariantAnnotation)
library(AnnotationHub)

ah <- AnnotationHub()
ah <- subset(ah, species == "Homo sapiens")


hist <- display(ah)
Human_gff = query(ah, c("Gencode", "gff", "human"))


library(rtracklayer)

#data = read.delim(file.choose(), sep="\t", skip=1)
#annotation = read.delim(file.choose(), sep="\t", skip=1)
#test = merge(data, annotation, by.x = "ID", by.y = "ORF")
#annotation_updated = import.gff(file.choose())

# Nieuwe versie van het script:

#library(rtracklayer)

data = read.delim(file.choose(), sep="\t", skip=1)
wcfs1_annotation = import.gff(file.choose())
data_annotated_wcfs1 = merge(data, wcfs1_annotation, by.x = "ID", by.y = "locus_tag")

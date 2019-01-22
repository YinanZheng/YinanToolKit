# Standard Enrichment analysis

### Install package if necessary

library(clusterProfiler)
library(DOSE)
library(AnnotationDbi)
library(org.Hs.eg.db)


### The input is the RefGene ID (NM_xxxxx or NR_xxxxx) from the Illumina Annotation file (i.e. the "UCSC_RefGene_Accession" column)

UCSC_RefGene_Accession <- NA

### 

## Split the ID and make a unique set
REGSEQ <- unique(unlist(strsplit(as.character(UCSC_RefGene_Accession), ";")))

## Translate ID to Entrez ID
geneSymbol <- select( org.Hs.eg.db, keys = REGSEQ,
                      columns = "ENTREZID",  keytype = "REFSEQ" )



### KEGG enrichment anlaysis

kk <- enrichKEGG(gene         = geneSymbol$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

write.csv(kk, "./Results/KEGG_enrich.csv")


#### GO enrichment anlaysis
ego <- enrichGO(gene         = geneSymbol$ENTREZID,
                OrgDb = 'org.Hs.eg.db', ont = "BP",
                pvalueCutoff = 0.05)
head(ego)

write.csv(ego, "./Results/GO_enrich.csv")


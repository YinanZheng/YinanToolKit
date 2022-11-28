# Manhattan plot

library(CMplot)

## Flat

## res_CC_plot is a dataframe with 4 columns: CpG ID, Chromosome, Position, p-value. (The column names does not matter) 

CMplot(res_CC_plot, plot.type="m", col=c("grey30"), LOG10=TRUE, threshold=c(5.965213e-08),
       threshold.lty=c(1), threshold.lwd=c(1), threshold.col=c("black"), amplify=TRUE,
       chr.den.col=NULL, signal.col=c("red"), cex = 0.8, signal.cex=c(0.8),signal.pch=c(19),
       file="jpg",memo="",dpi=600,file.output=TRUE,verbose=TRUE,width=14,height=6,chr.labels.angle=45)

## Note: Use those "highlight" parameter family to mark out the significant results. (Use CpG ID to specify which CpG to 
## highlight)


###################################################################################################################


# Enrichment plot

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(clusterProfiler)
library(missMethyl)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(ReactomePA)

annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

####

# cpg_list is a character vector of significant CpG IDs.

EntrezID<- getMappedEntrezIDs(sig.cpg = cpg_list, array.type = "EPIC")

## GO (BP)
ego <- enrichGO(gene          = EntrezID$sig.eg,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 0.2,
                readable      = TRUE)
data.frame(ego)
ego2 <- pairwise_termsim(ego)
p1 <- treeplot(ego2,nCluster = 3, method = "Wang")
p1

write.xlsx(data.frame(ego), file="GeneEnrichment_qval0.2.xlsx",
           sheetName="GO-BP", append=FALSE)

## GO (MF)
ego <- enrichGO(gene          = EntrezID$sig.eg,
                OrgDb         = org.Hs.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 0.2,
                readable      = TRUE)
data.frame(ego)
ego2 <- pairwise_termsim(ego)
p1 <- treeplot(ego2,nCluster = 3, method = "Wang")
p1

write.xlsx(data.frame(ego), file="GeneEnrichment_qval0.2.xlsx",
           sheetName="GO-MF", append=TRUE)


## KEGG
kk <- enrichKEGG(gene         = EntrezID$sig.eg,
                 organism     = 'hsa',
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 0.2)
head(kk)

kkx <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(kkx, showCategory = 10, cex_label_gene = 0.5)
p1

write.xlsx(data.frame(kk), file="GeneEnrichment_qval0.2.xlsx",
           sheetName="KEGG", append=TRUE)


## Reactome
x <- enrichPathway(gene=EntrezID$sig.eg,               
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)
head(x)
edox <- setReadable(x, 'org.Hs.eg.db', 'ENTREZID')

p1 <- cnetplot(edox,showCategory = 10, cex_label_gene = 0.5)
p1

write.xlsx(data.frame(edox), file="GeneEnrichment_qval0.2.xlsx",
           sheetName="Reactome", append=TRUE)


## Disease gene network
edo <- enrichDGN(EntrezID$sig.eg,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

p1 <- cnetplot(edox,showCategory = 10, cex_label_gene = 0.5)
p1

write.xlsx(data.frame(edox), file="GeneEnrichment_qval0.2.xlsx",
           sheetName="DGN", append=TRUE)

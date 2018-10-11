
## Produce descriptive table
descript <- function(dat, var, type = "continuous", digit = 1)
{
  d <- eval(parse(text = paste0("dat$", var)))
  if(type == "continuous")
  {
    res <- data.frame(Variable = var,
                      MeanOrFreq = mean(d,na.rm = T),
                      SDOrPerc = sd(d, na.rm = T),
                      Nmissing = sum(is.na(d)))
    res$Presentation <- paste0(round(res$MeanOrFreq, digit = digit), " (", round(res$SDOrPerc, digit = digit), ")")
  }
  if(type == "categorical")
  {
    temp <- table(d, useNA = "ifany")
    res <- data.frame(Variable = paste0(var, "-", names(temp)),
                      MeanOrFreq = as.numeric(temp),
                      SDOrPerc = as.numeric(temp)/length(d),
                      Nmissing = NA)
    res$Presentation <- paste0(res$MeanOrFreq, " (", round(res$SDOrPerc*100, digit = digit), ")")
  }
  res
}

## Remove outliers using IQR rule
IQRoutliers <- function(dat)
{
  iqr <- IQR(dat, na.rm = T)
  q <- quantile(dat, c(0.25, 0.75), na.rm = T)
  firstQ <- q[1]
  thirdQ <- q[2]
  outInd <- which(dat < firstQ - 3*iqr | dat > thirdQ + 3*iqr)
  if(length(outInd) > 0)
  {
    print(paste0(length(outInd), " outliers found!"))
    dat[outInd] <- NA
  } else {
    print("No outlier!")
  }
  return(dat)
}

## Impute missing with mean- row is variable (CpG methylation), column is sample
impute_w_rowMean <- function(dat){
  dat <- data.frame(dat, check.names = FALSE)
  for(i in 1:nrow(dat)){
    dat[i, is.na(dat[i, ])] <- mean(as.numeric(dat[i,]), na.rm = TRUE)
  }
  return(dat)
}

### lmCat toolkit function
lmCatCheck <- function(res, nlevel)
{
  cpglist <- unique(res$CpG)
  cpglist_bad <- cpglist[which(table(res$CpG) != nlevel)]
  if(length(cpglist_bad) == 0) message(paste0("All CpGs have ", nlevel, " groups tested.")) 
  else message(paste0("The following CpGs have less than ", nlevel, "group tested: ", paste0(cpglist_bad, collapse = ", ")))
}

lmCatSig <- function(res, sigcut, useFDR = TRUE)
{ 
  if(useFDR)
  {
    sigCpG <- na.omit(res[, c("CpG", "qvalue")])
    sigCpG <- subset(sigCpG, qvalue<sigcut)
    sigCpG <- as.character(sigCpG$CpG[order(sigCpG$qvalue)])
    mind <- match(res$CpG, sigCpG)
    res <- res[intersect(order(mind), which(!is.na(mind))), ]
  } else{
    sigCpG <- na.omit(res[, c("CpG", "TypeII_Pvalue")])
    sigCpG <- subset(sigCpG, TypeII_Pvalue<sigcut)
    sigCpG <- as.character(sigCpG$CpG[order(sigCpG$TypeII_Pvalue)])
    mind <- match(res$CpG, sigCpG)
    res <- res[intersect(order(mind), which(!is.na(mind))), ]
  }
  res
}

lmCatGetAnova <- function(res)
{
  res_anova <- res[which(res$Group == ""),]
  res_anova <- res_anova[, -c(2:6, 10:33)]
  res_anova
}

lmCatGetGroup <- function(res)
{
  res <- res[, -c(7:9, 34:35)]
  res
}


### Subset EWAS results by gene
geneInDB <- function(geneList, geneList_alias, TotalGeneList)
{
  if(all(geneList %in% TotalGeneList)) 
  {
    message("    All genes are covered!")
  } else {
    notinDB <- which(!geneList %in% TotalGeneList)
    notinDB_Gene <- geneList[notinDB]
    message(paste0("    ", paste0(notinDB_Gene, collapse = ", ")), " are not covered!")
    if(all(geneList_alias[notinDB] %in% TotalGeneList))
    {
      message(paste0("    But their alias, ", paste0(geneList_alias[notinDB], collapse = ", "), " are covered!"))
      geneList[notinDB] <- geneList_alias[notinDB]
    } else {
      notinDB_alias <- which(!geneList_alias[notinDB] %in% TotalGeneList)
      notinDB_alias_Gene <- geneList_alias[notinDB][notinDB_alias]
      notinDB_Gene <- geneList[notinDB][notinDB_alias]
      message(paste0("    The alias, ", paste0(notinDB_alias_Gene, collapse = ", "), " of ", paste0(geneList[notinDB][notinDB_alias], collapse = ", "), ", are still not covered!"))
      inDB_alias <- which(geneList_alias[notinDB] %in% TotalGeneList)
      geneList[notinDB][inDB_alias] <- geneList_alias[notinDB][inDB_alias]
    }
  }
  message(paste0("    Final working gene list: ", paste0(geneList, collapse = ", ")))
  return(geneList)
}

getGeneSubset <- function(res, geneList, 
                          array = c("EPIC","450K"), 
                          geneDB = c("UCSC", "GENCODE", "BOTH"),
                          promoter = FALSE,
                          FDR = TRUE)
{
  library(org.Hs.eg.db)
  library(DBI)
  
  array <- match.arg(array)
  geneDB <- match.arg(geneDB)
  
  dbCon <- org.Hs.eg_dbconn()
  sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
  geneList_alias <- aliasSymbol[match(geneList, aliasSymbol[,2]),5]
  
  if(length(geneList_alias) != length(geneList)) stop("Something wrong. Check gene symbol!")
  
  geneList <- toupper(geneList)
  geneList_alias <- toupper(geneList_alias)

  if(array == "EPIC")
  {
    message("Loading EPIC annotation...")
    library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
    annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
    annot <- annot[match(rownames(res), annot$Name),]
  }
  
  if(array == "450K")
  {
    message("Loading 450K annotation...")
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    annot <- annot[match(rownames(res), annot$Name),]
  }
  
  totalGeneList_UCSC <- unique(unlist(strsplit(annot$UCSC_RefGene_Name, ";")))
  totalGeneList_GENCODE <- unique(unlist(strsplit(annot$GencodeCompV12_NAME, ";")))
  
  promoter_UCSC <- grep("TSS1500|TSS200|5'UTR|1stExon", annot$UCSC_RefGene_Group)
  promoter_GENCODE <- grep("TSS1500|TSS200|5'UTR|1stExon", annot$GencodeCompV12_Group)
  
  if(geneDB == "UCSC")
  {
    message("Using UCSC gene DB...")
    geneList_UCSC_FINAL <- geneInDB(geneList, geneList_alias, totalGeneList_UCSC)
    geneList_grep <- paste0("\\b", paste0(geneList_UCSC_FINAL, collapse = "\\b|\\b"), "\\b")
    ind <- grep(geneList_grep, toupper(annot$UCSC_RefGene_Name))
    if(promoter) {message("Note: Promoter only."); ind <- intersect(ind, promoter_UCSC)}
  }
  if(geneDB == "GENCODE")
  {
    message("Using GENCODE gene DB...")
    geneList_GENCODE_FINAL <- geneInDB(geneList, geneList_alias, totalGeneList_GENCODE)
    geneList_grep <- paste0("\\b", paste0(geneList_GENCODE_FINAL, collapse = "\\b|\\b"), "\\b")
    ind <- grep(geneList_grep, toupper(annot$GencodeCompV12_NAME))
    if(promoter) {message("Note: Promoter only."); ind <- intersect(ind, promoter_GENCODE)}
  }
  if(geneDB == "BOTH")
  {
    message("Using UCSC gene DB...")
    geneList_UCSC_FINAL <- geneInDB(geneList, geneList_alias, totalGeneList_UCSC)
    geneList_UCSC_grep <- paste0("\\b", paste0(geneList_UCSC_FINAL, collapse = "\\b|\\b"), "\\b")
    
    message("Using GENCODE gene DB...")
    geneList_GENCODE_FINAL <- geneInDB(geneList, geneList_alias, totalGeneList_GENCODE)
    geneList_GENCODE_grep <- paste0("\\b", paste0(geneList_GENCODE_FINAL, collapse = "\\b|\\b"), "\\b")
    
    ind1 <- grep(geneList_UCSC_grep, toupper(annot$UCSC_RefGene_Name))
    ind2 <- grep(geneList_GENCODE_grep, toupper(annot$GencodeCompV12_NAME))
    ind <- union(ind1, ind2)
    if(promoter) {message("Note: Promoter only."); ind <- intersect(ind, union(promoter_UCSC, promoter_GENCODE))}
  }
  
  res_sub <- res[ind,]
  annot_sub <- annot[ind,]
  
  if(FDR)
  {
    res_sub$FDR <- p.adjust(res_sub$Pvalue)
  }
  
  res_sub <- cbind(res_sub, annot_sub)
  
  res_sub <- res_sub[order(res_sub$Pvalue),]
  
  return(res_sub)
}


function_list <- c("descript",
                   "IQRoutliers",
                   "impute_w_rowMean",
                   "DNAmAgeAccel",
                   "calibrateDNAmAge",
                   "NewDNAmAgeCleaner",
                   "lmCatCheck",
                   "lmCatSig",
                   "lmCatGetAnova",
                   "lmCatGetGroup",
                   "geneInDB",
                   "getGeneSubset",
                   "UCSCtoGRanges")

trash <- sapply(function_list, function(x) suppressWarnings(rm(x)))
                
message("The following tools have been succesfully loaded:")
trash <- sapply(function_list, function(x) message(x))

                
### Produce descriptive table
round_pad <- function(x, digits=0)
{
 format(round(x, digits), nsmall=digits)
}

descript <- function(dat, var, type = "continuous", digits = 1)
{
  d <- eval(parse(text = paste0("dat$", var)))
  if(type == "continuous")
  {
    res <- data.frame(Variable = var,
                      MeanOrFreq = mean(d,na.rm = T),
                      SDOrPerc = sd(d, na.rm = T),
                      Nmissing = sum(is.na(d)))
    res$Presentation <- paste0(round_pad(res$MeanOrFreq, digits = digits), " (", round_pad(res$SDOrPerc, digits = digits), ")")
  }
  if(type == "categorical")
  {
    temp <- table(d, useNA = "ifany")
    res <- data.frame(Variable = paste0(var, "-", names(temp)),
                      MeanOrFreq = as.numeric(temp),
                      SDOrPerc = as.numeric(temp)/length(d),
                      Nmissing = NA)
    res$Presentation <- paste0(res$MeanOrFreq, " (", round_pad(res$SDOrPerc*100, digits = digits), ")")
  }
  res
}

### Remove outliers using IQR rule
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

### Function to impute missing with mean- row is variable (CpG methylation), column is sample
impute_w_rowMean <- function(dat){
  dat <- data.frame(dat, check.names = FALSE)
  for(i in 1:nrow(dat)){
    dat[i, is.na(dat[i, ])] <- mean(as.numeric(dat[i,]), na.rm = TRUE)
  }
  return(dat)
}

### Function to compute Residule- note that the new column name is DNAmPhenoAgeAccel!
DNAmAgeAccel <- function(DNAmAge, Age, ID)
{
  dat <- data.frame(DNAmAge = DNAmAge, Age = Age)
  rownames(dat) <- ID
  DNAmPhenoAgeAccel <- resid(lm(DNAmAge ~ Age, data = dat))
  return(data.frame(SampleID = names(DNAmPhenoAgeAccel), DNAmPhenoAgeAccel = DNAmPhenoAgeAccel))
}

### Function to Calibrate mAge
calibrateDNAmAge <- function(DNAmAge, Age)
{
  dat <- data.frame(DNAmAge = DNAmAge, Age = Age)
  coe <- summary(lm(DNAmAge ~ Age, data = dat))$coefficient
  intercept <- coe[1,1]
  slope <- coe[2,1]
  DNAmAge_calibrated = (DNAmAge - intercept) / slope
  return(DNAmAge_calibrated)
}

### For new DNAm age calculator only. This function removes outliers for commonly used epigenetic age variables.
NewDNAmAgeCleaner <- function(DNAmAge_Output, filename)
{
  DNAmAge_Output$HorvathDNAmAge_clean <- DNAmAge_Output$DNAmAge
  message("Horvath DNAmAge:")
  DNAmAge_Output$HorvathDNAmAge_clean <- IQRoutliers(DNAmAge_Output$HorvathDNAmAge_clean)
  
  DNAmAge_Output$HannumDNAmAge_clean <- DNAmAge_Output$DNAmAgeHannum
  message("Hannum DNAmAge:")
  DNAmAge_Output$HannumDNAmAge_clean <- IQRoutliers(DNAmAge_Output$HannumDNAmAge_clean)
  
  
  ####Remove outliers for original PhenoAge and GrimAge
  DNAmAge_Output$DNAmPhenoAge_clean <- DNAmAge_Output$DNAmPhenoAge
  message("DNAmPhenoAge:")
  DNAmAge_Output$DNAmPhenoAge_clean <- IQRoutliers(DNAmAge_Output$DNAmPhenoAge_clean)
  
  DNAmAge_Output$DNAmGrimAge_clean <- DNAmAge_Output$DNAmGrimAge
  message("DNAmGrimAge:")
  DNAmAge_Output$DNAmGrimAge_clean <- IQRoutliers(DNAmAge_Output$DNAmGrimAge_clean)
  
  ####Remove outliers for original SkinAge
  DNAmAge_Output$DNAmSkinBloodAge_clean <- DNAmAge_Output$DNAmAgeSkinBloodClock
  message("DNAmSkinBloodAge:")
  DNAmAge_Output$DNAmSkinBloodAge_clean <- IQRoutliers(DNAmAge_Output$DNAmSkinBloodAge_clean)
  
  ####Remove outliers for original mPAI-1
  DNAmAge_Output$DNAmPAI1_clean <- DNAmAge_Output$DNAmPAI1
  message("DNAmPAI1:")
  DNAmAge_Output$DNAmPAI1_clean <- IQRoutliers(DNAmAge_Output$DNAmPAI1_clean)
  
  ####Removing outliers for IEAA and EEAA
  DNAmAge_Output$IEAA_clean <- DNAmAge_Output$IEAA
  message("IEAA:")
  DNAmAge_Output$IEAA_clean <- IQRoutliers(DNAmAge_Output$IEAA_clean)
  
  DNAmAge_Output$EEAA_clean <- DNAmAge_Output$EEAA
  message("EEAA:")
  DNAmAge_Output$EEAA_clean <- IQRoutliers(DNAmAge_Output$EEAA_clean)
  
  ####Calculate PAA, GAA, SBAA, using epigenetic age removed outliers
  ####Compute PAA(PhenoAge Acceleration)and GAA(GrimAge Acceleration)
  PAA <- DNAmAgeAccel(DNAmAge_Output$DNAmPhenoAge_clean, DNAmAge_Output$Age, DNAmAge_Output$SampleID)
  colnames(PAA)[colnames(PAA)=="DNAmPhenoAgeAccel"] <- "PAA_clean"
  
  GAA <- DNAmAgeAccel(DNAmAge_Output$DNAmGrimAge_clean, DNAmAge_Output$Age, DNAmAge_Output$SampleID)
  colnames(GAA)[colnames(GAA)=="DNAmPhenoAgeAccel"] <- "GAA_clean"
  
  ####Compute SBAA(SkinBloodAge Acceleration)
  SBAA <- DNAmAgeAccel(DNAmAge_Output$DNAmSkinBloodAge_clean, DNAmAge_Output$Age, DNAmAge_Output$SampleID)
  colnames(SBAA)[colnames(SBAA)=="DNAmPhenoAgeAccel"] <- "SBAA_clean"
  
  DNAmAge_Output <- merge(DNAmAge_Output, PAA, by ="SampleID", all = TRUE)
  DNAmAge_Output <- merge(DNAmAge_Output, GAA, by ="SampleID", all = TRUE)
  DNAmAge_Output <- merge(DNAmAge_Output, SBAA, by ="SampleID", all = TRUE)
  
  # ####Calibrate Horvath Age and Hannum Age
  # DNAmAge_Output$HorvathDNAmAge_Calibrated <- calibrateDNAmAge(DNAmAge_Output$HorvathDNAmAge_clean, DNAmAge_Output$Age) 
  # DNAmAge_Output$HannumDNAmAge_Calibrated <- calibrateDNAmAge(DNAmAge_Output$HannumDNAmAge_clean, DNAmAge_Output$Age) 
  # 
  # ####Calibrate PhenoAge and GrimAge
  # DNAmAge_Output$PhenoAge_Calibrated <- calibrateDNAmAge(DNAmAge_Output$DNAmPhenoAge_clean, DNAmAge_Output$Age) 
  # DNAmAge_Output$DNAmGrimAge_Calibrated <- calibrateDNAmAge(DNAmAge_Output$DNAmGrimAge_clean, DNAmAge_Output$Age) 
  
  write.csv(DNAmAge_Output, file = paste0(filename, ".csv"), row.names = FALSE)
  return(DNAmAge_Output)
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


### Subset EWAS results or Methylation data by gene
# Note: either EWAS or Methylation must of CpG ID as the row names!
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
  message(paste0("    Final working gene list (", length(geneList), "): ", paste0(geneList, collapse = ", ")))
  return(geneList)
}

getGeneSubset <- function(res, geneList, 
                          input = c("EWAS", "METHY"),
                          array = c("EPIC","450K"), 
                          geneDB = c("UCSC", "GENCODE", "BOTH"),
                          promoter = FALSE,
                          FDR = TRUE)
{
  # res <- as.data.frame(res, check.names = FALSE)
  library(org.Hs.eg.db)
  library(DBI)
  
  stopifnot(!is.null(rownames(res)) & "cg" %in% substring(rownames(res), 1,2))
  
  input <- match.arg(input)
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
    message("Loading EPIC annotation IlluminaHumanMethylationEPICanno.ilm10b4.hg19...")
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    annot <- annot[match(rownames(res), annot$Name),]
  }
  
  if(array == "450K")
  {
    message("Loading 450K annotation IlluminaHumanMethylation450kanno.ilmn12.hg19...")
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    annot <- annot[match(rownames(res), annot$Name),]
  }
  
  totalGeneList_UCSC <- toupper(unique(unlist(strsplit(annot$UCSC_RefGene_Name, ";"))))
  totalGeneList_GENCODE <- toupper(unique(unlist(strsplit(annot$GencodeCompV12_NAME, ";"))))
  
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
  
  res_sub <- as.data.frame(res[ind,], check.names = FALSE)
  annot_sub <- as.data.frame(annot[ind,], check.names = FALSE)
  
  if(input == "EWAS" & FDR)
  {
    message("Your input is EWAS results. Combining annotation data with subset results...")
    res_sub$FDR <- p.adjust(res_sub$Pvalue)
    res_sub <- cbind(res_sub, annot_sub)
    res_sub <- res_sub[order(res_sub$Pvalue),]
  }
  
  if(input == "METHY")
  {
    message("Your input is methyltion data. Preparing a list object to keep the subset of methylation and annotation data...")
    ord <- order(annot_sub$chr , annot_sub$pos)
    res_sub <- res_sub[ord, ]
    annot_sub <- annot_sub[ord, ]
    res_sub = list(methylation = res_sub, annotation = annot_sub)
  }
  
  return(res_sub)
}

### Convert a UCSC genomic region text (chr:start-end) to a GRanges object
UCSCtoGRanges <- function(text)
{
  library(GenomicRanges)
  split1 <- strsplit(text, ":")[[1]]
  chr <- split1[1]
  gr <- split1[2]
  split2 <- strsplit(gr, "-")[[1]]
  start <- as.integer(split2[1])
  end <- as.integer(split2[2])
  GR <- GRanges(seqnames = chr, IRanges(start = start, end = end))
  return(GR)
}

function_list <- c("round_pad",  
                   "descript",
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
                   "UCSCtoGRanges",
                   "corPlot",
                   "extractSNPdat",
                   "CpGAnnot_GENCODE",
                   "matrixHeatmap")

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

### Convert a UCSC genomic region text/text vector (chr:start-end or chr:pos) to a GRanges object
UCSCtoGRanges <- function(text)
{
  library(GenomicRanges)
  text <- gsub(",","",text)
  split1 <- do.call(rbind, strsplit(text, ":"))
  chr <- split1[,1]
  gr <- split1[,2]
  
  if(all(!grepl("-", text)))
  {
    message("Single position format...")
    start <- as.integer(gr)
    end <- as.integer(gr)
  } else {
    message("Range format...")
    split2 <- do.call(rbind, strsplit(gr, "-"))
    start <- as.integer(split2[,1])
    end <- as.integer(split2[,2])
  }
  
  GR <- GRanges(seqnames = chr, IRanges(start = start, end = end))
  return(GR)
}


## Make a quick scatterplot with regression line and correlation test
## x and y are column names in dat
corPlot <- function(dat, x, y, xlab = x, ylab = y, xlimit = NULL, ylimit = NULL, 
                    smoothtype = "lm", corlab_x = "left", corlab_y = "top", main = NULL)
{
  library(ggplot2)
  library(ggpubr)
  
  p <- eval(parse(text = paste0("ggplot(data = dat, aes(", x, ", ", y, "))")))
  p = p + theme_classic()
  
  p = p + geom_point(shape = 21, colour = "dark grey", fill = "#52478B", size = 2) + 
    theme(axis.text.x = element_text(color = "black", size=14), 
          axis.text.y = element_text(color = "black", size=14),
          axis.title = element_text(size=14))
  
  p = p + geom_smooth(method = smoothtype, colour = "black") + 
    xlab(xlab) + ylab(ylab)
  
  p = p + stat_cor(method = "pearson", label.x.npc = corlab_x, label.y.npc = corlab_y)
  
  
  
  if(!is.null(xlimit)) p = p + xlim(xlimit)
  if(!is.null(ylimit)) p = p + xlim(ylimit)
  if(!is.null(main)) p = p + ggtitle(main)
  
  return(p)
}

## Extract SNP from CARDIA genotype dataset in Quest
# race: "Black" or "White"
# library(BSgenome)
# available.SNPs()
# library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
# snpDB <- SNPlocs.Hsapiens.dbSNP144.GRCh37
                
extractSNPdat <- function(snpDB, rsids = NULL, seqnames = NULL, pos = NULL, race, infoOnly = TRUE)
{  
  if(all(is.na(seqnames)))
  {
    message("Locating SNP position...")
    my_snps <- snpsById(snpDB, rsids, ifnotfound='warning')
    my_snps <- sort(my_snps)
  } else {
    my_snps <- GRanges(seqnames = seqnames, ranges = IRanges(start = pos, end = pos), RefSNP_id = rsids)
  }

  chrList <- as.character(seqnames(my_snps))
  snpList <- as.character(my_snps$RefSNP_id)
  
  lastChr = ""
  
  for(i in 1:length(snpList))
  {
    rsid <- snpList[i]
    chr <- chrList[i]
    
    message("Processing ", rsid, " in chromosome ", chr)
    
    INFO_select <- NULL
    DOSAGE_select <- NULL
    GENOTYPE_select <- NULL
    
    if(chr != lastChr)
    {
      QuestPath <- file.path("/projects/b1096/CARDIA/GENOTYPE/Imputation", paste0(race, "_filter_RDS"))
      INFO_Path <- file.path(QuestPath, paste0(race, "_chr", chr, ".Minimac3_MI_1KGP_p3v5_typed_rsid_ALL.MAF0.01_Rsq0.3_INFO.rds"))
      DOSAGE_Path <- file.path(QuestPath, paste0(race, "_chr", chr, ".Minimac3_MI_1KGP_p3v5_typed_rsid_ALL.MAF0.01_Rsq0.3_DOSAGE.rds"))
      GENOTYPE_Path <- file.path(QuestPath, paste0(race, "_chr", chr, ".Minimac3_MI_1KGP_p3v5_typed_rsid_ALL.MAF0.01_Rsq0.3_GENOTYPE.rds"))
      
      INFO_dat <- readRDS(INFO_Path)
      
      if(!infoOnly)
      {
        DOSAGE_dat <- readRDS(DOSAGE_Path)
        GENOTYPE_dat <- readRDS(GENOTYPE_Path)
      }
      lastChr <- chr
    }
    
    SNPind <- which(INFO_dat$SNP_ID == rsid)
    
    if(length(SNPind) > 0)
    {
      INFO_select <- rbind(INFO_select, INFO_dat[SNPind,])
      
      if(!infoOnly)
      {
        DOSAGE_select <- rbind(DOSAGE_select, DOSAGE_dat[SNPind,])
        GENOTYPE_select <- rbind(GENOTYPE_select, GENOTYPE_dat[SNPind,])
      }
    } else {
      message("SNP ", rsid, " does not exist in the database!")
    }
  }
  return(list(INFO = INFO_select, DOSAGE = DOSAGE_select, GENOTYPE = GENOTYPE_select))
}

## GDM version (not imputed data only genotype data)
## Mother == TRUE then extract Mother's data, Mother == FALSE then extract Children's data
extractSNPdat <- function(snpDB, rsids = NULL, seqnames = NULL, pos = NULL, Mother = TRUE, infoOnly = TRUE)
{
  if(all(is.na(seqnames)))
  {
    message("Locating SNP position...")
    
    my_snps <- snpsById(snpDB, rsids, ifnotfound='warning')
    my_snps <- sort(my_snps)
  } else {
    my_snps <- GRanges(seqnames = seqnames, ranges = IRanges(start = pos, end = pos), RefSNP_id = rsids)
  }
  
  if(Mother) {message("Getting Mother's data..."); who = "MOTHER"} else {message("Getting Children's data..."); who = "CHILDREN"}
  
  chrList <- as.character(seqnames(my_snps))
  snpList <- as.character(my_snps$RefSNP_id)
  
  lastChr = ""
  
  for(i in 1:length(snpList))
  {
    rsid <- snpList[i]
    chr <- chrList[i]
    
    message("Processing ", rsid, " in chromosome ", chr)
    
    INFO_select <- NULL
    GENOTYPE_select <- NULL
    
    if(chr != lastChr)
    {
      QuestPath <- file.path("/projects/b1096/GDM/GWAS/Children_Mothers_Merged/RDS")
      INFO_Path <- file.path(QuestPath, paste0("chr", chr, "_rsid_INFO.rds"))
      GENOTYPE_Path <- file.path(QuestPath, paste0("chr", chr, "_", who, "_rsid_GENOTYPE.rds"))

      INFO_dat <- readRDS(INFO_Path)
      
      if(!infoOnly)
      {
        GENOTYPE_dat <- readRDS(GENOTYPE_Path)
      }
      lastChr <- chr
    }
    
    SNPind <- which(INFO_dat$SNP_ID == rsid)
    
    if(length(SNPind) > 0)
    {
      message("SNP ", rsid, " found in the database!")
      INFO_select <- rbind(INFO_select, INFO_dat[SNPind,])
      
      if(!infoOnly)
      {
        GENOTYPE_select <- rbind(GENOTYPE_select, GENOTYPE_dat[SNPind,])
      }
    } else {
      message("SNP ", rsid, " does not exist in the database!")
    }
  }
  return(list(INFO = INFO_select, GENOTYPE = GENOTYPE_select))
}                
                
## Annotate CpGs with GENCODE database                
# Download GENCODE gene annotation file GTF: https://www.gencodegenes.org/human/ https://www.gencodegenes.org/human/release_33lift37.html
# GENCODE_GTF <- rtracklayer::import("./Data/gencode.v33lift37.annotation.gtf.gz", format = "gtf")
# cpg: the list of CpG to be annotated.

# Outcome dataset -> Distance:  Negative distance = upstream of gene
CpGAnnot_GENCODE <- function(GENCODE_GTF, cpg, arrayType, win = 1000000)
{
   if(arrayType == "EPIC")
   {
      library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
      annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
      annot <- subset(annot, Name %in% cpg)
   }

   if(arrayType == "450K")
   {
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      annot <- subset(annot, Name %in% cpg)
   }

   message("Flanking region window size: ", win)
   
   annot_GR_extended <- GRanges(seqnames = annot$chr, IRanges(start = annot$pos-win, end = annot$pos+win), 
                                strand = as.character(annot$strand), CpG_ID = annot$Name, CpG_Pos = annot$pos)
   
   hits <- findOverlaps(annot_GR_extended, GENCODE_GTF, ignore.strand=TRUE)

   CpG_annot_full <- data.frame(CpG = annot_GR_extended$CpG_ID[queryHits(hits)],
                                CpG_Pos = annot_GR_extended$CpG_Pos[queryHits(hits)],
                                CpG_strand = strand(annot_GR_extended[queryHits(hits)]),
                                GENCODE_GTF[subjectHits(hits)])
   
   CpG_annot_full_reverse <- subset(CpG_annot_full, strand == "-")
   CpG_annot_full_forward <- subset(CpG_annot_full, strand == "+")
   
   CpG_annot_full_reverse$distance <- CpG_annot_full_reverse$end - CpG_annot_full_reverse$CpG_Pos ## Negative distance = upstream of gene
   CpG_annot_full_forward$distance <- CpG_annot_full_forward$CpG_Pos - CpG_annot_full_forward$start ## Negative distance = upstream of gene
   
   # CpG-Gene combination by strand
   CpG_annot_full_reverse_nearby <- CpG_annot_full_reverse %>% 
      group_by(CpG, gene_name) %>%
      filter(abs(distance) == min(abs(distance))) %>% select(CpG, CpG_Pos, CpG_strand, seqnames, start, end, strand, gene_name, gene_type, distance)
   
   CpG_annot_full_forward_nearby <- CpG_annot_full_forward %>% 
      group_by(CpG, gene_name) %>%
      filter(abs(distance) == min(abs(distance))) %>% select(CpG, CpG_Pos, CpG_strand, seqnames, start, end, strand, gene_name, gene_type, distance)
   
   CpG_annot_full_nearby <- rbind(CpG_annot_full_reverse_nearby, CpG_annot_full_forward_nearby)
   
   # by CpG
   CpG_annot_full_nearby <- CpG_annot_full_nearby %>% 
      group_by(CpG) %>% filter(abs(distance) == min(abs(distance))) %>%
      arrange(CpG)
   
   ## Organize results
   CpG_annot_full_nearby_summarise <- CpG_annot_full_nearby %>%
      group_by(CpG) %>%
      summarise(Chr = as.integer(gsub("chr", "", as.character(unique(seqnames)))), CpG_Pos = unique(CpG_Pos), CpG_strand = unique(CpG_strand), 
                Gene_strand = unique(strand), start = min(start), end = max(end), 
                gene_name = paste(unique(gene_name), collapse = ";"),
                gene_type = paste(unique(gene_type), collapse = ";"),
                Distance = paste(unique(distance), collapse = ";")) %>%
      arrange(Chr, CpG_Pos)
   
   CpG_annot_full_nearby_summarise$Abs_Distance <- abs(as.integer(CpG_annot_full_nearby_summarise$Distance))
   CpG_annot_full_nearby_summarise$UCSC_POS <- with(CpG_annot_full_nearby_summarise, paste0("chr", Chr, ":", CpG_Pos))
   
   n_cpg <- length(cpg)
   n_cpg_annot <- length(unique(CpG_annot_full_nearby_summarise$CpG))
   
   if(n_cpg == n_cpg_annot) message("All ", n_cpg, " CpGs are annotated!") else message("Only ", n_cpg_annot, " out of ", n_cpg, " are annotated!")
   return(CpG_annot_full_nearby_summarise)
}

                
                
# ## Visualize a matrix value with heatmap - it is useful to visualize correlation matrix
# data(iris)
# mat = as.matrix(iris[,1:4])
# res = matrixHeatmap(mat, legendName = "Correlation")
# res = matrixHeatmap(mat, showValue = TRUE, marksize = 5, legendName = "Correlation")
# res
# ## Save it as tiff file
# matrixHeatmap(mat, filename = "matrixHeatmap_irisExample")
matrixHeatmap<-function(mat, corrMat = FALSE, xlab = NULL, ylab = NULL, showValue = FALSE, marksize = 2, adjustMethod = "holm", legendName = "NULL", width = 5, height = 5, res = 300, filename = NULL)
{
  library(reshape2)
  library(RColorBrewer)
  library(ggplot2)
  library(psych)
  
  if(corrMat)
  {
    r = mat
    showValue = F
  } else {
    corTest = corr.test(mat, adjust = adjustMethod)
    r = corTest$r
    poriginal = padjust = corTest$p
    
    poriginal[upper.tri(poriginal)] <- t(poriginal)[upper.tri(poriginal)]
    padjust[lower.tri(padjust)] <- t(padjust)[lower.tri(padjust)]
    
    # For melt() to work seamlessly, myData has to be a matrix.
    longpadjust = melt(padjust)
  }
  
  longData <- melt(r)
  
  longData[,1] = factor(longData[,1], levels = rownames(r))
  longData[,2] = factor(longData[,2], levels = colnames(r))
  
  colnames(longData) = c("row","col","value")
  
  if(!corrMat)
  {
    colnames(longpadjust) = c("row","col","value")
    
    longData$mark = sprintf("%.2f", longData$value)
    twoStar = longpadjust$value<0.01
    oneStar = longpadjust$value<0.05 & longpadjust$value>=0.01
    if(sum(twoStar)>0) longData$mark[twoStar] = paste0(longData$mark[twoStar],"\n**")
    if(sum(oneStar)>0) longData$mark[oneStar] = paste0(longData$mark[oneStar],"\n*")
  }
  
  # Define palette
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  
  # Plot
  zp1 <- ggplot(longData,
                aes(x = col, 
                    y = row, fill = value))
  zp1 <- zp1 + xlab(xlab) + ylab(ylab)
  zp1 <- zp1 + geom_tile()
  if(showValue)   zp1 <- zp1 + geom_text(aes(fill = longData$value, label = longData$mark), size = marksize)
  zp1 <- zp1 + scale_fill_gradientn(colours = myPalette(100), values = seq(-1.2, 1.2, length = 100), name = legendName)
  zp1 <- zp1 + scale_x_discrete(expand = c(0, 0))
  zp1 <- zp1 + scale_y_discrete(expand = c(0, 0))
  zp1 <- zp1 + coord_equal()
  zp1 <- zp1 + theme_bw()
  zp1 <- zp1 + theme(legend.text = element_text(size=5), 
                     legend.title = element_text(size=9, face = "italic", vjust = 0.2))
  zp1 <- zp1 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  if(!is.null(filename))
  {
    tiff(file = paste0(filename,".tiff"), width = width, height = height, unit = "in", res = res, compress = "lzw")
    print(zp1)  
    dev.off()
  } else print(zp1)
  
  if(corrMat)   return(list(r=r))
  if(!corrMat)  return(list(r=r, p=poriginal, p.adj=padjust, adjust.method=adjustMethod))
}

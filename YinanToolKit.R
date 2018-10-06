
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


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
  dat <- data.frame(dat)
  for(i in 1:nrow(dat)){
    dat[i, is.na(dat[i, ])] <- mean(as.numeric(dat[i,]), na.rm = TRUE)
  }
  return(dat)
}

### vcftools example: extract data by posision
#vcftools --vcf chr18.dose.vcf --chr chr18 --from-bp 31598555 --to-bp 31598755 --recode --recode-INFO-all --out chr18_subset_rs76992529.snp
#vcftools --gzvcf chr6.dose.vcf.gz --chr chr6 --from-bp 26092813 --to-bp 26093013 --recode --recode-INFO-all --out chr6_subset_rs1800562.snp

setwd("C:\\Users\\yzk256\\Dropbox\\MESA")

readVCF <- function(f)
{
  tmp_vcf<-readLines(f)
  tmp_vcf_data<-read.table(f, stringsAsFactors = FALSE)
  
  # filter for the columns names
  tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
  vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
  names(tmp_vcf_data)<-vcf_names
  return(tmp_vcf_data)
}

vcf_geno <- function(x) sapply(strsplit(as.character(x), ":"),function(x) x[1])
vcf_dosage <- function(x) as.numeric(sapply(strsplit(as.character(x), ":"),function(x) x[2]))

genotypeTranslate <- function(genotype, REF, ALT)
{
  genotype <- as.character(genotype)
  genotype[genotype == "0|0"] <- paste0(REF, REF)
  genotype[genotype %in% c("1|0", "0|1")] <- paste0(REF, ALT)
  genotype[genotype == "1|1"] <- paste0(ALT, ALT)
  return(genotype)
}

reorganize <- function(RES, trait, model, ethnicity, study, array, filetype = c("txt", "csv", "rds"), gz = FALSE){
  filetype <- match.arg(filetype)
  RES_mod <- data.frame(ProbeID = rownames(RES), BETA = RES$Estimates, 
                        SE = RES$StdErr, P_VAL = RES$Pvalue,
                        N_samp = RES$Sample_Size, mean = RES$beta_Mean, SD = RES$beta_SD)
  out_file = paste(trait, model, ethnicity, study, array, format(Sys.time(),format="%m%d%Y"), sep = "_")
  
  if(filetype == "rds") {
    writeRDS(RES_mod, paste0(out_file, ".rds"))
  } else if (filetype == "csv") {
    if(gz){
      gz <- gzfile(paste0(out_file, ".csv.gz"), "w")
      write.csv(RES_mod, gz, quote=F, row.names=F)
      close(gz)
    } else {
      write.csv(RES_mod, paste0(out_file, ".csv"), quote=F, row.names=F)
    }
  } else {
    if(gz) {
      gz <- gzfile(paste0(out_file, ".txt.gz"), "w")
      write.table(RES_mod, gz, quote=F, row.names=F,sep="\t")
      close(gz)
    } else {
      write.table(RES_mod, paste0(out_file, ".txt"), quote=F, row.names=F,sep="\t")
    }
  }
  print(head(RES_mod))
}

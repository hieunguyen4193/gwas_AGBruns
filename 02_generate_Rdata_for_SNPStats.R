gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full_corrected_20240312"

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))
library(karyoploteR)

thres_call.freq <- 0.5
thres_hwe <- 0.05
outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "AGBruns_SNP_data_full_corrected_20240312"
for (thres_maf in c(0.01, 0.05)){
  
  path.to.maindir <- file.path(outdir, PROJECT)
  dir.create(file.path(path.to.maindir, "snpStats_RData"), showWarnings = FALSE, recursive = TRUE)
  all.GWAS.cases <- read.csv(file.path(path.to.project.src, "all_GWAS_cases.txt"), sep = "\t", header = FALSE)[["V1"]]
  
  for (group in all.GWAS.cases){
    
    print("---------------------------------------------------------------------")
    print(sprintf("Working on %s", group))
    print("---------------------------------------------------------------------")
    path.to.save.RData <- file.path(path.to.maindir, "snpStats_RData", sprintf("03_run_snpStats_workspace_CallFreq_%s_HWE_%s_MAF_%s.%s.RData", thres_call.freq, thres_hwe, thres_maf, group))
    
    if (file.exists(path.to.save.RData) == FALSE){
      print("Reading in ped data file...")
      path.to.snp.full.info1 <- file.path(path.to.maindir,  "01_output", "snp_table_full_info_plate1.map")
      path.to.snp.full.info2 <- file.path(path.to.maindir,  "01_output", "snp_table_full_info_plate2.map")
      
      snp.full.info1 <- read.csv(path.to.snp.full.info1, check.names = FALSE)
      snp.full.info2 <- read.csv(path.to.snp.full.info2, check.names = FALSE)
      
      path.to.pedfile <- file.path(path.to.maindir,  "01_output",
                                   sprintf("parallel_processing_%s", group), 
                                   "output_chunks", sprintf("%s.ped", group))
      path.to.info.file <- file.path(path.to.maindir,  "01_output", "snp_table_for_GWAS.map")
      
      peddf <- read.pedfile(path.to.pedfile, snps = path.to.info.file)
      print("Finish reading in ped data file. ")
      print("Preprocessing")
      snps.data <- peddf$genotypes
      snp.support <- peddf$map
      subject.support <- peddf$fam
      
      # replace NA value in SBP condtion by 0
      subject.support$affected[is.na(subject.support$affected)] <- 0
      
      sample.qc <- row.summary(snps.data)
      
      snp.full.info1 <- subset(snp.full.info1, select = -c(Index)) %>% subset(select = c(Name, Chr, Position, ChiTest100, Call.Freq, Minor.Freq))
      colnames(snp.full.info1) <- c(c("Name", "Chr", "Position"), to_vec(for (item in colnames(snp.full.info1)[4:6]) sprintf("%s_plate1", item)))
      snp.full.info2 <- subset(snp.full.info2, select = -c(Index)) %>% subset(select = c(Name, ChiTest100, Call.Freq, Minor.Freq))
      colnames(snp.full.info2) <- c("Name", "ChiTest100_plate2", "Call.Freq_plate2", "Minor.Freq_plate2")
      
      snp.full.info <- merge(snp.full.info1, snp.full.info2, by.x = "Name", by.y = "Name")
      snp.full.info <- subset(snp.full.info, snp.full.info$Chr %in% c("X", "Y", "MT", "XY") == FALSE)
      source(file.path(path.to.project.src, "00_preprocessing_SNP_full_info_table.R"))
      
      print("Precessing finished. Saving workspace...")
      save.image(file = path.to.save.RData)  
      print("Finish!")
    } 
  }
}


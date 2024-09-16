gc()
rm(list = ls())

source("/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full_corrected_20240312/02_generate_Rdata_for_SNPStats.R")

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "AGBruns_SNP_data_full_corrected_20240312"
path.to.save.HTML.dir <- file.path(outdir, PROJECT, "html_reports_20240312")
dir.create(path.to.save.HTML.dir, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full_corrected_20240312"
all.GWAS.cases <- read.csv(file.path(path.to.project.src, "all_GWAS_cases.txt"), sep = "\t", header = FALSE)[["V1"]]

#####----------------------------------------------------------------------#####
##### BASIC ASSOC TEST RESULTS
#####----------------------------------------------------------------------#####
path.to.Rmd.file <- "/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full_corrected_20240312/03_generate_GWAS_reports.Rmd"
for (maf in c(0.01, 0.05)){
  dir.create(file.path(path.to.save.HTML.dir, sprintf("GWAS_reports_snpStats_and_PLINK_maf_%s", maf)), showWarnings = FALSE, recursive = TRUE)
  for (group in all.GWAS.cases){
    print("#####------------------------------------------------#####")
    print(sprintf("Working on %s, MAF: %s", group, maf))
    print("#####------------------------------------------------#####")
    save.HTML.filename <- sprintf("GWAS_reports_snpStats_and_PLINK_maf_%s.%s.html", maf, group)
    dir.create(path.to.save.HTML.dir, showWarnings = FALSE)
    if (file.exists(file.path(path.to.save.HTML.dir, 
                              sprintf("GWAS_reports_snpStats_and_PLINK_maf_%s", maf),
                              save.HTML.filename)) == FALSE){
      rmarkdown::render(input = path.to.Rmd.file, 
                        params = list(
                          group = group,
                          maf = maf
                        ),
                        output_file = save.HTML.filename,
                        output_dir = file.path(path.to.save.HTML.dir, sprintf("GWAS_reports_snpStats_and_PLINK_maf_%s", maf))) 
    }
  }  
}

#####----------------------------------------------------------------------#####
##### DISEASE MODEL TEST RESULTS
#####----------------------------------------------------------------------#####
gc()
rm(list = ls())

source("/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full_corrected_20240312/02_generate_Rdata_for_SNPStats.R")

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "AGBruns_SNP_data_full_corrected_20240312"
path.to.save.HTML.dir <- file.path(outdir, PROJECT, "html_reports_20240312")
dir.create(path.to.save.HTML.dir, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full_corrected_20240312"
all.GWAS.cases <- read.csv(file.path(path.to.project.src, "all_GWAS_cases.txt"), sep = "\t", header = FALSE)[["V1"]]

path.to.Rmd.file <- "/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full_corrected_20240312/05_GWAS_with_models.Rmd"
for (disease.model in c("GENO", "TREND", "ALLELIC", "DOM", "REC")){
  for (maf in c(0.01, 0.05)){
    dir.create(file.path(path.to.save.HTML.dir, sprintf("GWAS_reports_disease_model_PLINK_maf_%s", maf)), showWarnings = FALSE, recursive = TRUE)
    for (group in all.GWAS.cases){
      print("#####------------------------------------------------#####")
      print(sprintf("Working on %s, MAF: %s", group, maf))
      print("#####------------------------------------------------#####")
      save.HTML.filename <- sprintf("GWAS_reports_disease_model_%s_PLINK_maf_%s.%s.html", disease.model, maf, group)
      dir.create(path.to.save.HTML.dir, showWarnings = FALSE)
      if (file.exists(file.path(path.to.save.HTML.dir, sprintf("GWAS_reports_disease_model_PLINK_maf_%s", maf), save.HTML.filename)) == FALSE){
        rmarkdown::render(input = path.to.Rmd.file, 
                          params = list(
                            group = group,
                            maf = maf,
                            disease.model = disease.model
                          ),
                          output_file = save.HTML.filename,
                          output_dir = file.path(path.to.save.HTML.dir, sprintf("GWAS_reports_disease_model_PLINK_maf_%s", maf))) 
      }
    }  
  }
}



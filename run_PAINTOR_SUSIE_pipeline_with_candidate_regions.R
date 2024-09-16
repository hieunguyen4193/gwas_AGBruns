gc()
rm(list = ls())

if ("argparse" %in% installed.packages() == FALSE){
  install.packages("argparse")
}

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full_corrected_20240312"
source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "AGBruns_SNP_data_full_corrected_20240312"

path.to.main.output <- file.path(outdir, PROJECT)

group <- "SBP_ever"
maf <- 0.01
path.to.main.src.project <- "/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full_corrected_20240312"
setwd(path.to.main.src.project)
path.to.postgwas.pipeline <- file.path(path.to.main.src.project, "04_PIPELINE_PAINTOR_AND_SUSIE.sh")
for (model in c("ALLELIC", "basic_assoc", "DOM", "GENO", "REC", "TREND")){
  path.to.candidate.region <- file.path(path.to.main.output, "candidate_regions", group, sprintf("MAF_%s", maf), sprintf("%s_MAF_%s_model_%s.tsv", group, maf, model))
  candidate.regiondf <- read.table(path.to.candidate.region)
  candidate.regiondf <- candidate.regiondf %>% rowwise() %>%
    mutate(V2 = ifelse(V2 <= 0, 0, V2))
  for (i in seq(1, nrow(candidate.regiondf))){
    tmpdf <- candidate.regiondf[i,]
    chrom <- tmpdf$V1 %>% str_replace("chr", "")
    region.start <- tmpdf$V2
    region.end <- tmpdf$V3

    path.to.pipeline.output <- file.path(path.to.main.output, "04_output", group, sprintf("MAF_%s", maf), model, sprintf("chr%s_%s_%s", chrom, region.start, region.end))
    dir.create(path.to.pipeline.output, showWarnings = FALSE, recursive = TRUE)
    
    if (file.exists(file.path(path.to.pipeline.output, 
                              sprintf("finished_region_%s_%s_%s.csv", chrom, region.start, region.end))) == FALSE){
      cmd <- sprintf("bash %s %s %s %s %s %s %s %s", path.to.postgwas.pipeline, maf, group, chrom, region.start, region.end, model, path.to.pipeline.output)
      print(sprintf("Working on MAF %s, Group %s, chrom %s, Start: %s, End: %s, disease model: %s", maf, group, chrom, region.start, region.end, model))
      system(cmd)
      write.csv(data.frame(status = c("finished_region_%s_%s_%s", chrom, region.start, region.end)), 
                file.path(path.to.pipeline.output, 
                          sprintf("finished_region_%s_%s_%s.csv", chrom, region.start, region.end)))      
    }
  }
}


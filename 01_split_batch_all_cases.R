gc()
rm(list = ls())

# stratum <- "SBP"

library(vroom)
library(comprehenr)
library(readxl)
library(writexl)
library(stringr)
library(tidyverse)
library(dplyr)

#####----------------------------------------------------------------------#####
##### PATHS AND CONFIGURATIONS
#####----------------------------------------------------------------------#####
path.to.main.dir <- "/home/hieunguyen/CRC1382"
PROJECT <- "AGBruns_SNP_data_full_corrected_20240312"
path.to.outdir <- "/home/hieunguyen/CRC1382/outdir"

path.to.storage <- "/home/hieunguyen/CRC1382/storage"
main.data.dir <- file.path(path.to.storage, "AGBruns_SNP_data_full")

path.to.main.output <- file.path(path.to.outdir, PROJECT)
path.to.00.output <- file.path(path.to.main.output, "00_output")
path.to.01.output <- file.path(path.to.main.output, "01_output")

dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- file.path("/home/hieunguyen/CRC1382/src_2023", PROJECT)
#####----------------------------------------------------------------------#####
##### PREPROCESSING INPUT SNP TABLES
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.01.output, "snp_table.map")) == FALSE){
  print("Generating snp_table.map")
  ##### reduce the SNP table
  all.snp.tables <- Sys.glob(file.path(main.data.dir, "A5057", "*", "*_SNP Table*"))
  
  snp.table.1 <- read.csv(all.snp.tables[[1]], sep = "\t")
  reduced.snp.table1 <- subset(snp.table.1, select = c("Index", "Name", "Chr", "Position",
                                                       "ChiTest100", "Het.Excess", "AA.Freq", "AB.Freq", "BB.Freq", 
                                                       "Minor.Freq", "Call.Freq",  "SNP"))
  
  snp.table.1 <- subset(snp.table.1, grepl("D", snp.table.1$SNP) == FALSE)
  snp.table.1 <- subset(snp.table.1, grepl("I", snp.table.1$SNP) == FALSE)
  
  snp.table.2 <- read.csv(all.snp.tables[[2]], sep = "\t")
  reduced.snp.table2 <- subset(snp.table.2, select = c("Index", "Name", "Chr", "Position",
                                                       "ChiTest100", "Het.Excess", "AA.Freq", "AB.Freq", "BB.Freq", 
                                                       "Minor.Freq", "Call.Freq",  "SNP"))
  
  snp.table.2 <- subset(snp.table.2, grepl("D", snp.table.2$SNP) == FALSE)
  snp.table.2 <- subset(snp.table.2, grepl("I", snp.table.2$SNP) == FALSE)
  
  ##### Generate a .map file
  mapdf <- subset(snp.table.1, select = c(Chr, Name))
  mapdf$Pos.morgans <- 0
  mapdf$coord <- snp.table.1$Position
  
  ##### UPDATE 16.10.2023
  # REMOVE CHROMOSOME X, XY, Y, MT from mapdf
  mapdf <- subset(mapdf, mapdf$Chr %in% c("MT", "Y", "X", "XY") == FALSE)
  
  write.csv(mapdf, file = file.path(path.to.01.output, "snp_table.map"))
  write.table(mapdf, file = file.path(path.to.01.output, "snp_table_for_GWAS.map"), sep = "\t", col.names = FALSE, row.names = FALSE)
  write.csv(snp.table.1, file = file.path(path.to.01.output, "snp_table_full_info_plate1.map"))
  write.csv(snp.table.2, file = file.path(path.to.01.output, "snp_table_full_info_plate2.map"))
  write.csv(subset(snp.table.1, select = c(Name, SNP)), file = file.path(path.to.01.output, "snp_table_full_info.map"))
} else {
  print("SNP_table.map existed!")
  mapdf <- read.csv(file.path(path.to.01.output, "snp_table.map"))
  if ("X" %in% colnames(mapdf)){
    mapdf <- subset(mapdf, select = -c(X))
  }
}

#####----------------------------------------------------------------------#####
##### PREPROCESSING INPUT DATA TABLES
#####----------------------------------------------------------------------#####
all.sample.tables <- Sys.glob(file.path(main.data.dir, "A5057", "*", "*_Samples Table*"))

sample.table <- data.frame()
for (file in all.sample.tables){
  tmp <- read.csv(file, sep = "\t")
  sample.table <- rbind(sample.table, tmp)
}

sample.table <- subset(sample.table, select = c("Unique.sample.A.number", "Sample.ID", "Gender.Est"))

colnames(sample.table) <- c("unique_ID", "Sample.ID", "Gender")

sample.table <- sample.table %>% rowwise() %>% mutate(Sample.ID = str_replace(Sample.ID, "-1", ".1"))

##### UPDATE 08.10.2023
full.meta.data <- read_excel(file.path(path.to.storage, "AGBruns_SNP_data_full", "20231007_GWAS_Clinical_Data.xlsx"))

all_gwas_cases <- read.csv(file.path(path.to.main.src, "all_GWAS_cases.txt"), sep = "\t", header = FALSE)[["V1"]]

for (stratum in all_gwas_cases){
# for (stratum in c("SBP_ever")){
  meta.data <- full.meta.data
  if (file.exists(file.path(path.to.01.output, sprintf("parallel_processing_%s", stratum))) == FALSE){
    path.to.parallel.workdir <- file.path(path.to.01.output, sprintf("parallel_processing_%s", stratum))
    dir.create(path.to.parallel.workdir, showWarnings = FALSE, recursive = TRUE)
    path.to.input.chunks <- file.path(path.to.parallel.workdir, "input_chunks")
    dir.create(path.to.input.chunks, showWarnings = FALSE, recursive = TRUE)
    
    meta.data <- meta.data[c("eurofins_ID", stratum)]
    
    colnames(meta.data) <- c("unique_ID", stratum)
    
    meta.data <- subset(meta.data, meta.data$unique_ID %in% unique(sample.table$unique_ID))
    merge.sample.table <- merge(sample.table, meta.data, by.x = "unique_ID", by.y = "unique_ID", all.x = TRUE)
    
    merge.sample.table <- merge.sample.table[!duplicated(merge.sample.table$unique_ID),]
    
    peddf <- merge.sample.table[,c("Sample.ID", "Gender", stratum)]
    
    ##### assuming missing values = no disease!!! ---> WRONG
    # peddf <- replace(peddf, is.na(peddf), 0)
    
    ##### UPDATE 13.03.2024: REMOVE NA 
    peddf <- na.omit(peddf)

    peddf$family_id <- to_vec(for (item in seq(1, nrow(peddf))) sprintf("fam_%s", item))
    peddf$father_id <- 0
    peddf$mother_id <- 0
    # peddf <- subset(peddf, select = c(family_id, Sample.ID, father_id, mother_id, Gender, SBP))
    peddf <- peddf[, c("family_id", "Sample.ID", "father_id", "mother_id", "Gender", stratum)]  
    
    all.full.tables <- Sys.glob(file.path(path.to.00.output, "*.reduced.txt"))
    
    full.table.df <- data.frame()
    for (file in all.full.tables){
      tmp.full.table <- read.csv(file, sep = ",") %>% column_to_rownames("X")
      tmp.full.table <- t(tmp.full.table)
      row.names(tmp.full.table) <- to_vec(for (item in row.names(tmp.full.table)) 
        str_replace_all(item, ".GType", ""))
      full.table.df <- rbind(full.table.df, tmp.full.table)
    }
    
    full.table.df <- full.table.df[, mapdf$Name]
    
    row.names(full.table.df) <- to_vec(for (item in row.names(full.table.df)) str_replace(item, "X", ""))
    
    peddf <- merge(peddf, full.table.df %>% rownames_to_column("Sample.ID"), by.x = "Sample.ID", by.y = "Sample.ID")
    
    peddf <- peddf[, c(c("Sample.ID", "family_id", "father_id", "mother_id", "Gender", stratum), colnames(peddf)[7:length(colnames(peddf))])]
    
    
    #####----------------------------------------------------------------------#####
    ##### Split the current peddf file to smaller chunks for faster processing
    #####----------------------------------------------------------------------#####
    all.mutations <- colnames(peddf)[7:ncol(peddf)]
    
    split.all.mutations <- split(all.mutations, ceiling(seq_along(all.mutations) / 10000))
    
    path.to.output.chunks <- file.path(path.to.parallel.workdir, "output_chunks")
    dir.create(path.to.output.chunks, showWarnings = FALSE, recursive = TRUE)
    
    for (group in names(split.all.mutations)){
      print(sprintf("Working on %s", group))
      selected.cols <- c(colnames(peddf)[1:6], split.all.mutations[[group]])
      tmp.peddf <- peddf[, selected.cols]
      
      write.csv(tmp.peddf, file.path(path.to.input.chunks, sprintf("chunk_%s.csv", group)))
    }
  }

}

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
library(hash)

#####----------------------------------------------------------------------#####
##### PATHS AND CONFIGURATIONS
#####----------------------------------------------------------------------#####

selected.mutations <- c("rs7085647", "rs1875051", "rs9987289")
save.name <- paste(selected.mutations, collapse = "_")

path.to.main.dir <- "/home/hieunguyen/CRC1382"
PROJECT <- "AGBruns_SNP_data_full_corrected_20240312"
path.to.outdir <- "/home/hieunguyen/CRC1382/outdir"

path.to.storage <- "/home/hieunguyen/CRC1382/storage"
main.data.dir <- file.path(path.to.storage, "AGBruns_SNP_data_full")

path.to.main.output <- file.path(path.to.outdir, PROJECT)
path.to.00.output <- file.path(path.to.main.output, "00_output")
path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.save.output <- file.path(path.to.main.output, "get_data_from_input_SNP")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

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
  if (file.exists(file.path(path.to.save.output, sprintf("%s_selected_mutation_%s", stratum, save.name))) == FALSE){
    print(sprintf("Working on case %s", stratum))
    meta.data <- full.meta.data
    meta.data <- meta.data[c("eurofins_ID", stratum)]
    
    colnames(meta.data) <- c("unique_ID", stratum)
    
    meta.data <- subset(meta.data, meta.data$unique_ID %in% unique(sample.table$unique_ID))
    merge.sample.table <- merge(sample.table, meta.data, by.x = "unique_ID", by.y = "unique_ID", all.x = TRUE)
    
    merge.sample.table <- merge.sample.table[!duplicated(merge.sample.table$unique_ID),]
    
    peddf <- merge.sample.table[,c("Sample.ID", "Gender", stratum)]
    
    ##### assuming missing values = no disease!!!
    # peddf <- replace(peddf, is.na(peddf), 0)
    peddf <- na.omit(peddf)
    peddf$family_id <- to_vec(for (item in seq(1, nrow(peddf))) sprintf("fam_%s", item))
    peddf$father_id <- 0
    peddf$mother_id <- 0
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
    
    all.mutations <- colnames(peddf)[7:ncol(peddf)]
    
    selected.cols <- c(colnames(peddf)[1:6], selected.mutations)
    tmp.peddf <- peddf[, selected.cols]
    write.csv(tmp.peddf, file.path(path.to.save.output, sprintf("%s_selected_mutation_%s", stratum, save.name)))
  } 
}
  
alldf <- hash()
for (stratum in all_gwas_cases){
  alldf[[stratum]] <- read.csv(file.path(path.to.save.output, sprintf("%s_selected_mutation_%s", stratum, save.name)))
}
  
maindf <- alldf$AlcoholCirrhosis[, c(2:7)]
for (stratum in all_gwas_cases){
  if (stratum != "AlcoholCirrhosis"){
    tmpdf <- alldf[[stratum]][, c(2, 7)]
    maindf <- merge(maindf, tmpdf, by.x = "Sample.ID", by.y = "Sample.ID")
  }
}

maindf <- merge(maindf, alldf$AlcoholCirrhosis[, c("Sample.ID", selected.mutations)], by.x = "Sample.ID", by.y = "Sample.ID")
convert.gt <- list(AA = 0, AB = 1, BB = 2)
for (mutation in selected.mutations){
  maindf[[sprintf("convert_%s", mutation)]] <- unlist(lapply(maindf[[mutation]], function(x){
    return(convert.gt[[x]])
  }))
}

for (mutation in selected.mutations){
  all.outputdf <- data.frame()
  for (stratum in all_gwas_cases){
    tmpdf <- maindf[, c("Sample.ID", stratum, mutation)]
    tmpdf.stratum <- tmpdf[tmpdf[[stratum]] == 1,]
    
    outputdf <- data.frame(data = c(stratum))
    colnames(outputdf) <- c(mutation)
    
    count.in.all.patients <- table(tmpdf[[mutation]])
    count.in.positive.patients <- table(tmpdf.stratum[[mutation]])
    for (item in c("AA", "AB", "BB")){
      if (item %in% names(count.in.all.patients) == FALSE){
        count.in.all.patients[[item]] <- 0
      }
    }
    
    for (item in c("AA", "AB", "BB")){
      if (item %in% names(count.in.positive.patients) == FALSE){
        count.in.positive.patients[[item]] <- 0
      }
    }
    
    outputdf[[sprintf("All patients (N = %s)", nrow(tmpdf))]] <- sprintf("%s (%s %%)", nrow(tmpdf.stratum), format(round(100 * nrow(tmpdf.stratum)/nrow(tmpdf), digits=2),nsmall=2))
    outputdf[[sprintf("AA (N = %s)", count.in.all.patients[["AA"]])]] <- sprintf("%s (%s %%)", count.in.positive.patients[["AA"]], format(round(100 * count.in.positive.patients[["AA"]]/count.in.all.patients[["AA"]], digits=2),nsmall=2))
    outputdf[[sprintf("AB (N = %s)", count.in.all.patients[["AB"]])]] <- sprintf("%s (%s %%)", count.in.positive.patients[["AB"]], format(round(100 * count.in.positive.patients[["AB"]]/count.in.all.patients[["AB"]], digits=2),nsmall=2))
    outputdf[[sprintf("BB (N = %s)", count.in.all.patients[["BB"]])]] <- sprintf("%s (%s %%)", count.in.positive.patients[["BB"]], format(round(100 * count.in.positive.patients[["BB"]]/count.in.all.patients[["BB"]] , digits=2),nsmall=2))
    all.outputdf <- rbind(all.outputdf, outputdf)  
  }
  
  writexl::write_xlsx(all.outputdf, file.path(path.to.save.output, sprintf("summary_table_SNP_%s.xlsx", mutation)), format_headers = TRUE)
}

writexl::write_xlsx(maindf, file.path(path.to.save.output, sprintf("raw_data_frame_SNPs_%s.xlsx", save.name)), format_headers = TRUE)

maindf.eurofins <- merge(maindf, subset(sample.table, select = c(unique_ID, Sample.ID)), by.x = "Sample.ID", by.y = "Sample.ID")
writexl::write_xlsx(maindf.eurofins, file.path(path.to.save.output, sprintf("raw_data_frame_SNPs_%s.with_eurofins.xlsx", save.name)), format_headers = TRUE)

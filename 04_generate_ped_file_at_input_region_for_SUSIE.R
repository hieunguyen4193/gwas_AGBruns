gc()
rm(list = ls())

library(vroom)
library(comprehenr)
library(readxl)
library(writexl)
library(stringr)
library(tidyverse)
library(dplyr)
library(argparse)
library(hash)

#####----------------------------------------------------------------------#####
##### ARGS PARSER
#####----------------------------------------------------------------------#####
parser <- ArgumentParser()

parser$add_argument("-s", "--stratum", type="character",
                    help="Choose the clinical condition")
parser$add_argument("-c", "--chrom", type="character",
                    help="Chromosome of the region")
parser$add_argument("-x", "--start", type="character",
                    help="Starting coordinate of the region")
parser$add_argument("-y", "--end", type="character",
                    help="Ending coordinate of the region")
parser$add_argument("-o", "--output_dir", type="character",
                    help="Path to save output")
args <- parser$parse_args()

stratum <- args$stratum
chrom <- args$chrom
region.start <- args$start
region.end <- args$end 
path.to.04.output <- args$output_dir

#####----------------------------------------------------------------------#####
##### PATHS AND CONFIGURATIONS
#####----------------------------------------------------------------------#####
PROJECT <- "AGBruns_SNP_data_full_corrected_20240312"
path.to.outdir <- "/home/hieunguyen/CRC1382/outdir"

path.to.storage <- "/home/hieunguyen/CRC1382/storage"
main.data.dir <- file.path(path.to.storage, "AGBruns_SNP_data_full")

path.to.main.output <- file.path(path.to.outdir, PROJECT)
path.to.00.output <- file.path(path.to.main.output, "00_output")
# path.to.04.output <- file.path(path.to.main.output, "04_output", stratum, sprintf("region_chr%s_%s_%s", chrom, region.start, region.end))

dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- file.path("/home/hieunguyen/CRC1382/src_2023", PROJECT)
#####----------------------------------------------------------------------#####
##### PREPROCESSING INPUT SNP TABLES
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.04.output, sprintf("region_chr%s_%s_%s.raw.map", chrom, region.start, region.end))) == FALSE){
  print(sprintf("Generating snp_table.map at region chromosome %s, start = %s, end = %s", chrom, region.start, region.end))
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
  mapdf <- subset(snp.table.1, select = c(Chr, Name, SNP))
  mapdf$Pos.morgans <- 0
  mapdf$coord <- snp.table.1$Position
  
  ##### UPDATE 16.10.2023
  # REMOVE CHROMOSOME X, XY, Y, MT from mapdf
  mapdf <- subset(mapdf, mapdf$Chr %in% c("MT", "Y", "X", "XY") == FALSE)
  
  mapdf <- subset(mapdf, (mapdf$Chr == chrom) & (mapdf$coord >= as.numeric(region.start) - 1) & (mapdf$coord <= as.numeric(region.end) + 1))
  
  write.csv(mapdf, file = file.path(path.to.04.output, sprintf("region_chr%s_%s_%s.raw.map", chrom, region.start, region.end)))
} else {
  print(sprintf("Found snp_table.map at region chromosome %s, start = %s, end = %s", chrom, region.start, region.end))
  print("reading in ...")
  mapdf <- read.csv(file.path(path.to.04.output, sprintf("region_chr%s_%s_%s.raw.map", chrom, region.start, region.end)))
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

meta.data <- full.meta.data

meta.data <- meta.data[c("eurofins_ID", stratum)]

colnames(meta.data) <- c("unique_ID", stratum)

meta.data <- subset(meta.data, meta.data$unique_ID %in% unique(sample.table$unique_ID))
merge.sample.table <- merge(sample.table, meta.data, by.x = "unique_ID", by.y = "unique_ID", all.x = TRUE)

merge.sample.table <- merge.sample.table[!duplicated(merge.sample.table$unique_ID),]

peddf <- merge.sample.table[,c("Sample.ID", "Gender", stratum)]

##### assuming missing values = no disease!!!
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

peddf <- subset(peddf, peddf[[stratum]] %in% c(0, 1))

peddf[, stratum] <- unlist(lapply(peddf[, stratum], function(x){
  if (x == 0){
    return(1)
  } else if (x == 1){
    return(2)
  }
}))

generate_variant_genotype_info <- function(mutation){
  gt <- subset(mapdf, mapdf$Name == mutation)$SNP
  allele1 <- str_replace_all(str_split(gt, "/")[[1]][[1]], "\\[", "")
  allele2 <- str_replace_all(str_split(gt, "/")[[1]][[2]], "\\]", "")
  
  convert_allele <- hash()
  convert_allele[["A"]] <- 1
  convert_allele[["C"]] <- 2
  convert_allele[["G"]] <- 3
  convert_allele[["T"]] <- 4
  
  convert_genotype_to_SNPCall <- function(allele1, allele2, genotype){
    allele1 <- convert_allele[[allele1]]
    allele2 <- convert_allele[[allele2]]
    
    if (genotype == "AA"){
      return(c(allele1, allele1))
    } else if(genotype == "BB"){
      return(c(allele2, allele2))
    } else if (genotype == "AB"){
      return(c(allele1, allele2))
    } else if (genotype == "NC"){
      return(c(0, 0))
    }
  }
  
  col1 <- unlist(lapply(peddf[[mutation]], function(x){
    return(convert_genotype_to_SNPCall(allele1, allele2, x)[[1]])
  }))
  
  col2 <- unlist(lapply(peddf[[mutation]], function(x){
    return(convert_genotype_to_SNPCall(allele1, allele2, x)[[2]])
  }))
  return(list(allele.col1 = col1, allele.col2 = col2))
}

all_mutation_names <- colnames(peddf)[7: ncol(peddf)]
output.peddf <- peddf[, 1:6]

for (i in seq(length(all_mutation_names))){
  mutation <- all_mutation_names[[i]]
  
  res <- generate_variant_genotype_info(mutation)
  
  col1 <- res$allele.col1
  col2 <- res$allele.col2
  
  output.peddf[[sprintf("%s_1", mutation)]] <- col1
  output.peddf[[sprintf("%s_2", mutation)]] <- col2
  
  if (i %% 100 == 0){
    print(sprintf("Working at the %s-th mutation", i))
  }
}

write.table(output.peddf, sep = "\t", file = file.path(path.to.04.output, sprintf("region_chr%s_%s_%s.ped", chrom, region.start, region.end)), col.names = FALSE, row.names = FALSE)

dim(output.peddf)

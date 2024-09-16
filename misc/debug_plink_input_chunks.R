gc()
rm(list = ls())


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
path.to.main.dir <- "/home/hieunguyen/CRC1382"

path.to.outdir <- file.path(path.to.main.dir, "outdir")

main.data.dir <- file.path(path.to.outdir, "AGBruns_SNP_data")

path.to.save.output <- file.path(path.to.outdir, "AGBruns_SNP_data_full")

path.to.01.output <- file.path(path.to.save.output, "01_output")

path_to_input_chunks <- "/home/hieunguyen/CRC1382/outdir/AGBruns_SNP_data_full/01_output/parallel_processing_history_of_cancer"

stratum <- "history_of_cancer"

input.chunks <- Sys.glob(file.path(path_to_input_chunks, "input_chunks" , "*"))

path.to.input.ped.chunk <- input.chunks[[1]]

mapdf <- read.csv(file.path(path.to.01.output, "snp_table_full_info.map")) 
#####----------------------------------------------------------------------#####
# MAIN ANALYSIS
#####----------------------------------------------------------------------#####
peddf <- read.csv(path.to.input.ped.chunk, check.names = FALSE, )
peddf <- peddf[, 2:ncol(peddf)]

# (1=unaff, 2=aff, 0=miss)
# Change stratum condition coding

peddf <- subset(peddf, peddf[[stratum]] %in% c(0, 1))

peddf[, stratum] <- unlist(lapply(peddf[, stratum], function(x){
  if (x == 0){
    return(1)
  } else if (x == 1){
    return(2)
  }
}))

##### TEST ONLY
# peddf <- peddf[, 1:20]
# colnames(peddf) <- to_vec(for (item in colnames(peddf)) str_replace(str_replace(item, "X", ""), "\\.", ":"))

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
official.peddf <- peddf[, 1:6]

for (i in seq(length(all_mutation_names))){
  mutation <- all_mutation_names[[i]]
  
  res <- generate_variant_genotype_info(mutation)
  
  col1 <- res$allele.col1
  col2 <- res$allele.col2
  
  official.peddf[[sprintf("%s_1", mutation)]] <- col1
  official.peddf[[sprintf("%s_2", mutation)]] <- col2
  
  if (i %% 100 == 0){
    print(sprintf("Working at the %s-th mutation", i))
  }
}
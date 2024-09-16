gc()
rm(list = ls())

library(stringr)
library(tidyverse)
library(dplyr)
library(comprehenr)

path.to.storage <- "/home/hieunguyen/CRC1382/storage"
path.to.maindir <- file.path(path.to.storage, "AGBruns_SNP_data_full", "A5057")
all.full.tables <- Sys.glob(file.path(path.to.maindir, "*", "*Full Data Table*"))

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "AGBruns_SNP_data_full_corrected_20240312"
path.to.main.output <- file.path(outdir, PROJECT)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.00.output <- file.path(path.to.main.output, "00_output")
dir.create(path.to.00.output, showWarnings = FALSE, recursive = TRUE)

for (file in all.full.tables){
  filename <- basename(file)
  full.table <- read.csv(file, sep = "\t")
  cols.to.keep <- c("Name", to_vec(for (item in colnames(full.table)) if(grepl("GType", item) == TRUE) item))
  write.csv(data.frame(data = cols.to.keep), file.path(path.to.00.output, str_replace(filename, ".txt", ".sampleIDonly.csv")))
  reduced.full.table <- full.table[, cols.to.keep] %>%#
    column_to_rownames("Name")
  write.csv(reduced.full.table, file.path(path.to.00.output, str_replace(filename, ".txt", ".reduced.txt")))
}


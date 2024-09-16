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
parser$add_argument("-m", "--maf", type="character",
                    help="Input Minor Allele Frequency")
parser$add_argument("-d", "--model", type="character",
                    help="Name of the association (disease) model")
parser$add_argument("-o", "--output_dir", type="character",
                    help="Path to save output html")

args <- parser$parse_args()

stratum <- args$stratum
chrom <- args$chrom
region.start <- args$start
region.end <- args$end
maf <- args$maf
model <- args$model
path.to.04.output <- args$output_dir

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full_corrected_20240312"

PROJECT <- "AGBruns_SNP_data_full_corrected_20240312"
path.to.outdir <- "/home/hieunguyen/CRC1382/outdir"

path.to.storage <- "/home/hieunguyen/CRC1382/storage"
main.data.dir <- file.path(path.to.storage, "AGBruns_SNP_data_full")

path.to.main.output <- file.path(path.to.outdir, PROJECT)
# path.to.04.output <- file.path(path.to.main.output, "04_output", stratum, sprintf("region_chr%s_%s_%s", chrom, region.start, region.end))

path.to.Rmd.file <- file.path(path.to.project.src, "04_SUSIE_and_PAINTOR_fine_mapping.Rmd")

save.HTML.filename <- sprintf("PostGWAS_SUSIE_%s_region_chr%s_%s_%s_MAF_%s.model_%s.html", stratum, chrom, region.start, region.end, maf, model)
path.to.save.HTML.dir <- file.path(path.to.04.output)

print(path.to.save.HTML.dir)
rmarkdown::render(input = path.to.Rmd.file, 
                    params = list(
                      maf = maf,
                      stratum = stratum,
                      chrom = chrom,
                      region.start = region.start,
                      region.end = region.end,
                      path.to.04.output = path.to.04.output
                    ),
                    output_file = save.HTML.filename,
                    output_dir = path.to.save.HTML.dir)  

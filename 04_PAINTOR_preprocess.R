gc()
rm(list = ls())

library(argparse)
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
                    help="Minor allele frequency")
parser$add_argument("-o", "--output_dir", type="character",
                    help="Path to output dir")
args <- parser$parse_args()

stratum <- args$stratum
chrom <- args$chrom
region.start <- args$start
region.end <- args$end
maf <- args$maf
path.to.04.output <- args$output_dir

PROJECT <- "AGBruns_SNP_data_full_corrected_20240312"
path.to.outdir <- "/home/hieunguyen/CRC1382/outdir"

path.to.project.src <- file.path("/home/hieunguyen/CRC1382/src_2023", PROJECT)

source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

path.to.storage <- "/home/hieunguyen/CRC1382/storage"
main.data.dir <- file.path(path.to.storage, "AGBruns_SNP_data_full")

path.to.main.output <- file.path(path.to.outdir, PROJECT)
path.to.00.output <- file.path(path.to.main.output, "00_output")
# path.to.04.output <- file.path(path.to.main.output, "04_output", stratum, sprintf("region_chr%s_%s_%s", chrom, region.start, region.end))
dir.create(file.path(path.to.04.output, "PAINTOR"), showWarnings = FALSE, recursive = TRUE)


mapdf <- read.csv(file.path(path.to.04.output, sprintf("region_chr%s_%s_%s.raw.map", chrom, region.start, region.end))) %>%
  subset(select = -c(X))
ld.matrix <- read.table(file.path(path.to.04.output, "plink.r.matrix.ld"), header = FALSE) 

colnames(ld.matrix) <- mapdf$Name
row.names(ld.matrix) <- mapdf$Name

path.to.plink.output <- file.path(path.to.outdir, PROJECT, "PLINK_results", sprintf("maf_%s", maf), sprintf("%s_results_maf_%s", stratum, maf))

assocdf <- read.table(file.path(path.to.plink.output, "plink.assoc"), header = TRUE) %>%
  subset(SNP %in% mapdf$Name) %>%
  rowwise() %>%
  mutate(logOR = log10(OR)) %>%
  mutate(signOR = ifelse(logOR >= 0, 1, -1)) %>% 
  mutate(zscore = signOR * qnorm(0.5 * P))

assocdf$CHR <- unlist(lapply(
  assocdf$SNP, function(x){
    return(subset(mapdf, mapdf$Name == x)$Chr)
  }
))

ld.matrix <- ld.matrix[assocdf$SNP, assocdf$SNP]

assocdf <- assocdf[, c("CHR", "BP", "SNP", "zscore")]

colnames(assocdf) <- c("CHR", "POS", "RSID", "ZSCORE.P1")

write.table(assocdf, row.names = FALSE, sep = " ", file = file.path(path.to.04.output, "PAINTOR", "Locus"), quote = FALSE)
write.table(ld.matrix, col.names = FALSE, row.names = FALSE, sep = " ", file = file.path(path.to.04.output, "PAINTOR", "Locus.LD"), quote = FALSE)
write.table(data.frame(data=c("Locus")), row.names = FALSE, col.names = FALSE, file = file.path(path.to.04.output, "PAINTOR", "input.file"), quote = FALSE)

fake.annotations <- to_vec(for (item in seq(1, nrow(assocdf)))1)
write.table(data.frame(fake_annotations=fake.annotations), row.names = FALSE, file = file.path(path.to.04.output, "PAINTOR", "Locus.annotations"), quote = FALSE)




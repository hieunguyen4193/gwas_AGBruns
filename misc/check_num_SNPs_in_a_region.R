chrom <- 16
region.start <- 50693606
region.end <- 50733075

path.to.region.data <- sprintf("/media/hieunguyen/HD0/outdir/CRC1382/AGBruns_SNP_data_full_removeXY_20231018/04_output/SBP_ever/region_chr%s_%s_%s", chrom, region.start, region.end)
#####----------------------------------------------------------------------#####
##### PATHS AND CONFIGURATIONS
#####----------------------------------------------------------------------#####
PROJECT <- "AGBruns_SNP_data_full_removeXY_20231018"
path.to.outdir <- "/media/hieunguyen/HD0/outdir/CRC1382"

path.to.storage <- "/home/hieunguyen/CRC1382/storage"
main.data.dir <- file.path(path.to.storage, "AGBruns_SNP_data_full")

path.to.main.output <- file.path(path.to.outdir, PROJECT)
path.to.00.output <- file.path(path.to.main.output, "00_output")
path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full"
#####----------------------------------------------------------------------#####
##### PREPROCESSING INPUT SNP TABLES
#####----------------------------------------------------------------------#####
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

mapdf.original <- mapdf
mapdf <- subset(mapdf, (mapdf$Chr == chrom) & (mapdf$coord >= as.numeric(region.start) - 1) & (mapdf$coord <= as.numeric(region.end) + 1))


peddf <- read.table(file.path(path.to.region.data, "region_chr16_50693606_50733075.ped"))

ld.matrix <- read.table(file.path(path.to.region.data, "plink.r2.matrix.ld"))
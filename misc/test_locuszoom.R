gc()
rm(list = ls())

maf <- 0.01
stratum <- "SBP_ever"
chrom <- 8
region.start <- "8000000"
region.end <- "9900000" 

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full"
source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

locuszoom.dir <- file.path(path.to.project.src, "LocusZooms")

source(file.path(locuszoom.dir, "functions", "locus_zoom.R"))

PROJECT <- "AGBruns_SNP_data_full_removeXY_20231018"
path.to.outdir <- "/media/hieunguyen/HD0/outdir/CRC1382"

path.to.storage <- "/home/hieunguyen/CRC1382/storage"
main.data.dir <- file.path(path.to.storage, "AGBruns_SNP_data_full")

path.to.main.output <- file.path(path.to.outdir, PROJECT)
path.to.00.output <- file.path(path.to.main.output, "00_output")
path.to.04.output <- file.path(path.to.main.output, "04_output", stratum, sprintf("region_chr%s_%s_%s", chrom, region.start, region.end))

Unique.genes <- read.delim(file.path(locuszoom.dir, "Gencode_GRCh37_Genes_UniqueList2021.txt"), stringsAsFactors = FALSE, header = TRUE)
#####----------------------------------------------------------------------#####
##### HELPER FUNCTIONS
#####----------------------------------------------------------------------#####
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}

#####----------------------------------------------------------------------#####
chosen.snp <- "rs9987289"

mapdf <- read.table(file.path(path.to.04.output, sprintf("region_chr%s_%s_%s.raw.map", chrom, region.start, region.end)), sep = ",", header = TRUE) %>%
  subset(select = -c(X))

path.to.ld.matrix <- file.path(path.to.04.output, "plink.r2.matrix.ld")
path.to.assoc.file <- file.path(path.to.outdir, PROJECT, "PLINK_results", sprintf("%s_results_maf_%s", stratum, maf), "plink.assoc")
assocdf <- read.table(path.to.assoc.file, stringsAsFactors = FALSE, header = TRUE)
ld.matrix <- read.table(path.to.ld.matrix, stringsAsFactors = FALSE, header = FALSE)

colnames(ld.matrix) <- mapdf$Name
row.names(ld.matrix) <- mapdf$Name

ld.matrix <- data.frame(R2 = ld.matrix[, c(chosen.snp)])
ld.matrix$SNP_A <- chosen.snp
ld.matrix$SNP_B <- mapdf$Name
ld.matrix <- subset(ld.matrix, ld.matrix$SNP_B != chosen.snp)

input.ld.matrix <- ld.matrix

assocdf <- assocdf %>%
  subset(SNP %in% mapdf$Name)

assocdf$CHR <- unlist(
  lapply(assocdf$SNP, function(x){
    return(subset(mapdf, mapdf$Name == x)$Chr)
  })
)

input.ld.matrix$CHR_A <- unlist(
  lapply(input.ld.matrix$SNP_A, function(x){
    return(subset(mapdf, mapdf$Name == x)$Chr)
  })
)

input.ld.matrix$CHR_B <- unlist(
  lapply(input.ld.matrix$SNP_B, function(x){
    return(subset(mapdf, mapdf$Name == x)$Chr)
  })
)

input.ld.matrix$BP_A <- unlist(
  lapply(input.ld.matrix$SNP_A, function(x){
    return(subset(mapdf, mapdf$Name == x)$BP)
  })
)

input.ld.matrix$BP_B <- unlist(
  lapply(input.ld.matrix$SNP_B, function(x){
    return(subset(mapdf, mapdf$Name == x)$BP)
  })
)

locus.zoom(data = assocdf,                                    
           region = c(chrom, as.numeric(region.start), as.numeric(region.end)),                             
           offset_bp = 0,                                                  
           ld.file = input.ld.matrix,                                           
           genes.data = Unique.genes,			                   
           plot.title = "Check",        
           file.name = file.path(path.to.project.src, sprintf("LocusZoom_%s_chr%s_%s_%s.jpg", stratum, chrom, region.start, region.end)),                                   
           secondary.label = TRUE,
           rsid.check = FALSE)   



